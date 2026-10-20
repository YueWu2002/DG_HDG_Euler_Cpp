#include "euler/hybrid.hpp"
namespace euler {
HybridOperator::HybridOperator(const Problem &p)
    : p_(p),volume_(p),faces_jacobian_(p),lu_(p.layout.block_size()),
      x_(static_cast<Eigen::Index>(p.layout.block_size())*p.local_trace_size()*p.layout.cells),
      schur_(static_cast<Eigen::Index>(p.local_trace_size())*p.local_trace_size()*p.layout.cells),
      y_(p.layout.block_size(),p.layout.cells),a_(p.layout.block_size(),p.layout.block_size()),
      local_rhs_(p.layout.block_size(),p.local_trace_size()+1),solved_(local_rhs_.rows(),local_rhs_.cols()),
      face_action_(p.local_trace_size(),p.local_trace_size()+1),
      volume_residual_(p.layout.size()),skeleton_residual_(p.trace_size()),rhs_(p.trace_size()),
      local_(p.layout.block_size()),faces_(p.local_trace_size()),face_result_(p.local_trace_size()),
      face_flux_(p.basis.nodes,NVAR) {}
void HybridOperator::gather(const Vec &input,int cell) const {
    const auto fs=p_.faces(cell); const int fw=p_.face_width();
    for(int s=0;s<NFACES;++s) faces_.segment(s*fw,fw)=input.segment(fs[s]*fw,fw);
}
void HybridOperator::scatter(Vec &output,int cell,const Vec &local) const {
    const auto fs=p_.faces(cell); const int fw=p_.face_width();
    for(int s=0;s<NFACES;++s) output.segment(fs[s]*fw,fw)+=local.segment(s*fw,fw);
}
void HybridOperator::residual(const Vec &state,const Vec &lambda,Vec &r,Vec &g) const {
    p_.hybrid_residual(state,lambda,r,g);
}
void HybridOperator::refresh(const Vec &state,const Vec &lambda) {
    residual(state,lambda,volume_residual_,skeleton_residual_);
    condense(state);
}
void HybridOperator::refresh_projected(const Vec &state,const Vec &r) {
    volume_residual_=r; skeleton_residual_.setZero(); condense(state);
}
void HybridOperator::assemble_cell(int cell) {
    const int n=p_.basis.nodes, fw=p_.face_width(), t=p_.local_trace_size();
    volume_.block(cell,a_);
    local_rhs_.leftCols(t).setZero();
    const auto fs=p_.faces(cell);
    for(int side=0;side<NFACES;++side) {
        const Vec &shape=side%2?p_.basis.right:p_.basis.left;
        for(int q=0;q<n;++q) {
            const auto derivative=faces_jacobian_.hdg(fs[side],q,side%2!=0);
            const double stabilization=p_.face_weight(fs[side],q)*faces_jacobian_.tau(fs[side],q,side%2!=0);
            for(int v=0;v<NVAR;++v) for(int i=0;i<n;++i)
                local_rhs_(p_.layout.local(v,p_.face_node(side,q,i)),side*fw+v*n+q)=-stabilization*shape[i];
            // Column-major assembly: contiguous x lines; y lines retain their n stride.
            for(int w=0;w<NVAR;++w) for(int j=0;j<n;++j) {
                auto column=a_.col(p_.layout.local(w,p_.face_node(side,q,j)));
                for(int v=0;v<NVAR;++v) {
                    const double scale=shape[j]*derivative(v,w);
                    if(side<2) column.segment(p_.layout.local(v,n*q),n).noalias()+=scale*shape;
                    else for(int i=0;i<n;++i) column[p_.layout.local(v,p_.face_node(side,q,i))]+=scale*shape[i];
                }
            }
        }
    }
    local_rhs_.col(t)=volume_residual_.segment(cell*p_.layout.block_size(),p_.layout.block_size());
}
void HybridOperator::condense(const Vec &state) {
    volume_.refresh(state); faces_jacobian_.refresh(state); rhs_=-skeleton_residual_;
    const int m=p_.layout.block_size(), t=p_.local_trace_size();
    for(int cell=0;cell<p_.layout.cells;++cell) {
        assemble_cell(cell);
        lu_.compute(a_);
        lu_.solve_multiple(local_rhs_,solved_);
        Eigen::Map<Mat>(x_.data()+static_cast<Eigen::Index>(cell)*m*t,m,t)=solved_.leftCols(t);
        y_.col(cell)=solved_.col(t);
        // C = weighted pointwise physical Jacobian after tensor face restriction.
        // This costs O(n^2 * t), rather than dense O(m*t^2).
        faces_jacobian_.cell_c_many(cell,solved_,face_action_);
        Eigen::Map<Mat>(schur_.data()+static_cast<Eigen::Index>(cell)*t*t,t,t)=-face_action_.leftCols(t);
        face_result_=face_action_.col(t);
        scatter(rhs_,cell,face_result_);
    }
}
void HybridOperator::apply(const Vec &direction,Vec &result) const {
    result.array()=faces_jacobian_.trace_diagonal().array()*direction.array();
    const int t=p_.local_trace_size();
    for(int cell=0;cell<p_.layout.cells;++cell) {
        gather(direction,cell);
        const Eigen::Map<const Mat> s(schur_.data()+static_cast<Eigen::Index>(cell)*t*t,t,t);
        face_result_.noalias()=s*faces_;
        scatter(result,cell,face_result_);
    }
}
void HybridOperator::block(int row,int col,BlockRef result) const {
    const int fw=p_.face_width(),t=p_.local_trace_size();
    result.setZero();
    if(row==col) result.diagonal()=faces_jacobian_.trace_diagonal().segment(row*fw,fw);
    const int axis=row/p_.layout.cells, first=row%p_.layout.cells;
    for(int adjacent=0;adjacent<2;++adjacent) {
        const int cell=adjacent==0?first:p_.neighbor(first,2*axis+1);
        const int row_side=2*axis+(adjacent==0?1:0);
        const auto fs=p_.faces(cell);
        const Eigen::Map<const Mat> s(schur_.data()+static_cast<Eigen::Index>(cell)*t*t,t,t);
        for(int side=0;side<NFACES;++side)
            if(fs[side]==col) result+=s.block(row_side*fw,side*fw,fw,fw);
    }
}
void HybridOperator::recover(const Vec &direction,Vec &result) const {
    const int m=p_.layout.block_size(),t=p_.local_trace_size();
    for(int cell=0;cell<p_.layout.cells;++cell) {
        gather(direction,cell);
        const Eigen::Map<const Mat> x(x_.data()+static_cast<Eigen::Index>(cell)*m*t,m,t);
        result.segment(cell*m,m).noalias()=-x*faces_;
        result.segment(cell*m,m)-=y_.col(cell);
    }
}
void HybridOperator::recovered_residual_from_schur(const Vec &error,Vec &result) const {
    result.setZero();
    const int n=p_.basis.nodes,fw=p_.face_width();
    for(int cell=0;cell<p_.layout.cells;++cell) {
        const auto fs=p_.faces(cell);
        for(int side=0;side<NFACES;++side) {
            const int face=fs[side];
            for(int v=0;v<NVAR;++v) for(int q=0;q<n;++q) {
                const double d=faces_jacobian_.trace_diagonal()[face*fw+v*n+q];
                const double tau=faces_jacobian_.tau(face,q,side%2!=0);
                face_flux_(q,v)=p_.face_weight(face,q)*tau*error[face*fw+v*n+q]/d;
            }
            p_.lift(side,face_flux_,p_.layout.cell(result,cell));
        }
    }
}
}
