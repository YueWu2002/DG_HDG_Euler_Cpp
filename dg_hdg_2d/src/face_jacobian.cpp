#include "euler/face_jacobian.hpp"
namespace euler {
FaceJacobian::FaceJacobian(const Problem &p)
    : p_(p), primitives_left_(p.basis.nodes,3*p.face_count()),
      primitives_right_(primitives_left_.rows(),primitives_left_.cols()),
      tau_left_(p.basis.nodes,p.face_count()),tau_right_(tau_left_.rows(),tau_left_.cols()),
      weights_(tau_left_.rows(),tau_left_.cols()),blend_left_(weights_.rows(),weights_.cols()),
      blend_right_(weights_.rows(),weights_.cols()),penalty_(weights_.rows(),weights_.cols()),
      diagonal_(p.trace_size()),left_(p.basis.nodes,NVAR),right_(p.basis.nodes,NVAR),
      flux_left_(p.basis.nodes,NVAR),flux_right_(p.basis.nodes,NVAR),
      samples_(p.face_width(),p.local_trace_size()+1),coefficients_(p.basis.nodes,NVAR*NVAR),
      pressure_direction_(p.basis.nodes) {
    for(int face=0;face<p.face_count();++face) for(int q=0;q<p.basis.nodes;++q)
        weights_(q,face)=p.face_weight(face,q);
}
void FaceJacobian::refresh(const Vec &state) {
    for(int axis=0;axis<NDIM;++axis) for(int cell=0;cell<p_.layout.cells;++cell) {
        const int face=axis*p_.layout.cells+cell, side=2*axis+1;
        p_.face_values(p_.layout.cell(state,cell),side,left_);
        p_.face_values(p_.layout.cell(state,p_.neighbor(cell,side)),side-1,right_);
        for(int q=0;q<p_.basis.nodes;++q) {
            const auto l=p_.example.gas.primitive(left_.row(q).transpose());
            const auto r=p_.example.gas.primitive(right_.row(q).transpose());
            primitives_left_(q,3*face)=l.velocity[0]; primitives_left_(q,3*face+1)=l.velocity[1];
            primitives_left_(q,3*face+2)=l.enthalpy;
            primitives_right_(q,3*face)=r.velocity[0]; primitives_right_(q,3*face+1)=r.velocity[1];
            primitives_right_(q,3*face+2)=r.enthalpy;
            const auto s=p_.speeds[face*p_.basis.nodes+q];
            tau_left_(q,face)=s.a; tau_right_(q,face)=s.b;
            blend_left_(q,face)=weights_(q,face)*s.b/(s.a+s.b);
            blend_right_(q,face)=weights_(q,face)*s.a/(s.a+s.b);
            penalty_(q,face)=weights_(q,face)*s.a*s.b/(s.a+s.b);
            for(int v=0;v<NVAR;++v)
                diagonal_[face*p_.face_width()+v*p_.basis.nodes+q]=-weights_(q,face)*(s.a+s.b);
        }
    }
}
void FaceJacobian::normal_action(const Mat &primitive,int face,ConstBlockRef direction,BlockRef result,int axis) const {
    const auto u=primitive.col(3*face).array(), v=primitive.col(3*face+1).array();
    const auto H=primitive.col(3*face+2).array();
    const auto dr=direction.col(Density).array(), dx=direction.col(MomentumX).array();
    const auto dy=direction.col(MomentumY).array(), de=direction.col(Energy).array();
    pressure_direction_.array()=(p_.example.gas.gamma()-1.)*(de-u*dx-v*dy+.5*(u.square()+v.square())*dr);
    if(axis==0) {
        result.col(Density)=direction.col(MomentumX);
        result.col(MomentumX).array()=2*u*dx-u.square()*dr+pressure_direction_.array();
        result.col(MomentumY).array()=u*dy+v*dx-u*v*dr;
        result.col(Energy).array()=H*(dx-u*dr)+u*(de+pressure_direction_.array());
    } else {
        result.col(Density)=direction.col(MomentumY);
        result.col(MomentumX).array()=u*dy+v*dx-u*v*dr;
        result.col(MomentumY).array()=2*v*dy-v.square()*dr+pressure_direction_.array();
        result.col(Energy).array()=H*(dy-v*dr)+v*(de+pressure_direction_.array());
    }
}
void FaceJacobian::hll_action(int face,ConstBlockRef l,ConstBlockRef r,BlockRef flux) const {
    const int axis=face/p_.layout.cells;
    normal_action(primitives_left_,face,l,flux_left_,axis);
    normal_action(primitives_right_,face,r,flux_right_,axis);
    for(int v=0;v<NVAR;++v)
        flux.col(v).array()=blend_left_.col(face).array()*flux_left_.col(v).array()+
            blend_right_.col(face).array()*flux_right_.col(v).array()+
            penalty_.col(face).array()*(l.col(v)-r.col(v)).array();
}
FluxJacobian FaceJacobian::physical(int face,int q,bool left) const {
    const Mat &a=left?primitives_left_:primitives_right_;
    return p_.example.gas.cached_jacobian(a(q,3*face),a(q,3*face+1),a(q,3*face+2),face/p_.layout.cells);
}
FluxJacobian FaceJacobian::hll(int face,int q,bool left) const {
    FluxJacobian result=physical(face,q,left)*(left?blend_left_(q,face):blend_right_(q,face));
    result.diagonal().array()+=left?penalty_(q,face):-penalty_(q,face);
    return result;
}
double FaceJacobian::tau(int face,int q,bool left) const {
    return left?tau_left_(q,face):tau_right_(q,face);
}
FluxJacobian FaceJacobian::hdg(int face,int q,bool left) const {
    FluxJacobian result=(left?1.:-1.)*physical(face,q,left);
    result.diagonal().array()+=tau(face,q,left);
    result*=weights_(q,face);
    return result;
}
void FaceJacobian::cell_c(int cell,Eigen::Ref<const Vec> direction,Eigen::Ref<Vec> result) const {
    Eigen::Map<const Mat> input(direction.data(),p_.layout.block_size(),1);
    Eigen::Map<Mat> output(result.data(),p_.local_trace_size(),1);
    cell_c_many(cell,input,output);
}
void FaceJacobian::cell_c_many(int cell,ConstBlockRef directions,BlockRef result) const {
    const int n=p_.basis.nodes, fw=p_.face_width(), columns=directions.cols();
    require(columns<=samples_.cols(),"Too many local RHS columns");
    const auto fs=p_.faces(cell);
    for(int side=0;side<NFACES;++side) {
        for(int j=0;j<columns;++j) {
            Eigen::Map<const Mat> nodal(directions.col(j).data(),p_.layout.plane(),NVAR);
            Eigen::Map<Mat> values(samples_.col(j).data(),n,NVAR);
            p_.face_values(nodal,side,values);
        }
        for(int q=0;q<n;++q) {
            const auto derivative=hdg(fs[side],q,side%2!=0);
            for(int w=0;w<NVAR;++w) for(int v=0;v<NVAR;++v)
                coefficients_(q,w*NVAR+v)=derivative(v,w);
        }
        for(int v=0;v<NVAR;++v) {
            auto output=result.block(side*fw+v*n,0,n,columns);
            output.array()=(samples_.block(0,0,n,columns).array().colwise()*coefficients_.col(v).array());
            for(int w=1;w<NVAR;++w)
                output.array()+=(samples_.block(w*n,0,n,columns).array().colwise()*coefficients_.col(w*NVAR+v).array());
        }
    }
}
}
