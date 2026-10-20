#include "euler/time_integrator.hpp"
#include <iostream>
#include <algorithm>
namespace {
using namespace euler;
void near(double e,double tolerance,const char *message) {
    if(!std::isfinite(e)||e>tolerance) { std::cerr<<message<<": "<<e<<" > "<<tolerance<<'\n'; throw std::runtime_error(message); }
}
Vec direction(int size) {
    Vec out(size); for(int i=0;i<size;++i) out[i]=.2*std::sin(.31*i)+.1*std::cos(.79*i); return out;
}
void jacobi_pool() {
    for(int width:{4,12,36,324}) {
        constexpr int blocks=3;
        Mat a=Mat::Zero(width,width);
        for(int i=0;i<width;++i) { a(i,(i+1)%width)=2.+.01*i; a(i,i)+=.1; }
        Vec exact=direction(blocks*width),rhs(exact.size()),solution(exact.size());
        for(int i=0;i<blocks;++i) rhs.segment(i*width,width).noalias()=a*exact.segment(i*width,width);
        BlockJacobi pc(blocks,width);
        Eigen::internal::set_is_malloc_allowed(false);
        pc.refresh([&](int,int,BlockRef out){out=a;});pc.apply(rhs,solution);
        Eigen::internal::set_is_malloc_allowed(true);
        near((solution-exact).norm(),2e-12,"Contiguous Jacobi LU and forced permutation");
    }
}
void physics() {
    Gas gas(3.);
    for(double rho:{.3,.7,1.3}) for(int axis=0;axis<2;++axis) {
        State u;u<<rho,.7*rho,-.2*rho,1.7+.265*rho;
        const auto primitive=gas.primitive(u);
        const auto cached=gas.cached_jacobian(primitive.velocity[0],primitive.velocity[1],primitive.enthalpy,axis);
        near((cached-gas.jacobian(u,axis)).norm(),1e-13,"Cached physical derivative");
        const State d=direction(NVAR); const double eps=2e-7;
        near(((gas.flux(u+eps*d,axis)-gas.flux(u-eps*d,axis))/(2*eps)-cached*d).norm(),1e-8,"Physics finite difference");
    }
}
void tensor_reference(int degree,bool overintegrate) {
    NodalBasis b(degree,degree+1+(overintegrate?2:0));
    const int n=b.nodes,Q=b.points; const double hx=.3,hy=.7;
    TensorProduct t(b,hx,hy);
    Mat fx(Q*Q,NVAR),fy(Q*Q,NVAR),out=Mat::Zero(n*n,NVAR),reference=out;
    for(int v=0;v<NVAR;++v) for(int q=0;q<Q*Q;++q) {
        fx(q,v)=std::sin(.23*q+.4*v);fy(q,v)=std::cos(.12*q-.7*v);
    }
    t.subtract_weak(fx,fy,out);
    for(int v=0;v<NVAR;++v) for(int j=0;j<n;++j) for(int i=0;i<n;++i)
        for(int qy=0;qy<Q;++qy) for(int qx=0;qx<Q;++qx)
            reference(i+n*j,v)-=.5*hy*b.weighted_derivatives(qx,i)*b.qw[qy]*b.values(qy,j)*fx(qx+Q*qy,v)
                +.5*hx*b.qw[qx]*b.values(qx,i)*b.weighted_derivatives(qy,j)*fy(qx+Q*qy,v);
    near((out-reference).norm(),3e-12,"Sum factorization against independent quadrature");
}
void operators(int degree,bool overintegrate) {
    Config c;c.nx=3;c.ny=4;c.degree=degree;c.quadrature=degree+1+(overintegrate?2:0);
    Problem p(c);Vec u(p.layout.size());p.initial(u);p.freeze(u);p.step_size=.003;
    DGOperator op(p); HybridOperator hybrid(p);
    Vec d=direction(u.size()),jv(u.size()),r(u.size()),rp(u.size()),rm(u.size()),up(u.size()),um(u.size());
    Vec lambda(p.trace_size()),g(p.trace_size()),hr(u.size()),dl=direction(p.trace_size()),sdl(p.trace_size()),sblock(p.trace_size());
    Vec du(u.size()),schur_error(p.trace_size()),check(u.size());
    Mat block(p.layout.block_size(),p.layout.block_size());Vec assembled=Vec::Zero(u.size());
    // First refresh/apply/block and all condensation work must be allocation-free.
    Eigen::internal::set_is_malloc_allowed(false);
    op.refresh(u);op.apply(d,jv);p.dg_residual(u,r);
    const double eps=2e-7;
    up=u+eps*d;um=u-eps*d;p.dg_residual(up,rp);p.dg_residual(um,rm);
    near((jv-(rp-rm)/(2*eps)).norm()/std::max(1.,jv.norm()),2e-8,"DG derivative finite difference");
    const auto graph=BlockGraph::cartesian(c.nx,c.ny); // std::vector allocations only
    for(int row=0;row<p.layout.cells;++row) for(int slot=graph.row_ptr[row];slot<graph.row_ptr[row+1];++slot) {
        const int col=graph.col_idx[slot];op.block(row,col,block);
        assembled.segment(row*block.rows(),block.rows()).noalias()+=block*d.segment(col*block.cols(),block.cols());
    }
    near((assembled-jv).norm()/std::max(1.,jv.norm()),2e-12,"DG matrix-free/block identity");
    p.trace(u,lambda);hybrid.residual(u,lambda,hr,g);
    near((hr-r).norm(),1e-12,"Projected DG/HDG residual identity");near(g.norm(),1e-12,"Projected skeleton equation");
    hybrid.refresh_projected(u,r);hybrid.apply(dl,sdl);hybrid.recover(dl,du);
    schur_error=sdl-hybrid.rhs();hybrid.recovered_residual_from_schur(schur_error,check);
    op.apply(du,assembled);assembled+=r;
    near((check-assembled).norm()/std::max(1.,assembled.norm()),3e-11,"Recovered DG residual identity");
    Eigen::internal::set_is_malloc_allowed(true);
    const auto tg=p.trace_graph(); Mat sb(p.face_width(),p.face_width());sblock.setZero();
    for(int row=0;row<tg.rows;++row) for(int slot=tg.row_ptr[row];slot<tg.row_ptr[row+1];++slot) {
        const int col=tg.col_idx[slot];hybrid.block(row,col,sb);
        sblock.segment(row*sb.rows(),sb.rows()).noalias()+=sb*dl.segment(col*sb.cols(),sb.cols());
    }
    near((sblock-sdl).norm()/std::max(1.,sdl.norm()),2e-12,"Schur matrix-free/block identity");
    std::cout<<"degree="<<degree<<" overintegration="<<overintegrate<<" operators passed\n";
}
void stepping() {
    for(const std::string integrator:{"be","ssp-dirk3"}) {
        Vec reference;
        for(Method method:{Method::DG,Method::HDG}) {
            Config c;c.nx=4;c.ny=4;c.degree=2;c.final_time=.01;c.cfl=2.;c.integrator=integrator;c.check_update=true;
            Problem p(c);Vec u(p.layout.size());p.initial(u);TimeIntegrator solver(p,method);const auto result=solver.advance(u);
            near(result.conservation_drift,2e-11,"Conservation in production stepping");
            if(reference.size()==0) reference=u;else near((u-reference).cwiseAbs().maxCoeff(),2e-10,"DG/HDG final solution agreement");
        }
    }
}
}
int main() {
    try {
        std::cout.setf(std::ios::unitbuf);
        physics();jacobi_pool();
        for(int k:{0,1,2,4,8}) for(bool over:{false,true}) { tensor_reference(k,over);operators(k,over); }
        stepping();std::cout<<"All 2D checks passed\n";
    } catch(const std::exception &e) {
        Eigen::internal::set_is_malloc_allowed(true);
        std::cerr<<"CHECK FAILED: "<<e.what()<<'\n';return 1;
    }
}
