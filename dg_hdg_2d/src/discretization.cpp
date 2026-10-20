#include "euler/discretization.hpp"
#include <algorithm>
namespace euler {
Problem::Problem(const Config &c)
    : config(c), basis(c.degree, c.quadrature), layout(c.nx, c.ny, basis.nodes), example(c.testcase, c.alpha_degrees),
      hx(example.lx / c.nx), hy(example.ly / c.ny), old(layout.size()), mass(layout.size()),
      inverse_mass(layout.size()), trace_mass(trace_size()), inverse_trace_mass(trace_size()),
      speeds(face_count() * basis.nodes), tensor_(basis, hx, hy),
      points_(basis.points * basis.points, NVAR), flux_x_(points_.rows(), NVAR),
      flux_y_(points_.rows(), NVAR), left_(basis.nodes, NVAR), right_(basis.nodes, NVAR),
      face_flux_(basis.nodes, NVAR), other_flux_(basis.nodes, NVAR) {
    validate(c);
    for (int cell = 0; cell < layout.cells; ++cell)
        for (int v = 0; v < NVAR; ++v)
            for (int j = 0; j < basis.nodes; ++j)
                for (int i = 0; i < basis.nodes; ++i)
                    mass[layout.index(cell, v, i + basis.nodes * j)] =
                        .25 * hx * hy * basis.w[i] * basis.w[j];
    for (int f = 0; f < face_count(); ++f)
        for (int v = 0; v < NVAR; ++v)
            for (int q = 0; q < basis.nodes; ++q)
                trace_mass[f * face_width() + v * basis.nodes + q] = face_weight(f, q);
    inverse_mass = mass.cwiseInverse();
    inverse_trace_mass = trace_mass.cwiseInverse();
}
int Problem::neighbor(int cell, int side) const {
    int x = cell % layout.nx, y = cell / layout.nx;
    if (side < 2) x = (x + (side == 0 ? layout.nx - 1 : 1)) % layout.nx;
    else y = (y + (side == 2 ? layout.ny - 1 : 1)) % layout.ny;
    return x + layout.nx * y;
}
std::array<int, NFACES> Problem::faces(int cell) const {
    return {neighbor(cell, 0), cell, layout.cells + neighbor(cell, 2), layout.cells + cell};
}
double Problem::face_weight(int face, int q) const {
    return .5 * (face < layout.cells ? hy : hx) * basis.w[q];
}
int Problem::face_node(int side, int q, int i) const {
    return side < 2 ? i + basis.nodes * q : q + basis.nodes * i;
}
State Problem::sample(const Vec &u, int cell, int side, int q) const {
    State state = State::Zero();
    const Vec &shape = side % 2 ? basis.right : basis.left;
    for (int v = 0; v < NVAR; ++v)
        for (int i = 0; i < basis.nodes; ++i)
            state[v] += shape[i] * u[layout.index(cell, v, face_node(side, q, i))];
    return state;
}
namespace {
template<int N> void restrict_face(ConstBlockRef nodal,int side,const Vec &shape,BlockRef values) {
    using Sheet=Eigen::Matrix<double,N,N>;
    using Line=Eigen::Matrix<double,N,1>;
    Eigen::Map<const Line> e(shape.data());
    for(int v=0;v<NVAR;++v) {
        Eigen::Map<const Sheet> sheet(nodal.col(v).data());
        Eigen::Map<Line> out(values.col(v).data());
        if(side<2) out.noalias()=sheet.transpose()*e;
        else out.noalias()=sheet*e;
    }
}
template<int N> void lift_face(int side,ConstBlockRef values,BlockRef nodal,const Vec &shape,double scale) {
    using Sheet=Eigen::Matrix<double,N,N>;
    using Line=Eigen::Matrix<double,N,1>;
    Eigen::Map<const Line> e(shape.data());
    for(int v=0;v<NVAR;++v) {
        Eigen::Map<Sheet> sheet(nodal.col(v).data());
        Eigen::Map<const Line> f(values.col(v).data());
        if(side<2) sheet.noalias()+=scale*e*f.transpose();
        else sheet.noalias()+=scale*f*e.transpose();
    }
}
}
void Problem::face_values(ConstBlockRef nodal,int side,BlockRef values) const {
    const Vec &shape=side%2?basis.right:basis.left;
    switch(basis.nodes) {
    case 1: restrict_face<1>(nodal,side,shape,values); return;
    case 2: restrict_face<2>(nodal,side,shape,values); return;
    case 3: restrict_face<3>(nodal,side,shape,values); return;
    case 4: restrict_face<4>(nodal,side,shape,values); return;
    case 5: restrict_face<5>(nodal,side,shape,values); return;
    case 6: restrict_face<6>(nodal,side,shape,values); return;
    case 7: restrict_face<7>(nodal,side,shape,values); return;
    case 8: restrict_face<8>(nodal,side,shape,values); return;
    case 9: restrict_face<9>(nodal,side,shape,values); return;
    }
    throw std::runtime_error("Unsupported face degree");
}
void Problem::lift(int side,ConstBlockRef values,BlockRef nodal,double scale) const {
    const Vec &shape=side%2?basis.right:basis.left;
    switch(basis.nodes) {
    case 1: lift_face<1>(side,values,nodal,shape,scale); return;
    case 2: lift_face<2>(side,values,nodal,shape,scale); return;
    case 3: lift_face<3>(side,values,nodal,shape,scale); return;
    case 4: lift_face<4>(side,values,nodal,shape,scale); return;
    case 5: lift_face<5>(side,values,nodal,shape,scale); return;
    case 6: lift_face<6>(side,values,nodal,shape,scale); return;
    case 7: lift_face<7>(side,values,nodal,shape,scale); return;
    case 8: lift_face<8>(side,values,nodal,shape,scale); return;
    case 9: lift_face<9>(side,values,nodal,shape,scale); return;
    }
    throw std::runtime_error("Unsupported face degree");
}
void Problem::evaluate_volume_flux(ConstBlockRef states) const {
    for (int q = 0; q < states.rows(); ++q) {
        State u = states.row(q).transpose();
        auto primitive = example.gas.primitive(u);
        flux_x_.row(q) = example.gas.flux(u, primitive, 0).transpose();
        flux_y_.row(q) = example.gas.flux(u, primitive, 1).transpose();
    }
}
void Problem::volume_residual(const Vec &u, Vec &r) const {
    r.array() = mass.array() * (u - old).array() / step_size;
    if (basis.nodes == 1) return;
    for (int cell = 0; cell < layout.cells; ++cell) {
        if (basis.points == basis.nodes) evaluate_volume_flux(layout.cell(u, cell));
        else {
            tensor_.interpolate(layout.cell(u, cell), points_);
            evaluate_volume_flux(points_);
        }
        tensor_.subtract_weak(flux_x_, flux_y_, layout.cell(r, cell));
    }
}
void Problem::dg_residual(const Vec &u, Vec &r) const {
    volume_residual(u, r);
    for (int axis = 0; axis < NDIM; ++axis)
        for (int cell = 0; cell < layout.cells; ++cell) {
            int side = 2 * axis + 1, next = neighbor(cell, side), face = axis * layout.cells + cell;
            face_values(layout.cell(u, cell), side, left_);
            face_values(layout.cell(u, next), side - 1, right_);
            for (int q = 0; q < basis.nodes; ++q) {
                State l = left_.row(q).transpose(), rr = right_.row(q).transpose();
                auto s = speeds[face * basis.nodes + q];
                face_flux_.row(q) = (face_weight(face, q) / (s.a + s.b)) *
                    (s.b * example.gas.flux(l, axis) + s.a * example.gas.flux(rr, axis) +
                     s.a * s.b * (l - rr)).transpose();
            }
            lift(side, face_flux_, layout.cell(r, cell));
            lift(side - 1, face_flux_, layout.cell(r, next), -1);
        }
}
void Problem::hybrid_residual(const Vec &u, const Vec &lambda, Vec &r, Vec &g) const {
    volume_residual(u, r);
    for (int axis = 0; axis < NDIM; ++axis)
        for (int cell = 0; cell < layout.cells; ++cell) {
            int side = 2 * axis + 1, next = neighbor(cell, side), face = axis * layout.cells + cell;
            face_values(layout.cell(u, cell), side, left_);
            face_values(layout.cell(u, next), side - 1, right_);
            Eigen::Map<const Mat> trace(lambda.data() + face * face_width(), basis.nodes, NVAR);
            Eigen::Map<Mat> skeleton(g.data() + face * face_width(), basis.nodes, NVAR);
            for (int q = 0; q < basis.nodes; ++q) {
                State l = left_.row(q).transpose(), rr = right_.row(q).transpose();
                State z = trace.row(q).transpose();
                auto s = speeds[face * basis.nodes + q];
                double weight = face_weight(face, q);
                face_flux_.row(q) = weight * (example.gas.flux(l, axis) + s.a * (l - z)).transpose();
                other_flux_.row(q) = weight * (example.gas.flux(rr, axis) + s.b * (z - rr)).transpose();
            }
            skeleton = face_flux_ - other_flux_;
            lift(side, face_flux_, layout.cell(r, cell));
            lift(side - 1, other_flux_, layout.cell(r, next), -1);
        }
}
double Problem::freeze(const Vec &u) {
    old = u;
    for (int axis = 0; axis < NDIM; ++axis)
        for (int cell = 0; cell < layout.cells; ++cell) {
            int side = 2 * axis + 1, next = neighbor(cell, side), face = axis * layout.cells + cell;
            face_values(layout.cell(u, cell), side, left_);
            face_values(layout.cell(u, next), side - 1, right_);
            for (int q = 0; q < basis.nodes; ++q)
                speeds[face * basis.nodes + q] = einfeldt(left_.row(q).transpose(), right_.row(q).transpose(),
                                                         example.gas, config.safety, axis);
        }
    double rate = 0;
    for (int cell = 0; cell < layout.cells; ++cell) {
        std::array<double, NDIM> alpha{};
        auto fs = faces(cell);
        for (int side = 0; side < NFACES; ++side)
            for (int q = 0; q < basis.nodes; ++q) {
                auto s = speeds[fs[side] * basis.nodes + q];
                alpha[side / 2] = std::max({alpha[side / 2], s.a, s.b});
            }
        tensor_.interpolate(layout.cell(u, cell), points_);
        for (int q = 0; q < points_.rows(); ++q) {
            State state = points_.row(q).transpose();
            auto primitive = example.gas.primitive(state);
            require(state[Density] > 0 && primitive.pressure > 0, "Invalid CFL volume state");
            double sound = std::sqrt(example.gas.gamma() * primitive.pressure * primitive.inverse_density);
            for (int axis = 0; axis < NDIM; ++axis)
                alpha[axis] = std::max(alpha[axis], std::abs(primitive.velocity[axis]) + sound);
        }
        rate = std::max(rate, alpha[0] / hx + alpha[1] / hy);
    }
    return rate;
}
void Problem::trace(const Vec &u, Vec &lambda) const {
    for (int axis = 0; axis < NDIM; ++axis)
        for (int cell = 0; cell < layout.cells; ++cell) {
            int side = 2 * axis + 1, next = neighbor(cell, side), face = axis * layout.cells + cell;
            face_values(layout.cell(u, cell), side, left_);
            face_values(layout.cell(u, next), side - 1, right_);
            Eigen::Map<Mat> trace(lambda.data() + face * face_width(), basis.nodes, NVAR);
            for (int q = 0; q < basis.nodes; ++q) {
                State l = left_.row(q).transpose(), r = right_.row(q).transpose();
                auto s = speeds[face * basis.nodes + q];
                trace.row(q) = ((s.a * l + s.b * r + example.gas.flux(l, axis) -
                                example.gas.flux(r, axis)) / (s.a + s.b)).transpose();
            }
        }
}
bool Problem::physical(ConstBlockRef states) const {
    if (!states.allFinite()) return false;
    for (int q = 0; q < states.rows(); ++q) {
        State state = states.row(q).transpose();
        if (state[Density] <= 1e-13 || !(example.gas.pressure(state) > 1e-13)) return false;
    }
    return true;
}
bool Problem::admissible(const Vec &u) const {
    if (u.size() != layout.size()) return false;
    for (int cell = 0; cell < layout.cells; ++cell) {
        if (!physical(layout.cell(u, cell))) return false;
        if (basis.points != basis.nodes) {
            tensor_.interpolate(layout.cell(u, cell), points_);
            if (!physical(points_)) return false;
        }
        for (int side = 0; side < NFACES; ++side) {
            face_values(layout.cell(u, cell), side, left_);
            if (!physical(left_)) return false;
        }
    }
    return true;
}
double Problem::residual_norm(const Vec &r) const {
    return std::sqrt((r.array().square() * inverse_mass.array()).sum() / (example.lx * example.ly));
}
State Problem::integral(const Vec &u) const {
    State result = State::Zero();
    for (int cell = 0; cell < layout.cells; ++cell)
        for (int v = 0; v < NVAR; ++v)
            result[v] += layout.cell(u, cell).col(v).dot(layout.cell(mass, cell).col(v));
    return result;
}
BlockGraph Problem::trace_graph() const {
    std::vector<std::vector<int>> adjacency(face_count());
    for (int cell = 0; cell < layout.cells; ++cell) {
        auto fs = faces(cell);
        for (int row : fs)
            for (int col : fs) adjacency[row].push_back(col);
    }
    return BlockGraph(adjacency);
}
} // namespace euler
