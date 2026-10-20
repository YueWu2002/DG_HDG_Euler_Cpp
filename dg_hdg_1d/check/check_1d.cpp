#include "check_1d.hpp"
#include "dg_reference.hpp"
#include "euler/time_integrator.hpp"
#include <iostream>
#include <vector>
namespace euler::check {
namespace {
void near(double error, double tolerance, const char *message) {
    if (!std::isfinite(error) || error > tolerance) {
        std::cerr << message << ": error=" << error << ", tolerance=" << tolerance << '\n';
        throw std::runtime_error(message);
    }
}
Vec deterministic(int n) {
    Vec x(n);
    for (int i = 0; i < n; ++i)
        x[i] = std::sin(.7 * i) + .2 * std::cos(1.9 * i);
    return x;
}
void tensor_and_basis() {
    Layout layout{5, 4};
    auto shape = layout.shape();
    require(shape == std::array<Eigen::Index, 5>{layout.nodes, 1, NVAR, layout.cells, 1},
            "Column-major tensor shape mismatch");
    Vec storage(layout.size());
    for (int c = 0; c < layout.cells; ++c)
        for (int v = 0; v < NVAR; ++v)
            for (int i = 0; i < layout.nodes; ++i)
                storage[layout.index(c, v, i)] = 100 * c + 10 * v + i;
    for (int c = 0; c < layout.cells; ++c)
        for (int v = 0; v < NVAR; ++v)
            for (int i = 0; i < layout.nodes; ++i)
                near(std::abs(layout.cell(storage, c)(i, v) - (100 * c + 10 * v + i)), 0,
                     "Tensor stride mismatch");
    auto tensor = layout.tensor(storage);
    for (int c = 0; c < layout.cells; ++c)
        for (int v = 0; v < NVAR; ++v)
            for (int i = 0; i < layout.nodes; ++i)
                near(std::abs(tensor(i, 0, v, c, 0) - storage[layout.index(c, v, i)]), 0,
                     "DOF tensor view mismatch");
    BlockMatrix matrix(BlockGraph::periodic(layout.cells), layout.block_size());
    for (int slot = 1; slot < static_cast<int>(matrix.graph.col_idx.size()); ++slot)
        require(matrix.block(slot).data() - matrix.block(slot - 1).data() ==
                    layout.block_size() * layout.block_size(),
                "Nonzero blocks must share one contiguous numeric pool");
    require(matrix.graph.shape() == std::array<Eigen::Index, 4>{layout.cells, 1, layout.cells, 1},
            "Sparse outer tensor shape mismatch");
    auto block = matrix.tensor_block(0, 0, 1, 0);
    for (int v = 0; v < NVAR; ++v)
        for (int w = 0; w < NVAR; ++w)
            for (int i = 0; i < layout.nodes; ++i)
                for (int j = 0; j < layout.nodes; ++j)
                    block(i, 0, v, j, 0, w) = 1000 * v + 100 * w + 10 * i + j;
    const BlockMatrix &constant = matrix;
    auto constant_block = constant.tensor_block(0, 0, 1, 0);
    for (int v = 0; v < NVAR; ++v)
        for (int w = 0; w < NVAR; ++w)
            for (int i = 0; i < layout.nodes; ++i)
                for (int j = 0; j < layout.nodes; ++j) {
                    double expected = 1000 * v + 100 * w + 10 * i + j;
                    near(std::abs(matrix.at(0, 1)(layout.local(v, i), layout.local(w, j)) -
                                  expected),
                         0, "Sparse block matrix view mismatch");
                    near(std::abs(constant_block(i, 0, v, j, 0, w) - expected), 0,
                         "Sparse const tensor view mismatch");
                }
    for (int degree : {0, 1, 3, 6}) {
        NodalBasis basis(degree, degree + 3);
        near(std::abs(basis.w.sum() - 2), 1e-13, "Gauss constant integration");
        for (int p = 0; p <= 2 * degree + 1; ++p) {
            double integral = (basis.w.array() * basis.x.array().pow(p)).sum();
            near(std::abs(integral - (p % 2 == 0 ? 2.0 / (p + 1) : 0)), 1e-13,
                 "Gauss polynomial integration");
        }
        for (int p = 0; p <= degree; ++p) {
            Vec nodal = basis.x.array().pow(p);
            Vec derivative = basis.weighted_derivatives * nodal;
            Vec exact =
                basis.qw.array() *
                (p == 0 ? Vec::Zero(basis.points) : Vec(p * basis.qx.array().pow(p - 1))).array();
            near((derivative - exact).norm(), 2e-12, "Lagrange derivative polynomial reproduction");
        }
    }
    std::cout << "tensor layout and reference element passed\n";
}
void dense_lu() {
    // Cyclic permutation forces multiple row pivots; compare with an independent
    // Eigen full-pivot solve, including in-place RHS and right solves.
    int n = NVAR * 3;
    Mat a = Mat::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        a(i, (i + 1) % n) = 2 + i;
        a(i, i) = .01;
    }
    Vec b = deterministic(n), x(n);
    DenseLU lu(n);
    lu.compute(a);
    Eigen::FullPivLU<Mat> reference(a);
    Vec exact = reference.solve(b);
    lu.solve(b, x);
    near((x - exact).norm(), 1e-13, "Pivoted forward/back substitution");
    lu.solve(b, b);
    near((b - exact).norm(), 1e-13, "Aliased LU solve");
    Mat right = Mat::Random(NVAR, n), original = right;
    lu.right_solve(right);
    near((right * a - original).norm(), 1e-12, "Transposed right LU solve");
    std::cout << "pivoted local LU passed\n";
}
void linear_solvers() {
    // Dense block graph: ILU(0) has no missing fill, hence must be exact even
    // after block COLAMD. Path graph tests a truly sparse exact block ILU.
    for (bool full : {false, true})
        for (bool colamd : {false, true}) {
            constexpr int count = 7, width = NVAR;
            std::vector<std::vector<int>> rows(count);
            for (int i = 0; i < count; ++i)
                for (int j = 0; j < count; ++j)
                    if (full || std::abs(i - j) <= 1)
                        rows[i].push_back(j);
            BlockGraph graph(rows);
            BlockMatrix matrix(graph, width);
            for (int i = 0; i < count; ++i)
                for (int p = graph.row_ptr[i]; p < graph.row_ptr[i + 1]; ++p) {
                    int j = graph.col_idx[p];
                    matrix.block(p) = .03 * Mat::Random(width, width);
                    if (i == j) {
                        matrix.block(p).setZero();
                        for (int v = 0; v < NVAR; ++v)
                            matrix.block(p)(v, (v + 1) % NVAR) = 3 + v;
                    }
                }
            Vec exact = deterministic(count * width), rhs(exact.size()), x(exact.size()),
                weights = Vec::Ones(exact.size());
            matrix.apply(exact, rhs);
            BlockPreconditioner preconditioner(graph, width, Preconditioner::BlockILU0, colamd);
            require(preconditioner.factors().graph.col_idx.size() == graph.col_idx.size(),
                    "ILU(0) must preserve the permuted block count");
            preconditioner.refresh([&](int r, int c, BlockRef b) { b = matrix.at(r, c); });
            preconditioner.apply(rhs, x);
            // Natural path ordering and full graphs need no fill; COLAMD on a path
            // may introduce fill, in which case GMRES still has to recover the truth.
            if (full || !colamd)
                near((x - exact).norm(), 1e-12, "Exact block ILU property");
            GMRES gmres(x.size(), 12, weights);
            auto result =
                gmres.solve(Action::bind(matrix), Action::bind(preconditioner), rhs, x, 100, 1e-12);
            near((x - exact).norm(), 1e-10, "Right-preconditioned GMRES solution");
            near(result.relative_residual, 1e-12, "True GMRES residual");
            // Reuse after a nonzero solve: callers must see the new equation's
            // true residual even when the zero RHS takes the early return.
            rhs.setZero();
            auto zero = gmres.solve(Action::bind(matrix), Action::bind(preconditioner),
                                    rhs, x, 100, 1e-12);
            require(zero.iterations == 0 && zero.relative_residual == 0 &&
                    x.norm() == 0 && gmres.residual().norm() == 0,
                    "Zero-RHS GMRES retained a stale residual");
            matrix.apply(exact, rhs);
            DirectSolver direct(graph, width);
            direct.refresh([&](int r, int c, BlockRef b) { b = matrix.at(r, c); });
            direct.solve(rhs, x, weights);
            near((x - exact).norm(), 1e-12, "Direct block reference");
        }
    std::cout << "block ILU(0), block COLAMD, GMRES, direct reference passed\n";
}
void jacobians_and_condensation() {
    for (int degree : {0, 1, 3})
        for (int extra : {0, 2}) {
            Config config;
            config.cells = 5;
            config.degree = degree;
            config.quadrature = degree + 1 + extra;
            Problem p(config);
            int size = p.layout.size();
            Vec u(size), r(size), rp(size), rm(size), direction = deterministic(size), action(size),
                                                      assembled_action(size);
            p.initial(u);
            p.freeze(u);
            p.step_size = .003;
            // Test away from U^n as well as at it; frozen speeds must not be differentiated.
            u += .0001 * direction;
            DGOperator op(p);
            VolumeJacobian volume_jacobian(p);
            op.refresh(u);
            volume_jacobian.refresh(u);
            op.apply(direction, action);
            reference::dg_action(p, u, direction, assembled_action);
            near((action - assembled_action).norm() / action.norm(), 1e-13,
                 "Cached DG action versus general quadrature reference");
            Mat optimized(p.layout.block_size(), p.layout.block_size());
            Mat expected(p.layout.block_size(), p.layout.block_size());
            for (int cell = 0; cell < config.cells; ++cell) {
                volume_jacobian.block(cell, optimized);
                reference::volume_block(p, u, cell, expected);
                near((optimized - expected).norm() / expected.norm(), 1e-13,
                     "Cached volume block versus general quadrature reference");
            }
            double epsilon = 2e-7;
            Vec plus = u + epsilon * direction, minus = u - epsilon * direction;
            p.dg_residual(plus, rp);
            p.dg_residual(minus, rm);
            near((action - (rp - rm) / (2 * epsilon)).norm() / action.norm(), 3e-9,
                 "DG analytic Jacobian finite difference");
            BlockMatrix blocks(BlockGraph::periodic(config.cells), p.layout.block_size());
            for (int i = 0; i < config.cells; ++i)
                for (int at = blocks.graph.row_ptr[i]; at < blocks.graph.row_ptr[i + 1]; ++at)
                    op.block(i, blocks.graph.col_idx[at], blocks.block(at));
            blocks.apply(direction, assembled_action);
            near((action - assembled_action).norm() / action.norm(), 1e-13,
                 "DG block assembly versus matrix-free action");
            Vec lambda(config.cells * NVAR);
            p.trace(u, lambda);
            p.dg_residual(u, r);
            HybridOperator hybrid(p);
            hybrid.refresh(u, lambda);
            DirectSolver dg(blocks.graph, p.layout.block_size()), hdg(blocks.graph, NVAR);
            dg.refresh([&](int row, int col, BlockRef b) { op.block(row, col, b); });
            hdg.refresh([&](int row, int col, BlockRef b) { hybrid.block(row, col, b); });
            Vec du(size), dl(lambda.size()), recovered(size), weights = Vec::Ones(lambda.size()),
                                                              rhs = -r;
            dg.solve(rhs, du, p.inverse_mass);
            hdg.solve(hybrid.rhs(), dl, weights);
            hybrid.recover(dl, recovered);
            near((du - recovered).norm(), 2e-12, "DG versus condensed HDG Newton correction");
            BlockMatrix schur(blocks.graph, NVAR);
            for (int row = 0; row < config.cells; ++row)
                for (int at = schur.graph.row_ptr[row]; at < schur.graph.row_ptr[row + 1]; ++at)
                    hybrid.block(row, schur.graph.col_idx[at], schur.block(at));
            Vec test = deterministic(lambda.size()), sx(lambda.size()), sy(lambda.size());
            hybrid.apply(test, sx);
            schur.apply(test, sy);
            near((sx - sy).norm() / sx.norm(), 1e-13, "Schur blocks versus matrix-free action");
            // Nonprojected trace: compare to a full augmented finite-difference
            // Jacobian. This catches signs in B,C,D and the nonzero skeleton RHS.
            lambda += .001 * test;
            hybrid.refresh(u, lambda);
            hdg.refresh([&](int row, int col, BlockRef b) { hybrid.block(row, col, b); });
            hdg.solve(hybrid.rhs(), dl, weights);
            hybrid.recover(dl, recovered);
            int augmented_size = size + lambda.size();
            Mat jacobian(augmented_size, augmented_size);
            Vec volume(size), skeleton(lambda.size()), plus_volume(size),
                plus_skeleton(lambda.size());
            Vec minus_volume(size), minus_skeleton(lambda.size()), augmented_rhs(augmented_size),
                correction(augmented_size);
            hybrid.residual(u, lambda, volume, skeleton);
            augmented_rhs.head(size) = -volume;
            augmented_rhs.tail(lambda.size()) = -skeleton;
            for (int j = 0; j < augmented_size; ++j) {
                Vec up = u, um = u, lp = lambda, lm = lambda;
                if (j < size) {
                    up[j] += epsilon;
                    um[j] -= epsilon;
                } else {
                    lp[j - size] += epsilon;
                    lm[j - size] -= epsilon;
                }
                hybrid.residual(up, lp, plus_volume, plus_skeleton);
                hybrid.residual(um, lm, minus_volume, minus_skeleton);
                jacobian.col(j).head(size) = (plus_volume - minus_volume) / (2 * epsilon);
                jacobian.col(j).tail(lambda.size()) =
                    (plus_skeleton - minus_skeleton) / (2 * epsilon);
            }
            correction.head(size) = recovered;
            correction.tail(lambda.size()) = dl;
            near((jacobian * correction - augmented_rhs).norm() / augmented_rhs.norm(), 1e-7,
                 "Augmented HDG finite-difference Newton equation");
        }
    std::cout << "DG derivatives, HDG condensation, augmented trace equations passed\n";
}
void time_runs() {
    for (int degree : {0, 1, 3})
        for (int count : {3, 8, 12}) {
            // A three-cell high-degree projection is not resolved enough for this
            // positive smooth-wave evolution. Exercise three-cell connectivity with
            // P0; use resolved meshes for higher-degree time integration.
            if (degree > 0 && count == 3)
                continue;
            Config c;
            c.cells = count;
            c.degree = degree;
            // Cover both collocation and overintegration in complete time loops.
            c.quadrature = count == 8 ? degree + 1 : degree + 3;
            c.final_time = .015;
            c.nonlinear_tolerance = 1e-11;
            c.linear_tolerance = 1e-12;
            Vec reference;
            for (LinearSolver linear : {LinearSolver::Direct, LinearSolver::GMRES})
                for (Preconditioner precond :
                     {Preconditioner::BlockILU0, Preconditioner::BlockJacobi})
                    for (bool colamd : {false, true})
                        for (Method method : {Method::DG, Method::HDG}) {
                            if (linear == LinearSolver::Direct &&
                                (colamd || precond == Preconditioner::BlockJacobi))
                                continue;
                            c.linear = linear;
                            c.preconditioner = precond;
                            c.colamd = colamd;
                            std::cout
                                << "time case: cells=" << count << " k=" << degree << " linear="
                                << (linear == LinearSolver::Direct ? "direct" : "gmres")
                                << " precond="
                                << (precond == Preconditioner::BlockILU0 ? "ilu0" : "jacobi")
                                << " colamd=" << colamd
                                << " method=" << (method == Method::DG ? "dg" : "hdg") << '\n';
                            Problem p(c);
                            Vec solution(p.layout.size());
                            p.initial(solution);
                            TimeIntegrator integrator(p, method);
                            // Guard is active over the FIRST step and every later step,
                            // including LU numeric updates, Arnoldi, and line searches.
                            auto result = integrator.advance(solution, true);
                            near(std::abs(result.final_time - c.final_time), 0,
                                 "Final-time clipping");
                            near(result.conservation_drift, 2e-11, "Periodic conservation");
                            require(p.admissible(solution), "Admissibility after Newton");
                            if (reference.size() == 0)
                                reference = solution;
                            else
                                near((solution - reference).cwiseAbs().maxCoeff(), 2e-10,
                                     "DG/HDG/direct/GMRES final states");
                        }
        }
    // Independent temporal accuracy check against the analytic smooth wave.
    double previous = 0;
    for (double cfl : {.4, .2, .1}) {
        Config c;
        c.integrator = "be";
        c.cells = 24;
        c.degree = 3;
        c.quadrature = 6;
        c.cfl = cfl;
        c.final_time = .04;
        Problem p(c);
        Vec solution(p.layout.size());
        p.initial(solution);
        TimeIntegrator integrator(p, Method::DG);
        integrator.advance(solution);
        double error = p.error(solution, c.final_time)[Density];
        if (previous > 0)
            require(error < .8 * previous,
                    "Backward Euler temporal refinement did not reduce error");
        previous = error;
    }
    // Large direct matrix exercises the factorization workspace, not only
    // tiny local LU. This also tests a first step with allocation prohibition.
    Config c;
    c.cells = 48;
    c.degree = 3;
    c.steps = 2;
    c.linear = LinearSolver::Direct;
    Problem p(c);
    Vec u(p.layout.size());
    p.initial(u);
    TimeIntegrator direct(p, Method::DG);
    direct.advance(u);
    std::cout << "allocation-free time loops, conservation, solver agreement, temporal refinement "
                 "passed\n";
}
} // namespace
void all() {
    tensor_and_basis();
    dense_lu();
    linear_solvers();
    jacobians_and_condensation();
    linearization_lifecycle();
    time_runs();
}
} // namespace euler::check
