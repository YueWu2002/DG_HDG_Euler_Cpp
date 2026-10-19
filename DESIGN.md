# 1D 实现设计

## 张量与块

目标全局 DOF 张量采用 Eigen 列主序 `[kx+1,ky+1,NVAR,Nx,Ny]`，第一维连续。
1D 固定 `Ny=1, ky=0`，因此逻辑形状为 `[k+1,1,NVAR,Nx,1]`。
`Layout::tensor` 提供真实的 Eigen TensorMap，`index/local/cell` 提供同一内存上的访问接口。
`cell` 是零拷贝的 Eigen `[nodes,NVAR]` 列主序矩阵视图。
变量数仅由 `inline constexpr int NVAR` 决定，分量使用 `Density/Momentum/Energy`。
物理公式里的系数 3、左右两个面等数学常数不属于变量个数。

DG 单块宽度 `m=NVAR*(k+1)`；单块矩阵的行列逻辑维度为
`[k+1,1,NVAR,k+1,1,NVAR]`。外层稀疏连接图坐标为
`[Nx_row,Ny_row,Nx_col,Ny_col]`，1D 为 `[Nx,1,Nx,1]`。
周期左右邻居构成块 CSR 图，保留整块结构零；不会把变量或节点独立交给 COLAMD。

BCSR 的 `row_ptr` 与 `col_idx` 只索引单元块；数值统一放在一个 Eigen 对齐的连续数值池。
`block(slot)` 映射池内一个列主序 `[m,m]` 矩阵，
`tensor_block(ix_row,iy_row,ix_col,iy_col)` 映射同一内存上的六维列主序张量。
块内偏移为 `i + n*v + m*(j+n*w)`，再加 `slot*m*m` 得到池内偏移。
不使用 `vector<MatrixXd>` 为各非零块分散分配内存；稀疏结构一经构造不可修改。

HDG 的 A 为 `[m,m]`，B 为 `[m,NFACES*NVAR]`，C 为 `[NFACES*NVAR,m]`；
局部面按左、右排序，面内按变量排序。全局 trace 列主序张量为 `[NVAR,Nx]`，
Schur 块宽度为 NVAR，仍按界面块图处理。

## 数学路径

空间离散沿用原项目：FastGL 生成 Gauss 节点/权重、Lagrange 基、对角质量、弱通量积分、
周期 Einfeldt/HLL 面通量。参考单元一次构造，默认体积分点数为 k+1，
可用 `--quad` 过积分；初始化使用高阶 L2 投影和守恒的凸正性缩放。

每步在 Uⁿ 冻结面波速，`rate=max_cell(alpha/h)`，alpha 同时覆盖体积分点声速与面波速。
`dt=CFL/rate`，不额外插入次数因子。终止时间模式截断最后一步；步数模式不截断。
自动选择的 dt 是算法必需的内部量，用户不能指定固定 dt。

DG 的 Newton 算子由 `DGOperator` 表示；它缓存体积分点和面的解析通量导数，
不存储全局 Jacobian。预条件器和直接参考通过 `DGOperator::block(row,col,out)`
获取需要的局部块，GMRES 通过 `DGOperator::apply(direction,result)` 计算 Jv。

HDG 保持独立 augmented trace，不在 Newton 更新后重投影：

```text
[ A  B ] [du] = -[r]
[ C  D ] [dl]   [g]
X = solve(A,B),  y = solve(A,r)
S dl = -g + C y,  S = D - C X
du = -y - X dl
```

A 按单元做选主元 LU。保存 `X=A⁻¹B` 的作用结果与 LU，不保存 A⁻¹。
Schur 的 `apply` 按单元 gather → X → C → scatter；预条件和直接路径可查询 Schur 的单块。
DG 用质量对偶 L2 RMS 残差；HDG 线搜索同时约束体积和 skeleton 残差。
两者都检查节点、体积分点和边界的密度及压力正性。

## 线性化模块与调用顺序

`discretization.hpp/.cpp` 中的 `Problem` 负责状态、残差、积分、正性等；
线性化单独放在以下模块，计算库不依赖 check 目录：

| 模块 | 职责 |
| --- | --- |
| `volume_jacobian.hpp/.cpp` | 体积分和时间项的系数缓存、matrix-free 作用、单元对角块 |
| `dg_operator.hpp/.cpp` | 复用体积分缓存，加上冻结 HLL 波速下的面导数；提供完整 DG 的 apply/block |
| `hybrid.hpp/.cpp` | 共用 VolumeJacobian 装填局部 A，再加 HDG 面项并进行静态凝聚 |

两个线性化对象的构造函数完成全部 Eigen 空间分配，必须在时间循环前构造。
`refresh(state)` 在已有工作区内更新缓存，不分配空间。对象不保存 state 的引用，
即使原 state 向量被原地修改或 swap，已有缓存仍表示上一次 refresh 的状态。
网格、基函数和质量矩阵在对象生命周期内必须保持不变。

DG 调用示例：

```cpp
DGOperator jacobian(problem);          // 时间循环前：分配全部缓存和工作区
// 每个时间步：problem.freeze(old_state)，然后设置 problem.step_size。
// 每次需要 Newton 修正、且残差尚未收敛时：
jacobian.refresh(newton_state);        // 同时快照当前 dt、体积分导数和面导数
preconditioner.refresh([&](int row, int col, BlockRef block) {
    jacobian.block(row, col, block);
});
// GMRES 内可重复调用，期间绝不 refresh：
jacobian.apply(direction, result);     // result 必须已分配，且不能与 direction 别名
```

`state`、`step_size` 或冻结波速发生变化后，下一次线性求解前必须重新 refresh。
`apply` 和 `block` 不会隐式刷新，也不比较向量指针或内容来猜测缓存是否有效；
首次 refresh 之前调用会报错。Newton 线搜索只计算候选状态的残差，不能在其中刷新
正在使用的线性化；接受候选状态后，在下一次 Newton 线性求解前刷新即可。
`TimeIntegrator` 已按这个顺序调用。直接参考和预条件器共用相同的缓存。
`HybridOperator::refresh` 内部调用自己的 `VolumeJacobian::refresh`，外部无需另行调用。

### 系数布局与 matrix-free 作用

令 `n=k+1`、`Q=basis.points`、`C=layout.cells`、`J_q=F'(U_q)`。
体积分只保存 `Q*NVAR*NVAR*C` 个 double，列主序逻辑张量为
`[Q,NVAR_flux,NVAR_state,C]`。同一个物理 Jacobian 分量在 q 方向连续，
从而用 Eigen array 运算计算 `z_q=J_q*x_q`。不展开成每单元 `n*n*NVAR*NVAR` 的矩阵。

`VolumeJacobian::apply` 先写入 `M*x/dt`，再逐单元计算
`x_q=V*x_nodal`、`z_q=J_q*x_q`、`result -= (D^w)^T*z`。
`(D^w)^T` 在构造时保存为连续列主序矩阵，避免反复访问原矩阵的跨步长行。
积分点状态、方向和 flux 工作区仅为一个单元的 `[Q,NVAR]`，反复复用。

DG 面导数分左右各存一个 `[NVAR_flux,NVAR_state,C]` 张量。
每次 refresh 每个面只计算一次 `G_left` 与 `G_right`；apply 使用
`G_left*x_left + G_right*x_right`，block 查询也直接复用它们。
面导数缓存约 `2*C*NVAR*NVAR` 个 double，与单元多项式阶数无关。

### 配点和超积分的两条路径

当前参考单元在 `Q==n` 时设置 `V=I`。此条件是此处配点特化的依据，
未来如果改变节点或积分规则，必须同时维护这一约束，不能只比较两个维数。

默认配点路径：

```text
U_q = U_nodal                         // 零拷贝读取，无插值
x_q = x_nodal                         // apply 中同样不插值
A[(v,i),(w,j)] = -D^w[j,i] * J_j[v,w]
A.diagonal += M/dt
```

因此体积分块装填从通用的 O(Q*n²*NVAR²) 简化为 O(n²*NVAR²)。
矩阵按 `w -> j -> v -> i` 的逻辑顺序写入，i 连续；所有体积分元素直接覆盖，
不先清零后反复累加。只有 n 个物理 Jacobian，完整块中的 n² 个节点耦合是
这些物理 Jacobian 与参考微分矩阵系数的组合。

超积分 `Q>n` 路径保留一般离散式。对每个 (v,w,j)，在已分配的 q 工作区中计算
`weights[q]=J_q[v,w]*V[q,j]`，再做连续输出列的矩阵向量积
`A[v_nodes,(w,j)]=-(D^w)^T*weights`。没有构造 Q×Q 对角矩阵或完整的三阶参考张量。
`Problem::quadrature_state` 的配点读取也已特化，残差、CFL 与正性检查共同受益。

本次优化保持现有逐非零块查询接口；没有把它描述成跨单元 SIMD 装填。
块内连续数组运算和 GEMV 可使用 Eigen 支持的向量化，物理 Jacobian 刷新仍逐积分点计算。
若后续实现跨单元 SIMD 或按单元批量填充块，需另外测量打包和重排序写入的代价。
工作区由对象复用，所以同一对象的 apply/block 不可被多个线程同时调用。

### 独立检查和小网格计时

原来的一般积分点求和实现保留在 `check/dg_reference.hpp/.cpp`，仅用于独立检查和
小网格基准，不链接进计算库。`check/linearization.cpp` 检查首次 refresh/apply/block
即禁止 Eigen 分配、同一 state 内存原地更新后快照不变，以及显式 refresh 后 dt、状态、
波速全部更新。原有有限差分、DG/HDG 凝聚一致性和时间推进检查继续执行。

`benchmark_dg` 比较无缓存的一般积分点参考实现与新实现，包括单次 Jv、所有体积分块、
refresh 及一次 refresh 加 20 次 Jv 的组合。网格仅 32 单元，k=0/1/3/6，
同时测试配点与超积分。它不是完整求解器的端到端加速测量；完整 Newton/GMRES 的
收益还取决于迭代次数、预条件器分解、缓存容量和硬件，不能直接套用单核函数倍数。

## Block ILU(0) 与 COLAMD

只对单元/界面连接图执行一次 COLAMD，得到块顺序；使用同一个顺序做对称行列置换，
不把 COLAMD 的列置换直接当作完整数值分解。数值更新保留置换后原块图的结构：

```text
Lij = right_solve(Ujj, Aij)  // Lij Ujj = Aij
Aik -= Lij Ujk              // 仅当 (i,k) 原本存在
```

对角块位置与每个下三角块所需的乘积更新位置在初始化时预计算；数值分解不查找稀疏图。
每个对角块在所有前序更新后做块内选主元 LU。右除通过 Uᵀ、Lᵀ、P 回代；
应用预条件器时做块级前代/后代和块内 LU 回代。没有显式逆、标量 ILU 或填充块。
block-Jacobi 只保存与分解对角块，使用同样的局部 LU 接口。

“matrix-free ILU”在此指预条件器直接向算子查询块，不要求组装标量全局稀疏矩阵；
不完全分解本身必然要保存块因子。空间复杂度为原图块数乘以块大小平方。

## GMRES、直接参考与空间分配

GMRES 是基于 Eigen 工作区实现的 restarted right-preconditioned GMRES。
两遍 modified Gram-Schmidt、Givens 旋转、真残差复核均使用预分配空间；
DG 使用 inverse-mass 加权内积，HDG trace 使用欧氏内积。
迭代解存储预条件后的 Krylov 向量，不将左预条件残差当作原方程残差。

直接路径从相同块接口装入预分配的全局稠密矩阵，允许跨单元行选主元。
它是小规模精度参照，避免将“不分配内存”的承诺建立在 Eigen SparseLU 的动态填充上。
不会在 GMRES 模式下分配此全局矩阵。

Problem 与 TimeIntegrator 构造完成所有 Eigen 空间分配，包括 Krylov 基、块因子、
LU 置换和回代空间、HDG 消元结果、Newton trial、残差和方向。
只构造所选求解器及方法的主要工作区。时间循环中不调用 resize、容器插入、
返回动态矩阵或动态向量的算子，也不写文件。
DenseLU 在预分配 Eigen 矩阵上做选主元、逐列秩一更新和显式回代，
避免通用 LU 的 GEMM packing 分配及 solve 表达式的隐含临时量。
所有检查程序和检查专用辅助函数均在 `dg_hdg_1d/check/`；运行时输入合法性和正性判定属于算法本身。

## 后续 2D 扩展边界

本次没有 2D 文件或求解实现。后续需要扩展几何 Layout、张量积体积算子、
四面 trace 及二维邻接图；不能只把 NVAR 从三改成四。
