# 1D / 2D 架构与性能设计

## 模块与边界

继续使用现有目录，不引入 MFEM 或外部 Krylov 求解器。
`dg_hdg_1d/src/physics.cpp` 实现 Euler/HLL；`discretization.cpp` 管残差、CFL、
投影及积分；`volume_jacobian.cpp` 缓存通量导数；`dg_operator.cpp` 提供 Jv 与单块；
`hybrid.cpp` 负责静态凝聚；`time_integrator.cpp` 实现 BE 与四阶段 DIRK3；
`main.cpp` 处理输入和时间循环之外的输出。
`common/src/linalg/` 是 Eigen 工作区上的右预条件 GMRES、块图和预条件器。
`dg_hdg_2d/` 采用相同职责划分，二维张量结构和存储策略见末节。

## 1D 布局和算子

状态 Eigen 列主序 `[k+1,1,NVAR,Nx,1]`，节点连续、其次变量、最后单元；
`Layout::cell` 为零拷贝 `[k+1,NVAR]` 视图。
质量只存对角向量。体积系数 `[Q,NVAR_flux,NVAR_state,Nx]` 的 Q 维连续，
Jv 通过点上物理 Jacobian 和一维微分收缩，避免每单元完整 Jacobian 乘法。
配点 Q=k+1 时跳过恒等插值。Eigen 连续 array 和 GEMV 使用 CPU 向量化；
未声称实现跨单元 SIMD。超积分保留一般插值路径。

DG 全局迭代不组装标量稀疏矩阵，只为 block-Jacobi 查询单元对角块。
HDG 局部选主元 LU 获得 X=A^-1 B、y=A^-1 r，缓存 -CX 的面耦合块；
Schur 作用只做 gather、块乘、scatter，不在每次 Krylov 乘法重复 C(Xv)。
1D 的 6×6 Schur 作用和 3×3 trace 三角回代使用固定大小 Eigen Map，
让小块运算内联，避免逐单元动态 BLAS 调用。
每次 Newton 重新投影 trace；DG 残差用于共同线搜索与停止判据。
独立的增广 trace 接口仍可用于凝聚等价性诊断。

GMRES 用加权两遍重正交与真残差复核，固定右预条件器允许将 Arnoldi 向量组合后
只作一次回代更新解，不必存储全部预条件 Krylov 向量。
block-Jacobi 为每个单元/面独立 LU，不需要 COLAMD；
block ILU(0) 保留块图和可选 COLAMD，没有标量填充。

## 缓存与分配

构造函数分配状态、阶段斜率、Krylov 基、块因子及线搜索空间。
refresh 在每次 Newton 线性化时更新快照；apply/block 复用该快照，
改变状态、步长或冻结波速后必须显式 refresh。同一算子的可变工作区不可并发使用。
时间推进从第一次调用开始禁止 Eigen 堆分配，Release 仍启用相应断言。
此保证只覆盖 Eigen 管理的分配，不宣称拦截标准库或系统全部分配。

## 计时

主指标是接受的 Newton 更新耗时，包括共同 DG 残差、线性化、预条件器建立、
凝聚、Krylov、真残差验收、恢复和线搜索。每次更新时间保存在预留容量的记录容器中（超出估算时仅标准库容器可扩容），
全部 CSV 在时间循环结束后写出，不把文件写出算进求解耗时。
速度运行 512 单元、BE 单步四次更新；第一次预热，后三次平均。
DG 和 HDG 各自预热，独立进程顺序运行，单线程；同时报告完整四次 time_loop_s。
固定更新次数不表示非线性收敛。收敛实验与性能实验严格分开。

## 构建与实验隔离

Intel/MKL LP64 sequential，-march=native，SIMD 编译选项，Release IPO，
不使用 fast-math。保留无 BLAS/LAPACK 的 Eigen 原生构建以作可移植验证。
Slurm 提交时复制源码到本次独立快照，计算节点单独构建，原目录后续编辑不改变在途作业。
不重编译或修改参考项目的在途作业依赖。

## 2D：从 1D 经验扩展

2D 已有独立 CMake 入口、main、时间推进、HDG 和正确性测试，1D 文件结构保持不变。
Eigen 列主序状态为 `[nx_node,ny_node,NVAR,Nx,Ny]`，x 节点连续。
Gauss 配点使用固定 n=1..9 的参考收缩与面限制/提升内核；超积分保留两方向一般收缩。
全局缓存体积 u/v/H 三场，面缓存两侧 u/v/H 与冻结 HLL 系数，Jv 不读取完整单元矩阵。

精确 DG 块 Jacobi 使用单个对角装填工作区和持久 LU 池，避免第二份全局对角矩阵。
HDG 不保留逐单元 A、B、C 或 LU：一个单元 m×m A、一个 LU 和 m×(t+1) 右端工作区复用。
X、y 和局部 Schur 为连续大数值池；面限制先对多个右端形成 fw×(t+1) 小工作区，
再沿连续的 Gauss 点轴施加 4×4 物理导数，避免稠密 C 的大量结构零乘法。

块 Jacobi 的 LU 系数和置换均连续保存，使用一套 LU 工作区逐块刷新。
DG 单元块和 HDG 面块均保留完整耦合并作精确 LU 求解，分别在 CSV 中标记为 dense、face-dense。
张量积优化保持代数等价，不使用方向乘积分裂、低秩截断或局部迭代替代完整块求解。

2D Slurm 使用独立源码快照、一个 CPU；Q2 精度扫描 16 GiB，速度扫描 100 GiB。速度 256²、k=8→0，四次更新首个预热，
不与 1D 的 512 单元结果混用，也不将降低局部求解精度包装为同一预条件器。
