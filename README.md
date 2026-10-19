# 1D Euler：分块 DG / HDG 重构

本目录是独立实现；`../DG_HDG_Euler_Cpp` 仅作为数学参照，不依赖或修改原目录。
本阶段只实现 1D 周期 Euler、隐式 Euler 时间离散、冻结 Einfeldt/HLL 波速，以及
代数等价的 DG 与 augmented-trace HDG。没有开始 2D 实现。

## 文件框架

```text
CMakeLists.txt                 继承原项目顶层结构，仅启用 1D
cmake/                        继承原项目编译选项、BLAS/LAPACK 配置
vendor/fastgl/                原项目 FastGL，始终参与构建与链接
dg_hdg_1d/CMakeLists.txt       1D 独立构建入口
dg_hdg_1d/include/euler/
  types.hpp                    constexpr NVAR、物理分量名称、张量 Layout
  config.hpp                   运行选项
  physics.hpp                  Gas、Euler 通量/解析导数、Einfeldt 波速
  simple_wave.hpp              算例持有 gamma，解析解与光滑时间范围
  nodal_basis.hpp              仅初始化时计算的参考单元，无 basis.cpp 包装层
  discretization.hpp           质量张量、DG 残差、状态插值与积分
  volume_jacobian.hpp          积分点 Jacobian 缓存、配点特化、体积分块与 Jv
  dg_operator.hpp              Newton 状态快照、缓存面导数、完整 DG 的 block/apply
  hybrid.hpp                   HDG 局部块、静态凝聚、Schur 算子及恢复
  time_integrator.hpp          预分配工作区、CFL、Newton、正性线搜索
  linalg/
    block_graph.hpp            单元/界面级 CSR 图、块矩阵、块图 COLAMD
    dense_lu.hpp               块内选主元 LU、左/右回代，均不求逆
    block_ilu.hpp              分块 ILU(0) / block-Jacobi
    gmres.hpp                  Eigen 工作区上的 matrix-free GMRES
    direct.hpp                 全局选主元 LU，供小规模准确参考
dg_hdg_1d/src/                对应模块的 .cpp；main.cpp 只处理运行与输出
  linalg/                     分块线性代数实现
dg_hdg_1d/check/
  main.cpp                     独立检查程序入口
  check_1d.hpp/.cpp             检查代码完全独立于计算库
  linearization.cpp            缓存生命周期、首次调用无 Eigen 分配检查
  dg_reference.hpp/.cpp        独立的一般积分点参考实现
  benchmark_dg.cpp             小网格的系数刷新、块装填与 Jv 计时
scripts/
  build.ps1                    没有 CMake 时使用 MinGW 的本地构建入口
DESIGN.md                      张量次序、方程、预条件及内存约束
```

## 构建与运行

依赖为 C++17、Eigen >= 3.4，以及随项目保留的 FastGL。
COLAMD 使用 Eigen 自带实现，不需要外部稀疏求解库。
原项目的 `DGH_BLAS_BACKEND`、`DGH_LAPACK_BACKEND` 和 `DGH_NATIVE_ARCH` 选项继续保留。
Eigen 可通过安装包提供 `Eigen3Config.cmake`，或把 Eigen 3.4.0 解压到
`.deps/eigen-3.4.0/`（此目录不纳入源代码）。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Windows 当前环境也可使用：

```powershell
./scripts/build.ps1
./build/check_1d.exe
./build/benchmark_dg.exe
./build/dg_hdg.exe --cells 32 --k 3 --cfl 0.2 --T 0.1 --method both
./build/dg_hdg.exe --cells 8 --k 2 --cfl 0.2 --T 0.03 --linear direct --prefix output/reference
./build/dg_hdg.exe --cells 32 --k 3 --steps 10 --precond block-jacobi --ordering natural
```

`--linear gmres`（默认）和 `--linear direct` 均支持 DG / HDG。
CMake 构建的程序位于 `build/dg_hdg_1d/`（多配置生成器可能另有 Release 子目录）。
直接路径使用全局稠密选主元 LU，**内存 O(n²)、分解 O(n³)，仅用于小规模参考**；
HDG 的 n 是界面 DOF 数。此选择保证直接路径也能在时间循环内禁止 Eigen 堆分配。
GMRES 路径没有全局标量稀疏矩阵或全局稠密矩阵。
DG 每次 Newton 线性求解前刷新一次积分点与面导数，所有 GMRES 作用及块查询复用该快照。
默认配点时不做恒等插值，体积分块直接按微分矩阵系数装填；`--quad` 超积分仍受支持。
模块调用顺序、缓存失效条件及辅助张量布局见 `DESIGN.md`。

时间步只由 `--cfl` 决定，可用 `--T` 或 `--steps` 停止，两者不能同时指定。
`--dt`、`--increment-state`、普通 `ilu`、外部稀疏求解后端均不存在。
`gamma` 只在 `SimpleWave` 中设置一次，不是状态、命令行或 CSV 字段。

输出为 `<prefix>_summary.csv` 和 `<prefix>_solution.csv`。
summary 记录 CFL、求解器、迭代统计、自动选出的步长范围、L2 误差和守恒漂移。
所有输出在时间循环结束后写出；`--help` 显示支持的参数。

## 验证

`check_1d` 检查张量步长、连续块数值池、ILU(0) 块数固定、Gauss 精度、基函数导数、强制换行选主元的左右 LU 回代、
block ILU 的无缺失填充情形、COLAMD、GMRES 真残差、DG Jacobian 有限差分、
DG/HDG Newton 修正等价性、Schur 块/算子一致性、非投影 trace 的完整 HDG 有限差分、
时间推进中的守恒、解析解时间收敛，以及 DG/HDG/direct/GMRES 的最终状态。

计算库启用 `EIGEN_RUNTIME_NO_MALLOC`；默认时间推进从第一步起禁止 Eigen 堆分配，
而非先执行一次迭代后再检查。Release 构建也保留这些断言。
这验证 Eigen 可控制的分配，不宣称拦截操作系统或 C++ 标准库的所有分配。

本次验证后已按存储空间要求清理临时依赖、下载归档、构建产物和测试输出。
后续重新编译时需提供 Eigen 与 CMake；FastGL 源码始终保留在 `vendor/fastgl/`。
