---
name: optimize
description: 优化本项目 Eigen Euler DG/HDG 1D/2D 的精确算子与局部求解，利用张量积、缓存友好布局、SIMD 和存储复用，并验证正确性及实测性能。
---

先读根目录 `AGENTS.md`，再查 `DESIGN.md`、`MATHEMATICS.md`、`MEMORY_ESTIMATE.md`。对照项目是 `../mfem_euler2d_cpp`；本项目目标是在保持协议一致的前提下以 Eigen 自行实现更高性能，不链接 MFEM，不改变 common/dg_hdg_1d/dg_hdg_2d 的架构。

## 目标与边界

以已有效的 1D 算法为核心，通过精确张量积扩展到 2D，优先利用质量对角、matrix-free、紧凑块存储、SIMD 与 cache friendliness。目前正式实验只使用完整局部块精确 LU 的 block-Jacobi 和固定右预条件 GMRES。

用户明确排除张量方向近似、ADI、近似张量逆及类似替代。变量系数 Euler 单元块通常不能精确分离：导数、限制和提升可 sum-factorize，不代表整块逆可写成一维逆的张量积。任何减少存储的结构假设必须由实际矩阵的代数结构证明。

1D k=8 的至少8倍 DG/HDG 加速是目标，不是每台 CPU 或2D的保证。不得通过放松容差、改求积、降低精度、改稳定化、改变工作负载或计时边界制造收益。

## 定位与实现

1. 依据实际阶段计时、RSS 或 profiling 选择瓶颈。保存基线版本、编译器、CPU、线程、网格、阶数、CFL、容差和 Krylov 次数；保留在途任务所需路径，另用独立构建。
2. 查 `common/src/linalg/gmres.cpp`、`common/include/euler/linalg/{block_jacobi,dense_lu}.hpp`，以及 `dg_hdg_{1d,2d}/src/{discretization,hybrid,time_integrator}.cpp`。2D 张量和通量内核在 `dg_hdg_2d/src/{tensor_product,volume_jacobian,face_jacobian}.cpp`。
3. 优先消除稠密块副本、重复分解/求解、热循环分配及 gather/scatter。使用连续块池、原位 pivoted LU、预分配 scratch、多 RHS 和稳定的非拥有视图。明确每个缓存依赖状态、dt、冻结波速还是仅依赖网格/阶数。
4. 张量收缩沿连续方向遍历，复用一维导数表与中间结果，减少大 Kronecker 矩阵和不必要转置。对齐、存储顺序、循环次序和小矩阵内核要通过测量选择，不能只凭 SIMD 标记宣称加速。
5. 采用流式单元凝聚、紧凑面耦合或省去缓存时，逐项核对 Schur 作用和恢复是否还需要原矩阵或 X=A⁻¹B，量化重算成本。不能覆盖仍被 matrix-free 作用读取的数据。

固定线性右预条件器允许在组合 Krylov 更新后再作用一次，避免保存第二套 Z 基；若预条件器在一次 GMRES 内变化，该推导失效。保留逆质量加权内积、真实残差重算和恢复后的共同 DG 验收。

共享 scratch 和分配守卫不自动具备线程安全；不能直接给单元循环加并行 pragma。保持正式单线程协议，若另研究线程扩展则单独标记实验。

## 内存估算

2D 令 q=k+1、E=N²、单元宽 m=4q²。DG 稠密 Jacobi 因子池约 `8 E m²` 字节，体积向量 `8 E m` 字节，GMRES(r) 基约 `(r+1)` 个向量。256²、k=8 时仅因子池约51.26 GiB。

用 `speed_test/2d/memory_estimate.py` 核对 HDG X、Schur、恢复、面因子和系数缓存；再算所有同时存活的阶段状态、临时矩阵、pivot、索引、库和分配器开销。matrix-free 不等于没有大预条件块；100 GiB 是该正式速度实验的参考申请，不是实际峰值或所有配置的保证。

## Release优化与Debug验证

- 正式性能构建使用Release：保留`-O3`、native指令、支持时的IPO及SIMD，以及Eigen/MKL BLAS/LAPACKE路径；通过`NDEBUG`和`EIGEN_NO_DEBUG`关闭Eigen运行时断言，不再对Release强制`-UNDEBUG`。不因此关闭项目`require(...)`的物理可容许性、求解失败和残差验收，也不自动启用fast-math或放宽数值协议。
- code_test及相关CTest使用独立Debug构建，保留Eigen尺寸/索引断言和`EIGEN_RUNTIME_NO_MALLOC`分配检查；`intel-build.sh`通过`BUILD_TYPE=Debug`选择真实Debug配置，目录名称不能代替配置。code_test入口核验CMakeCache并拒绝Release或无法核实的构建。
- 实验入口先用`scripts/build_debug_checks.sh <1d|2d> <运行目录>`在`build-debug/`构建、检查，再在`build/`独立构建Release用于accuracy/speed。保留两套配置、构建日志、Debug检查输出及对应二进制SHA256/provenance，不能将Debug检查冒充为Release二进制验证。
- Release中Eigen无分配断言同时失效，不能宣称该性能进程仍在运行时检查无分配；Debug测试通过也不证明所有输入或Release数值行为正确。关闭断言的实际加速需同配置、同设备对照测量，不能仅凭编译选项宣称收益。
- 两次构建和Debug检查的额外耗时须在提交前估算并单独申请；不改写在途快照，不把历史开启断言的测量与新Release测量混为同一构建配置。

## 验证与计时

先按 [debug](../debug/SKILL.md) 运行当前二进制的相关 CTest 与短验证；数学路径改动再按 [accuracy](../accuracy/SKILL.md) 检查精度和收敛阶。不能借优化任务修改正在运行的源码快照。

正式测速按 [speed](../speed/SKILL.md)：1D 512 单元、2D 256²，四次 Newton 更新、第一次预热、后三次平均，单线程 Slurm。额外重复或 profiling 单独记录，不替代该均值。对照同一方法优化前后及 DG/HDG 比值，同时报告 setup、迭代数、时间和 RSS，区分省内存与省计算的收益。

本地和 Slurm 新结果均落在对应 `running/`，结束后完整归档至 `results/`。更新受影响的设计说明，只有可追溯的完成测量才进入结果文档；保留所有 CSV 和设备日志。

新构建产物与测试入口遵循 [脚本目录约定](../../scripts/README.md)：可执行文件为各构建 `bin/*.exe`，静态库为构建目录内的 `*.a`；只忽略二进制产物，不忽略构建日志及配置。不改写在途快照或历史基线路径。
