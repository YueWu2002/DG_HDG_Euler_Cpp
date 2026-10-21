# 1D / 2D 验证协议

1D/2D 实际结果与作业状态统一见 [RESULTS.md](RESULTS.md)。
本文件规定测试，不把计划或历史结果写成已完成实验。

## 正确性

保留 `dg_hdg_1d/check/` 的独立体积分参考、有限差分、凝聚等价性、
LU、加权 GMRES 真残差及首次调用 Eigen 禁分配检查。
新增 DIRK3 配对推进检查和相同 CFL 下的时空同步加密；BE 检查显式指定积分器。
HDG 每次更新以恢复后的 DG 真残差验收，独立 Jv 诊断验证残差恢复公式。
日常生产路径检查入口为 `python3 code_test/1d/test_solver.py --exe <明确路径> --output code_test/running/<独立目录>`，
包含 DIRK3 smoke 和单步完全收敛 BE 的 DG/HDG 比较，并独立验证恢复残差恒等式。
默认生产构建只编译，不自动执行测试。

## 收敛

`accuracy_test/1d/convergence.py` 使用 32、64、128、256、512 单元，k=2，
SSP-DIRK(4,3)，CFL=2，T=0.2，HLL safety=1.1，
Newton 和线性相对容差各 1e-12，max Newton=30。
每个网格从投影初值开始分别完全收敛 DG 与 HDG。
记录三个分量 L2 误差、log2(e_h/e_h2)、守恒漂移、Newton/Krylov 数量、
实际 min/max/mean dt、DG/HDG 的积分 L2 差与最大节点差。
此为同步时空加密，k=2 对应预期渐近三阶。

## 性能

`speed_test/1d/benchmark.py` 在同一分配节点、一个 CPU 上顺序运行 DG/HDG，
512 单元、k=0..8、matrix-free、block-Jacobi、BE、CFL=5、steps=1，
固定四次 Newton 更新，第一更新预热，后三次平均。
保留每次更新的时间与线性/DG 真残差、完整循环时间、迭代数、资源信息。
报告 DG 平均/HDG 平均的比值；不以自由度之比或单核函数时间代替端到端测量。
只有实际测到 k=8 比值 >=8 才能宣称达到目标。

Slurm 精度与速度独立提交，不自动设置相互依赖。后续收敛批次按网格拆成独立任务，内存分别估算，具体入口规范见 accuracy skill。
新任务在对应顶层 running 下独立运行，结束后完整归档至 results，不覆盖先前实验。

## 2D 独立验证

见 `code_test/2d/README.md` 的独立数学检查，另运行旋转简单波、BE/DIRK3、
固定旋转算例与熵波的生产路径 DG/HDG 比较，不提供命令行旋转角输入。
收敛网格固定 32²、64²、128²、256²，Q2、DIRK3、CFL=2、T=0.2、容差均 1e-12。
速度在 Slurm 上运行 256²、k=8→0，100 GiB，四次 Newton 更新首个预热。
单元/面块均使用完整块的精确 LU，张量积运算通过独立求积和块/算子恒等式验证代数等价。
1D 历史实验与 2D 状态、完成程度分别记录于 `RESULTS.md` 对应章节。

## 本次已完成实验的残差解释

作业7640912和7635952的完成结果见 `RESULTS.md` 的 2D 章节。HDG Krylov 内部控制逆面质量加权 Schur 相对残差；结束后通过代数映射验收恢复后的 DG 残差。
验收阈值为 `max(10*ltol*||R||_(M^-1), allowance)`，并不是对恢复后的 DG 相对残差严格施加 ltol。1D HDG Krylov 内部的权重为全1，恢复后的 DG 验收采用相同阈值形式。
因此 `linear_relative_residual` 和 `dg_linear_relative_residual` 必须分别记录，不能宣称两方法 Krylov 使用完全相同的 DG 度量。独立 `--check-update` 检查用于验证恢复恒等式；正式 timing-only 运行不执行该独立参考求解。


## Debug验证与Release实验构建

Release以`NDEBUG`及`EIGEN_NO_DEBUG`关闭Eigen运行时断言，包括禁止Eigen堆分配的断言；项目`require(...)`物理/残差检查继续生效。Debug保留尺寸、索引和分配检查，不能据此保证未覆盖输入的正确性。
`intel-build.sh`支持`BUILD_TYPE=Debug|Release`（默认Release）；默认仅Debug开启BUILD_TESTING。code_test驱动要求可追溯到CMakeCache.txt的Debug二进制，拒绝Release或无法核实的构建。
共享`scripts/build_debug_checks.sh <1d|2d> <运行目录>`在`build-debug/`构建Debug，运行CTest及生产路径code_test；实验入口随后在独立`build/`构建Release，供accuracy/speed使用。保留Debug构建日志、CTest、code_test及provenance，也保留Release配置/二进制身份，不将两套二进制混称为同一验证对象。新增构建成本需计入提交前资源预估，不改写在途快照。

本地code_test默认在 `code_test/running/` 创建唯一目录，显式输出也限于此处；Slurm内嵌检查由 `DGH_CHECK_ROOT` 固定到所属任务 `checks/`，不拆成第二份档案。本地驱动保留 `exit_code.txt`，使用 `scripts/archive_finished.py --run-root <目录>` 在确认进程结束后归档到 `code_test/results/`。脚本布局说明见 [scripts/README.md](scripts/README.md)。
