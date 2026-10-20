# 1D / 2D 验证协议

1D 实际结果与作业状态见 [RESULTS.md](RESULTS.md)，2D 见 [RESULTS_2D.md](RESULTS_2D.md)。
本文件规定测试，不把计划或历史结果写成已完成实验。

## 正确性

保留 `dg_hdg_1d/check/` 的独立体积分参考、有限差分、凝聚等价性、
LU、加权 GMRES 真残差及首次调用 Eigen 禁分配检查。
新增 DIRK3 配对推进检查和相同 CFL 下的时空同步加密；BE 检查显式指定积分器。
HDG 每次更新以恢复后的 DG 真残差验收，独立 Jv 诊断验证残差恢复公式。
日常生产路径检查入口为 `python3 code_test/1d/test_solver.py --exe <明确路径> --output <独立目录>`，
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

Slurm 速度作业依赖收敛作业 afterok，保证先精度后性能。
所有结果写进各测试目录下的独立 job_ID，不覆盖先前实验。

## 2D 独立验证

见 `code_test/2d/README.md` 的独立数学检查，另运行旋转简单波、BE/DIRK3、
30° 旋转与熵波的生产路径 DG/HDG 比较。
收敛网格固定 32²、64²、128²、256²，Q2、DIRK3、CFL=2、T=0.2、容差均 1e-12。
速度在 Slurm 上运行 256²、k=8→0，100 GiB，四次 Newton 更新首个预热。
单元/面块均使用完整块的精确 LU，张量积运算通过独立求积和块/算子恒等式验证代数等价。
1D 历史实验仍单独保留；2D 状态和完成程度见 `RESULTS_2D.md`。
