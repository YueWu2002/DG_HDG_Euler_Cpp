# 2D 单线程速度测试

正式协议：256×256，k=8,7,6,...,0，BE、CFL=5、单个时间步固定四次 Newton 更新，
每个方法第一次预热，只统计后三次均值。右 block-Jacobi，Newton/线性容差 1e-12，
restart=20、max-linear=1500，HLL safety=1.1，Gauss 配点、alpha=45°。
DG/HDG 在同一计算节点、同一 CPU 依次运行，全部阶段建立和恢复成本均计入。

```sh
bash speed_test/2d/submit.sh
```

Slurm 申请 100 GiB、一个 Intel CPU、一天。在本次独立 source 快照内原地构建，
保存源文件哈希、二进制哈希、编译器、CPU、逐更新记录和每进程峰值 RSS。
原有 1D 和 MFEM 的在途作业目录不会被覆盖。

DG 单元块和 HDG 面块固定使用完整块的精确 LU 逆作用。
主路径不提供近似块逆或方向分裂选项；张量积只用于代数等价的算子计算。

```sh
python3 speed_test/2d/memory_estimate.py --cells 256 --k 8
```

256²、k=8 时，DG 精确 LU 数值池约 51.26 GiB；
HDG 的 X、局部 Schur、面 LU 约 22.78、10.13、1.27 GiB。
再加 Krylov 基和状态工作区，预算仍在 100 GiB 内；最终以资源文件的 RSS 为准。
不假定二维加速比与一维相同，也不假定 HDG 对所有次数更快。

完整实验结束后可自动验收并归档报告：

```sh
bash scripts/submit_2d_report.sh <accuracy submission/job_ID> <speed submission/job_ID>
```

此小作业依赖两个实验均成功，检查完整网格、次数、四次更新及均值计算，
只向新的报告目录写入，不覆盖原始作业。当前提交 ID 与目录见 [RESULTS_2D.md](../../RESULTS_2D.md)。

新任务统一保存到测试根目录的 `running/<维数>_submission_*/`；结束后运行 `python3 scripts/archive_finished.py` 整体归档到 `results/`。每次运行保留输出、后处理数据、CPU、参数和代码版本；详见根目录 `AGENTS.md`。旧的运行中目录不搬移。

本地运行同样使用顶层 `running/`；省略 `--output` 自动创建带维数、时间戳和随机后缀的目录。显式 `--output` 也必须位于相应测试的 `running/` 下。Slurm 仍通过 `submit.sh` 提交。
