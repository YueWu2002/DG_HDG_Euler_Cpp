# 1D 收敛测试

`bash accuracy_test/1d/submit.sh` 提交独立源码快照；一台 Intel 节点、一个 CPU、
60 GiB、一天时限，节点内 Intel/MKL/Eigen Release native 构建并运行 CTest。
随后逐网格运行 DG/HDG，32/64/128/256/512，k=2、DIRK3、CFL=2、T=0.2、
Newton/线性容差 1e-12、restart=20、max-linear=1500、max-newton=30、HLL safety=1.1。

输出在 `accuracy_test/running/1d_submission_时间_随机值/job_JOBID/`，data 下有每个网格的 summary、
Newton 记录、节点解、DG/HDG 积分 L2 差和最大节点差、convergence.csv。
收敛表逐网格更新；仅存在 SUCCESS 文件才表示五个网格全部成功。
编译日志、环境、源码清单和可执行文件哈希保留，历史结果不覆盖。

也可以显式指定可执行文件本地运行：

```sh
python3 accuracy_test/1d/convergence.py --exe build/bin/dg_hdg.exe \
  --output accuracy_test/running/local_unique
```

新任务统一保存到测试根目录的 `running/<维数>_submission_*/`；结束后运行 `python3 scripts/archive_finished.py` 整体归档到 `results/`。每次运行保留输出、后处理数据、CPU、参数和代码版本；详见根目录 `AGENTS.md`。旧的运行中目录不搬移。

本地运行同样使用顶层 `running/`；省略 `--output` 自动创建带维数、时间戳和随机后缀的目录。显式 `--output` 也必须位于相应测试的 `running/` 下。Slurm 仍通过 `submit.sh` 提交。
