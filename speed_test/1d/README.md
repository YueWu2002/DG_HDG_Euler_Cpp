# 1D Slurm 速度测试

```sh
bash speed_test/1d/submit.sh --dependency afterok:ACCURACY_JOB_ID
```

512 个单元，k=0..8，DG/HDG 独立进程顺序运行，单 CPU、单线程，
BE、CFL=5、一个时间步内固定四次 Newton 更新；每个方法第一个更新为预热，
只对后三个 iteration_s 求均值。Newton/线性容差 1e-12、restart=20、max-linear=1500，
matrix-free 和右 block-Jacobi。完整求解循环时间仍包含全部四次更新。

提交时冻结源码，在 Intel 计算节点构建，保存环境、编译日志、哈希、每次更新时间和
资源数据。data/speedup.csv 给出 DG_mean/HDG_mean；timings.csv 包含原始统计。
k=8 的 8 倍是性能目标，测试不会把未达到目标的结果包装成成功加速。
SUCCESS 仅说明实验执行完毕，不代表达到 8 倍，也不代表四步 Newton 完全收敛。

仅运行指定次数：

```sh
python3 speed_test/1d/benchmark.py --exe build/bin/dg_hdg.exe \
  --degrees 8 --output speed_test/running/local_unique
```

新任务统一保存到测试根目录的 `running/<维数>_submission_*/`；结束后运行 `python3 scripts/archive_finished.py` 整体归档到 `results/`。每次运行保留输出、后处理数据、CPU、参数和代码版本；详见根目录 `AGENTS.md`。旧的运行中目录不搬移。

本地运行同样使用顶层 `running/`；省略 `--output` 自动创建带维数、时间戳和随机后缀的目录。显式 `--output` 也必须位于相应测试的 `running/` 下。Slurm 仍通过 `submit.sh` 提交。

## restart 对照的独立任务提交

```sh
python3 scripts/speed_campaign.py prepare --dimension 1d --restart 100 --minimum-memory-gib 16
python3 scripts/speed_campaign.py submit --run-root <上一命令生成的目录>
```

该对照使用SHA256与旧基线一致的二进制，每个k和方法一个Slurm任务；本次1d指定原节点。计算任务各至少16 GiB，较大任务根据历史RSS加新增基存储和余量申请。计算结束后afterany依赖汇总，再用独立依赖任务从目录外整体归档。原始数据在各job_ID/data中，Slurm out/err在scheduler中，汇总表在postprocessing中，报告在report_archived中。所有信息在同一实验文件夹保留。只在Slurm计算节点运行求解器。

1D 自定义网格与同节点 restart 对照：`prepare --dimension 1d --restart 20|100 --cells 2048 --node <已核实节点> --minimum-memory-gib 16`。两套均覆盖 k=0..8 的 DG/HDG，复用同一基线二进制；1D 历史源码身份未知，以二进制完整 SHA256 为准。内存估算按单元数放大历史峰值并增加 Krylov 基存储。网格与历史基线不同时，不计算历史耗时比；两套新结果可按匹配参数比较。
