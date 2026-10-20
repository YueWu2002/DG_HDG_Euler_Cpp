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
python3 speed_test/1d/benchmark.py --exe build/dg_hdg_1d/dg_hdg \
  --degrees 8 --output speed_test/running/local_unique
```

新任务统一保存到测试根目录的 `running/<维数>_submission_*/`；结束后运行 `python3 scripts/archive_finished.py` 整体归档到 `results/`。每次运行保留输出、后处理数据、CPU、参数和代码版本；详见根目录 `AGENTS.md`。旧的运行中目录不搬移。

本地运行同样使用顶层 `running/`；省略 `--output` 自动创建带维数、时间戳和随机后缀的目录。显式 `--output` 也必须位于相应测试的 `running/` 下。Slurm 仍通过 `submit.sh` 提交。
