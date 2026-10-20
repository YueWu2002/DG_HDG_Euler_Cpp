# 2D 精度和收敛阶

`bash accuracy_test/2d/submit.sh` 提交 Intel 单 CPU、16 GiB、一天时限的独立源码快照。
在计算节点构建并执行 CTest、生产路径检查后，顺序完成 32²、64²、128²、256² 网格。
k=2，SSP-DIRK(4,3)、CFL=2、T=0.2，gamma=3、alpha=45°，
Newton/线性相对容差均为 1e-12，restart=20、max-linear=1500、max-newton=30，
Gauss 配点及 HLL safety=1.1。DG 和 HDG 均完全收敛，使用精确局部块 Jacobi。

计算域为 [-1/cos(alpha),1/cos(alpha)] × [-1/sin(alpha),1/sin(alpha)]，
与参考项目一致。输出四个守恒分量的 L2 误差、同步时空加密阶、守恒漂移，
DG/HDG 各分量积分 L2 差、联合 L2 差和最大节点差。
节点解文件默认在扫描中关闭，误差积分和比较仍完整执行。

```sh
python3 accuracy_test/2d/convergence.py --exe build-2d-intel/dg_hdg_2d/dg_hdg_2d \
  --output accuracy_test/running/local_unique
```

`--grids 32 64` 可运行独立的部分扫描；SUCCESS 仅说明所请求的网格完成，
完整四网格结论必须检查 convergence.csv 的 nx/ny。默认四网格扫描不跳过历史结果。
局部块求解固定为精确 LU，不提供近似块逆选项。

新任务统一保存到测试根目录的 `running/<维数>_submission_*/`；结束后运行 `python3 scripts/archive_finished.py` 整体归档到 `results/`。每次运行保留输出、后处理数据、CPU、参数和代码版本；详见根目录 `AGENTS.md`。旧的运行中目录不搬移。

本地运行同样使用顶层 `running/`；省略 `--output` 自动创建带维数、时间戳和随机后缀的目录。显式 `--output` 也必须位于相应测试的 `running/` 下。Slurm 仍通过 `submit.sh` 提交。
