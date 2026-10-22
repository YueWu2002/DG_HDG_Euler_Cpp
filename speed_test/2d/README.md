# 2D 单线程速度测试

默认入口协议：256×256，k=8,7,6,...,0，BE、CFL=5、单个时间步固定四次 Newton 更新，
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

新任务统一保存到测试根目录的 `running/<维数>_submission_*/`；结束后运行 `python3 scripts/archive_finished.py` 整体归档到 `results/`。每次运行保留输出、后处理数据、CPU、参数和代码版本；详见根目录 `AGENTS.md`。旧的运行中目录不搬移。

本地运行同样使用顶层 `running/`；省略 `--output` 自动创建带维数、时间戳和随机后缀的目录。显式 `--output` 也必须位于相应测试的 `running/` 下。Slurm 仍通过 `submit.sh` 提交。

## 统一后处理

```sh
python3 speed_test/2d/postprocess.py --run-root speed_test/results/2d_submission_20261019T205542_KhFC6z
```

`--run-root` 是实验目录。与另一类实验共用 `scripts/postprocess.py`：验收原始数据、重算指标，然后写入 `postprocessing/speed_summary.csv` 和唯一的 `report_archived_<日期>/REPORT.md`、`inputs.json`。CSV保留参数、设备、版本与任务信息；inputs.json记录相对输入路径、输入/输出SHA256及覆盖范围。报告只引用原始数据，不复制整套CSV和日志；重复执行复用已有报告目录，原始文件不变。

精度按原始summary核验并重算收敛阶；速度按四次更新、首个预热和后三次平均核验，并记录各方法总时长、RSS、Krylov数及残差。自定义部分扫描明确标记覆盖范围，不等于完整正式实验。

旧联合入口 `scripts/collect_2d.py` 只调用两类后处理，不再复制数据或创建联合报告；`--output` 已废弃。两类任务与目录规则详见各自SKILL。

## restart 对照的独立任务提交

```sh
python3 scripts/speed_campaign.py prepare --dimension 2d --restart 100 --minimum-memory-gib 16
python3 scripts/speed_campaign.py submit --run-root <上一命令生成的目录>
```

该对照使用SHA256与旧基线一致的二进制，每个k和方法一个Slurm任务；本次2d指定原节点。计算任务各至少16 GiB，较大任务根据历史RSS加新增基存储和余量申请。计算结束后afterany依赖汇总，再用独立依赖任务从目录外整体归档。原始数据在各job_ID/data中，Slurm out/err在scheduler中，汇总表在postprocessing中，报告在report_archived中。所有信息在同一实验文件夹保留。只在Slurm计算节点运行求解器。


## 新源码二维速度独立作业入口

`python3 scripts/prepare_snapshot_speed.py --node <探测节点> --cells <每方向单元数> --restart <GMRES重启长度> --factor <资源倍率>`，随后对生成目录执行 `python3 scripts/speed_campaign.py submit --run-root <目录>`。
入口优先用MEMORY_ESTIMATE.md的二维256²/GMRES(20)实测行估算，保存引用片段和文档哈希；内存先按单元数缩放，再加`8*全局线性自由度*(restart-20)`的新增基存储，最后加余量；时间按历史完整进程wall乘倍率，restart变化的迭代/正交化成本另标为不确定。倍率需覆盖单元数比例，不能当作任意新协议的通用时长保证。
所有(k,method)独立提交到同一指定节点，共享构建成功后各自调度；共享构建保留必要依赖，汇总/归档亦保留依赖。源码快照在该计算节点编译并检查，记录实际CPU及二进制身份；结束后保留全部档案并移除临时source/。可选节点探测作业记录在同批任务清单中。

新快照入口支持`--build-time HH:MM:SS`，用于共享Debug验证加Release构建作业，单独记录预算；与每个算例的求解时限分开。构建/验证预算与求解时限分别记录；已完成共享Debug检查加Release构建的实测见MEMORY_ESTIMATE.md，不能据此推断单独Debug峰值。

## 已完成结果与资源参考

512²/GMRES(30)已完成18项、9组配对，k=8 DG/HDG为5.2949；该对的运行节点不同。正式数据见[项目结果](../../RESULTS.md#二维512重启30新release)，逐阶RSS/wall及构建成本见[资源参考](../../MEMORY_ESTIMATE.md)。批次状态仅由[状态索引](../STATUS.md)导航，不重复保存报告或数据。
