---
name: speed
description: 运行或分析本项目 Eigen DG/HDG 1D/2D 速度与内存实验，核对四次 Newton 更新的计时、CPU、Slurm 状态和版本来源。
---

先读根目录 `AGENTS.md`。按需读 `speed_test/{1d,2d}/benchmark.py`、对应 `submit.sh`、`scripts/submit_experiment.sh`、`scripts/submit_2d.sh`、`scripts/experiment{,_2d}_job.sh`、`scripts/experiment.py` 与 `MEMORY_ESTIMATE.md`。

## 正式工作负载

- 1D：512 单元，k=0～8；2D：256×256，k=8、7、6、…、0。
- BE、CFL=5、一个时间步内固定四次 Newton 更新；第一次为预热，只取后三次 `iteration_s` 的算术平均。两方法都按该协议，不以整个进程 wall time 或 `time_loop_s/4` 替代。
- 精确 block-Jacobi、固定右预条件 GMRES、自然排序、线性与 Newton 容差 1e-12、restart=20、max-linear=1500、HLL safety=1.1。2D alpha=45°。读取脚本确认完整参数，用户的新要求优先。
- 使用 `intel-build.sh` 的 Intel/MKL/Eigen Release、native 等实际配置，单线程。正式速度任务提交 Slurm；本地运行只用于诊断。DG/HDG 在同一次资源分配中依次执行，记录各自 Krylov 数和全部计时边界。
- 1D k=8 的目标为 DG/HDG 至少 8 倍；这是目标，不能替代测量，也不能外推到 2D 或低阶。

## 在途任务与存储

通过每次运行的 `job_id.txt`、`squeue`、`scontrol`、`sacct` 和日志交叉确认。登录节点 `ps` 不能证明远程任务停止；无法确认时保护其依赖。

新任务放入 `speed_test/running/<独立运行目录>/`，包含源码快照、构建、环境、参数、版本及 `job_<id>/data/`。本地 benchmark 省略 `--output` 创建唯一 running 目录。旧的在途目录保持原位。独立优化构建不得覆盖在途任务的二进制或脚本；状态查询不重新提交、取消或重跑任务。

2D 256²/Q8 精确 block-Jacobi 的申请预算可从 100 GiB 起核对；这不是实测 RSS。使用 `python3 speed_test/2d/memory_estimate.py --cells 256 --k 8`，再计入未列出的工作区、库、阶段状态和构建开销。其他规模按工作集估算，不机械复制申请量。

## 读取与验收

- 每方法的 `<prefix>_iterations.csv` 恰有四次更新，`measured` 为 0/1/1/1，summary 的 measured_updates=3；重算后三次平均并与 mean_newton_s 对照。固定更新的最终误差不用于收敛结论。
- `timings.csv` 是方法计时和资源汇总；`speedup.csv` 是 DG/HDG 均值比。完整 1D/2D 扫描分别有18条方法记录、9条加速比。未完成配对不得形成比值。
- `*.resources.json` 记录单个求解器进程 wall/user/system 时间、peak_rss_kib 和退出码；`environment.txt`、provenance 与命令记录提供 CPU、参数和代码身份。运行中可用 sstat，结束后用 sacct；申请内存不是实际峰值。
- 空 CSV、缓冲日志或很久没有新行不能独立证明失败或停滞。说明最新可见阶段，不推断已完成的 Newton 次数。
- 任务结束后核对状态再运行 `python3 scripts/archive_finished.py`。保留失败数据，结果放在 `speed_test/results/<同名目录>/`。报告和后处理留在所属运行目录，跨运行报告记录每份输入的身份。

历史 `speed_test/results/1d_slurm_7622721/` 同时收纳 profile CSV；profile 的独立来源记录是 `profile_provenance.json`。不得将该目录里的 Slurm CPU/二进制身份自动套用给 profile，或声称两者已确认同版本。

只对匹配机器、线程、编译和数值工作负载的完成测量宣称性能收益。参数或硬件不同则分别展示。保留原始数据及设备日志；结果文档可记测量值，不能把历史值写成当前运行结论。

## 运行目录与归档开发规范

- accuracy_test 与 speed_test 均在顶层使用 `running/` 和 `results/`：即 `accuracy_test/{running,results}/`、`speed_test/{running,results}/`。`1d/`、`2d/` 继续保存测试驱动和说明，不再分散保存新结果。
- 每次运行使用独立文件夹，名称带维数、时间戳和唯一后缀，例如 `running/2d_submission_<时间戳>_<唯一后缀>/`。排队和运行中的任务都放在 `running/`。
- 每次运行的文件夹必须完整包含其计算输出、后处理数据、CPU 型号、实际参数、代码版本及源码哈希、构建环境、任务状态等；保留现有 `job_<id>/data/` 输出布局，不额外新设计算输出目录。报告放在对应运行文件夹内，跨运行报告必须记录输入运行的身份。
- 任务结束后整体移到 `results/`，保持文件夹内部结构。成功、失败、取消均应归档并记录状态；不得把部分结果标为成功，不得覆盖已有运行目录。使用 `python3 scripts/archive_finished.py` 核对 Slurm 并归档；该命令应在任务结束后执行，而不是从运行中的批处理脚本内搬移自身目录。
- 已经运行或排队的旧任务属于过渡例外：保留原有路径，禁止搬移或修改其文件、源码快照和日志。待任务完全结束并确认没有其他任务仍依赖其路径后，再整理到顶层 `results/`。
- 历史数据按已知运行批次整理到 `results/<运行标识>/`；保留所有可用 CSV 和设备日志，无法恢复的 CPU、代码版本或原始输出必须明确标记为未知或缺失，不得伪造。不同独立运行不得因日期相同而混为一次运行。
- `.gitignore` 必须允许 CSV、设备/环境日志、命令、版本清单和追溯元数据纳入 Git。只忽略构建产物和重复源码快照等；Git 忽略规则不代表允许删除磁盘上的原始运行数据。
- 新增本地测试、Slurm 提交、后处理或报告入口时，必须遵循上述目录和归档约定，并更新相关文档引用。既有在运行任务使用的旧入口快照保持不动。
- speed 的 1D/2D 本地测试入口省略 `--output` 时，必须自动创建 `running/<维数>_local_<时间戳>_<唯一后缀>/`；显式输出路径也必须在对应 `running/` 内。Slurm 测试输出位于 `running/<提交目录>/job_<id>/data/`。不得直接在 `results/` 内启动计算。
