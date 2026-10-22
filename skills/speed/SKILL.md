---
name: speed
description: 运行或分析本项目 Eigen DG/HDG 1D/2D 速度与内存实验，核对四次 Newton 更新的计时、CPU、Slurm 状态和版本来源。
---

先读根目录 `AGENTS.md`。按需读 `speed_test/{1d,2d}/benchmark.py`、对应 `submit.sh`、`scripts/submit_experiment.sh`、`scripts/submit_2d.sh`、`scripts/experiment{,_2d}_job.sh`、`scripts/experiment.py` 与 `MEMORY_ESTIMATE.md`。

## 正式工作负载

- 1D：512 单元，k=0～8；2D：256×256，k=8、7、6、…、0。
- BE、CFL=5、一个时间步内固定四次 Newton 更新；第一次为预热，只取后三次 `iteration_s` 的算术平均。两方法都按该协议，不以整个进程 wall time 或 `time_loop_s/4` 替代。
- 精确 block-Jacobi、固定右预条件 GMRES、自然排序、线性与 Newton 容差 1e-12、restart=20、max-linear=1500、HLL safety=1.1。2D alpha=45°。读取脚本确认完整参数，用户的新要求优先。
- 使用 `intel-build.sh` 的 Intel/MKL/Eigen Release、native 等实际配置，单线程。正式速度任务提交 Slurm；本地运行只用于诊断。默认配对比较可在同一资源分配中顺序执行；用户要求拆任务时按独立方法提交，保持节点要求并记录并发条件、各自 Krylov 数和全部计时边界。
- 1D k=8 的目标为 DG/HDG 至少 8 倍；这是目标，不能替代测量，也不能外推到 2D 或低阶。

## 提交前资源预估（必须执行）

- 计算时限必须先根据匹配历史数据和实际规模预估，再充分增大。Krylov迭代次数、restart、谱性质、网格加密及节点争用会改变成本，不能准确预测总耗时；不得只按单元数线性缩放或只加很小余量。根据不确定性明确选择充分的时间倍率，记录原始估计、放大依据、最终时限；用户指定双倍时限时，逐项对上一轮申请时限翻倍。超时数据只能作为耗时下界，不能作为成功运行需求。
- 新构建可执行文件统一位于各构建目录的 `bin/*.exe`，静态库为构建目录内的 `*.a`；批次运行副本位于该批次 `bin/solver.exe`。Debug与Release各自保留配置、日志及二进制身份。历史命令与旧路径不回写。


- 优先读取根目录 [MEMORY_ESTIMATE.md](../../MEMORY_ESTIMATE.md) 的匹配实测表与估算说明，作为资源预估依据；其中已记录规模、协议、实际RSS/wall和来源，覆盖充分时不必翻找或重新遍历历史数据，以减少重复读取和token开销。仅在缺少对应规模/配置、版本不匹配、需要更高精度或数据存疑时，沿该文档来源链接定向核查必要文件。新增有效测量后及时更新该参考表，未知信息仍明确标注。
- 提交任何新 Slurm 测试前，必须先预估该任务的峰值内存和总运行时长，再据此确定 `--mem`、`--time`；不得先提交统一的大额资源申请、事后才估算，也不得将最大算例申请机械套用到所有任务。
- 优先查找已完成且可比的历史运行，核对维数、单元总数、阶数、DG/HDG、积分器、步数或终止时间、容差、GMRES restart、预条件器、代码/二进制版本、CPU和线程。以实际 `peak_rss_kib`、进程 `wall_s` 及对应 Slurm accounting 为依据；申请量不作为实测值，失败/OOM或超时数据不能当作成功运行的完整需求。
- 内存按网格和自由度变化、局部精确块存储、Krylov基及辅助工作区缩放；时长按单步成本、实际步数、Newton/Krylov迭代变化及硬件差异估算。网格加密不必然使耗时仅按单元数线性增长。没有可比历史记录时，先用存储公式和已有小规模数据估算；必要时提交有完整provenance的小规模探测，明确不确定性，不在登录节点运行正式计算。
- 每个独立任务分别计算申请值，在估算峰值/总时长上增加明确且合理的余量，再按Slurm单位向上取整。总时长须涵盖该作业内的构建、检查、全部求解进程及辅助步骤；共享构建和后处理独立作业分别估算。顺序进程内存按同时驻留的峰值计算，并发进程按重叠工作集计算。余量或倍率按本次用户要求和证据确定，不把某次实验的4倍、5倍或固定时限当作通用默认。
- 在提交元数据中保存历史来源路径/任务号及版本、原始RSS和耗时、规模与协议差异、缩放公式/倍率、预计峰值内存/总时长、余量及最终申请值；提交前核对节点容量、账户配额和并发总申请。运行后分别记录各求解器进程的 `peak_rss_kib`、`wall_s`、退出码及Slurm资源统计，以更新后续估算；修改已排队或运行作业的资源须有用户授权，不覆盖原始证据。

- speed：按每个 `(维数,网格,k,method,restart)` 的实际工作负载估算。申请时长使用完整进程wall及作业额外开销，不能将仅后三次平均的 `mean_newton_s` 直接当成总时长；四次更新中含一次预热。

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

- 归档时遵循根目录 `AGENTS.md` 的“完整测试档案的保留”：将 Slurm `.out`、`.err` 和全部计算输出、后处理、日志、元数据及运行辅助文件一并保存在 `speed_test/results/<对应测试目录>/`；只允许在任务及依赖全部结束后删除特别大的可执行文件，并保留删除记录及二进制身份。不得仅归档 CSV 或清空构建目录。

- accuracy_test 与 speed_test 均在顶层使用 `running/` 和 `results/`：即 `accuracy_test/{running,results}/`、`speed_test/{running,results}/`。`1d/`、`2d/` 继续保存测试驱动和说明，不再分散保存新结果。
- 每次运行使用独立文件夹，名称带维数、时间戳和唯一后缀，例如 `running/2d_submission_<时间戳>_<唯一后缀>/`。排队和运行中的任务都放在 `running/`。
- 每次运行的文件夹必须完整包含其计算输出、后处理数据、CPU 型号、实际参数、代码版本及源码哈希、构建环境、任务状态等；保留现有 `job_<id>/data/` 输出布局，不额外新设计算输出目录。报告放在对应运行文件夹内，跨运行报告必须记录输入运行的身份。
- 任务结束后整体移到 `results/`，保持文件夹内部结构。成功、失败、取消均应归档并记录状态；不得把部分结果标为成功，不得覆盖已有运行目录。使用 `python3 scripts/archive_finished.py` 核对 Slurm 并归档；该命令应在任务结束后执行，而不是从运行中的批处理脚本内搬移自身目录。
- 已经运行或排队的旧任务属于过渡例外：保留原有路径，禁止搬移或修改其文件、源码快照和日志。待任务完全结束并确认没有其他任务仍依赖其路径后，再整理到顶层 `results/`。
- 历史数据按已知运行批次整理到 `results/<运行标识>/`；保留所有可用 CSV 和设备日志，无法恢复的 CPU、代码版本或原始输出必须明确标记为未知或缺失，不得伪造。不同独立运行不得因日期相同而混为一次运行。
- `.gitignore` 必须允许 CSV、设备/环境日志、命令、版本清单和追溯元数据纳入 Git。只忽略构建产物和重复源码快照等；Git 忽略规则不代表允许删除磁盘上的原始运行数据。
- 新增本地测试、Slurm 提交、后处理或报告入口时，必须遵循上述目录和归档约定，并更新相关文档引用。既有在运行任务使用的旧入口快照保持不动。
- speed 的 1D/2D 本地测试入口省略 `--output` 时，必须自动创建 `running/<维数>_local_<时间戳>_<唯一后缀>/`；显式输出路径也必须在对应 `running/` 内。Slurm 测试输出位于 `running/<提交目录>/job_<id>/data/`。不得直接在 `results/` 内启动计算。

源码副本只作为任务排队、构建和运行期间的临时输入。新提交通过 `retain_source_after_completion=false` 标记，确认任务与依赖结束后由归档工具移除 `source/`；保留 `source.sha256`、版本、配置及命令。不得提前删除在途源码。

## 统一的2D后处理流程

与另一类实验共用 `scripts/postprocess.py`，通过本技能对应的 `speed_test/2d/postprocess.py --run-root <实验目录>` 调用。参数 `--run-root` 指向包含 job_id.txt 的实验目录，不指向 job 或 data。本地执行后处理不启动求解器、不提交新计算。

- 原始文件唯一保存在 `job_<id>/data/`。后处理 CSV 放在同级 `postprocessing/`，正式报告及输入清单放在唯一的 `report_archived_<日期>/`。反复执行复用已有报告目录，不因归档或日期变化复制报告与数据。
- accuracy 生成 `convergence_summary.csv`；speed 生成 `speed_summary.csv`。两者都保留参数、CPU、主机、任务号、源码与二进制身份，原始数据不改动。旧任务缺少统一版本号时，根据保留的基准提交和源码清单哈希计算来源标识，并在 inputs.json 中记录实际依据，不推断丢失内容。
- 先检查 SUCCESS、数据完整性、配对、数值协议和身份，再重算科学指标。accuracy 核验各网格 summary 和差值文件，重算实际网格比的收敛阶；speed 核验原始 summary、四次更新及0/1/1/1标记、后三次平均、资源记录与配对比值。不同来源不兼容时拒绝汇总，不用缺失值冒充零。
- `inputs.json` 使用相对实验目录的输入路径及 SHA256，记录后处理器哈希、输出哈希、实际网格/阶数及 full_formal_set。部分自定义扫参可以生成报告，但不能宣称完整正式覆盖。报告的 SUCCESS 仅表示所声明覆盖通过验收。
- 报告通过相对链接引用原始与后处理文件，不 copytree 原始目录。采用新结构后只维护当前 REPORT.md 和 inputs.json，不创建 history/ 或 before_unification_* 备份；必要的独有来源信息合并到现行元数据。原始实验数据仍完整保留，删除重复数据文件必须先验证与保留的原件逐文件哈希一致。
- 报告明确说明实际残差口径、时间单位、预热与统计方式。检查失败时不覆盖已有正式结果；不能通过后处理修改求解器输出。
- `scripts/collect_2d.py` 仅为兼容入口，分别调用两类后处理；不再生成第三份联合数据档案。`--output` 已废弃，改用各实验的标准目录。需要 Slurm 执行时 `scripts/submit_2d_report.sh` 保存后处理任务自身日志，计算结果仍写入两次实验各自的标准后处理目录。

数值算例配置与求解器参数分开：后续结果CSV、报告和后处理参数表不记录 `alpha_degrees`、`gamma`，不把它们作为求解器配置列或配对协议字段。保持实际算例计算不变；目前暂不增加案例identifier。历史原始数据不回写。

旋转角属于算例内部配置，求解器不接受 `--alpha-deg`（或 `alpha_deg`）命令行输入；新测试命令不传该参数，直接使用 Example 内的固定配置。旧的 command.json 是历史执行证据，不修改。

## 按独立算例拆分 Slurm 与自动后处理

- 用户要求并行扫参时，每个独立计算任务单独提交：速度按 `(维数,k,method)` 拆分，精度按单网格拆分（可保留网格内DG/HDG配对）。计算任务之间不设依赖；只在用户要求复现指定机器时设置 nodelist，精度默认无需相同节点。
- 按历史单进程峰值、当前工作集和新增Krylov基向量估算内存并留余量。更改restart时增加 `8*N*(restart_new-restart_old)` 字节及小工作区；本次用户要求每个任务至少16 GiB，计算、汇总、归档均遵守。
- 全部计算任务之后提交轻量 `afterany` 汇总任务，核对每个进程退出状态、原始CSV、四次更新及预热规则；失败项明确保留，不补零，不将失败批次标记SUCCESS。汇总CSV和报告使用统一的 `postprocessing/`、唯一 `report_archived*` 布局，不复制原始数据。
- 再以 `afterany` 依赖汇总任务提交归档任务；归档在计算与汇总均结束后，从实验目录外执行，将整个运行目录从 running 移到 results。归档保留所有scheduler out/err、命令、CPU、参数、版本、资源与失败记录；不要从仍在计算的任务中搬移其工作目录。
- restart对照实验优先复用与基线SHA256一致的二进制，保持物理、容差、预条件、正交化和计时协议不变；旧源码身份无法恢复时明确未知，记录二进制身份，不能把当前源码版本冒充基线版本。
- 分任务同节点并发会改变节点争用条件，即使求解器只改restart，也要在报告中注明；相同node并不保证与旧顺序运行的系统负载相同。

restart=100 的独立任务入口为 `scripts/speed_campaign.py`：先 `prepare --dimension 1d|2d --restart 100 --minimum-memory-gib 16` 生成不可变驱动与基线二进制副本，再 `submit --run-root <生成的目录>`。提交阶段只做文件准备和sbatch，不在登录节点执行求解器。每维18个计算任务、1个汇总任务、1个归档任务；submitted_jobs.csv与job_ids.json记录完整依赖链。该入口保留本次原节点和历史内存依据，不将其节点硬编码规则推广到其他实验。

collect与archive任务不限制node、CPU型号或Intel特性，可在任意可用计算节点执行；原节点/指定CPU的约束仅用于实际数值计算任务。后处理仍记录实际设备信息，不对设备型号作验收限制。

1D 自定义网格与同节点 restart 对照：`prepare --dimension 1d --restart 20|100 --cells 2048 --node <已核实节点> --minimum-memory-gib 16`。两套均覆盖 k=0..8 的 DG/HDG，复用同一基线二进制；1D 历史源码身份未知，以二进制完整 SHA256 为准。内存估算按单元数放大历史峰值并增加 Krylov 基存储。网格与历史基线不同时，不计算历史耗时比；两套新结果可按匹配参数比较。


## 新源码二维速度独立作业入口

`python3 scripts/prepare_snapshot_speed.py --node <探测节点> --cells <每方向单元数> --restart <GMRES重启长度> --factor <资源倍率> --time-multiplier <时间放大倍率>`，随后对生成目录执行 `python3 scripts/speed_campaign.py submit --run-root <目录>`。
入口优先用MEMORY_ESTIMATE.md的二维256²/GMRES(20)实测行估算，保存引用片段和文档哈希；内存先按单元数缩放，再加`8*全局线性自由度*(restart-20)`的新增基存储，最后加余量；时间先按历史完整进程wall乘资源倍率并向上取整，再乘独立的 `--time-multiplier`（默认2），为Krylov迭代变化充分放大；实际倍率与最终时限写入提交元数据。restart变化的迭代/正交化成本仍标为不确定。倍率需覆盖单元数比例，不能当作任意新协议的通用时长保证。
所有(k,method)独立提交到同一指定节点，共享构建成功后各自调度；共享构建保留必要依赖，汇总/归档亦保留依赖。源码快照在该计算节点编译并检查，记录实际CPU及二进制身份；结束后保留全部档案并移除临时source/。可选节点探测作业记录在同批任务清单中。


## Debug验证与Release实验构建

Release以`NDEBUG`及`EIGEN_NO_DEBUG`关闭Eigen运行时断言，包括禁止Eigen堆分配的断言；项目`require(...)`物理/残差检查继续生效。Debug保留尺寸、索引和分配检查，不能据此保证未覆盖输入的正确性。
`intel-build.sh`支持`BUILD_TYPE=Debug|Release`（默认Release）；默认仅Debug开启BUILD_TESTING。code_test驱动要求可追溯到CMakeCache.txt的Debug二进制，拒绝Release或无法核实的构建。
共享`scripts/build_debug_checks.sh <1d|2d> <运行目录>`在`build-debug/`构建Debug，运行CTest及生产路径code_test；实验入口随后在独立`build/`构建Release，供accuracy/speed使用。保留Debug构建日志、CTest、code_test及provenance，也保留Release配置/二进制身份，不将两套二进制混称为同一验证对象。新增构建成本需计入提交前资源预估，不改写在途快照。

新快照入口支持`--build-time HH:MM:SS`，用于共享Debug验证加Release构建作业，单独记录预算；与每个算例的求解时限分开。构建/验证时限与各算例时限分别记录；512²/GMRES(30)批次已有共享Debug检查加Release构建实测，见MEMORY_ESTIMATE.md。该数据不代表单独Debug峰值，后续预算保留编译与负载余量。

共享脚本的目录约定见 [scripts/README.md](../../scripts/README.md)。本地驱动写入进程身份与退出状态；显式使用 `scripts/archive_finished.py --run-root <目录>` 归档已停止的本地运行。源码快照不打包bin、构建目录及二进制；Git只忽略产物，保留配置、日志和完整来源。

指定CPU型号时，新快照入口使用 `--required-cpu "<精确CPU型号>"`；构建开始后先核验实际CPU与节点，再执行Debug和Release构建。`build_identity.json`保存Release构建相对目录，使`bin/solver.exe`运行副本在归档前后均可追溯CMake配置。

已完成512²/GMRES(30)批次的18项RSS/wall和共享构建成本已加入[资源参考](../../MEMORY_ESTIMATE.md)。后续匹配规模优先查该表；k=8 DG/HDG跨节点，速度比保留节点与并发负载限制。批次结束后将STATUS.md收敛为正式报告索引，并更新RESULTS.md、README.md和资源参考，避免保留过期running路径。
