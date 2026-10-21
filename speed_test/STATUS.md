# 速度测试进度

此页保存尚未整批完成的临时测量与状态，不是正式实验报告。原始档案保留在运行目录；完整批次验收后以唯一正式报告为准。

<a id="二维512重启30新release"></a>
### 部分完成速度：512×512、GMRES(30)，Release

**状态与验收（2026-10-21T19:29:51.466227+00:00）。** 已完成17/18项方法测量，形成8组完整DG/HDG配对；完成作业7687732～7687748均为COMPLETED/0:0，逐项SUCCESS、进程退出码、四次更新及0/1/1/1预热标记已核对，后三次均值已重算。k=8 DG作业7718103已在node2318运行，尚无完整测量，不填零、不计算该阶速度比；汇总7718104和归档7718105等待依赖，本批不标为整体完成。

**环境与构建。** 已完成17项均在node2347，Intel(R) Xeon(R) Platinum 8268 CPU @ 2.90GHz，每任务1核、单线程；各方法独立作业，可在同节点并发。共享构建7687730成功，Debug的CTest为1/1通过，生产路径BE、DIRK3、entropy检查通过；独立Release使用IntelLLVM/icpx 2025.2.0，`-O3 -DNDEBUG`，Eigen断言关闭。Debug检查与Release计时为不同二进制，不混称为同一验证对象。

**协议。** BE、CFL=5、512×512、k=0～8，单步固定四次Newton更新；首个预热，仅后三次`iteration_s`取平均。精确右block-Jacobi、natural、ntol=ltol=1e-12、restart=30、max-linear=1500、HLL safety=1.1、timing-only。单方法进程wall和Slurm总时长均不替代平均更新时间。

| k | DG s/更新 | HDG s/更新 | DG/HDG | DG Krylov总数 | HDG Krylov总数 |
|---:|---:|---:|---:|---:|---:|
| 0 | [23.457089](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687747/data/k0_dg_summary.csv) | [30.489270](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687748/data/k0_hdg_summary.csv) | 0.7694 | 442 | 489 |
| 1 | [74.986093](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687745/data/k1_dg_summary.csv) | [86.743751](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687746/data/k1_hdg_summary.csv) | 0.8645 | 422 | 449 |
| 2 | [241.207543](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687743/data/k2_dg_summary.csv) | [251.747531](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687744/data/k2_hdg_summary.csv) | 0.9581 | 570 | 611 |
| 3 | [520.504391](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687741/data/k3_dg_summary.csv) | [421.557828](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687742/data/k3_hdg_summary.csv) | 1.2347 | 602 | 566 |
| 4 | [1414.225684](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687739/data/k4_dg_summary.csv) | [448.307715](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687740/data/k4_hdg_summary.csv) | 3.1546 | 719 | 647 |
| 5 | [3170.293415](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687737/data/k5_dg_summary.csv) | [927.543823](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687738/data/k5_hdg_summary.csv) | 3.4179 | 810 | 686 |
| 6 | [4784.598631](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687735/data/k6_dg_summary.csv) | [1146.095157](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687736/data/k6_hdg_summary.csv) | 4.1747 | 902 | 713 |
| 7 | [6505.288222](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687733/data/k7_dg_summary.csv) | [1900.703343](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687734/data/k7_hdg_summary.csv) | 3.4226 | 989 | 809 |
| 8 | —（未完成） | [2336.131207](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687732/data/k8_hdg_summary.csv) | — | — | 796 |

**测量结论。** k=0～2时DG更快；k=3起已完成配对中HDG更快。k=4～7的DG/HDG依次为3.1546、3.4179、4.1747、3.4226；k=8 HDG为2336.131207 s/更新，缺DG配对。Krylov总数包含四次更新，时间平均仅覆盖后三次。已完成进程峰值RSS最高约154.960 GiB（k=7 DG）；k=8 HDG为153.556 GiB，均为进程实测而非申请值。

**残差口径。** 17项各自线性系统的最大相对残差为9.98548e-13；DG为逆体积质量加权残差，HDG为逆面质量加权Schur残差。HDG回代后的DG相对残差最高7.52330e-12，不能据此标称两方法均满足1e-12的DG相对残差。timing-only误差范数字段nan为预期，固定四次更新不用于精度或非线性收敛结论。

#### 与256×256、GMRES(20)历史测量对照

旧批采用同型号8268 CPU，在node2301顺序运行；新批在node2347分方法并发。网格单元数增加4倍，restart从20变30，代码/二进制及Release配置身份也不同。下表为实测描述，不把耗时变化全部归因于restart或网格，也不据此宣称某项代码优化收益。耗时倍率为同方法的512²/G30平均更新时间除以256²/G20平均更新时间。

| k | 256²/G20 DG/HDG | 512²/G30 DG/HDG | DG耗时倍率 | HDG耗时倍率 |
|---:|---:|---:|---:|---:|
| 0 | 1.0474 | 0.7694 | 4.0071 | 5.4553 |
| 1 | 0.9005 | 0.8645 | 5.3916 | 5.6161 |
| 2 | 1.1073 | 0.9581 | 5.8092 | 6.7138 |
| 3 | 1.5823 | 1.2347 | 5.9385 | 7.6101 |
| 4 | 2.2169 | 3.1546 | 6.9072 | 4.8541 |
| 5 | 2.5633 | 3.4179 | 8.5518 | 6.4136 |
| 6 | 3.4049 | 4.1747 | 6.9619 | 5.6782 |
| 7 | 3.7243 | 3.4226 | 5.5083 | 5.9939 |
| 8 | 5.3598 | — | — | 5.8337 |

**对照结论。** k=4～6的HDG相对优势比旧批增大，k=7略减小；k=2由HDG略快变为DG略快。当前17项的512²/G30相对256²/G20耗时倍率为4.0071～8.5518，不能按单元数4倍简单预测。k=5 DG为8.5518倍；k=8 HDG为5.8337倍。

**版本与来源。** 新批版本为`ba0b3fcb4135+src.fc044fefaca984cb`；完整源码清单SHA256为`fc044fefaca984cbe489da8f0be56afc5ebd191791a35464565b724e6fc85b37`，Release可执行文件SHA256为`5429a2f3a7556476bfb88ebd9241bc0b1abe70315af34d26f7a1f03f42715e79`。旧批速度版本标识为`5383fb8b2864+src.2fe268224441020d`，完整源码清单SHA256为`2fe268224441020d0ddb5dbe685ea03f61e0d083e86ecfb30061002e94c7a665`，二进制SHA256为`812e1075efc2854e89cfc02de8426e23ed757eb1abd497fc104ce2771e98b533`；不同哈希范围不能仅凭摘要判断实现差异。

新批：[配置与资源估算](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/submission.json)、[源码清单](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/source.sha256)、[Release身份](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/build_identity.json)、[构建配置](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687730/build/CMakeCache.txt)、[CTest日志](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687730/ctest.log)、[code_test日志](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/job_7687730/code_test.log)、[任务清单](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/submitted_jobs.csv)。表中时间链接指向原始summary；同目录保留逐更新CSV、完整命令、CPU/环境、二进制身份及资源文件。

k=8 DG实际运行主机为node2318.oscar.ccv.brown.edu，已核对CPU为Intel(R) Xeon(R) Platinum 8268 CPU @ 2.90GHz，执行文件SHA256与本批Release一致。调度使用`48core&intel&cascade&edr`，不固定节点，运行时仍精确核验8268型号；其最终实际节点、CPU和二进制身份以启动记录为准。来源：[调度条件变更](../speed_test/running/2d_n512_restart30_20261020T220633_rt2c3wx6/node_policy_change_20261021/change.json)。

旧批：[原始计时和资源](../speed_test/results/2d_submission_20261019T205542_KhFC6z/job_7635952/data/timings.csv)、[历史正式报告](../speed_test/results/2d_submission_20261019T205542_KhFC6z/report_archived_20261020/REPORT.md)。本节仅整理结果并引用原件；在途目录保持原位，不复制原始数据或创建第二份批次报告，待整个批次完成并通过验收后，再由既有汇总/归档流程生成正式报告并更新结果页。
