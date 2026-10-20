# 1D N=2048 GMRES(20) 速度报告

18 个计算任务全部 COMPLETED/0:0；汇总、归档成功，完整覆盖 k=0..8 的 DG/HDG。
node1848，Intel Xeon Platinum 8268 @ 2.90GHz，单线程；BE、CFL=5、右 block-Jacobi、natural、ntol=ltol=1e-12、max-linear=1500、HLL safety=1.1。
一个时间步固定四次 Newton 更新，首个预热，仅后三次 iteration_s 算术平均；不以固定更新后的误差宣称非线性收敛。
复用历史 Intel 二进制，本批未重新编译；源码身份和原始编译节点未知，用户已接受本批例外。二进制 SHA256：`da18d8c6e359e544259dd911f8e43ae08f757ed6b5dc5d83a01c5ce54357f8f2`。
每算例独立 Slurm 作业，同节点并发但启动时段和争用不同；仅报告本次测量，不将全部差异归因于 restart。与 N=512 历史结果不计算耗时比。

| k | DG s/更新 | HDG s/更新 | DG/HDG | DG Krylov总数 | HDG Krylov总数 |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.025673 | 0.039709 | 0.6465 | 487 | 482 |
| 1 | 0.113190 | 0.042731 | 2.6489 | 425 | 416 |
| 2 | 0.367583 | 0.046697 | 7.8717 | 487 | 409 |
| 3 | 0.312750 | 0.051437 | 6.0803 | 428 | 432 |
| 4 | 0.625003 | 0.055633 | 11.2344 | 506 | 412 |
| 5 | 0.621798 | 0.064852 | 9.5880 | 436 | 438 |
| 6 | 0.807785 | 0.071251 | 11.3372 | 516 | 412 |
| 7 | 0.861830 | 0.091169 | 9.4531 | 442 | 443 |
| 8 | 1.176961 | 0.092542 | 12.7182 | 523 | 415 |

Krylov总数包含全部四次更新（含预热）。linear_relative_residual 为各自线性系统残差；HDG 恢复后的 dg_linear_relative_residual 单列，不能将两者视为同一残差。
[完整时间/资源](../postprocessing/timings.csv)；[方法比值](../postprocessing/speed_summary.csv)；[原始输入清单](inputs.json)；[本次复核](verification.json)。
原始数据唯一保存在 job_<id>/data/；scheduler/ 保留 out/err。历史命令中的 running 路径按 archive.json 映射到当前 results 路径，不回写原始证据。
