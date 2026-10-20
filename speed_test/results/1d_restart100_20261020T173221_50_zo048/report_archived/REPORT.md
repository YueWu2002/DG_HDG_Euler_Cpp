# 1D N=2048 GMRES(100) 速度报告

18 个计算任务全部 COMPLETED/0:0；汇总、归档成功，完整覆盖 k=0..8 的 DG/HDG。
node1848，Intel Xeon Platinum 8268 @ 2.90GHz，单线程；BE、CFL=5、右 block-Jacobi、natural、ntol=ltol=1e-12、max-linear=1500、HLL safety=1.1。
一个时间步固定四次 Newton 更新，首个预热，仅后三次 iteration_s 算术平均；不以固定更新后的误差宣称非线性收敛。
复用历史 Intel 二进制，本批未重新编译；源码身份和原始编译节点未知，用户已接受本批例外。二进制 SHA256：`da18d8c6e359e544259dd911f8e43ae08f757ed6b5dc5d83a01c5ce54357f8f2`。
每算例独立 Slurm 作业，同节点并发但启动时段和争用不同；仅报告本次测量，不将全部差异归因于 restart。与 N=512 历史结果不计算耗时比。

| k | DG s/更新 | HDG s/更新 | DG/HDG | DG Krylov总数 | HDG Krylov总数 |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.074716 | 0.127305 | 0.5869 | 483 | 477 |
| 1 | 0.270156 | 0.130379 | 2.0721 | 422 | 413 |
| 2 | 0.630409 | 0.138213 | 4.5611 | 469 | 406 |
| 3 | 0.724235 | 0.125587 | 5.7668 | 425 | 427 |
| 4 | 1.072600 | 0.132629 | 8.0872 | 487 | 409 |
| 5 | 1.253979 | 0.144229 | 8.6944 | 433 | 435 |
| 6 | 1.574890 | 0.139067 | 11.3247 | 497 | 411 |
| 7 | 1.675802 | 0.169274 | 9.8999 | 438 | 438 |
| 8 | 1.977131 | 0.125909 | 15.7029 | 502 | 413 |

Krylov总数包含全部四次更新（含预热）。linear_relative_residual 为各自线性系统残差；HDG 恢复后的 dg_linear_relative_residual 单列，不能将两者视为同一残差。
[完整时间/资源](../postprocessing/timings.csv)；[方法比值](../postprocessing/speed_summary.csv)；[原始输入清单](inputs.json)；[本次复核](verification.json)。
原始数据唯一保存在 job_<id>/data/；scheduler/ 保留 out/err。历史命令中的 running 路径按 archive.json 映射到当前 results 路径，不回写原始证据。
