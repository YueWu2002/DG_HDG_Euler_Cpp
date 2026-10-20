# 1D / 2D 内存估算

令 C 为单元数、n=k+1、m=3n、t=6、R 为 GMRES restart。
DG 状态和质量为 O(Cm)，通量导数 O(9Cn)，
block-Jacobi 对角块和 LU 为 O(Cm²)，Krylov 基 O(CmR)。
没有 O((Cm)²) 全局矩阵，`--linear direct` 的小规模参考例外。

HDG 体积局部凝聚保存 A/LU、B/C、X、y 与 t×t Schur 块，
存储 O(C(m²+mt+t²))；全局 trace 大小为 3C，
Krylov 存储 O(3CR)。缓存局部 Schur 增加 36C 个 double，
512 单元约 144 KiB，避免每次乘法重复体积维度收缩。
实际 RSS 还包括 Eigen、动态库、工作区、计时记录和运行环境，须以实验资源数据为准。
二维高阶实验的独立估算见下节。

## 2D，256×256，k=8，100 GiB 预算

n=9、m=324、单面宽 fw=36、局部面宽 t=144、C=65536。
主数值池（GiB，未含库与分配器开销）：

| 数值池 | GiB |
|---|---:|
| DG 精确 LU（只存一份） | 51.26 |
| HDG X | 22.78 |
| HDG 局部 Schur | 10.13 |
| HDG 面块 Jacobi LU | 1.27 |
| DG GMRES 基，restart=20 | 3.32 |
| HDG GMRES 基，restart=20 | 0.74 |

此外为状态、阶段斜率、质量/逆质量、Krylov 辅助向量、置换和少量局部工作区。
DG/HDG 使用独立进程顺序运行，两个方法的主要池不同时驻留。
估计精确 DG 总量约 60 GiB、HDG 约 40 GiB；保留 100 GiB 作业预算，记录实际峰值 RSS。
`python3 speed_test/2d/memory_estimate.py` 可重算主数值池。
2D 的精确局部凝聚只复用一套 m×m A/LU，不按 C 复制。
