---
name: debug
description: 定位和修复本项目 Eigen Euler DG/HDG 1D/2D 的构建、残差、Newton/GMRES、张量索引和更新一致性问题，并运行针对性验证。
---

先读根目录 `AGENTS.md`，按问题查 `DESIGN.md`、`MATHEMATICS.md` 和 `VALIDATION.md`。保留 Eigen、现有 common/dg_hdg_1d/dg_hdg_2d 架构，不引入 MFEM 依赖。对照实现位于 `../mfem_euler2d_cpp`。

## 源码定位

| 问题 | 路径 |
|---|---|
| CLI、输出与统计 | `dg_hdg_{1d,2d}/src/config.cpp`、`main.cpp` |
| DIRK/BE、CFL、Newton、线搜索和残差验收 | `dg_hdg_{1d,2d}/src/time_integrator.cpp` |
| DG 残差及作用 | `dg_hdg_{1d,2d}/src/discretization.cpp` |
| HDG 凝聚、恢复、迹系统 | `dg_hdg_{1d,2d}/src/hybrid.cpp` |
| 2D 体积、面与张量布局 | `dg_hdg_2d/src/{volume_jacobian,face_jacobian,tensor_product}.cpp` |
| 2D 初始化与误差 | `dg_hdg_2d/src/{example,initial_and_error,physics}.cpp` |
| 线性求解与局部因子 | `common/src/linalg/gmres.cpp`、`common/include/euler/linalg/{gmres,block_jacobi,dense_lu}.hpp` |

先取得实际失败命令、源码和二进制哈希、设备、日志、退出状态。只读诊断不顺便启动正式扫描；实现修复时使用独立构建与输出目录，不修改在途快照。

## 数值不变量

GMRES 为自实现的重启右预条件求解器，内积采用逆质量权重。查看构造处的 DG/迹权重与最终物理 DG 残差验收；不能把 Euclidean norm、Arnoldi 估计或迹残差单独当作共同 DG 残差。HDG 恢复更新必须满足原 DG 线性化方程。零右端提前返回也要保持残差缓存正确。

核对 HLL 波速冻结、阶段 dt、Newton 状态和缓存失效时机。每个隐式阶段冻结波速，不能在 Newton/线搜索内擅自重估而改变方程。检查质量逆缩放、体积/面权重、周期连接、张量展平与变量顺序、边界面符号、pivot 和 view 生命周期。

不要用放宽容差、减少求积、近似局部逆或取消独立残差检查掩盖错误。初始化不满足物理可容许性时先检查投影、网格和算例，不直接断言迭代失败。CSV 末行可能未写完整；更新前残差不是最终收敛残差。

## 构建与验证

用 `intel-build.sh`，显式设置独立 BUILD_DIR、BUILD_1D、BUILD_2D 和 BUILD_JOBS；目录名含 debug 不代表调试编译。保持 BLAS/OpenMP 单线程，先运行相关 CTest，再运行对应 `code_test/<维数>/test_solver.py --exe <本次二进制> --output <独立检查目录>`。

- 1D 常规验证：DIRK smoke 与完全收敛的单步 BE，使用 `--check-update` 对照 DG/HDG 更新及最终解。
- 2D 常规验证：沿用现有驱动的 dirk3、be、entropy 三项（旋转角由算例内部固定，不提供命令行角度输入），并检查解差与守恒；不直接照搬 MFEM 的文件字段或两项测试假设。
- 需要严格单步 BE 诊断时明确使用 `--integrator be --steps 1`，按收敛准则停止，不设置测速的固定 Newton 次数。

`*_difference.csv` 是方法最终解差；更新检查验证的是线性化更新，两者不可互替。高阶索引、缓存、恢复等改动增加必要的定向检查即可，不扩成长耗时的日常扫描。精度归 accuracy_test，性能归 speed_test，code_test 保持短小验证。

报告修复、依据、实际运行的检查与限制；明确区分当前源码构建的二进制和旧二进制。缺失的历史元信息保留为未知。


## Debug验证与Release实验构建

Release以`NDEBUG`及`EIGEN_NO_DEBUG`关闭Eigen运行时断言，包括禁止Eigen堆分配的断言；项目`require(...)`物理/残差检查继续生效。Debug保留尺寸、索引和分配检查，不能据此保证未覆盖输入的正确性。
`intel-build.sh`支持`BUILD_TYPE=Debug|Release`（默认Release）；默认仅Debug开启BUILD_TESTING。code_test驱动要求可追溯到CMakeCache.txt的Debug二进制，拒绝Release或无法核实的构建。
共享`scripts/build_debug_checks.sh <1d|2d> <运行目录>`在`build-debug/`构建Debug，运行CTest及生产路径code_test；实验入口随后在独立`build/`构建Release，供accuracy/speed使用。保留Debug构建日志、CTest、code_test及provenance，也保留Release配置/二进制身份，不将两套二进制混称为同一验证对象。新增构建成本需计入提交前资源预估，不改写在途快照。
