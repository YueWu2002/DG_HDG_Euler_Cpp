# 2D 正确性

`ctest --test-dir build-2d-intel --output-on-failure` 检查：

- 原始物理 Jacobian、压缩缓存重构与有限差分。
- k=0/1/2/4/8，配点与超积分下的独立求积和张量收缩。
- 完整 DG Jv 的有限差分、块矩阵/算子一致性。
- 代数凝聚、Schur 块/作用及恢复后 DG 真残差恒等式。
- 首次刷新、Jv、装填与凝聚时禁止 Eigen 堆分配。
- 连续块 LU 池的置换/回代顺序，包括 k=8 对应块宽 324 的强制行置换。
- 小网格低阶 BE/DIRK3 的 DG/HDG 最终解与守恒，均使用精确块 Jacobi。

生产路径入口为 `python3 code_test/2d/test_solver.py --exe <exe> --output <独立目录>`，
旋转简单波采用算例内部固定角度，另检查二维熵波；不再提供命令行旋转角输入。构建本身不执行测试。


## Debug验证与Release实验构建

Release以`NDEBUG`及`EIGEN_NO_DEBUG`关闭Eigen运行时断言，包括禁止Eigen堆分配的断言；项目`require(...)`物理/残差检查继续生效。Debug保留尺寸、索引和分配检查，不能据此保证未覆盖输入的正确性。
`intel-build.sh`支持`BUILD_TYPE=Debug|Release`（默认Release）；默认仅Debug开启BUILD_TESTING。code_test驱动要求可追溯到CMakeCache.txt的Debug二进制，拒绝Release或无法核实的构建。
共享`scripts/build_debug_checks.sh <1d|2d> <运行目录>`在`build-debug/`构建Debug，运行CTest及生产路径code_test；实验入口随后在独立`build/`构建Release，供accuracy/speed使用。保留Debug构建日志、CTest、code_test及provenance，也保留Release配置/二进制身份，不将两套二进制混称为同一验证对象。新增构建成本需计入提交前资源预估，不改写在途快照。
