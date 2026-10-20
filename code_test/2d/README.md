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
另检查旋转角 30° 和二维熵波。构建本身不执行测试。
