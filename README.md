# Eigen Euler：1D / 2D DG 与凝聚 HDG

保留 `common/`、`dg_hdg_1d/`、`dg_hdg_2d/`、`cmake/`、`vendor/fastgl/` 架构，
以现有 1D 算法为核心，使用 Eigen 和 FastGL，不依赖 MFEM。
1D 与 2D 分别构建独立可执行文件和测试；公共线性代数复用 Eigen。
数值协议参照 [mfem_euler2d_cpp](https://github.com/YueWu2002/mfem_euler2d_cpp)，
将旋转二维简单波限制到一维周期区间 [-1,1]，gamma=3。

- [数学定义](MATHEMATICS.md)：简单波、Gauss DG、DIRK3、凝聚及残差标准。
- [架构与性能边界](DESIGN.md)：Eigen 布局、matrix-free、局部块与预分配。
- [验证协议](VALIDATION.md)、[实验结果（1D/2D）](RESULTS.md)、[内存估算](MEMORY_ESTIMATE.md)。
- `accuracy_test/1d/`：独立 1D 收敛测试；`speed_test/1d/`：独立 1D Slurm 速度测试。

## 构建

C++17、CMake >=3.22、Eigen >=3.4。构建和配置不会运行测试。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DDGH_NATIVE_ARCH=ON
cmake --build build -j 4
ctest --test-dir build --output-on-failure
python3 code_test/1d/test_solver.py --exe build/bin/dg_hdg.exe \
  --output code_test/running/1d_local_unique
```

OSCAR 使用 `bash intel-build.sh`（默认构建 1D），加载 Intel 编译器、MKL、Eigen 和 CMake，
单线程、Release、native、SIMD、IPO；实际编译信息写入构建日志。
Release关闭Eigen断言（`NDEBUG`/`EIGEN_NO_DEBUG`），保留项目`require(...)`检查。
`BUILD_TYPE=Debug BUILD_DIR=<独立目录> bash intel-build.sh`启用断言及Eigen分配检查；`code_test`仅接受可核实的Debug构建。
实验入口先在`build-debug/`构建并运行CTest和code_test，再在`build/`独立构建Release用于accuracy/speed；两套日志和二进制身份分别保留。
速度作业在分配到的 Intel 计算节点编译并运行，避免跨 CPU 的 native 指令差异。

## 1D 运行

```sh
# 测试驱动自动创建 running/ 下的独立目录并记录命令、设备及版本。
BUILD_TYPE=Release BUILD_DIR=build-release bash intel-build.sh
python3 accuracy_test/1d/convergence.py --exe build-release/bin/dg_hdg.exe
python3 speed_test/1d/benchmark.py --exe build-release/bin/dg_hdg.exe --degrees 8
```

默认 DIRK3、block-Jacobi、自然块顺序；`--integrator be` 保留后向 Euler。
`--linear direct` 为小规模稠密参考，主路径始终 matrix-free；
block ILU(0) 继续可用。输出 summary、逐 Newton 更新计时和可选节点解。
`--timing-only` 跳过误差积分与解文件，仍保留残差验收及正性检查。

## 1D 实验顺序

```sh
bash accuracy_test/1d/submit.sh
# 精度与速度独立提交，不自动设置依赖
bash speed_test/1d/submit.sh
```

收敛网格 32/64/128/256/512，k=2、DIRK3、CFL=2、T=0.2。
速度网格 512，BE、CFL=5、每次运行四次 Newton 更新，首个预热，后三个均值。
脚本保留命令、源码快照、编译日志、CPU/编译器信息和可执行文件哈希，
每次作业独立输出，不覆盖旧实验。8 倍目标不能替代测量结果。

## 2D 张量积实现

```sh
BUILD_TYPE=Debug BUILD_1D=OFF BUILD_2D=ON BUILD_DIR=build-2d-debug bash intel-build.sh
ctest --test-dir build-2d-debug --output-on-failure
python3 code_test/2d/test_solver.py --exe build-2d-debug/bin/dg_hdg_2d.exe \
  --output code_test/running/2d_local_unique
BUILD_TYPE=Release BUILD_1D=OFF BUILD_2D=ON BUILD_DIR=build-2d-intel bash intel-build.sh
python3 accuracy_test/2d/convergence.py --exe build-2d-intel/bin/dg_hdg_2d.exe
```

2D 采用旋转周期简单波和参考项目相同的中心矩形计算域，算例内部固定 alpha=45°、gamma=3，不提供旋转角命令行输入。
体积 Jv 为方向一维收缩，缓存仅 u/v/H 三场；HDG 采用单元流式精确凝聚和多右端回代，
面耦合按限制到面、点上小矩阵、提升分解，避免完整 C 矩阵乘法。
块 Jacobi 对每个完整单元/面块执行精确 LU 求解，只保存一份 LU，复用装填矩阵。
张量积用于精确算子收缩、面限制/提升和局部装填，保持完整的方向耦合。

[2D 正确性](code_test/2d/README.md)、[收敛](accuracy_test/2d/README.md)、
[256×256 Slurm 速度测试与 100 GiB 内存预算](speed_test/2d/README.md)。

```sh
bash accuracy_test/2d/submit.sh
bash speed_test/2d/submit.sh
```

精度网格 N=32/64/128/256；默认速度入口为256²、k=8→0。另已完成512²、GMRES(30)独立作业批次。

实验代码版本、参数与 CPU 的追溯规则见 [AGENTS.md](AGENTS.md)。

已完成2D实验：Q2收敛至N=256，密度阶2.8993；256²/GMRES(20)、k=8测得DG/HDG为5.3598倍。最新512²/GMRES(30)完整18项、9组配对，k=8测得5.2949倍，峰值RSS为239.0270/153.5564 GiB；k=8两方法在不同节点运行。实际残差验收口径及完整数据见 [RESULTS.md](RESULTS.md#2d)。

新构建产物统一为 `<构建目录>/bin/*.exe` 和 `<构建目录>/*.a`；Debug、Release分别保存。速度批次运行副本为 `<批次>/bin/solver.exe`，历史运行路径保持原始记录。

脚本与目录约定见 [scripts/README.md](scripts/README.md)。本地 code_test 省略 `--output` 自动创建 `code_test/running/` 下的独立目录；结束后用 `python3 scripts/archive_finished.py --run-root <目录>` 显式归档。
