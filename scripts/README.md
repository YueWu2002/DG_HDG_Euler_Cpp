# 脚本产物与运行目录

新构建统一布局（根项目、1D/2D独立构建、Debug和Release均适用）：

```text
<build>/CMakeCache.txt
<build>/compile_commands.json
<build>/lib*.a
<build>/bin/<solver-or-check>.exe
```

`BUILD_DIR`选择独立构建根；Debug与Release必须使用不同目录。CMake管理全部求解器、CTest检查程序和benchmark目标的命名。测试驱动的`--exe`要求`bin/*.exe`；code_test额外要求可追溯到Debug的CMakeCache。

- `accuracy_test/{running,results}/`：精度实验。
- `speed_test/{running,results}/`：速度实验。
- `code_test/{running,results}/`：独立本地验证。
- Slurm实验保留`<submission>/job_<id>/{build,build-debug,data,checks}/`、版本清单和scheduler日志；独立速度campaign的执行副本位于`<submission>/bin/solver.exe`。
- 内嵌code_test只使用所属运行目录的`checks/`，由共享构建入口设置`DGH_CHECK_ROOT`，不另建重复档案。

本地驱动省略`--output`自动创建唯一running目录，显式路径也必须位于对应running中。驱动保存CPU、命令、二进制哈希、构建配置、原始CSV、日志、资源和退出状态；失败同样保留。只有SUCCESS与退出状态均满足验收时才能宣称完成。

`python3 scripts/archive_finished.py --run-root <目录>`显式归档已停止的本地运行；要求本机记录的PID已不存在、退出状态存在且无Slurm路径依赖。PID仍存在或无法核实则保留原位。普通Slurm任务通过job_id、squeue及sacct核验后整体归档；独立速度campaign继续由自身afterany归档任务处理。归档不清空构建目录、不删除日志、配置或唯一数据。仅临时source在任务及依赖结束、清单存在且策略允许时移除。

快照打包排除running/results、bin、build/build-*、CMakeFiles及.exe/.a/.o；保留源码、构建脚本和逐文件SHA256。历史基线无后缀可执行文件的读取路径属于原始证据，不重命名；生成的新副本使用bin/solver.exe。历史输出和在途快照不批量修改。

`.gitignore`只排除编译产物、缓存和重复源码快照，不排除整棵build/bin目录；构建配置、日志、CSV、环境、patch、out/err及元数据均可纳入Git。忽略规则不授权删除磁盘文件。

不启动求解器或Slurm的回归检查：

```sh
python3 -m unittest discover -s scripts/tests -v
```

新源码速度批次可用 `prepare_snapshot_speed.py --required-cpu "<精确CPU型号>"` 在编译前核验指定节点。运行副本通过 `build_identity.json` 中相对路径关联原Release构建配置，整体归档后仍可解析；不把提交节点硬件当作运行节点信息。

## Git归档与同步

执行者应按根目录 [AGENTS.md](../AGENTS.md) 的“Git归档与同步检查”主动核对文件、跟踪状态和忽略规则，无需单独运行审查脚本。具体保留、忽略、索引处理及原字节保存要求统一维护在AGENTS.md。
