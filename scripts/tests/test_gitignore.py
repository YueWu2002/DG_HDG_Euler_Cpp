"""Guard retention of provenance while excluding snapshots and compiler output."""
from pathlib import Path
import subprocess
import unittest


class GitIgnoreTests(unittest.TestCase):
    def test_archive_retention_matrix(self):
        project = Path(__file__).resolve().parents[2]
        paths = {
            'job_1/build/CMakeCache.txt': False,
            'job_1/build/compile_commands.json': False,
            'job_1/build/CMakeFiles/CMakeConfigureLog.yaml': False,
            'job_1/build/CMakeFiles/3.31.8/CMakeCXXCompiler.cmake': False,
            'job_1/build/CMakeFiles/3.31.8/CMakeSystem.cmake': False,
            'job_1/build/Testing/Temporary/LastTest.log': False,
            'job_1/build/CMakeFiles/euler_2d.dir/build.make': True,
            'job_1/build/CMakeFiles/euler_2d.dir/link.txt': True,
            'job_1/build/Makefile': True,
            'job_1/build/dg_hdg_2d/dg_hdg_2d': True,
            'job_1/build/dg_hdg_2d/check_2d': True,
            'job_1/build/libeuler_2d.a': True,
            'job_1/build/bin/dg_hdg_2d.exe': True,
            'job_1/build-debug/CMakeFiles/CTest.log': False,
            'job_1/data/summary.csv': False,
            'scheduler/job.err': False,
            'source/CMakeLists.txt': True,
            'source.sha256': False,
            'driver.sha256': False,
        }
        for tree in ['accuracy_test/results/a', 'speed_test/results/a', 'code_test/running/a']:
            for name, expected in paths.items():
                with self.subTest(tree=tree, path=name):
                    result = subprocess.run(['git', '-C', str(project), 'check-ignore',
                                             '--no-index', '-q', tree+'/'+name])
                    self.assertIn(result.returncode, [0, 1])
                    self.assertEqual(result.returncode == 0, expected)
