"""Filesystem-only regression checks; never launch a solver or Slurm job."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import experiment
import archive_finished
import provenance


class ArtifactLayoutTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.exe = self.root/'build/bin/solver.exe'
        self.exe.parent.mkdir(parents=True)
        self.exe.write_bytes(b'layout fixture only')
        (self.root/'build/CMakeCache.txt').write_text('CMAKE_BUILD_TYPE:STRING=Debug\n')
        experiment._OUTPUT = None
        self.addCleanup(setattr, experiment, '_OUTPUT', None)

    def test_rejects_old_executable_layout(self):
        for exe in [self.root/'solver.exe', self.root/'build/bin/solver']:
            with self.assertRaisesRegex(RuntimeError, 'bin/.*exe'):
                experiment.require_solver_layout(exe)
        self.assertEqual(experiment.require_solver_layout(self.exe), self.exe)

    def test_checks_keep_debug_requirement(self):
        (self.root/'build/CMakeCache.txt').write_text('CMAKE_BUILD_TYPE:STRING=Release\n')
        with self.assertRaisesRegex(RuntimeError, 'requires a Debug'):
            experiment.prepare_checks(None, self.exe, '1d')

    def test_invalid_output_does_not_create_directory(self):
        running = self.root/'speed_test/running'
        outside = self.root/'speed_test/results/new'
        with patch.dict('os.environ', {'DGH_RESULT_ROOT': str(running/'submission')}):
            with self.assertRaisesRegex(RuntimeError, 'inside'):
                experiment.prepare_experiment(outside, self.exe, 'speed', '2d')
        self.assertFalse(outside.exists())

    def test_embedded_checks_remain_in_parent_run(self):
        parent = self.root/'speed_test/running/submission/job_1'
        with patch.dict('os.environ', {'DGH_CHECK_ROOT': str(parent)}), patch.object(experiment, 'prepare') as prepare:
            experiment.prepare_checks(None, self.exe, '2d')
            prepare.assert_called_once_with(parent/'checks', self.exe)
            with self.assertRaisesRegex(RuntimeError, 'checks'):
                experiment.prepare_checks(parent/'elsewhere', self.exe, '2d')

    def test_nonempty_output_is_never_overwritten(self):
        output = self.root/'existing'
        output.mkdir()
        (output/'keep.csv').write_text('unique evidence')
        with self.assertRaisesRegex(RuntimeError, 'empty'):
            experiment.prepare(output, self.exe)
        self.assertEqual((output/'keep.csv').read_text(), 'unique evidence')

    def test_campaign_copy_keeps_build_configuration_after_move(self):
        campaign = self.root/'campaign'
        exe = campaign/'bin/solver.exe'
        exe.parent.mkdir(parents=True)
        exe.write_bytes(b'fixture')
        build = campaign/'job_1/build'
        build.mkdir(parents=True)
        (build/'CMakeCache.txt').write_text('CMAKE_BUILD_TYPE:STRING=Release\n')
        (campaign/'build_identity.json').write_text(json.dumps(dict(build_directory='job_1/build')))
        moved = self.root/'archived'
        campaign.rename(moved)
        with patch.object(provenance, 'source_identity', return_value={}):
            info = provenance.runtime(moved/'bin/solver.exe')
        self.assertIn('CMAKE_BUILD_TYPE:STRING=Release', info['cmake_cache'])

    def test_exit_status_survives_failure(self):
        experiment._OUTPUT = self.root
        def fail():
            raise RuntimeError('fixture failure')
        with self.assertRaisesRegex(RuntimeError, 'fixture failure'):
            experiment.run_main(fail)
        self.assertEqual((self.root/'exit_code.txt').read_text(), '1\n')
        experiment.run_main(lambda: None)
        self.assertEqual((self.root/'exit_code.txt').read_text(), '0\n')

    def test_local_archive_requires_explicit_selection_and_dead_process(self):
        folder = self.root/'running/local'
        folder.mkdir(parents=True)
        (folder/'local_process.json').write_text(json.dumps(dict(pid=123, host=archive_finished.platform.node())))
        (folder/'exit_code.txt').write_text('1\n')
        archive_finished.archive(self.root)
        self.assertTrue(folder.exists())
        with patch.object(archive_finished.os, 'kill', return_value=None):
            archive_finished.archive(self.root, folder)
        self.assertTrue(folder.exists())
        response = unittest.mock.Mock(stdout='', returncode=0)
        with patch.object(archive_finished.os, 'kill', side_effect=ProcessLookupError), patch.object(archive_finished.subprocess, 'run', return_value=response):
            archive_finished.archive(self.root, folder)
        self.assertFalse(folder.exists())
        self.assertEqual((self.root/'results/local/exit_code.txt').read_text(), '1\n')


if __name__ == '__main__':
    unittest.main()
