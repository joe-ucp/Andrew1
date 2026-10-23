"""A Git checkout must retain strict scientific-file integrity checks."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import release


class CheckoutIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.package = self.root / 'package'
        self.evidence = self.package / 'evidence/card001/results/external_summary.csv'
        self.evidence.parent.mkdir(parents=True)
        self.evidence.write_text('policy,events,warned\nM15,8648,5315\n', encoding='utf-8')
        files = ['PACKAGE_SHA256.json', 'RELEASE_ALLOWLIST.json',
                 'evidence/card001/results/external_summary.csv']
        (self.package / 'RELEASE_ALLOWLIST.json').write_text(
            json.dumps({'files': files}), encoding='utf-8')
        hashes = {name: release.sha(self.package / name)
                  for name in files if name != 'PACKAGE_SHA256.json'}
        (self.package / 'PACKAGE_SHA256.json').write_text(json.dumps(hashes), encoding='utf-8')
        self.package_patch = patch.object(release, 'PACKAGE', self.package)
        self.package_patch.start()
        self.addCleanup(self.package_patch.stop)

    def git_directory(self):
        self.git('init', '--quiet')
        return self.package / '.git'

    def git(self, *arguments):
        return subprocess.run(
            ['git', '-c', 'user.name=Checkout test fixture',
             '-c', 'user.email=checkout-test@example.invalid',
             '-C', str(self.package), *arguments],
            capture_output=True, text=True, check=True)

    def ignore_private_directories(self):
        (self.package / '.git/info/exclude').write_text('/_UCP/\n/_private/\n', encoding='utf-8')

    def test_root_git_directory_is_excluded(self):
        metadata = self.git_directory()
        result = release.verify()
        (metadata / 'description').write_text('Updated local Git description\n', encoding='utf-8')
        self.assertEqual(release.verify(), result)
        self.assertEqual(result['files'], 3)
        self.assertEqual(result['hashes_checked'], 2)

    def test_root_worktree_git_file_is_excluded(self):
        gitdir = self.root / 'shared-metadata/worktrees/review'
        gitdir.parent.mkdir(parents=True)
        self.git('init', '--quiet', '--separate-git-dir', str(gitdir))
        self.assertEqual(release.verify()['status'], 'PASS')

    def test_invalid_root_git_file_is_rejected(self):
        (self.package / '.git').write_text('private notes disguised as metadata\n', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'Invalid root Git metadata'):
            release.verify()

    def test_missing_evidence_still_fails_in_git_checkout(self):
        self.git_directory()
        self.evidence.unlink()
        with self.assertRaisesRegex(ValueError, 'Allowlist mismatch: missing='):
            release.verify()

    def test_changed_evidence_still_fails_in_git_checkout(self):
        self.git_directory()
        self.evidence.write_text('policy,events,warned\nM15,8648,5316\n', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'Changed package file: evidence/'):
            release.verify()

    def test_nested_git_directory_is_rejected_even_when_empty(self):
        self.git_directory()
        (self.package / 'evidence/.git').mkdir()
        with self.assertRaisesRegex(ValueError, 'Forbidden release member'):
            release.verify()

    def test_nested_git_file_is_rejected(self):
        (self.package / 'evidence/.git').write_text('gitdir: ../elsewhere\n', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'Forbidden release member'):
            release.verify()

    def test_root_private_directories_are_excluded_without_traversal(self):
        self.git_directory()
        self.ignore_private_directories()
        for name in ['_UCP', '_private']:
            folder = self.package / name / 'historical/.git'
            folder.mkdir(parents=True)
            (folder / 'private-record.txt').write_text('Private history stays local', encoding='utf-8')
        result = release.verify()
        self.assertEqual(result['root_private_directories_excluded'], ['_UCP', '_private'])
        self.assertEqual(result['files'], 3)
        self.assertEqual(result['git_privacy']['status'], 'PASS')

    def test_unignored_root_private_directory_is_rejected(self):
        self.git_directory()
        (self.package / '_private').mkdir()
        with self.assertRaisesRegex(ValueError, 'Root private directory must be ignored'):
            release.verify()

    def test_nested_private_directory_is_rejected_even_when_empty(self):
        (self.package / 'evidence/_private').mkdir()
        with self.assertRaisesRegex(ValueError, 'Forbidden release member'):
            release.verify()

    def test_nested_workflow_directory_is_rejected_even_when_empty(self):
        (self.package / 'evidence/_UCP').mkdir()
        with self.assertRaisesRegex(ValueError, 'Forbidden release member'):
            release.verify()

    def test_root_private_exclusion_does_not_hide_regular_files(self):
        (self.package / '_private').write_text('Not a directory', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'Root private exclusion must be a directory'):
            release.verify()

    def test_staged_private_file_is_rejected_even_if_ignored(self):
        self.git_directory()
        self.ignore_private_directories()
        secret = self.package / '_private/review.txt'
        secret.parent.mkdir()
        secret.write_text('Private notes', encoding='utf-8')
        self.git('add', '-f', '--', '_private/review.txt')
        with self.assertRaisesRegex(ValueError, 'Private paths in Git index'):
            release.verify()

    def test_private_history_is_rejected_after_private_file_is_deleted(self):
        self.git_directory()
        self.ignore_private_directories()
        secret = self.package / '_UCP/workflow.txt'
        secret.parent.mkdir()
        secret.write_text('Private workflow', encoding='utf-8')
        self.git('add', '-f', '--', '_UCP/workflow.txt')
        self.git('commit', '--quiet', '-m', 'Fixture with historical private file')
        self.git('rm', '--quiet', '--', '_UCP/workflow.txt')
        self.git('commit', '--quiet', '-m', 'Fixture removes private file')
        with self.assertRaisesRegex(ValueError, 'Private paths in Git history'):
            release.verify()

    def test_private_history_on_another_branch_is_rejected(self):
        self.git_directory()
        self.ignore_private_directories()
        self.git('add', '--', 'evidence/card001/results/external_summary.csv')
        self.git('commit', '--quiet', '-m', 'Fixture public baseline')
        public_branch = self.git('branch', '--show-current').stdout.strip()
        self.git('checkout', '--quiet', '-b', 'private-fixture')
        secret = self.package / '_private/review.txt'
        secret.parent.mkdir()
        secret.write_text('Private branch material', encoding='utf-8')
        self.git('add', '-f', '--', '_private/review.txt')
        self.git('commit', '--quiet', '-m', 'Fixture private branch')
        self.git('checkout', '--quiet', public_branch)
        with self.assertRaisesRegex(ValueError, 'Private paths in Git history'):
            release.verify()

    def test_non_allowlisted_index_path_is_rejected_after_worktree_deletion(self):
        self.git_directory()
        note = self.package / 'internal-review.txt'
        note.write_text('Unlisted local material', encoding='utf-8')
        self.git('add', '--', note.name)
        note.unlink()
        with self.assertRaisesRegex(ValueError, 'Non-allowlisted paths in Git index'):
            release.verify()

    def test_non_allowlisted_history_is_rejected_after_deletion(self):
        self.git_directory()
        note = self.package / 'old-review.txt'
        note.write_text('Historical internal material', encoding='utf-8')
        self.git('add', '--', note.name)
        self.git('commit', '--quiet', '-m', 'Fixture with non-public history')
        self.git('rm', '--quiet', '--', note.name)
        self.git('commit', '--quiet', '-m', 'Fixture removes non-public file')
        with self.assertRaisesRegex(ValueError, 'Non-allowlisted paths in Git history'):
            release.verify()

    def test_valid_public_index_and_history_are_accepted(self):
        self.git_directory()
        self.git('add', '--', 'evidence/card001/results/external_summary.csv')
        self.git('commit', '--quiet', '-m', 'Fixture public evidence')
        result = release.verify()
        self.assertEqual(result['git_privacy']['index_paths_checked'], 1)
        self.assertEqual(result['git_privacy']['commits_checked'], 1)

    def test_unlisted_private_file_is_rejected(self):
        (self.package / 'private-review.txt').write_text('Internal review notes', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'Allowlist mismatch:.*extra='):
            release.verify()

    def test_root_git_symlink_is_rejected(self):
        outside = self.root / 'outside'
        outside.mkdir()
        try:
            (self.package / '.git').symlink_to(outside, target_is_directory=True)
        except OSError as error:
            self.skipTest(f'Creating filesystem symlinks is unavailable: {error}')
        with self.assertRaisesRegex(ValueError, 'Forbidden release member'):
            release.verify()

    @unittest.skipUnless(sys.platform == 'win32', 'Windows junction-specific check')
    def test_root_git_junction_is_rejected(self):
        outside = self.root / 'outside'
        outside.mkdir()
        result = subprocess.run(
            ['cmd', '/c', 'mklink', '/J', str(self.package / '.git'), str(outside)],
            capture_output=True, text=True, check=False, timeout=10)
        if result.returncode:
            self.skipTest('Creating a Windows junction is unavailable: ' + result.stderr.strip())
        self.assertTrue((self.package / '.git').is_junction())
        with self.assertRaisesRegex(ValueError, 'Forbidden release member'):
            release.verify()

    @unittest.skipUnless(sys.platform == 'win32', 'Windows junction-specific check')
    def test_root_private_junction_is_rejected(self):
        outside = self.root / 'outside'
        outside.mkdir()
        result = subprocess.run(
            ['cmd', '/c', 'mklink', '/J', str(self.package / '_private'), str(outside)],
            capture_output=True, text=True, check=False, timeout=10)
        if result.returncode:
            self.skipTest('Creating a Windows junction is unavailable: ' + result.stderr.strip())
        with self.assertRaisesRegex(ValueError, 'Forbidden release member'):
            release.verify()

    def test_evidence_symlink_is_rejected_even_with_matching_bytes(self):
        outside = self.root / 'outside.csv'
        outside.write_bytes(self.evidence.read_bytes())
        self.evidence.unlink()
        try:
            self.evidence.symlink_to(outside)
        except OSError as error:
            self.skipTest(f'Creating filesystem symlinks is unavailable: {error}')
        with self.assertRaisesRegex(ValueError, 'Forbidden release member'):
            release.verify()

    def test_all_release_commands_accept_root_git_metadata(self):
        self.git_directory()
        actual_run = subprocess.run
        for command in ['verify', 'tests', 'reproduce']:
            with self.subTest(command=command):
                arguments = ['release.py', command]
                if command != 'verify':
                    arguments += ['--output-root', str(self.root / command)]
                if command == 'reproduce':
                    arguments += ['--inputs-root', str(self.root / 'inputs')]
                # Exercise real verification and directory guards; do not launch
                # analyses or require participant data for a checkout-safety test.
                analysis_commands = []
                def dispatch(arguments, **kwargs):
                    if arguments[0] == 'git':
                        return actual_run(arguments, **kwargs)
                    analysis_commands.append(arguments)
                with patch.object(sys, 'argv', arguments), \
                     patch.object(release, 'check_environment', return_value={}), \
                     patch.object(release, 'check_inputs', return_value={}), \
                     patch.object(release.subprocess, 'run', side_effect=dispatch), \
                     contextlib.redirect_stdout(io.StringIO()):
                    release.main()
                self.assertEqual(len(analysis_commands), int(command != 'verify'))


if __name__ == '__main__':
    unittest.main()
