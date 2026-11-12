"""Check Docker arguments without accessing the Docker daemon."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

BENCH = Path(__file__).resolve().parents[1]


class LauncherTests(unittest.TestCase):
    def test_reference_and_output_paths_with_spaces(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            refdir = root / 'reference with spaces'
            refdir.mkdir()
            ref = refdir / 'hg19 genome.fa'
            ref.write_text('>chr22\nA\n')
            Path(str(ref) + '.fai').write_text('chr22\t1\t7\t1\t2\n')
            output = root / 'results with spaces'
            capture = root / 'args.json'
            docker = root / 'docker'
            docker.write_text('#!/usr/bin/env python3\nimport json, os, sys\n'
                              'open(os.environ["CAPTURE"], "w").write(json.dumps(sys.argv[1:]))\n')
            docker.chmod(0o755)
            env = {**os.environ, 'PATH': str(root) + os.pathsep + os.environ['PATH'], 'CAPTURE': str(capture)}
            result = subprocess.run(['bash', str(BENCH / 'run.sh'), '--ref', str(ref),
                                     '--out', str(output), '--tools', 'NanoSV eLaRodON'], env=env,
                                    text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            args = json.loads(capture.read_text())
            self.assertIn(f'type=bind,source={refdir},target=/reference,readonly', args)
            self.assertIn(f'type=bind,source={output},target=/data', args)
            self.assertIn('/reference/hg19 genome.fa', args)
            self.assertIn('NanoSV eLaRodON', args)
            self.assertNotIn('--bam', args)
            self.assertTrue(output.is_dir())


if __name__ == '__main__':
    unittest.main()
