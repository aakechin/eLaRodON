"""Runner integration tests: real VCF normalization/sort, simulated Truvari executable.

The stand-in checks inputs and FILTER policy, not Truvari's scientific matching.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

BENCH = Path(__file__).resolve().parents[1]


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.scripts = self.root / 'scripts'
        shutil.copytree(BENCH / 'scripts', self.scripts)
        self.ref = self.root / 'hg19.fa'
        self.ref.write_text(''.join(f'>chr{c}\n' + 'A' * 2000 + '\n' for c in list(range(1, 23)) + ['X', 'Y']))
        subprocess.run(['samtools', 'faidx', str(self.ref)], check=True)
        truth = self.scripts / 'truth/na12878_golden_set.vcf'
        truth.write_text('22\t100\tgold\t' + 'A' * 251 + '\tA\t60\tPASS\t.\tGT\t0/1\n')
        self.out = self.root / 'out'
        (self.out / 'vcf').mkdir(parents=True)
        (self.out / 'vcf/NanoSV.vcf').write_text(
            '##cmdline=test\n##fileformat=VCFv4.2\n'
            '##contig=<ID=chr22,length=2000>\n'
            '##FILTER=<ID=MapQual,Description="Low mapping quality">\n'
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            'chr22\t100\tnano\tA\t<DEL>\t60\tMapQual\tSVTYPE=DEL;END=350;SVLEN=-250\n')
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        self.executable('micromamba', '#!/bin/bash\nshift 3\nexec "$@"\n')
        self.executable('truvari', '''#!/usr/bin/env python3
import json, pathlib, sys
import pysam
args = sys.argv[1:]
assert '--passonly' not in args, 'caller FILTER must not exclude records'
out = pathlib.Path(args[args.index('-o')+1]); out.mkdir()
vcf = pysam.VariantFile(args[args.index('-c')+1])
records = list(vcf)
assert records, 'filtered input must survive normalization'
if 'NanoSV' in str(vcf.filename):
    assert 'MapQual' in records[0].filter
(out/'summary.json').write_text(json.dumps({'TP-base':0,'FN':1,'TP-comp':0,'FP':len(records), 'recall':0,'precision':0,'f1':0}))
''')
        self.env = {**os.environ, 'PATH': str(self.bin) + os.pathsep + os.environ['PATH']}

    def executable(self, name, text):
        p = self.bin / name
        p.write_text(text)
        p.chmod(0o755)

    def run_runner(self, *extra):
        return subprocess.run(['bash', str(self.scripts / 'run_benchmark.sh'),
            '--ref', str(self.ref), '--out', str(self.out), '--steps', 'truvari',
            '--tools', 'NanoSV', *extra], env=self.env, text=True, capture_output=True)

    def test_two_truth_sets_keep_filtered_calls(self):
        p = self.run_runner()
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        for name in ('golden_set', 'full_svaba'):
            summary = self.out / 'truvari' / name / 'truvari_summary.tsv'
            self.assertTrue(summary.is_file(), p.stdout + p.stderr)
            self.assertIn('NanoSV', summary.read_text())
        manifest = json.loads((self.out / 'truvari/run_manifest.json').read_text())
        self.assertEqual(manifest['filter_policy'], 'all_records')
        self.assertEqual(len(manifest['truth_sets']), 2)

    def test_missing_selected_vcf_fails(self):
        (self.out / 'vcf/NanoSV.vcf').unlink()
        p = self.run_runner()
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn('NanoSV', p.stdout + p.stderr)

    def test_truvari_rejects_non_200_threshold(self):
        p = self.run_runner('--min-length', '50')
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn('200', p.stdout + p.stderr)

    def test_filtered_elarodon_call_survives_normalization(self):
        source = self.out / 'vcf/NanoSV.vcf'
        text = source.read_text().replace('SVLEN=-250', 'SVLEN=-250;J1=R;J2=L')
        (self.out / 'vcf/eLaRodON.vcf').write_text(text)
        p = self.run_runner('--tools', 'eLaRodON')
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        import pysam
        records = list(pysam.VariantFile(str(self.out / 'truvari/prepared/eLaRodON.vcf.gz')))
        self.assertEqual(len(records), 1)

    def test_nanovar_pass_only_output_is_rejected(self):
        shutil.copy(self.out / 'vcf/NanoSV.vcf', self.out / 'vcf/NanoVar.vcf')
        p = self.run_runner('--tools', 'NanoVar')
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('total.vcf', p.stdout + p.stderr)

    def test_existing_comparison_is_preserved(self):
        target = self.out / 'truvari'
        target.mkdir()
        marker = target / 'previous.txt'
        marker.write_text('keep')
        p = self.run_runner()
        self.assertNotEqual(p.returncode, 0)
        self.assertEqual(marker.read_text(), 'keep')

    def test_truvari_failure_propagates(self):
        self.executable('truvari', '#!/bin/bash\nexit 23\n')
        p = self.run_runner()
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)

    def test_call_copy_failure_cannot_report_success_or_reuse_stale_vcf(self):
        bam = self.root / 'input.bam'
        bam.touch()
        Path(str(bam) + '.bai').touch()
        stale = self.out / 'vcf/Sniffles2.vcf'
        stale.write_text('old calls')
        self.executable('cp', '#!/bin/bash\nexit 7\n')
        self.executable('sniffles', '#!/bin/bash\nwhile [ "$#" -gt 0 ]; do\n'
                        'if [ "$1" = -v ]; then printf "new calls" > "$2"; fi\nshift\ndone\n')
        p = self.run_runner('--tools', 'Sniffles2', '--steps', 'call', '--bam', str(bam))
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertFalse(stale.exists())
        self.assertNotIn('records ->', p.stdout)

    def test_svim_traceback_cannot_report_success(self):
        bam = self.root / 'input.bam'
        bam.touch()
        Path(str(bam) + '.bai').touch()
        self.executable('svim', '''#!/bin/bash
printf '##fileformat=VCFv4.2\\n' > "$2/variants.vcf"
printf 'Traceback (most recent call last):\\nAttributeError: plot failed\\n' >&2
exit 0
''')
        p = self.run_runner('--tools', 'SVIM', '--steps', 'call', '--bam', str(bam))
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertFalse((self.out / 'vcf/SVIM.vcf').exists())


if __name__ == '__main__':
    unittest.main()
