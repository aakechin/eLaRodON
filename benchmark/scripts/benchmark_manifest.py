"""Record exactly which inputs and comparison policy were used for both Truvari runs."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys


def file_info(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return {'path': str(path), 'sha256': digest.hexdigest(), 'bytes': path.stat().st_size}


def main():
    directory, reference, chromosome, golden, full, *tools = sys.argv[1:]
    root = Path(directory)
    try:
        version = importlib.metadata.version('truvari')
    except importlib.metadata.PackageNotFoundError:
        version = 'not-installed-as-python-package'
    commit_file = Path('/opt/benchmark/environments/elarodon.commit.txt')
    if commit_file.is_file():
        commit = commit_file.read_text().strip()
    else:
        result = subprocess.run(['git', '-C', '/opt/eLaRodON', 'rev-parse', 'HEAD'],
                                capture_output=True, text=True)
        commit = result.stdout.strip() or None
    data = {
        'filter_policy': 'all_records', 'chromosome': chromosome,
        'reference': file_info(reference), 'reference_fai': file_info(reference + '.fai'),
        'truvari_version': version, 'elarodon_commit': commit,
        'parameters': {'pctseq': 0, 'pctsize': 0.7, 'refdist': 500, 'bnddist': 100,
                       'sizemin': 200, 'sizefilt': 200, 'sizemax': 300000000},
        'truth_sets': {'golden_set': file_info(root / 'truth/golden_set.typed.vcf.gz'),
                       'full_svaba': file_info(full)},
        'golden_set_source': file_info(golden),
        'calls': {tool: file_info(root / f'prepared/{tool}.vcf.gz') for tool in tools},
        'original_calls': {tool: file_info(root.parent / f'vcf/{tool}.vcf') for tool in tools},
        'runner': file_info(Path(__file__).with_name('run_benchmark.sh')),
    }
    (root / 'run_manifest.json').write_text(json.dumps(data, indent=2) + '\n')


if __name__ == '__main__':
    main()
