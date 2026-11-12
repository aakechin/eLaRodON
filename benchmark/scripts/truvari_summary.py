"""
Table of Truvari results: one row per <dir>/<tool>/summary.json.

Usage: python3 truvari_summary.py <truvari dir>
"""
import glob
import json
import os
import sys

FIELDS = ['TP-base', 'FN', 'TP-comp', 'FP', 'recall', 'precision', 'f1']

print('\t'.join(['tool'] + FIELDS))
for path in sorted(glob.glob(os.path.join(sys.argv[1], '*', 'summary.json'))):
    summary = json.load(open(path))
    values = [summary.get(f) for f in FIELDS]
    print('\t'.join([os.path.basename(os.path.dirname(path))] +
                    ['NA' if v is None else f'{v:.3f}' if isinstance(v, float) else str(v) for v in values]))
