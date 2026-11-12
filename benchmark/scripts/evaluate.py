"""
Compare SV calls of eLaRodON and other tools with the NA12878 golden set (SvABA records of at least 200 bp
detected by at least one tool), as in Figure 5B of the article, and draw precision-recall and ROC curves.

Matching rules (the article, with the corrections described in the answer to Reviewer 1, comment 1):
  * calls shorter than 200 bp, eLaRodON INVTD and BND_INVTD are not compared; FILTER is not used;
    junctions between chromosomes have no length and are compared whatever the length;
  * a call matches a golden-set record of a compatible type on the same chromosome(s) if both breakpoints
    (for INS: position and length) differ by at most 10% of the call length (100 bp for calls without length
    and for junctions between chromosomes);
  * an eLaRodON TRL between chromosomes (fragment c2:[END, END+SVLEN-1] inserted at c1:POS) has two junctions,
    c1:POS-c2:END and c1:POS-c2:END+SVLEN-1; the call matches if either of them does;
  * DUP is treated as TD; BND records are typed by ALT (']p]t' TD, 't[p[' DEL, other brackets BND_INV,
    another chromosome BND_TRL);
  * if several calls match one golden-set record, only the first one enters the curves;
  * ROC: golden-set records missed by a tool are added as positives with QUAL = -1.
Precision-recall: at QUAL threshold t, precision = matched calls with QUAL >= t / calls with QUAL >= t,
recall = golden-set records matched by these calls / all golden-set records.

Usage:
  python3 evaluate.py --truth golden_set.vcf --out results/evaluation [--chrom chr22] \
      --vcf eLaRodON=elarodon.vcf --vcf Sniffles2=sniffles.vcf ...
Tool names: eLaRodON, NanoSV, NanoVar, SVIM, Sniffles2, Sniffles2_mosaic, Severus.
"""
import argparse
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

MIN_LENGTH = 200
COLORS = {'eLaRodON': '#d62728', 'NanoSV': '#1f77b4', 'NanoVar': '#ff7f0e', 'SVIM': '#2ca02c',
          'Sniffles2': '#9467bd', 'Sniffles2_mosaic': '#8c564b', 'Severus': '#17becf'}


def digits(s):
    return re.sub(r'\D', '', s)


def chr_name(c):
    try:
        return 'chr' + str(int(c))
    except ValueError:
        return 'chr' + c if c in ('X', 'Y') else c


def info_of(col):
    return dict(item.split('=', 1) for item in col[7].split(';') if '=' in item)


def bracket_type(alt):
    """']..' -> TD, '..[' -> DEL, other brackets -> BND_INV; returns the type and the mate chromosome."""
    if alt[0] == ']':
        return 'TD', alt[1:alt.find(':')]
    if alt[-1] == '[':
        return 'DEL', alt[alt.find('[') + 1:alt.find(':')]
    if alt[0] == '[':
        return 'BND_INV', alt[1:alt.find(':')]
    if alt[-1] == ']':
        return 'BND_INV', alt[alt.find(']') + 1:alt.find(':')]
    raise ValueError(f'unexpected ALT {alt}')


def vcf_records(path):
    for line in open(path):
        if not line.startswith('#'):
            yield line.rstrip('\n').split('\t')


def read_truth(path):
    """SvABA records; a record whose second end lies before the first gives a second row with swapped ends."""
    rows = []
    for col in vcf_records(path):
        info = info_of(col)
        alt, ref = col[4], col[3]
        c1 = c2 = chr_name(col[0])
        p1 = int(col[1])
        if not ('DBSNP' in col[7] or 'LOD' in col[7]):
            p2 = int(digits(alt.split(':')[1]))
            svtype, c2 = bracket_type(alt)
            c2 = chr_name(c2)
            if c2 != c1:
                svtype, svlen = 'BND_TRL', 0
            else:
                svlen = p2 - p1 + 1
                if p2 < p1 and svtype in ('DEL', 'TD'):      # right end of a DEL or TD: type of the left end
                    svtype = {'DEL': 'TD', 'TD': 'DEL'}[svtype]
        else:
            svtype = 'INS' if len(alt) - len(ref) >= 0 else 'DEL'
            svlen = int(info['SPAN'])
            p2 = p1 + abs(svlen) if svtype == 'DEL' else p1      # deletion: second end = POS + length
        rec = dict(ID=col[2], CHR1=c1, CHR2=c2, POS1=p1, POS2=p2, SVLEN=abs(svlen), SVTYPE=svtype)
        rows.append(dict(rec))
        if p2 < p1:
            rec.update(CHR1=c2, CHR2=c1, POS1=p2, POS2=p1)
        rows.append(rec)
    return rows


# ---------- parsers: (svtype, chrom1, pos1, chrom2, pos2, svlen, qual) or None to skip ----------

def parse_elarodon(col):
    info = info_of(col)
    svtype = info['SVTYPE'][1:-1]
    if svtype in ('BND_INVTD', 'INVTD'):
        return None
    p1 = int(col[1])
    return (svtype, col[0], p1, info.get('CHROM2', col[0]), int(info.get('END', p1)),
            abs(int(info.get('SVLEN', 0))), float(col[5]))


def parse_nanosv(col):
    try:
        qual = float(col[5])
    except ValueError:
        return None
    info = info_of(col)
    c1, p1, alt = col[0], int(col[1]), col[4]
    if '<' not in alt:
        p2 = int(digits(alt.split(':')[1]))
        svtype, c2 = bracket_type(alt)
        svlen = 0 if c2 != c1 else abs(int(info.get('SVLEN', 0)))
        if c2 != c1:
            svtype = 'BND_TRL'
    else:
        svtype = alt[1:-1]
        svlen = abs(int(info.get('SVLEN', 0)))
        c2 = info.get('CHROM2', c1)
        p2 = int(info.get('END', p1))
    return ('TD' if svtype == 'DUP' else svtype), c1, p1, c2, p2, svlen, qual


def parse_nanovar(col):
    info = info_of(col)
    p1 = int(col[1])
    svtype = 'TD' if info['SVTYPE'] == 'DUP' else info['SVTYPE']
    if svtype == 'BND':
        # interchromosomal BND: the mate is only in ALT (CHR2, END = POS + 1)
        m = re.search(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]', col[4])
        if m and m.group(1) != col[0]:
            return 'BND_TRL', col[0], p1, m.group(1), int(m.group(2)), 0, float(col[5])
    return (svtype, col[0], p1, info.get('CHROM2', col[0]), int(info.get('END', p1)),
            abs(int(digits(info['SVLEN']))), float(col[5]))


def bnd_by_alt(c1, p1, alt):
    p2 = int(digits(alt.split(':')[1]))
    svtype, c2 = bracket_type(alt)
    if c2 != c1:
        return 'BND_TRL', c2, p2, 0
    return svtype, c2, p2, p2 - p1


def parse_svim(col):
    info = info_of(col)
    c1, p1, alt = col[0], int(col[1]), col[4]
    c2 = info.get('CHROM2', c1)
    p2 = int(info.get('END', p1))
    svlen = abs(int(info.get('SVLEN', p2 - p1 + 1)))
    svtype = info['SVTYPE']
    if svtype == 'BND':
        svtype, c2, p2, svlen = bnd_by_alt(c1, p1, alt)
    if 'DUP' in svtype:
        svtype = 'TD'
    return svtype, c1, p1, c2, p2, svlen, float(col[5])


def parse_sniffles(col):
    info = info_of(col)
    c1, p1, alt = col[0], int(col[1]), col[4]
    svtype = info['SVTYPE']
    c2 = info.get('CHROM2', c1)
    p2 = int(info.get('END', p1))
    svlen = abs(int(info.get('SVLEN', 0)))
    if svtype == 'BND':
        svtype, c2, p2, svlen = bnd_by_alt(c1, p1, alt)
    return ('TD' if svtype == 'DUP' else svtype), c1, p1, c2, p2, svlen, float(col[5])


def parse_severus(col):
    info = info_of(col)
    c1, p1, alt = col[0], int(col[1]), col[4]
    c2 = info.get('CHROM2', c1)
    p2 = int(info.get('END', p1))
    try:
        svlen = abs(int(digits(info['SVLEN']))) if 'SVLEN' in info else abs(p2 - p1)
    except ValueError:
        svlen = abs(p2 - p1)
    svtype = info['SVTYPE']
    if svtype == 'BND':
        m = re.search(r'[\[\]]([^:\]\[]+):(\d+)[\[\]]', alt)
        if m:
            c2, p2 = m.group(1), int(m.group(2))
            svlen = abs(p2 - p1)
        else:
            svlen = 0
        if alt[0] == ']':
            svtype = 'TD'
        elif alt[-1] == '[':
            svtype = 'DEL'
        elif alt[0] == '[' or alt[-1] == ']':
            svtype = 'BND_INV'
        if c2 != c1:
            svtype, svlen = 'BND_TRL', 0
    return ('TD' if svtype == 'DUP' else svtype), c1, p1, c2, p2, svlen, float(col[5])


PARSERS = {'eLaRodON': parse_elarodon, 'NanoSV': parse_nanosv, 'NanoVar': parse_nanovar, 'SVIM': parse_svim,
           'Sniffles2': parse_sniffles, 'Sniffles2_mosaic': parse_sniffles, 'Severus': parse_severus}


def find_match(call, truth_by_chr1):
    """First golden-set row (in file order) that matches the call, or None."""
    svtype, c1, p1, c2, p2, svlen, _ = call
    err = round(svlen * 0.1) if svlen > 0 else 100
    if c1 != c2:
        err = 100      # junction between chromosomes: 100 bp at both ends, whatever the call length
    for t in truth_by_chr1.get(c1, []):
        if t['CHR2'] != c2 or not (t['SVTYPE'] in svtype or svtype in t['SVTYPE']):
            continue
        if t['SVTYPE'] == 'INS':
            if abs(t['POS1'] - p1) <= err and abs(t['SVLEN'] - svlen) <= err:
                return t
        elif abs(t['POS1'] - p1) <= err and abs(t['POS2'] - p2) <= err:
            return t
    return None


def evaluate(tool, path, truth_by_chr1, chrom):
    """Calls entering the curves: list of (QUAL, golden-set ID or None)."""
    parse = PARSERS[tool]
    used = set()
    calls = []
    for col in vcf_records(path):
        call = parse(col)
        # interchromosomal records have no length and are not filtered by it
        if call is None or (call[5] < MIN_LENGTH and call[1] == call[3]):
            continue
        if chrom and chrom not in (call[1], call[3]):
            continue
        t = find_match(call, truth_by_chr1)
        if t is None and call[0] == 'TRL' and call[1] != call[3]:
            t = find_match(call[:4] + (call[4] + call[5] - 1,) + call[5:], truth_by_chr1)
        if t is None:
            calls.append((call[6], None))
        elif t['ID'] not in used:
            used.add(t['ID'])
            calls.append((call[6], t['ID']))
    return calls


def pr_curve(calls, n_truth):
    """Thresholds, precision and recall, one point per distinct QUAL from the highest."""
    calls = sorted(calls, key=lambda c: -c[0])
    thresholds, precision, recall = [], [], []
    kept = matched = 0
    found = set()
    i = 0
    while i < len(calls):
        qual = calls[i][0]
        while i < len(calls) and calls[i][0] == qual:
            kept += 1
            if calls[i][1] is not None:
                matched += 1
                found.add(calls[i][1])
            i += 1
        thresholds.append(qual)
        precision.append(matched / kept)
        recall.append(len(found) / n_truth)
    return thresholds, precision, recall


def average_precision(precision, recall):
    """Sum of (R_k - R_k-1) * P_k, as in sklearn.metrics.average_precision_score."""
    ap, prev_r = 0.0, 0.0
    for p, r in zip(precision, recall):
        ap += (r - prev_r) * p
        prev_r = r
    return ap


def roc_data(calls, n_truth):
    found = sum(1 for _, t in calls if t is not None)
    scores = [q for q, _ in calls] + [-1.0] * (n_truth - found)
    labels = [int(t is not None) for _, t in calls] + [1] * (n_truth - found)
    return scores, labels


def summarize(tool, calls, n_truth):
    t, p, r = pr_curve(calls, n_truth)
    f1 = [2 * a * b / (a + b) if a + b else 0.0 for a, b in zip(p, r)]
    k = max(range(len(f1)), key=f1.__getitem__)
    scores, labels = roc_data(calls, n_truth)
    return dict(tool=tool, calls=len(calls), golden_set=n_truth, found=round(r[-1] * n_truth),
                recall_all_calls=r[-1], precision_all_calls=p[-1], average_precision=average_precision(p, r),
                best_f1=f1[k], best_f1_qual=t[k], precision_at_best_f1=p[k], recall_at_best_f1=r[k],
                roc_auc=roc_auc_score(labels, scores))


def main():
    par = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    par.add_argument('--truth', required=True, help='golden-set VCF (SvABA records)')
    par.add_argument('--vcf', action='append', required=True, metavar='TOOL=PATH')
    par.add_argument('--chrom', default='', help='compare only calls with a breakpoint on this chromosome')
    par.add_argument('--out', required=True, help='output directory')
    par.add_argument('--min-length', type=int, default=200,
                     help='minimal call length, bp (default 200; use the golden set built with the same threshold, '
                          'truth/na12878_golden_set_min<L>.vcf)')
    args = par.parse_args()
    os.makedirs(args.out, exist_ok=True)
    global MIN_LENGTH
    MIN_LENGTH = args.min_length

    truth = read_truth(args.truth)
    if args.chrom:
        truth = [t for t in truth if args.chrom in (t['CHR1'], t['CHR2'])]
    truth_by_chr1 = {}
    for t in truth:
        truth_by_chr1.setdefault(t['CHR1'], []).append(t)
    n_truth = len({t['ID'] for t in truth})

    rows, points, curves = [], [], []
    for item in args.vcf:
        tool, path = item.split('=', 1)
        if not os.path.exists(path):
            print(f'{tool}: no VCF ({path}), skipped')
            continue
        calls = evaluate(tool, path, truth_by_chr1, args.chrom)
        if not calls:
            print(f'{tool}: no calls of at least {MIN_LENGTH} bp, skipped')
            continue
        s = summarize(tool, calls, n_truth)
        rows.append(s)
        t, p, r = pr_curve(calls, n_truth)
        points.append(pd.DataFrame({'tool': tool, 'qual': t, 'precision': p, 'recall': r}))
        fpr, tpr, _ = roc_curve(*reversed(roc_data(calls, n_truth)))
        curves.append((tool, s, p, r, fpr, tpr))
        print(f'{tool}: {s["calls"]} calls, found {s["found"]} of {n_truth}, AP {s["average_precision"]:.3f}, '
              f'best F1 {s["best_f1"]:.3f}, AUC {s["roc_auc"]:.3f}', flush=True)

    # one plot per file, legend under the axes
    title = f'NA12878, {args.chrom or "all chromosomes"}, golden set: {n_truth}'
    for name, x, y in (('roc_curves.png', 'False positive rate', 'True positive rate'),
                       ('pr_curves.png', 'Recall', 'Precision')):
        fig, ax = plt.subplots(figsize=(6.5, 7))
        for tool, s, p, r, fpr, tpr in curves:
            if name == 'roc_curves.png':
                ax.plot(fpr, tpr, lw=1.8, color=COLORS[tool], label=f'{tool} (AUC = {s["roc_auc"]:.2f})')
            else:
                # the precision of point k holds on the recall interval that ends at point k
                ax.plot([0] + r, [p[0]] + p, drawstyle='steps-pre', lw=1.8, color=COLORS[tool],
                        label=f'{tool} (AP = {s["average_precision"]:.2f}, F1 = {s["best_f1"]:.2f})')
                ax.plot(s['recall_at_best_f1'], s['precision_at_best_f1'], 'o', color=COLORS[tool], markersize=6,
                        markeredgecolor='black', markeredgewidth=0.6)
        if name == 'roc_curves.png':
            ax.plot([0, 1], [0, 1], color='grey', lw=1, linestyle=':')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.3)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=9, frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, name), dpi=300)
        plt.close(fig)

    pd.DataFrame(rows).to_csv(os.path.join(args.out, 'metrics.tsv'), sep='\t', index=False)
    pd.concat(points).to_csv(os.path.join(args.out, 'pr_points.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
