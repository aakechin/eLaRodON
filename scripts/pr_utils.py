"""
Precision-recall curve over QUAL thresholds.

A call is a pair (qual, truth_ids): truth_ids is the set of golden-set variants the call
matched (empty set for a false call). At threshold t all calls with QUAL >= t are kept:
    precision = matched calls / kept calls
    recall    = distinct golden-set variants matched by kept calls / all golden-set variants
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# fixed colour per program, the same on every panel
COLORS = {
    'eLaRodON': '#d62728', 'NanoSV': '#1f77b4', 'NanoVar': '#ff7f0e', 'SVIM': '#2ca02c',
    'Sniffles2': '#9467bd', 'Sniffles2_mosaic': '#8c564b', 'Severus': '#17becf',
}
# low-support runs: colour of the default run, dashed line
LOWSUPPORT_BASE = {'Sniffles2_mosaic_germline': 'Sniffles2', 'NanoVar_mincov1': 'NanoVar',
                   'NanoSV_cluster1': 'NanoSV', 'Severus_minsupport1': 'Severus'}


def pr_curve(calls, n_truth):
    """Returns lists of thresholds, precision and recall (one point per distinct QUAL)."""
    calls = sorted(calls, key=lambda c: -c[0])
    thresholds, precision, recall = [], [], []
    kept = matched_calls = 0
    found = set()
    i = 0
    while i < len(calls):
        qual = calls[i][0]
        while i < len(calls) and calls[i][0] == qual:
            kept += 1
            if calls[i][1]:
                matched_calls += 1
                found |= calls[i][1]
            i += 1
        thresholds.append(qual)
        precision.append(matched_calls / kept)
        recall.append(len(found) / n_truth)
    return thresholds, precision, recall


def average_precision(precision, recall):
    """Sum of (R_k - R_k-1) * P_k, as in sklearn.metrics.average_precision_score."""
    ap, prev_r = 0.0, 0.0
    for p, r in zip(precision, recall):
        ap += (r - prev_r) * p
        prev_r = r
    return ap


def best_f1(thresholds, precision, recall):
    best = (0.0, None, None, None)
    for t, p, r in zip(thresholds, precision, recall):
        f1 = 2 * p * r / (p + r) if p + r else 0.0
        if f1 > best[0]:
            best = (f1, t, p, r)
    return best


def summarize(program, calls, n_truth):
    t, p, r = pr_curve(calls, n_truth)
    f1, thr, p_at, r_at = best_f1(t, p, r)
    return dict(program=program, calls=len(calls), matched_calls=sum(1 for c in calls if c[1]),
                golden_set=n_truth, found=round(r[-1] * n_truth), recall_all_calls=r[-1],
                precision_all_calls=p[-1], average_precision=average_precision(p, r),
                best_f1=f1, best_f1_qual=thr, precision_at_best_f1=p_at, recall_at_best_f1=r_at), (p, r)


def curve_table(program, calls, n_truth):
    t, p, r = pr_curve(calls, n_truth)
    return pd.DataFrame({'program': program, 'qual': t, 'precision': p, 'recall': r})


def plot_pr(ax, curves, title):
    """curves: table with program, precision, recall (one point per QUAL)."""
    for program, c in curves.groupby('program', sort=False):
        ap = average_precision(c['precision'], c['recall'])
        # precision of point k holds on the recall interval that ends at point k
        ax.plot([0] + list(c['recall']), [c['precision'].iloc[0]] + list(c['precision']),
                drawstyle='steps-pre', lw=2, color=COLORS[LOWSUPPORT_BASE.get(program, program)],
                linestyle='--' if program in LOWSUPPORT_BASE else '-',
                label=f'{program.replace("_", " ")} (AP={ap:.2f})')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(loc='upper right', fontsize=9)
