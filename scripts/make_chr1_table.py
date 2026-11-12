"""
NA12878 (rel7, 40x, hg19): chr1 calls of each program with a +/- column
(+ = the call matches a variant of the SvABA golden set).

Parsing and matching rules are those of rel7/common_plot_ROC.py (Figure 5B), with two fixes:
  1. DUP is renamed to TD for Sniffles2, NanoVar, SVIM and NanoSV
     (the original had `svtype == 'TD'` or no renaming, so DUP never matched SvABA TD);
  2. each SvABA variant missed by a program is counted once in the ROC (the original counted it twice);
  3. Sniffles2 BND records (no SVLEN) are parsed from ALT, as SVIM BND records.
With --as-published the fixes are off; AUC then equals rel7/results_ROC/metrics_table_rel7_svaba_200.csv.

Rules kept as in the original: calls < 200 bp are not compared; eLaRodON INVTD and BND_INVTD are
not compared; FILTER is not used;
if several calls match one SvABA variant, only the first enters the ROC (all of them get "+" in the table).
Sniffles2 mosaic is parsed like Sniffles2.

Usage: python3 make_chr1_table.py <output.xlsx> [--as-published]
"""
import os
import re
import sys

import pandas as pd
from sklearn.metrics import roc_auc_score

REL7 = '/media/regina/MR_GREY/human/rel7'
SVABA = f'{REL7}/results_svaba_compare/matched_svaba_variants.vcf'
PROGRAMS = [
    ('eLaRodON', f'{REL7}/vcf_hg19/filtered_elarodon/rel7_on_hg19_minimap2_mapont_sorted.40_all_LGRS.vcf'),
    ('NanoSV', f'{REL7}/vcf_hg19/filtered_nanosv/nanosv_rel7_hg19_40.vcf'),
    ('NanoVar', f'{REL7}/vcf_hg19/filtered_nanovar/rel7_on_hg19_minimap2_mapont_sorted.40.nanovar.pass.vcf'),
    ('SVIM', f'{REL7}/vcf_hg19/filtered_svim/SVIM_rel7_hg19_40.vcf'),
    ('Sniffles2', f'{REL7}/vcf_hg19/filtered_sniffles/sniffles_rel7_40_hg19.vcf'),
    ('Sniffles2_mosaic', f'{REL7}/vcf_hg19/filtered_sniffles.mosaic/sniffles_rel7_mosaic_40_hg19.vcf'),
    ('Severus', f'{REL7}/vcf_hg19/filtered_severus/severus_all.vcf'),
]
# NA12878_REVISION=1: VCF files of the revision runs (new eLaRodON, NanoVar 1.8.3, Sniffles2 2.7.1 with
# --tandem-repeats) instead of the article ones
REVISION_DIR = '/media/regina/MR_GREY/human/rel7/New_files/na12878_revision'
if os.environ.get('NA12878_REVISION') == '1':
    REVISION_VCFS = {'eLaRodON': f'{REVISION_DIR}/elarodon_revision/rel7_on_hg19_minimap2_mapont_sorted.40_all_LGRS.vcf',
                     'NanoVar': f'{REVISION_DIR}/nanovar183/rel7_on_hg19_minimap2_mapont_sorted.40.nanovar.pass.vcf',
                     'Sniffles2': f'{REVISION_DIR}/sniffles271_tr/sniffles271_tr_rel7_40_hg19.vcf'}
    PROGRAMS = [(p, REVISION_VCFS.get(p, path)) for p, path in PROGRAMS]
    SVABA = '/media/regina/MR_GREY/human/rel7/New_files/goldset_rebuild/svaba_compare_200/matched_svaba_variants.vcf'      # rebuilt reference set (goldset_rebuild/join_all_fixed.py)
TABLE_CHROM = 'chr1'
EXCEL_CELL_MAX = 32767
VCF_COLUMNS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
FIXES = '--as-published' not in sys.argv


def digits(s):
    return re.sub(r'\D', '', s)


def chr_name(c):
    try:
        return 'chr' + str(int(c))
    except ValueError:
        return 'chr' + c if c in ('X', 'Y') else c


def info_of(col):
    return dict(item.split('=') for item in col[7].split(';') if '=' in item)


def bracket_type(alt):
    """SvABA/NanoSV/SVIM rule: ']..' -> TD, '..[' -> DEL, other brackets -> BND_INV."""
    if alt[0] == ']':
        return 'TD', alt[1:alt.find(':')]
    if alt[-1] == '[':
        return 'DEL', alt[alt.find('[') + 1:alt.find(':')]
    if alt[0] == '[':
        return 'BND_INV', alt[1:alt.find(':')]
    if alt[-1] == ']':
        return 'BND_INV', alt[alt.find(']') + 1:alt.find(':')]
    raise ValueError(f'unexpected ALT {alt}')


def read_svaba(path):
    """Each SvABA record gives two rows: as written and with ends swapped if POS2 < POS1."""
    rows = []
    for line in open(path):
        if line.startswith('#'):
            continue
        col = line.rstrip('\n').split('\t')
        info = info_of(col)
        alt, ref = col[4], col[3]
        c1 = chr_name(col[0])
        c2 = c1
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
        svlen = abs(svlen)
        rec = dict(ID=col[2], CHR1=c1, CHR2=c2, POS1=p1, POS2=p2, SVLEN=svlen, SVTYPE=svtype,
                   key=(c1, c2, p1, p2, svlen, svtype), line=line.rstrip('\n'))
        rows.append(dict(rec))
        if p2 < p1:
            rec.update(CHR1=c2, CHR2=c1, POS1=p2, POS2=p1)
        rows.append(rec)
    return rows


# ---------- parsers: return (svtype, chrom1, pos1, chrom2, pos2, svlen, qual) or a reason to skip ----------

def parse_elarodon(col):
    info = info_of(col)
    svtype = info['SVTYPE'][1:-1]
    if svtype in ('BND_INVTD', 'INVTD'):
        return 'INVTD/BND_INVTD excluded'
    p1 = int(col[1])
    return (svtype, col[0], p1, info.get('CHROM2', col[0]), int(info.get('END', p1)),
            abs(int(info.get('SVLEN', 0))), float(col[5]))


def parse_nanosv(col):
    try:
        qual = float(col[5])
    except ValueError:
        return 'QUAL is not a number'
    info = info_of(col)
    c1, p1, alt = col[0], int(col[1]), col[4]
    if '<' not in alt:
        p2 = int(digits(alt.split(':')[1]))
        svtype, c2 = bracket_type(alt)
        if c2 != c1:
            svtype, svlen = 'BND_TRL', 0
        else:
            svlen = abs(int(info.get('SVLEN', 0)))
    else:
        svtype = alt[1:-1]
        svlen = abs(int(info.get('SVLEN', 0)))
        c2 = info.get('CHROM2', c1)
        p2 = int(info.get('END', p1))
    if FIXES and svtype == 'DUP':
        svtype = 'TD'
    return svtype, c1, p1, c2, p2, svlen, qual


def parse_nanovar(col):
    info = info_of(col)
    p1 = int(col[1])
    svtype = 'TD' if FIXES and info['SVTYPE'] == 'DUP' else info['SVTYPE']
    if FIXES and svtype == 'BND':
        # interchromosomal BND: the mate is only in ALT (CHR2, END = POS + 1)
        m = re.search(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]', col[4])
        if m and m.group(1) != col[0]:
            return 'BND_TRL', col[0], p1, m.group(1), int(m.group(2)), 0, float(col[5])
    return (svtype, col[0], p1, info.get('CHROM2', col[0]), int(info.get('END', p1)),
            abs(int(digits(info['SVLEN']))), float(col[5]))


def parse_svim(col):
    info = info_of(col)
    c1, p1, alt = col[0], int(col[1]), col[4]
    c2 = info.get('CHROM2', c1)
    p2 = int(info.get('END', p1))
    svlen = abs(int(info.get('SVLEN', p2 - p1 + 1)))
    svtype = info['SVTYPE']
    if svtype == 'BND':
        p2 = int(digits(alt.split(':')[1]))
        svtype, c2 = bracket_type(alt)
        if c2 != c1:
            svtype, svlen = 'BND_TRL', 0
        else:
            svlen = p2 - p1
    if FIXES and 'DUP' in svtype:
        svtype = 'TD'
    return svtype, c1, p1, c2, p2, svlen, float(col[5])


def parse_sniffles(col):
    info = info_of(col)
    c1, p1, alt = col[0], int(col[1]), col[4]
    svtype = info['SVTYPE']
    c2 = info.get('CHROM2', c1)
    p2 = int(info.get('END', p1))
    svlen = abs(int(info.get('SVLEN', 0)))
    if FIXES and svtype == 'BND':
        # same rule as for SVIM BND
        p2 = int(digits(alt.split(':')[1]))
        svtype, c2 = bracket_type(alt)
        if c2 != c1:
            svtype, svlen = 'BND_TRL', 0
        else:
            svlen = p2 - p1
    if FIXES and svtype == 'DUP':
        svtype = 'TD'
    return svtype, c1, p1, c2, p2, svlen, float(col[5])


def parse_severus(col):
    info = info_of(col)
    c1, p1, alt = col[0], int(col[1]), col[4]
    c2 = info.get('CHROM2', c1)
    p2 = int(info.get('END', p1))
    if 'SVLEN' in info:
        try:
            svlen = abs(int(digits(info['SVLEN'])))
        except ValueError:
            svlen = abs(p2 - p1)
    else:
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
    if svtype == 'DUP':
        svtype = 'TD'
    return svtype, c1, p1, c2, p2, svlen, float(col[5])


PARSERS = {'eLaRodON': parse_elarodon, 'NanoSV': parse_nanosv, 'NanoVar': parse_nanovar,
           'SVIM': parse_svim, 'Sniffles2': parse_sniffles, 'Sniffles2_mosaic': parse_sniffles,
           'Severus': parse_severus}


def find_match(call, truth_by_chr1):
    """First SvABA row (in file order) that matches the call, or None."""
    svtype, c1, p1, c2, p2, svlen, _ = call
    err = round(svlen * 0.1) if svlen > 0 else 100
    if FIXES and c1 != c2:
        err = 100      # junction between chromosomes: 100 bp at both ends, whatever the call length
    for t in truth_by_chr1.get(c1, []):
        if t['CHR2'] != c2:
            continue
        if not (t['SVTYPE'] in svtype or svtype in t['SVTYPE']):
            continue
        if t['SVTYPE'] == 'INS':
            if abs(t['POS1'] - p1) <= err and abs(t['SVLEN'] - svlen) <= err:
                return t
        elif abs(t['POS1'] - p1) <= err and abs(t['POS2'] - p2) <= err:
            return t
    return None


def truth_key(t):
    return f"{t['SVTYPE']}_{t['CHR1']}_{t['POS1']}_{t['CHR2']}_{t['POS2']}_{t['SVLEN']}"


def evaluate(program, path, truth, truth_by_chr1):
    parse = PARSERS[program]
    record_keys = {t['key'] for t in truth}   # SvABA records as written in the file
    used_keys = set()      # SvABA records already matched by an earlier call
    found_keys = set()     # matched SvABA rows of calls that entered the ROC (Real_Info)
    roc_scores, roc_labels = [], []
    roc_calls = []         # (QUAL, {SvABA ID} or empty set) for precision-recall
    table = []
    for line in open(path):
        if line.startswith('#'):
            continue
        col = line.rstrip('\n').split('\t')
        call = parse(col)
        # interchromosomal records have no length and are not filtered by it
        if isinstance(call, str) or (call[5] < 200 and not (FIXES and call[1] != call[3])):
            continue
        t = find_match(call, truth_by_chr1)
        if t is None and FIXES and call[0] == 'TRL' and call[1] != call[3]:
            # eLaRodON TRL: fragment c2:[END, END+SVLEN-1] inserted at c1:POS, the second junction is at its far end
            t = find_match(call[:4] + (call[4] + call[5] - 1,) + call[5:], truth_by_chr1)
        in_roc = True
        if t is not None:
            row_key = (t['CHR1'], t['CHR2'], t['POS1'], t['POS2'], t['SVLEN'], t['SVTYPE'])
            reverse_key = (t['CHR2'], t['CHR1'], t['POS2'], t['POS1'], t['SVLEN'], t['SVTYPE'])
            key = row_key if row_key in record_keys else reverse_key if reverse_key in record_keys else None
            if key is None or key in used_keys:
                in_roc = False
            if key is not None:
                used_keys.add(key)
        if in_roc:
            roc_scores.append(call[6])
            roc_labels.append(int(t is not None))
            roc_calls.append((call[6], {t['ID']} if t is not None else set()))
            if t is not None:
                found_keys.add(truth_key(t))
        if col[0] == TABLE_CHROM:
            table.append(dict(call=call, col=col, line=line.rstrip('\n'), truth=t, in_roc=in_roc))
    calls_in_roc = len(roc_scores)
    tp_calls = sum(roc_labels)
    # SvABA variants not found by the program enter the ROC as positives with QUAL = -1.
    # Each variant has two rows in `truth`; the original script counted every missed row.
    found_ids = {t['ID'] for t in truth if truth_key(t) in found_keys}
    if FIXES:
        missed = len({t['ID'] for t in truth} - found_ids)
    else:
        missed = sum(1 for t in truth if truth_key(t) not in found_keys)
    roc_scores += [-1.0] * missed
    roc_labels += [1] * missed
    stats = dict(program=program, auc=roc_auc_score(roc_labels, roc_scores),
                 calls_compared=calls_in_roc, true_calls=tp_calls, svaba_missed=missed)
    return table, found_ids, stats, roc_calls


def excel_cell(value):
    if len(value) <= EXCEL_CELL_MAX:
        return value
    note = f' [truncated, full length {len(value)}]'
    return value[:EXCEL_CELL_MAX - len(note)] + note


def program_sheet(table):
    rows = []
    for r in table:
        row = {name: excel_cell(v) for name, v in zip(VCF_COLUMNS, r['col'])}
        row['SvABA'] = '+' if r['truth'] else '-'
        rows.append(row)
    df = pd.DataFrame(rows)
    for name in ('POS',):
        df[name] = df[name].astype(int)
    return df


def main(out_xlsx):
    truth = read_svaba(SVABA)
    truth_by_chr1 = {}
    for t in truth:
        truth_by_chr1.setdefault(t['CHR1'], []).append(t)

    all_stats = []
    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as xw:
        for program, path in PROGRAMS:
            table, _, stats, _ = evaluate(program, path, truth, truth_by_chr1)
            program_sheet(table).to_excel(xw, sheet_name=program, index=False)
            all_stats.append(stats)
            print(stats, flush=True)
        for ws in xw.book.worksheets:
            ws.freeze_panes = 'A2'
            for cells in ws.columns:
                ws.column_dimensions[cells[0].column_letter].width = 60 if cells[0].value in ('INFO', 'REF', 'ALT') else 14
    pd.DataFrame(all_stats).to_csv(out_xlsx.replace('.xlsx', '_auc.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main(sys.argv[1])
