"""
One insertion of a mobile element (IS6110 etc.) is reported by eLaRodON as several TRL,
because the inserted part of each read is aligned to a different copy of the element in the reference.
Here such TRL are merged into one record:
  - first breakpoints (CHROM, POS) form a chain with steps <= 30 bp (default -join of eLaRodON);
  - second breakpoint (CHROM2, END) of every member lies in a copy of the same mobile element family
    (mobile_element features of the GenBank record, +-50 bp; family = name before the copy number, IS6110-12 -> IS6110).
The merged record is the member with the highest QUAL, with SR = sum of SR, DP = max DP and QUAL rescaled
by the two eLaRodON QUAL terms that depend on SR (define_type_create_vcf_LRs.py):
    alt_cov = min(1.27*SR/(2.512+SR), 1),  vaf_prob = min(1.218*VAF/(0.095+VAF), 1),  VAF = SR/DP
    QUAL = min(100, QUAL_best * alt_cov(SR_sum)*vaf_prob(SR_sum/DP) / (alt_cov(SR_best)*vaf_prob(SR_best/DP)))
INFO gets MERGED=<number of TRL>. Other records are copied unchanged.

Usage: python3 merge_trl_mobile.py <in_root> <out_root> <mobile_bed>
       (in_root/<sample>/*_all_LGRS.vcf -> out_root/<sample>/<same name>)
"""
import glob
import os
import re
import sys

JOIN_DIST = 30
ELEMENT_PAD = 50


def read_elements(path):
    elements = []
    for line in open(path):
        chrom, start, end, name, kind = line.rstrip('\n').split('\t')
        if kind == 'mobile_element':
            elements.append((chrom, int(start), int(end), re.sub(r"[-'].*$", '', name)))
    return elements


def element_family(elements, chrom, pos):
    for c, start, end, family in elements:
        if c == chrom and start - ELEMENT_PAD <= pos <= end + ELEMENT_PAD:
            return family
    return None


def alt_cov(sr):
    return min(1.27 * sr / (2.512 + sr), 1)


def vaf_prob(vaf):
    return min(1.218 * vaf / (0.095 + vaf), 1)


def info_of(cols):
    return dict(x.split('=', 1) for x in cols[7].split(';') if '=' in x)


def set_info(cols, key, value):
    parts = cols[7].split(';')
    for i, x in enumerate(parts):
        if x.startswith(key + '='):
            parts[i] = f'{key}={value}'
            break
    else:
        parts.append(f'{key}={value}')
    cols[7] = ';'.join(parts)


def merge_group(group):
    """group: list of VCF column lists, all TRL. Returns one column list."""
    best = max(group, key=lambda c: float(c[5]))
    best_info = info_of(best)
    sr_best = int(best_info['SR'])
    sr_sum = sum(int(info_of(c)['SR']) for c in group)
    depth = max(int(info_of(c)['DP']) for c in group)
    scale = alt_cov(sr_sum) * vaf_prob(min(sr_sum / depth, 1)) / (alt_cov(sr_best) * vaf_prob(min(sr_best / depth, 1)))
    merged = list(best)
    merged[5] = str(round(min(100, float(best[5]) * scale), 2))
    set_info(merged, 'SR', sr_sum)
    set_info(merged, 'DP', depth)
    set_info(merged, 'VAF', round(min(sr_sum / depth, 1), 2))
    set_info(merged, 'MERGED', len(group))
    return merged


def process(vcf_in, vcf_out, elements):
    header, records, trl = [], [], []
    for line in open(vcf_in):
        if line.startswith('#'):
            header.append(line)
            continue
        cols = line.rstrip('\n').split('\t')
        info = info_of(cols)
        family = None
        if info['SVTYPE'] == '<TRL>':
            family = element_family(elements, info['CHROM2'], int(info['END']))
        if family is None:
            records.append(cols)
        else:
            trl.append((cols[0], family, int(cols[1]), cols))
    n_in = len(trl)
    trl.sort(key=lambda x: (x[0], x[1], x[2]))
    groups = []
    for chrom, family, pos, cols in trl:
        last = groups[-1] if groups else None
        if last and last['key'] == (chrom, family) and pos - last['pos'] <= JOIN_DIST:
            last['cols'].append(cols)
            last['pos'] = pos
        else:
            groups.append({'key': (chrom, family), 'pos': pos, 'cols': [cols]})
    for g in groups:
        records.append(g['cols'][0] if len(g['cols']) == 1 else merge_group(g['cols']))
    records.sort(key=lambda c: (c[0], int(c[1])))
    with open(vcf_out, 'w') as out:
        for line in header:
            if line.startswith('#CHROM'):
                out.write('##INFO=<ID=MERGED,Number=1,Type=Integer,Description="Number of TRL merged into this record '
                          '(second breakpoints in copies of one mobile element; merge_trl_mobile.py)">\n')
            out.write(line)
        for cols in records:
            out.write('\t'.join(cols) + '\n')
    return n_in, len(groups)


def main(in_root, out_root, bed):
    elements = read_elements(bed)
    total_in = total_out = 0
    for vcf_in in sorted(glob.glob(os.path.join(in_root, '*', '*_all_LGRS.vcf'))):
        sample = os.path.basename(os.path.dirname(vcf_in))
        os.makedirs(os.path.join(out_root, sample), exist_ok=True)
        n_in, n_out = process(vcf_in, os.path.join(out_root, sample, os.path.basename(vcf_in)), elements)
        total_in += n_in
        total_out += n_out
        print(sample, 'TRL in mobile elements:', n_in, '->', n_out)
    print('total:', total_in, '->', total_out)


if __name__ == '__main__':
    main(*sys.argv[1:4])
