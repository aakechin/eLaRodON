#!/usr/bin/env python3
"""
Convert eLaRodON or SvABA VCF into a standard VCF for third-party SV benchmarking
tools (Truvari, SURVIVOR, Minda, SVanalyzer, witty.er).

Both callers are reduced to the same internal model and written in the same way,
so that format differences do not create false events:
  * junction  - one adjacency (chr1:pos1, side1) <-> (chr2:pos2, side2);
                side R = the joined piece lies left of pos (VCF "t[p[" / "t]p]"),
                side L = the joined piece lies right of pos ("]p]t" / "[p[t")
  * insertion - chr:pos + length (no sequence)
  * inversion - two-boundary eLaRodON INV (interval + its two junctions)

Output modes:
  typed    - intrachromosomal junction -> DEL (R,L) / DUP (L,R) / INV (R,R or L,L),
             interchromosomal junction -> BND pair, INS -> INS
  breakend - every junction -> BND pair (two records with MATEID), INS -> INS

Redundant records describing the same event are collapsed (highest QUAL kept), see
README in docs/vcf_normalization_for_benchmark.md.
"""

import argparse
import bisect
import collections
import gzip
import os
import re
import subprocess
import sys

import pysam

CANONICAL = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
CHR_ORDER = {c: i for i, c in enumerate(CANONICAL)}


def norm_chrom(c):
    if c.startswith('chr'):
        return c
    return 'chr' + c


# ---------------------------------------------------------------- data model

class Event:
    """One unit for collapsing/writing.
    kind: 'INT' (typed interval DEL/DUP/INV), 'INS', 'JNC' (junction for BND pair)"""
    __slots__ = ('kind', 'svtype', 'c1', 'p1', 's1', 'c2', 'p2', 's2', 'length',
                 'qual', 'sr', 'ids', 'orig_types', 'orig_svlen', 'tdrn', 'nmerged')

    def __init__(self, kind, svtype, c1, p1, c2=None, p2=None, s1='', s2='',
                 length=0, qual=0.0, sr=0, rid='', orig_type='', orig_svlen=None,
                 tdrn=None):
        self.kind = kind
        self.svtype = svtype
        self.c1, self.p1, self.s1 = c1, p1, s1
        self.c2, self.p2, self.s2 = c2, p2, s2
        self.length = length
        self.qual = qual
        self.sr = sr
        self.ids = [rid]
        self.orig_types = collections.Counter([orig_type])
        self.orig_svlen = orig_svlen
        self.tdrn = tdrn
        self.nmerged = 1


def junction_events(c1, p1, s1, c2, p2, s2, mode, **kw):
    """Turn one junction into the event written for the chosen mode."""
    # canonical order: (c1, p1) <= (c2, p2)
    if (CHR_ORDER[c1], p1) > (CHR_ORDER[c2], p2):
        c1, p1, s1, c2, p2, s2 = c2, p2, s2, c1, p1, s1
    if mode == 'typed' and c1 == c2:
        if (s1, s2) == ('R', 'L'):
            t = 'DEL'
        elif (s1, s2) == ('L', 'R'):
            t = 'DUP'
        else:
            t = 'INV'
        return Event('INT', t, c1, p1, c1, p2, s1, s2, length=p2 - p1, **kw)
    length = p2 - p1 if c1 == c2 else 0
    return Event('JNC', 'BND', c1, p1, c2, p2, s1, s2, length=length, **kw)


# ---------------------------------------------------------------- eLaRodON

def parse_info(info):
    d = {}
    for item in info.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            d[k] = v
        else:
            d[item] = True
    return d


def read_elarodon(path, args, rep):
    """Return list of raw records (dicts) after filters, with link absorption."""
    recs = []
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            col = line.rstrip('\n').split('\t')
            info = parse_info(col[7])
            t = info['SVTYPE'].strip('<>')
            rep['input', t, col[6]] += 1
            if args.pass_only and col[6] != 'PASS':
                rep['drop_filter', t, col[6]] += 1
                continue
            if t in ('INVTD', 'BND_INVTD') and not args.keep_invtd:
                rep['drop_invtd', t, ''] += 1
                continue
            c1 = col[0]
            c2 = info.get('CHROM2', c1)
            if c1 not in CHR_ORDER or c2 not in CHR_ORDER:
                rep['drop_noncanonical_chrom', t, ''] += 1
                continue
            pos = int(col[1])
            end = int(info.get('END', pos))
            svlen = abs(int(info['SVLEN'])) if 'SVLEN' in info else None
            # length that the standard record will have (for min-length filter)
            if t == 'INS' or t == 'TRL':
                std_len = svlen
            elif c1 != c2:
                std_len = None  # interchromosomal: no length, always kept
            else:
                std_len = abs(end - pos)
            if std_len is not None and std_len < args.min_len:
                rep['drop_minlen', t, ''] += 1
                continue
            sblr = info.get('SBLR', '""').strip('"')
            recs.append({
                'id': col[2], 'type': t, 'c1': c1, 'pos': pos, 'c2': c2, 'end': end,
                'svlen': svlen, 'qual': float(col[5]), 'sr': int(info.get('SR', 0)),
                'j1': info.get('J1', ''), 'j2': info.get('J2', ''),
                's2': info.get('S2', ''), 'd1': int(info['D1']) if 'D1' in info else None,
                'd2': int(info['D2']) if 'D2' in info else None,
                'tdrn': int(info['TDRN']) if 'TDRN' in info else None,
                'links': [x for x in sblr.split(',') if x] if sblr else [],
            })
    if args.absorb_linked:
        recs = absorb_linked(recs, args, rep)
    return recs


TWO_BOUNDARY = {'DEL', 'TD', 'INV', 'TRL', 'INVTD'}


def absorb_linked(recs, args, rep):
    """Drop records that eLaRodON itself linked (SBLR tag) to a better-resolved
    representation of the same event:
      - BND_* (one boundary) linked to a two-boundary record or to an INS;
      - INS linked to TD or TRL (CIGAR insertion = the duplicated/inserted fragment).
    Links are used in both directions; the partner must pass all filters and lie
    within args.link_dist bp."""
    by_id = {r['id']: r for r in recs}
    partners = collections.defaultdict(set)
    for r in recs:
        for x in r['links']:
            if x in by_id and x != r['id']:
                partners[r['id']].add(x)
                partners[x].add(r['id'])

    def near(a, b):
        pa = {(a['c1'], a['pos']), (a['c2'], a['end'])}
        pb = {(b['c1'], b['pos']), (b['c2'], b['end'])}
        return any(ca == cb and abs(xa - xb) <= args.link_dist
                   for ca, xa in pa for cb, xb in pb)

    drop = set()
    for r in recs:
        t = r['type']
        for pid in partners.get(r['id'], ()):
            p = by_id[pid]
            if not near(r, p):
                continue
            if args.trl == 'ins_intra' and t.startswith('BND_') and r['c1'] != r['c2']:
                # junction between chromosomes is kept as a junction (as in SvABA and other callers):
                # a linked INS (sequence from another chromosome) is dropped, a linked interchromosomal TRL
                # gives the same junction itself
                if p['type'] == 'INS':
                    drop.add(pid)
                    rep['drop_linked', 'INS', 'into_' + t] += 1
                elif p['type'] == 'TRL' and p['c1'] != p['c2']:
                    drop.add(r['id'])
                    rep['drop_linked', t, 'into_TRL'] += 1
                    break
                continue
            if t.startswith('BND_') and (p['type'] in TWO_BOUNDARY or p['type'] == 'INS'):
                drop.add(r['id'])
                rep['drop_linked', t, 'into_' + p['type']] += 1
                break
            if t == 'INS' and p['type'] in ('TD', 'TRL'):
                drop.add(r['id'])
                rep['drop_linked', t, 'into_' + p['type']] += 1
                break
    return [r for r in recs if r['id'] not in drop]


def elarodon_events(r, args):
    """Convert one filtered eLaRodON record into standard events."""
    t = r['type']
    kw = dict(qual=r['qual'], sr=r['sr'], rid=r['id'], orig_type=t,
              orig_svlen=r['svlen'], tdrn=r['tdrn'])
    if t == 'INS':
        return [Event('INS', 'INS', r['c1'], r['pos'], length=r['svlen'], **kw)]
    if t == 'INV':
        # read: -->(1) ends at POS-D1 ; <--(2) = [POS, END] ; -->(3) starts at END+D2
        ja = (r['c1'], r['pos'] - r['d1'], 'R', r['c1'], r['end'], 'R')
        jb = (r['c1'], r['pos'], 'L', r['c1'], r['end'] + r['d2'], 'L')
        if args.mode == 'typed':
            return [Event('INT', 'INV', r['c1'], r['pos'], r['c1'], r['end'],
                          length=r['end'] - r['pos'], **kw)]
        return [junction_events(*ja, args.mode, **kw), junction_events(*jb, args.mode, **kw)]
    if t == 'TRL':
        # fragment c2:[END, END+SVLEN-1] inserted at c1:POS; D1<0 = target-site duplication
        if args.mode == 'typed' and (args.trl == 'ins' or (args.trl == 'ins_intra' and r['c1'] == r['c2'])):
            return [Event('INS', 'INS', r['c1'], r['pos'], length=r['svlen'], **kw)]
        d1, ln = r['d1'], r['svlen']
        if d1 >= 0:
            a1, a2 = r['pos'], r['pos'] + d1 - 1
        else:
            a2, a1 = r['pos'], r['pos'] - d1 - 1
        if r['s2'] == '+':
            ja = (r['c1'], a1, 'R', r['c2'], r['end'], 'L')
            jb = (r['c2'], r['end'] + ln - 1, 'R', r['c1'], a2, 'L')
        else:
            ja = (r['c1'], a1, 'R', r['c2'], r['end'] + ln - 1, 'R')
            jb = (r['c2'], r['end'], 'L', r['c1'], a2, 'L')
        return [junction_events(*ja, args.mode, **kw), junction_events(*jb, args.mode, **kw)]
    if t == 'INVTD':
        return []  # only with --keep-invtd; two-boundary INVTD is not representable simply
    # DEL, TD and all BND_*: a single junction given by J1/J2 (sides relative to reference)
    return [junction_events(r['c1'], r['pos'], r['j1'], r['c2'], r['end'], r['j2'],
                            args.mode, **kw)]


# ---------------------------------------------------------------- SvABA

ALT_RE = re.compile(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]')


def bnd_sides(alt):
    """VCF breakend ALT -> (side at POS, side at mate)."""
    if alt.endswith('['):
        return 'R', 'L'   # t[p[
    if alt.endswith(']'):
        return 'R', 'R'   # t]p]
    if alt.startswith(']'):
        return 'L', 'R'   # ]p]t
    return 'L', 'L'       # [p[t


def svaba_events(paths, args, rep):
    events = []
    seen_junctions = set()
    for path in paths:
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                col = line.rstrip('\n').split('\t')
                info = parse_info(col[7])
                c1 = norm_chrom(col[0])
                pos = int(col[1])
                qual = float(col[5]) if col[5] != '.' else 0.0
                if info.get('SVTYPE') == 'BND':
                    t = 'BND'
                    rep['input', t, col[6]] += 1
                    m = ALT_RE.search(col[4])
                    c2, p2 = norm_chrom(m.group(1)), int(m.group(2))
                    if c1 not in CHR_ORDER or c2 not in CHR_ORDER:
                        rep['drop_noncanonical_chrom', t, ''] += 1
                        continue
                    if args.pass_only and col[6] != 'PASS':
                        rep['drop_filter', t, col[6]] += 1
                        continue
                    s1, s2 = bnd_sides(col[4])
                    key = tuple(sorted([(CHR_ORDER[c1], pos, s1), (CHR_ORDER[c2], p2, s2)]))
                    if key in seen_junctions:
                        rep['drop_mate_duplicate', t, ''] += 1
                        continue
                    seen_junctions.add(key)
                    if c1 == c2 and abs(p2 - pos) < args.min_len:
                        rep['drop_minlen', t, ''] += 1
                        continue
                    events.append(junction_events(c1, pos, s1, c2, p2, s2, args.mode,
                                                  qual=qual, rid=col[2], orig_type='BND'))
                else:
                    diff = len(col[4]) - len(col[3])
                    t = 'INS' if diff > 0 else 'DEL'
                    rep['input', t, col[6]] += 1
                    if c1 not in CHR_ORDER:
                        rep['drop_noncanonical_chrom', t, ''] += 1
                        continue
                    if args.pass_only and col[6] != 'PASS':
                        rep['drop_filter', t, col[6]] += 1
                        continue
                    if abs(diff) < args.min_len:
                        rep['drop_minlen', t, ''] += 1
                        continue
                    if t == 'INS':
                        events.append(Event('INS', 'INS', c1, pos, length=diff, qual=qual,
                                            rid=col[2], orig_type='INS', orig_svlen=diff))
                    else:
                        events.append(junction_events(c1, pos, 'R', c1, pos + abs(diff) + 1, 'L',
                                                      args.mode, qual=qual, rid=col[2],
                                                      orig_type='DEL', orig_svlen=abs(diff)))
    return events


# ---------------------------------------------------------------- collapse

def collapse(events, args, rep):
    """Greedy clustering, highest QUAL first: an event joins an existing cluster if
    both ends are within args.window bp of the cluster seed (INS: position within
    window and size ratio >= args.pctsize). Seed (= max QUAL) is written out."""
    groups = collections.defaultdict(list)
    for e in events:
        if e.kind == 'INS':
            key = ('INS', e.c1)
        elif e.kind == 'INT':
            key = ('INT', e.svtype, e.c1)
        else:
            key = ('JNC', e.c1, e.s1, e.c2, e.s2)
        groups[key].append(e)
    out = []
    w = args.window
    for key, evs in groups.items():
        evs.sort(key=lambda e: (-e.qual, -e.sr, e.p1))
        seed_pos, seeds = [], []
        for e in evs:
            i = bisect.bisect_left(seed_pos, e.p1 - w)
            hit = None
            while i < len(seed_pos) and seed_pos[i] <= e.p1 + w:
                s = seeds[i]
                if e.kind == 'INS':
                    ok = min(s.length, e.length) / max(s.length, e.length) >= args.pctsize
                else:
                    ok = abs(s.p2 - e.p2) <= w
                if ok:
                    hit = s
                    break
                i += 1
            if hit is None:
                j = bisect.bisect_left(seed_pos, e.p1)
                seed_pos.insert(j, e.p1)
                seeds.insert(j, e)
            else:
                hit.nmerged += 1
                hit.sr += e.sr
                hit.ids.extend(e.ids)
                hit.orig_types.update(e.orig_types)
                rep['collapsed', key[0] if key[0] != 'INT' else key[1], ''] += 1
        out.extend(seeds)
    return out


# ---------------------------------------------------------------- writing

def bnd_alt(base, side_here, side_mate, c_mate, p_mate):
    p = f'{c_mate}:{p_mate}'
    if side_here == 'R':
        return f'{base}[{p}[' if side_mate == 'L' else f'{base}]{p}]'
    return f']{p}]{base}' if side_mate == 'R' else f'[{p}[{base}'


def write_vcf(events, args, source):
    fa = pysam.FastaFile(args.ref)

    def base(c, p):
        b = fa.fetch(c, max(p - 1, 0), max(p, 1)).upper()
        return b if b in ('A', 'C', 'G', 'T') else 'N'

    lines = []
    for n, e in enumerate(events, 1):
        common = [f'NMERGED={e.nmerged}',
                  'ORIG_TYPES=' + ','.join(f'{k}:{v}' for k, v in sorted(e.orig_types.items())),
                  'ORIG_IDS=' + ','.join(e.ids[:args.max_ids]) + (',...' if len(e.ids) > args.max_ids else '')]
        if e.sr:
            common.append(f'SR={e.sr}')
        if e.orig_svlen is not None:
            common.append(f'ORIG_SVLEN={e.orig_svlen}')
        if e.tdrn is not None:
            common.append(f'TDRN={e.tdrn}')
        q = f'{e.qual:g}'
        vid = f'{source}_{n}'
        if e.kind == 'INS':
            info = [f'SVTYPE=INS', f'SVLEN={e.length}', f'END={e.p1}'] + common
            lines.append((e.c1, e.p1, [e.c1, str(e.p1), vid, base(e.c1, e.p1), '<INS>', q, 'PASS', ';'.join(info), 'GT', '0/1']))
        elif e.kind == 'INT':
            svlen = -e.length if e.svtype == 'DEL' else e.length
            info = [f'SVTYPE={e.svtype}', f'SVLEN={svlen}', f'END={e.p2}'] + common
            lines.append((e.c1, e.p1, [e.c1, str(e.p1), vid, base(e.c1, e.p1), f'<{e.svtype}>', q, 'PASS', ';'.join(info), 'GT', '0/1']))
        else:
            id1, id2 = vid + '_1', vid + '_2'
            b1, b2 = base(e.c1, e.p1), base(e.c2, e.p2)
            info1 = ['SVTYPE=BND', f'MATEID={id2}', f'CHR2={e.c2}', f'END2={e.p2}'] + common
            info2 = ['SVTYPE=BND', f'MATEID={id1}', f'CHR2={e.c1}', f'END2={e.p1}'] + common
            lines.append((e.c1, e.p1, [e.c1, str(e.p1), id1, b1, bnd_alt(b1, e.s1, e.s2, e.c2, e.p2), q, 'PASS', ';'.join(info1), 'GT', '0/1']))
            lines.append((e.c2, e.p2, [e.c2, str(e.p2), id2, b2, bnd_alt(b2, e.s2, e.s1, e.c1, e.p1), q, 'PASS', ';'.join(info2), 'GT', '0/1']))
    lines.sort(key=lambda x: (CHR_ORDER[x[0]], x[1]))

    header = ['##fileformat=VCFv4.2',
              f'##source=sv_to_std_vcf.py ({source}; mode={args.mode}; trl={args.trl}; '
              f'min_len={args.min_len}; pass_only={args.pass_only}; window={args.window}; '
              f'absorb_linked={args.absorb_linked})',
              f'##reference={os.path.abspath(args.ref)}']
    for c in CANONICAL:
        header.append(f'##contig=<ID={c},length={fa.get_reference_length(c)}>')
    header += [
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=DUP,Description="Duplication (tandem)">',
        '##ALT=<ID=INV,Description="Inversion">',
        '##ALT=<ID=INS,Description="Insertion (sequence not reported)">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant (negative for DEL)">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
        '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">',
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome of mate breakend">',
        '##INFO=<ID=END2,Number=1,Type=Integer,Description="Position of mate breakend">',
        '##INFO=<ID=NMERGED,Number=1,Type=Integer,Description="Number of original records collapsed into this one">',
        '##INFO=<ID=ORIG_TYPES,Number=.,Type=String,Description="Original types of collapsed records (type:count)">',
        '##INFO=<ID=ORIG_IDS,Number=.,Type=String,Description="Original IDs of collapsed records (representative first)">',
        '##INFO=<ID=SR,Number=1,Type=Integer,Description="Sum of supporting reads over collapsed records">',
        '##INFO=<ID=ORIG_SVLEN,Number=1,Type=Integer,Description="Original SVLEN of representative record">',
        '##INFO=<ID=TDRN,Number=1,Type=Integer,Description="eLaRodON number of tandem repeats">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (placeholder)">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + args.sample,
    ]
    out = args.output
    plain = out[:-3] if out.endswith('.gz') else out
    with open(plain, 'w') as f:
        f.write('\n'.join(header) + '\n')
        for _, _, col in lines:
            f.write('\t'.join(col) + '\n')
    if out.endswith('.gz'):
        subprocess.run(['bgzip', '-f', plain], check=True)
        subprocess.run(['tabix', '-f', '-p', 'vcf', out], check=True)
    return len(lines)


def write_report(rep, events_before, events_after, n_lines, path):
    with open(path, 'w') as f:
        f.write('step\ttype\tdetail\tcount\n')
        for (step, t, d), v in sorted(rep.items()):
            f.write(f'{step}\t{t}\t{d}\t{v}\n')
        for label, evs in (('events_before_collapse', events_before), ('events_after_collapse', events_after)):
            c = collections.Counter(e.svtype for e in evs)
            for t, v in sorted(c.items()):
                f.write(f'{label}\t{t}\t\t{v}\n')
        f.write(f'vcf_records\tall\t\t{n_lines}\n')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('caller', choices=['elarodon', 'svaba'])
    p.add_argument('-i', '--input', nargs='+', required=True, help='input VCF(s); for svaba sv and/or indel VCF')
    p.add_argument('-o', '--output', required=True, help='output VCF (.vcf.gz -> bgzip + tabix)')
    p.add_argument('-r', '--ref', required=True, help='reference FASTA (UCSC hg19, chr-names) with .fai')
    p.add_argument('--mode', choices=['typed', 'breakend'], default='typed')
    p.add_argument('--trl', choices=['bnd', 'ins', 'ins_intra'], default='bnd',
                   help='eLaRodON TRL in typed mode: two junctions (bnd), insertion (ins), or insertion within '
                        'a chromosome and two junctions between chromosomes (ins_intra; an interchromosomal '
                        'BND linked to an INS is then kept instead of the INS). Default: bnd')
    p.add_argument('--min-len', type=int, default=200)
    p.add_argument('--pass-only', action='store_true')
    p.add_argument('--keep-invtd', action='store_true', help='keep BND_INVTD (INVTD is always dropped)')
    p.add_argument('--no-absorb-linked', dest='absorb_linked', action='store_false')
    p.add_argument('--link-dist', type=int, default=1000)
    p.add_argument('--no-collapse', dest='collapse', action='store_false')
    p.add_argument('--window', type=int, default=200, help='collapse window, bp')
    p.add_argument('--pctsize', type=float, default=0.7, help='INS size ratio for collapse')
    p.add_argument('--sample', default='NA12878')
    p.add_argument('--max-ids', type=int, default=20)
    p.add_argument('--report', help='TSV report (default: <output>.report.tsv)')
    args = p.parse_args()

    rep = collections.Counter()
    if args.caller == 'elarodon':
        recs = []
        for path in args.input:
            recs.extend(read_elarodon(path, args, rep))
        events = []
        for r in recs:
            events.extend(elarodon_events(r, args))
    else:
        events = svaba_events(args.input, args, rep)
    before = list(events)
    if args.collapse:
        events = collapse(events, args, rep)
    n = write_vcf(events, args, args.caller)
    write_report(rep, before, events, n, args.report or args.output + '.report.tsv')
    print(f'{args.caller}: {len(before)} events -> {len(events)} after collapse -> {n} VCF records: {args.output}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
