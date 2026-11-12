# Rebuild of the NA12878 reference set (29.10): copy of length_thresholds/join_all_minlen.py with two fixes and
# the VCF files of the latest tool versions (eLaRodON, NanoVar 1.8.3, Sniffles2 2.7.1 --tandem-repeats).
# Fix 1: deletions of na12878.svaba.indel.vcf got END = POS, so no call could match their second breakpoint;
#        now END = POS + length.
# Fix 2: the second end of an SV of na12878.svaba.sv.vcf was typed by its own ALT, so the right end of a DEL
#        became TD (and of a TD became DEL) and both ends stayed in the set; now the right end gets the type of
#        the left one, and the two ends collapse into one record.
# Fix 3 (29.10 evening): interchromosomal junctions (BND_TRL) have no length and were dropped by the length
#        filter in the reference set and in the calls; now the filter does not apply to any interchromosomal record. Sniffles2 and
#        NanoVar write the mate of an interchromosomal BND only in ALT (CHR2, no END), so the mate was lost;
#        now it is taken from ALT, as for SVIM, NanoSV and Severus. An interchromosomal eLaRodON TRL (fragment
#        c2:[END, END+SVLEN-1] inserted at c1:POS) has two junctions; both are compared.
# Golden set rule is unchanged: a SvABA record >= MIN_LEN is kept if at least one of the six tools found it.
# Usage: python3 make_na12878_reference_set.py 200 results --inputs inputs.json
import sys
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description='Rebuild the SvABA reference set from caller VCF directories.')
parser.add_argument('min_length', type=int)
parser.add_argument('output')
parser.add_argument('--inputs', required=True, help='JSON mapping caller keys to input directories')
args = parser.parse_args()
MIN_LEN = args.min_length
OUT_DIR = args.output
if MIN_LEN < 1:
    parser.error('min_length must be positive')
config_path = Path(args.inputs).resolve()
with config_path.open() as config_file:
    input_directories = json.load(config_file)
required_inputs = ('svaba_path', 'elarodon_path', 'svim_path', 'sniffles_path',
                   'nanovar_path', 'nanosv_path', 'severus_path')
for key in required_inputs:
    if key not in input_directories:
        parser.error('Missing input key: ' + key)
    directory = Path(input_directories[key])
    if not directory.is_absolute():
        directory = config_path.parent / directory
    if not directory.is_dir():
        parser.error('Input directory not found: ' + str(directory))
    input_directories[key] = str(directory.resolve())
import re
import os
from multiprocessing import Pool
from collections import defaultdict
import csv
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict


def remove_non_numeric(input_string):
    return re.sub(r'\D', '', input_string)

def interchrom_mate(alt_field, chrom1):
    """(chrom2, pos2) of a breakend ALT whose mate lies on another chromosome, else None."""
    m = re.search(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]', alt_field)
    if m and m.group(1) != chrom1:
        return m.group(1), int(m.group(2))
    return None


def process_vcf_file(file_path, file_type):
    data = defaultdict(list)
    data_orriginal_lines = {}

    original_lrs = set()

    with open(file_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            info_field = columns[7]
            info_dict = dict(item.split('=') for item in info_field.split(';') if '=' in item)
            alt_field = columns[4]
            ref_field = columns[3]

            try:
               chrom1 = 'chr' + str(int(columns[0]))
            except ValueError:
                if columns[0] in ('X', 'Y'):
                   chrom1 = 'chr' + columns[0]
                else:
                   chrom1 = columns[0]
            chrom2 = chrom1
            pos1 = int(columns[1])
            svtype, svlen, pos2 = '', 0, pos1
            second_junction = None      # interchromosomal eLaRodON TRL: (chrom, pos) of the fragment's far end

            if file_type == 'svaba':
                if 'svaba.sv' in file_path:
                    pos2 = int(remove_non_numeric(columns[4].split(':')[1]))
                    if alt_field[0] == ']':
                        svtype = 'TD'
                        chrom2 = alt_field[1:alt_field.find(':')]
                        try:
                            chrom2 = 'chr' + str(int(chrom2))
                        except ValueError:
                            if chrom2 in ('X', 'Y'):
                                chrom2 = 'chr' + chrom2
                            else:
                                chrom2 = chrom2
                    elif alt_field[-1] == '[':
                        svtype = 'DEL'
                        chrom2 = alt_field[alt_field.find('[') + 1:alt_field.find(':')]
                        try:
                            chrom2 = 'chr' + str(int(chrom2))
                        except ValueError:
                            if chrom2 in ('X', 'Y'):
                                chrom2 = 'chr' + chrom2
                            else:
                                chrom2 = chrom2
                    elif alt_field[0] == '[' or alt_field[-1] == ']':
                        svtype = 'BND_INV'
                        if alt_field[0] == '[':
                            chrom2 = alt_field[1:alt_field.find(':')]
                            try:
                                chrom2 = 'chr' + str(int(chrom2))
                            except ValueError:
                                if chrom2 in ('X', 'Y'):
                                    chrom2 = 'chr' + chrom2
                                else:
                                    chrom2 = chrom2
                        elif alt_field[-1] == ']':
                            chrom2 = alt_field[alt_field.find(']') + 1:alt_field.find(':')]
                            try:
                                chrom2 = 'chr' + str(int(chrom2))
                            except ValueError:
                                if chrom2 in ('X', 'Y'):
                                    chrom2 = 'chr' + chrom2
                                else:
                                    chrom2 = chrom2
                    if chrom2 != chrom1:
                        svtype = 'BND_TRL'
                        svlen = 0
                    else:
                        svlen = pos2 - pos1 + 1
                        if pos2 < pos1 and svtype in ('DEL', 'TD'):       # fix 2: right end of the SV
                            svtype = {'DEL': 'TD', 'TD': 'DEL'}[svtype]
                elif 'svaba.indel' in file_path:
                    ins_del_len = len(alt_field) - len(ref_field)
                    svtype = 'INS' if ins_del_len >= 0 else 'DEL'
                    svlen = int(info_dict['SPAN'])
                    pos2 = pos1 + abs(svlen) if svtype == 'DEL' else pos1      # fix 1
                    svlen = abs(svlen)

            elif file_type == 'elarodon':
                try:
                    svlen = int(info_dict['SVLEN'])
                except KeyError:
                    svlen = 0
                svtype = info_dict['SVTYPE'][1:-1]
                if svtype == 'INS':
                    pos2 = pos1
                else:
                    try:
                        chrom2 = info_dict['CHROM2']
                    except KeyError:
                        print(info_dict)
                        exit()
                    pos2 = int(info_dict['END'])
                    if svtype == 'TRL' and chrom2 != chrom1:
                        second_junction = (chrom2, pos2 + abs(svlen) - 1)

            elif file_type == 'svim':
                svtype = columns[2].split('.')[1]
                if svtype != 'BND':
                    pos2 = int(info_dict['END'])
                    try:
                        svlen = abs(int(info_dict['SVLEN']))
                    except KeyError:
                        svlen = pos2 - pos1 + 1
                else:
                    pos2 = int(remove_non_numeric(columns[4].split(':')[1]))
                    if alt_field[0] == ']':
                        svtype = 'TD'
                        chrom2 = alt_field[1:alt_field.find(':')]
                    elif alt_field[-1] == '[':
                        svtype = 'DEL'
                        chrom2 = alt_field[alt_field.find('[') + 1:alt_field.find(':')]
                    elif alt_field[0] == '[' or alt_field[-1] == ']':
                        svtype = 'BND_INV'
                        if alt_field[0] == '[':
                            chrom2 = alt_field[1:alt_field.find(':')]
                        elif alt_field[-1] == ']':
                            chrom2 = alt_field[alt_field.find(']') + 1:alt_field.find(':')]
                    if chrom2 !=chrom1 :
                        svtype = 'BND_TRL'
                        svlen = 0
                    else:
                        svlen = pos2 - pos1
            
            elif file_type == 'sniffles':
                chrom1 = columns[0]
                pos1 = int(columns[1])
                chrom2 = info_dict.get('CHROM2', chrom1)
                pos2 = int(info_dict.get('END', pos1))
                svlen = abs(int((info_dict.get('SVLEN', 0))))
                svtype = info_dict.get('SVTYPE')
                mate = interchrom_mate(alt_field, chrom1)       # fix 3
                if svtype == 'BND' and mate:
                    (chrom2, pos2), svtype, svlen = mate, 'BND_TRL', 0
            
            elif file_type == 'nanovar':
                chrom1 = columns[0]
                pos1 = int(columns[1])
                chrom2 = info_dict.get('CHROM2', chrom1)
                pos2 = int(info_dict.get('END', pos1))
                svlen = abs(int(remove_non_numeric(info_dict.get('SVLEN', 0))))
                svtype = info_dict.get('SVTYPE')
                mate = interchrom_mate(alt_field, chrom1)       # fix 3
                if svtype == 'BND' and mate:
                    (chrom2, pos2), svtype, svlen = mate, 'BND_TRL', 0
            
            elif file_type == 'nanosv':
                chrom1 = columns[0]
                pos1 = int(columns[1])
                alt_field = columns[4]
                if '<' not in alt_field:
                    pos2 = int(remove_non_numeric(columns[4].split(':')[1]))
                    if alt_field[0] == ']':
                        svtype = 'TD'
                        chrom2 = alt_field[1:alt_field.find(':')]
                    elif alt_field[-1] == '[':
                        svtype = 'DEL'
                        chrom2 = alt_field[alt_field.find('[') + 1:alt_field.find(':')]
                    elif alt_field[0] == '[' or alt_field[-1] == ']':
                        svtype = 'BND_INV'
                        if alt_field[0] == '[':
                            chrom2 = alt_field[1:alt_field.find(':')]
                        elif alt_field[-1] == ']':
                            chrom2 = alt_field[alt_field.find(']') + 1:alt_field.find(':')]
                    if chrom2 != chrom1:
                        svtype = 'BND_TRL'
                        svlen = 0
                    else:
                        svlen = abs(int((info_dict.get('SVLEN', 0))))
                else:
                    svtype = alt_field[1:-1]
                    svlen = abs(int((info_dict.get('SVLEN', 0))))
                    chrom2 = info_dict.get('CHROM2', chrom1)
                    pos2 = int(info_dict.get('END', pos1))
                if svtype in ('DUP_TANDEM', 'DUP_INT'):
                    svtype = 'TD'
            
            elif file_type == 'severus':
                chrom1 = columns[0]
                pos1 = int(columns[1])
                chrom2 = info_dict.get('CHROM2', chrom1)
                pos2 = int(info_dict.get('END', pos1))
                svlen = None
                if 'SVLEN' in info_dict and info_dict.get('SVLEN') is not None:
                    try:
                        svlen = abs(int(remove_non_numeric(info_dict.get('SVLEN'))))
                    except Exception:
                        svlen = abs(pos2 - pos1)
                else:
                    svlen = abs(pos2 - pos1)

                svtype = info_dict.get('SVTYPE')
                if svtype == 'BND':
                    mate_pos = None
                    mate_chrom = chrom1
                    m = re.search(r'[\[\]]([^:\]\[]+):(\d+)[\[\]]', alt_field)
                    if m:
                        mate_chrom = m.group(1)
                        mate_pos = int(m.group(2))
                    # если удалось вытащить координату/хром, используем её
                    if mate_pos is not None:
                        pos2 = mate_pos
                        chrom2 = mate_chrom
                        svlen = abs(pos2 - pos1)
                    else:
                        # иначе считаем как бессмысленную BND — присваиваем нулевую длину
                        svlen = 0
                    # уточняем тип по форме ALT (по аналогии с вашими severus-правилами)
                    if alt_field and alt_field[0] == ']':
                        svtype = 'TD'
                    elif alt_field and alt_field[-1] == '[':
                        svtype = 'DEL'
                    elif alt_field and (alt_field[0] == '[' or alt_field[-1] == ']'):
                        svtype = 'BND_INV'
                    # транс-хромосомные (если хромы разные) — помечаем как BND_TRL и длину = 0
                    if chrom2 != chrom1:
                        svtype = 'BND_TRL'
                        svlen = 0
                # нормализуем DUP -> TD
                if svtype == 'DUP':
                    svtype = 'TD'
            
            if second_junction:
                c_far, p_far = second_junction
                c_ins, p_ins = chrom1, pos1
            if pos2 < pos1:
                chrom1 , chrom2 = chrom2, chrom1 
                pos1, pos2 = pos2, pos1                
            
            svlen = abs(svlen)

            if (svlen < MIN_LEN and chrom1 == chrom2) or 'INVTD' in svtype:     # fix 3: no length filter between chromosomes
                continue

            if second_junction:
                # second junction of an interchromosomal TRL, ordered by position like all records
                a, b = (c_ins, p_ins), (c_far, p_far)
                if b[1] < a[1]:
                    a, b = b, a
                data[(a[0], b[0])].append((a[0], b[0], a[1], b[1], svlen, svtype, alt_field))
            
            if (chrom1, chrom2, pos1, pos2, svtype) not in original_lrs and (chrom2, chrom1, pos2, pos1, svtype) not in original_lrs:
                data_orriginal_lines[(chrom1, chrom2, pos1, pos2, svtype)] = line
                data_orriginal_lines[(chrom2, chrom1, pos2, pos1, svtype)] = line

            # Ensure the entry has the correct format
            entry = (chrom1, chrom2, pos1, pos2, svlen, svtype, alt_field)
            
            if (chrom1, chrom2, pos1, pos2, svtype) not in original_lrs and (chrom2, chrom1, pos2, pos1, svtype) not in original_lrs:
                original_lrs.add((chrom1, chrom2, pos1, pos2, svtype))
                original_lrs.add((chrom2, chrom1, pos2, pos1, svtype))
            else:
                continue

            if len(entry) != 7:
                print(f"Unexpected entry format: {entry}")
                continue

            data[(chrom1, chrom2)].append(entry)

    return data, data_orriginal_lines

def process_folder(folder_path, file_type):
    file_paths = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.vcf')]
    with Pool(processes=12) as pool:
        results = pool.starmap(process_vcf_file, [(fp, file_type) for fp in file_paths])

    combined_data = defaultdict(list)
    all_original_lines = {}
    for result, all_original_lines_temp in results:
        for key, value in result.items():
            combined_data[key].extend(value)
        all_original_lines.update(all_original_lines_temp)

    return combined_data, all_original_lines

def find_matches(svaba_batch, elarodon_data, svim_data, sniffles_data, nanovar_data, nanosv_data, severus_data):
    matching_rows = []
    matched_svaba_lines = []
    unmatched_svaba_lines = []

    elarodon_matches_big = []
    svim_matches_big = []
    nanovar_matches_big = []
    sniffles_matches_big = []
    nanosv_matches_big = []
    severus_matches_big = []

    for chr_pair, list_svaba in svaba_batch.items():
        for svaba_entry in list_svaba:
            if len(svaba_entry) != 7:
                print(f"Unexpected entry format: {svaba_entry}")
                continue

            chr_key = (svaba_entry[0], svaba_entry[1])
            pos1, pos2, svlen, svtype, alt = svaba_entry[2:]
            elarodon_matches = []
            svim_matches = []
            sniffles_matches = []
            nanovar_matches = []
            nanosv_matches = []
            severus_matches = []

            err = round(svlen * 0.1) if svlen > 0 else 100

            if svtype == 'INS':
                elarodon_matches = [entry for entry in elarodon_data.get(chr_key, []) if entry[5] == svtype and abs(entry[2] - pos1) <= err and abs(abs(entry[4]) - abs(svlen)) <= err]
                # if elarodon_matches == []:
                svim_matches = [entry for entry in svim_data.get(chr_key, []) if entry[5] == svtype and abs(entry[2] - pos1) <= err and abs(abs(entry[4]) - abs(svlen)) <= err]
                sniffles_matches = [entry for entry in sniffles_data.get(chr_key, []) if entry[5] == svtype and abs(entry[2] - pos1) <= err and abs(abs(entry[4]) - abs(svlen)) <= err]
                nanovar_matches = [entry for entry in nanovar_data.get(chr_key, []) if entry[5] == svtype and abs(entry[2] - pos1) <= err and abs(abs(entry[4]) - abs(svlen)) <= err]
                nanosv_matches = [entry for entry in nanosv_data.get(chr_key, []) if entry[5] == svtype and abs(entry[2] - pos1) <= err and abs(abs(entry[4]) - abs(svlen)) <= err]
                severus_matches = [entry for entry in severus_data.get(chr_key, []) if entry[5] == svtype and abs(entry[2] - pos1) <= err and abs(abs(entry[4]) - abs(svlen)) <= err]
            else:
                elarodon_matches = [entry for entry in elarodon_data.get(chr_key, []) if (entry[5] in svtype or svtype in entry[5]) and abs(entry[2] - pos1) <= err and abs(entry[3] - pos2) <= err]
                # if elarodon_matches == []:
                svim_matches = [entry for entry in svim_data.get(chr_key, []) if (entry[5] in svtype or svtype in entry[5]) and abs(entry[2] - pos1) <= err and abs(entry[3] - pos2) <= err]
                sniffles_matches = [entry for entry in sniffles_data.get(chr_key, []) if (entry[5] in svtype or svtype in entry[5]) and abs(entry[2] - pos1) <= err and abs(entry[3] - pos2) <= err]
                nanovar_matches = [entry for entry in nanovar_data.get(chr_key, []) if (entry[5] in svtype or svtype in entry[5]) and abs(entry[2] - pos1) <= err and abs(entry[3] - pos2) <= err]
                nanosv_matches = [entry for entry in nanosv_data.get(chr_key, []) if (entry[5] in svtype or svtype in entry[5]) and abs(entry[2] - pos1) <= err and abs(entry[3] - pos2) <= err]
                severus_matches = [entry for entry in severus_data.get(chr_key, []) if (entry[5] in svtype or svtype in entry[5]) and abs(entry[2] - pos1) <= err and abs(entry[3] - pos2) <= err]
            # else:
            #     elarodon_matches = [entry for entry in elarodon_data.get(chr_key, []) if abs(entry[2] - pos1) <= 200 and abs(entry[3] - pos2) <= 200]
            #     # if elarodon_matches == []:
            #     svim_matches = [entry for entry in svim_data.get(chr_key, []) if abs(entry[2] - pos1) <= 200 and abs(entry[3] - pos2) <= 200]

            match_found = False
            matched_row_elarodon = None
            matched_row_svim = None
            matched_row_sniffles = None
            matched_row_nanovar = None
            matched_row_nanosv = None
            matched_row_severus = None

            if elarodon_matches:
                match_found = True
                matched_row_elarodon = elarodon_matches[0]
                elarodon_matches_big.extend(elarodon_matches[:1])

            if svim_matches: # not match_found and
                match_found = True
                matched_row_svim = svim_matches[0]
                svim_matches_big.extend(svim_matches[:1])
            
            if sniffles_matches: # not match_found and
                match_found = True
                matched_row_sniffles = sniffles_matches[0]
                sniffles_matches_big.extend(sniffles_matches[:1])
            
            if nanovar_matches: # not match_found and
                match_found = True
                matched_row_nanovar = nanovar_matches[0]
                nanovar_matches_big.extend(nanovar_matches[:1])
            
            if nanosv_matches: # not match_found and
                match_found = True
                matched_row_nanosv = nanosv_matches[0]
                nanosv_matches_big.extend(nanosv_matches[:1])
            
            if severus_matches: # not match_found and
                match_found = True
                matched_row_severus = severus_matches[0]
                severus_matches_big.extend(severus_matches[:1])

            if match_found:
                new_row = {
                    'CHR1': chr_key[0],
                    'CHR2': chr_key[1],
                    'POS1': pos1,
                    'POS2': pos2,
                    'SVLEN': svlen,
                    'SVTYPE': svtype,
                    'ALT': alt
                }
                if matched_row_elarodon:
                    new_row.update({f'{k}_elarodon': v for k, v in zip(['CHR1', 'CHR2', 'POS1', 'POS2', 'SVLEN', 'SVTYPE', 'ALT'], matched_row_elarodon)})
                if matched_row_svim:
                    new_row.update({f'{k}_svim': v for k, v in zip(['CHR1', 'CHR2', 'POS1', 'POS2', 'SVLEN', 'SVTYPE', 'ALT'], matched_row_svim)})
                if matched_row_sniffles:
                    new_row.update({f'{k}_sniffles': v for k, v in zip(['CHR1', 'CHR2', 'POS1', 'POS2', 'SVLEN', 'SVTYPE', 'ALT'], matched_row_sniffles)})
                if matched_row_nanovar:
                    new_row.update({f'{k}_nanovar': v for k, v in zip(['CHR1', 'CHR2', 'POS1', 'POS2', 'SVLEN', 'SVTYPE', 'ALT'], matched_row_nanovar)})
                if matched_row_nanosv:
                    new_row.update({f'{k}_nanosv': v for k, v in zip(['CHR1', 'CHR2', 'POS1', 'POS2', 'SVLEN', 'SVTYPE', 'ALT'], matched_row_nanosv)})
                if matched_row_severus:
                    new_row.update({f'{k}_severus': v for k, v in zip(['CHR1', 'CHR2', 'POS1', 'POS2', 'SVLEN', 'SVTYPE', 'ALT'], matched_row_severus)})
                matching_rows.append(new_row)
                try:
                    matched_svaba_lines.append(svaba_original_lines[(chr_key[0], chr_key[1], pos1, pos2, svtype)])
                except KeyError:
                    pass
            else:
                try:
                    unmatched_svaba_lines.append(svaba_original_lines[(chr_key[0], chr_key[1], pos1, pos2, svtype)])
                except KeyError:
                    pass

    return matching_rows, matched_svaba_lines, unmatched_svaba_lines, elarodon_matches_big, svim_matches_big, sniffles_matches_big, nanovar_matches_big, nanosv_matches_big, severus_matches_big

def compare_data(svaba_data, elarodon_data, svim_data, sniffles_data, nanovar_data, nanosv_data, severus_data, original_lines, batch_size=1000):
    matching_rows = []
    matched_svaba_lines = []
    unmatched_svaba_lines = []
    # total_entries = sum(len(entries) for entries in svaba_data.values())
    original_elarodon = set()
    original_svim = set()
    original_sniffles = set()
    original_nanovar = set()
    original_nanosv = set()
    original_severus = set()

    # batches = []
    # current_batch = []

    # for entries in svaba_data.values():
    #     for entry in entries:
    #         current_batch.append(entry)
    #         if len(current_batch) >= batch_size:
    #             batches.append((current_batch,))
    #             current_batch = []
    # if current_batch:
    #     batches.append((current_batch, ))

    # with Pool(processes=12) as pool:
    #     for index, (batch_result, matched_svaba, unmatched_svaba, elarodon_match, svim_match) in enumerate(pool.starmap(find_matches, [(batch, elarodon_data, svim_data) for batch in batches])):

    batch_result, matched_svaba_lines, unmatched_svaba, elarodon_match, svim_match, sniffles_match, nanovar_match, nanosv_match, severus_match = find_matches(svaba_data, elarodon_data, svim_data, sniffles_data, nanovar_data, nanosv_data, severus_data)     
    matching_rows.extend(batch_result)
    # print(matched_svaba)
    unmatched_svaba_lines.extend(unmatched_svaba)
    # print('elarodon_match', elarodon_match)
    for m in elarodon_match:
        original_elarodon.add(m)
    for m in svim_match:
        original_svim.add(m)
    for m in sniffles_match:
        original_sniffles.add(m)
    for m in nanovar_match:
        original_nanovar.add(m)
    for m in nanosv_match:
        original_nanosv.add(m)
    for m in severus_match:
        original_severus.add(m)
    
    print(f'\nTotal number of matches in original VCF eLaRodON lines: {len(original_elarodon)}')
    print(f'Total number of matches in original VCF SVIM lines: {len(original_svim)}')
    print(f'Total number of matches in original VCF SNIFFLES lines: {len(original_sniffles)}')
    print(f'Total number of matches in original VCF NANOVAR lines: {len(original_nanovar)}')
    print(f'Total number of matches in original VCF NANOSV lines: {len(original_nanosv)}')
    print(f'Total number of matches in original VCF SEVERUS lines: {len(original_severus)}')
    # input()

    return matching_rows, matched_svaba_lines, unmatched_svaba_lines

def save_to_csv(data, file_path):
    if isinstance(data, dict):
        rows = []
        for (chr1, chr2), entries in data.items():
            for entry in entries:
                row = {
                    'CHR1': chr1,
                    'CHR2': chr2,
                    'POS1': entry[2],
                    'POS2': entry[3],
                    'SVLEN': entry[4],
                    'SVTYPE': entry[5],
                    'ALT': entry[6]
                }
                rows.append(row)
    else:
        rows = data

    if rows:
        fieldnames = sorted(set(key for row in rows for key in row.keys()), key=lambda x: (x.endswith('svaba'), x.endswith('elarodon'), x.endswith('svim'), x))

        with open(file_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

def save_matched_svaba_to_vcf(matched_svaba_lines, output_path):
    with open(output_path, 'w') as output_vcf:
        # Write matched entries
        for line in matched_svaba_lines:
            # print(line)
            # input()
            output_vcf.write(f"{line}")

def get_length_category(svlen):
    # if svlen <= 50:
    #     return "0-50"
    # elif 50 < svlen <= 200:
    #     return "50-200"
    if 200 <= svlen <= 1000:
        return "200-1000"
    elif 1000 < svlen <= 5000:
        return "1000-5000"
    elif 5000 < svlen <= 10000:
        return "5000-10000"
    else:
        return ">10000"
    
def analyze_unmatched_stats(unmatched_lines):
    type_stats = defaultdict(int)
    len_categories = ["200-1000", "1000-5000", "5000-10000", ">10000"]
    len_stats = {category: 0 for category in len_categories}
    type_len_stats = defaultdict(lambda: {category: 0 for category in len_categories})
    
    for line in unmatched_lines:
        if not line.startswith('#'):
            columns = line.strip().split('\t')
            info_field = columns[7]
            info_dict = dict(item.split('=') for item in info_field.split(';') if '=' in item)
            alt_field = columns[4]
            ref_field = columns[3]
            
            # Определяем тип перестройки (аналогично process_vcf_file)
            svtype = ''
            svlen = 0

            pos1 = int(columns[1])
            
            # Для SVABA
            if 'DBSNP' not in info_field and 'LOD' not in info_field:
                pos2 = int(remove_non_numeric(columns[4].split(':')[1]))
                if alt_field[0] == ']':
                    svtype = 'TD'
                elif alt_field[-1] == '[':
                    svtype = 'DEL'
                elif alt_field[0] == '[' or alt_field[-1] == ']':
                    svtype = 'BND_INV'
                svlen = pos2 - int(columns[1]) + 1
            else:
                ins_del_len = len(alt_field) - len(ref_field)
                svtype = 'INS' if ins_del_len >= 0 else 'DEL'
                svlen = int(info_dict['SPAN'])
                pos2 = pos1 - svlen
                if pos2 < pos1:
                    pos2 = pos1
            
            if pos2 < pos1:
                pos1, pos2 = pos2, pos1
            
            svlen = abs(svlen)

            if svlen < MIN_LEN or 'INVTD' in svtype:
                continue
            
            # Определяем категорию длины
            len_category = get_length_category(svlen)
            
            # Обновляем статистику
            type_stats[svtype] += 1
            len_stats[len_category] += 1
            type_len_stats[svtype][len_category] += 1
    
    return type_stats, len_stats, type_len_stats

def print_unmatched_stats(type_stats, len_stats, type_len_stats):
    print("\n=== Statistics for unmatched SVABA variants ===")
    
    # Общая статистика по типам
    print("\nBy SV type (total counts):")
    for svtype, count in sorted(type_stats.items()):
        print(f"{svtype}: {count} variants ({count/sum(type_stats.values()):.1%})")
    
    # Статистика по категориям длин
    print("\nBy length categories:")
    total = sum(len_stats.values())
    for category in ["200-1000", "1000-5000", "5000-10000", ">10000"]:
        count = len_stats[category]
        if count > 0:
            print(f"{category.ljust(8)} bp: {count} variants ({count/total:.1%})")
    
    # Детальная статистика по типам и длинам
    print("\nDetailed statistics by SV type and length:")
    for svtype in sorted(type_len_stats.keys()):
        print(f"\n{svtype}:")
        total = type_stats[svtype]
        for category in ["200-1000", "1000-5000", "5000-10000", ">10000"]:
            count = type_len_stats[svtype][category]
            if count > 0:
                print(f"  {category.ljust(8)} bp: {count} variants ({count/total:.1%})")
    
    print("\n=== End of statistics ===")

def print_histogram(data, title, max_width=50):
    if not data:
        return
    
    print(f"\n{title}")
    max_value = max(data.values())
    
    for key, value in sorted(data.items()):
        bar_length = int((value / max_value) * max_width) if max_value > 0 else 0
        bar = '#' * bar_length
        print(f"{str(key).ljust(10)}: {bar} {value}")

def filter_svim_only_matches(compared_results):
    """Filter results to find Svaba variants matched only in SVIM but not in eLaRodON"""
    svim_only_matches = []
    
    for row in compared_results:
        # Проверяем, что есть совпадение с SVIM, но нет с eLaRodON
        has_svim_match = any(key.endswith('_svim') for key in row.keys())
        has_elarodon_match = any(key.endswith('_elarodon') for key in row.keys())
        
        if has_svim_match and not has_elarodon_match:
            # Создаем новую строку только с нужными полями (без _elarodon)
            new_row = {k: v for k, v in row.items() if not k.endswith('_elarodon')}
            svim_only_matches.append(new_row)
    
    return svim_only_matches

def filter_elarodon_only_matches(compared_results):
    """Filter results to find Svaba variants matched only in eLaRodON but not in SVIM"""
    elarodon_only_matches = []
    
    for row in compared_results:
        has_svim_match = any(key.endswith('_svim') for key in row.keys())
        has_elarodon_match = any(key.endswith('_elarodon') for key in row.keys())
        
        if has_elarodon_match and not has_svim_match:
            # Создаем новую строку только с нужными полями (без _svim)
            new_row = {k: v for k, v in row.items() if not k.endswith('_svim')}
            elarodon_only_matches.append(new_row)
    
    return elarodon_only_matches


os.makedirs(OUT_DIR, exist_ok=True)
vcfs = input_directories

global_svaba_variants, svaba_original_lines = process_folder(vcfs['svaba_path'], 'svaba')
global_elarodon_variants, _ = process_folder(vcfs['elarodon_path'], 'elarodon')
global_svim_variants, _ = process_folder(vcfs['svim_path'], 'svim')
global_sniffles_variants, _ = process_folder(vcfs['sniffles_path'], 'sniffles')
global_nanovar_variants, _ = process_folder(vcfs['nanovar_path'], 'nanovar')
global_nanosv_variants, _ = process_folder(vcfs['nanosv_path'], 'nanosv')
global_severus_variants, _ = process_folder(vcfs['severus_path'], 'severus')

save_to_csv(global_svaba_variants, OUT_DIR + '/global_svaba_variants.csv')
print('SVABA saved...')
save_to_csv(global_elarodon_variants, OUT_DIR + '/global_elarodon_variants.csv')
print('ELARODON saved...')
save_to_csv(global_svim_variants, OUT_DIR + '/global_svim_variants.csv')
print('SVIM saved...')
save_to_csv(global_sniffles_variants, OUT_DIR + '/global_sniffles_variants.csv')
print('SNIFFLES saved...')
save_to_csv(global_nanovar_variants, OUT_DIR + '/global_nanovar_variants.csv')
print('NANOVAR saved...')
save_to_csv(global_nanosv_variants, OUT_DIR + '/global_nanosv_variants.csv')
print('NANOSV saved...')
save_to_csv(global_severus_variants, OUT_DIR + '/global_severus_variants.csv')
print('NANOSV saved...')

print('Start finding SVABA...')
compared_results, matched_svaba_lines, unmatched_svaba_lines = compare_data(
    global_svaba_variants, 
    global_elarodon_variants, 
    global_svim_variants, 
    global_sniffles_variants,
    global_nanovar_variants,
    global_nanosv_variants,
    global_severus_variants,
    svaba_original_lines
)

print(f"\nTotal SVABA variants: {len(svaba_original_lines)}")
print(f"\nMatched variants: {len(matched_svaba_lines)}")
print(f"Unmatched variants: {len(unmatched_svaba_lines)}")

if unmatched_svaba_lines:
    type_stats, len_stats, type_len_stats = analyze_unmatched_stats(unmatched_svaba_lines)
    print_unmatched_stats(type_stats, len_stats, type_len_stats)
else:
    print("No unmatched variants found")

print('Comparison was finished! Saving...')
save_to_csv(compared_results, OUT_DIR + '/compared_svaba_variants.csv')

# Save matched SVABA entries to a new VCF file
save_matched_svaba_to_vcf(matched_svaba_lines, OUT_DIR + '/matched_svaba_variants.vcf')
print('Matched SVABA entries saved to VCF...')

# Save unmatched SVABA entries to a separate VCF file
save_matched_svaba_to_vcf(unmatched_svaba_lines, OUT_DIR + '/unmatched_svaba_variants.vcf')
print('Unmatched SVABA entries saved to VCF...')



# svim_only_matches = filter_svim_only_matches(compared_results)
# print(f"\nFound {len(svim_only_matches)} variants matched in Svaba and SVIM but not in eLaRodON")

# if svim_only_matches:
#     svim_only_file = OUT_DIR + '/svaba_svim_only_matches.csv'
#     save_to_csv(svim_only_matches, svim_only_file)
#     print(f"Svaba-SVIM only matches saved to {svim_only_file}")




# elarodon_only_matches = filter_elarodon_only_matches(compared_results)
# print(f"\nFound {len(elarodon_only_matches)} variants matched in Svaba and eLaRodON but not in SVIM")

# # Сохраняем их в отдельный файл
# if elarodon_only_matches:
#     elarodon_only_file = OUT_DIR + '/svaba_elarodon_only_matches.csv'
#     save_to_csv(elarodon_only_matches, elarodon_only_file)
#     print(f"Svaba-eLaRodON only matches saved to {elarodon_only_file}")


# print_histogram(len_stats, "Length categories distribution")
# for svtype in sorted(type_len_stats.keys()):
#     print_histogram(type_len_stats[svtype], f"Length distribution for {svtype}")