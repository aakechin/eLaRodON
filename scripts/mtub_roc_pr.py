"""
Figure 5A (M. tuberculosis, 38 genomes): copy of mtub/sens_compare/quality_ROC_all.py
with precision-recall curves added.
Sniffles2 mosaic (--mosaic, Sniffles2 2.6.0) is added and parsed like Sniffles2.
Fixes (off with --as-published): DUP is renamed to TD for Sniffles2, NanoVar, SVIM and NanoSV.
With --as-published the ROC AUC equals mtub/sens_compare/plots_ROC/metrics_table_mtub.csv.

With --lowsupport the low-support runs (rel7/New_files/mtub_lowsupport) are added as separate programs.
With fixes on, a VCF whose name has another sample ID is skipped (svim/ERR9030471 also holds svim_ERR9040471.vcf).

Usage: python3 mtub_roc_pr.py <output_dir> [--as-published] [--lowsupport]
       [--vcf-root=DIR] [--mosaic-root=DIR] [--lowsupport-root=DIR]  (other VCF folders, same layout)
"""
import os
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pr_utils

FIXES = '--as-published' not in sys.argv
# Sniffles2 --mosaic runs (not in the published figure)
MOSAIC_ROOT = '/media/regina/MR_GREY/human/rel7/New_files/mtub_sniffles_mosaic'
LOWSUPPORT = '--lowsupport' in sys.argv
LOWSUPPORT_ROOT = '/media/regina/MR_GREY/human/rel7/New_files/mtub_lowsupport'


def path_arg(name, default):
    values = [a.split('=', 1)[1] for a in sys.argv if a.startswith(name + '=')]
    return values[0] if values else default


MOSAIC_ROOT = path_arg('--mosaic-root', MOSAIC_ROOT)
LOWSUPPORT_ROOT = path_arg('--lowsupport-root', LOWSUPPORT_ROOT)
from sklearn.metrics import roc_curve, auc, roc_auc_score
import re
from sklearn.metrics import precision_recall_curve, f1_score
import numpy as np

def remove_non_numeric(input_string):
    numeric_string = re.sub(r'\D', '', input_string)
    return numeric_string

def unmatched_variants_to_csv(data, output):
    with open(output, 'w+') as out_csv:
        for row in data:
            out_csv.write('\t'.join(row)+'\n')

# Автоматическое создание списков образцов для каждого инструмента
def create_samples_list(tool_name, root=None):
    samples = []
    tool_dir = os.path.join(root or vcf_root, tool_name)
    
    if not os.path.exists(tool_dir):
        print(f"Предупреждение: папка инструмента {tool_name} не найдена: {tool_dir}")
        return samples
    
    # Перебираем все подпапки в папке инструмента
    for sample_dir in os.listdir(tool_dir):
        sample_path = os.path.join(tool_dir, sample_dir)
        # print('sample_path', sample_path)
        if not os.path.isdir(sample_path):
            continue
            
        # Ищем соответствующий CSV файл
        sample_id = sample_dir
        print(sample_id)
        csv_path = None
        
        # Ищем файлы CSV, начинающиеся с ID образца
        for fname in os.listdir(csv_root):
            if fname.startswith(sample_id) and fname.endswith('.csv'):
                csv_path = os.path.join(csv_root, fname)
                print('csv_path', csv_path)
                break
        
        if csv_path and os.path.exists(csv_path):
            samples.append({
                'csv_path': csv_path,
                'vcf_folder': sample_path
            })
        else:
            print(f"Предупреждение: для образца {sample_id} не найден CSV файл в {csv_root}")
    print(samples)
    # input()
    return samples
    
def elarodon_process_vcf_for_dataframe(real_variants, vcf_path):
    data = []
    unmatched_real_variants = []
    matched_indexes = []
    real_sv_matched = []
    vcf_rearrangements = 0 
    matches = 0
    mismatches = 0
    # print(real_variants.to_dict())
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            info_field = columns[7]  # INFO field
            chrom = columns[0]
            qual_field = float(columns[5])  # Use the QUAL value directly
            info_dict = dict(item.split('=') for item in info_field.split(';') if '=' in item)
            # Извлекаем SR (количество поддерживающих ридов)
            sr_value = info_dict.get('SR', '0')
            try:
                sr_value = int(sr_value)
            except ValueError:
                sr_value = None

            pos = int(columns[1])
            # if pos == 3048413:
            #     print('STOP')
            #     print(info_dict)
            #     input()
            svlen = abs(int(info_dict.get('SVLEN', 0)))
            pos2 = int(info_dict.get('END', pos))
            svtype = info_dict.get('SVTYPE')[1:-1]

            if svlen < 200 or 'INVTD' in svtype:
                continue
            
            svlen_diff = 200

            if svlen > 0 and svlen < 10000:
                svlen_diff = svlen * 0.1
            elif svlen >= 10000:
                svlen_diff = svlen * 0.05

            found = False
            real_info_list = []
            real_info = ''

            vcf_rearrangements += 1

            # if pos == 3048413 and svtype == 'TRL':
            #     print('STOP')
            #     print(info_dict)
            #     print('pos', pos)
            #     print('pos2', pos2)
            #     print('svlen', svlen)
                # input()
            # c = 0
            for index, real_row in real_variants.iterrows():
                # if pos:
                # print(real_row)
                # c += 1
                # print(c)
                # input()
                real_pos = real_row['POS']
                # print("real_row['SVLEN']", real_row['SVLEN'])
                real_svlen = abs(real_row['SVLEN'])
                
                real_svtype = real_row['SVTYPE'].replace(" ", "")
                
                if 'INVTD' in svtype:
                    # print(real_row)
                    continue

                if 'BND' in svtype:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'INV':
                    if ((svtype in real_svtype or real_svtype in svtype) and abs(real_pos - pos) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or (real_svtype == 'INS_DEL' and abs(real_pos - pos) <= 200):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TD':
                    if (real_svtype == 'INS' and abs(real_pos - pos2) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or ((svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TRL':
                    # if pos == 3048413:
                    #     print(real_row)
                    # if pos == 3048413 and real_pos == 3890779:
                    #     print()
                    #     print(real_row)
                    #     print(info_dict)
                    #     print('pos', pos)
                    #     print('pos2', pos2)
                    #     print('svlen', svlen)

                    if (((svtype in real_svtype or real_svtype in svtype) or real_svtype == 'INS') and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200)) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                else:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
            # if pos == 3048413:
            #     print('c', c)
            if found:
                for real_info in real_info_list:
                    # print(real_info)
                    real_info = real_info.split('_')
                    real_svtype = real_info[0]
                    try:
                        real_pos = int(real_info[1])
                        real_svlen = int(real_info[2])
                    except ValueError:
                        real_svtype = real_info[0]+'_'+real_info[1]
                        real_pos = int(real_info[2])
                        real_svlen = int(real_info[3])
                    if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
                        real_sv_matched.append((real_pos, real_svlen, real_svtype))
                        matches += 1
                real_info = '_OR_'.join(real_info_list)
            else:
                mismatches += 1

            vcf_info = f"{chrom}_{pos}_{pos2}_{svlen}"
            data.append([
                        svtype,              # тип реаранжировки
                        int(found),          # обнаружено?
                        qual_field,          # QUAL
                        vcf_info,            # VCF_Info
                        real_info,           # Real_Info
                        sr_value             # ← новй столбец
                    ])
    for index, real_row in real_variants.iterrows():
        if index not in matched_indexes:
            unmatched_real_variants.append((str(real_row['POS']), str(real_row['SVTYPE']), str(real_row['SVLEN'])))

    for index, real_row in real_variants.iterrows():
        real_pos = real_row['POS']
        real_svlen = abs(real_row['SVLEN'])
        real_svtype = real_row['SVTYPE'].replace(" ", "")
        if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
            real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
            data.append(['ND', 1, -1, '', real_info])
    print()
    print(f'VCF path: {vcf_path}')
    print(f"Total rearrangements in vcf: {vcf_rearrangements}")
    print(f"Matches with real rearrangements: {matches}")
    print(f"Mismatches with real rearrangements: {mismatches}")

    return pd.DataFrame(
                        data,
                        columns=['Rearrangement', 'Present', 'QUAL', 'VCF_Info', 'Real_Info', 'SR']
                        ), unmatched_real_variants

def nanosv_process_vcf_for_dataframe(real_variants, vcf_path):
    data = []
    unmatched_real_variants = []
    matched_indexes = []
    real_sv_matched = []
    vcf_rearrangements = 0 
    matches = 0
    mismatches = 0
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            info_field = columns[7]  # INFO field
            
            try:
                qual_field = float(columns[5])  # Use the QUAL value directly
            except ValueError:
                return None
            
            info_dict = dict(item.split('=') for item in info_field.split(';') if '=' in item)
            
            chrom1 = columns[0]
            pos = int(columns[1])

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
                pos2 = int(info_dict.get('END', pos))

            if FIXES and svtype == 'DUP':
                svtype = 'TD'
            if svlen < 200:
                continue

            svlen_diff = 200

            if svlen > 0 and svlen < 10000:
                svlen_diff = svlen * 0.1
            elif svlen >= 10000:
                svlen_diff = svlen * 0.05

            found = False
            real_info = ""
            real_info_list = []

            vcf_rearrangements += 1

            for index, real_row in real_variants.iterrows():
                real_pos = real_row['POS']
                real_svlen = abs(real_row['SVLEN'])
                real_svtype = real_row['SVTYPE'].replace(" ", "")
                if real_svtype == 'INVTD':
                    # print(real_row)
                    continue

                if 'BND' in svtype:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'INV':
                    if ((svtype in real_svtype or real_svtype in svtype) and abs(real_pos - pos) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or (real_svtype == 'INS_DEL' and abs(real_pos - pos) <= 200):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TD':
                    if (real_svtype == 'INS' and abs(real_pos - pos2) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or ((svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TRL':

                    if (((svtype in real_svtype or real_svtype in svtype) or real_svtype == 'INS') and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200)) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                else:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
            if found:
                for real_info in real_info_list:
                    # print(real_info)
                    real_info = real_info.split('_')
                    real_svtype = real_info[0]
                    try:
                        real_pos = int(real_info[1])
                        real_svlen = int(real_info[2])
                    except ValueError:
                        real_svtype = real_info[0]+'_'+real_info[1]
                        real_pos = int(real_info[2])
                        real_svlen = int(real_info[3])
                    if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
                        real_sv_matched.append((real_pos, real_svlen, real_svtype))
                        matches += 1
                real_info = '_OR_'.join(real_info_list)
            else:
                mismatches += 1

            vcf_info = f"{chrom1}_{pos}_{pos2}_{svlen}"
            data.append([svtype, int(found), qual_field, vcf_info, real_info])
    
    for index, real_row in real_variants.iterrows():
        if index not in matched_indexes:
            unmatched_real_variants.append((str(real_row['POS']), str(real_row['SVTYPE']), str(real_row['SVLEN'])))

    for index, real_row in real_variants.iterrows():
        real_pos = real_row['POS']
        real_svlen = abs(real_row['SVLEN'])
        real_svtype = real_row['SVTYPE'].replace(' ', '')
        if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
            real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
            data.append(['ND', 1, -1, '', real_info])
    print()
    print(f'VCF path: {vcf_path}')
    print(f"Total rearrangements in vcf: {vcf_rearrangements}")
    print(f"Matches with real rearrangements: {matches}")
    print(f"Mismatches with real rearrangements: {mismatches}")

    return pd.DataFrame(data, columns=['Rearrangement', 'Present', 'QUAL', 'VCF_Info', 'Real_Info']), unmatched_real_variants

def nanovar_process_vcf_for_dataframe(real_variants, vcf_path):
    data = []
    unmatched_real_variants = []
    matched_indexes = []
    real_sv_matched = []
    vcf_rearrangements = 0 
    matches = 0
    mismatches = 0
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            info_field = columns[7]  # INFO field
            
            qual_field = float(columns[5])  # Use the QUAL value directly
            info_dict = dict(item.split('=') for item in info_field.split(';') if '=' in item)
            
            chrom = columns[0]
            pos = int(columns[1])
            chrom2 = info_dict.get('CHROM2', chrom)
            pos2 = int(info_dict.get('END', pos))
            svlen = abs(int(remove_non_numeric(info_dict.get('SVLEN', 0))))
            svtype = info_dict.get('SVTYPE')
            if svtype == 'DUP':
                svtype = 'TD' if FIXES else svtype

            if svlen < 200:
                continue

            svlen_diff = 200

            if svlen > 0 and svlen < 10000:
                svlen_diff = svlen * 0.1
            elif svlen >= 10000:
                svlen_diff = svlen * 0.05

            found = False
            real_info = ""
            real_info_list = []

            vcf_rearrangements += 1

            for index, real_row in real_variants.iterrows():
                real_pos = real_row['POS']
                real_svlen = abs(real_row['SVLEN'])
                real_svtype = real_row['SVTYPE'].replace(" ", "")
                if real_svtype == 'INVTD':
                    # print(real_row)
                    continue

                if 'BND' in svtype:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'INV':
                    if ((svtype in real_svtype or real_svtype in svtype) and abs(real_pos - pos) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or (real_svtype == 'INS_DEL' and abs(real_pos - pos) <= 200):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TD':
                    if (real_svtype == 'INS' and abs(real_pos - pos2) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or ((svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TRL':

                    if (((svtype in real_svtype or real_svtype in svtype) or real_svtype == 'INS') and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200)) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                else:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
            if found:
                for real_info in real_info_list:
                    # print(real_info)
                    real_info = real_info.split('_')
                    real_svtype = real_info[0]
                    try:
                        real_pos = int(real_info[1])
                        real_svlen = int(real_info[2])
                    except ValueError:
                        real_svtype = real_info[0]+'_'+real_info[1]
                        real_pos = int(real_info[2])
                        real_svlen = int(real_info[3])
                    if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
                        real_sv_matched.append((real_pos, real_svlen, real_svtype))
                        matches += 1
                real_info = '_OR_'.join(real_info_list)
            else:
                mismatches += 1

            vcf_info = f"{chrom}_{pos}_{pos2}_{svlen}"
            data.append([svtype, int(found), qual_field, vcf_info, real_info])
    
    for index, real_row in real_variants.iterrows():
        if index not in matched_indexes:
            unmatched_real_variants.append((str(real_row['POS']), str(real_row['SVTYPE']), str(real_row['SVLEN'])))

    for index, real_row in real_variants.iterrows():
        real_pos = real_row['POS']
        real_svlen = abs(real_row['SVLEN'])
        real_svtype = real_row['SVTYPE'].replace(' ', '')
        if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
            real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
            data.append(['ND', 1, -1, '', real_info])
    print()
    print(f'VCF path: {vcf_path}')
    print(f"Total rearrangements in vcf: {vcf_rearrangements}")
    print(f"Matches with real rearrangements: {matches}")
    print(f"Mismatches with real rearrangements: {mismatches}")

    return pd.DataFrame(data, columns=['Rearrangement', 'Present', 'QUAL', 'VCF_Info', 'Real_Info']), unmatched_real_variants

def sniffles_process_vcf_for_dataframe(real_variants, vcf_path):
    data = []
    unmatched_real_variants = []
    matched_indexes = []
    real_sv_matched = []
    vcf_rearrangements = 0 
    matches = 0
    mismatches = 0
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            info_field = columns[7]  # INFO field
            
            qual_field = float(columns[5])  # Use the QUAL value directly
            info_dict = dict(item.split('=') for item in info_field.split(';') if '=' in item)
            
            chrom = columns[0]
            pos = int(columns[1])
            chrom2 = info_dict.get('CHROM2', chrom)
            pos2 = int(info_dict.get('END', pos))
            svlen = abs(int((info_dict.get('SVLEN', 0))))
            svtype = info_dict.get('SVTYPE')
            if svtype == 'DUP':
                svtype = 'TD' if FIXES else svtype
        

            if svlen < 200:
                continue

            svlen_diff = 200

            if svlen > 0 and svlen < 10000:
                svlen_diff = svlen * 0.1
            elif svlen >= 10000:
                svlen_diff = svlen * 0.05

            found = False
            real_info = ""
            real_info_list = []

            vcf_rearrangements += 1

            for index, real_row in real_variants.iterrows():
                real_pos = real_row['POS']
                real_svlen = abs(real_row['SVLEN'])
                real_svtype = real_row['SVTYPE'].replace(" ", "")
                if real_svtype == 'INVTD':
                    # print(real_row)
                    continue

                if 'BND' in svtype:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'INV':
                    if ((svtype in real_svtype or real_svtype in svtype) and abs(real_pos - pos) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or (real_svtype == 'INS_DEL' and abs(real_pos - pos) <= 200):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TD':
                    if (real_svtype == 'INS' and abs(real_pos - pos2) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or ((svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TRL':

                    if (((svtype in real_svtype or real_svtype in svtype) or real_svtype == 'INS') and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200)) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                else:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
            if found:
                for real_info in real_info_list:
                    # print(real_info)
                    real_info = real_info.split('_')
                    real_svtype = real_info[0]
                    try:
                        real_pos = int(real_info[1])
                        real_svlen = int(real_info[2])
                    except ValueError:
                        real_svtype = real_info[0]+'_'+real_info[1]
                        real_pos = int(real_info[2])
                        real_svlen = int(real_info[3])
                    if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
                        real_sv_matched.append((real_pos, real_svlen, real_svtype))
                        matches += 1
                real_info = '_OR_'.join(real_info_list)
            else:
                mismatches += 1

            vcf_info = f"{chrom}_{pos}_{pos2}_{svlen}"
            data.append([svtype, int(found), qual_field, vcf_info, real_info])
    
    for index, real_row in real_variants.iterrows():
        if index not in matched_indexes:
            unmatched_real_variants.append((str(real_row['POS']), str(real_row['SVTYPE']), str(real_row['SVLEN'])))

    for index, real_row in real_variants.iterrows():
        real_pos = real_row['POS']
        real_svlen = abs(real_row['SVLEN'])
        real_svtype = real_row['SVTYPE'].replace(' ', '')
        if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
            real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
            data.append(['ND', 1, -1, '', real_info])
    print()
    print(f'VCF path: {vcf_path}')
    print(f"Total rearrangements in vcf: {vcf_rearrangements}")
    print(f"Matches with real rearrangements: {matches}")
    print(f"Mismatches with real rearrangements: {mismatches}")

    return pd.DataFrame(data, columns=['Rearrangement', 'Present', 'QUAL', 'VCF_Info', 'Real_Info']), unmatched_real_variants

def svim_process_vcf_for_dataframe(real_variants, vcf_path):
    data = []
    unmatched_real_variants = []
    matched_indexes = []
    real_sv_matched = []
    vcf_rearrangements = 0 
    matches = 0
    mismatches = 0
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            info_field = columns[7]  # INFO field
            alt_field = columns[4]
            
            qual_field = float(columns[5])  # Use the QUAL value directly
            info_dict = dict(item.split('=') for item in info_field.split(';') if '=' in item)
            
            chrom1 = columns[0]
            pos = int(columns[1])
            chrom2 = info_dict.get('CHROM2', chrom1)
            pos2 = int(info_dict.get('END', pos))
            svlen = abs(int((info_dict.get('SVLEN', pos2 - pos + 1))))
            svtype = info_dict.get('SVTYPE')
            
            if svtype == 'BND':
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
                    svlen = pos2 - pos

            if 'DUP' in svtype:
                svtype = 'TD' if FIXES else svtype
            
            if svlen < 200:
                continue

            svlen_diff = 200

            if svlen > 0 and svlen < 10000:
                svlen_diff = svlen * 0.1
            elif svlen >= 10000:
                svlen_diff = svlen * 0.05

            found = False
            real_info = ""
            real_info_list = []

            vcf_rearrangements += 1

            for index, real_row in real_variants.iterrows():
                real_pos = real_row['POS']
                real_svlen = abs(real_row['SVLEN'])
                real_svtype = real_row['SVTYPE'].replace(" ", "")
                if real_svtype == 'INVTD':
                    # print(real_row)
                    continue

                if 'BND' in svtype:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'INV':
                    if ((svtype in real_svtype or real_svtype in svtype) and abs(real_pos - pos) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or (real_svtype == 'INS_DEL' and abs(real_pos - pos) <= 200):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TD':
                    if (real_svtype == 'INS' and abs(real_pos - pos2) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or ((svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TRL':

                    if (((svtype in real_svtype or real_svtype in svtype) or real_svtype == 'INS') and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200)) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                else:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
            if found:
                for real_info in real_info_list:
                    # print(real_info)
                    real_info = real_info.split('_')
                    real_svtype = real_info[0]
                    try:
                        real_pos = int(real_info[1])
                        real_svlen = int(real_info[2])
                    except ValueError:
                        real_svtype = real_info[0]+'_'+real_info[1]
                        real_pos = int(real_info[2])
                        real_svlen = int(real_info[3])
                    if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
                        real_sv_matched.append((real_pos, real_svlen, real_svtype))
                        matches += 1
                real_info = '_OR_'.join(real_info_list)
            else:
                mismatches += 1

            vcf_info = f"{chrom1}_{pos}_{pos2}_{svlen}"
            data.append([svtype, int(found), qual_field, vcf_info, real_info])
    
    for index, real_row in real_variants.iterrows():
        if index not in matched_indexes:
            unmatched_real_variants.append((str(real_row['POS']), str(real_row['SVTYPE']), str(real_row['SVLEN'])))

    for index, real_row in real_variants.iterrows():
        real_pos = real_row['POS']
        real_svlen = abs(real_row['SVLEN'])
        real_svtype = real_row['SVTYPE'].replace(' ', '')
        if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
            real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
            data.append(['ND', 1, -1, '', real_info])
    print()
    print(f'VCF path: {vcf_path}')
    print(f"Total rearrangements in vcf: {vcf_rearrangements}")
    print(f"Matches with real rearrangements: {matches}")
    print(f"Mismatches with real rearrangements: {mismatches}")

    return pd.DataFrame(data, columns=['Rearrangement', 'Present', 'QUAL', 'VCF_Info', 'Real_Info']), unmatched_real_variants

def severus_process_vcf_for_dataframe(real_variants, vcf_path):
    data = []
    unmatched_real_variants = []
    matched_indexes = []
    real_sv_matched = []
    vcf_rearrangements = 0 
    matches = 0
    mismatches = 0

    # print(vcf_path)
    # input()
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            # print(line)
            # input()

            if line.startswith('#'):
                continue
            
            columns = line.strip().split('\t')
            info_field = columns[7]
            alt_field = columns[4]

            # print(columns)

            qual_field = float(columns[5])  # используем QUAL напрямую
            # парсим INFO в словарь (только пары ключ=значение)
            info_dict = dict(item.split('=') for item in info_field.split(';') if '=' in item)

            chrom1 = columns[0]
            pos = int(columns[1])
            chrom2 = info_dict.get('CHROM2', chrom1)
            pos2 = int(info_dict.get('END', pos))

            # print(chrom1, pos, chrom2, pos2)
            # input()

            # SVLEN: если есть в INFO — безопасно приводим к int через remove_non_numeric,
            # иначе считаем по координатам (END - POS)
            svlen = None
            if 'SVLEN' in info_dict and info_dict.get('SVLEN') is not None:
                try:
                    svlen = abs(int(remove_non_numeric(info_dict.get('SVLEN'))))
                except Exception:
                    svlen = abs(pos2 - pos)
            else:
                svlen = abs(pos2 - pos)

            svtype = info_dict.get('SVTYPE')

            # Обработка BND (разбор ALT вида ]chrX:pos] или [chrX:pos[ и т.п.)
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
                    svlen = abs(pos2 - pos)
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
            
            

            if svlen < 200:
                continue

            svlen_diff = 200

            if svlen > 0 and svlen < 10000:
                svlen_diff = svlen * 0.1
            elif svlen >= 10000:
                svlen_diff = svlen * 0.05

            found = False
            # print(pos, pos2, svlen)
            real_info = ""
            real_info_list = []

            vcf_rearrangements += 1

            for index, real_row in real_variants.iterrows():
                real_pos = real_row['POS']
                real_svlen = abs(real_row['SVLEN'])
                real_svtype = real_row['SVTYPE'].replace(" ", "")
                if real_svtype == 'INVTD':
                    # print(real_row)
                    continue

                if 'BND' in svtype:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'INV':
                    if ((svtype in real_svtype or real_svtype in svtype) and abs(real_pos - pos) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or (real_svtype == 'INS_DEL' and abs(real_pos - pos) <= 200):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TD':
                    if (real_svtype == 'INS' and abs(real_pos - pos2) <= 200 and abs(real_svlen - svlen) <= svlen_diff) or ((svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff):
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                elif svtype == 'TRL':

                    if (((svtype in real_svtype or real_svtype in svtype) or real_svtype == 'INS') and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200)) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
                else:
                    if (svtype in real_svtype or real_svtype in svtype) and (abs(real_pos - pos) <= 200 or abs(real_pos - pos2) <= 200) and abs(real_svlen - svlen) <= svlen_diff:
                        found = True
                        real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
                        real_info_list.append(real_info)
                        matched_indexes.append(index)
                        # break
            if found:
                # print('found')
                # input()
                for real_info in real_info_list:
                    # print(real_info)
                    real_info = real_info.split('_')
                    real_svtype = real_info[0]
                    try:
                        real_pos = int(real_info[1])
                        real_svlen = int(real_info[2])
                    except ValueError:
                        real_svtype = real_info[0]+'_'+real_info[1]
                        real_pos = int(real_info[2])
                        real_svlen = int(real_info[3])
                    if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
                        real_sv_matched.append((real_pos, real_svlen, real_svtype))
                        matches += 1
                real_info = '_OR_'.join(real_info_list)
            else:
                mismatches += 1

            vcf_info = f"{chrom1}_{pos}_{pos2}_{svlen}"
            data.append([svtype, int(found), qual_field, vcf_info, real_info])
    
    for index, real_row in real_variants.iterrows():
        if index not in matched_indexes:
            unmatched_real_variants.append((str(real_row['POS']), str(real_row['SVTYPE']), str(real_row['SVLEN'])))

    for index, real_row in real_variants.iterrows():
        real_pos = real_row['POS']
        real_svlen = abs(real_row['SVLEN'])
        real_svtype = real_row['SVTYPE'].replace(' ', '')
        if (real_pos, real_svlen, real_svtype) not in real_sv_matched:
            real_info = f"{real_svtype}_{real_pos}_{real_svlen}"
            data.append(['ND', 1, -1, '', real_info])
    print()
    print(f'VCF path: {vcf_path}')
    print(f"Total rearrangements in vcf: {vcf_rearrangements}")
    print(f"Matches with real rearrangements: {matches}")
    print(f"Mismatches with real rearrangements: {mismatches}")

    return pd.DataFrame(data, columns=['Rearrangement', 'Present', 'QUAL', 'VCF_Info', 'Real_Info']), unmatched_real_variants

def aggregate_data_for_samples(samples, process_func, tool_name):
    """Общая функция для агрегации данных по всем инструментам"""
    all_data = []
    real_lrs = 0
    
    for sample in samples:
        unmatched_real_variants_list = [['POS', 'SVTYPE', 'SVLEN']]
        # Читаем CSV файл с реальными вариантами
        real_variants_global = pd.read_csv(sample['csv_path'])
        
        # Обрабатываем SVLEN для разных типов вариантов
        def process_svlen(row):
            if row['SVTYPE'] == 'INS_DEL':
                # print(row['SVLEN'])
                # Для INS_DEL извлекаем оба значения и берем минимальное по модулю
                ins_len, del_len = map(float, row['SVLEN'].split('_'))
                return int(min(abs(ins_len), abs(del_len)))
            else:
                # Для остальных типов просто берем абсолютное значение
                return int(abs(float(row['SVLEN'])))
        
        def process_pos(row):
            return int(abs(float(row['POS'])))
        
        # print(sample['csv_path'])
        # print(real_variants_global)

        # Создаем временный столбец с обработанными длинами
        real_variants_global['SVLEN'] = real_variants_global.apply(process_svlen, axis=1)
        real_variants_global['POS'] = real_variants_global.apply(process_pos, axis=1)
        
        # Фильтруем варианты по длине (>= 200)
        filtered_variants = real_variants_global[real_variants_global['SVLEN'] >= 200]
        real_lrs += len(filtered_variants)
        
        print(f'\nReal LRs in {sample["csv_path"]}: {len(filtered_variants)}')
        name_of_sample = sample.get('name', os.path.basename(sample['vcf_folder']))
        
        # Обрабатываем все VCF файлы в папке
        for vcf_file in os.listdir(sample['vcf_folder']):
            if vcf_file.endswith('.vcf'):
                if FIXES and 'ERR' in vcf_file and name_of_sample not in vcf_file:
                    print(f'skipped, another sample: {vcf_file}')
                    continue
                vcf_path = os.path.join(sample['vcf_folder'], vcf_file)
                print(f"Обработка VCF: {vcf_file}")
                
                # print('filtered_variants\n', filtered_variants)
                # input()
                # Обработка VCF с помощью специфичной для инструмента функции
                result = process_func(filtered_variants, vcf_path)
                
                # Обработка возвращаемого значения
                if isinstance(result, tuple):
                    df, unmatched_real_variants = result
                    unmatched_real_variants_list.extend(unmatched_real_variants)
                else:
                    df = result
                
                # print(df)
                # input()
                if df is not None:
                    df['Sample'] = name_of_sample
                all_data.append(df)
        
        # Сохраняем несопоставленные варианты (если есть)
        if unmatched_real_variants_list:
            output_path = f'{output_dir}/unmatched_real_variants_{tool_name}_{name_of_sample}.csv'
            # pd.DataFrame(unmatched_real_variants_list).to_csv(output_path, index=False)
            with open(output_path, 'w') as unmatch_file:
                for um_line in unmatched_real_variants_list:
                    unmatch_file.write(','.join(um_line)+'\n')
            print(f"Сохранены несопоставленные варианты: {output_path}")
    
    print(f'\nОбщее количество реальных LRs: {real_lrs}')
    return pd.concat(all_data, ignore_index=True)

# Для elarodon
def elarodon_aggregate_data_for_samples(samples):
    return aggregate_data_for_samples(samples, elarodon_process_vcf_for_dataframe, 'eLaRodON')

# Для NanoSV
def nanosv_aggregate_data_for_samples(samples):
    return aggregate_data_for_samples(samples, nanosv_process_vcf_for_dataframe, 'NanoSV')

# Для Nanovar
def nanovar_aggregate_data_for_samples(samples):
    return aggregate_data_for_samples(samples, nanovar_process_vcf_for_dataframe, 'Nanovar')

# Для Sniffles
def sniffles_aggregate_data_for_samples(samples):
    return aggregate_data_for_samples(samples, sniffles_process_vcf_for_dataframe, 'Sniffles')

# Для SVIM
def svim_aggregate_data_for_samples(samples):
    return aggregate_data_for_samples(samples, svim_process_vcf_for_dataframe, 'SVIM')

def sniffles_mosaic_aggregate_data_for_samples(samples):
    return aggregate_data_for_samples(samples, sniffles_process_vcf_for_dataframe, 'Sniffles2_mosaic')

def severus_aggregate_data_for_samples(samples):
    return aggregate_data_for_samples(samples, severus_process_vcf_for_dataframe, 'Severus')

def plot_combined_roc_curve(dataframes, titles, output_path):
    plt.figure(figsize=(12, 8))

    for df, title in zip(dataframes, titles):
        y_true = df['Present'].astype(int)  # Ensure y_true is binary (0 or 1)
        y_scores = df['QUAL']

        fpr, tpr, thrs = roc_curve(y_true, y_scores)
        roc_auc = roc_auc_score(y_true, y_scores)

        plt.plot(fpr, tpr, lw=2, label=f'{title}')

    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=23, fontweight='bold')
    plt.ylabel('True Positive Rate', fontsize=23, fontweight='bold')
    plt.xticks(fontsize=23, fontweight='bold')
    plt.yticks(fontsize=23, fontweight='bold')
    
    # Добавляем жирный шрифт в легенду
    plt.legend(
        loc="upper left", 
        fontsize=23,
        prop={'weight': 'bold', 'size': 23}  # Жирный шрифт и размер
    )
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=600)
    plt.close()


# Конфигурация путей
csv_root = '/media/regina/MR_GREY/mtub/sens_compare/csv_checked_manually_INSDEL'
vcf_root = path_arg('--vcf-root', '/media/regina/MR_GREY/mtub/sens_compare/vcf_all')
output_dir = sys.argv[1]
os.makedirs(output_dir, exist_ok=True)

# Создаем списки для каждого инструмента
samples_el = create_samples_list('elarodon')
samples_nsv = create_samples_list('nanosv')
samples_nv = create_samples_list('nanovar')
samples_sf = create_samples_list('sniffles')
samples_sv = create_samples_list('svim')
samples_sev = create_samples_list('severus')
samples_sfm = create_samples_list('', root=MOSAIC_ROOT)


# Обработка каждого инструмента
tool_datasets = [
    {'samples': samples_el, 'name': 'eLaRodON', 'func': elarodon_aggregate_data_for_samples, 'output': 'aggregated_df_eLaRodON'},
    {'samples': samples_nsv, 'name': 'NanoSV', 'func': nanosv_aggregate_data_for_samples, 'output': 'aggregated_df_NanoSV'},
    {'samples': samples_nv, 'name': 'Nanovar', 'func': nanovar_aggregate_data_for_samples, 'output': 'aggregated_df_Nanovar'},
    {'samples': samples_sv, 'name': 'SVIM', 'func': svim_aggregate_data_for_samples, 'output': 'aggregated_df_SVIM'},
    {'samples': samples_sf, 'name': 'Sniffles2', 'func': sniffles_aggregate_data_for_samples, 'output': 'aggregated_df_Sniffles2'},
    {'samples': samples_sfm, 'name': 'Sniffles2_mosaic', 'func': sniffles_mosaic_aggregate_data_for_samples, 'output': 'aggregated_df_Sniffles2_mosaic'},
    {'samples': samples_sev, 'name': 'Severus', 'func': severus_aggregate_data_for_samples, 'output': 'aggregated_df_Severus'}
]

if LOWSUPPORT:
    def lowsupport_samples(tool_dir, sub=''):
        # same as create_samples_list, the VCF may lie in a subfolder (Severus: all_SVs)
        samples = []
        for s in create_samples_list(tool_dir, root=LOWSUPPORT_ROOT):
            samples.append({'csv_path': s['csv_path'], 'vcf_folder': os.path.join(s['vcf_folder'], sub),
                            'name': os.path.basename(s['vcf_folder'])})
        return samples
    for name, tool_dir, sub, func in [
            ('Sniffles2_mosaic_germline', 'sniffles_mosaic_germline', '', sniffles_process_vcf_for_dataframe),
            ('NanoVar_mincov1', 'nanovar_mincov1', '', nanovar_process_vcf_for_dataframe),
            ('NanoSV_cluster1', 'nanosv_cluster1', '', nanosv_process_vcf_for_dataframe),
            ('Severus_minsupport1', 'severus_minsupport1', 'all_SVs', severus_process_vcf_for_dataframe)]:
        tool_datasets.append({'samples': lowsupport_samples(tool_dir, sub), 'name': name,
                              'func': lambda smp, f=func, n=name: aggregate_data_for_samples(smp, f, n),
                              'output': 'aggregated_df_' + name})

# R1.2 test versions of eLaRodON: --elarodon-variants baseline,fix_td  (folders elarodon_test_fixes/out_<variant>/)
ELARODON_TEST_ROOT = '/home/regina/work/elarodon_test_fixes'   # MR_GREY is full; scripts in rel7/New_files/elarodon_test_fixes
variants_arg = [a for a in sys.argv if a.startswith('--elarodon-variants=')]
if variants_arg:
    for variant in variants_arg[0].split('=', 1)[1].split(','):
        name = 'eLaRodON_' + variant
        tool_datasets.append({'samples': create_samples_list('out_' + variant, root=ELARODON_TEST_ROOT), 'name': name,
                              'func': lambda smp, n=name: aggregate_data_for_samples(smp, elarodon_process_vcf_for_dataframe, n),
                              'output': 'aggregated_df_' + name})

df_list = []
df_titles = []

for dataset in tool_datasets:
    if not dataset['samples']:
        print(f"Пропуск {dataset['name']}: нет образцов для обработки")
        continue
    
    print('-'*len(dataset['name']))
    print(dataset['name'])
    print('-'*len(dataset['name']))
    
    # Агрегируем данные
    print(dataset['samples'])
    df = dataset['func'](dataset['samples'])
    
    # Сохраняем результаты
    output_path = os.path.join(output_dir, dataset['output'] + '.csv')
    df.to_csv(output_path, index=False)
    print(f"Результаты сохранены в {output_path}\n")
    
    # Добавляем в список для ROC-кривой
    df_list.append(df)
    df_titles.append(dataset['name'])

# Построение комбинированного графика ROC
# plot_combined_roc_curve(
#     df_list,
#     df_titles,
#     output_path='/media/regina/KINGSTON/mtub/sens_compare/plots_ROC/combined_ROC_curves_mtub.png'
# )

metrics = []
for df, title in zip(df_list, df_titles):
    # ensure columns exist and types are correct
    y_true = df['Present'].astype(int).values
    # если QUAL может быть пустым — подставляем 0
    y_scores = df['QUAL'].fillna(0).astype(float).values

    # безопасно считаем AUC (если в y_true присутствует только один класс, roc_auc_score упадёт)
    try:
        auc_val = roc_auc_score(y_true, y_scores)
    except Exception:
        auc_val = float('nan')

    # пытаемся найти порог, максимизирующий F1, через precision_recall_curve
    # try:
    precision, recall, pr_thresholds = precision_recall_curve(y_true, y_scores)
    # precision_recall_curve возвращает len(precision) = len(recall) = len(thresholds)+1
    # Соответственно для порогов используем precision[1:], recall[1:]
    if len(pr_thresholds) > 0:
        f1s_at_thr = 2 * precision[1:] * recall[1:] / (precision[1:] + recall[1:] + 1e-12)
        best_idx = int(np.nanargmax(f1s_at_thr))
        best_f1 = float(f1s_at_thr[best_idx])
        best_threshold = float(pr_thresholds[best_idx])
    else:
        # нет порогов (возможно все оценки одинаковы) — делаем fallback:
        # пробуем порог = медиана или 0.5
        fallback_thr = np.median(y_scores) if np.nanmedian(y_scores) != 0 else 0.5
        y_pred = (y_scores >= fallback_thr).astype(int)
        best_f1 = float(f1_score(y_true, y_pred))
        best_threshold = float(fallback_thr)
    # except Exception:
    #     # в случае ошибки — ставим NaN
    #     best_f1 = float('nan')
    #     best_threshold = float('nan')

    metrics.append({
        'program': title,
        'auc': auc_val,
        'f1': best_f1,
        'best_threshold': best_threshold
    })

metrics_df = pd.DataFrame(metrics)

# Сохраняем таблицу (путь можно поменять)
out_metrics_csv = os.path.join(output_dir, 'metrics_table_mtub.csv')
metrics_df.to_csv(out_metrics_csv, index=False)

print("Metrics table saved to:", out_metrics_csv)
print(metrics_df)

# Теперь обновим подписи для графика, чтобы показать AUC и F1 прямо в легенде
df_titles_with_metrics = []
for row in metrics_df.itertuples(index=False):
    try:
        auc_display = f"{row.auc:.2f}" if not np.isnan(row.auc) else "nan"
    except Exception:
        auc_display = "nan"
    try:
        f1_display = f"{row.f1:.2f}" if not np.isnan(row.f1) else "nan"
    except Exception:
        f1_display = "nan"
    df_titles_with_metrics.append(f"{row.program} (AUC={auc_display}, F1={f1_display})")

# Строим график с обновлёнными подписями
plot_combined_roc_curve(
    df_list,
    df_titles_with_metrics,
    output_path=os.path.join(output_dir, 'combined_ROC_curves_mtub_200.png')
)


# ---------- precision-recall ----------
def golden_set_ids():
    """(sample, 'TYPE_POS_LEN') for golden-set variants >= 200 bp, as the matching above sees them."""
    ids = set()
    for fname in os.listdir(csv_root):
        if not fname.endswith('.csv'):
            continue
        sample = fname[:-len('.csv')]
        for _, row in pd.read_csv(os.path.join(csv_root, fname)).iterrows():
            svtype = row['SVTYPE'].replace(' ', '')
            if row['SVTYPE'] == 'INS_DEL':
                ins_len, del_len = map(float, row['SVLEN'].split('_'))
                svlen = int(min(abs(ins_len), abs(del_len)))
            else:
                svlen = int(abs(float(row['SVLEN'])))
            if svlen >= 200:
                ids.add((sample, f"{svtype}_{int(abs(float(row['POS'])))}_{svlen}"))
    return ids


truth_ids = golden_set_ids()
print('Golden set variants:', len(truth_ids))
pr_rows, pr_curves = [], []
for df, title in zip(df_list, df_titles):
    program = 'NanoVar' if title == 'Nanovar' else title
    calls = []
    for _, row in df[df['Rearrangement'] != 'ND'].iterrows():
        matched = set()
        if row['Present'] == 1:
            matched = {(row['Sample'], x) for x in row['Real_Info'].split('_OR_')}
            assert matched <= truth_ids, matched - truth_ids
        calls.append((float(row['QUAL']), matched))
    summary, (p, r) = pr_utils.summarize(program, calls, len(truth_ids))
    pr_rows.append(summary)
    pr_curves.append(pr_utils.curve_table(program, calls, len(truth_ids)))
pd.DataFrame(pr_rows).to_csv(os.path.join(output_dir, 'pr_metrics_mtub.tsv'), sep='\t', index=False)
print(pd.DataFrame(pr_rows).to_string())
pd.concat(pr_curves).to_csv(os.path.join(output_dir, 'pr_points_mtub.tsv'), sep='\t', index=False)
