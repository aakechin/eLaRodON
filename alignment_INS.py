import argparse
import subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import shutil
from statistics import median, mean
import pysam
import sys
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import logging

class InsertionProcessor:
    def __init__(self, log_path, ins_file, lr_file, ref_genome, name, sam_file=None, threads=6, not_remove_trash=False):
        self.log_path = log_path
        logging.basicConfig(level=logging.DEBUG, filename=self.log_path, filemode="a",
                            format="%(asctime)s %(funcName)s %(lineno)d %(message)s")
        self.ins_file = ins_file
        self.lr_file = lr_file
        self.ref_genome = ref_genome
        self.name = name
        self.sam_file = sam_file
        self.insertions = {}
        self.tandem_duplications = []
        self.tandem_duplications_ids = []
        self.threads = threads
        self.not_remove_trash = not_remove_trash

    def read_insertions(self):
        print('Reading INS-file...')
        logging.info('Reading INS-file...')
        df = pd.read_csv(self.ins_file, delimiter='\t')
        for index, row in df.iterrows():
            ins_id = index
            if ',' in row['Read_Names']:
                row['Read_Number'] = str(len(row['Read_Names'].split(',')))
            self.insertions[ins_id] = row

    def write_fasta(self, output_file):
        print('Writing FASTA-file...')
        logging.info('Writing FASTA-file...')
        logging.info('output_file: '+output_file)
        with open(output_file, 'w') as fasta_file:
            for ins_id, data in self.insertions.items():
                # print('data', data)
                # input()
                seq_id = str(ins_id)+'*'+data['Chrom']+'*'+str(data['Pos'])+'*'+str(data['Insertion_Length'])
                sequence = SeqRecord(Seq(data['Sequence']), id=seq_id, description='')
                SeqIO.write(sequence, fasta_file, 'fasta')
        # input('STOP FASTA')

    def map_with_minimap2(self, fasta_file):
        print('Mapping...')
        logging.info('Mapping...')
        if self.sam_file and os.path.exists(self.sam_file):
            print("Using existing SAM file...")
            logging.info('"Using existing SAM file...')
            return
        output_sam = self.name + "_mapping_results.sam"
        logging.info('output_sam: '+output_sam)
        subprocess.run([
            'minimap2', '-a', self.ref_genome, fasta_file],
            stdout=open(output_sam, 'w')
        )
        logging.info(' '.join(['minimap2', '-a', self.ref_genome, fasta_file]))
        self.sam_file = output_sam
        # input('STOP mapping')

    def analyze_mapping_results(self):
        print('Searching tandem duplications...')
        logging.info('Searching tandem duplications...')
        with open(self.sam_file, 'r') as sam_file:
            for line in sam_file:
                if line.startswith('@'):
                    continue
                parts = line.split('\t')
                # print('parts', parts)
                # input()
                read_id = parts[0]
                ins_id, chrom, coordinate, length = read_id.split('*')
                ins_id = int(ins_id)
                coordinate = int(coordinate)
                length = int(length)

                map_chrom = parts[2]
                map_pos = int(parts[3])
                if map_chrom == chrom and coordinate - length - 50 < map_pos < coordinate + length + 50:
                    self.tandem_duplications.append((ins_id, map_pos, length))
                    self.tandem_duplications_ids.append(ins_id)
        # print('self.tandem_duplications', self.tandem_duplications)
        # print()
        # input('STOP SAM')

    def merge_tandem_duplications(self):
        print('Merging tandem duplications...')
        logging.info('Merging tandem duplications...')
        self.merged_insertions = defaultdict(list)
        for ins_id_1, map_pos_1, length_1 in self.tandem_duplications:
            for ins_id_2, map_pos_2, length_2 in self.tandem_duplications:
                if abs(map_pos_1 - map_pos_2) < 50 and abs(length_1 - length_2):
                    self.merged_insertions[ins_id_1].append((ins_id_2, map_pos_2, length_2))
                    self.merged_insertions[ins_id_2].append((ins_id_1, map_pos_1, length_1))
        # print('self.merged_insertions', self.merged_insertions)
        # input('STOP merging')

    def copy_file(self, source_path, destination_path):
        try:
            shutil.copy(source_path, destination_path)
            # print(f"File copied from {source_path} to {destination_path}")
        except Exception as e:
            print("Error occurred while copying file: " + e)
            logging.error("Error occurred while copying file: " + e)
    
    def get_sequence_from_genome(self, chrom_name, temp_start, temp_end, diff, flag):
        # attempts = 0
        seq = None
        while (seq is None):
            if flag == 'L':
                try:
                    seq = str(ff.fetch(region=str(chrom_name)+':'+str(temp_start-diff)+'-'+str(temp_end)))
                except:
                    seq = None
                    if diff > 5:
                        diff -= 5
                        if diff < 5:
                            diff = 5
                            continue
            elif flag == 'R':
                try:
                    seq = str(ff.fetch(region=str(chrom_name)+':'+str(temp_start)+'-'+str(temp_end+diff)))
                except:
                    seq = None
                    if diff > 5:
                        diff -= 5
                        if diff < 5:
                            diff = 5
                            continue
            # if attempts == 4:
            #     diff = 2
            # attempts += 1
            if seq is None:
                ff = pysam.FastaFile(self.ref_genome, filepath_index=self.ref_genome+'.fai')
                if diff <= 5:
                    diff -= 1
                    if diff >= 1:
                        continue
                    else:
                        # print()
                        # print('Error! Unable to extract the sequence from the reference genome!')
                        # print(chrom_name, temp_start, temp_end, diff, flag)
                        return ''
                        # sys.exit()
        # if flag == 'R':
        #     print(chrom_name, temp_start, temp_end, diff, flag)
        #     print(seq)
        return (seq.upper())
    
    def process_extraction(self, data_ins):

        left_s = self.get_sequence_from_genome(data_ins['Chrom'], int(data_ins['map_pos']), int(data_ins['map_pos']), 30, 'L')
        right_s = self.get_sequence_from_genome(data_ins['Chrom'], int(data_ins['Pos']), int(data_ins['Pos']), 30, 'R')

        return (data_ins['ID'], left_s, right_s)

    def extract_sequences_before_after(self):
        print('Extracting sequences before and after insertions...')
        logging.info('Extracting sequences before and after insertions...')
        rows = []
        for ins_id, map_pos, length in self.tandem_duplications:  
            self.insertions[ins_id]['ID'] = ins_id
            self.insertions[ins_id]['map_pos'] = ins_id
            rows.append(self.insertions[ins_id])

        # create pool
        num_processes = self.threads
        pool = Pool(processes=num_processes)

        # give single rows to pool
        if len(rows) < 1000:
            chunk_size_val = 100
        elif len(rows) < 10000:
            chunk_size_val = 1000
        else:
            chunk_size_val = 10000

        results = list(tqdm(pool.imap_unordered(self.process_extraction, rows, chunksize=chunk_size_val), total=len(rows), desc="Processing insertions")) #, colour='blue'))
        # results = pool.map(self.process_row_INS, rows)

        pool.close()
        pool.join()

        for result in results:
            id = result[0]
            seq_bef = result[1]
            seq_aft = result[2]
            self.insertions[id]['SBB'] = seq_bef
            self.insertions[id]['SAB'] = seq_aft


    def update_large_rearrangements_file(self):
        print('Updating LRs data...')
        logging.info('Updating LRs data...')
        uniq = []
        with open(self.lr_file[:-4]+'_mapped_ins.csv', 'a') as lr_file:
            for ins_id, map_pos, length in self.tandem_duplications:
                if ins_id not in self.merged_insertions:
                    ins_info = self.insertions[ins_id]
                    # print('ins_info', ins_info)
                    # input()
                    lr_file.write('\t'.join([str(ins_info['Chrom']), # Chrom1
                                            str(map_pos)+','+str(ins_info['Pos']), # Pos1
                                            '+', # Strand1
                                            'D', # Junction_Side1
                                            str(ins_info['Chrom']), # Chrom2
                                            str(map_pos)+','+str(map_pos+length), # Pos2
                                            '+', # Strand2
                                            'D', # Jinction_Side2
                                            str(ins_info['Read_Number']), # Read_Number
                                            str(int(ins_info['Total_coverage']) - int(ins_info['Read_Number'])), # Ref_Coverage
                                            str(ins_info['MQ_median']), # MQ_median
                                            '99', # Read_Len1
                                            '99', # Read_Len2
                                            str(ins_info['Read_Muts']), # Read_Muts
                                            '0', # Inter_Joints
                                            '1', # Scarless_Joints
                                            '0', # NewSeq_Joints
                                            '1', # Cigar_Num
                                            '0', # Inversion_Num
                                            '0', # Start_End_Locations
                                            '0', # Number_of_parts
                                            ','.join(ins_info['SBB']), # Seqs_before_BND
                                            ','.join(ins_info['SAB']), # Seqs_after_BND
                                            ins_info['Read_Names'] # Read_names
                                            ]) + '\n')
                    uniq.append(ins_id)
                else:
                    ins_ids = [ins_id]
                    ins_ids.extend([i for i in self.merged_insertions[ins_id]])

                    ins_info_main = [self.insertions[ins_id]]
                    # print('ins_info_main', ins_info_main)
                    # input()
                    ins_info_main.extend([self.insertions[id[0]] for id in self.merged_insertions[ins_id]])

                    map_poses = [map_pos]
                    map_poses.extend([pos[1] for pos in self.merged_insertions[ins_id]])

                    map_lens = [length]
                    map_lens.extend([l[2] for l in self.merged_insertions[ins_id]])
                    
                    poses = [x['Pos'] for x in ins_info_main]
                    chroms = [x['Chrom'] for x in ins_info_main]
                    rn = [int(x['Read_Number']) for x in ins_info_main]
                    total_covs = median([int(x['Total_coverage']) for x in ins_info_main])
                    mq = median([float(x['MQ_median']) for x in ins_info_main])
                    # input()
                    muts = ','.join([str(x['Read_Muts']) for x in ins_info_main])
                    # print('ins_info_main', ins_info_main)
                    # print([str(x['Read_Muts']) for x in ins_info_main])
                    # print('muts', muts)
                    rnames = ','.join([str(x['Read_Names']) for x in ins_info_main])

                    left_seqs = [x['SBB'] for x in ins_info_main]      
                    right_seqs = [x['SAB'] for x in ins_info_main]            
                    

                    if ins_id not in uniq:
                        lr_file.write('\t'.join([str(ins_info_main[0]['Chrom']), # Chrom1
                                            str(min(map_poses))+','+str(max(poses)), # Pos1
                                            '+', # Strand1
                                            'D', # Junction_Side1
                                            str(ins_info_main[0]['Chrom']), # Chrom2
                                            str(round(mean(map_poses)))+','+str(round(mean(map_poses))+round(mean(map_lens))), # Pos2
                                            '+', # Strand2
                                            'D', # Jinction_Side2
                                            str(sum(rn)), # Read_Number
                                            str(max(total_covs - sum(rn), 0)), # Ref_Coverage
                                            str(mq), # MQ_median
                                            '99', # Read_Len1
                                            '99', # Read_Len2
                                            str(muts), # Read_Muts
                                            '0', # Inter_Joints
                                            '1', # Scarless_Joints
                                            '0', # NewSeq_Joints
                                            str(len(poses)), # Cigar_Num
                                            '0', # Inversion_Num
                                            '0', # Start_End_Locations
                                            '0', # Number_of_parts
                                            ','.join(left_seqs),
                                            ','.join(right_seqs),
                                            rnames # Read_names
                                            ]) + '\n')
                    for id in ins_ids:
                        uniq.append(id)  
        # input('STOP writing LRs')     

    def update_original_csv(self, output_file):
        print('Updating INS data...')
        logging.info('Updating INS data...')
        df = pd.read_csv(self.ins_file, delimiter='\t')
        filtered_df = df[~df.index.isin(self.tandem_duplications_ids)]
        filtered_df.to_csv(output_file, sep='\t', index=False)

    # def main(self):
    #     processor = InsertionProcessor(self.ins_file, self.lr_file, self.ref_genome, self.name, self.sam_file, self.threads, self.not_remove_trash)
    #     processor.read_insertions()
    #     fasta_file = self.name + '_insertions.fasta'
    #     processor.write_fasta(fasta_file)
    #     processor.map_with_minimap2(fasta_file)
    #     processor.analyze_mapping_results()
    #     processor.merge_tandem_duplications()
    #     processor.copy_file(self.lr_file, self.lr_file[:-4]+'_mapped_ins.csv')
    #     processor.extract_sequences_before_after()
    #     processor.update_large_rearrangements_file()
    #     updated_csv = self.ins_file[:-4]+'_mapped_ins.csv'
    #     processor.update_original_csv(updated_csv)

    #     print('Analysis of INS was done!')

# if __name__ == '__main__':
#     main()
