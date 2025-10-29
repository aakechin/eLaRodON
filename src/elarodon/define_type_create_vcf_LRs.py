import pandas as pd
import pysam
import sys
import math
from statistics import median, variance
from datetime import datetime
import argparse
from collections import OrderedDict
from sys import stdout
import subprocess as sp
from copy import copy, deepcopy
import os
from Bio import Seq
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from shutil import rmtree
import toml
import gzip
import logging


class AnalyzeLR():

    def __init__(self, args):
        self.path_to_junc_file = args.juncFile
        self.path_to_ins_file = args.insFile
        self.log_path = args.workDir+'/elarodon.log'
        logging.basicConfig(level=logging.DEBUG, filename=self.log_path, filemode="a",
                            format="%(asctime)s %(funcName)s %(lineno)d %(message)s")
        # if os.path.isfile(self.path_to_junc_file[:-4] + '_mapped_ins.csv'):
        #     self.path_to_junc_file = self.path_to_junc_file[:-4] + '_mapped_ins.csv'
        # if os.path.isfile(self.path_to_ins_file[:-4] + '_mapped_ins.csv'):
        #     self.path_to_ins_file = self.path_to_ins_file[:-4] + '_mapped_ins.csv'
        self.path_to_bam_file = args.bamFile
        self.reference_genome = args.refGen
        self.path_vcf_out = args.outVCF
        self.path_vcfanno = args.VCFanno
        # self.path_lua = args.lua_file
        # self.path_toml = args.toml_file
        temp_path = args.workDir +'/supplementary/'
        self.path_bed = args.bedFile
        self.threads = args.threads
        self.name_sample = self.path_to_junc_file[self.path_to_junc_file.rfind('/')+1:self.path_to_junc_file.find('.')]
        self.temp_dir = temp_path+'temp_'+self.path_vcf_out.split('/')[-1][:-4]
        self.temp_vcf_file_pos1 = self.temp_dir+'/temp1_' + self.path_vcf_out.split('/')[-1]
        self.temp_vcf_file_pos2 = self.temp_dir+'/temp2_' + self.path_vcf_out.split('/')[-1]
        self.temp_vcf_ann_1 = self.temp_dir+'/res_ann_1_' + self.path_vcf_out.split('/')[-1]
        self.temp_vcf_ann_2 = self.temp_dir+'/res_ann_2_' + self.path_vcf_out.split('/')[-1]
        self.path_lua = self.temp_dir+'custom.lua'
        self.path_toml = self.temp_dir+'conf.toml'
        self.not_remove_trash = args.notRemoveTrashAnno

    # show percentage of work
    def show_perc_work(self, done_work, all_work):
        perc_done_work = round((done_work/all_work)*100, 2)
        stdout.write('\r'+str(perc_done_work)+'%')
        stdout.flush()

    # to create complementary sequence
    def complement_seq(self, nuc):
        return (str(Seq.Seq(nuc).complement())).upper()
    
    # read the bam file and saving contigs 
    def read_bam_file(self):
        samfile = pysam.AlignmentFile(self.path_to_bam_file, "rb")
        iter = samfile.fetch()
        print()
        print('Reading BAM-file...')
        logging.info('Reading BAM-file...')
        print()
        for x in iter:
            header = x.header
            names = header.references
            lengths = header.lengths
            self.contigs_dict = dict(zip(names, lengths))
            break
        return self.contigs_dict
    
    # sort chromosome names (contigs) for correct recording into vcf-file
    def sort_chrom_names(self):
        general_names_int = []
        general_names_str = []
        additional_names = []
        for name, len in self.contigs_dict.items():
            if '.' in name or '_' in name:
                additional_names.append(name)
            else:
                name = "".join(c for c in name if c.isdigit() or c in ['X', 'Y', 'M'])
                try:
                    general_names_int.append(int(name))
                except:
                    general_names_str.append(name)
        # general_names_int contains the main chromosomes (contigs) names include numbers (chr1, chr2, chr10)
        # general_names_str contains the main chromosomes (contigs) names include letters (chrX, chrY, chrM)
        # additional_names contains specific chromosome (contigs) names that include letters and other characters (e.g. J02459.1)
        general_names_int = list(sorted(general_names_int))
        additional_names = list(sorted(additional_names))
        self.general_dict_LR = OrderedDict()
        for name in general_names_int:
            self.general_dict_LR['chr'+str(name)] = []
        for name in general_names_str:
            self.general_dict_LR['chr'+str(name)] = []
        for name in additional_names:
            self.general_dict_LR[name] = []
        # print(self.general_dict_LR)
        for name, len in self.contigs_dict.items():
            self.general_dict_LR[name].append(len)
            self.general_dict_LR[name].append([])
    
    # write a new title for the vcf file
    def build_new_header(self):
        self.new_header = []
        header_list_1 = ['##fileformat=VCFv4.2',
                        '##fileDate='+str(datetime.strftime(datetime.now(), "%d-%m-%Y")), '##source=NAME_OF_SOURCE', 
                        '##source_reads='+self.path_to_bam_file, '##reference='+self.reference_genome, 
                        '##mapping=-', '##phasing=none', '##depth_of_coverage=-']
        header_list_2 = ['##ALT=<ID=<INS>,Description="Insertion of novel sequence relative to the reference">',
                        '##ALT=<ID=<INV>,Description="Inversion relative to the reference">', 
                        '##ALT=<ID=<BND_INV>,Description="Inversion relative to the reference (unconfirmed by multiple reads)">', 
                        '##ALT=<ID=<DEL>,Description="Deletion relative to the reference">',
                        '##ALT=<ID=<BND_DEL>,Description="Deletion relative to the reference (unconfirmed by multiple reads)">',
                        '##ALT=<ID=<TD>,Description="Tandem duplication">',
                        '##ALT=<ID=<BND_TD>,Description="Possible tandem duplication (one of the fragments is remoted)">',
                        '##ALT=<ID=<INVTD>,Description="Inverted tandem duplication">',
                        '##ALT=<ID=<BND_INVTD>,Description="Inverted tandem duplication (unconfirmed by multiple reads)">',
                        '##ALT=<ID=<TRL>,Description="Translocation">',
                        '##ALT=<ID=<BND_TRL>,Description="Two pieces of DNA adjacent to the breakpoint positions are directly connected to each other in the sample">',
                        
                        '##FILTER=<ID=FAIL,Description="None of genomic rearrangement types">',
                        '##FILTER=<ID=PASS,Description="Type of genomic rearrangement was determined">',

                        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of genomic rearrangement (structural variant)">', 
                        '##INFO=<ID=CHROM2,Number=1,Type=String,Description="Second chromosome involved in the genomic rearrangement">',
                        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the genomic rearrangement">', 
                        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Exact length of the insertion/deletion or approximate lenght of the tandem duplication/inversion">',
                        '##INFO=<ID=TDRN,Number=1,Type=Integer,Description="Number of tandem repeats">',
                        
                        '##INFO=<ID=J1,Number=1,Type=String,Description="On which side (relative to pos1) is the piece of read (L - left, R - right, ND - not defined)">',
                        '##INFO=<ID=J2,Number=1,Type=String,Description="On which side (relative to pos2) is the piece of read (L - left, R - right, ND - not defined)">',
                        '##INFO=<ID=S1,Number=1,Type=String,Description="Strand of the read with pos1">',
                        '##INFO=<ID=S2,Number=1,Type=String,Description="Strand of the read with pos2">',
                        '##INFO=<ID=D1,Number=1,Type=Integer,Description="First distance between breakpoints">', 
                        '##INFO=<ID=D2,Number=1,Type=Integer,Description="Second distance between breakpoints">', 

                        '##INFO=<ID=RL1,Number=1,Type=Integer,Description="First read length">',
                        '##INFO=<ID=RL2,Number=1,Type=Integer,Description="Second read length">',
                        '##INFO=<ID=SR,Number=1,Type=Integer,Description="Number of supporting reads">',
                        '##INFO=<ID=MQ,Number=1,Type=Float,Description="Measure of the confidence that a read actually comes from the position it is aligned to by the mapping algorithm">',
                        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth (total coverage)">',
                        '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">',

                        '##INFO=<ID=MH,Number=1,Type=Integer,Description="Length of the microhomology sequence at the ends of the rearrangement">', 
                        '##INFO=<ID=MHS,Number=1,Type=String,Description="Microhomology sequence at the ends of the rearrangement">', 
                        '##INFO=<ID=HOM,Number=1,Type=Integer,Description="Length of the homeology sequence at the ends of the rearrangement">', 
                        '##INFO=<ID=HOMS,Number=1,Type=String,Description="Homeology sequence at the ends of the rearrangement">',

                        '##INFO=<ID=MUTM,Number=1,Type=Float,Description="Median number of percentage errors in mapping read fragments">',
                        '##INFO=<ID=MUTV,Number=1,Type=Float,Description="Variance of percentage errors in mapping read fragments">',
                        
                        '##INFO=<ID=SBLR,Number=1,Type=String,Description="Possible second boundry of LGR">',
                        '##INFO=<ID=SBTRL,Number=1,Type=String,Description="Possible translocation: its length">',
                        '##INFO=<ID=SBINV,Number=1,Type=String,Description="Possible inversion: its length">',
                        '##INFO=<ID=SBTD,Number=1,Type=String,Description="Possible tandem duplication: its length">',

                        '##INFO=<ID=LERD,Number=1,Type=Integer,Description="Distance from the left position to the nearest external repeat (or mobile element)">',
                        '##INFO=<ID=LERN,Number=1,Type=Integer,Description="Name of the external repeat (or mobile element) on the left side">',
                        '##INFO=<ID=LIRD,Number=1,Type=Integer,Description="Distance from the left position to the nearest internal repeat (or mobile element)">',
                        '##INFO=<ID=LIRN,Number=1,Type=Integer,Description="Name of the internal repeat (or mobile element) on the left side">',
                        '##INFO=<ID=LCRN,Number=1,Type=Integer,Description="Name of the repeat (or mobile element) that contains the left coordinate">',
                        '##INFO=<ID=RERD,Number=1,Type=Integer,Description="Distance from the right position to the nearest external repeat (or mobile element)">',
                        '##INFO=<ID=RERN,Number=1,Type=Integer,Description="Name of the external repeat (or mobile element) on the right side">',
                        '##INFO=<ID=RIRD,Number=1,Type=Integer,Description="Distance from the right position to the nearest internal repeat (or mobile element)">',
                        '##INFO=<ID=RIRN,Number=1,Type=Integer,Description="Name of the internal repeat (or mobile element) on the right side">',
                        '##INFO=<ID=RCRN,Number=1,Type=Integer,Description="Name of the repeat (or mobile element) that contains the right coordinate">',        
                        
                        '##INFO=<ID=ISN,Number=1,Type=Integer,Description="Number of reads with overlapping fragments of reads (InterSeq number)">',
                        '##INFO=<ID=ISNMFS,Number=1,Type=Integer,Description="Length of the most frequent sequence of reads with overlapping fragments of reads (InterSeq number)">',
                        '##INFO=<ID=NSN,Number=1,Type=Integer,Description="Number of reads with unmapped fragment between them (NewSeq number)">',
                        '##INFO=<ID=NSNMFS,Number=1,Type=Integer,Description="Length of the most frequent sequence of reads with unmapped fragment between them (NewSeq number)">',
                        '##INFO=<ID=SLN,Number=1,Type=Integer,Description="Number of reads with fragments that fit exactly together (ScarLess number)">',
                        '##INFO=<ID=CN,Number=1,Type=Integer,Description="Number of reads with CIGAR">',
                        '##INFO=<ID=IVN,Number=1,Type=Integer,Description="Number of reads with inversions">',
                        '##INFO=<ID=NSP,Number=1,Type=Integer,Description="Value of new sequence homopolymeric tract pattern">',

                        '##INFO=<ID=MNF,Number=1,Type=Integer,Description="Median number of rearrangement fragments">',
                        '##INFO=<ID=PFF,Number=1,Type=Integer,Description="Percentage of first fragments">',
                        '##INFO=<ID=PLF,Number=1,Type=Integer,Description="Percentage of last fragments">',
                        
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
                        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele (reading depth, total coverage)">',
                        '##FORMAT=<ID=VF,Number=R,Type=Float,Description="Variant allele frequency">' 
                        ]
        for line in header_list_1:
            try:
                self.new_header.append(line)
            except ValueError:
                print('New header error! New_line: ', line)
                logging.error('New header error! New_line: '+line)
                break
        for name_contig, len_contig in self.general_dict_LR.items():
            len_contig = len_contig[0]
            line = '##contig=<ID='+str(name_contig)+',length='+str(len_contig)+'>'
            try:
                self.new_header.append(line)
            except ValueError:
                print('New header error! New_line: ', line)
                logging.error('New header error! New_line: '+line)
                break
        for line in header_list_2:
            try:
                self.new_header.append(line)
            except ValueError:
                print('New header error! New_line: ', line)
                logging.error('New header error! New_line: '+line)
                break
        self.new_header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+str(self.name_sample))
        
    # get sequences from reference genome 
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
                ff = pysam.FastaFile(self.reference_genome, filepath_index=self.reference_genome+'.fai')
                if diff <= 5:
                    diff -= 1
                    if diff >= 1:
                        continue
                    else:
                        # print()
                        # print('Error! Unable to extract the sequence from the reference genome!')
                        # print(chrom_name, temp_start, temp_end, diff, flag)
                        # sys.exit()
                        logging.error('Unable to extract the sequence from the reference genome! '+str(chrom_name)+' '+str(temp_start)+' '+str(temp_end)+' '+str(diff)+' '+str(flag))
                        return ''
        return (seq.upper())
    
    # find the best alignment
    def best_alignments(self, als):
        max_micro_seq = ''
        max_hom_seq = ''
        for num_al in range(len(als)):
            al = als[num_al]
            if al == ['']:
                continue
            # print('als', als)
            # print('al', al)
            # input()
            if len(al[0]) == 0:
                # micro
                if (num_al <= 2) and len(al[1]) > len(max_micro_seq) and len(al[1]) <= 3:
                    max_micro_seq = al[1]
                # homeo 
                temp_al = deepcopy(al)
                if len(temp_al) % 2 == 1:
                    temp_al = temp_al[:-1]
                while len(temp_al) > 1:
                    coml_areas = ''
                    nocompl_areas = ''
                    sequence = ''
                    for i in range(len(temp_al)):
                        if i == 0:
                            continue
                        if i % 2 == 1:
                            coml_areas += temp_al[i]
                        elif i % 2 == 0:
                            nocompl_areas += temp_al[i]
                        sequence += temp_al[i]
                    if len(coml_areas) / len(sequence) >= 0.8:
                        if len(sequence) > len(max_hom_seq) and len(sequence) > 3:
                            max_hom_seq = sequence
                    # print('temp_al', temp_al)
                    # print('len(temp_al)', len(temp_al))
                    # print('coml_areas', coml_areas)
                    # print('nocompl_areas', nocompl_areas)
                    # print('sequence', sequence)
                    # print('len(coml_areas) / len(sequence)', len(coml_areas) / len(sequence))
                    temp_al = temp_al[:-2]
                    # print('new_temp_al', temp_al)
                    # input()
                    if temp_al == 1:
                        break
        # print('!!!!!!')
        # print('max_micro_seq', max_micro_seq)
        # print('max_hom_seq', max_hom_seq)
        # print('!!!!!!')
        return (max_micro_seq, max_hom_seq)
    
    # local alignment between two sequences
    def local_alignment(self, seq1, seq2, side):
        align_result = []
        if side == 'R':
            seq1 = seq1[::-1]
            seq2 = seq2[::-1]
        align_0 = ['']
        align_1 = ['']
        align_1_rev = ['']
        align_2 = ['']
        align_2_rev = ['']
        align_3 = ['']
        align_3_rev = ['']
        for i in range(min(len(seq1), len(seq2))):
            i_1 = i + 1
            i_2 = i + 2
            i_3 = i + 3
            # no slide
            if seq1[i] != seq2[i] and len(align_0) % 2 == 1:
                align_0[-1] = align_0[-1] + seq1[i]
            elif seq1[i] != seq2[i] and len(align_0) % 2 == 0:
                align_0.append('')
                align_0[-1] = align_0[-1] + seq1[i]
            elif seq1[i] == seq2[i] and len(align_0) % 2 == 0:
                align_0[-1] = align_0[-1] + seq1[i]
            elif seq1[i] == seq2[i] and len(align_0) % 2 == 1:
                align_0.append('')
                align_0[-1] = align_0[-1] + seq1[i]
            # slide (1 nuc)
            if i_1 < min(len(seq1), len(seq2)):
                # slided seq2
                if seq1[i] != seq2[i_1] and len(align_1) % 2 == 1:
                    align_1[-1] = align_1[-1] + seq1[i]
                elif seq1[i] != seq2[i_1] and len(align_1) % 2 == 0:
                    align_1.append('')
                    align_1[-1] = align_1[-1] + seq1[i]
                elif seq1[i] == seq2[i_1] and len(align_1) % 2 == 0:
                    align_1[-1] = align_1[-1] + seq1[i]
                elif seq1[i] == seq2[i_1] and len(align_1) % 2 == 1:
                    align_1.append('')
                    align_1[-1] = align_1[-1] + seq1[i]
                # slided seq1
                if seq2[i] != seq1[i_1] and len(align_1_rev) % 2 == 1:
                    align_1_rev[-1] = align_1_rev[-1] + seq2[i]
                elif seq2[i] != seq1[i_1] and len(align_1_rev) % 2 == 0:
                    align_1_rev.append('')
                    align_1_rev[-1] = align_1_rev[-1] + seq2[i]
                elif seq2[i] == seq1[i_1] and len(align_1_rev) % 2 == 0:
                    align_1_rev[-1] = align_1_rev[-1] + seq2[i]
                elif seq2[i] == seq1[i_1] and len(align_1_rev) % 2 == 1:
                    align_1_rev.append('')
                    align_1_rev[-1] = align_1_rev[-1] + seq2[i]
            # slide (2 nuc)
            if i_2 < min(len(seq1), len(seq2)):
                # slided seq2
                if seq1[i] != seq2[i_2] and len(align_2) % 2 == 1:
                    align_2[-1] = align_2[-1] + seq1[i]
                elif seq1[i] != seq2[i_2] and len(align_2) % 2 == 0:
                    align_2.append('')
                    align_2[-1] = align_2[-1] + seq1[i]
                elif seq1[i] == seq2[i_2] and len(align_2) % 2 == 0:
                    align_2[-1] = align_2[-1] + seq1[i]
                elif seq1[i] == seq2[i_2] and len(align_2) % 2 == 1:
                    align_2.append('')
                    align_2[-1] = align_2[-1] + seq1[i]
                # slided seq1
                if seq2[i] != seq1[i_2] and len(align_2_rev) % 2 == 1:
                    align_2_rev[-1] = align_2_rev[-1] + seq2[i]
                elif seq2[i] != seq1[i_2] and len(align_2_rev) % 2 == 0:
                    align_2_rev.append('')
                    align_2_rev[-1] = align_2_rev[-1] + seq2[i]
                elif seq2[i] == seq1[i_2] and len(align_2_rev) % 2 == 0:
                    align_2_rev[-1] = align_2_rev[-1] + seq2[i]
                elif seq2[i] == seq1[i_2] and len(align_2_rev) % 2 == 1:
                    align_2_rev.append('')
                    align_2_rev[-1] = align_2_rev[-1] + seq2[i]
            # slide (3 nuc)
            if i_3 < min(len(seq1), len(seq2)):
                # slided seq2
                if seq1[i] != seq2[i_1] and len(align_3) % 2 == 1:
                    align_3[-1] = align_3[-1] + seq1[i]
                elif seq1[i] != seq2[i_1] and len(align_3) % 2 == 0:
                    align_3.append('')
                    align_3[-1] = align_3[-1] + seq1[i]
                elif seq1[i] == seq2[i_1] and len(align_3) % 2 == 0:
                    align_3[-1] = align_3[-1] + seq1[i]
                elif seq1[i] == seq2[i_1] and len(align_3) % 2 == 1:
                    align_3.append('')
                    align_3[-1] = align_3[-1] + seq1[i]
                # slided seq1
                if seq2[i] != seq1[i_1] and len(align_3_rev) % 2 == 1:
                    align_3_rev[-1] = align_3_rev[-1] + seq2[i]
                elif seq2[i] != seq1[i_1] and len(align_3_rev) % 2 == 0:
                    align_3_rev.append('')
                    align_3_rev[-1] = align_3_rev[-1] + seq2[i]
                elif seq2[i] == seq1[i_1] and len(align_3_rev) % 2 == 0:
                    align_3_rev[-1] = align_3_rev[-1] + seq2[i]
                elif seq2[i] == seq1[i_1] and len(align_3_rev) % 2 == 1:
                    align_3_rev.append('')
                    align_3_rev[-1] = align_3_rev[-1] + seq2[i]
        align_result = (align_0, align_1, align_1_rev, align_2, align_2_rev, align_3, align_3_rev)
        return align_result

    # search microhomology
    def find_microhomology(self, microseq_1_list, microseq_2_list, flag):
        microseq_max = ''
        homeoseq_max = ''
        # scores for the alignment
        if len(microseq_1_list) == 1 and len(microseq_2_list) == 1:
            microseq_1 = microseq_1_list[0]
            microseq_2 = microseq_2_list[0]
            alignments = self.local_alignment(microseq_1, microseq_2, flag)
            microseq_max, homeoseq_max = self.best_alignments(alignments)   
        elif len(microseq_1_list) != 1 and len(microseq_2_list) == 1:
            microseq_2 = microseq_2_list[0]
            for microseq_1 in microseq_1_list:
                alignments = self.local_alignment(microseq_1, microseq_2, flag)
                microseq, homeoseq = self.best_alignments(alignments)
                if len(microseq) > len(microseq_max):
                    microseq_max = microseq
                if len(homeoseq) > len(homeoseq_max):
                    homeoseq_max = homeoseq
        elif len(microseq_1_list) == 1 and len(microseq_2_list) != 1:
            microseq_1 = microseq_1_list[0]
            for microseq_2 in microseq_2_list:
                alignments = self.local_alignment(microseq_1, microseq_2, flag)
                microseq, homeoseq = self.best_alignments(alignments)   
                if len(microseq) > len(microseq_max):
                    microseq_max = microseq
                if len(homeoseq) > len(homeoseq_max):
                    homeoseq_max = homeoseq
        else:
            print('Error! There are several sequences in each list!')
            logging.error('Error! There are several sequences in each list!')
            sys.exit()
        return (microseq_max, homeoseq_max)

    # process one row from junctions file
    def process_row_LR(self, row):
        
        # save all information about genomic rearrangemnt
        info_LR = {}

        false_LGR = False

        info_LR['Chrom'] = row['Chrom1']
        info_LR['Chrom2'] = row['Chrom2']
        if row['Junction_Side1'] == 'L' or row['Junction_Side1'] == 'R':
            info_LR['Pos'] = int(row['Pos1']) + 1
            info_LR['Pos2'] = int(row['Pos2']) + 1
        info_LR['Strand1'] = row['Strand1']
        info_LR['Strand2'] = row['Strand2']
        info_LR['Junc_1'] = row['Junction_Side1']
        info_LR['Junc_2'] = row['Junction_Side2']
        if ',' in str(row['Read_Len1']):
            if 'D' in info_LR['Junc_1']:
                list_rl1 = row['Read_Len1'].split(',')
                info_LR['Read_Len1'] = int(round(median([int(x) for i,x in enumerate(list_rl1) if i % 2 == 0])))
                info_LR['Diff_reads'] = int(round(median([int(x) for i,x in enumerate(list_rl1) if i % 2 == 1]),0))
            else:
                info_LR['Read_Len1'] = max([int(x) for x in row['Read_Len1'].split(',')])
        else:
            info_LR['Read_Len1'] = int(row['Read_Len1'])
        info_LR['Read_Len2'] = int(row['Read_Len2'])
        info_LR['Read_Number'] = int(row['Read_Number'])
        
        line_num = int(row['Line_Num'])

        # save the length of InterSeq
        inter_seq_lens = []
        inter_seqs = []
        if row['Inter_Joints'] != 0 and row['Inter_Joints'] != '0':
            if ',' in row['Inter_Joints']:
                for inter_seq_len in str(row['Inter_Joints']).split(','):
                    try:
                        inter_seq_lens.append(int(inter_seq_len))
                    except ValueError:
                        inter_seq_lens.append(len(inter_seq_len))
                        inter_seqs.append(inter_seq_len)
            else:
                try:
                    inter_seq_lens.append(int(row['Inter_Joints']))
                except ValueError:
                    inter_seq_lens.append(len(row['Inter_Joints']))
                    inter_seqs.append(row['Inter_Joints'])
            info_LR['Inter_Most_Freq_Seq'] = int(round(median(inter_seq_lens)))
        else:
            info_LR['Inter_Most_Freq_Seq'] = 0
        
        if len(inter_seqs) > 0:
            info_LR['Inter_Seq'] = inter_seqs[0]
        else:
            info_LR['Inter_Seq'] = ''

        # save the length of NewSeq
        new_seq_lens = []
        if row['NewSeq_Joints'] != 0 and row['NewSeq_Joints'] != '0':
            if ',' in row['NewSeq_Joints']:
                for new_seq_len in str(row['NewSeq_Joints']).split(','):
                    try:
                        new_seq_lens.append(int(new_seq_len))
                    except ValueError:
                        new_seq_lens.append(len(new_seq_len))
            else:
                try:
                    new_seq_lens.append(int(row['NewSeq_Joints']))
                except ValueError:
                    new_seq_lens.append(len(row['NewSeq_Joints']))
            info_LR['NewSeq_Most_Freq_Seq'] = int(round(median(new_seq_lens)))
        else:
            info_LR['NewSeq_Most_Freq_Seq'] = 0

        # save the percent of errors
        try:
            read_muts = [float(x) for x in row['Read_Muts'].split(',')]
        except ValueError:
            print('ERROR of getting read muts:', row)
            print(row['Read_Muts'])
            print(row['Read_Muts'].split(','))
            logging.error('ERROR of getting read muts:'+str(row))
            sys.exit()
        info_LR['Read_Muts_median'] = int(round(median(read_muts),0))
        if len(read_muts) == 1:
            info_LR['Read_Muts_variance'] = 0
        else:
            info_LR['Read_Muts_variance'] = int(round(variance(read_muts),0))

        # additional parameters about line number, coverage and quality
        info_LR['LR_num'] = int(row['Line_Num'])
        # print("row['Ref_Coverage']", row['Ref_Coverage'])
        # print("row['Read_Number']", row['Read_Number'])
        info_LR['Ref_Coverage'] = max(0, int(row['Ref_Coverage']) - row['Read_Number'])
        # print("info_LR['Ref_Coverage']", info_LR['Ref_Coverage'])
        info_LR['MQ_median'] = int(row['MQ_median'])

        # ISN, SLN, NSN, CN, IVN
        info_LR['Inter_Num'] = len(inter_seq_lens)
        info_LR['NewSeq_Num'] = len(new_seq_lens)
        if ',' in str(row['Scarless_Joints']):
            info_LR['Scarless_Num'] = len(str(row['Scarless_Joints']).split(','))
        else:
            if row['Scarless_Joints'] == 0 or row['Scarless_Joints'] == '0':
                info_LR['Scarless_Num'] = 0  
            else:
                try:
                    info_LR['Scarless_Num'] = int(row['Scarless_Joints'])
                except ValueError:
                    info_LR['Scarless_Num'] = len(row['Scarless_Joints'])
        try:
            info_LR['Cigar_Num'] = int(row['Cigar_Num'])
        except ValueError:
            print('ERROR!', row)
            logging.error('ERROR!', row)
            sys.exit()
        info_LR['Inversion_Num'] = int(row['Inversion_Num'])

        # define sequences before and after break point
        try:
            if ',' in row['Seqs_after_BND']:
                int_rseq = row['Seqs_after_BND'].split(',')
            else:
                int_rseq = [row['Seqs_after_BND']]
        except TypeError:
            int_rseq = ['']
        try:
            if ',' in row['Seqs_before_BND']:
                int_lseq = row['Seqs_before_BND'].split(',')
            else:
                int_lseq = [row['Seqs_before_BND']]
        except TypeError:
            int_lseq = ['']

        # define numbers if reads are first/last/middle
        try:
            if type(row['Start_End_Locations']) == int:
                info_LR['New_Seq_Pattern'] = int(row['Start_End_Locations'])
            else:
                info_LR['New_Seq_Pattern'] = max([int(x) for x in str(row['Start_End_Locations']).split(',')])
            if type(row['Number_of_parts']) == int:
                info_LR['Number_of_parts'] = [int(row['Number_of_parts'])]
            else:
                info_LR['Number_of_parts'] = [int(x) for x in str(row['Number_of_parts']).split(',')]
            info_LR['MNF'] = round(median(info_LR['Number_of_parts']), 1)
            info_LR['PFF'] = round(info_LR['New_Seq_Pattern'].count(1) / len(info_LR['New_Seq_Pattern']), 1)
            info_LR['PLF'] = round(info_LR['New_Seq_Pattern'].count(2) / len(info_LR['New_Seq_Pattern']), 1)
        except AttributeError:
            info_LR['MNF'] = 0
            info_LR['PFF'] = 0
            info_LR['PLF'] = 0
        
        determined_type_num = 0

        # determine confirmed TDs, INVs, TRLs, DELs 
        # translocations
        if info_LR['Junc_1'] == 'T':
            info_LR['Type'] = 'TRL'
            # print('row', row)
            # where LGR was translocated
            try:
                info_LR['T_pos1_1'] = int(row['Pos1'].split(',')[0]) + 1
                info_LR['T_pos1_2'] = int(row['Pos1'].split(',')[1]) + 1
            except IndexError:
                print('Error in the record of coordinates for T!')
                print('row.Pos1:', row['Pos1'])
                logging('Error in the record of coordinates for T!')
                logging(str(row))
                sys.exit()
            # from where this piece was extracted
            try:
                info_LR['T_pos2_1'] = int(row['Pos2'].split(',')[0]) + 1
                info_LR['T_pos2_2'] = int(row['Pos2'].split(',')[1]) + 1
            except IndexError:
                print('Error in the record of coordinates for T!')
                print('row.Pos2:', row['Pos2'])
                logging('Error in the record of coordinates for T!')
                logging(str(row))
                sys.exit()
            # junction sides
            info_LR['Junc_1'] = 'R'
            # define length of TRL
            info_LR['LR_len'] = abs(info_LR['T_pos2_1'] - info_LR['T_pos2_2']) + 1
            # pos - where this piece was added
            info_LR['Pos'] = min(info_LR['T_pos1_1'], info_LR['T_pos1_2'])
            # pos2 - from were this piece was extracted
            info_LR['Pos2'] = min(info_LR['T_pos2_1'], info_LR['T_pos2_2'])
            # if pos2_1 > pos2_2 then it means "-" strand
            if info_LR['T_pos2_1'] > info_LR['T_pos2_2']:
                info_LR['Strand2'] = '-'
                info_LR['Junc_2'] = 'R'
            else:
                info_LR['Strand2'] = '+'
                info_LR['Junc_2'] = 'L'
            # if pos2_1 > pos2_2 then it means "-" strand

            # define distance between coords where this piece was added
            # if Dist1 < 0 then this is TRL with TD
            if info_LR['T_pos1_2'] - info_LR['T_pos1_1'] >= 0:
                info_LR['Dist1'] = info_LR['T_pos1_2'] - info_LR['T_pos1_1'] + 1
            else:
                info_LR['Dist1'] = info_LR['T_pos1_2'] - info_LR['T_pos1_1'] - 1
            # print('info_LR', info_LR)
            # input()
            determined_type_num += 1

        # inversion
        elif info_LR['Junc_1'] == 'I':
            info_LR['Type'] = 'INV'
            # pos1 (1 and 2) = gap between --> (1) and <-- (2)
            try:
                info_LR['I_pos1'] = int(row['Pos1'].split(',')[0]) + 1
                info_LR['I_pos2_1'] = int(row['Pos1'].split(',')[1]) + 1
            except IndexError:
                print('Error in the record of coordinates for I!')
                print('row.Pos1:', row['Pos1'])
                logging('Error in the record of coordinates for I!')
                logging(str(row))
                sys.exit()
            # pos2 (1 and 2) = gap between <-- (2) and --> (3)
            try:
                info_LR['I_pos2_2'] = int(row['Pos2'].split(',')[0]) + 1
                info_LR['I_pos3'] = int(row['Pos2'].split(',')[1]) + 1
            except IndexError:
                print('Error in the record of coordinates for I!')
                print('row.Pos2:', row['Pos2'])
                logging('Error in the record of coordinates for I!')
                logging(str(row))
                sys.exit()
            # junction sides
            info_LR['Junc_1'] = 'L'
            info_LR['Junc_2'] = 'R'
            info_LR['Strand1'] = '-'
            info_LR['Strand2'] = '-'
            # pos and pos2 = coordinates of inverted fragment of read
            info_LR['Pos'] = info_LR['I_pos2_1']
            info_LR['Pos2'] = info_LR['I_pos2_2']
            # # may be useless information
            # info_LR['INV_start'] = info_LR['Pos']
            # info_LR['INV_end'] = info_LR['Pos2']
            # length of LGR is length of interted fragment
            info_LR['LR_len'] = info_LR['Pos2'] - info_LR['Pos'] + 1
            # define distance of the first and the second gaps
            info_LR['Dist1'] = info_LR['I_pos2_1'] - info_LR['I_pos1']
            info_LR['Dist2'] = info_LR['I_pos3'] - info_LR['I_pos2_2']
            determined_type_num += 1

        # tandem duplication
        elif 'D' in info_LR['Junc_1']:
            info_LR['Type'] = 'TD'
            # pos1 (1 and 2) = junction between --> (1) and --> (2)
            try:
                info_LR['D_pos1'] = int(row['Pos1'].split(',')[0]) + 1
                info_LR['D_pos2_1'] = int(row['Pos1'].split(',')[1]) + 1
            except IndexError:
                print('Error in the record of coordinates for D!')
                print('row.Pos1:', row['Pos1'])
                logging('Error in the record of coordinates for D!')
                logging(str(row))
                sys.exit()
            # pos2 (1 and 2) = junction between --> (2) and --> (3)
            try:
                info_LR['D_pos2_2'] = int(row['Pos2'].split(',')[0]) + 1
                info_LR['D_pos3'] = int(row['Pos2'].split(',')[1]) + 1
            except IndexError:
                print('Error in the record of coordinates for D!')
                print('row.Pos2:', row['Pos2'])
                logging('Error in the record of coordinates for D!')
                logging(str(row))
                sys.exit()
            if len(info_LR['Junc_1'])>1:
                repeat_num=int(info_LR['Junc_1'].replace('D',''))+1
            else:
                repeat_num=2
            info_LR['TDRN']=repeat_num
            # junction sides
            info_LR['Junc_1'] = 'L'
            info_LR['Junc_2'] = 'R'
            info_LR['Strand1'] = '+'
            info_LR['Strand2'] = '+'
            # pos and pos2 = coordinates (start and end) of duplicated fragment of read
            info_LR['Pos'] = info_LR['D_pos2_1']
            info_LR['Pos2'] = max(info_LR['D_pos1'],info_LR['D_pos2_2'])
            # need to delete InterSeqLen (overlapping piece of fragments)
            if info_LR['Pos'] + info_LR['Inter_Most_Freq_Seq'] < info_LR['Pos2'] - 50:
                info_LR['Pos'] += info_LR['Inter_Most_Freq_Seq']
                # print('info_LR', info_LR)
            else:
                # print('FALSE TD', info_LR)
                false_LGR = True
            # length of LGR is length of duplicated fragment (at least 2 repeats of fragment)
            info_LR['LR_len'] = (info_LR['Pos2'] - info_LR['Pos'] + 1) * repeat_num
            # in this case we need to consider NewSeq as a part of TD length
            info_LR['LR_len'] += int(round(info_LR['NewSeq_Most_Freq_Seq'] * (info_LR['NewSeq_Num'] / info_LR['Read_Number'])))
            # minus difference between starts of fragments 
            if 'Diff_reads' in info_LR:
                info_LR['LR_len'] -= int(round(info_LR['Diff_reads']))
            # define distance between breaks (--> (1) and --> (3))
            info_LR['Dist1'] = info_LR['D_pos3'] - info_LR['D_pos1'] + 1
            # print('info_LR', sorted(info_LR.items()))
            determined_type_num += 1

        # inverted tandem duplication
        elif info_LR['Junc_1'] == 'V':
            info_LR['Type'] = 'INVTD'
            # pos1 (1 and 2) = junction between --> (1) and <-- (2)
            try:
                info_LR['V_pos1'] = int(row['Pos1'].split(',')[0]) + 1 # R+
                info_LR['V_pos2_2'] = int(row['Pos1'].split(',')[1]) + 1 # R-
            except IndexError:
                print('Error in the record of coordinates for V!')
                print('row.Pos1:', row['Pos1'])
                logging('Error in the record of coordinates for V!')
                logging(str(row))
                sys.exit()
            # pos2 (1 and 2) = junction between <-- (2) and --> (3)
            try:
                info_LR['V_pos2_1'] = int(row['Pos2'].split(',')[0]) + 1 # L-
                info_LR['V_pos3'] = int(row['Pos2'].split(',')[1]) + 1 # L+
            except IndexError:
                print('Error in the record of coordinates for V!')
                print('row.Pos2:', row['Pos2'])
                logging('Error in the record of coordinates for V!')
                logging(str(row))
                sys.exit()
            # junction sides
            info_LR['Junc_1'] = 'L'
            info_LR['Junc_2'] = 'R'
            info_LR['Strand1'] = '-'
            info_LR['Strand2'] = '-'
            # pos and pos2 = coordinates (start and end) of inverted and duplicated fragment of read
            info_LR['Pos'] = info_LR['V_pos2_1']
            info_LR['Pos2'] = info_LR['V_pos2_2']
            # length of LGR is length of duplicated fragment
            info_LR['LR_len'] = info_LR['Pos2'] - info_LR['Pos'] + 1
            # define distance between breaks (--> (1) and --> (3))
            info_LR['Dist1'] = info_LR['V_pos3'] - info_LR['V_pos1'] + 1
            determined_type_num += 1

        # not exact types of LGRs
        elif info_LR['Chrom'] == info_LR['Chrom2']:
            if info_LR['Strand1'] == info_LR['Strand2']: # ++ or --
                if info_LR['Junc_1'] == 'L' and info_LR['Junc_2'] == 'R': # L+R+ L-R-
                    if ((info_LR['Strand1'] == '+' and info_LR['Pos'] > info_LR['Pos2'] - info_LR['Read_Len2'] and 
                        abs((info_LR['Pos'] + info_LR['Read_Len1']) - info_LR['Pos2']) <= 100) or 
                        (info_LR['Strand1'] == '-' and info_LR['Pos'] + info_LR['Read_Len1'] < info_LR['Pos2'] and
                        abs((info_LR['Pos2'] - info_LR['Read_Len2']) - info_LR['Pos']) <= 100)):
                        info_LR['Type'] = 'TD'
                        # for this situation we need to add a length of interseq
                        if info_LR['Pos'] + int(info_LR['Inter_Most_Freq_Seq']) < info_LR['Pos2'] - 50:
                            info_LR['Pos'] += int(info_LR['Inter_Most_Freq_Seq'])
                        else:
                            false_LGR = True
                        info_LR['TDRN']=1
                        info_LR['LR_len'] = info_LR['Pos2'] - info_LR['Pos'] + 1 + info_LR['NewSeq_Most_Freq_Seq']
                        # it is neccessary that pos2 > pos1 
                        if info_LR['LR_len'] < 0:
                            print('Error! The rearrangement length is less than zero!')
                            print('Information dict:', info_LR)
                            logging.error('Error! The rearrangement length is less than zero!')
                            logging.error('Information dict: '+str(info_LR))
                            sys.exit()
                    else:
                        # may be TD or TRL
                        info_LR['Type'] = 'BND_TD'
                        # for this situation we need to add a length of interseq
                        if info_LR['Pos'] + int(info_LR['Inter_Most_Freq_Seq']) < info_LR['Pos2'] - 50:
                            info_LR['Pos'] += int(info_LR['Inter_Most_Freq_Seq'])
                        else:
                            false_LGR = True
                        info_LR['LR_len'] = info_LR['Pos2'] - info_LR['Pos'] + 1
                        # it is neccessary that pos2 > pos1 
                        if info_LR['LR_len'] < 0:
                            print('Error! The rearrangement length is less than zero!')
                            print('Information dict:', info_LR)
                            logging.error('Error! The rearrangement length is less than zero!')
                            logging.error('Information dict: '+str(info_LR))
                            sys.exit()
                    determined_type_num += 1
                elif info_LR['Junc_1'] == 'R' and info_LR['Junc_2'] == 'L': # R+L+ R-L-
                    # may be DEL
                    if info_LR['Cigar_Num'] != 0:
                        # more possibility that it is DEL
                        info_LR['Type'] = 'DEL'
                    else:
                        info_LR['Type'] = 'BND_DEL'
                    info_LR['LR_len'] = info_LR['Pos2'] - info_LR['Pos'] + 1
                    # it is neccessary that pos2 > pos1 
                    if info_LR['LR_len'] < 0:
                        print('Error! The rearrangement length is less than zero!')
                        print('Information dict:', info_LR)
                        logging.error('Error! The rearrangement length is less than zero!')
                        logging.error('Information dict: '+str(info_LR))
                        sys.exit()
                    determined_type_num += 1
            else: # +- or -+
                if info_LR['Junc_1'] == info_LR['Junc_2']: # RR or LL --> R+R- R-R+ L-L+ L+L-
                    if (info_LR['Junc_1'] == 'R' and info_LR['Pos2'] - info_LR['Read_Len2'] + 1 >= info_LR['Pos'] or 
                        info_LR['Junc_1'] == 'L' and info_LR['Pos'] + info_LR['Read_Len1'] -1 <= info_LR['Pos2']):
                        info_LR['Type'] = 'BND_INV'
                        # # may be useless
                        # if info_LR['Junc_1'] == 'R':
                        #     info_LR['INV_type'] = '+'
                        # else:
                        #     info_LR['INV_type'] = '-'
                        # info_LR['INV_start'] = info_LR['Pos']
                        # info_LR['INV_end'] = info_LR['Pos2']
                        info_LR['LR_len'] = info_LR['Pos2'] - info_LR['Pos'] + 1
                        if info_LR['LR_len'] < 0:
                            print('Error! The rearrangement length is less than zero!')
                            print('Information dict:', info_LR)
                            logging.error('Error! The rearrangement length is less than zero!')
                            logging.error('Information dict: '+str(info_LR))
                            sys.exit()
                        determined_type_num += 1
                    else: # also R+R- R-R+ L-L+ L+L- but with overlap
                        info_LR['Type'] = 'BND_INVTD'
                        # # may be useless
                        # info_LR['INVTD_start'] = info_LR['Pos']
                        # info_LR['INVTD_end'] = info_LR['Pos2']
                        if info_LR['Junc_1'] == 'L':
                            info_LR['LR_len'] = info_LR['Pos2'] + info_LR['Read_Len2'] - info_LR['Pos'] + 1
                        else:                            
                            info_LR['LR_len'] = info_LR['Pos2'] - info_LR['Pos'] + info_LR['Read_Len1'] + 1
                        if info_LR['LR_len'] < 0:
                            print('Error! The rearrangement length is less than zero!')
                            print('Information dict:', info_LR)
                            logging.error('Error! The rearrangement length is less than zero!')
                            logging.error('Information dict: '+str(info_LR))
                            sys.exit()
                        determined_type_num += 1
                        

        # BND between two chromosomes
        elif info_LR['Chrom'] != info_LR['Chrom2']:
            # may be translocation
            info_LR['Type'] = 'BND_TRL'
            # may be useless
            # first read direction from the position (UP - upstream, DOWN - downstream)
            # if info_LR['Junc_1'] == 'R':
            #     info_LR['SD_1'] = 'UP'
            # else:
            #     info_LR['SD_1'] = 'DOWN'
            # # second read direction from the position (UP - upstream, DOWN - downstream)
            # if info_LR['Junc_2'] == 'R':
            #     info_LR['SD_2'] = 'UP'
            # else:
            #     info_LR['SD_2'] = 'DOWN'
            determined_type_num += 1

        if false_LGR:
            # print('info_LR', info_LR)
            return {}

        if 'LR_len' in info_LR.keys():
            if info_LR['LR_len'] < 50:
                return {}

        # ID of LGR
        info_LR['ID'] = str(line_num) + info_LR['Junc_1'] + info_LR['Junc_2']

        # check if there if several types of the rearrangements
        if determined_type_num > 1:
            print('Error! There are several types of the rearrangement!')
            print('Information dict:', info_LR)
            logging.error('Error! There are several types of the rearrangement!')
            logging.error('Information dict: '+str(info_LR))
            sys.exit()

        # check microhomology for all types of LGRs
        info_LR['Micro'] = 0
        info_LR['MicroSeq'] = '""'
        info_LR['Homeo'] = 0
        info_LR['HomeoSeq'] = '""'

        # debag
        if info_LR['Inter_Num'] == 0 or (info_LR['Inter_Num'] != 0 and info_LR['Inter_Most_Freq_Seq'] == 0):
            if 'chr' in info_LR['Chrom'] and 'chr' in info_LR['Chrom2']:
                # for left checking
                # R+L+ L+R+
                if info_LR['Strand1'] == '+' and info_LR['Strand2'] == '+':
                    if info_LR['Junc_1'] == 'R' and info_LR['Junc_2'] == 'L':
                        seq_1_l = int_lseq
                        seq_2_l = (self.get_sequence_from_genome(info_LR['Chrom2'], info_LR['Pos2'], info_LR['Pos2'], 34, 'L'),)
                    elif info_LR['Junc_1'] == 'L' and info_LR['Junc_2'] == 'R':
                        seq_1_l = (self.get_sequence_from_genome(info_LR['Chrom'], info_LR['Pos'], info_LR['Pos'], 34, 'L'),)
                        seq_2_l = int_rseq
                # R-L- L-R-
                elif info_LR['Strand1'] == '-' and info_LR['Strand2'] == '-':
                    if info_LR['Junc_1'] == 'R' and info_LR['Junc_2'] == 'L':
                        seq_1_l = (self.get_sequence_from_genome(info_LR['Chrom'], info_LR['Pos'], info_LR['Pos'], 34, 'R'),)
                        seq_2_l = int_rseq
                    elif info_LR['Junc_1'] == 'L' and info_LR['Junc_2'] == 'R':
                        seq_1_l = int_lseq
                        seq_2_l = (self.get_sequence_from_genome(info_LR['Chrom2'], info_LR['Pos2'], info_LR['Pos2'], 34, 'L'),)
                # R+R- L+L-
                elif info_LR['Strand1'] == '+' and info_LR['Strand2'] == '-':
                    if info_LR['Junc_1'] == 'R' and info_LR['Junc_2'] == 'R':
                        seq_1_l = int_lseq
                        seq_2_l = (self.get_sequence_from_genome(info_LR['Chrom2'], info_LR['Pos2'], info_LR['Pos2'], 34, 'R'),)
                    elif info_LR['Junc_1'] == 'L' and info_LR['Junc_2'] == 'L':
                        seq_1_l = (self.complement_seq(self.get_sequence_from_genome(info_LR['Chrom'], info_LR['Pos'], info_LR['Pos'], 34, 'L')),)
                        seq_2_l = int_rseq
                 # L-L+ R-R+
                elif info_LR['Strand1'] == '-' and info_LR['Strand2'] == '+':
                    if info_LR['Junc_1'] == 'L' and info_LR['Junc_2'] == 'L':
                        seq_1_l = int_lseq
                        seq_2_l = (self.get_sequence_from_genome(info_LR['Chrom2'], info_LR['Pos2'], info_LR['Pos2'], 34, 'L'),)
                    elif info_LR['Junc_1'] == 'R' and info_LR['Junc_2'] == 'R':
                        seq_1_l = (self.complement_seq(self.get_sequence_from_genome(info_LR['Chrom'], info_LR['Pos'], info_LR['Pos'], 34, 'R')),)
                        seq_2_l = int_rseq
                try:
                    left_homology_info = self.find_microhomology(seq_1_l, seq_2_l, 'L')
                except UnboundLocalError:
                    print('info_LR', info_LR)
                    sys.exit()
                # for right checking
                # R+L+ L+R+
                if info_LR['Strand1'] == '+' and info_LR['Strand2'] == '+':
                    if info_LR['Junc_1'] == 'R' and info_LR['Junc_2'] == 'L':
                        seq_1_r = (self.get_sequence_from_genome(info_LR['Chrom'], info_LR['Pos'], info_LR['Pos'], 34, 'R'),)
                        seq_2_r = int_rseq
                    elif info_LR['Junc_1'] == 'L' and info_LR['Junc_2'] == 'R':
                        seq_1_r = int_lseq
                        seq_2_r = (self.get_sequence_from_genome(info_LR['Chrom2'], info_LR['Pos2'], info_LR['Pos2'], 34, 'R'),)
                # R-L- L-R-
                elif info_LR['Strand1'] == '-' and info_LR['Strand2'] == '-':
                    if info_LR['Junc_1'] == 'R' and info_LR['Junc_2'] == 'L':
                        seq_1_r = int_lseq
                        seq_2_r = (self.get_sequence_from_genome(info_LR['Chrom2'], info_LR['Pos2'], info_LR['Pos2'], 34, 'L'),)
                    elif info_LR['Junc_1'] == 'L' and info_LR['Junc_2'] == 'R':
                        seq_1_r = (self.get_sequence_from_genome(info_LR['Chrom'], info_LR['Pos'], info_LR['Pos'], 34, 'L'),)
                        seq_2_r = int_rseq
                # R+R- L+L-
                elif info_LR['Strand1'] == '+' and info_LR['Strand2'] == '-':
                    if info_LR['Junc_1'] == 'R' and info_LR['Junc_2'] == 'R':
                        seq_1_r = (self.complement_seq(self.get_sequence_from_genome(info_LR['Chrom'], info_LR['Pos'], info_LR['Pos'], 34, 'R')),)
                        seq_2_r = int_rseq
                    elif info_LR['Junc_1'] == 'L' and info_LR['Junc_2'] == 'L':
                        seq_1_r = int_lseq
                        seq_2_r = (self.get_sequence_from_genome(info_LR['Chrom2'], info_LR['Pos2'], info_LR['Pos2'], 34, 'L'),)
                # L-L+ R-R+
                elif info_LR['Strand1'] == '-' and info_LR['Strand2'] == '+':
                    if info_LR['Junc_1'] == 'L' and info_LR['Junc_2'] == 'L':
                        seq_1_r = (self.complement_seq(self.get_sequence_from_genome(info_LR['Chrom'], info_LR['Pos'], info_LR['Pos'], 34, 'L')),)
                        seq_2_r = int_rseq
                    elif info_LR['Junc_1'] == 'R' and info_LR['Junc_2'] == 'R':
                        seq_1_r = int_lseq
                        seq_2_r = (self.get_sequence_from_genome(info_LR['Chrom2'], info_LR['Pos2'], info_LR['Pos2'], 34, 'R'),)
                try:
                    right_homology_info = self.find_microhomology(seq_1_r, seq_2_r, 'R')
                except UnboundLocalError:
                    print('info_LR', info_LR)
                    sys.exit()
                microseq_l = left_homology_info[0]
                microseq_r = right_homology_info[0]
                if microseq_l != '' and microseq_r == '':
                    info_LR['Micro'] = len(microseq_l)
                    info_LR['MicroSeq'] = microseq_l
                elif microseq_l == '' and microseq_r != '':
                    info_LR['Micro'] = len(microseq_r)
                    info_LR['MicroSeq'] = microseq_r
                elif microseq_l == '' and microseq_r == '':
                    pass
                else:
                    if len(microseq_r) >= len(microseq_l):
                        info_LR['Micro'] = len(microseq_r)
                        info_LR['MicroSeq'] = microseq_r
                    else:
                        info_LR['Micro'] = len(microseq_l)
                        info_LR['MicroSeq'] = microseq_l
                # homeology
                homeoseq_l = left_homology_info[1]
                homeoseq_r = right_homology_info[1]
                if homeoseq_l != '' and homeoseq_r == '':
                    info_LR['Homeo'] = len(homeoseq_l)
                    info_LR['HomeoSeq'] = homeoseq_l
                elif homeoseq_l == '' and homeoseq_r != '':
                    info_LR['Homeo'] = len(homeoseq_r)
                    info_LR['HomeoSeq'] = homeoseq_r
                elif homeoseq_l == '' and homeoseq_r == '':
                    pass
                else:
                    if len(homeoseq_r) >= len(homeoseq_l):
                        info_LR['Homeo'] = len(homeoseq_r)
                        info_LR['HomeoSeq'] = homeoseq_r
                    else:
                        info_LR['Homeo'] = len(homeoseq_l)
                        info_LR['HomeoSeq'] = homeoseq_l
        elif info_LR['Inter_Most_Freq_Seq'] != 0:
            if info_LR['Inter_Most_Freq_Seq'] <= 3:
                info_LR['Micro'] = info_LR['Inter_Most_Freq_Seq']
                info_LR['MicroSeq'] = info_LR['Inter_Seq']
            else:
                info_LR['Homeo'] = info_LR['Inter_Most_Freq_Seq']
                if info_LR['Inter_Seq'] != '':
                    info_LR['HomeoSeq'] = len(info_LR['Inter_Seq'])
                else:
                    info_LR['HomeoSeq'] = info_LR['Inter_Most_Freq_Seq']

        # handling of possible errors in determining the LGR
        if determined_type_num == 0:
            print('Error! Failed to determine the type of the rearrangement!')
            print('Information dict:', info_LR)
            logging.error('Error! Failed to determine the type of the rearrangement!')
            logging.error('Information dict: '+str(info_LR))
            sys.exit()

        return info_LR

    # read the main csv-file with coordinates of rearrangements (all except insertions)
    def read_junc_file(self):
        print('Reading CSV-files with information about genomic rearrangements...')
        logging.info('Reading CSV-files with information about genomic rearrangements...')
        junc_file = pd.read_csv(self.path_to_junc_file, delimiter='\t')
        junc_file['Line_Num'] = range(1, len(junc_file) + 1)
        
        # delete lines with the number of reads less than 1
        junc_file = junc_file[junc_file['Read_Number'] >= 1]
        
        # DataFrame to the list of dicts
        rows = [row.to_dict() for _, row in junc_file.iterrows()]
        
        # create general data
        self.general_dict_LR_pos2 = deepcopy(self.general_dict_LR)
        
        # save chroms (their indexes will be numbers)
        self.chroms_to_numbers = []

        # create a dict with coordinates for checking LGRs located nearly
        self.types_secondary_boundry = {}
        self.sr_by_ids_secondary_boundry = {}

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
        # else:
        #     chunk_size_val = 100000
        # for row in rows:
        #     results = self.process_row_LR(row)
            
        results = list(tqdm(pool.imap_unordered(self.process_row_LR, rows, chunksize=chunk_size_val), total=len(rows), desc="Processing rows")) #, colour='blue'))
        
        pool.close()
        pool.join()
        
        # junside + strand
        # self.junc_strand_int = {'R+':11,
        #                         'R-':10,
        #                         'L+':1,
        #                         'L-':0}
        # join results
        for result in results:
            if result == {}:
                continue
            # save the information into a large dictionary for further sorting by chromosomes
            if result['Chrom'] not in self.general_dict_LR:
                self.general_dict_LR[result['Chrom']] = [0, []]
            self.general_dict_LR[result['Chrom']][1].append(result)

            # and fill another dictionary for checking repeats
            if result['Chrom2'] not in self.general_dict_LR_pos2:
                self.general_dict_LR_pos2[result['Chrom2']] = [0, []]
            self.general_dict_LR_pos2[result['Chrom2']][1].append(result)

            if result['Chrom'] not in self.chroms_to_numbers:
                self.chroms_to_numbers.append(result['Chrom'])
            if result['Chrom2'] not in self.chroms_to_numbers:
                self.chroms_to_numbers.append(result['Chrom2'])
            
            # save only such LGRs which can be secondary boundries
            if result['Type'] in ('BND_TD', 'BND_DEL', 'BND_INV', 'DEL', 'INV', 'TD'):
                type_sec = result['Type']
                if type_sec not in self.types_secondary_boundry:
                    # create internal lists
                    self.types_secondary_boundry[type_sec] = [
                        [],  # pos1
                        [],  # pos2
                        [],  # LR_num
                        [],  # chr1
                        [],  # chr2
                        [],  # juncside + strand int 
                        # [],  # len
                        [],  # (readlen1, readlen2)
                    ]
                self.types_secondary_boundry[type_sec][0].append(result['Pos'])
                self.types_secondary_boundry[type_sec][1].append(result['Pos2'])
                self.types_secondary_boundry[type_sec][2].append(result['LR_num'])
                self.types_secondary_boundry[type_sec][3].append(self.chroms_to_numbers.index(result['Chrom']))
                self.types_secondary_boundry[type_sec][4].append(self.chroms_to_numbers.index(result['Chrom2']))
                self.types_secondary_boundry[type_sec][5].append((result['Junc_1'], result['Strand1'], result['Junc_2'], result['Strand2']))
                self.types_secondary_boundry[type_sec][6].append((result['Read_Len1'], result['Read_Len2']))
                # if type_sec == 'BND_INV':
                #     print(self.types_secondary_boundry[type_sec])
                #     input()

            # save LR ID for defining number of reads

            if result['Type'] in ('BND_TD', 'BND_INV', 'INV', 'TRL', 'BND_DEL', 'TD', 'DEL'):
                self.sr_by_ids_secondary_boundry[result['LR_num']] = result['Read_Number']
        
        # print('self.types_secondary_boundry', self.types_secondary_boundry)
        # print('self.sr_by_ids_secondary_boundry', self.sr_by_ids_secondary_boundry)
        # print('self.chroms_to_numbers', self.chroms_to_numbers)
        # input()

        # delete data
        results = []
    
    # process every row in INS file
    def process_row_INS(self, row):
        
        # start_time = monotonic()
        # save information about insertions
        info_INS = {}
        info_INS['Type'] = 'INS'
        info_INS['Chrom'] = row['Chrom']
        info_INS['Pos'] = int(row['Pos']) + 1
        info_INS['LR_len'] = int(row['Insertion_Length'])
        info_INS['Read_Number'] = int(row['Read_Number'])
        info_INS['Total_coverage'] = int(row['Total_coverage'])
        info_INS['MQ_median'] = row['MQ_median']
        info_INS['ID'] = int(row['Line_Num'])
        
        if info_INS['LR_len'] < 50:
            return {}
        # end_time = monotonic()
        # print('1: ', end_time-start_time)
        
        # start_time = monotonic()
        # save the percent of errors
        read_muts = [int(round(float(x), 0)) for x in row['Read_Muts'].split(',')]
        info_INS['Read_Muts_median'] = median(read_muts)
        if len(read_muts) > 1:
            info_INS['Read_Muts_variance'] = int(round(variance(read_muts),0))
        else:
            info_INS['Read_Muts_variance'] = 0
        
        # end_time = monotonic()
        # print('2:', end_time-start_time)
        
        return info_INS
        
    # read a csv-file with information about insertions
    def read_file_INS(self):
        print()
        print('CSV-file reading with information about insertions...')
        logging.info('CSV-file reading with information about insertions...')
        ins_file = pd.read_csv(self.path_to_ins_file, delimiter='\t')
        ins_file['Line_Num'] = range(1, len(ins_file) + 1)
        # delete lines with the number of reads less than 1
        ins_file = ins_file[ins_file['Read_Number']>=1]
        
        # DataFrame to the list of dicts
        rows = [row.to_dict() for _, row in ins_file.iterrows()]
        # print('len(rows)', len(rows))

        # create additional parameters
        self.ins_chrom_dict = {}

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
        # else:
        #     chunk_size_val = 100000
        results = list(tqdm(pool.imap_unordered(self.process_row_INS, rows, chunksize=chunk_size_val), total=len(rows), desc="Processing rows")) #, colour='blue'))
        # results = pool.map(self.process_row_INS, rows)

        pool.close()
        pool.join()

        for result in results:

            if result == {}:
                continue

            if result['Chrom'] not in self.chroms_to_numbers:
                self.chroms_to_numbers.append(result['Chrom'])

            # save the information into a large dictionary for further sorting by chromosomes
            if result['Chrom'] not in self.general_dict_LR:
                self.general_dict_LR[result['Chrom']] = [0, []]
            self.general_dict_LR[result['Chrom']][1].append(result)
            # save information about insertions for checking connection with TRL
            if self.chroms_to_numbers.index(result['Chrom']) not in self.ins_chrom_dict.keys():
                self.ins_chrom_dict[self.chroms_to_numbers.index(result['Chrom'])]=[[],
                                                                                    [],
                                                                                    []]

            self.ins_chrom_dict[self.chroms_to_numbers.index(result['Chrom'])][0].append(result['Pos'])
            self.ins_chrom_dict[self.chroms_to_numbers.index(result['Chrom'])][1].append(result['LR_len'])
            self.ins_chrom_dict[self.chroms_to_numbers.index(result['Chrom'])][2].append(result['ID'])
            
        # print(self.ins_chrom_dict[self.chroms_to_numbers.index(result['Chrom'])][2])
        # input()
        
        # print('self.ins_chrom_dict', self.ins_chrom_dict)
        # input()

    # create temporary title for vcf-file
    def build_temp_header(self):
        self.temp_header = []
        header_list_1 = ['##fileformat=VCFv4.2',
                        '##fileDate='+str(datetime.strftime(datetime.now(), "%d-%m-%Y")), '##source=NAME_OF_SOURCE', 
                        '##source_reads='+self.path_to_bam_file, '##reference='+self.reference_genome, 
                        '##mapping=-', '##phasing=none', '##depth_of_coverage=-']
        header_list_2 = ['##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the rearrangement">',
                         '##INFO=<ID=POS1,Number=1,Type=Integer,Description="First break point">',
                         '##INFO=<ID=POS2,Number=1,Type=Integer,Description="Second break point">',
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
                        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele (reading depth, total coverage)">',
                        '##FORMAT=<ID=VF,Number=R,Type=Float,Description="Variant allele frequency">'
                        ]
        for line in header_list_1:
            try:
                self.temp_header.append(line)
            except ValueError:
                print('Error in the temporary header! New_line: ', line)
                logging.error('Error in the temporary header! New_line: '+line)
                sys.exit()
        for name_contig, len_contig in self.general_dict_LR.items():
            len_contig = len_contig[0]
            line = '##contig=<ID='+str(name_contig)+',length='+str(len_contig)+'>'
            try:
                self.temp_header.append(line)
            except ValueError:
                print('Error in the temporary header! New_line: ', line)
                logging.error('Error in the temporary header! New_line: '+line)
                sys.exit()
        for line in header_list_2:
            try:
                self.temp_header.append(line)
            except ValueError:
                print('Error in the temporary header! New_line: ', line)
                logging.error('Error in the temporary header! New_line: '+line)
                sys,exit()
        self.temp_header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+str(self.name_sample))

    # fill two temporary vcf-files
    def build_temp_vcf_data(self):
        # with first position
        with open(self.temp_vcf_file_pos1, 'a') as temp_vcf:
            for chrom, chrom_info in self.general_dict_LR.items():
                dicts_LGRs = chrom_info[1]
                for dict_LR_info in dicts_LGRs:
                    if dict_LR_info['Type'] != 'INS':
                        # for position_1 in range(dict_LR_info['Pos'] - 100, dict_LR_info['Pos'] + 101):
                        #     if position_1 < 1:
                        #         continue
                        new = '\t'.join([
                            str(dict_LR_info['Chrom']),
                            str(dict_LR_info['Pos'] - 1),
                            str(dict_LR_info['ID'])+'1',
                            'N',
                            str(dict_LR_info['Type']),
                            '0',
                            'PASS',
                            'SVLEN=1;POS1='+str(dict_LR_info['Pos'] - 1)+';POS2='+str(dict_LR_info['Pos2'] - 1),
                            'GT:DP:AD:VF',
                            '0'+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+'0.0'
                        ])
                        new = new + '\n'
                        temp_vcf.write(new)
        # with second position
        with open(self.temp_vcf_file_pos2, 'a') as temp_vcf:
            for chrom, chrom_info in self.general_dict_LR_pos2.items():
                dicts_LGRs = chrom_info[1]
                for dict_LR_info in dicts_LGRs:
                    if dict_LR_info['Type'] != 'INS':
                        # for position_2 in range(dict_LR_info['Pos2'] - 100, dict_LR_info['Pos2'] + 101):
                        #     if position_2 < 1:
                        #         continue
                        new = '\t'.join([
                            str(dict_LR_info['Chrom2']),
                            str(dict_LR_info['Pos2'] - 1),
                            str(dict_LR_info['ID'])+'2',
                            'N',
                            str(dict_LR_info['Type']),
                            '0',
                            'PASS',
                            'SVLEN=1;POS1='+str(dict_LR_info['Pos'] - 1)+';POS2='+str(dict_LR_info['Pos2'] - 1),
                            'GT:DP:AD:VF',
                            '0'+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+'0.0'
                        ])
                        new = new + '\n'
                        temp_vcf.write(new)

    # adjust bed-file (+- 200)
    def adjust_bed_coordinates(self, input_bed_path, output_bed_path):
        if input_bed_path[-3:]=='.gz':
            infile=gzip.open(input_bed_path,'rt')
        else:
            infile=open(input_bed_path,'rt')
        outfile=open(output_bed_path,'w')
        for line in infile:
            fields = line.strip().split('\t')
            chromosome = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name_el = fields[3]
            new_start = str(max(1, start - 200))
            new_end = str(end + 200)
            outfile.write(chromosome+'\t'+new_start+'\t'+new_end+'\t'+name_el+'\n')
       
    # create lua
    def create_lua_file(self, file_path):
        lua_content = """
    function mean(vals)
        local sum=0
        for i=1,#vals do
            sum = sum + vals[i]
        end
        return sum / #vals
    end

    function loc(chrom, start, stop)
        return chrom .. ":" .. start .. "-" .. stop
    end

    CLINVAR_LOOKUP = {}
    CLINVAR_LOOKUP['0'] = 'unknown'
    CLINVAR_LOOKUP['1'] = 'germline'
    CLINVAR_LOOKUP['2'] = 'somatic'
    CLINVAR_LOOKUP['4'] = 'inherited'
    CLINVAR_LOOKUP['8'] = 'paternal'
    CLINVAR_LOOKUP['16'] = 'maternal'
    CLINVAR_LOOKUP['32'] = 'de-novo'
    CLINVAR_LOOKUP['64'] = 'biparental'
    CLINVAR_LOOKUP['128'] = 'uniparental'
    CLINVAR_LOOKUP['256'] = 'not-tested'
    CLINVAR_LOOKUP['512'] = 'tested-inconclusive'
    CLINVAR_LOOKUP['1073741824'] = 'other'

    CLINVAR_SIG = {}
    CLINVAR_SIG['0'] = 'uncertain'
    CLINVAR_SIG['1'] = 'not-provided'
    CLINVAR_SIG['2'] = 'benign'
    CLINVAR_SIG['3'] = 'likely-benign'
    CLINVAR_SIG['4'] = 'likely-pathogenic'
    CLINVAR_SIG['5'] = 'pathogenic'
    CLINVAR_SIG['6'] = 'drug-response'
    CLINVAR_SIG['7'] = 'histocompatibility'
    CLINVAR_SIG['255'] = 'other'
    CLINVAR_SIG['.'] = '.'

    function intotbl(ud)
        local tbl = {}
        for i=1,#ud do
            tbl[i] = ud[i]
        end
        return tbl
    end

    -- from lua-users wiki
    function split(str, sep)
            local sep, fields = sep or ":", {}
            local pattern = string.format("([^%s]+)", sep)
            str:gsub(pattern, function(c) fields[#fields+1] = c end)
            return fields
    end

    function contains(str, tok)
        return string.find(str, tok) ~= nil
    end

    function div2(a, b)
        if(a == 0) then return "0.0" end
        return string.format("%.9f", (a + 0) / b)
    end

    function ratio(vals)
        vals = vals[1] -- get 2 values per element. ref and alt counts.
        if vals[2] == 0 then return "0.0" end
        return string.format("%.9f", vals[2] / (vals[1] + vals[2]))
    end

    function clinvar_sig(vals)
        local t = type(vals)
        -- just a single-value
        if(t == "string" or t == "number") and not contains(vals, "|") then
            return CLINVAR_SIG[vals]
        elseif t ~= "table" then
            if not contains(t, "userdata") then
                if t == "string" then
                    vals = split(vals, ",")
                else
                    vals = {vals}
                end
            else
                vals = intotbl(vals)
            end
        end
        local ret = {}
        for i=1,#vals do
            if not contains(vals[i], "|") then
                ret[#ret+1] = CLINVAR_SIG[vals[i]]
            else
                local invals = split(vals[i], "|")
                local inret = {}
                for j=1,#invals do
                    inret[#inret+1] = CLINVAR_SIG[invals[j]]
                end
                ret[#ret+1] = join(inret, "|")
            end
        end
        return join(ret, ",")
    end

    join = table.concat

    function check_clinvar_aaf(clinvar_sig, max_aaf_all, aaf_cutoff)
        -- didn't find an aaf for this so can't be common
        if max_aaf_all == nil or clinvar_sig == nil then
            return false
        end
        if type(clinvar_sig) ~= "string" then
            clinvar_sig = join(clinvar_sig, ",")
        end
        if false == contains(clinvar_sig, "pathogenic") then
            return false
        end
        if type(max_aaf_all) ~= "table" then
            return max_aaf_all > aaf_cutoff
        end
        for i, aaf in pairs(max_aaf_all) do
            if aaf > aaf_cutoff then
                return true
            end
        end
        return false
    end

    function setid(...)
        local t = {...}
        local res = {}
        local seen = {}
        for i, v in pairs(t) do
            if v ~= "." and v ~= nil and v ~= "" then
                if seen[v] == nil then
                    res[#res+1] = string.gsub(v, ",", ";")
                    seen[v] = true
                end
            end
        end
        return table.concat(res, ";")
    end
    """

        with open(file_path, 'w') as lua_file:
            lua_file.write(lua_content)
    
    # sorting and indexing bed file
    def process_bed_file(self, input_bed_path):
        sorted_bed_path = input_bed_path.replace('.bed', '_sorted.bed')
        sort_command = "bedtools sort -i "+input_bed_path+" > "+sorted_bed_path
        sp.run(sort_command, shell=True, check=True)
        self.compressed_bed_path = sorted_bed_path.replace('.bed', '.bed.gz')
        compress_command = "bgzip -c "+sorted_bed_path+" > "+self.compressed_bed_path
        sp.run(compress_command, shell=True, check=True)
        tabix_command = "tabix -p bed "+self.compressed_bed_path
        sp.run(tabix_command, shell=True, check=True)
    
    # create toml
    def create_toml_file(self, bed_file_path, toml_file_path):
        toml_content = {
            "annotation": [
                {
                    "file": bed_file_path,
                    "columns": [1, 2, 3, 4],
                    "names": ["chrom", "first_position", "last_position", "name"],
                    "ops": ["self", "self", "self", "self"]
                }
            ]
        }
        with open(toml_file_path, 'w') as toml_file:
            toml.dump(toml_content, toml_file)
    
    # to read joined bed-file
    def process_vcfanno_row(self, row):
        chrom_lr, pos1_lr, name_sides, x3, type, x5, x6, info, x7, x8 = row
        result_list = [0, [], # 0 - LERD, 1 - LERN
                        0, [], # 2 - LIRD, 3 - LIRN
                        [], # 4 - LCRN
                        0, [], # 5 - RERD, 6 - RERN
                        0, [], # 7 - RIRD, 8 - RIRN
                        []] # 9 - RCRN
        lr_num = int(name_sides[:-3])
        pos_num = int(name_sides[-1])
        if pos_num == 1:
            junc_pos1 = name_sides[-3]
            junc_pos2 = 0
        elif pos_num == 2:
            junc_pos2 = name_sides[-2]
            junc_pos1 = 0
        info_dict = dict(item.split('=') for item in info.split(';') if '=' in item)
        pos1_lr = int(info_dict['POS1'])
        pos2_lr = int(info_dict['POS2'])
        try:
            if ',' not in info_dict['first_position']:
                pos1_el_list = [int(info_dict['first_position'])+200]
                pos2_el_list = [int(info_dict['last_position'])-200]
                name_el_list = [info_dict['name']]
            else:
                pos1_el_list = [int(x)+200 for x in info_dict['first_position'].split(',')]
                pos2_el_list = [int(x)-200 for x in info_dict['last_position'].split(',')]
                name_el_list = [x for x in info_dict['name'].split(',')]
        except KeyError:
            return [lr_num] + result_list
        # print('name_sides', name_sides)
        # print('info', info)
        # print('name_el_list', name_el_list)
        # print('junc_pos1', junc_pos1)
        # print('junc_pos2', junc_pos2)
        # print('pos1_lr', pos1_lr)
        # print('pos2_lr', pos2_lr)
        # print('pos1_el_list', pos1_el_list)
        # print('pos2_el_list', pos2_el_list)
        # input()
        for i in range(len(name_el_list)):
            pos1_el = pos1_el_list[i]
            pos2_el = pos2_el_list[i]
            name_el = name_el_list[i]
            if junc_pos1 != 0:
                dist_pos1_right = pos1_el - pos1_lr + 1 
                if 0 < dist_pos1_right <= 200:
                    # external
                    if junc_pos1 == 'R':
                        result_list[0] += 1
                        result_list[1].append(name_el)
                    # internal
                    elif junc_pos1 == 'L':
                        result_list[2] += 1
                        result_list[3].append(name_el)
                elif dist_pos1_right < 0:
                    dist_pos1_right = pos2_el - pos1_lr + 1
                    dist_pos1_left = pos1_lr - pos1_el + 1
                    if dist_pos1_right > 0:
                        result_list[4].append(name_el)
                    else:
                        dist_pos1_left = pos1_lr - pos2_el + 1
                        if 0 < dist_pos1_left <= 200:
                            # external
                            if junc_pos1 == 'L':
                                result_list[0] += 1
                                result_list[1].append(name_el)
                            # internal
                            elif junc_pos1 == 'R':
                                result_list[2] += 1
                                result_list[3].append(name_el)         
            elif junc_pos2 != 0:
                dist_pos2_left = pos2_lr - pos2_el + 1
                if 0 < dist_pos2_left <= 200:
                    # external
                    if junc_pos2 == 'L':
                        result_list[5] += 1
                        result_list[6].append(name_el)
                    # internal
                    elif junc_pos2 == 'R':
                        result_list[7] += 1
                        result_list[8].append(name_el)
                elif dist_pos2_left < 0:
                    dist_pos2_right = pos2_el - pos2_lr + 1
                    dist_pos2_left = pos2_lr - pos1_el + 1
                    if dist_pos2_left > 0:
                        result_list[9].append(name_el)
                    else:
                        dist_pos2_right = pos1_el - pos2_lr + 1
                        if 0 < dist_pos2_right <= 200:
                            # internal
                            if junc_pos2 == 'L':
                                result_list[7] += 1
                                result_list[8].append(name_el)
                            # external
                            elif junc_pos2 == 'R':
                                result_list[5] += 1
                                result_list[6].append(name_el)

        return [lr_num] + result_list

    # to read results of analysis genomic elements
    def read_results_vcfanno(self):
        # analyze first file
        rows = []
        with open(self.temp_vcf_ann_1, 'r') as result_vcfanno:
            for row in result_vcfanno:
                if row.startswith('#'):
                    continue
                rows.append(row.strip().split('\t'))
        if len(rows) < 1000:
            chunk_size_val = 100
        elif len(rows) < 10000:
            chunk_size_val = 1000
        else:
            chunk_size_val = 10000
        # save data
        self.repeats_from_id = {}
        # create pool
        num_processes = self.threads
        pool = Pool(processes=num_processes)
        # for row_1 in rows:
        #     results = self.process_vcfanno_row(row_1)
        results_1 = list(tqdm(pool.imap_unordered(self.process_vcfanno_row, rows, chunksize=chunk_size_val), total=len(rows), desc="Processing vcfanno file 1")) #, colour='blue'))
        for result in results_1:
            id_LR = result[0]
            genomic_elements_info = result[1:]
            if id_LR not in self.repeats_from_id:
                self.repeats_from_id[id_LR] = [0, [],
                                               0, [],
                                               [],
                                               0, [],
                                               0, [],
                                               []]
            for ind in (1, 3, 4, 6, 8, 9):
                for name_el in genomic_elements_info[ind]:
                    # if name_el not in self.repeats_from_id[id_LR][ind]:
                    self.repeats_from_id[id_LR][ind].append(name_el)
                    if ind != 4 and ind != 9:
                        self.repeats_from_id[id_LR][ind-1] += 1
        results_1 = []
        pool.close()
        pool.join()
        # analyze second file
        rows = []
        with open(self.temp_vcf_ann_2, 'r') as result_vcfanno:
            for row in result_vcfanno:
                if row.startswith('#'):
                    continue
                rows.append(row.strip().split('\t'))
        if len(rows) < 1000:
            chunk_size_val = 100
        elif len(rows) < 10000:
            chunk_size_val = 1000
        else:
            chunk_size_val = 10000
        # create pool
        num_processes = self.threads
        pool = Pool(processes=num_processes)
        results_2 = list(tqdm(pool.imap_unordered(self.process_vcfanno_row, rows, chunksize=chunk_size_val), total=len(rows), desc="Processing vcfanno file 2")) #, colour='blue'))
        for result in results_2:
            id_LR = result[0]
            genomic_elements_info = result[1:]
            if id_LR not in self.repeats_from_id:
                self.repeats_from_id[id_LR] = [0, [],
                                               0, [],
                                               [],
                                               0, [],
                                               0, [],
                                               []]
            for ind in (1, 3, 4, 6, 8, 9):
                for name_el in genomic_elements_info[ind]:
                    # if name_el not in self.repeats_from_id[id_LR][ind]:
                    self.repeats_from_id[id_LR][ind].append(name_el)
                    if ind != 4 and ind != 9:
                        self.repeats_from_id[id_LR][ind-1] += 1 
        results_2 = []
        pool.close()
        pool.join()

    # to find possible secondary boundry
    def secondary_boundry(self, chrom_num_lr):
        
        current_chrom, current_num_d = chrom_num_lr
        dict_LR_info = self.general_dict_LR[current_chrom][1][current_num_d]
        if dict_LR_info['Type'] == 'INS':
            return (chrom_num_lr, )
        
        sblr_set = set()
        sbtrl_set = set()
        sbinv_set = set()
        sbtd_set = set()
        sb_lens = {}
        pos1_current = dict_LR_info['Pos']
        pos2_current = dict_LR_info['Pos2']
        current_type = dict_LR_info['Type']
        current_lr_num = dict_LR_info['LR_num']

        if current_type == 'BND_DEL':
            if 'BND_TD' in self.types_secondary_boundry.keys():
                for ind_second_boundry, second_boundry_pos1 in enumerate(self.types_secondary_boundry['BND_TD'][0]):
                    second_boundry_pos2 = self.types_secondary_boundry['BND_TD'][1][ind_second_boundry]
                    second_boundry_lr_num = self.types_secondary_boundry['BND_TD'][2][ind_second_boundry]
                    second_boundry_chrom1 = self.chroms_to_numbers[self.types_secondary_boundry['BND_TD'][3][ind_second_boundry]]
                    second_boundry_chrom2 = self.chroms_to_numbers[self.types_secondary_boundry['BND_TD'][4][ind_second_boundry]]
                    second_boundry_rl1 = self.types_secondary_boundry['BND_TD'][6][ind_second_boundry][0]
                    second_boundry_rl2 = self.types_secondary_boundry['BND_TD'][6][ind_second_boundry][1]
                    second_boundry_junc1, second_boundry_strand1, second_boundry_junc2, second_boundry_strand2 = self.types_secondary_boundry['BND_TD'][5][ind_second_boundry]
                    if set([dict_LR_info['Chrom'], dict_LR_info['Chrom2']]) != set([second_boundry_chrom1, second_boundry_chrom2]):
                        continue

                    if current_lr_num == second_boundry_lr_num:
                        continue

                    if ((abs(pos1_current - second_boundry_pos1) <= 200) or
                        (abs(pos2_current - second_boundry_pos2) <= 200)):
                        # sblr_set.add(second_boundry_lr_num)
                        strand_check = False
                        # len of trl
                        if (abs(pos1_current - second_boundry_pos1) <= 200):
                            # dict_LR_info['Strand1'] == second_boundry_strand1 and
                            if (pos2_current < second_boundry_pos2 and pos1_current < second_boundry_pos1 and
                                pos2_current < second_boundry_pos2 - second_boundry_rl2 and pos2_current + dict_LR_info['Read_Len2'] < second_boundry_pos2):
                                # print('BND_DEL', dict_LR_info)
                                # print('BND_TD', second_boundry_lr_num, 
                                #       second_boundry_chrom1,
                                #       second_boundry_pos1,
                                #       second_boundry_rl1,
                                #       second_boundry_junc1,
                                #       second_boundry_strand1,
                                #       second_boundry_chrom2,
                                #       second_boundry_pos2,
                                #       second_boundry_rl2,
                                #       second_boundry_junc2,
                                #       second_boundry_strand2)
                                # input() 
                                trl_len = second_boundry_pos2 - pos2_current + 1
                                strand_check = True
                            # elif dict_LR_info['Strand1'] != second_boundry_strand1 and pos2_current > second_boundry_pos2:
                            #    trl_len = pos2_current - second_boundry_pos2
                            #    strand_check = True
                        elif (abs(pos2_current - second_boundry_pos2) <= 200):
                            if (pos1_current > second_boundry_pos1 and pos2_current > second_boundry_pos2 and
                                second_boundry_pos1 + second_boundry_rl1 < pos1_current and pos1_current - dict_LR_info['Read_Len1'] > second_boundry_pos1): 
                                # print('BND_DEL', dict_LR_info)
                                # print('BND_TD', second_boundry_lr_num, 
                                #       second_boundry_chrom1,
                                #       second_boundry_pos1,
                                #       second_boundry_rl1,
                                #       second_boundry_junc1,
                                #       second_boundry_strand1,
                                #       second_boundry_chrom2,
                                #       second_boundry_pos2,
                                #       second_boundry_rl2,
                                #       second_boundry_junc2,
                                #       second_boundry_strand2)
                                # input()                               
                                trl_len = second_boundry_pos1 - pos1_current + 1
                                strand_check = True
                        if strand_check:
                            sbtrl_set.add(second_boundry_lr_num)
                            sb_lens[second_boundry_lr_num] = trl_len
                            sb_lens[current_lr_num] = trl_len

            if 'DEL' in self.types_secondary_boundry.keys():
                for ind_second_boundry, second_boundry_pos1 in enumerate(self.types_secondary_boundry['DEL'][0]):
                    second_boundry_pos2 = self.types_secondary_boundry['DEL'][1][ind_second_boundry]
                    second_boundry_lr_num = self.types_secondary_boundry['DEL'][2][ind_second_boundry]
                    second_boundry_chrom1 = self.chroms_to_numbers[self.types_secondary_boundry['DEL'][3][ind_second_boundry]]
                    second_boundry_chrom2 = self.chroms_to_numbers[self.types_secondary_boundry['DEL'][4][ind_second_boundry]]
                    second_boundry_junc1, second_boundry_strand1, second_boundry_junc2, second_boundry_strand2 = self.types_secondary_boundry['DEL'][5][ind_second_boundry]
                    if set([dict_LR_info['Chrom'], dict_LR_info['Chrom2']]) != set([second_boundry_chrom1, second_boundry_chrom2]):
                        continue
                    
                    if current_lr_num == second_boundry_lr_num:
                        continue

                    if (abs(pos1_current - second_boundry_pos1) <= 200) and (abs(pos2_current - second_boundry_pos2) <= 200):
                        sblr_set.add(second_boundry_lr_num)

        elif current_type == 'BND_INV':

            current_lr_len = dict_LR_info['LR_len']
            
            if not (dict_LR_info['Junc_1'] == 'R' and
                dict_LR_info['Junc_2'] == 'R' and
                dict_LR_info['Strand1'] != dict_LR_info['Strand2']):
                return (chrom_num_lr, )
            
            if 'BND_INV' in self.types_secondary_boundry.keys():
                for ind_second_boundry, second_boundry_pos1 in enumerate(self.types_secondary_boundry['BND_INV'][0]):
                    second_boundry_pos2 = self.types_secondary_boundry['BND_INV'][1][ind_second_boundry]
                    second_boundry_lr_num = self.types_secondary_boundry['BND_INV'][2][ind_second_boundry]
                    second_boundry_chrom1 = self.chroms_to_numbers[self.types_secondary_boundry['BND_INV'][3][ind_second_boundry]]
                    second_boundry_chrom2 = self.chroms_to_numbers[self.types_secondary_boundry['BND_INV'][4][ind_second_boundry]]
                    second_boundry_rl1 = self.types_secondary_boundry['BND_INV'][6][ind_second_boundry][0]
                    second_boundry_rl2 = self.types_secondary_boundry['BND_INV'][6][ind_second_boundry][1]
                    second_boundry_junc1, second_boundry_strand1, second_boundry_junc2, second_boundry_strand2 = self.types_secondary_boundry['BND_INV'][5][ind_second_boundry]

                    if set([dict_LR_info['Chrom'], dict_LR_info['Chrom2']]) != set([second_boundry_chrom1, second_boundry_chrom2]):
                        continue

                    if not (second_boundry_junc1 == 'L' and
                        second_boundry_junc2 == 'L' and
                        second_boundry_strand1 != second_boundry_strand2):
                        continue

                    if current_lr_num == second_boundry_lr_num:
                        continue

                    strand_check = False
                    if dict_LR_info['Strand1'] == second_boundry_strand1:
                        if (abs(pos1_current - second_boundry_pos1) <= 200 and second_boundry_pos2 + second_boundry_rl2 < pos2_current - dict_LR_info['Read_Len2']):
                            # len trl
                            trl_len = pos2_current - second_boundry_pos2 + 1
                            strand_check = True
                        elif (abs(pos2_current - second_boundry_pos2) <= 200 and second_boundry_pos1 + second_boundry_rl1 < pos1_current - dict_LR_info['Read_Len1']):
                            # len trl
                            trl_len = pos1_current - second_boundry_pos1 + 1
                            strand_check = True
                        if strand_check:
                            sbtrl_set.add(second_boundry_lr_num)
                            sb_lens[second_boundry_lr_num] = trl_len
                            sb_lens[current_lr_num] = trl_len
                    else:
                        if (pos1_current < second_boundry_pos1 and pos2_current < second_boundry_pos2 and 
                           abs(pos1_current - second_boundry_pos1) <= 200 and abs(pos2_current - second_boundry_pos2) <= 200):
                            inv_len = second_boundry_pos2 - pos1_current + 1
                            sbinv_set.add(second_boundry_lr_num)
                            sb_lens[second_boundry_lr_num] = inv_len
                            sb_lens[current_lr_num] = inv_len
                        
            if 'INV' in self.types_secondary_boundry.keys():
                for ind_second_boundry, second_boundry_pos1 in enumerate(self.types_secondary_boundry['INV'][0]):
                    second_boundry_pos2 = self.types_secondary_boundry['INV'][1][ind_second_boundry]
                    second_boundry_lr_num = self.types_secondary_boundry['INV'][2][ind_second_boundry]
                    second_boundry_chrom1 = self.chroms_to_numbers[self.types_secondary_boundry['INV'][3][ind_second_boundry]]
                    second_boundry_chrom2 = self.chroms_to_numbers[self.types_secondary_boundry['INV'][4][ind_second_boundry]]
                    second_junc_strand = self.types_secondary_boundry['INV'][5][ind_second_boundry]
                
                    if set([dict_LR_info['Chrom'], dict_LR_info['Chrom2']]) != set([second_boundry_chrom1, second_boundry_chrom2]):
                        continue

                    # if not(dict_LR_info['Strand1'] == second_boundry_strand1):
                    #     continue

                    if current_lr_num == second_boundry_lr_num:
                        continue

                    if (dict_LR_info['Strand1'] != second_boundry_strand1 and pos1_current < second_boundry_pos1 and pos2_current < second_boundry_pos2 and 
                       abs(pos1_current - second_boundry_pos1) <= 200 and abs(pos2_current - second_boundry_pos2) <= 200):
                        sblr_set.add(second_boundry_lr_num)

        elif current_type == 'TD':

            if 'TD' in self.types_secondary_boundry.keys():
                for ind_second_boundry, second_boundry_pos1 in enumerate(self.types_secondary_boundry['TD'][0]):
                    second_boundry_pos2 = self.types_secondary_boundry['TD'][1][ind_second_boundry]
                    second_boundry_lr_num = self.types_secondary_boundry['TD'][2][ind_second_boundry]
                    second_boundry_chrom1 = self.chroms_to_numbers[self.types_secondary_boundry['TD'][3][ind_second_boundry]]
                    second_boundry_chrom2 = self.chroms_to_numbers[self.types_secondary_boundry['TD'][4][ind_second_boundry]]
                    second_boundry_junc1, second_boundry_strand1, second_boundry_junc2, second_boundry_strand2 = self.types_secondary_boundry['TD'][5][ind_second_boundry]
                    if set([dict_LR_info['Chrom'], dict_LR_info['Chrom2']]) != set([second_boundry_chrom1, second_boundry_chrom2]):
                        continue

                    if current_lr_num == second_boundry_lr_num:
                        continue
                    
                    if not(second_boundry_strand1 == dict_LR_info['Strand1']):
                        continue

                    if (abs(pos1_current - second_boundry_pos1) <= 200 and abs(pos2_current - second_boundry_pos2) <= 200):
                        sblr_set.add(second_boundry_lr_num)
            
            if 'BND_TD' in self.types_secondary_boundry.keys():
                for ind_second_boundry, second_boundry_pos1 in enumerate(self.types_secondary_boundry['BND_TD'][0]):
                    second_boundry_pos2 = self.types_secondary_boundry['BND_TD'][1][ind_second_boundry]
                    second_boundry_lr_num = self.types_secondary_boundry['BND_TD'][2][ind_second_boundry]
                    second_boundry_chrom1 = self.chroms_to_numbers[self.types_secondary_boundry['BND_TD'][3][ind_second_boundry]]
                    second_boundry_chrom2 = self.chroms_to_numbers[self.types_secondary_boundry['BND_TD'][4][ind_second_boundry]]
                    second_boundry_rl1 = self.types_secondary_boundry['BND_TD'][6][ind_second_boundry][0]
                    second_boundry_rl2 = self.types_secondary_boundry['BND_TD'][6][ind_second_boundry][1]
                    second_boundry_junc1, second_boundry_strand1, second_boundry_junc2, second_boundry_strand2 = self.types_secondary_boundry['BND_TD'][5][ind_second_boundry]

                    if set([dict_LR_info['Chrom'], dict_LR_info['Chrom2']]) != set([second_boundry_chrom1, second_boundry_chrom2]):
                        continue

                    if current_lr_num == second_boundry_lr_num:
                        continue

                    if not(second_boundry_strand1 == dict_LR_info['Strand1']):
                        continue

                    if abs(pos1_current - second_boundry_pos1) <= 200 and abs(pos2_current - second_boundry_pos2) <= 200:
                        # sblr_set.add(second_boundry_lr_num)
                        td_len = max(pos2_current, second_boundry_pos2) - min(pos1_current, second_boundry_pos1) + 1
                        sbtd_set.add(second_boundry_lr_num)
                        sb_lens[second_boundry_lr_num] = td_len
                        sb_lens[current_lr_num] = td_len
            
        if current_type == 'TRL':

            if dict_LR_info['Strand1'] == dict_LR_info['Strand2']:

                if 'BND_TD' in self.types_secondary_boundry.keys():
                    for ind_second_boundry, second_boundry_pos1 in enumerate(self.types_secondary_boundry['BND_TD'][0]):
                        second_boundry_pos2 = self.types_secondary_boundry['BND_TD'][1][ind_second_boundry]
                        second_boundry_lr_num = self.types_secondary_boundry['BND_TD'][2][ind_second_boundry]
                        second_boundry_chrom1 = self.chroms_to_numbers[self.types_secondary_boundry['BND_TD'][3][ind_second_boundry]]
                        second_boundry_chrom2 = self.chroms_to_numbers[self.types_secondary_boundry['BND_TD'][4][ind_second_boundry]]
                        second_boundry_rl1 = self.types_secondary_boundry['BND_TD'][6][ind_second_boundry][0]
                        second_boundry_rl2 = self.types_secondary_boundry['BND_TD'][6][ind_second_boundry][1]
                        second_boundry_junc1, second_boundry_strand1, second_boundry_junc2, second_boundry_strand2 = self.types_secondary_boundry['BND_TD'][5][ind_second_boundry]

                        if set([dict_LR_info['Chrom'], dict_LR_info['Chrom2']]) != set([second_boundry_chrom1, second_boundry_chrom2]):
                            continue

                        if current_lr_num == second_boundry_lr_num:
                            continue

                        if (second_boundry_pos2 < max(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and 
                            second_boundry_pos2 - second_boundry_rl2 > min(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and
                            min(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2']) < second_boundry_pos1 < max(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2'])):
                            sblr_set.add(second_boundry_lr_num)
                        elif (second_boundry_pos1 > min(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and 
                            second_boundry_pos1 + second_boundry_rl1 < max(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and
                            min(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2']) < second_boundry_pos2 < max(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2'])):
                            sblr_set.add(second_boundry_lr_num)
                
                if 'BND_DEL' in self.types_secondary_boundry.keys():
                    for ind_second_boundry, second_boundry_pos1 in enumerate(self.types_secondary_boundry['BND_DEL'][0]):
                        second_boundry_pos2 = self.types_secondary_boundry['BND_DEL'][1][ind_second_boundry]
                        second_boundry_lr_num = self.types_secondary_boundry['BND_DEL'][2][ind_second_boundry]
                        second_boundry_chrom1 = self.chroms_to_numbers[self.types_secondary_boundry['BND_DEL'][3][ind_second_boundry]]
                        second_boundry_chrom2 = self.chroms_to_numbers[self.types_secondary_boundry['BND_DEL'][4][ind_second_boundry]]
                        second_boundry_rl1 = self.types_secondary_boundry['BND_DEL'][6][ind_second_boundry][0]
                        second_boundry_rl2 = self.types_secondary_boundry['BND_DEL'][6][ind_second_boundry][1]
                        second_boundry_junc1, second_boundry_strand1, second_boundry_junc2, second_boundry_strand2 = self.types_secondary_boundry['BND_DEL'][5][ind_second_boundry]

                        if set([dict_LR_info['Chrom'], dict_LR_info['Chrom2']]) != set([second_boundry_chrom1, second_boundry_chrom2]):
                            continue

                        if current_lr_num == second_boundry_lr_num:
                            continue

                        if (second_boundry_pos2 > min(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and 
                            second_boundry_pos2 + second_boundry_rl2 < max(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and
                            min(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2']) < second_boundry_pos1 < max(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2'])):
                            sblr_set.add(second_boundry_lr_num)
                        elif (second_boundry_pos1 < max(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and 
                            second_boundry_pos1 - second_boundry_rl1 > min(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and
                            min(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2']) < second_boundry_pos2 < max(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2'])):
                            sblr_set.add(second_boundry_lr_num)
            
            elif dict_LR_info['Strand1'] != dict_LR_info['Strand2']:

                if 'BND_INV' in self.types_secondary_boundry.keys():
                    for ind_second_boundry, second_boundry_pos1 in enumerate(self.types_secondary_boundry['BND_INV'][0]):
                        second_boundry_pos2 = self.types_secondary_boundry['BND_INV'][1][ind_second_boundry]
                        second_boundry_lr_num = self.types_secondary_boundry['BND_INV'][2][ind_second_boundry]
                        second_boundry_chrom1 = self.chroms_to_numbers[self.types_secondary_boundry['BND_INV'][3][ind_second_boundry]]
                        second_boundry_chrom2 = self.chroms_to_numbers[self.types_secondary_boundry['BND_INV'][4][ind_second_boundry]]
                        second_boundry_rl1 = self.types_secondary_boundry['BND_INV'][6][ind_second_boundry][0]
                        second_boundry_rl2 = self.types_secondary_boundry['BND_INV'][6][ind_second_boundry][1]
                        second_boundry_junc1, second_boundry_strand1, second_boundry_junc2, second_boundry_strand2 = self.types_secondary_boundry['BND_INV'][5][ind_second_boundry]

                        if set([dict_LR_info['Chrom'], dict_LR_info['Chrom2']]) != set([second_boundry_chrom1, second_boundry_chrom2]):
                            continue

                        if current_lr_num == second_boundry_lr_num:
                            continue
                        
                        if second_boundry_strand2 == '-':
                            # L-
                            if (second_boundry_pos2 > min(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and 
                                second_boundry_pos2 + second_boundry_rl2 < max(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and
                                min(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2']) < second_boundry_pos1 < max(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2'])):
                                sblr_set.add(second_boundry_lr_num)
                            # R -
                            elif (second_boundry_pos2 < max(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and 
                                second_boundry_pos2 - second_boundry_rl2 > min(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and
                                min(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2']) < second_boundry_pos1 < max(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2'])):
                                sblr_set.add(second_boundry_lr_num)
                        elif second_boundry_strand1 == '-':
                            # R -
                            if (second_boundry_pos1 < max(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and 
                                second_boundry_pos1 - second_boundry_rl1 > min(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and
                                min(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2']) < second_boundry_pos2 < max(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2'])):
                                sblr_set.add(second_boundry_lr_num)
                            # L -
                            elif (second_boundry_pos1 > min(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and 
                                second_boundry_pos1 + second_boundry_rl1 < max(dict_LR_info['T_pos2_1'], dict_LR_info['T_pos2_2']) and
                                min(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2']) < second_boundry_pos2 < max(dict_LR_info['T_pos1_1'], dict_LR_info['T_pos1_2'])):
                                sblr_set.add(second_boundry_lr_num)
                            

        return ((current_chrom, current_num_d, current_lr_num), sblr_set, sbtrl_set, sbtd_set, sbinv_set, sb_lens)
        
    # to find connected INS and TRL
    def connected_INS(self, chrom_num_lr):

        current_chrom, current_num_d = chrom_num_lr
        dict_LR_info = self.general_dict_LR[current_chrom][1][current_num_d]
        cins_set = set()
        current_lr_num = dict_LR_info['LR_num']
        index_current_chrom = self.chroms_to_numbers.index(dict_LR_info['Chrom'])
        index_current_chrom2 = self.chroms_to_numbers.index(dict_LR_info['Chrom2'])
                
        if dict_LR_info['Type'] == 'TRL':
            
            if index_current_chrom not in self.ins_chrom_dict.keys():
                # print('self.ins_chrom_dict.keys()', self.ins_chrom_dict.keys())
                # print('self.chroms_to_numbers', self.chroms_to_numbers)
                # print('index_current_chrom', index_current_chrom)
                # input()
                return (chrom_num_lr, )
            trl_len = dict_LR_info['LR_len'] + dict_LR_info['NewSeq_Most_Freq_Seq']            

            for ind_ins, ins_pos in enumerate(self.ins_chrom_dict[index_current_chrom][0]):
                ins_len = self.ins_chrom_dict[index_current_chrom][1][ind_ins]
                ins_name = self.ins_chrom_dict[index_current_chrom][2][ind_ins]
                
                if ((abs(ins_pos - dict_LR_info['T_pos1_1']) <= 200 or abs(ins_pos - dict_LR_info['T_pos1_2']) <= 200) and
                    (abs(ins_len - trl_len) / ins_len <= 0.1)):
                    cins_set.add(ins_name)
                # else:
                #     print(ins_pos, ins_len, ins_name)
        
        elif dict_LR_info['Type'] == 'BND_TRL':
            if index_current_chrom in self.ins_chrom_dict.keys():
                # We check for chrom1 and then for chrom2
                for ind_ins, ins_pos in enumerate(self.ins_chrom_dict[index_current_chrom][0]):
                    ins_len = self.ins_chrom_dict[index_current_chrom][1][ind_ins]
                    ins_name = self.ins_chrom_dict[index_current_chrom][2][ind_ins]
                    
                    if (abs(ins_pos - dict_LR_info['Pos']) <= 200 and
                        # the length of the insertion should be non less than the read fragment length from the chrom2
                        ins_len >= dict_LR_info['Read_Len2'] - 0.1 * ins_len):
                        cins_set.add(ins_name)
                    # else:
                    #     print(ins_pos, ins_len, ins_name)

            if index_current_chrom2 in self.ins_chrom_dict.keys():
                # Check if oppositely fragment was translocated from chrom1 to chrom2
                for ind_ins, ins_pos in enumerate(self.ins_chrom_dict[index_current_chrom2][0]):
                    ins_len = self.ins_chrom_dict[index_current_chrom2][1][ind_ins]
                    ins_name = self.ins_chrom_dict[index_current_chrom2][2][ind_ins]
                    
                    if (abs(ins_pos - dict_LR_info['Pos2']) <= 200 and
                        # the length of the insertion should be non less than the read fragment length from the chrom2
                        ins_len >= dict_LR_info['Read_Len1'] - 0.1 * ins_len):
                        cins_set.add(ins_name)
                    # else:
                    #     print(ins_pos, ins_len, ins_name)
        
        if dict_LR_info['Type'] == 'TD':
            
            if index_current_chrom not in self.ins_chrom_dict.keys():
                # print('self.ins_chrom_dict.keys()', self.ins_chrom_dict.keys())
                # print('self.chroms_to_numbers', self.chroms_to_numbers)
                # print('index_current_chrom', index_current_chrom)
                # input()
                return (chrom_num_lr, )
            
            td_len = dict_LR_info['LR_len']          

            for ind_ins, ins_pos in enumerate(self.ins_chrom_dict[index_current_chrom][0]):
                ins_len = self.ins_chrom_dict[index_current_chrom][1][ind_ins]
                ins_name = self.ins_chrom_dict[index_current_chrom][2][ind_ins]
                
                if ((abs(ins_pos - dict_LR_info['Pos']) <= 200 or abs(ins_pos - dict_LR_info['Pos2']) <= 200) and
                    (abs(ins_len - td_len) / ins_len <= 0.1)):
                    cins_set.add(ins_name)
                # else:
                #     print(ins_pos, ins_len, ins_name)

        return ((current_chrom, current_num_d, current_lr_num), cins_set)

    # to find connected LGRs
    def define_boundries(self):
        print()
        print('Trying to find secondary boundries and connected insertions for LGRs...')  
        logging.info('Trying to find secondary boundries and connected insertions for LGRs...') 
        # for saving "second" reads for LGRs
        self.secondary_LGRs_dict = {}
        self.secondary_TRLs_dict = {}
        self.secondary_TDs_dict = {}
        self.secondary_INVs_dict = {}
        self.connected_INS_dict = {}
        self.sb_lens_dict = {}

        # Split the data into chunks
        chunks_secondary_b = []
        for chrom, chrom_info in self.general_dict_LR.items():
            dicts_LR_info = chrom_info[1]
            for num_d, dict_LR_info in enumerate(dicts_LR_info):
                # if dict_LR_info['Type'] == 'TTD':
                #     print('TTD 1851')
                chunks_secondary_b.append((chrom, num_d))

        if len(chunks_secondary_b) < 1000:
            chunk_size_val = 100
        elif len(chunks_secondary_b) < 10000:
            chunk_size_val = 1000
        else:
            chunk_size_val = 10000
        # else:
        #     chunk_size_val = 100000

        # create pool
        num_processes = self.threads
        pool = Pool(processes=num_processes)
        
        # define possible secondary boundries
        results_secondary_boundary = list(tqdm(pool.imap_unordered(self.secondary_boundry, chunks_secondary_b, chunksize=chunk_size_val), total=len(chunks_secondary_b), desc="Processing LGRs"))
        
        # Split the data into chunks
        chunks_connected_ins = []
        for chrom, chrom_info in self.general_dict_LR.items():
            dicts_LR_info = chrom_info[1]
            for num_d, dict_LR_info in enumerate(dicts_LR_info):
                if dict_LR_info['Type'] in ('TRL', 'BND_TRL', 'TD'):
                    chunks_connected_ins.append((chrom, num_d))

        if len(chunks_connected_ins) < 1000:
            chunk_size_val = 100
        elif len(chunks_connected_ins) < 10000:
            chunk_size_val = 1000
        else:
            chunk_size_val = 10000
        # else:
        #     chunk_size_val = 100000
        
        # define possible connected TRL and INS
        results_connected_ins = list(tqdm(pool.imap_unordered(self.connected_INS, chunks_connected_ins, chunksize=chunk_size_val), total=len(chunks_connected_ins), desc="Processing TRLs and INS"))

        pool.close()
        pool.join()

        # combine the results
        for result in results_secondary_boundary:
            if len(result) == 1:
                chrom_name = result[0][0]
                num_d_lr = result[0][1]
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBLR'] = set()
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBTRL'] = set()
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBTD'] = set()
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBINV'] = set()
                continue
            chrom_name, num_d_lr, lr_num = result[0]
            if result[1] == set():
                chrom_name = result[0][0]
                num_d_lr = result[0][1]
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBLR'] = set()
            else:
                names_secondary_lrs = result[1]
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBLR'] = names_secondary_lrs
                for name_sec_lr in names_secondary_lrs:
                    if name_sec_lr not in self.secondary_LGRs_dict:
                        self.secondary_LGRs_dict[name_sec_lr] = set()
                    self.secondary_LGRs_dict[name_sec_lr].add(lr_num)
            if result[2] == set():
                chrom_name = result[0][0]
                num_d_lr = result[0][1]
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBTRL'] = set()
            else:
                names_secondary_trls = result[2]
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBTRL'] = names_secondary_trls
                for name_sec_trl in names_secondary_trls:
                    if name_sec_trl not in self.secondary_TRLs_dict:
                        self.secondary_TRLs_dict[name_sec_trl] = set()
                    self.secondary_TRLs_dict[name_sec_trl].add(lr_num)
            if result[3] == set():
                chrom_name = result[0][0]
                num_d_lr = result[0][1]
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBTD'] = set()
            else:
                names_secondary_tds = result[3]
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBTD'] = names_secondary_tds
                for name_sec_td in names_secondary_tds:
                    if name_sec_td not in self.secondary_TDs_dict:
                        self.secondary_TDs_dict[name_sec_td] = set()
                    self.secondary_TDs_dict[name_sec_td].add(lr_num)
            if result[4] == set():
                chrom_name = result[0][0]
                num_d_lr = result[0][1]
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBINV'] = set()
            else:
                names_secondary_invs = result[4]
                self.general_dict_LR[chrom_name][1][num_d_lr]['SBINV'] = names_secondary_invs
                for name_sec_inv in names_secondary_invs:
                    if name_sec_inv not in self.secondary_INVs_dict:
                        self.secondary_INVs_dict[name_sec_inv] = set()
                    self.secondary_INVs_dict[name_sec_inv].add(lr_num)
            for key, val in result[5].items():
                self.sb_lens_dict[key] = val      
        
        for result in results_connected_ins:
            if len(result) == 1 or result[1] == set():
                chrom_name = result[0][0]
                num_d_lr = result[0][1]
                self.general_dict_LR[chrom_name][1][num_d_lr]['CINS'] = set()
                continue
            chrom_name, num_d_lr, lr_num = result[0]
            names_connected_ins = result[1]
            self.general_dict_LR[chrom_name][1][num_d_lr]['CINS'] = names_connected_ins
            for name_ins in names_connected_ins:
                if name_ins not in self.connected_INS_dict:
                    self.connected_INS_dict[name_ins] = set()
                self.connected_INS_dict[name_ins].add(lr_num)
                                 
    # write all the information for each rearrangement into a vcf file
    def build_new_vcf_data(self):
        # # may be useless
        # inversion_before_trans = [] # list for saving possible inversion or inverted tandem duplications starts
        print()
        print('Creating data for VCF-file...')
        logging.info('Creating data for VCF-file...')
        with open(self.path_vcf_out, 'a') as vcf_file:
            # done_work = 0
            # all_work = sum([len(x[1][1]) for x in self.general_dict_LR.items()])
            # print(len(new_data))
            for key, val in tqdm(self.general_dict_LR.items(),  desc="Writing LGRs"): #, colour='blue'):
                # key = chrom 
                # val = list of rearrangements information
                val = val[1]
                for dict_LR_info in val:
                    pass_no = 'PASS'
                    lr_type_from_dict = '<'+dict_LR_info['Type']+'>'
                    if dict_LR_info['Type'] != 'INS':
                         # define repeats (or mobile elements) near to boundries of LGR
                        # print(self.repeats_from_id[dict_LR_info['LR_num']])
                        
                        prev_len = len(self.repeats_from_id[dict_LR_info['LR_num']][1])
                        now_len = len(set(self.repeats_from_id[dict_LR_info['LR_num']][1]))
                        diff_len = prev_len - now_len
                        dict_LR_info['LERD'] = self.repeats_from_id[dict_LR_info['LR_num']][0] - diff_len
                        dict_LR_info['LERN'] = ','.join(list(set(self.repeats_from_id[dict_LR_info['LR_num']][1])))
                        if dict_LR_info['LERN'] == '':
                            dict_LR_info['LERN'] = '""'

                        prev_len = len(self.repeats_from_id[dict_LR_info['LR_num']][3])
                        now_len = len(set(self.repeats_from_id[dict_LR_info['LR_num']][3]))
                        diff_len = prev_len - now_len
                        dict_LR_info['LIRD'] = self.repeats_from_id[dict_LR_info['LR_num']][2] - diff_len
                        dict_LR_info['LIRN'] = ','.join(list(set(self.repeats_from_id[dict_LR_info['LR_num']][3])))
                        if dict_LR_info['LIRN'] == '':
                            dict_LR_info['LIRN'] = '""'
                        
                        dict_LR_info['LCRN'] = ','.join(list(set(self.repeats_from_id[dict_LR_info['LR_num']][4])))
                        if dict_LR_info['LCRN'] == '':
                            dict_LR_info['LCRN'] = '""'
                        
                        prev_len = len(self.repeats_from_id[dict_LR_info['LR_num']][6])
                        now_len = len(set(self.repeats_from_id[dict_LR_info['LR_num']][6]))
                        diff_len = prev_len - now_len
                        dict_LR_info['RERD'] = self.repeats_from_id[dict_LR_info['LR_num']][5] - diff_len
                        dict_LR_info['RERN'] = ','.join(list(set(self.repeats_from_id[dict_LR_info['LR_num']][6])))
                        if dict_LR_info['RERN'] == '':
                            dict_LR_info['RERN'] = '""'
                        
                        prev_len = len(self.repeats_from_id[dict_LR_info['LR_num']][8])
                        now_len = len(set(self.repeats_from_id[dict_LR_info['LR_num']][8]))
                        diff_len = prev_len - now_len
                        dict_LR_info['RIRD'] = self.repeats_from_id[dict_LR_info['LR_num']][7] - diff_len
                        dict_LR_info['RIRN'] = ','.join(list(set(self.repeats_from_id[dict_LR_info['LR_num']][8])))
                        if dict_LR_info['RIRN'] == '':
                            dict_LR_info['RIRN'] = '""'
                        
                        dict_LR_info['RCRN'] = ','.join(list(set(self.repeats_from_id[dict_LR_info['LR_num']][9])))
                        if dict_LR_info['RCRN'] == '':
                            dict_LR_info['RCRN'] = '""'

                        # calculate QUAL
                        alt_cov_to_prob = min((1.27*dict_LR_info['Read_Number'])/(2.512+dict_LR_info['Read_Number']), 1)
                        vaf = dict_LR_info['Read_Number'] / (dict_LR_info['Read_Number'] + dict_LR_info['Ref_Coverage'])
                        vaf_to_prob = min((1.218*vaf)/(0.095+vaf), 1)
                        if dict_LR_info['MQ_median'] >= 60:
                            mq_to_prob = 1
                        else:
                            mq_median = dict_LR_info['MQ_median']
                            mq_to_prob = min(0.000282*(mq_median**2)-0.003384, 1)
                        qual = alt_cov_to_prob * vaf_to_prob * mq_to_prob

                        # if dict_LR_info['Type'] == 'DEL':
                        #     print()
                        #     print('qual', qual)
                        #     print('alt_cov_to_prob', alt_cov_to_prob)
                        #     print('vaf_to_prob', vaf_to_prob)
                        #     print('mq_to_prob', mq_to_prob)


                        # make a correction: if the distance between the reads is not zero, it reduces the quality
                        # if dict_LR_info['NewSeq_Num'] != 0:
                        #     if dict_LR_info['Inter_Num'] != 0 or dict_LR_info['Scarless_Num'] != 0:
                        #         coef = 0.2 / (dict_LR_info['NewSeq_Num'] / ((dict_LR_info['NewSeq_Num']) + (dict_LR_info['Inter_Num']) + (dict_LR_info['Scarless_Num'])))
                        #     else:
                        #         coef = 0.2
                        #     qual *= coef

                        # make correction: length of rearrangement should be proportional to min length of reads
                        if 'LR_len' in dict_LR_info.keys():
                            read_len_min = min(dict_LR_info['Read_Len1'], dict_LR_info['Read_Len2'])
                            if dict_LR_info['LR_len'] < 1000:
                                p_len = 1 / (1 + math.exp(-0.059096 * (read_len_min - 100))) - 0.00291
                            elif 1000 <= dict_LR_info['LR_len'] < 10000:
                                p_len = 1 / (1 + math.exp(-0.019699 * (read_len_min - 300))) - 0.00291
                            elif 10000 <= dict_LR_info['LR_len'] < 100000:
                                p_len = 1 / (1 + math.exp(-0.011819 * (read_len_min - 500))) - 0.00291
                            elif 100000 <= dict_LR_info['LR_len']:
                                p_len = 1 / (1 + math.exp(-0.005910 * (read_len_min - 1000))) - 0.00291
                            p_len = min(p_len, 1)
                            qual *= p_len
                        
                        # if dict_LR_info['Type'] == 'DEL':
                        #     print()
                        #     print('qual', qual)
                        #     print('p_len', p_len)

                        # correction for BND (TRL)
                        elif dict_LR_info['Type'] == 'BND_TRL': #or dict_LR_info['Type'] == 'TRL':
                            read_len_min = min(dict_LR_info['Read_Len1'], dict_LR_info['Read_Len2'])
                            p_len = 1 / (1 + math.exp(-0.005910 * (read_len_min - 1000))) - 0.00291
                            p_len = min(p_len, 1)
                            qual *= p_len

                        # correction for new sequence homopolymeric pattern
                        if dict_LR_info['New_Seq_Pattern']>0:
                            if type(dict_LR_info['NewSeq_Most_Freq_Seq'])==int:
                                val = dict_LR_info['New_Seq_Pattern']-dict_LR_info['NewSeq_Most_Freq_Seq']/100
                            else:
                                val = dict_LR_info['New_Seq_Pattern']-len(dict_LR_info['NewSeq_Most_Freq_Seq'])/100
                            if val>0:
                                qual/=2                        
                        
                        # increase qual for rearrangements with two boundries
                        if 'BND' not in dict_LR_info['Type']:
                            qual *= 1.5
                        
                        qual = min(100,round(qual*100, 2))

                        # if dict_LR_info['Type'] == 'DEL':
                        #     print()
                        #     print('qual', qual)
                        #     input()


                        vaf = round(vaf, 2) 
                        if qual <= 0:
                            qual = 0.0
                        if dict_LR_info['Ref_Coverage'] not in ['0', 0]:
                            genotype = '0/1'
                        else:
                            genotype = '1/1'
                        
                        # # may be useless
                        # # define possible inversion or inverted tandem duplications starts
                        # maybe_start = ''
                        # if dict_LR_info['Type'] == 'INV' and dict_LR_info['INV_type']=='+':
                        #     inv_type = 'START'
                        #     inversion_before_trans.append((dict_LR_info['INV_start'], dict_LR_info['INV_end'], qual))
                        # elif dict_LR_info['Type'] == 'INV' and dict_LR_info['INV_type']=='-':
                        #     inv_type = 'END'
                        #     inversion_before_trans = list(filter(lambda x: abs(x[1]-dict_LR_info['INV_start'])<= 50, inversion_before_trans))
                        #     if inversion_before_trans != []:
                        #         inversion_before_trans = list(sorted(inversion_before_trans, key=lambda x: x[2], reverse=True))
                        #         maybe_start = inversion_before_trans[0][0]
                        #     inversion_before_trans = []
                        # elif dict_LR_info['Type'] == 'INVTD' and inversion_before_trans != []:
                        #     inversion_before_trans_INVTD = list(filter(lambda x: abs(x[1]-dict_LR_info['INVTD_start'])<= 50, inversion_before_trans))
                        #     if inversion_before_trans_INVTD != []:
                        #         inversion_before_trans_INVTD = list(sorted(inversion_before_trans_INVTD, key=lambda x: x[2], reverse=True))
                        #         maybe_start = inversion_before_trans_INVTD[0][0]

                        
                        if dict_LR_info['LR_num'] in self.secondary_LGRs_dict.keys():
                            # print("self.secondary_LGRs_dict[dict_LR_info['LR_num']]", self.secondary_LGRs_dict[dict_LR_info['LR_num']])
                            # print("dict_LR_info['SBLR']", dict_LR_info['SBLR'])
                            if dict_LR_info['SBLR'] != set():
                                # print("dict_LR_info['SBLR']", dict_LR_info['SBLR'])
                                # print("self.secondary_LGRs_dict[dict_LR_info['LR_num']]", self.secondary_LGRs_dict[dict_LR_info['LR_num']])
                                for x in self.secondary_LGRs_dict[dict_LR_info['LR_num']]:
                                    dict_LR_info['SBLR'].add(x)
                                # print("dict_LR_info['SBLR']", dict_LR_info['SBLR'])
                                # print()
                            else:
                                dict_LR_info['SBLR'] = self.secondary_LGRs_dict[dict_LR_info['LR_num']]
                        
                        if dict_LR_info['LR_num'] in self.secondary_TRLs_dict.keys():
                            if dict_LR_info['SBTRL'] != set():
                                for x in self.secondary_TRLs_dict[dict_LR_info['LR_num']]:
                                    dict_LR_info['SBTRL'].add(x)
                            else:
                                dict_LR_info['SBTRL'] = self.secondary_TRLs_dict[dict_LR_info['LR_num']]
                        
                        if dict_LR_info['LR_num'] in self.secondary_TDs_dict.keys():
                            if dict_LR_info['SBTD'] != set():
                                for x in self.secondary_TDs_dict[dict_LR_info['LR_num']]:
                                    dict_LR_info['SBTD'].add(x)
                            else:
                                dict_LR_info['SBTD'] = self.secondary_TDs_dict[dict_LR_info['LR_num']]
                        
                        if dict_LR_info['LR_num'] in self.secondary_INVs_dict.keys():
                            if dict_LR_info['SBINV'] != set():
                                for x in self.secondary_INVs_dict[dict_LR_info['LR_num']]:
                                    dict_LR_info['SBINV'].add(x)
                            else:
                                dict_LR_info['SBINV'] = self.secondary_INVs_dict[dict_LR_info['LR_num']]

                        if dict_LR_info['Type'] == 'TRL' or dict_LR_info['Type'] == 'BND_TRL' or dict_LR_info['Type'] == 'TD':
                            if dict_LR_info['CINS'] != set():
                                # print("dict_LR_info['CINS']", dict_LR_info['CINS'])
                                for x in list(dict_LR_info['CINS']):
                                    dict_LR_info['SBLR'].add('INS'+str(x))
                                # print("dict_LR_info['SBLR']", dict_LR_info['SBLR'])

                        # PASS or FAIL
                        if qual == 0:
                            pass_no = 'FAIL'

                        if dict_LR_info['SBLR'] != set():
                            dict_LR_info['SBLR'] = ['LR'+str(x) if 'INS' not in str(x) else str(x) for x in list(dict_LR_info['SBLR'])]
                            # print("dict_LR_info['SBLR']", dict_LR_info['SBLR'])
                        else:
                            dict_LR_info['SBLR'] = ['""']
                        # trl
                        if dict_LR_info['SBTRL'] != set():
                            dict_LR_info['SBTRL'] = ['LR'+str(x)+':'+str(self.sb_lens_dict[x]) if 'INS' not in str(x) else str(x) for x in list(dict_LR_info['SBTRL'])]
                        else:
                            dict_LR_info['SBTRL'] = ['""']
                        # td
                        if dict_LR_info['SBTD'] != set():
                            dict_LR_info['SBTD'] = ['LR'+str(x)+':'+str(self.sb_lens_dict[x]) if 'INS' not in str(x) else str(x) for x in list(dict_LR_info['SBTD'])]
                        else:
                            dict_LR_info['SBTD'] = ['""']
                        # inv
                        if dict_LR_info['SBINV'] != set():
                            dict_LR_info['SBINV'] = ['LR'+str(x)+':'+str(self.sb_lens_dict[x]) if 'INS' not in str(x) else str(x) for x in list(dict_LR_info['SBINV'])]
                            # print("dict_LR_info['SBLR']", dict_LR_info['SBLR'])
                        else:
                            dict_LR_info['SBINV'] = ['""']

                        # # false positive
                        # if (dict_LR_info['Type'] == 'BND_INV' and dict_LR_info['LR_len'] < 1000 or
                        #     dict_LR_info['Type'] == 'BND_TD' and dict_LR_info['LR_len'] < 1000):
                        #     pass_no = 'FP'
                        # # true positive
                        # elif (dict_LR_info['Type'] in ('TD', 'TRL', 'INV') and dict_LR_info['Read_Muts_median'] < 10 
                        #     and dict_LR_info['SBLR'] != ['""'] and dict_LR_info['Read_Number'] == 1):
                        #     pass_no = 'TP'
                            
                        # debag
                        if dict_LR_info['Pos'] > dict_LR_info['Pos2'] and dict_LR_info['Type'] != 'TRL' and dict_LR_info['Type'] != 'BND_TRL':
                            print()
                            print("Warning! dict_LR_info['Pos'] > dict_LR_info['Pos2']")
                            print('dict_LR_info: ', dict_LR_info)
                            logging.warning("Warning! dict_LR_info['Pos'] > dict_LR_info['Pos2']")
                            logging.warning('dict_LR_info: '+dict_LR_info)
                            print()
                            # input()                              

                        # create new strings
                        if dict_LR_info['Type'] == 'BND_INV':

                            new = '\t'.join([
                                str(dict_LR_info['Chrom']),
                                str(dict_LR_info['Pos']),
                                'LR'+str(dict_LR_info['LR_num']),
                                'N',
                                str(lr_type_from_dict),
                                str(qual),
                                str(pass_no),
                                ('SVTYPE='+str(lr_type_from_dict)+';CHROM2='+str(dict_LR_info['Chrom2'])+';END='+str(dict_LR_info['Pos2'])+';SVLEN='+str(dict_LR_info['LR_len'])
                                +';J1='+str(dict_LR_info['Junc_1'])+';J2='+str(dict_LR_info['Junc_2'])+';S1='+str(dict_LR_info['Strand1'])+';S2='+str(dict_LR_info['Strand2'])
                                +';RL1='+str(dict_LR_info['Read_Len1'])+';RL2='+str(dict_LR_info['Read_Len2'])+';SR='+str(dict_LR_info['Read_Number'])+';MQ='+str(dict_LR_info['MQ_median'])
                                +';DP='+str(dict_LR_info['Ref_Coverage']+dict_LR_info['Read_Number'])+';VAF='+str(vaf)
                                +';MH='+str(dict_LR_info['Micro'])+';MHS='+str(dict_LR_info['MicroSeq'])+';HOM='+str(dict_LR_info['Homeo'])+';HOMS='+str(dict_LR_info['HomeoSeq'])
                                +';MUTM='+str(dict_LR_info['Read_Muts_median'])+';MUTV='+str(dict_LR_info['Read_Muts_variance'])
                                +';SBLR='+str(','.join(dict_LR_info['SBLR']))+';SBTRL='+str(','.join(dict_LR_info['SBTRL']))+';SBTD='+str(','.join(dict_LR_info['SBTD']))+';SBINV='+str(','.join(dict_LR_info['SBINV']))
                                +';LERD='+str(dict_LR_info['LERD'])+';LERN='+str(dict_LR_info['LERN'])+';LIRD='+str(dict_LR_info['LIRD'])+';LIRN='+str(dict_LR_info['LIRN'])+';LCRN='+str(dict_LR_info['LCRN'])
                                +';RERD='+str(dict_LR_info['RERD'])+';RERN='+str(dict_LR_info['RERN'])+';RIRD='+str(dict_LR_info['RIRD'])+';RIRN='+str(dict_LR_info['RIRN'])+';RCRN='+str(dict_LR_info['RCRN'])
                                +';ISN='+str(dict_LR_info['Inter_Num'])+';ISNMFS='+str(dict_LR_info['Inter_Most_Freq_Seq'])+';NSN='+str(dict_LR_info['NewSeq_Num'])
                                +';NSNMFS='+str(dict_LR_info['NewSeq_Most_Freq_Seq'])+';SLN='+str(dict_LR_info['Scarless_Num'])
                                +';CN='+str(dict_LR_info['Cigar_Num'])+';IVN='+str(dict_LR_info['Inversion_Num'])+';NSP='+str(dict_LR_info['New_Seq_Pattern'])
                                +';MNF='+str((dict_LR_info['MNF']))+';PFF='+str((dict_LR_info['PFF']))+';PLF='+str((dict_LR_info['PLF']))
                                ),
                                'GT:DP:AD:VF',
                                str(genotype)+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+str(vaf)
                            ])
                            new = new + '\n'
                        
                        elif dict_LR_info['Type'] == 'INV':

                            new = '\t'.join([
                                str(dict_LR_info['Chrom']),
                                str(dict_LR_info['Pos']),
                                'LR'+str(dict_LR_info['LR_num']),
                                'N',
                                str(lr_type_from_dict),
                                str(qual),
                                str(pass_no),
                                ('SVTYPE='+str(lr_type_from_dict)+';CHROM2='+str(dict_LR_info['Chrom2'])+';END='+str(dict_LR_info['Pos2'])+';SVLEN='+str(dict_LR_info['LR_len'])
                                +';J1='+str(dict_LR_info['Junc_1'])+';J2='+str(dict_LR_info['Junc_2'])+';S1='+str(dict_LR_info['Strand1'])+';S2='+str(dict_LR_info['Strand2'])
                                +';D1='+str(dict_LR_info['Dist1'])+';D2='+str(dict_LR_info['Dist2'])
                                +';RL1='+str(dict_LR_info['Read_Len1'])+';RL2='+str(dict_LR_info['Read_Len2'])+';SR='+str(dict_LR_info['Read_Number'])+';MQ='+str(dict_LR_info['MQ_median'])
                                +';DP='+str(dict_LR_info['Ref_Coverage']+dict_LR_info['Read_Number'])+';VAF='+str(vaf)
                                +';MH='+str(dict_LR_info['Micro'])+';MHS='+str(dict_LR_info['MicroSeq'])+';HOM='+str(dict_LR_info['Homeo'])+';HOMS='+str(dict_LR_info['HomeoSeq'])
                                +';MUTM='+str(dict_LR_info['Read_Muts_median'])+';MUTV='+str(dict_LR_info['Read_Muts_variance'])
                                +';SBLR='+str(','.join(dict_LR_info['SBLR']))+';SBTRL='+str(','.join(dict_LR_info['SBTRL']))+';SBTD='+str(','.join(dict_LR_info['SBTD']))+';SBINV='+str(','.join(dict_LR_info['SBINV']))
                                +';LERD='+str(dict_LR_info['LERD'])+';LERN='+str(dict_LR_info['LERN'])+';LIRD='+str(dict_LR_info['LIRD'])+';LIRN='+str(dict_LR_info['LIRN'])+';LCRN='+str(dict_LR_info['LCRN'])
                                +';RERD='+str(dict_LR_info['RERD'])+';RERN='+str(dict_LR_info['RERN'])+';RIRD='+str(dict_LR_info['RIRD'])+';RIRN='+str(dict_LR_info['RIRN'])+';RCRN='+str(dict_LR_info['RCRN'])
                                +';ISN='+str(dict_LR_info['Inter_Num'])+';ISNMFS='+str(dict_LR_info['Inter_Most_Freq_Seq'])+';NSN='+str(dict_LR_info['NewSeq_Num'])
                                +';NSNMFS='+str(dict_LR_info['NewSeq_Most_Freq_Seq'])+';SLN='+str(dict_LR_info['Scarless_Num'])
                                +';CN='+str(dict_LR_info['Cigar_Num'])+';IVN='+str(dict_LR_info['Inversion_Num'])+';NSP='+str(dict_LR_info['New_Seq_Pattern'])
                                +';MNF='+str((dict_LR_info['MNF']))+';PFF='+str((dict_LR_info['PFF']))+';PLF='+str((dict_LR_info['PLF']))
                                ),
                                'GT:DP:AD:VF',
                                str(genotype)+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+str(vaf)
                            ])
                            new = new + '\n'

                        elif dict_LR_info['Type'] == 'BND_INVTD':
                            
                            new = '\t'.join([
                                str(dict_LR_info['Chrom']),
                                str(dict_LR_info['Pos']),
                                'LR'+str(dict_LR_info['LR_num']),
                                'N',
                                str(lr_type_from_dict),
                                str(qual),
                                str(pass_no),
                                ('SVTYPE='+str(lr_type_from_dict)+';CHROM2='+str(dict_LR_info['Chrom2'])+';END='+str(dict_LR_info['Pos2'])+';SVLEN='+str(dict_LR_info['LR_len'])
                                +';J1='+str(dict_LR_info['Junc_1'])+';J2='+str(dict_LR_info['Junc_2'])+';S1='+str(dict_LR_info['Strand1'])+';S2='+str(dict_LR_info['Strand2'])
                                +';RL1='+str(dict_LR_info['Read_Len1'])+';RL2='+str(dict_LR_info['Read_Len2'])+';SR='+str(dict_LR_info['Read_Number'])+';MQ='+str(dict_LR_info['MQ_median'])
                                +';DP='+str(dict_LR_info['Ref_Coverage']+dict_LR_info['Read_Number'])+';VAF='+str(vaf)
                                +';MH='+str(dict_LR_info['Micro'])+';MHS='+str(dict_LR_info['MicroSeq'])+';HOM='+str(dict_LR_info['Homeo'])+';HOMS='+str(dict_LR_info['HomeoSeq'])
                                +';MUTM='+str(dict_LR_info['Read_Muts_median'])+';MUTV='+str(dict_LR_info['Read_Muts_variance'])
                                +';SBLR='+str(','.join(dict_LR_info['SBLR']))+';SBTRL='+str(','.join(dict_LR_info['SBTRL']))+';SBTD='+str(','.join(dict_LR_info['SBTD']))+';SBINV='+str(','.join(dict_LR_info['SBINV']))
                                +';LERD='+str(dict_LR_info['LERD'])+';LERN='+str(dict_LR_info['LERN'])+';LIRD='+str(dict_LR_info['LIRD'])+';LIRN='+str(dict_LR_info['LIRN'])+';LCRN='+str(dict_LR_info['LCRN'])
                                +';RERD='+str(dict_LR_info['RERD'])+';RERN='+str(dict_LR_info['RERN'])+';RIRD='+str(dict_LR_info['RIRD'])+';RIRN='+str(dict_LR_info['RIRN'])+';RCRN='+str(dict_LR_info['RCRN'])
                                +';ISN='+str(dict_LR_info['Inter_Num'])+';ISNMFS='+str(dict_LR_info['Inter_Most_Freq_Seq'])+';NSN='+str(dict_LR_info['NewSeq_Num'])
                                +';NSNMFS='+str(dict_LR_info['NewSeq_Most_Freq_Seq'])+';SLN='+str(dict_LR_info['Scarless_Num'])
                                +';CN='+str(dict_LR_info['Cigar_Num'])+';IVN='+str(dict_LR_info['Inversion_Num'])+';NSP='+str(dict_LR_info['New_Seq_Pattern'])
                                +';MNF='+str((dict_LR_info['MNF']))+';PFF='+str((dict_LR_info['PFF']))+';PLF='+str((dict_LR_info['PLF']))
                                ),
                                'GT:DP:AD:VF',
                                str(genotype)+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+str(vaf)
                            ])
                            new = new + '\n'
                        
                        elif dict_LR_info['Type'] == 'INVTD':
                            
                            new = '\t'.join([
                                str(dict_LR_info['Chrom']),
                                str(dict_LR_info['Pos']),
                                'LR'+str(dict_LR_info['LR_num']),
                                'N',
                                str(lr_type_from_dict),
                                str(qual),
                                str(pass_no),
                                ('SVTYPE='+str(lr_type_from_dict)+';CHROM2='+str(dict_LR_info['Chrom2'])+';END='+str(dict_LR_info['Pos2'])+';SVLEN='+str(dict_LR_info['LR_len'])
                                +';J1='+str(dict_LR_info['Junc_1'])+';J2='+str(dict_LR_info['Junc_2'])+';S1='+str(dict_LR_info['Strand1'])+';S2='+str(dict_LR_info['Strand2'])
                                +';D1='+str(dict_LR_info['Dist1'])
                                +';RL1='+str(dict_LR_info['Read_Len1'])+';RL2='+str(dict_LR_info['Read_Len2'])+';SR='+str(dict_LR_info['Read_Number'])+';MQ='+str(dict_LR_info['MQ_median'])
                                +';DP='+str(dict_LR_info['Ref_Coverage']+dict_LR_info['Read_Number'])+';VAF='+str(vaf)
                                +';MH='+str(dict_LR_info['Micro'])+';MHS='+str(dict_LR_info['MicroSeq'])+';HOM='+str(dict_LR_info['Homeo'])+';HOMS='+str(dict_LR_info['HomeoSeq'])
                                +';MUTM='+str(dict_LR_info['Read_Muts_median'])+';MUTV='+str(dict_LR_info['Read_Muts_variance'])
                                +';SBLR='+str(','.join(dict_LR_info['SBLR']))+';SBTRL='+str(','.join(dict_LR_info['SBTRL']))+';SBTD='+str(','.join(dict_LR_info['SBTD']))+';SBINV='+str(','.join(dict_LR_info['SBINV']))
                                +';LERD='+str(dict_LR_info['LERD'])+';LERN='+str(dict_LR_info['LERN'])+';LIRD='+str(dict_LR_info['LIRD'])+';LIRN='+str(dict_LR_info['LIRN'])+';LCRN='+str(dict_LR_info['LCRN'])
                                +';RERD='+str(dict_LR_info['RERD'])+';RERN='+str(dict_LR_info['RERN'])+';RIRD='+str(dict_LR_info['RIRD'])+';RIRN='+str(dict_LR_info['RIRN'])+';RCRN='+str(dict_LR_info['RCRN'])
                                +';ISN='+str(dict_LR_info['Inter_Num'])+';ISNMFS='+str(dict_LR_info['Inter_Most_Freq_Seq'])+';NSN='+str(dict_LR_info['NewSeq_Num'])
                                +';NSNMFS='+str(dict_LR_info['NewSeq_Most_Freq_Seq'])+';SLN='+str(dict_LR_info['Scarless_Num'])
                                +';CN='+str(dict_LR_info['Cigar_Num'])+';IVN='+str(dict_LR_info['Inversion_Num'])+';NSP='+str(dict_LR_info['New_Seq_Pattern'])
                                +';MNF='+str((dict_LR_info['MNF']))+';PFF='+str((dict_LR_info['PFF']))+';PLF='+str((dict_LR_info['PLF']))
                                ),
                                'GT:DP:AD:VF',
                                str(genotype)+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+str(vaf)
                            ])
                            new = new + '\n'
                        
                        elif dict_LR_info['Type'] == 'BND_DEL' or dict_LR_info['Type'] == 'DEL':
                            new = '\t'.join([
                                str(dict_LR_info['Chrom']),
                                str(dict_LR_info['Pos']),
                                'LR'+str(dict_LR_info['LR_num']),
                                'N',
                                str(lr_type_from_dict),
                                str(qual),
                                str(pass_no),
                                ('SVTYPE='+str(lr_type_from_dict)+';CHROM2='+str(dict_LR_info['Chrom2'])+';END='+str(dict_LR_info['Pos2'])+';SVLEN='+str(dict_LR_info['LR_len'])
                                +';J1='+str(dict_LR_info['Junc_1'])+';J2='+str(dict_LR_info['Junc_2'])+';S1='+str(dict_LR_info['Strand1'])+';S2='+str(dict_LR_info['Strand2'])
                                +';RL1='+str(dict_LR_info['Read_Len1'])+';RL2='+str(dict_LR_info['Read_Len2'])+';SR='+str(dict_LR_info['Read_Number'])+';MQ='+str(dict_LR_info['MQ_median'])
                                +';DP='+str(dict_LR_info['Ref_Coverage']+dict_LR_info['Read_Number'])+';VAF='+str(vaf)
                                +';MH='+str(dict_LR_info['Micro'])+';MHS='+str(dict_LR_info['MicroSeq'])+';HOM='+str(dict_LR_info['Homeo'])+';HOMS='+str(dict_LR_info['HomeoSeq'])
                                +';MUTM='+str(dict_LR_info['Read_Muts_median'])+';MUTV='+str(dict_LR_info['Read_Muts_variance'])
                                +';SBLR='+str(','.join(dict_LR_info['SBLR']))+';SBTRL='+str(','.join(dict_LR_info['SBTRL']))+';SBTD='+str(','.join(dict_LR_info['SBTD']))+';SBINV='+str(','.join(dict_LR_info['SBINV']))
                                +';LERD='+str(dict_LR_info['LERD'])+';LERN='+str(dict_LR_info['LERN'])+';LIRD='+str(dict_LR_info['LIRD'])+';LIRN='+str(dict_LR_info['LIRN'])+';LCRN='+str(dict_LR_info['LCRN'])
                                +';RERD='+str(dict_LR_info['RERD'])+';RERN='+str(dict_LR_info['RERN'])+';RIRD='+str(dict_LR_info['RIRD'])+';RIRN='+str(dict_LR_info['RIRN'])+';RCRN='+str(dict_LR_info['RCRN'])
                                +';ISN='+str(dict_LR_info['Inter_Num'])+';ISNMFS='+str(dict_LR_info['Inter_Most_Freq_Seq'])+';NSN='+str(dict_LR_info['NewSeq_Num'])
                                +';NSNMFS='+str(dict_LR_info['NewSeq_Most_Freq_Seq'])+';SLN='+str(dict_LR_info['Scarless_Num'])
                                +';CN='+str(dict_LR_info['Cigar_Num'])+';IVN='+str(dict_LR_info['Inversion_Num'])+';NSP='+str(dict_LR_info['New_Seq_Pattern'])
                                +';MNF='+str((dict_LR_info['MNF']))+';PFF='+str((dict_LR_info['PFF']))+';PLF='+str((dict_LR_info['PLF']))
                                ),
                                'GT:DP:AD:VF',
                                str(genotype)+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+str(vaf)
                            ])
                            new = new + '\n'

                        elif dict_LR_info['Type'] == 'BND_TRL':
                            new = '\t'.join([
                                str(dict_LR_info['Chrom']),
                                str(dict_LR_info['Pos']),
                                'LR'+str(dict_LR_info['LR_num']),
                                'N',
                                str(lr_type_from_dict),
                                str(qual),
                                str(pass_no),
                                ('SVTYPE='+str(lr_type_from_dict)+';CHROM2='+str(dict_LR_info['Chrom2'])+';END='+str(dict_LR_info['Pos2'])
                                +';J1='+str(dict_LR_info['Junc_1'])+';J2='+str(dict_LR_info['Junc_2'])+';S1='+str(dict_LR_info['Strand1'])+';S2='+str(dict_LR_info['Strand2'])
                                +';RL1='+str(dict_LR_info['Read_Len1'])+';RL2='+str(dict_LR_info['Read_Len2'])+';SR='+str(dict_LR_info['Read_Number'])+';MQ='+str(dict_LR_info['MQ_median'])
                                +';DP='+str(dict_LR_info['Ref_Coverage']+dict_LR_info['Read_Number'])+';VAF='+str(vaf)
                                +';MH='+str(dict_LR_info['Micro'])+';MHS='+str(dict_LR_info['MicroSeq'])+';HOM='+str(dict_LR_info['Homeo'])+';HOMS='+str(dict_LR_info['HomeoSeq'])
                                +';MUTM='+str(dict_LR_info['Read_Muts_median'])+';MUTV='+str(dict_LR_info['Read_Muts_variance'])
                                +';SBLR='+str(','.join(dict_LR_info['SBLR']))+';SBTRL='+str(','.join(dict_LR_info['SBTRL']))+';SBTD='+str(','.join(dict_LR_info['SBTD']))+';SBINV='+str(','.join(dict_LR_info['SBINV']))
                                +';LERD='+str(dict_LR_info['LERD'])+';LERN='+str(dict_LR_info['LERN'])+';LIRD='+str(dict_LR_info['LIRD'])+';LIRN='+str(dict_LR_info['LIRN'])+';LCRN='+str(dict_LR_info['LCRN'])
                                +';RERD='+str(dict_LR_info['RERD'])+';RERN='+str(dict_LR_info['RERN'])+';RIRD='+str(dict_LR_info['RIRD'])+';RIRN='+str(dict_LR_info['RIRN'])+';RCRN='+str(dict_LR_info['RCRN'])
                                +';ISN='+str(dict_LR_info['Inter_Num'])+';ISNMFS='+str(dict_LR_info['Inter_Most_Freq_Seq'])+';NSN='+str(dict_LR_info['NewSeq_Num'])
                                +';NSNMFS='+str(dict_LR_info['NewSeq_Most_Freq_Seq'])+';SLN='+str(dict_LR_info['Scarless_Num'])
                                +';CN='+str(dict_LR_info['Cigar_Num'])+';IVN='+str(dict_LR_info['Inversion_Num'])+';NSP='+str(dict_LR_info['New_Seq_Pattern'])
                                +';MNF='+str((dict_LR_info['MNF']))+';PFF='+str((dict_LR_info['PFF']))+';PLF='+str((dict_LR_info['PLF']))
                                ),
                                'GT:DP:AD:VF',
                                str(genotype)+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+str(vaf)
                            ])
                            new = new + '\n'
                        
                        elif dict_LR_info['Type'] == 'TRL':

                            new = '\t'.join([
                                str(dict_LR_info['Chrom']),
                                str(dict_LR_info['Pos']),
                                'LR'+str(dict_LR_info['LR_num']),
                                'N',
                                str(lr_type_from_dict),
                                str(qual),
                                str(pass_no),
                                ('SVTYPE='+str(lr_type_from_dict)+';CHROM2='+str(dict_LR_info['Chrom2'])+';END='+str(dict_LR_info['Pos2'])+';SVLEN='+str(dict_LR_info['LR_len'])
                                +';J1='+str(dict_LR_info['Junc_1'])+';J2='+str(dict_LR_info['Junc_2'])+';S1='+str(dict_LR_info['Strand1'])+';S2='+str(dict_LR_info['Strand2'])
                                +';D1='+str(dict_LR_info['Dist1'])
                                +';RL1='+str(dict_LR_info['Read_Len1'])+';RL2='+str(dict_LR_info['Read_Len2'])+';SR='+str(dict_LR_info['Read_Number'])+';MQ='+str(dict_LR_info['MQ_median'])
                                +';DP='+str(dict_LR_info['Ref_Coverage']+dict_LR_info['Read_Number'])+';VAF='+str(vaf)
                                +';MH='+str(dict_LR_info['Micro'])+';MHS='+str(dict_LR_info['MicroSeq'])+';HOM='+str(dict_LR_info['Homeo'])+';HOMS='+str(dict_LR_info['HomeoSeq'])
                                +';MUTM='+str(dict_LR_info['Read_Muts_median'])+';MUTV='+str(dict_LR_info['Read_Muts_variance'])
                                +';SBLR='+str(','.join(dict_LR_info['SBLR']))+';SBTRL='+str(','.join(dict_LR_info['SBTRL']))+';SBTD='+str(','.join(dict_LR_info['SBTD']))+';SBINV='+str(','.join(dict_LR_info['SBINV']))
                                +';LERD='+str(dict_LR_info['LERD'])+';LERN='+str(dict_LR_info['LERN'])+';LIRD='+str(dict_LR_info['LIRD'])+';LIRN='+str(dict_LR_info['LIRN'])+';LCRN='+str(dict_LR_info['LCRN'])
                                +';RERD='+str(dict_LR_info['RERD'])+';RERN='+str(dict_LR_info['RERN'])+';RIRD='+str(dict_LR_info['RIRD'])+';RIRN='+str(dict_LR_info['RIRN'])+';RCRN='+str(dict_LR_info['RCRN'])
                                +';ISN='+str(dict_LR_info['Inter_Num'])+';ISNMFS='+str(dict_LR_info['Inter_Most_Freq_Seq'])+';NSN='+str(dict_LR_info['NewSeq_Num'])
                                +';NSNMFS='+str(dict_LR_info['NewSeq_Most_Freq_Seq'])+';SLN='+str(dict_LR_info['Scarless_Num'])
                                +';CN='+str(dict_LR_info['Cigar_Num'])+';IVN='+str(dict_LR_info['Inversion_Num'])+';NSP='+str(dict_LR_info['New_Seq_Pattern'])
                                +';MNF='+str((dict_LR_info['MNF']))+';PFF='+str((dict_LR_info['PFF']))+';PLF='+str((dict_LR_info['PLF']))
                                ),
                                'GT:DP:AD:VF',
                                str(genotype)+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+str(vaf)
                            ])
                            new = new + '\n'
                        
                        elif dict_LR_info['Type'] == 'BND_TD': #or dict_LR_info['Type'] == 'TD':
                            
                            new = '\t'.join([
                                str(dict_LR_info['Chrom']),
                                str(dict_LR_info['Pos']),
                                'LR'+str(dict_LR_info['LR_num']),
                                'N',
                                str(lr_type_from_dict),
                                str(qual),
                                str(pass_no),
                                ('SVTYPE='+str(lr_type_from_dict)+';CHROM2='+str(dict_LR_info['Chrom2'])+';END='+str(dict_LR_info['Pos2'])+';SVLEN='+str(dict_LR_info['LR_len'])
                                +';J1='+str(dict_LR_info['Junc_1'])+';J2='+str(dict_LR_info['Junc_2'])+';S1='+str(dict_LR_info['Strand1'])+';S2='+str(dict_LR_info['Strand2'])
                                +';RL1='+str(dict_LR_info['Read_Len1'])+';RL2='+str(dict_LR_info['Read_Len2'])+';SR='+str(dict_LR_info['Read_Number'])+';MQ='+str(dict_LR_info['MQ_median'])
                                +';DP='+str(dict_LR_info['Ref_Coverage']+dict_LR_info['Read_Number'])+';VAF='+str(vaf)
                                +';MH='+str(dict_LR_info['Micro'])+';MHS='+str(dict_LR_info['MicroSeq'])+';HOM='+str(dict_LR_info['Homeo'])+';HOMS='+str(dict_LR_info['HomeoSeq'])
                                +';MUTM='+str(dict_LR_info['Read_Muts_median'])+';MUTV='+str(dict_LR_info['Read_Muts_variance'])
                                +';SBLR='+str(','.join(dict_LR_info['SBLR']))+';SBTRL='+str(','.join(dict_LR_info['SBTRL']))+';SBTD='+str(','.join(dict_LR_info['SBTD']))+';SBINV='+str(','.join(dict_LR_info['SBINV']))
                                +';LERD='+str(dict_LR_info['LERD'])+';LERN='+str(dict_LR_info['LERN'])+';LIRD='+str(dict_LR_info['LIRD'])+';LIRN='+str(dict_LR_info['LIRN'])+';LCRN='+str(dict_LR_info['LCRN'])
                                +';RERD='+str(dict_LR_info['RERD'])+';RERN='+str(dict_LR_info['RERN'])+';RIRD='+str(dict_LR_info['RIRD'])+';RIRN='+str(dict_LR_info['RIRN'])+';RCRN='+str(dict_LR_info['RCRN'])
                                +';ISN='+str(dict_LR_info['Inter_Num'])+';ISNMFS='+str(dict_LR_info['Inter_Most_Freq_Seq'])+';NSN='+str(dict_LR_info['NewSeq_Num'])
                                +';NSNMFS='+str(dict_LR_info['NewSeq_Most_Freq_Seq'])+';SLN='+str(dict_LR_info['Scarless_Num'])
                                +';CN='+str(dict_LR_info['Cigar_Num'])+';IVN='+str(dict_LR_info['Inversion_Num'])+';NSP='+str(dict_LR_info['New_Seq_Pattern'])
                                +';MNF='+str((dict_LR_info['MNF']))+';PFF='+str((dict_LR_info['PFF']))+';PLF='+str((dict_LR_info['PLF']))
                                ),
                                'GT:DP:AD:VF',
                                str(genotype)+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+str(vaf)
                            ])
                            new = new + '\n'

                        elif dict_LR_info['Type'] == 'TD':
                            
                            new = '\t'.join([
                                str(dict_LR_info['Chrom']),
                                str(dict_LR_info['Pos']),
                                'LR'+str(dict_LR_info['LR_num']),
                                'N',
                                str(lr_type_from_dict),
                                str(qual),
                                str(pass_no),
                                ('SVTYPE='+str(lr_type_from_dict)+';CHROM2='+str(dict_LR_info['Chrom2'])+';END='+str(dict_LR_info['Pos2'])+';SVLEN='+str(dict_LR_info['LR_len'])+';TDRN='+str(dict_LR_info['TDRN'])
                                +';J1='+str(dict_LR_info['Junc_1'])+';J2='+str(dict_LR_info['Junc_2'])+';S1='+str(dict_LR_info['Strand1'])+';S2='+str(dict_LR_info['Strand2'])
                                +';RL1='+str(dict_LR_info['Read_Len1'])+';RL2='+str(dict_LR_info['Read_Len2'])+';SR='+str(dict_LR_info['Read_Number'])+';MQ='+str(dict_LR_info['MQ_median'])
                                +';DP='+str(dict_LR_info['Ref_Coverage']+dict_LR_info['Read_Number'])+';VAF='+str(vaf)
                                +';MH='+str(dict_LR_info['Micro'])+';MHS='+str(dict_LR_info['MicroSeq'])+';HOM='+str(dict_LR_info['Homeo'])+';HOMS='+str(dict_LR_info['HomeoSeq'])
                                +';MUTM='+str(dict_LR_info['Read_Muts_median'])+';MUTV='+str(dict_LR_info['Read_Muts_variance'])
                                +';SBLR='+str(','.join(dict_LR_info['SBLR']))+';SBTRL='+str(','.join(dict_LR_info['SBTRL']))+';SBTD='+str(','.join(dict_LR_info['SBTD']))+';SBINV='+str(','.join(dict_LR_info['SBINV']))
                                +';LERD='+str(dict_LR_info['LERD'])+';LERN='+str(dict_LR_info['LERN'])+';LIRD='+str(dict_LR_info['LIRD'])+';LIRN='+str(dict_LR_info['LIRN'])+';LCRN='+str(dict_LR_info['LCRN'])
                                +';RERD='+str(dict_LR_info['RERD'])+';RERN='+str(dict_LR_info['RERN'])+';RIRD='+str(dict_LR_info['RIRD'])+';RIRN='+str(dict_LR_info['RIRN'])+';RCRN='+str(dict_LR_info['RCRN'])
                                +';ISN='+str(dict_LR_info['Inter_Num'])+';ISNMFS='+str(dict_LR_info['Inter_Most_Freq_Seq'])+';NSN='+str(dict_LR_info['NewSeq_Num'])
                                +';NSNMFS='+str(dict_LR_info['NewSeq_Most_Freq_Seq'])+';SLN='+str(dict_LR_info['Scarless_Num'])
                                +';CN='+str(dict_LR_info['Cigar_Num'])+';IVN='+str(dict_LR_info['Inversion_Num'])+';NSP='+str(dict_LR_info['New_Seq_Pattern'])
                                +';MNF='+str((dict_LR_info['MNF']))+';PFF='+str((dict_LR_info['PFF']))+';PLF='+str((dict_LR_info['PLF']))
                                ),
                                'GT:DP:AD:VF',
                                str(genotype)+':'+str(dict_LR_info['Read_Number']+dict_LR_info['Ref_Coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Ref_Coverage'])+':'+str(vaf)
                            ])
                            new = new + '\n'

                        else:
                            print()
                            print('Unknown rearrangement:', dict_LR_info['Type'])
                            print('dict_LR_info: ', dict_LR_info)
                            logging.error('Unknown rearrangement:'+str(dict_LR_info['Type']))
                            logging.error('dict_LR_info: '+str(dict_LR_info))
                            sys.exit()
                        
                        vcf_file.write(new)
                    
                    # similar calculations for insertions
                    else:
                        alt_cov_to_prob = (1.27*dict_LR_info['Read_Number'])/(2.512+dict_LR_info['Read_Number'])
                        vaf = dict_LR_info['Read_Number'] / dict_LR_info['Total_coverage']
                        vaf_to_prob = (1.218*vaf)/(0.095+vaf)
                        if dict_LR_info['MQ_median'] >= 60:
                            mq_to_prob = 1
                        else:
                            mq_median = dict_LR_info['MQ_median']
                            mq_to_prob =  0.000282*(mq_median**2)-0.003384
                        qual = alt_cov_to_prob * vaf_to_prob * mq_to_prob
                        qual = round(qual*100, 2) 
                        vaf = round(vaf, 2) 
                        if qual <= 0:
                            qual = 0.0
                        if vaf>=0.9:
                            genotype = '1/1'                            
                        else:
                            genotype = '0/1'
                        
                        if dict_LR_info['ID'] in self.connected_INS_dict.keys():
                            dict_LR_info['SBLR'] = self.connected_INS_dict[dict_LR_info['ID']]
                            dict_LR_info['SBLR'] = ['LR'+str(x) for x in list(dict_LR_info['SBLR'])]
                            # print("dict_LR_info['SBLR']", dict_LR_info['SBLR'])
                        else:
                            dict_LR_info['SBLR'] = ['""']
        

                        new = '\t'.join([
                            str(dict_LR_info['Chrom']),
                            str(dict_LR_info['Pos']),
                            'INS'+str(dict_LR_info['ID']),
                            'N',
                            str(lr_type_from_dict),
                            str(qual),
                            str(pass_no),
                            ('SVTYPE='+str(lr_type_from_dict)+';SVLEN='+str(dict_LR_info['LR_len'])+';SR='+str(dict_LR_info['Read_Number'])
                            +';MQ='+str(dict_LR_info['MQ_median'])+';DP='+str(dict_LR_info['Total_coverage'])+';VAF='+str(vaf)
                            +';MUTM='+str((dict_LR_info['Read_Muts_median']))+';MUTV='+str((dict_LR_info['Read_Muts_variance']))
                            +';SBLR='+str(','.join(dict_LR_info['SBLR']))
                            ),
                            'GT:DP:AD:VF',
                            str(genotype)+':'+str(dict_LR_info['Total_coverage'])+':'+str(dict_LR_info['Read_Number'])+','+str(dict_LR_info['Total_coverage']-dict_LR_info['Read_Number'])+':'+str(vaf)
                        ])
                        new = new + '\n'
                        vcf_file.write(new)
                    
                    # done_work += 1
                    # self.show_perc_work(done_work, all_work)
    
    # general function
    def main_func(self):
        # read the bam-file with all contigs --> save contigs names
        self.read_bam_file()
        self.sort_chrom_names()
        # create new header for vcf-file
        self.build_new_header() 
        # read the file with all found genomic rearrangements
        ff = pysam.FastaFile(self.reference_genome, filepath_index=self.reference_genome+'.fai') 
        # start_time = monotonic()
        self.read_junc_file()
        ff.close()
        self.read_file_INS()
        print('Sorting data for VCF-files...')
        logging.info('Sorting data for VCF-files...')
        # sort all rearrangement by chr + coordinate
        for key, val in self.general_dict_LR.items():
            self.general_dict_LR[key][1] = list(sorted(val[1], key=lambda x: x['Pos']))
        # sort second dictionary for pos 2
        for key, val in self.general_dict_LR_pos2.items():
            self.general_dict_LR_pos2[key][1] = list(sorted(val[1], key=lambda x: x['Pos2']))
        # write header to vcf-file
        with open(self.path_vcf_out, 'w', encoding='utf-8') as vcf_file:
            for line in self.new_header:
                vcf_file.write(line+'\n')
        # create two temporary files for repeats annotation
        if not os.path.isdir(self.temp_dir):
            os.mkdir(self.temp_dir)
        print()
        print('Creating temporary VCF-files for annotating genomic elements...')
        logging.info('Creating temporary VCF-files for annotating genomic elements...')
        self.build_temp_header()
        with open(self.temp_vcf_file_pos1, 'w', encoding='utf-8') as temp_1_vcf:
            for line in self.temp_header:
                temp_1_vcf.write(line+'\n')
        with open(self.temp_vcf_file_pos2, 'w', encoding='utf-8') as temp_2_vcf:
            for line in self.temp_header:
                temp_2_vcf.write(line+'\n')
        self.build_temp_vcf_data()
        # adjust bed file
        print('Preparing files for vcfanno...')
        logging.info('Preparing files for vcfanno...')
        self.path_bed_adjusted = self.path_bed[:-4] + '_adjusted.bed'
        self.adjust_bed_coordinates(self.path_bed, self.path_bed_adjusted)
        self.create_lua_file(self.path_lua)
        self.process_bed_file(self.path_bed_adjusted)
        self.create_toml_file(self.compressed_bed_path, self.path_toml)
        # delete data from dict with pos2
        self.general_dict_LR_pos2 = {}
        # annotate repeats near the ends of rearrangements
        print('Using vcfanno...')
        logging.info('Using vcfanno...')
        sp.check_output(self.path_vcfanno+" -lua "+self.path_lua+" "+self.path_toml+" "+self.temp_vcf_file_pos1+" > "+self.temp_vcf_ann_1, shell=True, stderr=sp.STDOUT).decode('utf-8')
        sp.check_output(self.path_vcfanno+" -lua "+self.path_lua+" "+self.path_toml+" "+self.temp_vcf_file_pos2+" > "+self.temp_vcf_ann_2, shell=True, stderr=sp.STDOUT).decode('utf-8')
        logging.info(self.path_vcfanno+" -lua "+self.path_lua+" "+self.path_toml+" "+self.temp_vcf_file_pos1+" > "+self.temp_vcf_ann_1)
        logging.info(self.path_vcfanno+" -lua "+self.path_lua+" "+self.path_toml+" "+self.temp_vcf_file_pos2+" > "+self.temp_vcf_ann_2)
        # read annotation results
        print('Saving the results of rearrangement analysis for overlapping genomic elements...')
        logging.info('Saving the results of rearrangement analysis for overlapping genomic elements...')
        self.read_results_vcfanno()
        if self.not_remove_trash == False:
            if os.path.isdir(self.temp_dir):
                rmtree(self.temp_dir)
                logging.info('rmtree(self.temp_dir) '+self.temp_dir)
            os.remove(self.path_toml)
            logging.info('os.remove(self.path_toml) '+self.path_toml)
            os.remove(self.path_lua)
            logging.info('os.remove(self.path_lua) '+self.path_lua)
        # define secondary boundries and connectes INS
        for key, val in self.types_secondary_boundry.items():
            types_secondary_boundry_list = list(zip(val[0], val[1], val[2], val[3], val[4], val[5], val[6]))
            types_secondary_boundry_list.sort(key=lambda x: x[0])
            self.types_secondary_boundry[key] = [[x[0] for x in types_secondary_boundry_list],
                                        [x[1] for x in types_secondary_boundry_list],
                                        [x[2] for x in types_secondary_boundry_list],
                                        [x[3] for x in types_secondary_boundry_list],
                                        [x[4] for x in types_secondary_boundry_list],
                                        [x[5] for x in types_secondary_boundry_list],
                                        [x[6] for x in types_secondary_boundry_list]]
        # print(self.ins_chrom_dict)
        for key, val in self.ins_chrom_dict.items():
            ins_chrom_list = list(zip(val[0], val[1], val[2]))
            ins_chrom_list.sort(key=lambda x: x[0])
            self.ins_chrom_dict[key] = [[x[0] for x in ins_chrom_list],
                                        [x[1] for x in ins_chrom_list],
                                        [x[2] for x in ins_chrom_list]]
        # print(self.ins_chrom_dict)
        self.define_boundries()
        # write data to vcf-file
        self.build_new_vcf_data()
        print()
        text_final = 'VCF-file was created successfully!'
        logging.info('VCF-file was created successfully!')
        print('-'*len(text_final))
        print(text_final)
        print('-'*len(text_final))

# if __name__ == '__main__':

#     parser = argparse.ArgumentParser()
#     parser.add_argument("-junc", "--junc_file", type=str, required=True,
#                         help="full path to the file with all found genomic rearrangements")
#     parser.add_argument("-ins", "--ins_file", type=str, required=True,
#                         help="full path to the file with all found insertions")
#     parser.add_argument("-bam", "--bam_file", type=str, required=True,
#                         help="full path to the bam-file")
#     parser.add_argument("-ref", "--reference", type=str, required=True,
#                         help="full path to the reference genome file")
#     parser.add_argument("-out_vcf", "--out_vcf_file", type=str, required=True,
#                         help="full path to the output vcf-file")
#     parser.add_argument("-vcfanno", "--vcf_anno", type=str, required=True,
#                         help="full path to the vcfanno binary file")
#     # parser.add_argument("-toml", "--toml_file", type=str, required=True,
#     #                     help="full path to the configuration toml-file")
#     # parser.add_argument("-lua", "--lua_file", type=str, required=True,
#     #                     help="full path to the custom lua-file")
#     parser.add_argument("-bed", "--bed_file", type=str, required=True,
#                         help="full path to the bed-file with repeats and other genome elements for program annotation")
#     parser.add_argument("-th", "--threads", type=int, required=False,
#                         help="number of processes for pool",
#                         default=cpu_count())
#     parser.add_argument("-nrt", "--not_remove_trash", action='store_true', required=False,
#                         help="use this parameter if you don't want to remove temporary files (.lua, .toml, directories with vcfanno results)",
#                         default=False)
    
#     args = parser.parse_args()

#     analyze_LR = AnalyzeLR(args)
#     analyze_LR.main_func()
