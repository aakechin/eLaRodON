from read_bam_file import BamReader
from ONTLRcaller import ONTLRCaller
from joinONTLRs import JoinLR
from alignment_INS import InsertionProcessor
from define_type_create_vcf_LRs import AnalyzeLR
import argparse
from multiprocessing import cpu_count
import glob
import os
import sys


if __name__ == '__main__':

    par = argparse.ArgumentParser(description='This script calls large rearrangements in ONT data')

    # ONTLRcaller
    par.add_argument('--workdir', '-dir',
                    dest='workDir', type=str,
                    help='full path to working directory',
                    required=True)
    par.add_argument('--bam-file', '-bam',
                    dest='bamFile', type=str,
                    help='full path to BAM file with NGS reads',
                    required=True)
    par.add_argument('--divide-chroms', '-div',
                    dest='div_chroms', action='store_true', default=False,
                    help='Divide chromosome to analyze. Default: all chromosomes',
                    required=False)
    par.add_argument('--div-length', '-dvlen',
                    dest='len_division', type=int, default=0,
                    help='Length of regions for division chromosome to analyze',
                    required=False)
    par.add_argument('--minimal-length', '-len',
                    dest='minVarLen', type=int,
                    help='minimal acceptable length of variant (for InDels). Default: 50',
                    required=False, default=50)
    par.add_argument('--minimal-clipped-length', '-clip',
                    dest='minClipLen', type=int,
                    help='minimal acceptable length of a clipped read part. Default: 100',
                    required=False, default=100)
    par.add_argument('--dist-to-join-trl', '-dist',
                    dest='distToJoinTrl', type=int,
                    help='maximal acceptable distance between two coordinates. Default: 1000',
                    required=False, default=1000)
    par.add_argument('--threads', '-th',
                    dest='threads', type=int,
                    help='number of threads to call LRs. Default: 4',
                    required=False, default=4)
    
    # joinONTLRs
    par.add_argument('--input-files', '-in',
                    dest='inFiles', type=str,
                    help='regular expression for CSV-files with called LRs (if you have changed files names)',
                    required=False)
    par.add_argument('--maximal-distance-join', '-join',
                    dest='maxDistToJoin', type=int,
                    help='maximal acceptable distance of two neighboring fusions to join them. Default: 30',
                    required=False, default=30)
    par.add_argument('--output-lrs', '-lrs',
                    dest='outLrsFile', type=str,
                    help='CSV output file for other large rearrangements (if you want special names)',
                    required=False)
    par.add_argument('--output-ins', '-ins',
                    dest='outInsFile', type=str,
                    help='CSV output file for insertions (if you want special names)',
                    required=False)
    
    # alignment INS
    par.add_argument('--ref-genome', '-ref',
                    dest='refGen', help='path to the reference genome',
                    required=True)
    par.add_argument('--name', '-n', help='Name of INS alignment run', 
                    dest='nameINS', required=False)
    par.add_argument('--sam_file', '-sam',
                    dest='samFile', help='path to the SAM-file if mapping has already been done',
                    required=False)
    par.add_argument("--not_remove_trash_align", "-nrt_ins",
                    dest='notRemoveTrashAlign', action='store_true', required=False,
                    help="Use this parameter if you don't want to remove temporary files after INS alignment (.fasta)",
                    default=False)
    
    # define type
    par.add_argument("--junc_file", "-juncf", 
                    dest='juncFile', type=str, required=False,
                    help="full path to the file with all found genomic rearrangements after the alignment of insertions")
    par.add_argument("--ins_file", "-insf",
                    dest='insFile', type=str, required=False,
                    help="full path to the file with all found insertions after the alignment of insertions")
    par.add_argument("--out_vcf", "-out",
                    dest='outVCF', type=str, required=False,
                    help="full path to the output vcf-file")
    par.add_argument("--vcf_anno", "-vcfanno",
                    dest='VCFanno', type=str, required=True,
                    help="full path to the vcfanno binary file")
    par.add_argument("--bed_file", "-bed",
                    dest='bedFile', type=str, required=True,
                    help="full path to the bed-file with repeats and other genome elements for program annotation")
    par.add_argument("--not_remove_trash_anno", "-nrt_anno",
                    dest='notRemoveTrashAnno', action='store_true', required=False,
                    help="use this parameter if you don't want to remove temporary files (.lua, .toml, directories with vcfanno results)",
                    default=False)
    
    args=par.parse_args()

    def check_path(path, arg_name):
        if not os.path.exists(path):
            raise FileNotFoundError(f"ERROR: {arg_name} does not exist: {path}")
        return path

    try:
        check_path(args.bamFile, "--bam-file")
        check_path(args.refGen, "--ref-genome")
        check_path(args.VCFanno, "--vcf_anno")
        check_path(args.bedFile, "--bed_file")

        if not os.path.isdir(args.workDir):
            os.makedirs(args.workDir, exist_ok=True)
            os.makedirs(os.path.join(args.workDir, 'supplemetary'), exist_ok=True)
            print(f"NOTE: Created working directory at {args.workDir}")
        
        if args.div_chroms and args.len_division <= 0:
            print('ERROR: Set the fragment size value for chromosome division!')
            print('Use --div-length or -dvlen with positive integer value')
            sys.exit(1)
            
    except FileNotFoundError as e:
        print(str(e))
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: {str(e)}")
        sys.exit(1)

    bam_file_name = args.bamFile.split('/')[-1][:-4]
    bam_results_path = args.workDir +'/supplemetary/' + bam_file_name 

    if not hasattr(args, 'inFiles') or args.inFiles is None:
        args.inFiles = bam_results_path+'.junction_stat.*.*ions.csv'
    if not hasattr(args, 'outLrsFile') or args.outLrsFile is None:
        args.outLrsFile = bam_results_path+'.junction_stat.LRs_join100.csv'
    if not hasattr(args, 'outInsFile') or args.outInsFile is None:
        args.outInsFile = bam_results_path+'.junction_stat.INS_join100.csv'
    if not hasattr(args, 'nameINS') or args.nameINS is None:
        args.nameINS = bam_file_name
    if not hasattr(args, 'juncFile') or args.juncFile is None:
        args.juncFile = bam_results_path +'.junction_stat.LRs_join100_mapped_ins.csv'
    if not hasattr(args, 'insFile') or args.insFile is None:
        args.insFile = bam_results_path +'.junction_stat.INS_join100_mapped_ins.csv'
    if not hasattr(args, 'outVCF') or args.outVCF is None:
        args.outVCF = args.workDir + '/' + bam_file_name + '_all_LGRS.vcf'

    # ONTLRcaller
    print()
    print('=== ONTLRcaller ===')
    
    bamreader = BamReader(args.bamFile, args.threads)
    chromosomes_dict = bamreader.read_bam()
    # we analyze only main chomosomes for human genome
    chom_list = set(chromosomes_dict.keys())
    chrom_intersection_v1 = chom_list & {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'}
    chrom_intersection_v2 = chom_list & {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M'}
    if len(chrom_intersection_v1) == 25:
        chromosomes_dict = {k: v for k, v in chromosomes_dict.items() if k in chrom_intersection_v1}
    elif len(chrom_intersection_v2) == 25:
        chromosomes_dict = {k: v for k, v in chromosomes_dict.items() if k in chrom_intersection_v2}
    for chrom, chrom_len in chromosomes_dict.items():
        if args.div_chroms == False:
            ontc=ONTLRCaller(args.bamFile,
                            chrom,
                            args.start,
                            args.end,
                            args.threads,
                            args.minVarLen,
                            args.minClipLen,
                            args.distToJoinTrl)
            ontc.readBamFile()
        else:
            n = args.len_division
            for start in range(1, chrom_len, n):
                end = min(start + n, chrom_len)
                ontc=ONTLRCaller(args.bamFile,
                                chrom,
                                start,
                                end,
                                args.threads,
                                args.minVarLen,
                                args.minClipLen,
                                args.distToJoinTrl)
                ontc.readBamFile()

    #joinONTLRs
    print()
    print('=== joinONTLRs ===')
    ds = glob.glob(args.inFiles)
    if len(ds)==0:
        print('ERROR (1)! No files were chosen:')
        print(args.inFiles)
        exit(1)
    print('The number of files chosen:', len(ds))
    print('Reading input files with', args.threads, 'threads...')
    jl = JoinLR(ds, args)
    print('The total number of fusions:', sum(len(item) for item in jl.allFusions.values()))
    print('The total number of insertions:', sum(len(item) for item in jl.allInsertions.values()))
    jl.joinAllSimilarLRs()
    jl.writeFusionToOutput(args.outLrsFile)
    jl.writeInsertionsToOutput(args.outInsFile)
    print()

    # alignment INS
    print()
    print('=== alignment_INS ===')
    alignINS_processor = InsertionProcessor(args.outInsFile, 
                                            args.outLrsFile, 
                                            args.refGen, 
                                            args.nameINS, 
                                            args.samFile,
                                            args.threads, 
                                            args.notRemoveTrashAlign)
    alignINS_processor.read_insertions()
    fasta_file = args.workDir +'/supplemetary/'+ args.nameINS + '_insertions.fasta'
    alignINS_processor.write_fasta(fasta_file)
    alignINS_processor.map_with_minimap2(fasta_file)
    alignINS_processor.analyze_mapping_results()
    alignINS_processor.merge_tandem_duplications()
    alignINS_processor.copy_file(args.outLrsFile, args.outLrsFile[:-4]+'_mapped_ins.csv')
    alignINS_processor.extract_sequences_before_after()
    alignINS_processor.update_large_rearrangements_file()
    updated_csv = args.outInsFile[:-4]+'_mapped_ins.csv'
    alignINS_processor.update_original_csv(updated_csv)
    print('Analysis of INS was done!')

    # define types of LGRs
    print()
    print('=== define_type_create_vcf_LRs ===')
    analyze_LR = AnalyzeLR(args=args)
    analyze_LR.main_func()

