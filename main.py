from ONTLRcaller import ONTLRCaller
from joinONTLRs import JoinLR
from alignment_INS import InsertionProcessor
from define_type_create_vcf_LRs import AnalyzeLR
import argparse
from multiprocessing import cpu_count
import glob
import os
# import ONTLRcaller_old as ontlr


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
    par.add_argument('--chromosome', '-chr',
                    dest='chrom', type=str,
                    help='chromosome to analyze. Default: all chromosomes',
                    required=False)
    par.add_argument('--start', '-start',
                    dest='start', type=int,
                    help='start coordinate for extracting reads. Default: whole chromosome or whole genome',
                    required=False, default=None)
    par.add_argument('--end', '-end',
                    dest='end', type=int,
                    help='end coordinate for extracting reads. Default: whole chromosome or whole genome',
                    required=False, default=None)
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

    bam_file_name = args.bamFile.split('/')[-1][:-4]
    if not os.path.isdir(args.workDir):
        os.mkdir(args.workDir)
    bam_results_path = args.workDir + '/' + bam_file_name

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
        args.outVCF = bam_results_path+'_all_LGRS.vcf'

    # print(args)

    # ONTLRcaller
    print()
    print('=== ONTLRcaller ===')
    ontc=ONTLRCaller(args.bamFile,
                     args.chrom,
                     args.start,
                     args.end,
                     args.threads,
                     args.minVarLen,
                     args.minClipLen,
                     args.distToJoinTrl)
    ontc.readBamFile()
    # without class
    # ontlr.readBamFile(args.bamFile,
    #                  args.chrom,
    #                  args.start,
    #                  args.end,
    #                  args.threads,
    #                  args.minVarLen,
    #                  args.minClipLen,
    #                  args.distToJoinTrl)
    print()
    print('ONTLRcaller finished')

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
    fasta_file = args.nameINS + '_insertions.fasta'
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

