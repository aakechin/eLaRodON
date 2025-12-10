import pysam
import logging

class BamReader():

    def __init__(self, bamFile, ths, log_path):

        self.bamFile = bamFile
        self.threads = ths
        self.log_path = log_path
        logging.basicConfig(level=logging.DEBUG, filename=self.log_path, filemode="a",
                format="%(asctime)s %(funcName)s %(lineno)d %(message)s")

    def read_bam(self):

        bam = pysam.AlignmentFile(self.bamFile, 'rb', threads=self.threads)
        chromosomes = bam.references
        lengths = bam.lengths
        chrom_lengths = dict(zip(chromosomes, lengths))
        logging.info('BAM file was read: '+str(self.bamFile))
        logging.info('Chromosomes: '+str(chromosomes))
        logging.info('Lengths: '+str(lengths))
        
        return chrom_lengths

        
        

