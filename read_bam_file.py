import pysam

class BamReader():

    def __init__(self, bamFile, ths):
        self.bamFile = bamFile
        self.threads = ths

    def read_bam(self):

        bam = pysam.AlignmentFile(self.bamFile, 'rb', threads=self.threads)
        chromosomes = bam.references
        lengths = bam.lengths
        chrom_lengths = dict(zip(chromosomes, lengths))
        
        return chrom_lengths

        
        

