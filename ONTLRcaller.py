# This script calls LRs in BAM file of ONT data
# 
import re
import sys
import argparse
import pysam
import math
from multiprocessing import Pool,Queue,Manager,Process,set_start_method
from multiprocessing.pool import ThreadPool
from threading import Thread
import statistics as stat
from tqdm import tqdm
import logging

class ONTLRCaller():

    def __init__(self,bamFile,log_path,
                 chrom=None,start=None,end=None,
                 threads=8,minVarLen=50,
                 minClipLen=100,distToJoinTrl=1000, workdir=''):
        if ((start==None or
             end==None) and 
            start!=end):
            print('ERROR (2)! If the start or the end for extracting reads is specified, another argument should also be defined.')
            print(start,end)
            exit(2)
        if bamFile[-2:] != 'gz':
            self.bam_file_name = bamFile.split('/')[-1][:-4]
        else:
            self.bam_file_name = bamFile.split('/')[-1][:-7]
        self.bamFile=bamFile
        self.log_path = log_path
        logging.basicConfig(level=logging.DEBUG, filename=self.log_path, filemode="a",
                            format="%(asctime)s %(funcName)s %(lineno)d %(message)s")
        self.workdir = workdir
        #self.outFile=outFile
        self.chrom=chrom
        self.start=start
        self.end=end
        self.threads=threads
        self.minVarLen=minVarLen
        self.minClipLen=minClipLen
        self.distToJoinTrl=distToJoinTrl
        # Stores coverage of genome regions splitted by minClipLen (e.g. 30 bp)
        # {chr1:{pos1:[cov1_primary,cov1_second],
        #        pos2:[cov1_primary,cov1_second]...}}
        self.genomeCoverage={}
        self.readToSecs={}
        #self.q=Queue()

    def readBamFile(self):
        global outFiles
        print('Reading BAM-file...')
        chromosomes=[self.chrom]       
        for chrom in chromosomes:
            if chrom not in self.genomeCoverage.keys():
                self.genomeCoverage[chrom]=[]
            # Collect all reads of one chromosome
            # {read_id:[minVarLen,minClipLen,read,sup_read1,sup_read2... ]}
            reads=[]
            # Stores secondary reads
            # Stores stat of read part intersections
            # [[start<end],[start==end],[start>end]]
            # For each variant, we store intersection length and its sequence
            # [intersection/distance length,sequence]
            # For start==end, we store +-1 bases
            bam=pysam.AlignmentFile(self.bamFile,threads=self.threads)
            if (self.start!=None and
                self.end!=None):
                regionForName='_'.join([str(chrom),
                                        str(self.start),
                                        str(self.end)])
            else:
                regionForName=str(chrom)
            for read in tqdm(bam.fetch(chrom,
                                  self.start,
                                  self.end)):
                # Calculate positions rounded to 30
                readStart=self.minVarLen*math.ceil(read.pos/self.minVarLen)
                readEnd=self.minVarLen*math.floor((read.pos+read.alen)/self.minVarLen)
                covValueNum=0
                # If read is secondary
                if read.is_secondary:
                    if read.query_name not in self.readToSecs.keys():
                        self.readToSecs[read.query_name]=[]
                    self.readToSecs[read.query_name].append([str(read.reference_name),
                                                        int(read.pos),
                                                        int(read.alen),
                                                        read.cigartuples,
                                                        read.is_reverse])
                    continue
                for pos in range(readStart,readEnd+self.minVarLen,self.minVarLen):
                    simplePos=int(math.floor(pos/self.minVarLen))
                    # If the position to add is more than the length of the list
                    if simplePos>=len(self.genomeCoverage[chrom]):
                        # We add all skipped elements
                        for i in range(len(self.genomeCoverage[chrom])-1,
                                       simplePos+1):
                            self.genomeCoverage[chrom].append(0)#[0,0]
                    self.genomeCoverage[chrom][simplePos]+=1 # [covValueNum]
                if (read.is_supplementary or
                    read.is_secondary):
                    continue
                # We take only primary mappings because they contain whole sequence
                # We save only necessary information about the read
                # [cigartuples,read_length,position,is_reverse,SA_tag]
                if read.has_tag('NM'):
                    nmTag=read.get_tag('NM')
                else:
                    nmTag=0
                reads.append([self.minVarLen,
                              self.minClipLen,
                              [],
                              [read.cigartuples,
                               read.infer_read_length(),
                               read.pos,
                               read.is_reverse,
                               read.query_sequence,
                               chrom,
                               read.alen,
                               read.query_name,
                               [read.mapping_quality],
                               nmTag]])
                if read.has_tag('SA'):
                    reads[-1][3].append(str(read.get_tag('SA')))
                else:
                    reads[-1][3].append([])
            if len(reads)==0:
                continue
            print('Analyzing reads for region '+regionForName+'...')
            logging.info('Analyzing reads for region '+regionForName+'...')
            # Create pool of processes
            p=ThreadPool(self.threads)
            # Write to output
            outFileName='.'.join([self.bam_file_name,
                                'junction_stat',regionForName,'xls'])
            outFileName = self.workdir + outFileName
            outFiles=[]
            outFileName1=outFileName[:outFileName.rfind('.')]+'.fusions.csv'
            outFiles.append(open(outFileName1,'w'))
            outFiles[-1].write('\t'.join(['Chrom1','Pos1','Strand1','Junction_Side1',
                                            'Chrom2','Pos2','Strand2','Junction_Side2',
                                            'MQ','Ref_Cov','Read_Len1','Read_Len2','Read_Muts',
                                            'Joint_Type','Joint_Length','Joint_Sequence',
                                            'First_Last_Part','Num_of_parts','Seq_before_BND','Seq_after_BND',
                                            'Read_Name'])+'\n')
            outFileName2=outFileName[:outFileName.rfind('.')]+'.insertions.csv'
            outFiles.append(open(outFileName2,'w'))
            outFiles[-1].write('\t'.join(['Chrom1','Pos1','Insertion_Length',
                                            'Sequence','Read_Number','Total_Cov',
                                            'MQ','Read_Muts1','Read_Name'])+'\n')
            outFileName3=outFileName[:outFileName.rfind('.')]+'.new_sequences.fa'
            outFileName4=outFileName[:outFileName.rfind('.')]+'.junctions.fa'
            outFiles.append(open(outFileName3,'w'))
            outFiles.append(open(outFileName4,'w'))
            for i,read in enumerate(reads):
                if read[3][7] in self.readToSecs.keys():
                    reads[i][2]=self.readToSecs[read[3][7]]
            results=[]
            m=Manager()
            self.thsize=0
            results=[]
            for res in tqdm(p.imap_unordered(self.callLRs,reads),
                            total=len(reads)):
                results.append(res)
            print('Writing fusions...')
            logging.info('Writing fusions...')
            for result in tqdm(results,
                                total=len(results)):
                t2=Thread(target=self.writeOutput,
                            args=(result,),daemon=True)
                t2.start()
                t2.join()
            for outFile in outFiles:
                outFile.close()

    def showPercWork(self,done,allWork):
        percDoneWork=round((done/allWork)*100,2)
        sys.stdout.write("\r"+str(percDoneWork)+"%")
        sys.stdout.flush()

    def makeAllRevCompl(self,seqsList):
        seqsListRC=[]
        for seq in seqsList:
            try:
                int(seq)
                seqsListRC.append(seq)
            except ValueError:
                seqsListRC.append(self.revComplement(seq))
        return(seqsListRC)

    def revComplement(self,seq):
        nucToRC={'A':'T','T':'A','G':'C','C':'G',
                 'a':'t','t':'a','g':'c','c':'g'}
        seqRC=[]
        for n in seq:
            seqRC.append(nucToRC[n])
        return(''.join(seqRC))        

    def getInDelsFromRead(self,cigar,chrom,curRefPos,
                          readSeq,mapQual,NM,
                          name,minVarLen):
        deletions=[]
        insertions=[]
        ciTypeStrToInt={'M':0,'I':1,'D':2,
                        'N':3,'S':4,'H':5,
                        'P':6,'X':7,'=':9}                    
        # If cigar is provided as string, convert it to tuple
        if type(cigar)==str:
            cigarTuples=[]
            cigarPat=re.compile(r'(\d+)([MIDNSHPX=]+)')
            for ciLen,ciType in cigarPat.findall(cigar):
                cigarTuples.append((ciTypeStrToInt[ciType],
                                    int(ciLen)))
        else:
            cigarTuples=cigar
        # Extract DELs and INSs from CIGAR
        curReadPos=0
        alignedLenPart1=0
        for ciType,ciLen in cigarTuples:
            if (ciType in [1,2] and
                ciLen>=minVarLen):
                # If it is insertion
                if ciType==1:
                    insertedSeq=readSeq[curReadPos:curReadPos+ciLen]
                    if (len(insertions)>0 and
                        abs(insertions[-1][1]-curRefPos)<200):
                        insertions[-1][2]+=ciLen
                        insertions[-1][5]+=insertedSeq
                    else:
                        insertions.append([chrom,
                                           curRefPos,
                                           ciLen,
                                           mapQual,NM,
                                           insertedSeq,
                                           name])
                elif ciType==2:
                    deletions.append([chrom,
                                      curRefPos,
                                      curRefPos+ciLen-1,
                                      mapQual,NM,
                                      curReadPos,
                                      alignedLenPart1])
            if ciType in [0,2,3]:
                curRefPos+=ciLen
                alignedLenPart1+=ciLen
            if ciType in [0,1,4]:
                curReadPos+=ciLen
        return(deletions,
               insertions)

    # The following function calls large rearrangements from read object of pysam
    def callLRs(self,readInfo):
        self.thsize+=1
        minVarLen=readInfo[0]
        minClipLen=readInfo[1]
        # Contains supplementary alignment positions in the read
        # [[16,181],[180,1152]...]
        readSupPoses=[]
        # We take the first aligment met and check for some fusions
        # [cigartuples,read_length,position,
        #  is_reverse,seq,chrom,alen,name,[MQ],NM,SA_tag]
        mainReadPart=readInfo[3]
        secAligns=readInfo[2]
        # Get main read part strand
        if mainReadPart[3]:
            mainReadPartStrand='-'
        else:
            mainReadPartStrand='+'
        # Deletions are stored as [start,end]
        deletions=[]
        # Insertions are stored as [position,length,seq]
        insertions=[]
        # Stores CIGAR regular expression
        cigarPat=re.compile(r'(\d+)([MIDNSHPX=]+)')
        if mainReadPart[0]==None:
            return([],mainReadPart,
                   deletions,insertions)
        # Extract DELs and INSs from CIGAR
        inDels=self.getInDelsFromRead(mainReadPart[0],
                                 mainReadPart[5],
                                 mainReadPart[2]+1,
                                 mainReadPart[4],
                                 mainReadPart[8],
                                 mainReadPart[9], # NM
                                 mainReadPart[7],
                                 minVarLen)
        deletions.extend(inDels[0])
        insertions.extend(inDels[1])
        # If the left or the right part of the cigar contains clipped bases
        if (((mainReadPart[0][0][0] in [4,5] and
              # and it is at least of minimal accepted length
              mainReadPart[0][0][1]>=minClipLen) or
             (mainReadPart[0][-1][0] in [4,5] and
              # and it is at least of minimal accepted length
              mainReadPart[0][-1][1]>=minClipLen)) and
            # One more IMPORTANT condition is the presence of supplementary alignments for the read
            # IMPORTANT! We do not consider bases clipped at read ends as LRs
            len(mainReadPart[10])>0):
            # Stores left or right side by which this part is connected to preavious and next parts, respectively
            sidesForConn=['N','N']
            # Determine coordinates of the main alignment in the read
            # Start is a number of clipped bases +1
            if mainReadPart[0][0][0] in [4,5]:
                start=mainReadPart[0][0][1]+1
                # if it is on '-' strand
                if mainReadPart[3]:
                    sidesForConn[0]='L'
                else:
                    sidesForConn[0]='L'
            else:
                start=1
                
            # End is a start+aligned_length-1
            if mainReadPart[0][-1][0] in [4,5]:
                end=mainReadPart[1]-mainReadPart[0][-1][1]
                # if it is on '-' strand
                if mainReadPart[3]:
                    sidesForConn[1]='R'
                else:
                    sidesForConn[1]='R'
            else:
                end=mainReadPart[1]
            # Add this part to read parts
            # startInReadSeq,endInReadSeq,chrom,pos,alen,NM,strand
            readSupPoses.append([start,end,
                                 mainReadPart[5],
                                 mainReadPart[2]+1,
                                 mainReadPart[6],
                                 int(mainReadPart[9]),
                                 mainReadPartStrand,
                                 sidesForConn])
            # Then go through all supplementary alignments and determine their coordinates in the read
            for readSup in mainReadPart[10].split(';')[:-1]:
                # Get supplementary info
                supChrom,supCoord,supStrand,supCigar,supMQ,supNM=readSup.split(',')
                # We check if some secondary mapping with better alignment has been saved
                mainReadPart[8].append(int(supMQ))
                # Stores left or right side by which this part is connected to previous and next parts, respectively
                sidesForConn=['N','N']
                # If strand of supplementary is the same as in the main part
                if mainReadPartStrand==supStrand:
                    # Then its CIGAR coordinates correspond main read part
                    cigarMatches=cigarPat.findall(supCigar)
                    # We determine this supplement alignment position in the read
                    # If its 1st CIGAR part is S
                    if cigarMatches[0][1]=='S':
                        # Then the start is the number of clipped bases +1
                        start=int(cigarMatches[0][0])+1
                        # If supplementary is on '+' strand
                        if supStrand=='+':
                            if mainReadPartStrand=='+':
                                sidesForConn[0]='L'
                            else:
                                sidesForConn[1]='L'
                        else:
                            if mainReadPartStrand=='+':
                                sidesForConn[1]='L'
                            else:
                                sidesForConn[0]='L'
                    else:
                        # Then the start is 1
                        start=1
                    # If the last CIGAR part is S
                    if cigarMatches[-1][1]=='S':
                        # Then the end is a (read_length - the_number_of_clipped_bases) 
                        end=mainReadPart[1]-int(cigarMatches[-1][0])
                        # If supplementary is on '+' strand
                        if supStrand=='+':
                            if mainReadPartStrand=='+':
                                sidesForConn[1]='R'
                            else:
                                sidesForConn[0]='R'
                        else:
                            if mainReadPartStrand=='+':
                                sidesForConn[0]='R'
                            else:
                                sidesForConn[1]='R'
                    else:
                        # Then end is a read length
                        end=mainReadPart[1]
                    # Calculate aligned part length
                    alen=0
                    for cigarPart in cigarMatches:
                        if cigarPart[1] in ['M','D']:
                            alen+=int(cigarPart[0])
                    if alen>=minClipLen:
                        for secReadPart in secAligns:
                            distToMain=abs(secReadPart[1]-mainReadPart[2])
                            distToSup=abs(secReadPart[1]-int(supCoord))
                            mainToSupDist=abs(int(supCoord)-mainReadPart[2])
                            # Compare values of the length of aligned parts for SA
                            if secReadPart[2]/alen>=0.8:
                                # Compare chromosomes
                                if (secReadPart[0]==mainReadPart[5] and
                                    (supChrom!=mainReadPart[5] or
                                     (distToMain+1000<min(10000,mainToSupDist) and
                                      distToSup>distToMain))):
                                    return([],mainReadPart,
                                           deletions,insertions)
                                    # If chromosome is the same like in the primary alignment and it's closer to it than current sup_align
                                    # Then, now we need to check that current secondary alignment refers to the same part of read like sup_align
                            # Compare values of the length of aligned parts for main mapping
                            if secReadPart[2]/mainReadPart[6]>=0.8:
                                # Compare chromosomes
                                if (secReadPart[0]==supChrom and
                                    (supChrom!=mainReadPart[5] or
                                     distToSup+1000<min(10000,mainToSupDist))):
                                    return([],mainReadPart,
                                           deletions,insertions)
                                    # If chromosome is the same like in the supplementary alignment and it's closer to it than current primary_al
                                    # Then, now we need to check that current secondary alignment refers to the same part of read like sup_align
                        # Add this read part to read parts
                        # startInReadSeq,endInReadSeq,chrom,pos,alen,strand
                        readSupPoses.append([start,end,
                                             supChrom,int(supCoord),
                                             alen,int(supNM),supStrand,
                                             sidesForConn])
                # If strand of supplementary is the opposite to the main part strand
                else:
                    # Then its CIGAR coordinates correspond main read part but reversed
                    cigarMatches=cigarPat.findall(supCigar)
                    # We determine this supplement alignment position in the read
                    # If its last CIGAR part is S
                    if cigarMatches[-1][1]=='S':
                        # Then the start is the number of clipped bases +1
                        start=int(cigarMatches[-1][0])+1
                        # If supplementary is on '+' strand
                        if supStrand=='+':
                            if mainReadPartStrand=='+':
                                sidesForConn[1]='R'
                            else:
                                sidesForConn[0]='R'
                        else:
                            if mainReadPartStrand=='+':
                                sidesForConn[0]='R'
                            else:
                                sidesForConn[1]='R'
                    else:
                        # Then the start is 1
                        start=1
                    # If the 1st CIGAR part is S
                    if cigarMatches[0][1]=='S':
                        # Then the end is a (read_length - the_number_of_clipped_bases) 
                        end=mainReadPart[1]-int(cigarMatches[0][0])
                        # If supplementary is on '+' strand
                        if supStrand=='+':
                            if mainReadPartStrand=='+':
                                sidesForConn[0]='L'
                            else:
                                sidesForConn[1]='L'
                        else:
                            if mainReadPartStrand=='+':
                                sidesForConn[1]='L'
                            else:
                                sidesForConn[0]='L'
                    else:
                        # Then end is a read length
                        end=mainReadPart[1]
                    # Calculate aligned part length
                    alen=0
                    for cigarPart in cigarMatches:
                        if cigarPart[1] in ['M','D']:
                            alen+=int(cigarPart[0])
                    if alen>=minClipLen:
                        for secReadPart in secAligns:
                            distToMain=abs(secReadPart[1]-mainReadPart[2])
                            distToSup=abs(secReadPart[1]-int(supCoord))
                            mainToSupDist=abs(int(supCoord)-mainReadPart[2])
                            # Compare values of the length of aligned parts for SA
                            if secReadPart[2]/alen>=0.8:
                                # Compare chromosomes
                                if (secReadPart[0]==mainReadPart[5] and
                                    (supChrom!=mainReadPart[5] or
                                     (distToMain+1000<min(10000,mainToSupDist) and
                                      distToSup>distToMain))):
                                    return([],mainReadPart,
                                           deletions,insertions)
                                    # If chromosome is the same like in the primary alignment and it's closer to it than current sup_align
                                    # Then, now we need to check that current secondary alignment refers to the same part of read like sup_align
                            # Compare values of the length of aligned parts for main mapping
                            if secReadPart[2]/mainReadPart[6]>=0.8:
                                # Compare chromosomes
                                if (secReadPart[0]==supChrom and
                                    (supChrom!=mainReadPart[5] or
                                     distToSup+1000<min(10000,mainToSupDist))):
                                    return([],mainReadPart,
                                           deletions,insertions)
                                    # If chromosome is the same like in the primary alignment and it's closer to it than current sup_align
                                    # Then, now we need to check that current secondary alignment refers to the same part of read like sup_align
                        # Add this read part to read parts
                        # startInReadSeq,endInReadSeq,chrom,pos,alen,strand
                        readSupPoses.append([start,end,
                                             supChrom,int(supCoord),
                                             alen,int(supNM),supStrand,
                                             sidesForConn])
        return(readSupPoses,mainReadPart,
               deletions,insertions)

    def countHomoNucNum(self,seq,nuc='T',
                        minLen=4):
        p=re.compile(nuc+'+')
        n=0
        sumLen=0
        for m in p.findall(seq):
            if len(m)>=minLen:
                n+=1
                sumLen+=len(m)
        return(n,sumLen)

    def countDifNucHomoTracts(self,seq):
        nucToHomoSumLen={'A':0,'T':0,'G':0,'C':0}
        nucToHomoTractNum={'A':0,'T':0,'G':0,'C':0}
        for nuc in nucToHomoSumLen.keys():
            n,sumLen=self.countHomoNucNum(seq,nuc=nuc)
            nucToHomoSumLen[nuc]=sumLen*100/len(seq)
            nucToHomoTractNum[nuc]=n #*100/len(seq)
        return(max(nucToHomoSumLen.values()),sum(nucToHomoTractNum.values()))

    def writeOutput(self,out):
        global outFiles
        self.thsize+=1
        fusions=[]
        try:
            readParts,(cigarTuples,readLen,readPos,readReverse,readSeq,readChrom,readAlignedLen,readName,mapQual,mapNM,readSups),deletions,insertions=out#out.get()
        except ValueError:
            print('ERROR',out)
            exit()
        # Process saved read parts and their sequences
        # First, sort read parts
        readParts.sort()
        fusion=None
        # Stores number of BND saved
        bndSaved=0
        homoTractSumValue=0
        # Go through all read parts and compare positions of the neighboring parts
        for i,readPart in enumerate(readParts[1:]):
            fusion=[]
            # Add info from readParts
            # startInReadSeq,endInReadSeq,chrom,pos,alen,strand
            # Add chrom
            # 'Chrom1'
            fusion.append(readParts[i][2])
            # 'Pos1'
            # For the previous part, we take end genome position regardless of its strand
            if readParts[i][7][1]=='R':
                fusion.append(readParts[i][3]+readParts[i][4]-1)
            elif readParts[i][7][1]=='L':
                fusion.append(readParts[i][3])
            else:
                break
            # 'Strand1','Junction_Side1'
            fusion.append(readParts[i][6])
            fusion.append(readParts[i][7][1])
            # 'Chrom2'
            fusion.append(readPart[2])
            # 'Pos2'
            if readPart[7][0]=='R':
                fusion.append(readPart[3]+readPart[4]-1)
            elif readPart[7][0]=='L':
                fusion.append(readPart[3])
            else:
                break
            # 'Strand2','Junction_Side2'
            fusion.append(readPart[6])
            fusion.append(readPart[7][0])
            # MQ
            fusion.append(stat.median(mapQual))
            # 'Read_Len1','Read_Len2'
            fusion.append(readParts[i][4])
            fusion.append(readPart[4])
            # NM
            errPerc=[str(int(round(int(readParts[i][5])*100/int(readParts[i][4]),0))),
                        str(int(round(int(readPart[5])*100/int(readPart[4]),0)))]
            fusion.append(','.join(errPerc))
            # Get junction sequence
            junctionSeq=readSeq[readParts[i][0]:readPart[1]]
            junctionPos1=readParts[i][1]-readParts[i][0]
            junctionPos2=readPart[0]-readParts[i][0]
            homoTractSumValue=0
            # Choose variant how we will write junction seq
            # If there is intersection or scarless junction
            #         -------<----->-----    ---||---
            # Or sequence between them is short (less 50 bp)
            #    --------->~~~~~~~~~~~~~<-----------------
            if (readPart[0]<=readParts[i][1]+1 or
                readPart[0]-readParts[i][1]<8):
                outFiles[3].write('_'.join(['>'+readName,
                                            str(fusion[0]),
                                            str(fusion[1]),
                                            str(fusion[4]),
                                            str(fusion[5]),
                                            'not_long_dist',
                                            str(junctionPos1),
                                            str(junctionPos2)])+'\n')
                outFiles[3].write(junctionSeq+'\n')
            else:
                # Else we write junction of left part with distseq
                outFiles[3].write('_'.join(['>'+readName,
                                            'long_dist1',
                                            str(fusion[0]),
                                            str(fusion[1]),
                                            str(fusion[4]),
                                            str(fusion[5]),
                                            str(junctionPos1),
                                            str(junctionPos2)])+'\n')
                outFiles[3].write(junctionSeq+'\n')
            # Then we write whole left part and whole right part
            # If start of the next less than end of the previous (intersection)
            if readPart[0]<=readParts[i][1]:
                interLen=readParts[i][1]-readPart[0]+1
                seq=readSeq[readPart[0]-1:readParts[i][1]+1-1]
                fusion.extend(['inter',
                                interLen,
                                seq])
                
            # If start==end+1
            elif readPart[0]-readParts[i][1]==1:
                interLen=1
                seq=readSeq[readPart[0]-1-1:readPart[0]+1-1]
                fusion.extend(['scarless',
                                interLen,
                                seq])
            # If start of the next more than end of the previous (distance)
            elif readPart[0]>readParts[i][1]:
                interLen=readPart[0]-readParts[i][1]-1
                seq=readSeq[readParts[i][1]+1-1:readPart[0]-1]
                homoTractSumLen,homoTractSumValue=self.countDifNucHomoTracts(seq)
                outFiles[2].write('_'.join(['>'+readName,
                                            str(fusion[0]),
                                            str(fusion[1]),
                                            str(fusion[4]),
                                            str(fusion[5]),
                                            str(junctionPos1),
                                            str(junctionPos2),
                                            str(readParts[i][1]+1),
                                            str(readPart[0]-1),
                                            str(homoTractSumValue)])+'\n')
                outFiles[2].write(seq+'\n')
                if interLen>6:
                    seq=''
                fusion.extend(['newseq',
                                interLen,
                                seq])
            # Add seqBeforeBND and seqAfterBND
            seqBefore=readSeq[max(0,readParts[i][1]-30):readParts[i][1]]
            seqAfter=readSeq[readPart[0]-1:readPart[0]-1+30]
            # Add First_Last_Part, Number_of_parts, seq_before_BND, seq_after_BND
            if i==0:
                firstLastPart=1
            elif i==len(readParts)-2:
                firstLastPart=2
            else:
                firstLastPart=0
            fusions.append(fusion+[homoTractSumValue,
                                    len(readParts),
                                    seqBefore,
                                    seqAfter,
                                    readName])
        # If number of saved BNDs is more than 2
        # then we need to try determine the 2nd coordinate of the LGR
        if len(fusions)>=2:
            fusionsJoint=[]
            for i in range(1,len(fusions)):
                # Go through all fusion boundaries and compare their
                # strands (+/-), positions (>/<) and side of junction (L/R)
                pos1,pos2,pos3,pos4=[fusions[i-1][1],
                                        fusions[i-1][5],
                                        fusions[i][1],
                                        fusions[i][5]]
                # Joint strand+size lists
                ss1,ss2,ss3,ss4=[fusions[i-1][2]+fusions[i-1][3],
                                    fusions[i-1][6]+fusions[i-1][7],
                                    fusions[i][2]+fusions[i][3],
                                    fusions[i][6]+fusions[i][7]]
                # Determine end of the 1st fragment
                if ss1[1]=='L':
                    fragEnd1=pos1+fusions[i-1][9]
                elif ss1[1]=='R':
                    fragEnd1=pos1-fusions[i-1][9]
                # Determine end of the 1st fragment
                if ss2[1]=='L':
                    fragEnd2=pos4+fusions[i-1][10]
                elif ss2[1]=='R':
                    fragEnd2=pos4-fusions[i-1][10]
                ssOrder=[ss1,ss2,ss3,ss4]
                trlSides=[['+R','+L','+R','+L'],
                            ['+R','-R','-L','+L'],
                            ['-L','+L','+R','-R'],
                            ['-L','-R','-L','-R']]
                tdSides=[['+R','+L','+R','+L'],
                            ['-L','-R','-L','-R']]
                invSides=[['-L','+L','+R','-R'],
                            ['+R','-R','-L','+L']]
                invTdSides=[['-L','+L','+R','-R'],
                            ['+R','-R','-L','+L'],
                            ['+R','-R','-L','-R'],
                            #['+R','+L','+R','-R'],
                            ['-L','+L','+R','+L']]
                # Get sequences before and after BND
                posList=[pos1,pos2,pos3,pos4]
                errPerc=','.join([fusions[i-1][11],
                                    fusions[i][11]])
                jointType=[fusions[i-1][12],
                            fusions[i][12]]
                jointLen=[str(fusions[i-1][13]),
                            str(fusions[i][13])]
                jointSeq=[str(fusions[i-1][14]),
                            str(fusions[i][14])]
                firstLastPart=[str(fusions[i-1][15]),
                                str(fusions[i][15])]
                seqsNearBND=[str(fusions[i-1][17]),
                                str(fusions[i][17]),
                                str(fusions[i-1][18]),
                                str(fusions[i][18])]
                if ss1[1]=='L':
                    if ss1[0]=='-':
                        posList=[str(pos4),str(pos1),str(pos3),str(pos2)]
                    else:
                        posList=[str(pos4),str(pos1),str(pos3),str(pos2)]
                    jointType=jointType[::-1]
                    jointLen=jointLen[::-1]
                    jointSeq=jointSeq[::-1]
                    firstLastPart=firstLastPart[::-1]
                    seqsBeforeBNDs=seqsNearBND[::-1][:2]
                    seqsAfterBNDs=seqsNearBND[::-1][2:]
                else:
                    if ss1[0]=='-':
                        posList=[str(pos1),str(pos4),str(pos2),str(pos3)]
                    else:
                        posList=[str(pos1),str(pos4),str(pos2),str(pos3)]
                    seqsBeforeBNDs=seqsNearBND[:2]
                    seqsAfterBNDs=seqsNearBND[2:]
                if ss1[0]=='-':
                    jointSeq=self.makeAllRevCompl(jointSeq)
                    seqsBeforeBNDs=self.makeAllRevCompl(seqsBeforeBNDs)
                    seqsAfterBNDs=self.makeAllRevCompl(seqsAfterBNDs)
                # Translocation
                twoJointLR=False
                # Check order of strands and sides
                if ((ssOrder in trlSides or
                        ssOrder[::-1] in trlSides) and
                    # Check if position of the 2nd part is out of the 1st fragment location
                    fusions[i-1][0]==fusions[i][4] and # chrom of 1st and 4th positions are the same
                    abs(pos4-pos1)<=self.distToJoinTrl and # distance between two positions flanking position into which TRL was inserted should be not so large
                    (fusions[i-1][0]!=fusions[i-1][4] or # chrom for 2nd and 3rd position are another
                        (pos2<min(pos1,fragEnd1,pos4,fragEnd2) and pos3<min(pos1,fragEnd1,pos4,fragEnd2) or # or those positions are out of range of both fragments
                                pos2>max(pos1,fragEnd1,pos4,fragEnd2) and pos3>max(pos1,fragEnd1,pos4,fragEnd2)))):
                    # We need to write it in the new form
                    # chrom1 - chromosome of position of genome where new fragment was added
                    # pos1 - position of genome where new fragment was added
                    # strand1 - it is always +, because we don't need to write it in opposite direction
                    # side1 - it was L or R, for translocations it will be T
                    # chrom2 - chromosome of the position of genome from which fragment was transfered
                    # pos2,pos3 - positions of genome from which fragment was transfered
                    # mq - mean mapping quality
                    # fraglen1 - length of one fragment around site of new fragment insertion
                    # fraglen2 - length of the 2nd fragment around site of new fragment insertion
                    # jointType- type of joinment, listed with separation by comma for two junctions
                    # jointLen - length of the joinment, listed with separation by comma for two junctions
                    # jointSeq - sequence of the joinment, listed with separation by comma for two junctions
                    # firstLastPart - is it first or the last part of the read (1 - 1st, 2 - last, 0 - other)
                    # partsNum - number of read parts
                    # seqBeforeBND - sequence before juntions, listed with separation by comma for two junctions
                    # seqAfterBND - sequence after juntions, listed with separation by comma for two junctions
                    # readName - Read_Name
                    mq=str(round((fusions[i-1][8]+fusions[i][8])/2,0))
                    fusionsJoint.append([fusions[i-1][0],','.join(posList[:2]),'+','T',
                                            fusions[i-1][4],','.join(posList[2:]),
                                            '+','T',mq,
                                            fusions[i-1][9], # fraglen1
                                            fusions[i-1][10], # fraglen2
                                            errPerc, # errPercents
                                            ','.join(jointType), # jointType
                                            ','.join(jointLen), # jointLen
                                            ','.join(jointSeq), # jointSeq
                                            # number was deleted from fusions
                                            ','.join(firstLastPart), # firstLastPart
                                            fusions[i-1][16], # partsNum
                                            ','.join(seqsBeforeBNDs), # seqBeforeBND
                                            ','.join(seqsAfterBNDs), # seqAfterBND
                                            fusions[i][19]])
                    ############################################################
                    ############### Check that some sequence was duplicated around TRL, e.g. chr1:1929454,1929364
                    # chr1:8,912,128-8,912,409
                    # chr1    8912269,8912173 +       T       chr15   43014532,43015654       +       
                    # T       49.0    6       216     1123    inter,inter     94,608
                    ############################################################
                    twoJointLR=True
                elif ((ssOrder in tdSides or
                        ssOrder[::-1] in tdSides) and
                        fusions[i-1][0]==fusions[i-1][4]==fusions[i][0]==fusions[i][4] and
                        abs(pos1-pos3)<100 and
                        abs(pos2-pos4)<100):
                    mq=str(round((fusions[i-1][8]+fusions[i][8])/2,0))
                    posList=[str(pos1),str(pos2),str(pos3),str(pos4)]
                    if ss1[1]=='L':
                        posList=posList[::-1]
                    # Calculate distance of the third fragment to the end of the repeat
                    # If it is far from it, we should subtract this distanca from the length of the TD
                    lenBias=max(0,int(posList[0])-int(fusions[i-1][9])-int(posList[1]),
                                int(posList[2])-(int(posList[3])+int(fusions[i][10])))
                    # If such fragment did not reach the end the repeat we need take only the small repeat length
                    # ++++++-------------++++++ - ref
                    # ++++++-------------++++++-------------++++++-------------++++++ - TD
                    #           ---------++++++ - third fragment
                    if lenBias>max([fusions[i-1][13],
                                    fusions[i][13]]):
                        lenBias=max([fusions[i-1][13],
                                        fusions[i][13]])
                    # Check, if the previous fusionJoint was also TD
                    # Then, it means that it is the extension of the previous TD
                    if (len(fusionsJoint)>0 and
                        fusionsJoint[-1][3][-1]=='D'):
                        # And then we work with the previous fusionJoint, but not this one
                        # Check, if this extension of previous TD has been already done
                        if len(fusionsJoint[-1][3])>1:
                            tdSize=int(fusionsJoint[-1][3].replace('D',''))
                            fusionsJoint[-1][3]=str(tdSize+1)+'D'
                            fusionsJoint[-1][7]=str(tdSize+1)+'D'
                        else:
                            fusionsJoint[-1][3]='2D'
                            fusionsJoint[-1][7]='2D'
                        fusionsJoint[-1][10]=fusions[i][10]
                        prevLenBias=int(fusionsJoint[-1][9].split(',')[1])
                        if prevLenBias<lenBias:
                            fusionsJoint[-1][9]=','.join([fusionsJoint[-1][9].split(',')[0],
                                                            str(lenBias)])
                        fusionsJoint[-1][11]+=','+errPerc
                        fusionsJoint[-1][12]+=','+','.join(jointType)
                        fusionsJoint[-1][13]+=','+','.join(jointLen)
                        fusionsJoint[-1][14]+=','+','.join(jointSeq)
                        fusionsJoint[-1][15]+=','+','.join(firstLastPart)
                        fusionsJoint[-1][17]+=','+','.join(seqsBeforeBNDs)
                        fusionsJoint[-1][18]+=','+','.join(seqsAfterBNDs)
                    else:
                        fusionsJoint.append([fusions[i-1][0],','.join(posList[:2]),'+','D',
                                                fusions[i-1][4],','.join(posList[2:]),'+','D',mq,
                                                ','.join([str(fusions[i-1][9]),
                                                        str(lenBias)]), # fraglen1 of fusion1
                                                fusions[i][10], # fraglen2 of fusion2
                                                errPerc, # errPercents
                                                ','.join(jointType), # jointType
                                                ','.join(jointLen), # jointLen
                                                ','.join(jointSeq), # jointSeq
                                                # number was deleted from fusions
                                                ','.join(firstLastPart), # firstLastPart
                                                fusions[i-1][16], # partsNum
                                                ','.join(seqsBeforeBNDs), # seqBeforeBND
                                                ','.join(seqsAfterBNDs), # seqAfterBND
                                                fusions[i][19]])
                    twoJointLR=True
                elif ((ssOrder in invSides or
                        ssOrder[::-1] in invSides) and
                        fusions[i-1][0]==fusions[i-1][4]==fusions[i][0]==fusions[i][4] and
                        pos1<=pos2<=pos4 and
                        pos1<=pos3<=pos4 and
                        abs(pos1-pos2)>=self.minVarLen and
                        abs(pos3-pos4)>=self.minVarLen):
                    mq=str(round((fusions[i-1][8]+fusions[i][8])/2,0))
                    fusionsJoint.append([fusions[i-1][0],','.join([str(min(pos1,pos4)),
                                                                    str(min(pos2,pos3))]),
                                            '+','I',
                                            fusions[i-1][4],','.join([str(max(pos2,pos3)),
                                                                    str(max(pos1,pos4))]),
                                            '+','I',mq,
                                            fusions[i-1][9], # fraglen1
                                            fusions[i-1][10], # fraglen2
                                            errPerc, # errPercents
                                            ','.join(jointType), # jointType
                                            ','.join(jointLen), # jointLen
                                            ','.join(jointSeq), # jointSeq
                                            # number was deleted from fusions
                                            ','.join(firstLastPart), # firstLastPart
                                            fusions[i-1][16], # partsNum
                                            ','.join(seqsBeforeBNDs), # seqBeforeBND
                                            ','.join(seqsAfterBNDs), # seqAfterBND
                                            fusions[i][19]])
                    # Add detection of deletions near the inversion
                    # If the distance between two fragments is more than the LR minimal length
                    # We save it as deletion
                    if abs(pos1-pos3)>=self.minVarLen:
                        fusionsJoint.append([fusions[i-1][0],str(min(pos1,pos3)),'+','R',
                                                fusions[i-1][4],str(max(pos2,pos3)),'+','L',mq,
                                                fusions[i-1][9], # fraglen1
                                                fusions[i-1][10], # fraglen2
                                                errPerc, # errPercents
                                                'inversion',
                                                0,'',
                                                # number was deleted from fusions
                                                0,
                                                fusions[i-1][15], # partsNum
                                                '', # seqBeforeBND
                                                '', # seqAfterBND
                                                fusions[i][18]])
                    if abs(pos2-pos4)>=self.minVarLen:
                        fusionsJoint.append([fusions[i-1][0],str(min(pos2,pos4)),'+','R',
                                                fusions[i-1][4],str(max(pos2,pos4)),'+','L',mq,
                                                fusions[i][9], # fraglen1
                                                fusions[i][10], # fraglen2
                                                errPerc, # errPercents
                                                'inversion',
                                                0,'',
                                                # number was deleted from fusions
                                                0,
                                                fusions[i-1][16], # partsNum
                                                '', # seqBeforeBND
                                                '', # seqAfterBND
                                                fusions[i][19]])
                    twoJointLR=True
                elif ((ssOrder in invTdSides or
                        ssOrder[::-1] in invTdSides) and
                        fusions[i-1][0]==fusions[i-1][4]==fusions[i][0]==fusions[i][4] and
                        abs(pos1-pos4)<100 and
                        (abs(pos1-pos2)<100 or
                        abs(pos3-pos4)<100)):
                    mq=str(round((fusions[i-1][8]+fusions[i][8])/2,0))
                    fusionsJoint.append([fusions[i-1][0],','.join(posList[:2]),
                                            '+','V',
                                            fusions[i-1][4],','.join(posList[2:]),
                                            '+','V',mq,
                                            fusions[i-1][9], # fraglen1
                                            fusions[i-1][10], # fraglen2
                                            errPerc, # errPercents
                                            ','.join(jointType), # jointType
                                            ','.join(jointLen), # jointLen
                                            ','.join(jointSeq), # jointSeq
                                            # number was deleted from fusions
                                            ','.join(firstLastPart), # firstLastPart
                                            fusions[i-1][16], # partsNum
                                            ','.join(seqsBeforeBNDs), # seqBeforeBND
                                            ','.join(seqsAfterBNDs), # seqAfterBND
                                            fusions[i][19]])
                    twoJointLR=True
                else:
                    fusionsJoint.append(fusions[i-1])
                    if i==len(fusions)-1:
                        fusionsJoint.append(fusions[i])
            if twoJointLR==True:
                fusions=fusionsJoint[:]
        # Write deletions and insertions
        for deletion in deletions:
            errPerc=int(round(deletion[4]*100/readAlignedLen,0))
            fusion=[deletion[0],str(deletion[1]),'+','R',
                    deletion[0],str(deletion[2]),'+','L',
                    stat.median(deletion[3]),
                    deletion[6],readAlignedLen-deletion[6]-(deletion[2]-deletion[1]+1),errPerc,'cigar',0,'']
            fusions.append(fusion+[0,
                                    len(readParts),
                                    readSeq[max(0,deletion[5]-30):deletion[5]],
                                    readSeq[deletion[5]:deletion[5]+30],
                                    readName])
            anyFusionsFound=True
        for fusion in fusions:
            # [chr1,pos1,strand1,side1,
            #  chr2,pos2,strand2,side2,
            #  MQ,refCov,readLen1,readLen2,
            #  type(intersection,equal, or distance),
            #  length,sequence,number,
            #  firstLastPart,partsNum,
            #  seqBeforeBND,seqAfterBND]
            out=[]
            for fusionInfoNum in range(8):
                out.append(str(fusion[fusionInfoNum])) 

            ###########################################################################################
            ####### TO DO: Analyze coverage statistics to get additional characteristics for DELs #####
            ###########################################################################################
            # MQ
            out.append(str(fusion[8]))
            pos1=int(str(fusion[1]).split(',')[0])
            if (fusion[3]=='T' or
                'D' in fusion[3] or
                fusion[3]=='V'):
                pos2=pos1
            elif fusion[3]=='I':
                pos2=int(str(fusion[5]).split(',')[1])
            else:
                pos2=int(fusion[5])
            # Get reference coverage for the left coordinate
            if (fusion[3]=='R' or
                fusion[3]=='T' or
                'D' in fusion[3] or
                fusion[3]=='I' or
                fusion[3]=='V'):
                leftCoord=self.minVarLen*math.ceil(pos1/self.minVarLen)-self.minVarLen
            elif fusion[3]=='L':
                leftCoord=self.minVarLen*math.floor(pos1/self.minVarLen)+self.minVarLen
            # If fusion's chromosome is not the same as the current input
            leftCoordSimple=int(math.floor(leftCoord/self.minVarLen))
            if (fusion[0] in self.genomeCoverage.keys() and
                leftCoordSimple<len(self.genomeCoverage[fusion[0]])):
                leftCov=self.genomeCoverage[fusion[0]][leftCoordSimple]#[0]
            else:
                leftCov=0
            # Get reference coverage for the right coordinate
            if fusion[7]=='R':
                rightCoord=self.minVarLen*math.ceil(pos2/self.minVarLen)-self.minVarLen
            elif (fusion[7]=='L' or
                    fusion[7]=='T' or
                    'D' in fusion[7] or
                    fusion[7]=='I' or
                    fusion[7]=='V'):
                rightCoord=self.minVarLen*math.floor(pos2/self.minVarLen)+self.minVarLen
            # If fusion's chromosome is not the same as the current input
            rightCoordSimple=int(math.floor(rightCoord/self.minVarLen))
            if (fusion[4] in self.genomeCoverage.keys() and
                rightCoordSimple<len(self.genomeCoverage[fusion[4]])):
                rightCov=self.genomeCoverage[fusion[4]][rightCoordSimple]#[0]
            else:
                rightCov=0
            # For intrachromosomal we get maximal
            if fusion[0]==fusion[4]:
                out.append(str(max(leftCov,rightCov)))
            # For interchromosomal, we get sum
            else:
                out.append(str(max([leftCov,rightCov])))
            try:
                out.append(str(fusion[9]))
            except TypeError:
                print('ERROR!',fusion)
                print(fusion[9])
                print(fusion[10])
                exit()
            for fusionInfoNum in range(10,20):
                out.append(str(fusion[fusionInfoNum]))
            outFiles[0].write('\t'.join(out)+'\n')
        # chrom,curRefPos,ciLen,mapQual,insertedSeq,readName
        for insertion in insertions:
            out=[]
            out.extend([insertion[0],
                        str(insertion[1])])
            leftCoord=self.minVarLen*math.ceil(insertion[1]/self.minVarLen)-self.minVarLen
            leftCoordSimple=int(math.floor(leftCoord/self.minVarLen))
            rightCoord=self.minVarLen*math.floor(insertion[1]/self.minVarLen)+self.minVarLen
            rightCoordSimple=int(math.floor(rightCoord/self.minVarLen))
            if (readChrom in self.genomeCoverage.keys() and
                leftCoordSimple<len(self.genomeCoverage[readChrom])):
                leftCov=self.genomeCoverage[readChrom][leftCoordSimple]#[0]
            else:
                leftCov=0
            if (readChrom in self.genomeCoverage.keys() and
                rightCoordSimple<len(self.genomeCoverage[readChrom])):
                rightCov=self.genomeCoverage[readChrom][rightCoordSimple]#[0]
            else:
                rightCov=0
            errPerc=int(round(insertion[4]*100/readAlignedLen,0))
            out.extend([str(insertion[2]),
                        insertion[5],
                        '1',str(max([leftCov,
                                        rightCov])),
                        str(stat.median(insertion[3])),
                        str(errPerc),
                        insertion[6]])
            outFiles[1].write('\t'.join(out)+'\n')
        self.thsize-=1

if __name__=='__main__':        

    # Section of reading arguments
    par=argparse.ArgumentParser(description='This script calls large rearrangements in ONT data')
    par.add_argument('--bam-file','-bam',
                     dest='bamFile',type=str,
                     help='BAM file with NGS reads',
                     required=True)
    par.add_argument('--chromosome','-chr',
                     dest='chrom',type=str,
                     help='chromosome to analyze. Default: all chromosomes',
                     required=False)
    par.add_argument('--start','-start',
                     dest='start',type=int,
                     help='Start coordinate for extracting reads. Default: whole chromosome or whole genome',
                     required=False,default=None)
    par.add_argument('--end','-end',
                     dest='end',type=int,
                     help='End coordinate for extracting reads. Default: whole chromosome or whole genome',
                     required=False,default=None)
    par.add_argument('--minimal-length','-len',
                     dest='minVarLen',type=int,
                     help='minimal acceptable length of variant (for InDels). Default: 50',
                     required=False,default=50)
    par.add_argument('--minimal-clipped-length','-clip',
                     dest='minClipLen',type=int,
                     help='minimal acceptable length of a clipped read part. Default: 100',
                     required=False,default=100)
    par.add_argument('--dist-to-join-trl','-dist',
                     dest='distToJoinTrl',type=int,
                     help='maximal acceptable distance between two coordinates. Default: 1000',
                     required=False,default=1000)
    par.add_argument('--threads','-th',
                     dest='threads',type=int,
                     help='number of threads to call LRs. Default: 4',
                     required=False,default=4)
    #par.add_argument('--output-file','-out',
    #                 dest='outFile',type=str,
    #                 help='VCF-file for output',
    #                 required=False)
    args=par.parse_args()
    
    ontc=ONTLRCaller(args.bamFile,args.chrom,
                     args.start,args.end,args.threads,
                     args.minVarLen,args.minClipLen,
                     args.distToJoinTrl)
    ontc.readBamFile()
    
    
