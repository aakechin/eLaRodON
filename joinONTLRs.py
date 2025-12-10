# This script joins called LRs

# import argparse
from multiprocessing import Value,Pool
from typing import Type
from xml.dom.minidom import Attr
import pandas
import glob
import sys
import statistics
import traceback
from tqdm import tqdm
import logging

class JoinLR():

    def __init__(self,ds,
                 args):
        self.th=args.threads
        self.log_path = args.workDir+'/elarodon.log'
        logging.basicConfig(level=logging.DEBUG, filename=self.log_path, filemode="a",
                            format="%(asctime)s %(funcName)s %(lineno)d %(message)s")
        p=Pool(self.th)
        self.maxDistToJoin=args.maxDistToJoin
        self.allFusions={}
        self.allInsertions={}
        for res in tqdm(p.imap_unordered(self.processFile,ds,
                                         chunksize=10),
                        total=len(ds)):
            fileFusions,fileInsertions=res
            if len(fileFusions)>0:
                for chrom,fusionToInfo in fileFusions.items():
                    if chrom not in self.allFusions.keys():
                        self.allFusions[chrom]=fusionToInfo
                    else:
                        for fusion,info in fusionToInfo.items():
                            if fusion not in self.allFusions[chrom].keys():
                                self.allFusions[chrom][fusion]=info
                            else:
                                for i in [0,2,4,5,6,7,8,9,10]:
                                    self.allFusions[chrom][fusion][i].extend(info[i])
                                    # Error percents
                                for jointType,jointSeqToNum in enumerate(info[1]):
                                    for jointSeq,num in jointSeqToNum.items():
                                        if jointSeq not in self.allFusions[chrom][fusion][1][jointType].keys():
                                            self.allFusions[chrom][fusion][1][jointType][jointSeq]=num
                                        else:
                                            self.allFusions[chrom][fusion][1][jointType][jointSeq]+=num
            if len(fileInsertions)>0:
                for chrom1_start,info in fileInsertions.items():
                    if chrom1_start not in self.allInsertions.keys():
                        self.allInsertions[chrom1_start]=info
                    else:
                        for insInfo in info:
                            self.allInsertions[chrom1_start].append(insInfo)

    def processFile(self,d):
        allFusions={}
        allInsertions={}
        if 'fusions.csv' in d:
            try:
                data=pandas.read_csv(d,sep='\t',
                                     header=0,
                                     keep_default_na=False,
                                     low_memory=False,
                                     encoding='utf-8')
            except UnicodeDecodeError:
                print('ERROR (6)! Could not read the following file (UnicodeDecodeError):')
                print(d)
                exit(6)
            for i,row in data.iterrows():
                row = row.to_list()
                fusion=row[:8]
                # Compare chromosomes
                fusionOrdered=fusion[:]
                try:
                    chrom1=int(fusion[0].replace('chr',''))
                except TypeError:
                    chrom1=fusion[0]
                except ValueError:
                    chrom1=fusion[0]
                try:
                    chrom2=int(fusion[4].replace('chr',''))
                except TypeError:
                    chrom2=fusion[4]
                except ValueError:
                    chrom2=fusion[4]
                except AttributeError:
                    print('AttribureError:',fusion[4])
                    print(d)
                    print(row)
                    print(fusion)
                    exit()
                readLen1=row[10]
                readLen2=row[11]
                # Change the order only if it is not T, D, or I
                if (row[3] not in ['T','D','I','V'] and
                    'D' not in row[3]):
                    if ((type(chrom1)==int and
                         type(chrom2)==int) or
                        (type(chrom1)==str and
                          type(chrom2)==str)):
                        # If chrom2 is less than chrom1
                        if chrom2<chrom1:
                            # We change their order
                            fusionOrdered=[*fusion[4:8],
                                           *fusion[0:4]]
                            readLen1=row[11]
                            readLen2=row[10]
                        # If pos2<pos1, we also change the order
                        elif (chrom1==chrom2 and 
                              checkPosOrder(fusion[1],
                                            fusion[5],
                                            fusion[3])):
                            fusionOrdered=[*fusion[4:8],
                                           *fusion[0:4]]
                            readLen1=row[11]
                            readLen2=row[10]
                    # If the 1st chromosome is str and 2nd is int
                    elif (type(chrom1)==str and
                          type(chrom2)==int):
                        # We also change their order
                        fusionOrdered=[*fusion[4:8],
                                       *fusion[0:4]]
                        readLen1=row[11]
                        readLen2=row[10]
                errPercs=str(row[12]).split(',')
                # Start_End_Locations,Number_of_parts
                startEndLoc=row[16]
                partsNum=int(row[17])
                # SeqBefore,SeqAfter
                seqBefore=row[18]
                seqAfter=row[19]
                fusionOrdered=tuple(fusionOrdered)
                chrom=fusionOrdered[0]
                if type(fusionOrdered[1])==int:
                    pos=fusionOrdered[1]
                else:
                    pos=int(fusionOrdered[1].split(',')[0])
                # Now save it by the 1st chromosome
                if chrom not in allFusions.keys():
                    allFusions[chrom]={}
                if fusionOrdered not in allFusions[chrom].keys():
                    # [reads_names,[jointType1,jointType2,jointType3],
                    # [MQs],refCov,[readLens1],[readLens2],[errPercs]
                    # [First_Last_Part],[Num_of_Parts],[Seq_before],[Seq,after]]
                    allFusions[chrom][fusionOrdered]=[[],[{},{},{},{},{}],
                                                           [],row[9],[],[],[],
                                                           [],[],[],[]]
                if row[20] not in allFusions[chrom][fusionOrdered][0]:
                    # read_name
                    allFusions[chrom][fusionOrdered][0].append(row[20])
                    # MQ
                    allFusions[chrom][fusionOrdered][2].append(row[8])
                    # Read lengths
                    allFusions[chrom][fusionOrdered][4].append(str(readLen1))
                    allFusions[chrom][fusionOrdered][5].append(int(readLen2))
                    # Error percents
                    allFusions[chrom][fusionOrdered][6].extend(errPercs)
                    # Parts of reads (numbers)
                    allFusions[chrom][fusionOrdered][7].append(str(startEndLoc))
                    allFusions[chrom][fusionOrdered][8].append(str(partsNum))
                    # Seq before and after BND
                    allFusions[chrom][fusionOrdered][9].append(seqBefore)
                    allFusions[chrom][fusionOrdered][10].append(seqAfter)
                    # If the sequence of jointType is absent
                    if (row[15]=='' or
                        type(row[15])==float):
                        # Then we save its length
                        jointSeq=str(row[14])
                    else:
                        jointSeq=row[15]
                    # Determine, to which list we should save the joint info
                    for jointType,jointSeq in zip(str(row[13]).split(','),
                                                  str(row[14]).split(',')):
                        jointNum=['inter','scarless','newseq','cigar','inversion'].index(jointType)
                        if jointSeq not in allFusions[chrom][fusionOrdered][1][jointNum].keys():
                            allFusions[chrom][fusionOrdered][1][jointNum][jointSeq]=0
                        allFusions[chrom][fusionOrdered][1][jointNum][jointSeq]+=1
            del(data)
            return(allFusions,allInsertions)
        if 'insertions.csv' in d:
            data=pandas.read_csv(d,sep='\t')
            sorted_data = data.sort_values(by='Pos1')

            miniGroups = []
            miniGroupInds = {}

            miniGroupSeqs = {}
            
            for i,row in sorted_data.iterrows():
                row = row.to_list()
                
                ins_added = False 
                # print(row)

                # row
                # 0 - chrom, 1 - pos, 2 - len ins, 3 - sequence, 4 - read number, 5 - total cov, 6 - MQ, 7 - Read Muts1, 8 - read name 
                
                if miniGroups != []:
                    for ind_last_group, last_ins in sorted(miniGroupInds.items(), key=lambda x: x[0], reverse=True):
                        # print('ind_last_group', ind_last_group)
                        if (ins_added == False and 
                            row[0] == last_ins[0] and
                            abs(int(row[1]) - int(last_ins[1])) <= 100 and
                            abs(int(row[2]) - int(last_ins[2])) <= 100):
                            
                            # minigroup
                            # 0 - chrom, 1 - pos, 2 - len, 3 - total cov, 4 - MQ, 5 - Read_Muts1, 6 - read name

                            miniGroups[ind_last_group][0].append((row[0], # chrom 0
                                                                int(row[1]), # pos 1
                                                                int(row[2]), # len 2
                                                                int(row[5]), # tot cov 3
                                                                int(row[6]), # MQ 4
                                                                float(row[7]), # RM1 5
                                                                row[8])) # RN 6
                            miniGroupInds[ind_last_group] = [row[0], int(row[1]), int(row[2])] # add new last ins in group
                            
                            if row[3] not in miniGroups[ind_last_group][1]:
                                miniGroups[ind_last_group][1][row[3]] = 0                    
                            miniGroups[ind_last_group][1][row[3]] += 1
                            # print('We add INS to ', ind_last_group, ' group')
                            # print(miniGroups[ind_last_group])
                            # print(miniGroupInds[ind_last_group])
                            # input()
                            ins_added = True
                    
                else:
                    # minigroups and minigroup sequences
                    miniGroups.append([[], {}])
                    miniGroups[-1][0].append((row[0], # chrom 0
                                            int(row[1]), # pos 1
                                            int(row[2]), # len 2
                                            int(row[5]), # tot cov 3
                                            int(row[6]), # MQ 4
                                            float(row[7]), # RM1 5
                                            row[8])) # RN 6
                    miniGroupInds[len(miniGroups)-1] = [row[0], int(row[1]), int(row[2])]

                    if row[3] not in miniGroups[-1][1]:
                        miniGroups[-1][1][row[3]] = 0                    
                    miniGroups[-1][1][row[3]] += 1
                    ins_added = True
                
                if ins_added == False:

                    miniGroups.append([[], {}])
                    miniGroups[-1][0].append((row[0], # chrom 0
                                            int(row[1]), # pos 1
                                            int(row[2]), # len 2
                                            int(row[5]), # tot cov 3
                                            int(row[6]), # MQ 4
                                            float(row[7]), # RM1 5
                                            row[8])) # RN 6
                    miniGroupInds[len(miniGroups)-1] = [row[0], int(row[1]), int(row[2])]

                    if row[3] not in miniGroups[-1][1]:
                        miniGroups[-1][1][row[3]] = 0                    
                    miniGroups[-1][1][row[3]] += 1
                
                # print('ins_added', ins_added)
                # print('miniGroups', miniGroups)
                # input()

            for miniGroup in miniGroups:

                # print('miniGroup', miniGroup)
                # input()

                chrom1_start = miniGroup[0][0][0]
                pos1_start = miniGroup[0][0][1]
                total_cov_middle = round(sum([x[3] for x in miniGroup[0]]) / len(miniGroup[0]))      
                len_ins_middle = round(sum([x[2] for x in miniGroup[0]]) / len(miniGroup[0]))
                read_muts = [str(x[5]) for x in miniGroup[0]]
                read_names = [x[6] for x in miniGroup[0]]
                mqs = [x[4] for x in miniGroup[0]]

                miniGroupSeqs = miniGroup[1]
                
                if chrom1_start not in allInsertions.keys():
                    allInsertions[chrom1_start]=[]
                allInsertions[chrom1_start].append([pos1_start,
                                                    read_names, 
                                                    miniGroupSeqs,
                                                    total_cov_middle,
                                                    mqs,
                                                    read_muts])
        del(miniGroups)
        del(data)
        return(allFusions,allInsertions)

    def joinAllSimilarLRs(self):
        print('Joining similar fusions...')
        p=Pool(self.th)
        self.allFusionsJoined1={}
        for res in tqdm(p.imap_unordered(self.joinChromSimilarLRs,
                                         self.allFusions.keys(),
                                         chunksize=10),
                        total=len(self.allFusions.keys())):
            self.allFusionsJoined1[res[0]]=res[1]
        

    def joinChromSimilarLRs(self,chrom):
        chromFusionsJoined1={}
        # Sort by 1st coordinate
        fusionsSorted=sorted(self.allFusions[chrom],
                             key=lambda item:[item[0],
                                              int(str(item[1]).split(',')[0]), # pos
                                              item[4],
                                              int(str(item[5]).split(',')[0])])
        # Stores fusion IDs that have been already joined to another fusion
        alreadyJoinedFusionIDs=set()
        # Join fusions by searching for the first coordinate
        for i,fusion in enumerate(fusionsSorted):
            if (i==len(fusionsSorted)-1 or
                i in alreadyJoinedFusionIDs):
                continue
            joinedNum=self.allFusions[chrom][fusion]
            joinedFusion=fusion[:]
            fusionNum=i+1
            # Compare coordinates of the current fusion and several further
            while(True):
                if fusionNum==len(fusionsSorted)-1:
                    break
                nextFusion=fusionsSorted[fusionNum]
                # If chromosomes are the same
                if (fusion[0]==nextFusion[0] and
                    fusion[4]==nextFusion[4] and
                    # And strands are corresponding
                    # [+,+] vs [+,+] or [-,-] | [-,-] vs [+,+] or [-,-]
                    # [+,-] vs [+,-] or [-,+]
                    (((fusion[2]==fusion[6] and
                       nextFusion[2]==nextFusion[6])) or
                     (fusion[2]!=fusion[6] and
                      nextFusion[2]!=nextFusion[6])) and
                    (fusion[3]==nextFusion[3] and
                     fusion[7]==nextFusion[7]) and
                    # And the coordinates are similar
                    compareCoords(fusion[1],fusion[5],
                                  nextFusion[1],nextFusion[5],
                                  self.maxDistToJoin)):
                    # Join total numbers for the joint fusions
                    joinedNum[0].extend(self.allFusions[chrom][nextFusion][0])
                    joinedNum[2].extend(self.allFusions[chrom][nextFusion][2])
                    joinedNum[4].extend(self.allFusions[chrom][nextFusion][4])
                    joinedNum[5].extend(self.allFusions[chrom][nextFusion][5])
                    joinedNum[6].extend(self.allFusions[chrom][nextFusion][6])
                    joinedNum[7].extend(self.allFusions[chrom][nextFusion][7])
                    joinedNum[8].extend(self.allFusions[chrom][nextFusion][8])
                    joinedNum[9].extend(self.allFusions[chrom][nextFusion][9])
                    joinedNum[9].extend(self.allFusions[chrom][nextFusion][10])
                    # Join numbers for different joint types
                    for typeNum in range(5):
                        # Go through all sequences of the nextFusion
                        for nextFusionSeq,seqNum in self.allFusions[chrom][nextFusion][1][typeNum].items():
                            if nextFusionSeq not in joinedNum[1][typeNum].keys():
                                joinedNum[1][typeNum][nextFusionSeq]=0
                            joinedNum[1][typeNum][nextFusionSeq]+=seqNum
                    alreadyJoinedFusionIDs.add(fusionNum)
                elif (fusion[0]!=nextFusion[0] or
                      abs(int(str(fusion[1]).split(',')[0])-int(str(nextFusion[1]).split(',')[0]))>self.maxDistToJoin):
                    break
                fusionNum+=1
            if joinedFusion not in chromFusionsJoined1.keys():
                chromFusionsJoined1[joinedFusion]=joinedNum
            else:
                chromFusionsJoined1[joinedFusion][0].extend(joinedNum[0])
                for typeNum in range(5):
                    # Go through all sequences of the nextFusion
                    for nextFusionSeq,seqNum in joinedNum[1][typeNum].items():
                        if nextFusionSeq not in chromFusionsJoined1[joinedFusion][1][typeNum].keys():
                            chromFusionsJoined1[joinedFusion][1][typeNum][nextFusionSeq]=0
                        chromFusionsJoined1[joinedFusion][1][typeNum][nextFusionSeq]+=seqNum
                chromFusionsJoined1[joinedFusion][2].extend(joinedNum[2])
        return(chrom,chromFusionsJoined1)

    def writeFusionToOutput(self,outFile):
        print('Writing fusions to output...')
        rFile=open(outFile,'w')
        rFile.write('\t'.join(['Chrom1','Pos1','Strand1','Junction_Side1',
                               'Chrom2','Pos2','Strand2','Junction_Side2',
                               'Read_Number','Ref_Coverage','MQ_median',
                               'Read_Len1','Read_Len2','Read_Muts',
                               'Inter_Joints', 
                               'Scarless_Joints',
                               'NewSeq_Joints',
                               'Cigar_Num','Inversion_Num',
                               'Start_End_Locations','Number_of_parts',
                               'Seqs_before_BND','Seqs_after_BND',
                               'Read_Names'])+'\n')
        sumLen=sum(len(item.keys()) for item in self.allFusionsJoined1.values())
        showPercWork(0,sumLen)
        fusioNum=0
        # Stores GC-content values for LRs of different groups:
        # inter_lambda, equal_lambda, distance_lambda, inter_dif_chroms, equal_dif_chroms, dist_dif_chroms, 
        # inter_same_chrom, equal_same_chrom, dist_same_chrom
        gcContents=[[[],[],[]],
                    [[],[],[]],
                    [[],[],[]]]
        # Stores GC-content values for LRs with different coverage values:
        # [cov==1,cov==2,cov==3...cov==9,cov>=10]
        # It doesn't include LRs with distance
        gcContentForCovs=[[],[],[],
                          [],[],[],
                          [],[],[],[]]
        for chrom in tqdm(sorted(self.allFusionsJoined1.keys()),
                          total=len(self.allFusionsJoined1.keys())):
            for fusion,numbers in self.allFusionsJoined1[chrom].items():
                # Stores length of the joint or its sequence
                jointTypeInfo=[]
                for typeNum in range(5):
                    if typeNum>2:
                        jointTypeInfo.append(str(sum(numbers[1][typeNum].values())))
                        continue
                    jointSeqLens=[]
                    for jointSeqLen,num in numbers[1][typeNum].items():
                        jointSeqLens.extend([str(jointSeqLen)]*num)
                    if len(jointSeqLens)==0:
                        jointTypeInfo.append('0')
                    else:
                        jointTypeInfo.append(','.join(jointSeqLens))
                try:
                    out='\t'.join([fusion[0],str(fusion[1]),
                                   fusion[2],fusion[3],
                                   fusion[4],str(fusion[5]),
                                   fusion[6],fusion[7],
                                   str(len(numbers[0])),
                                   str(int(numbers[3])),
                                   str(statistics.median(numbers[2])),
                                   ','.join(numbers[4]),
                                   str(max(numbers[5])),
                                   ','.join(numbers[6]),
                                   *jointTypeInfo,
                                   ','.join(numbers[7]),
                                   ','.join(numbers[8]),
                                   ','.join(numbers[9]),
                                   ','.join(numbers[10]),
                                   ','.join(numbers[0])])
                except TypeError:
                    print('ERROR (2)! Incorrect type of data for output:')
                    print(jointTypeInfo)
                    print(fusion)
                    print(numbers)
                    exit(2)
                rFile.write(out+'\n')
        rFile.close()
        print()

    def writeInsertionsToOutput(self,outFile):
        # Insertions
        print('Writing insertions to output...')
        rFile=open(outFile,'w')
        rFile.write('\t'.join(['Chrom','Pos','Insertion_Length','Sequence',
                               'Read_Number','Total_coverage','MQ_median',
                               'Read_Names','Read_Muts'])+'\n')
        for chrom,posToInsertions in tqdm(sorted(self.allInsertions.items()),
                                          total=len(self.allInsertions.keys())):
            for insInfo in sorted(posToInsertions, key=lambda x: x[0]):
                pos = insInfo[0]
                insInfo = insInfo[1:]
                outLineList=[str(chrom),
                             str(pos)]
                # Get the most frequent seq
                # [altCov,{insSeqs:1},totalCov,[MQs]]
                for seq,num in sorted(insInfo[1].items(),
                                      key=lambda item:item[1],
                                      reverse=True):
                    outLineList.append(str(len(seq)))
                    outLineList.append(seq)
                    break
                outLineList.extend([str(len(insInfo[0])),
                                    str(insInfo[2]),
                                    str(statistics.median(insInfo[3])),
                                    ','.join(insInfo[0]),
                                    ','.join(insInfo[4])])
                rFile.write('\t'.join(outLineList)+'\n')
        rFile.close()

def showPercWork(done,allWork):
    percDoneWork=round((done/allWork)*100,2)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

def getGC(seq):
    cNum=seq.count('C')
    gNum=seq.count('G')
    gcContent=round((cNum+gNum)*100/len(seq),0)
    return(gcContent)

def compareCoords(pos1,pos2,
                  nextPos1,nextPos2,
                  maxDist):
    pos1=str(pos1).split(',')
    pos2=str(pos2).split(',')
    nextPos1=str(nextPos1).split(',')
    nextPos2=str(nextPos2).split(',')
    if (len(pos1)!=len(nextPos1) or
        len(pos2)!=len(nextPos2)):
        return(False)
    for p1,np1 in zip(pos1,nextPos1):
        if abs(int(p1)-int(np1))>maxDist:
            return(False)
    for p2,np2 in zip(pos2,nextPos2):
        if abs(int(p2)-int(np2))>maxDist:
            return(False)
    return(True)

# Returns True if pos2<pos1
def checkPosOrder(pos1,pos2,juncType='L'):
    if juncType!='L' and juncType!='R':
        return(False)
    pos1=str(pos1).split(',')
    pos2=str(pos2).split(',')
    if int(pos2[0])<int(pos1[0]):
        return(True)
    return(False)

# if __name__=='__main__':

#     # Section of reading arguments
#     par=argparse.ArgumentParser(description='This script joined called LRs')
#     par.add_argument('--input-files','-in',
#                      dest='inFiles',type=str,
#                      help='regular expression for CSV-files with called LRs',
#                      required=True)
#     par.add_argument('--maximal-distance-join','-join',
#                      dest='maxDistToJoin',type=int,
#                      help='maximal acceptable distance of two neighboring fusions to join them. Default: 30',
#                      required=False,default=30)
#     par.add_argument('--output-insertions-file','-outi',
#                      dest='outInsFile',type=str,
#                      help='CSV output file for insertions',
#                      required=True)
#     par.add_argument('--output-file','-out',
#                      dest='outFile',type=str,
#                      help='CSV file for other large rearrangements',
#                      required=True)
#     par.add_argument('--threads','-th',
#                      dest='threads',type=int,
#                      help='number of threads to process files. Default: 8',
#                      required=False,default=8)
#     args=par.parse_args()

#     ds=glob.glob(args.inFiles)
#     if len(ds)==0:
#         print('ERROR (1)! No files were chosen:')
#         print(args.inFiles)
#         exit(1)
#     print('The number of files chosen:',len(ds))

#     print('Reading input files with',args.threads,'threads...')
#     jl=JoinLR(ds,args)
    
#     print('The total number of fusions:',sum(len(item) for item in jl.allFusions.values()))
#     print('The total number of insertions:',sum(len(item) for item in jl.allInsertions.values()))

#     jl.joinAllSimilarLRs()
#     jl.writeFusionToOutput(args.outFile)
#     jl.writeInsertionsToOutput(args.outInsFile)

##    print('Writing GC-content of different LRs...')
##    # Output-file for GC-content of different LRs
##    # Dist lambda | equal lambda | inter lambda...
##    rFileGC=open(args.outFile[:args.outFile.rfind('.')]+'.GC_content.csv','w')
##    rFileGC.write('\t'.join(['Chromosomes',
##                             'GC-content',
##                             'Joint Type'])+'\n')
##    # Go through lambda, dif_chroms, same_chroms
##    for i,chromType in enumerate(['Lambda','Dif_Chrom','Same_Chrom']):
##        for j,readOverType in enumerate(['Inter','Equal','Dist']):
##            for value in gcContents[i][j]:
##                rFileGC.write('\t'.join([chromType,
##                                         str(value),
##                                         readOverType])+'\n')
##    rFileGC.close()
##    print('Writing GC-content for LRs with different coverage values...')
##    # Output-file for GC-content for LRs with different coverage values
##    # Coverage | GC-content
##    rFileGC=open(args.outFile[:args.outFile.rfind('.')]+'.GC_content_vs_cov.csv','w')
##    rFileGC.write('\t'.join(['Coverage',
##                             'GC-content'])+'\n')
##    # Go through all coverage values
##    for i in range(10):
##        for value in gcContentForCovs[i]:
##            if i==9:
##                rFileGC.write('\t'.join(['>=10',
##                                         str(value)])+'\n')
##            else:
##                rFileGC.write('\t'.join([str(i+1),
##                                         str(value)])+'\n')
##    rFileGC.close()
