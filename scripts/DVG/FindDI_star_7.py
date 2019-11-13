#finding DIs from the star alignment information

import os
import glob
import argparse
import sys
import subprocess
from subprocess import call
import pandas as pd
import numpy as np
import re

"""
To run:
python3 FindDI_star_4.py --strain CAL09 --file_direct N_files


Only considering two possibilities:
MNM (one gap) and MNMNM (two gaps)
Insertions are not considered when finding gap start and stop positions

Mismatches/Deletions/Insertions are only accepted to the total mismatch allowed (--total-mismatch).
This will quantify the number of mismatch before and after the gap - if it is
greater than what we allow- then it isn't considered in our candidate list.

INPUT: 'N' files from split read alignment
OUTPUT:
1 file with coordinates using the CDS
1 file with coordinates using the full length genome
1 file where two gaps were called
"""

parser = argparse.ArgumentParser()
parser.add_argument('--strain','-s',required=True,help='Indicate strain- necessary for size calculations can only be CAL09,PR8,H3N2')
parser.add_argument('--ref','-r',help='Give directory where reference is located') #args.ref
parser.add_argument('--file_direct','-f',required=True,help='Give directory path for where N_files are located')
parser.add_argument('--align_length','-l',type=int,default=10,help='Give the length that each portion of the read must be (default 15 bp)')
parser.add_argument('--gap_size','-g',type=int,default=1,help='Give the minimum size of the gap (default 1)')
parser.add_argument('--total_mismatch','-m',type=int,default=10,help='Only 3 mismatches allowed in a read called as a DVG this includes mismatches/insertions/deletions')
args = parser.parse_args()

directory_list = ['{0}'.format(args.file_direct)]
for direc in directory_list:
    if not os.path.exists(direc):
        os.makedirs(direc)

def GrabLengths(strain):
    """
    INPUT: Strain name
    OUTPUT: lengths for full length and cds for each segment
    """
    if strain.upper() == 'CAL09':
        Full_seg = {"PB2":2280,"PB1":2274,"PA":2151 ,"HA":1701, "NP":1497,"NA":1410,"MP":982 ,"NS":863}
        CDS_seg = {"PB2":2280,"PB1":2274,"PA":2151,'HA':1701,'NP':1497,'NA':1410,'MP':982,'NS':838}

    elif strain.upper() == 'PR8':
        Full_seg = {"PB2":2341,"PB1":2341,"PA":2233 ,"HA":1778, "NP":1565,"NA":1423,"MP":1027 ,"NS":890}
        CDS_seg = {"PB2":2280,"PB1":2274,"PA":2151,'HA':1701,'NP':1497,'NA':1410,'MP':982,'NS':838}

    elif strain.upper()== 'H3N2':
        Full_seg = {"PB2":2341,"PB1":2341,"PA":2233 ,"HA":1778, "NP":1565,"NA":1413,"MP":1027 ,"NS":890}
        CDS_seg = {"PB2":2280,"PB1":2274,"PA":2151,'HA':1701,'NP':1497,'NA':1410,'MP':982,'NS':838}

    elif strain.upper() == 'YAMA': #added 11/26/2018
        Full_seg = {"YAMA_PB2":2313,"YAMA_PB1":2259,"YAMA_PA":2181 ,"YAMA_HA":1755, "YAMA_NP":1683,"YAMA_NA":1401,"YAMA_MP":1076 ,"YAMA_NS":1024}
        CDS_seg = {"YAMA_PB2":2313,"YAMA_PB1":2259,"YAMA_PA":2181,'YAMA_HA':1755,'YAMA_NP':1683,'YAMA_NA':1401,'YAMA_MP':1076,'YAMA_NS':1024}

    elif strain.upper() == 'VICT': #added 11/26/2018
        Full_seg = {"VICT_PB2":2313,"VICT_PB1":2259,"VICT_PA":2181 ,"VICT_HA":1758, "VICT_NP":1683,"VICT_NA":1401,"VICT_MP":1076 ,"VICT_NS":1027}
        CDS_seg = {"VICT_PB2":2313,"VICT_PB1":2259,"VICT_PA":2181,'VICT_HA':1758,'VICT_NP':1693,'VICT_NA':1401,'VICT_MP':1076,'VICT_NS':1027}

    elif strain.upper() == 'ZIKA': #added 5/30/2019
        Full_seg = {"MR766":10260}
        CDS_seg = {"MR766":10260}

    return Full_seg, CDS_seg

def ORF_check(CDS_estimated_size):
    """
    INPUT: the estimated size using the CDS coordinates
    OUTPUT: Whether the estimated size is divisible by 3 and therefore maintaining open read frame
    """
    if CDS_estimated_size%3 == 0: #if the remainder is zero - ORF maintained
        ORF_flag = 'T'
    elif CDS_estimated_size %3 != 0: #elif the remainder is not zero
        ORF_flag = "F"
    return ORF_flag

def EstSize(segment,gap_size,strain):#change if CAL09, PR8, H3N2
    """
    INPUT: Segment, the measured gap size (from N_file), strain
    OUTPUT: The estimated size of the DVG using both the CDS and full segment coordinates
    """
    Full_seg, CDS_seg = GrabLengths(strain)
    Full_length = Full_seg[segment]
    CDS_length = CDS_seg[segment]
    Full_size = Full_length - int(gap_size)
    CDS_size = CDS_length - int(gap_size)
    ORF_flag = ORF_check(CDS_size)
    return ORF_flag,Full_length,CDS_length,Full_size,CDS_size

def MNM(start,M1,N,M2):
    """
    INPUT: Left-most mapping position for the read or 'start', length of 'match' before the gap,
    the size or length of the gap (from alignment file), and the the length of the alignment 'match'
    after the gap. 'Matches' or M are anything outside of insertion and deletion, so the
    base can be the same or different from the reference to be considered a 'match'.

    OUTPUT: The gap location information
    """
    gap_start = start + M1
    gap_end = gap_start + (N-1)
    gap = '{0}_{1}'.format(gap_start,gap_end)
    return gap,gap_start,gap_end

def MNMNM(start,M1,N1,M2,N2,M3):
    """
    INPUT: Alignment information from cigar string.
    A 'match' or M can be the reference nucleotide or not.
    However, it won't be an 'N', deletion, or insertion

    OUTPUT: Where the gaps start and end
    """
    gap_start1 = start + M1
    gap_end1 = gap_start1 + (N1 -1)
    after_gap1 = gap_start1 + N1
    gap1 = '{0}_{1}'.format(gap_start1,gap_end1)
    gap_start2 = after_gap1 + M2
    gap_end2 = gap_start2 + (N2 -1)
    gap2 = '{0}_{1}'.format(gap_start2,gap_end2)
    return gap1,gap2,gap_start1,gap_end1,gap_start2,gap_end2


#how is name decided??
def printer1(outfile,name,segment,readname,gap,gap_start,gap_end,
    estimated_size,length,segment_length,ORF_flag):
    """
    INPUT: info
    OUTPUT: print out info to file
    """
    print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}'.format(name,segment,readname,
    gap,gap_start,gap_end,estimated_size,length,segment_length,ORF_flag), end="\n",file=outfile)

def printer2(outfile,name,segment,readname,gap1,gap2,gap_start1,gap_end1,
        gap_start2,gap_end2,align_start,cigar):
    """
    INPUT: info for two gaps
    OUTPUT: print info for two gaps to a specific file
    """
    print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}'.format(name,segment,readname,gap1,gap2,
    gap_start1,gap_end1,gap_start2,gap_end2,align_start,cigar), end="\n",file=outfile)

Full_seg, CDS_seg = GrabLengths(args.strain)
HEADER1 = 'name,segment,readname,gap,gap_start,gap_end,gap_size,estimated_length,segment_size,ORF_flag'
HEADER2 = 'name,segment,readname,gap1,gap2,gap_start1,gap_end1,after_gap1,gap_start2,gap_end2,align_start,cigar'
path = args.file_direct

#three outfiles generated
outfile = open("CandidateDI_OneGap_FullLength.csv",'w')
print(HEADER1, end ="\n",file=outfile)

outfile3 = open("CandidateDI_OneGap_CDS.csv",'w')
print(HEADER1, end ="\n",file=outfile3)

outfile2 = open("CandidateDI_TwoGap.csv",'w')
print(HEADER2, end ="\n",file=outfile2)

cigar_dict = {}
for infile in glob.glob( os.path.join(path, '*.split.txt' ) ):#will go through each specified csv file
    print(infile)
    #name = infile.split('/')[-1].split('.')[0].split('_')[0] #split file if needed
    name = infile.split('/')[-1].split('.')[0] #split file if needed
    df = pd.read_csv(infile, sep = '\t', keep_default_na = False)
    df.columns= ['readname','segment','left_pos','cigar']

    for index, row in df.iterrows():
        cigar = row['cigar']
        start = row['left_pos']
        segment = row['segment']
        readname = row['readname']

        break_cigar = re.findall(r'(\d+)(\w)',cigar) #break cigar into tuples
        #print(break_cigar)

        #count number of deletions and insertions
        indel = [int(x[0]) for x in break_cigar if x[1] == 'I' or x[1] == 'D']
        total_indel = sum(indel)
        if total_indel > args.total_mismatch:
            continue

        elif total_indel <= args.total_mismatch:
            no_soft = [i for i in break_cigar if i[1] != 'S' and i[1] !='I'] #doesn't include soft clips or insertions in alignment
            tup = all(int(i[0]) >= args.align_length for i in no_soft if i[1] =='M' or i[1]=='N') #make sure gap and alignment is greater than aligned length output true or false

            if tup == bool('True'): #if match/mismatch is greater than length specified
                Cigar_len = len(no_soft)
                Where_N = [x for x, y in enumerate(no_soft) if y[1] == 'N']
                #print(Where_N)

                if len(Where_N) == 1: #if only one gap was called
                    N_idx = Where_N[0] #find the position where the gap is
                    N = int(no_soft[N_idx][0]) #the length of the gap

                    #incase there are deletions or insertions- identify everything before and after 'N'
                    Before_idx = range(0,N_idx) #find the index positions of everything before the gap in cigar
                    After_idx = range(N_idx+1,Cigar_len) #find index positions of everything after gap includes deletions

                    five_length = [] #how many indx positions are there before the 'N'
                    for x in Before_idx:
                        five_length.append(int(no_soft[x][0]))
                    five_length=(sum(five_length)) #sum up all of the positions with a match,deletion,insertion before 'n'

                    three_length=[] #how many idx positions are there after 'N'
                    for x in After_idx:
                        three_length.append(int(no_soft[x][0]))
                    three_length = sum(three_length) #sum up all of the matches,deletions, insertions after 'N'

                    M1 = five_length #new 'match' before N
                    M2 = three_length #new "match" length after 'N'

                    #calculate the following useing the function MNM for the pattern
                    gap,gap_start,gap_end = MNM(start,M1,N,M2)

                    if N > args.gap_size:
                        Conserved_Five = Full_seg[segment] - CDS_seg[segment]
                        gap_start_FULL = gap_start + Conserved_Five
                        gap_end_FULL = gap_end + Conserved_Five
                        gap_FULL = '{0}_{1}'.format(gap_start_FULL,gap_end_FULL)
                        ORF_flag,FULL_length,CDS_length,FULL_size,CDS_size = EstSize(segment,N,args.strain)

                        printer1(outfile3,name,segment,readname,gap,gap_start,gap_end,N,CDS_size,CDS_length,ORF_flag)
                        printer1(outfile,name,segment,readname,gap_FULL,gap_start_FULL,gap_end_FULL,N,FULL_size,FULL_length,ORF_flag)


                if len(Where_N) == 2: #if there are two gaps called
                    N_idx1 = Where_N[0] #find the first gap in the pattern
                    N1 = int(no_soft[N_idx1][0]) #N1 is the length of first gap
                    N_idx2 = Where_N[1] #find the second gap in the pattern
                    N2 = int(no_soft[N_idx2][0]) #length of second gap
                    match1 = range(0,N_idx1) #'match positions before N1
                    match2 = range(N_idx1 + 1, N_idx2) #match positions between N1 and
                    match3 = range(N_idx2 + 1, Cigar_len)

                    match_1 = []
                    for x in match1:
                        match_1.append(int(no_soft[x][0]))
                    match_1=(sum(match_1))

                    match_2 =[]
                    for x in match2:
                        match_2.append(int(no_soft[x][0]))
                    match_2 = sum(match_2)

                    match_3 = []
                    for x in match3:
                        match_3.append(int(no_soft[x][0]))
                    match_3 = sum(match_3)

                    gap1,gap2,gap_start1,gap_end1,gap_start2,gap_end2, = MNMNM(start,
                    match_1,N1,match_2,N2,match_3)

                    if N1 > args.gap_size and N2 > args.gap_size:
                        printer2(outfile2,name,segment,readname,gap1,gap2,
                        gap_start1,gap_end1,gap_start2,gap_end2,
                        start,cigar)

                elif len(Where_N) >= 3:
                    continue

outfile.close()
outfile2.close()
outfile3.close()
