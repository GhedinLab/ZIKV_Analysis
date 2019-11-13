#create consensus fasta
"""
python ConsensusFasta.py --ref reference/vic/victoria.flub.fasta --var FILES/ --strain vic
python ConsensusFasta.py --ref reference/yam/yamagata.flub.fasta --var FILES/yam --strain yam
"""
import os
import glob
import pandas as pd
import numpy as np
import math
import scipy
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--ref','-r',required=True,help='Indicate path and reference fasta file') #args.ref
parser.add_argument('--var','-v',required=True,help='Indicate path to variant files') #args.ref
parser.add_argument('--strain','-s',required=True,help='Indicate strain') #args.ref

args = parser.parse_args()

def read_fasta(fp): #function to read each segment of fasta, from tim readreport
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def open_fasta(filename): #opens each fasta- makes a dictionary of segment name from tim readreport
    # filename = '../FILES/reference/'+someref
    segdict = {}
    with open(filename) as fp:
        for name, seq in read_fasta(fp):
            segdict[name[1:]] = seq
    return segdict


#########################################################################################################################

ref = args.ref
path = args.var
STRAIN = args.strain
MINORFREQ_CUTOFF = 0.03 #the minor freq cutoff
freq = 3 #percentage for the minor freq cutoff
COVERAGE_CUTOFF = 200 #the lowest coverage that you will allow - this is usually acceptable for the major nucleotide
refdict = open_fasta(ref)


for SEGMENT in refdict:
    print SEGMENT
    outf = open('{0}.{1}.{2}.fasta'.format(SEGMENT,COVERAGE_CUTOFF,STRAIN), 'w')
    for infile in glob.glob( os.path.join(path, '*{0}.{1}.0.01.snplist.csv'.format(STRAIN.upper(),SEGMENT)) ):#will go through each specified csv file
        consensus_seq = []
        df = pd.read_csv(infile) #pandas will read the csv file create dataframe
        name = list(set(df['name']))[0]
        fullname = '{0} {1}'.format(name, SEGMENT)

        print(fullname)
        for index, row in df.iterrows():
            if row['totalcount'] >= COVERAGE_CUTOFF:
                consensus_seq.append(row['major'])

            elif row['totalcount'] < COVERAGE_CUTOFF:
                consensus_seq.append('N')

        consensus = ''.join(consensus_seq)
        outf.write('>' + fullname + '\n' + consensus + '\n')
        print(consensus)
    outf.close()
