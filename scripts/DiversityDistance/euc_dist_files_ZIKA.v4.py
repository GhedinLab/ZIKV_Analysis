"""
run: python3 euc_dist_files_ZIKA.v4.py
needs: snplist files, reference, metadata

Last update: 07/03/2019 must add specified region
"""

import os
import glob
import pandas as pd
import numpy as np
import math
import scipy


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

def CheckFreq(NtFreq,MINORFREQ_CUTOFF):
    """INPUT: Frequency of nucleotide at that position, minor freq cutoff
    OUTPUT: Returns the frequency of that nucleotide at the ntpos to be added to a dict
    """
    if NtFreq >= MINORFREQ_CUTOFF:
        return NtFreq
    elif NtFreq < MINORFREQ_CUTOFF:
        return 0

def printer(outfile,name,segment,ntpos,A,C,G,T,mouse,transmission,tissue,sex,ID):
    """
    INPUT: info
    OUTPUT: print out info to file
    """
    print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}'.format(name,segment,ntpos,A,C,G,T,mouse,transmission,tissue,sex,ID),end="\n",file=outfile)

################################################################################
#The following may need to be changed depending on thresholds you are using
STRAIN = 'ZIKA'
path = '../../Variants/AllSnplist/'
refdict = open_fasta('../../reference/zika/MR766.ZIKA.CDS.fas')
mdata = '../../MetadataZika/Zika_Metadata_v4.csv'
MINORFREQ_CUTOFF = 0.03
freq = 3
COVERAGE_CUTOFF = 200
MINCOVERAGE_CUTOFF = 50
pcr = 'pcr1'
region = range(0,3883) #pcr1 region last full codon
#included 0 so we can tell if the alignment was off


################################################################################
HEADER = "sample,segment,ntpos,A,C,G,T,mouse,transmission,tissue,sex,ID"

#The following is specific to the metadata that you are working with!
dfm = pd.read_csv(mdata, keep_default_na = False)
namelist=[]
mouse_dict = {}
transmission_dict = {}
tissue_dict = {}
sex_dict = {}
id_dict = {}
for index, row in dfm.iterrows():
    sampname= row['name']
    namelist.append(sampname)
    mouse_dict[sampname] = row['mouse_id']
    transmission_dict[sampname] = row['transmission']
    tissue_dict[sampname] = row['tissue']
    sex_dict[sampname] = row['sex']
    id_dict[sampname] = row['id']
print(namelist)

###############################################################################

nucleotides = ['A', 'T', 'C', 'G']
gotthezeros = []
N_samples = set()
for SEGMENT in refdict:
    print(SEGMENT)
    #generate distance file to pipe into R
    outfile = open("../../Distance/PCR1/{0}.ntfreq.euc.{1}.{2}.{3}.{4}.csv".format(SEGMENT,STRAIN,MINCOVERAGE_CUTOFF,freq,pcr),'w')
    print(HEADER, end ="\n",file=outfile)
    #iterate through individual snplist files
    for infile in glob.glob( os.path.join(path, '*snplist.csv')):#will go through each specified csv file
        df = pd.read_csv(infile) #pandas will read the csv file create dataframe
        df = df[df['ntpos'].isin(region)]
        #print(df['ntpos'])  print just to make sure region pulled is accurate
        for index, row in df.iterrows(): #iterate through rows of snp files

            if row['ntpos'] == 0: #if there is low coverage it will add a position 0 in the beginning
                gotthezeros.append('%s_%s' %(row['name'],row['segment']))

            if row['name'] in namelist: #if the name is in our namelist (double check)
                if int(row['totalcount']) >= COVERAGE_CUTOFF:
                    #calculate the frequency of nucleotides at each position
                    A_freq = (float(row['A'])/(float(row['totalcount'])))
                    C_freq = (float(row['C'])/(float(row['totalcount'])))
                    G_freq = (float(row['G'])/(float(row['totalcount'])))
                    T_freq = (float(row['T'])/(float(row['totalcount'])))

                    A_check =CheckFreq(A_freq,MINORFREQ_CUTOFF) #check to make sure it passes our cutoff
                    C_check = CheckFreq(C_freq,MINORFREQ_CUTOFF)
                    G_check = CheckFreq(G_freq,MINORFREQ_CUTOFF)
                    T_check = CheckFreq(T_freq,MINORFREQ_CUTOFF)

                    Nt_dict = {'A':A_check,'C':C_check,'G':G_check,'T':T_check}


                #if the totalcount or coverage at the position is too low just add major freq
                #is this biasing things???
                elif int(row['totalcount']) < COVERAGE_CUTOFF and int(row['totalcount']) >= MINCOVERAGE_CUTOFF:
                    if row['major']== 'T':
                        Nt_dict = {'A':0,'C':0,'G':0,'T':1}
                    elif row['major']== 'A':
                        Nt_dict = {'A':1,'C':0,'G':0,'T':0}
                    elif row['major']== 'G':
                        Nt_dict = {'A':0,'C':0,'G':1,'T':0}
                    elif row['major']== 'C':
                        Nt_dict =  {'A':0,'C':1,'G':0,'T':0}
                    elif row['major'] =='-' or row['major'] =='N': #need to change so it just removes this nt position all together
                        Nt_dict = {'A':"N",'C':"N",'G':"N",'T':"N"}


                elif int(row['totalcount']) < MINCOVERAGE_CUTOFF:
                    Nt_dict = {'A':"N",'C':"N",'G':"N",'T':"N"}
                    N_samples.add(row['name'])

                #print(Nt_dict)
                printer(outfile,row['name'],row['segment'],row['ntpos'],Nt_dict['A'], Nt_dict['C'],Nt_dict['G'],Nt_dict['T'],
                mouse_dict[row['name']], transmission_dict[row['name']],tissue_dict[row['name']],sex_dict[row['name']],id_dict[row['name']])
print(gotthezeros)
print(N_samples)
