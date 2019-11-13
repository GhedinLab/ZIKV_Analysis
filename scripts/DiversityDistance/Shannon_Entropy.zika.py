#assumes that any position without a minor doesn't have any minor (aka will be 0 because freq is 1)
#uses log base 2 because we are only comparing major and minor
#will normalize to the size of the segment
#input minor and major frequency data
#output:
#1. Shannon by segment
#2. Shannon by entire genome
"""
Output a file that contains site specific calculations
"""


import os
import glob
import pandas as pd
import numpy as np
import math
import scipy

savedirect='../../Diversity/PCR1/'

region = range(1,3883) #PCR1
genlength = len(range(1,3883)) #only pcr1 product - ends with that last full codon doesn't include last number
#region = range(1,10261)
#genlength = len(range(1,10261)) #doesnt include last number
length_dict = {"MR766":genlength}
pcr = 'pcr1'

total_genome= 10260

MINCOVERAGE_CUTOFF= 200
MINORFREQ_CUTOFF = 0.03


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

def printer(outfile,sample, segment, segment_sum, idx_normalized,\
                tissue,transmission,sex,type_sample,mouse):
    print>>outfile, '{0},{1},{2},{3},{4},{5},{6},{7},{8}'.format(sample, segment, segment_sum, idx_normalized,
                    tissue,transmission,sex,type_sample,mouse)

def PositionPrinter(outfile,sample,ntpos,segment,majorfreq,minorfreq,shannon,\
                tissue,transmission,sex,type_sample,mouse):
    print>>outfile, '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}'.format(sample,ntpos,segment,
                    majorfreq,minorfreq,shannon,
                    tissue,transmission,sex,type_sample,mouse)

def shannon(freq):
    return(-(float(freq))*(math.log((float(freq)), 2))) #calculate the shannon index then sum over all variants

def shannon_mean(dataframe,SEGMENT,sample,outfile,outfile2): #08/13/2019
    #print(SEGMENT)
    min_list = []
    maj_list = []
    site_list = []


    for index, row in dataframe.iterrows():
        maj_list.append(shannon(row['majorfreq']))
        min_list.append(shannon(row['minorfreq']))
        #print to a file individual ntpos
        sdi_ntpos = (shannon(row['majorfreq']) + shannon(row['minorfreq']))

        PositionPrinter(outfile2,sample,row['ntpos'],SEGMENT,row['majorfreq'],
        row['minorfreq'],sdi_ntpos,
        tissue_dict[sample], transmission_dict[sample],
        sex_dict[sample], type_dict[sample],mouse_dict[sample])

    site_list = [x + y for x, y in zip (min_list, maj_list)] #sum to get sdi at ntpos
    segment_sum = (sum(site_list)) #sum across segment, no minor is 0

    #idx_normalized = segment_sum/length_dict[SEGMENT] #normalized by segment length


    idx_normalized = float(segment_sum)/(float(length_dict[SEGMENT])/float(1000)) #normalized by segment length in kb

    printer(outfile,sample,SEGMENT,segment_sum, idx_normalized,
            tissue_dict[sample], transmission_dict[sample],
            sex_dict[sample], type_dict[sample],mouse_dict[sample])


    return(segment_sum,idx_normalized)

mdata = '../../MetadataZika/Zika_Metadata_v4.csv'
dfm = pd.read_csv(mdata)
zikadict = open_fasta('../../reference/zika/MR766.ZIKA.CDS.fas')

tissue_dict={}
transmission_dict = {}
sex_dict={}
type_dict={}
mouse_dict= {}
id_dict={}
for index, row in dfm.iterrows():
    tissue_dict[row['name']] = row['tissue']
    transmission_dict[row['name']] = row['transmission']
    sex_dict[row['name']] = row['sex']
    type_dict[row['name']] = row['type']
    mouse_dict[row['name']] = row['mouse_id']
    id_dict[row['name']] = row['id']

nucleotides = ['A', 'G', 'T', 'C'] #don't want deletions but should use


outfile = open("{2}ShannonDiversity_STAR.zika.perkb.{0}.{1}.{3}.csv".format(MINCOVERAGE_CUTOFF,MINORFREQ_CUTOFF,savedirect,pcr),'w')
print>>outfile,"sample,segment,shannon_sum,normalized_perkb,tissue,transmission,sex,type,mouse"

outfile2 = open("{2}ShannonDiversity_STAR.zika.ntpos.{0}.{1}.{3}.csv".format(MINCOVERAGE_CUTOFF,MINORFREQ_CUTOFF,savedirect,pcr),'w')
print>>outfile2,"sample,ntpos,segment,majorfreq,minorfreq,shannon,\
                tissue,transmission,sex,type_sample,mouse"

filename='../../Variants/Linegraphs/ZikaVariants.{0}.{1}.csv'.format(MINCOVERAGE_CUTOFF,MINORFREQ_CUTOFF)
df = pd.read_csv(filename,keep_default_na=False)
names = list(set(df['sample']))

for name in names:
    #print(name)
    segment_sums=[]
    segment_idxs = []
    df_names = df[df['sample'] == name]
    for SEGMENT in zikadict:
        print(SEGMENT)
        df_seg = df_names[df_names['segment'] ==SEGMENT]
        df_seg = df_seg[(df_seg['ntpos'].isin(region)) &
        (df_seg['binocheck']== bool("TRUE"))]
        #df_seg = df_names[df_names['segment'] ==SEGMENT]
        seg_sum,segment_idx = shannon_mean(df_seg,SEGMENT,name,outfile,outfile2)
        segment_idxs.append(segment_idx)
        segment_sums.append(seg_sum)

    segments_summed = sum(segment_sums)
    normalized_segments = segments_summed/(total_genome/1000)

    printer(outfile, name, 'FULL', segments_summed, normalized_segments,
            tissue_dict[name], transmission_dict[name],
            sex_dict[name], type_dict[name],mouse_dict[name])
