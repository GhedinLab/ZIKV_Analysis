"""
Updated 09/30/2019 KJ
"""

import os
import glob
import pandas as pd
import numpy as np

aminoacid = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}


def read_fasta(fp):
    """
    INPUT: Fasta file reference
    """
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def open_fasta(filename):
    """
    INPUT: fasta file reference
    OUTPUT: dictionary with segment name as key and sequence of segment as value
    """
    segdict = {}
    with open(filename) as fp:
        for name, seq in read_fasta(fp):
            segdict[name[1:]] = seq
    return segdict

def printerMajor(outfile,name,segment,ntpos,nt,majmin,freq,aa,codon,binocheck,totalcount,aapos):
    """
    INPUT: Major nucleotide information
    OUTPUT: Major information written to specified file
    """
    freq = str(freq)
    print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}'.format(name,segment,ntpos,nt,majmin,
    freq,aa,codon,'',binocheck,totalcount,aapos),end="\n",file=outfile)

def getindex(codon):
    """
    Identifying where the nucleotide falls in codon
    INPUT: major codon from snplist file
    OUTPUT: Index location of nucleotide
    """
    return [i for i, c in enumerate(codon) if c.isupper()]

def printerMinor(outfile,name,segment,ntpos,nt,majmin,freq,aa,codon,binocheck,totalcount,aapos):
    if(codon): #if codon exists
        if '-' in codon: #if there is a deletion
            #don't print any information for codon or amino acid
            print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}'.format(name,segment,
            ntpos,nt,majmin,freq,'','','',binocheck,totalcount,aapos),end = "\n",file=outfile)


        elif 'n' in codon:
            #don't print any information for codon or amino acid
            print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}'.format(name,segment,
            ntpos,nt,majmin,freq,'','','',binocheck,totalcount,aapos),end="\n",file=outfile)


        else:
            upperIndex = getindex(codon)[0] #identify where nucleotide is
            minorcodon = list(codon) #make codon a list of three nucleotides
            minorcodon[upperIndex] = nt.upper() #put minor nt in index position and make uppercase
            minorcodon = ''.join(minorcodon) #join into string (no longer list)


            if '-' in minorcodon: #if there is a deletion for minor nt
                #print without codon or amino acid
                print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}'.format(name,segment,
                ntpos,nt,majmin,freq,'','','',binocheck,totalcount,aapos),end="\n",file=outfile)


            else:
                #grab amino acid infromation
                minoraa = aminoacid[minorcodon.lower()]
                freq = str(freq)
                if minoraa == aa:
                    #if minor aa is the same as the major amino acid
                    nonsyn = 'syn'
                elif minoraa != aa:
                    nonsyn = 'nonsyn'
                print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}'.format(name,segment,ntpos,nt,majmin,freq,
                                minoraa,minorcodon,nonsyn,binocheck,totalcount,aapos),
                                end="\n",file=outfile)

    else: #if there isn't a codon (aka not divisible by 3) print the following without amino acid info
        print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}'.format(name,segment,ntpos,nt,
        majmin,freq,'','','',binocheck,totalcount,aapos),end="\n",file=outfile)




################################################################

STRAIN = 'ZIKA' #strain using, h3n2, h1n1, zika, etc.
path = '../../Variants/AllSnplist/' #where the snplist files are
savedirect = "../../Variants/Linegraphs/" #where save file
ref_file = open_fasta('../../reference/zika/MR766.ZIKA.CDS.fas')
HEADER = "sample,segment,ntpos,nt,majmin,freq,aa,codon,nonsyn,binocheck,totalcount,aapos"
MINORFREQ_CUTOFF = 0.03 #minimum frequency cutoff
COVERAGE_CUTOFF = 200 #need this coverage for minor variant
MINCOVERAGE = 200 #need this coverage for 'major' nucleotide
NT_list = ['A','T','G','C'] #doesn't include deletions
samptype = 'ZIKA'


######################################################################################################
outfile3 = open("{0}{1}.all.listVariants.{2}.csv".format(savedirect,STRAIN,MINORFREQ_CUTOFF),'w')
print(HEADER,end="\n",file=outfile3)
#print>>outfile3,HEADER

for SEGMENT in ref_file:
    print(SEGMENT)
    all_list = [] #collecting all of the minor ntpos in the collective group of csvs
    for infile in glob.glob( os.path.join(path, '*{0}.{1}.0.01.snplist.csv'.format(STRAIN,SEGMENT)) ):
        df = pd.read_csv(infile) #pandas will read the csv file create dataframe

        dfmin = df[(df.minor.isin(NT_list)) & (df.minorfreq >= MINORFREQ_CUTOFF) &
                    (df.totalcount >= COVERAGE_CUTOFF) & (df.binocheck ==bool('True'))]
                    #filter for minor variants
        ntpos_list = list(dfmin.ntpos)
        all_list = all_list + ntpos_list #concatenate lists
    all_list = list(set(all_list)) #remove duplicate positions

    outfile = open("{0}{1}.nobino.all.linegraph.{2}.{3}.{4}.{5}.csv".format(savedirect,SEGMENT,
    samptype,STRAIN,COVERAGE_CUTOFF,MINORFREQ_CUTOFF),'w')
    print(HEADER,end="\n",file=outfile)

    outfile2 = open("{0}{1}.all.listVariants.{2}.{3}.{4}.{5}.csv".format(savedirect,SEGMENT,
    samptype,STRAIN,COVERAGE_CUTOFF,MINORFREQ_CUTOFF),'w')
    print(HEADER,end="\n",file=outfile2)


    for infile in glob.glob( os.path.join(path, '*{0}.0.01.snplist.csv'.format(SEGMENT)) ):#will go through each specified csv file
        df = pd.read_csv(infile, keep_default_na=False)
        dfgraball = df[(df.ntpos.isin(all_list)) & (
                        df.totalcount>=MINCOVERAGE) &
                        df.major.isin(NT_list)]



        dfmin1 = dfgraball[(dfgraball.minor.isin(NT_list)) &
                            (dfgraball.totalcount >= COVERAGE_CUTOFF)]
        dfmin1['minorfreq'] = dfmin1['minorfreq'].astype(float)
        dfmin=dfmin1[dfmin1.minorfreq>=MINORFREQ_CUTOFF]

        for index, row in dfgraball.iterrows():
            printerMajor(outfile,row['name'],SEGMENT,row['ntpos'],row['major'],'major',
            row['majorfreq'],row['majoraa'],row['majorcodon'],
            row['binocheck'],row['totalcount'],row['aapos'])

        for index, row in dfmin.iterrows():
            print(row['minorfreq'])
            printerMinor(outfile,row['name'],SEGMENT,row['ntpos'],row['minor'],'minor',
            row['minorfreq'],row['majoraa'],row['majorcodon'],
            row['binocheck'],row['totalcount'],row['aapos'])

            printerMajor(outfile2,row['name'],SEGMENT,row['ntpos'], row['major'],
            'major',row['majorfreq'],row['majoraa'],row['majorcodon'],
            row['binocheck'],row['totalcount'],row['aapos'])

            printerMinor(outfile2,row['name'],SEGMENT,row['ntpos'], row['minor'],'minor',
            row['minorfreq'],row['majoraa'],row['majorcodon'],
            row['binocheck'],row['totalcount'],row['aapos'])

            printerMajor(outfile3,row['name'],SEGMENT,row['ntpos'], row['major'],'major',
            row['majorfreq'],row['majoraa'],row['majorcodon'],
            row['binocheck'],row['totalcount'],row['aapos'])

            printerMinor(outfile3,row['name'],SEGMENT,row['ntpos'], row['minor'],'minor',
            row['minorfreq'],row['majoraa'],row['majorcodon'],
            row['binocheck'],row['totalcount'],row['aapos'])


outfile.close()
outfile2.close()
outfile3.close()
