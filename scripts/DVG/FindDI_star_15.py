import os
import argparse
import pandas as pd
from scipy.stats.distributions import binom
import re

"""
Notes: Bam files cannot have '.' within the names
Notes: Needs to be ran with full pipeline ()
LAST EDIT: 07/10/2020
For SMART1 pacbio data (h3n2 did not have segments
named 'ha' but rather just 'ha')

To run:
python3 FindDI_star_4.py --strain CAL09 --file_direct N_files

NOTE:
Before running this script- make separate 'N' files:
module load samtools/intel/1.6
for i in $(ls ./bamfiles/rmdup_bams); do if [[ $i == *rmd.star.bam ]];
then samtools view ./bamfiles/rmdup_bams/$i | awk '$6 ~ /N/' | cut -f 1,3,4,6 |
uniq > $i.split.txt;fi;done
And move all 'N files to their own dierectory'

Only considering two possibilities:
MNM (one gap) and MNMNM (two gaps)
Insertions are not considered when finding gap start and stop positions

Mismatches/Deletions/Insertions are only accepted to the total mismatch allowed
(--total-mismatch).

This will quantify the number of mismatch before and after the gap - if it is
greater than what we allow- then it isn't considered in our candidate list.


Suggestions:
Alignments should be done to the coding sequence for the reference file
If there are samples aligned to multiple references- it is best to run it for
each strain aligned to (such as H1N1 and H3N2)


INPUT: 'N' files from split read alignment
OUTPUT:
1 file with coordinates using the CDS
1 file with coordinates using the full length genome
1 file where two gaps were called


Required Input files:
Split reads csv file
Idxstats file

To run:
python3 FindDI_star_11.py --strain H3N2 --ref ../../reference/h3n2/A_NewYork_238_2005_H3N2_refseq.FULL.fasta --file ../../h3n2.SplitReads.csv --idx ../../h3n2.IdxStats.csv
python3 FindDI_star_13.py --strain H3N2 --ref ../../reference/h3n2/A_NewYork_238_2005_H3N2_refseq.FULL.fasta --file ../h3n2.SplitReads.csv --idx ../H3N2.FeatureCounts.csv


python3 FindDI_star_13.py --strain H3N2 --ref ../../reference/h3n2/A_NewYork_238_2005_H3N2_refseq.FULL.fasta --file ../h3n2.SplitReads.csv --idx ../H3N2.FeatureCounts.csv
python3 FindDI_star_13.py --strain PR8 --ref ../../reference/PR8/Influenza_A_H1N1_PR8_refseq.FULL.fasta --file ../../pr8_control/pr8.SplitReads.csv --idx ../../pr8_control/pr8.FeaturesOutput.csv --save_dir ../../pr8_control

python3 FindDI_star_13.py --strain H3N2 --ref ../../reference/h3n2/A_NewYork_238_2005_H3N2_refseq.FULL.fasta --file ../../simulations/h3n2_sim/h3n2.SplitReads.csv --idx ../../simulations/h3n2_sim/h3n2.FeaturesOutput.csv --save_dir ../../simulations/h3n2_sim
python3 FindDI_star_13.py --strain H3N2 --ref ../../reference/h3n2/A_NewYork_238_2005_H3N2_refseq.FULL.fasta --file ../../simulations/h3n2_sim/h3n2.SplitReads.csv --idx ../../simulations/h3n2_sim/h3n2.FeaturesOutput.csv --save_dir ../../simulations/h3n2_sim --align_length 10

python3 FindDI_star_14.py --strain COV19 --ref /home/kate/Lab/COV19/MergedRuns/reference/SARS-COV-2.fasta --file /home/kate/Lab/COV19/MergedRuns/DVG/cov19.SplitReads.csv --idx /home/kate/Lab/COV19/MergedRuns/DVG/cov19.FeaturesOutput.csv --save_dir /home/kate/Lab/COV19/MergedRuns/DVG --align_length 25 --total_mismatch 1 --bam_path ./merged_bams/

Outputs several files in case they want to be checked later on:

python3 FindDI_star_15.py --strain 1 --bam_path ../SortedSubBams/ --ref ../../reference/zika/MR766.ZIKA.CDS.fas --file ../../SubVariants/SubDVG/1.SplitReads.csv --idx ../../SubVariants/SubDVG/1.FeaturesOutput.csv --save_dir ../../SubVariants/SubDVG --align_length 25 --total_mismatch  1


Should include an arguement indicating which sample is the 'control', this
can then use the highest percent abundance frequency as its threshold for the binocheck

"""

parser = argparse.ArgumentParser()
parser.add_argument('--strain','-s',required=True,help='Indicate strain- necessary for size calculations can only be CAL09,PR8,H3N2')
parser.add_argument('--ref','-r',required=True,help='Give directory where reference is located')
parser.add_argument('--file','-f',required=True,help='split txt file')
parser.add_argument('--idx','-i',required=True,help='split txt file')
#default variables below this, but can be input by user
parser.add_argument('--save_dir','-sd',default='.',help='indicate save directory, dont use / at end')
parser.add_argument('--bam_path','-bp',default='./bamfiles/rmdup_bams/',help='Indicate bamfile path (beginning of feature counts)')
parser.add_argument('--align_length','-l',type=int,default=10,help='Give the length that each portion of the read must be (default 10 bp)')
parser.add_argument('--gap_size','-g',type=int,default=1,help='Give the minimum size of the gap (default 1)')
parser.add_argument('--total_mismatch','-m',type=int,default=50,help='Only 3 mismatches allowed in a read called as a DVG this includes mismatches/insertions/deletions')
args = parser.parse_args()

############################### FUNCTIONS ##################################
def ORF_check(CDS_estimated_size):
    """
    INPUT: the estimated size using the CDS coordinates
    OUTPUT: Whether the estimated size is divisible by 3 and therefore maintaining open read frame
    """
    if CDS_estimated_size%3 == 0: #if the remainder is zero - ORF maintained
        ORF_flag = 'T'
    elif CDS_estimated_size%3 != 0: #elif the remainder is not zero
        ORF_flag = "F"
    return ORF_flag

def EstSize(segment_length,gap_size):#change if CAL09, PR8, H3N2
    """
    INPUT: Segment, the measured gap size (from N_file), strain
    OUTPUT: The estimated size of the DVG using both the CDS and full segment coordinates
    """
    CDS_size = int(segment_length) - int(gap_size)
    ORF_flag = ORF_check(CDS_size)
    return ORF_flag,CDS_size

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

def printer1(outfile,name,segment,readname,gap,gap_start,gap_end,
    estimated_size,length,segment_length,ORF_flag,totalreads,
    readflags,strain):
    """
    INPUT: info
    OUTPUT: print out info to file
    """
    print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12}'.format(name,segment,readname,
    gap,gap_start,gap_end,estimated_size,length,segment_length,
    ORF_flag,totalreads,readflags,strain), end="\n",file=outfile)

def printer2(outfile,name,segment,readname,gap1,gap2,gap_start1,gap_end1,
        gap_start2,gap_end2,gap_size1,gap_size2,align_start,cigar,totalreads,
        readflags,strain,segment_length):
    """
    INPUT: info for two gaps
    OUTPUT: print info for two gaps to a specific file
    """

    print('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16}'.format(name,segment,
    readname,gap1,gap2,gap_start1,gap_end1,gap_start2,gap_end2,gap_size1,gap_size2,align_start,cigar,
    totalreads,readflags,strain,segment_length), end="\n",file=outfile)

def PrepFreq(filename, GapNumber): #df
    """
    INPUT: The candidate gap dataframe file ({0}/{1}.CandidateDI_OneGap.N{2}.Mis{3}.csv)
    OUTPUT: Outputs dataframe to be written as a csv file used for the grouping scripts CW wrote
    """
    if GapNumber == 1:
        # columns to be used to group and merge by
        groups = ['segment','gap_start','gap_end','gap','gap_size','estimated_length','ORF_flag']

        df = pd.read_csv(filename, sep = ',', keep_default_na = False) # read file

        # count number of reads that match particular dvg
        dvg_freq = df.groupby(groups).size().reset_index(name="freq")

        # drop read information to count number of samples with dvg
        df2 = df.drop(['readname','totalreads','readflags'],axis=1)

        # drop any duplicates to not double count
        df2 = df2.drop_duplicates()

        # count number of samples with dvg type
        samp_freq = df2.groupby(groups).size().reset_index(name='number_of_samples')

        # merge the two dataframes
        df_merge = pd.merge(dvg_freq,samp_freq,how='left',left_on=groups,right_on=groups)

        return df_merge # return dataframe to be used in grouping script

    if GapNumber == 2:
        # columns used to group and merge with
        groups = ['segment','gap_start1','gap_end1','gap1','gap_start2','gap_end2','gap2','gap_size1','gap_size2']

        # read file
        df = pd.read_csv(filename, sep = ',', keep_default_na = False)

        # count the number of reads that match particular dvg across samples
        dvg_freq = df.groupby(groups).size().reset_index(name="freq")

        # drop read information to count number of samples with given dvg
        df2 = df.drop(['readname','totalreads','readflags'],axis=1)

        # drop any duplicates
        df2 = df2.drop_duplicates()

        # count number of samples with DVG types
        samp_freq = df2.groupby(groups).size().reset_index(name='number_of_samples')

        # merge the two dataframes to be used in grouping script
        df_merge = pd.merge(dvg_freq,samp_freq,how='left',left_on=groups, right_on=groups)

        return df_merge

def MergeReadsGroups(ReadsFile, GroupedFile, GapNumber):
    """
    INPUT: "Reads file" which is the file that contains all read infromation.
    "Grouped File" which contains the grouping information. Gap number - whether
    we are looking at 1 gap or 2 gaps.

    OUTPUT: Outputs a dataframe that has read info along with grouped read info.
    To be used to calculate counts, rpkm, percentage, etc.
    """

    if GapNumber == 1:
        # read the candidate reads file
        df = pd.read_csv(ReadsFile, sep = ',', keep_default_na = False)

        # read the grouped/new gap files from grouping script
        d1g =  pd.read_csv(GroupedFile, sep = ',', keep_default_na = False)

        # merge the two dataframes
        d1g_m = pd.merge(df,d1g,how='left',on=['segment','gap','gap_start','gap_end','gap_size','estimated_length','ORF_flag'])

        # drop any duplicates that were generated by merge
        d1g_m = d1g_m.drop_duplicates()

        # rename the freq column to keep for later use
        d1g_m = d1g_m.rename(columns = {"freq":"FreqAcrossSamples"})

        # return to be used for rpkm, percentages, etc.
        return d1g_m

    if GapNumber == 2:
        # read the candidate two gap file
        df2 = pd.read_csv(ReadsFile, sep = ',', keep_default_na = False)

        # read the grouped New gap file made using grouping script
        d2g = pd.read_csv(GroupedFile, sep = ',', keep_default_na = False)

        # merge the two dataframes
        d2g_m = pd.merge(df2,d2g,how='left',on=['segment','gap1','gap_start1','gap_end1','gap2','gap_start2','gap_end2','gap_size1','gap_size2'])

        # drop any duplicates potentially generated from merging
        d2g_m = d2g_m.drop_duplicates()

        # rename freq column to be used later on and not get confused
        d2g_m = d2g_m.rename(columns = {"freq":"FreqAcrossSamples"})

        # return to be used for counting, rpkm, percentage, etc.
        return d2g_m

def CountGroupedDVGs(df, GapNumber):
    """
    INPUT: Take in grouping/read merged dataframe from MergedReadGroups, count the
    new gap information for each sample

    OUTPUT: new dataframe that has the freq info for total number of gaps and
    Count information for each individual dvg type within the samples. This will
    then be used as input into the binomial check. Counts calculated with
    "New" gaps generated by grouping script
    """
    if GapNumber == 1:
        # count num reads / dvg
        dvg_freq = df.groupby(['name','segment','NewStart','NewEnd','NewGap']).size().reset_index(name="DVG_freq")

        # count num reads /sample segment
        total_dvg = df.groupby(['name','segment']).size().reset_index(name="SegTotalDVG")

        # merge the dvg freq and total dvg counts into df
        dvg_c = pd.merge(dvg_freq,total_dvg,on=['name','segment'])

        # merge the dvg freq, total dvg counts with the full dataframe with read info
        dvg_f = pd.merge(df, dvg_c, on=['name','segment','NewStart','NewEnd','NewGap'])

        # drop any duplicates generated during merging
        dvg_f.drop_duplicates()

        # return for rpkm, percentage calcs
        return dvg_f

    if GapNumber == 2:
        # count num reads / dvg
        dvg_freq = df.groupby(['name','segment','NewStart1','NewEnd1','NewGap1','NewStart2','NewEnd2','NewGap2']).size().reset_index(name="DVG_freq")

        # count num reads/sample segment
        total_dvg = df.groupby(['name','segment']).size().reset_index(name="SegTotalDVG")

        # merge the dvg counts and total dvg counts
        dvg_c = pd.merge(dvg_freq,total_dvg,on=['name','segment'])

        # merge count information with all data which includes read info
        dvg_f = pd.merge(df, dvg_c, on=['name','segment','NewStart1','NewEnd1','NewGap1',"NewStart2","NewEnd2","NewGap2"])

        # drop any duplicates potentially generated during merge
        dvg_f.drop_duplicates()

        # return dataframe to be used for rpkm, binocheck, percentage, etc.
        return dvg_f

def binocheck(df, GapNumber, cutoff = 0.00125):
    """
    INPUT: Dataframe with count information for the 'new' gap info generated
    using the grouping script and CountGroupedDVG function above. Also includes
    whether the dataframe is a 1 gap df or 2 gap df. And a cutoff used to determine
    the significance.

    FUTURE: Collect multiple controls for DVG cutoff - compute distribution for
    e (error rate, or p) then use a T-test to determine whether the # of DVGs
    being identified through pipeline is stat sig compared to our e distribution.

    Can only make claims one way - which is confirming that there are DVGs within the
    sample. If trying to prove that there aren't any DVGs - you can't do that.
    The alpha would need to change and be less than 0.05 (but not sure what) -
    a distribution of e would be needed.

    cutoff of 0.0025 was generated using the full template pr8 with
    superscript 3. this can be changed based off of a

    # these flags indicate "Mapped within the insert size and correct orientation"
    # info from: https://www.samformat.info/sam-format-flag
    # AllowedFlags = [99, 147, 83, 163]
    """

    alpha = 0.05 # sig threshold
    ForwardFlags = [147,83] # Mate reverse strand first, mate reverse strand second
    ReverseFlags = [99, 163] # Read reverse strand second, read reverse strand first

    fdf = df[df.readflags.isin(ForwardFlags)] # subset forward data
    rdf = df[df.readflags.isin(ReverseFlags)] # subset reverse data

    if GapNumber == 1:
        # group by name, segment and dvg information to count number of rows which
        # corresponds to # forward or # reverse
        groups = ['name','segment','NewStart','NewEnd','NewGap','DVG_freq','SegTotalDVG','totalreads']
        s_fdf = fdf.groupby(groups).size().reset_index(name="Forward_freq")
        s_rdf = rdf.groupby(groups).size().reset_index(name="Reverse_freq")

        # merge the two dataframes and keep everything (outer)
        m_df = pd.merge(s_fdf, s_rdf, on = groups, how='outer')
        m_df = m_df.drop_duplicates()

        # after we merge fill any na information with a 0
        # this will indicate 0 forward or 0 reverse reads cover that specific DVG
        m_df['Forward_freq'] = m_df['Forward_freq'].fillna(0)
        m_df['Reverse_freq'] = m_df['Reverse_freq'].fillna(0)


        m_df['pforward'] = 1 - binom.cdf(m_df['Forward_freq'], m_df['totalreads'], cutoff)
        m_df['preverse'] = 1 - binom.cdf(m_df['Reverse_freq'], m_df['totalreads'], cutoff)

        #m_df.loc[(m_df['pforward'] <= alpha/2) & (m_df['preverse'] <= alpha/2), "binom"] = 'pass'
        m_df.loc[(m_df['pforward'] <= alpha) & (m_df['preverse'] <= alpha), "binom"] = 'pass'
        m_df.loc[(m_df['pforward'] <= alpha) & (m_df['preverse'] > alpha), "binom"] = 'reverse_fail'
        m_df.loc[(m_df['pforward'] > alpha) & (m_df['preverse'] <= alpha), "binom"] = 'forward_fail'
        m_df.loc[(m_df['pforward'] > alpha) & (m_df['preverse'] > alpha), "binom"] = 'both_fail'

        # merge with all read data to have
        f_df = pd.merge(df, m_df, how='left', on = groups)
        f_df = f_df.drop_duplicates()

        return f_df

    if GapNumber == 2:
        # group by name, segment and dvg information to count number of rows which
        # corresponds to # forward or # reverse
        groups = ['name','segment','NewStart1','NewEnd1','NewGap1','NewStart2','NewEnd2','NewGap2','DVG_freq','SegTotalDVG','totalreads']
        s_fdf = fdf.groupby(groups).size().reset_index(name="Forward_freq")
        s_rdf = rdf.groupby(groups).size().reset_index(name="Reverse_freq")

        # merge the two dataframes and keep everything (outer)
        m_df = pd.merge(s_fdf, s_rdf, on = groups, how='outer')
        m_df = m_df.drop_duplicates()

        # after we merge fill any na information with a 0
        # this will indicate 0 forward or 0 reverse reads cover that specific DVG
        m_df['Forward_freq'] = m_df['Forward_freq'].fillna(0)
        m_df['Reverse_freq'] = m_df['Reverse_freq'].fillna(0)

        # binom check
        m_df['pforward'] = 1 - binom.cdf(m_df['Forward_freq'], m_df['totalreads'], cutoff)
        m_df['preverse'] = 1 - binom.cdf(m_df['Reverse_freq'], m_df['totalreads'], cutoff)

        m_df.loc[(m_df['pforward'] <= alpha) & (m_df['preverse'] <= alpha), "binom"] = 'pass'
        m_df.loc[(m_df['pforward'] <= alpha) & (m_df['preverse'] > alpha), "binom"] = 'reverse_fail'
        m_df.loc[(m_df['pforward'] > alpha) & (m_df['preverse'] <= alpha), "binom"] = 'forward_fail'
        m_df.loc[(m_df['pforward'] > alpha) & (m_df['preverse'] > alpha), "binom"] = 'both_fail'

        # merge with the entire read data to have
        f_df = pd.merge(df, m_df, how='left', on = groups)
        f_df = f_df.drop_duplicates()

        return f_df

def rpkm_and_percent(df):
    """
    INPUT: Dataframe with read counts calculated using the new gap information
    from the grouping script.

    OUTPUT: DVG specific RPKM, total gapped read RPKM, DVG type percent abundance,
    Total gapped-read percent abundance, DVG type percentage out of total gapped reads.

    All calculated using segment specific information and read counts. All are
    normalized by the segment size.
    """
    #print(df.columns)
    df['DVG_RPKM'] = (df['DVG_freq'])/((df['segment_size']/1000) * (df['totalreads']/1000000))
    df['Segment_RPKM'] = (df['SegTotalDVG'])/((df['segment_size']/1000) * (df['totalreads']/1000000))

    df['DVG_PercentAbundancePerkb'] = ((df['DVG_freq'])/((df['segment_size']/1000) * (df['totalreads']))) * 100
    df['Segment_PercentAbundancePerkb'] = ((df['SegTotalDVG'])/((df['segment_size']/1000) * (df['totalreads']))) *100

    df['DVG_PercentTotalDVG'] = ((df['DVG_freq'])/((df['segment_size']/1000) * (df['SegTotalDVG']))) *100

    return df

def ReduceDF(df, GapNumber):
    """
    INPUT: Read dataframe
    OUTPUT: A reduced dataframe, with info that we care about for figures
    """
    if GapNumber == 1:
        dropcols = ['readname', 'gap', 'gap_start', 'gap_end',
       'gap_size', 'estimated_length','ORF_flag',
       'readflags','Forward_freq', 'Reverse_freq', 'pforward',
       'preverse','FreqAcrossSamples','number_of_samples']

        df = df.drop(dropcols, axis = 1)

        df = df.drop_duplicates()

        return df

    if GapNumber == 2:
        dropcols = ['readname','align_start','cigar','readflags',
        'Forward_freq','Reverse_freq','pforward','preverse','gap_start1','gap_end1',
        'gap1','gap_start2','gap_end2','gap2','gap_size1','gap_size2','FreqAcrossSamples','number_of_samples']

        df = df.drop(dropcols,axis=1)

        df = df.drop_duplicates()

        return df

def CleanFeatureCounts(df,bam_dir,strain):
    idcols = list(df.columns[0:6])
    sampcols = list(df.columns[6:])
    #print(df)
    fdf = pd.melt(df,id_vars=idcols, value_vars=sampcols,var_name="bam_name",value_name='totalreads')
    #data.composers.str.split('\s+').str[0]

    # add name, strain in bam must be lower case, bam dir default is ./bamfiles/rmdup_bams
    # make sure that bam directory used as input is similar
    fdf['name'] = fdf.bam_name.str.split(bam_dir).str[1].str.split('.{0}.'.format(strain.lower())).str[0]
    # clean up and change column names
    #print(fdf)
    fdf = fdf[['Chr','Length','name','totalreads']]

    #rename columns to fit in script
    fdf = fdf.rename(columns={'Chr':'segment','Length':'SegmentLength','name':'name','totalreads':'totalreads'})

    #print(fdf)

    return fdf


###############################################################################

# read in split read file
infile = args.file
fdf = pd.read_csv(infile, sep = ',', keep_default_na = False)
# read in index stat file
#print(fdf)

idxstats = args.idx
#print(idxstats)
idf = pd.read_csv(idxstats, sep = ',', keep_default_na = False)
idf = CleanFeatureCounts(idf, args.bam_path, args.strain)

# merge the two files to combine information (will use this to work on)
df = pd.merge(fdf,idf, how='left', left_on = ['segment','name'], right_on= ['segment','name'])

print(df)

# Open one gap file to write to
Reads1N = "{0}/{1}.CandidateDVG_OneGap.N{2}.Mis{3}.M{4}.csv".format(args.save_dir,
args.strain,args.gap_size,args.total_mismatch,args.align_length)
OneGapFile = open(Reads1N,'w')
HEADER1 = 'name,segment,readname,gap,gap_start,gap_end,gap_size,estimated_length,segment_size,ORF_flag,totalreads,readflags,strain'
print(HEADER1, end ="\n",file=OneGapFile)

# Open two gap file to write to
Reads2N = "{0}/{1}.CandidateDVG_TwoGap.N{2}.Mis{3}.M{4}.csv".format(args.save_dir,
args.strain,args.gap_size,args.total_mismatch,args.align_length)
TwoGapFile = open(Reads2N,'w')
HEADER2 = 'name,segment,readname,gap1,gap2,gap_start1,gap_end1,gap_start2,gap_end2,gap_size1,gap_size2,align_start,cigar,totalreads,readflags,strain,segment_size'
print(HEADER2, end ="\n",file=TwoGapFile)

# Iterate through rows to generate gap infromation for one and two gapped reads
one_count = 0
two_count = 0
for index, row in df.iterrows():
    # set different variables from row
    cigar = row['cigar']
    start = row['left_pos']
    segment = row['segment']
    readname = row['readname']
    seg_length = row['SegmentLength']
    name = row['name']
    totalreads = row['totalreads']
    #num_unmapped = row['num_unmapped']
    readflag = row['flags']

    break_cigar = re.findall(r'(\d+)(\w)',cigar) # break cigar into tuples

    # count number of deletions and insertions
    indel = [int(x[0]) for x in break_cigar if x[1] == 'I' or x[1] == 'D']
    total_indel = sum(indel) # sum total deletions and insertions

    if total_indel > args.total_mismatch:
        # if total deletions/insertions are greater than set cutoff ignore and continue
        # default set to 50
        continue

    elif total_indel <= args.total_mismatch:
        # else if the number is below cutoff calculate dvg information
        no_soft = [i for i in break_cigar if i[1] != 'S' and i[1] !='I'] # doesn't include soft clips or insertions in alignment
        tup = all(int(i[0]) >= args.align_length for i in no_soft if i[1] =='M' or i[1]=='N') # make sure gap and alignment is greater than aligned length output true or false
        Cigar_len = len(no_soft)

        # determine the number of 'N'
        Where_N = [x for x, y in enumerate(no_soft) if y[1] == 'N']

        if (tup == bool('True') and len(Where_N) == 1): #if match/mismatch is greater than length specified
            #Cigar_len = len(no_soft)
            # determine the number of 'N'
            #Where_N = [x for x, y in enumerate(no_soft) if y[1] == 'N']

            #if len(Where_N) == 1: # if only one gap was called
            N_idx = Where_N[0] # find the position where the gap is
            N = int(no_soft[N_idx][0]) #the length of the gap

                #incase there are deletions or insertions- identify everything before and after 'N'
            Before_idx = range(0,N_idx) #find the index positions of everything before the gap in cigar
            After_idx = range(N_idx+1,Cigar_len) #find index positions of everything after gap includes deletions

            five_length = [] # how many indx positions are there before the 'N'
            for x in Before_idx:
                five_length.append(int(no_soft[x][0]))
            five_length=(sum(five_length)) # sum up all of the positions with a match,deletion,insertion before 'n'


            three_length=[] # how many idx positions are there after 'N'
            for x in After_idx:
                three_length.append(int(no_soft[x][0]))
            three_length = sum(three_length) # sum up all of the matches,deletions, insertions after 'N'

            M1 = five_length #new 'match' before N
            M2 = three_length #new "match" length after 'N'

                # calculate the following useing the function MNM for the pattern
            gap,gap_start,gap_end = MNM(start,M1,N,M2)

            if N > args.gap_size:
                one_count += 1
                ORF_flag,est_length = EstSize(seg_length, N)
                printer1(OneGapFile,name,segment,readname,gap,gap_start,
                gap_end,N,est_length,seg_length,ORF_flag,totalreads,
                readflag,args.strain)


        if (tup == bool('True') and len(Where_N) == 2) :  # if there are two gaps called
            N_idx1 = Where_N[0] #  find the first gap in the pattern
            N1 = int(no_soft[N_idx1][0]) #  N1 is the length of first gap
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
                two_count += 1
                printer2(TwoGapFile,name,segment,readname,gap1,gap2,
                gap_start1,gap_end1,gap_start2,gap_end2, N1, N2,
                start,cigar,totalreads,
                readflag,args.strain,seg_length)

        elif len(Where_N) >= 3: # else ignore if greater than two gaps
            continue

TwoGapFile.close()
OneGapFile.close()
###############################################################################

#PREPPING THE FILES TO BE PUT THROUGH CHANGS GROUPING CODE
# ONE GAP FILE # Reads1N
if one_count > 0:
    # freq file to be used in grouping
    Freqs1N = "{0}/{1}.DVG.freqs.OneGap.N{2}.Mis{3}.M{4}.csv".format(args.save_dir,
    args.strain,args.gap_size,args.total_mismatch,args.align_length)

    # groups files which will be output from grouping
    Group1N = "{0}/{1}.DVG.grouped.OneGap.N{2}.Mis{3}.M{4}.csv".format(args.save_dir,
    args.strain,args.gap_size,args.total_mismatch,args.align_length)

    # final file that will have all information (excessive) for reads and sample
    All1N = "{0}/{1}.DVG.AllInfo.OneGap.N{2}.Mis{3}.M{4}.csv".format(args.save_dir,
    args.strain,args.gap_size,args.total_mismatch,args.align_length)

    # final file with information needed for figures only
    Final1N = "{0}/{1}.DVG.FINAL.OneGap.N{2}.Mis{3}.M{4}.csv".format(args.save_dir,
    args.strain,args.gap_size,args.total_mismatch,args.align_length)

    # generate freq files for OneGap DVGs- to be used for grouping
    df_merge = PrepFreq(Reads1N, 1)

    # Output freq information to be used for the one gap grouping script
    df_merge.to_csv(Freqs1N,index=False)

    # run grouping script for one gap samples
    os.system("python3 DI_boundary_grouping_kj_v4.py --input_file {0} --output_file {1}".format(Freqs1N, Group1N))

    # Merging read information with NewGap information
    d1g_m = MergeReadsGroups(Reads1N, Group1N, 1)

    # Count the New Gap information for each sample and segment
    d1c = CountGroupedDVGs(d1g_m, 1)

    # Determine whether forward/reverse reads have gap - uses new gap information
    b1 = binocheck(d1c, 1)

    # pass df counts, and frequency/rpkm etc.
    # Calculate rpkm and percent read data using new gap information
    f1 = rpkm_and_percent(b1)

    # write all information to csv files
    f1.to_csv(All1N,index=False)

    #write reduced csv files to be used with figures
    f1r = ReduceDF(f1, 1)
    f1r.to_csv(Final1N, index=False)

else:
    print("No reads with one gap")

if two_count > 0:
    # pass df counts, and frequency/rpkm etc.
    # Calculate rpkm and percent read data using new gap information
    f1 = rpkm_and_percent(b1)

    # write all information to csv files
    f1.to_csv(All1N,index=False)

    #write reduced csv files to be used with figures
    f1r = ReduceDF(f1, 1)
    f1r.to_csv(Final1N, index=False)

    # Name for two gap freq files used for grouping script
    Freqs2N = "{0}/{1}.DVG.freqs.TwoGap.N{2}.Mis{3}.M{4}.csv".format(args.save_dir,
    args.strain,args.gap_size,args.total_mismatch,args.align_length)

    # file name for two gap grouped data- output of the grouping script
    Group2N = "{0}/{1}.DVG.grouped.TwoGap.N{2}.Mis{3}.M{4}.csv".format(args.save_dir,
    args.strain,args.gap_size,args.total_mismatch,args.align_length)

    # Final csv file for two gap data, includes all read info (excessive)
    All2N = "{0}/{1}.DVG.All.TwoGap.N{2}.Mis{3}.M{4}.csv".format(args.save_dir,
    args.strain,args.gap_size,args.total_mismatch,args.align_length)

    # Final csv file with only info needed for figures
    Final2N = "{0}/{1}.DVG.FINAL.TwoGap.N{2}.Mis{3}.M{4}.csv".format(args.save_dir,
    args.strain,args.gap_size,args.total_mismatch,args.align_length)

    # Generate freq files for the two gapped data
    df_merge2 = PrepFreq(Reads2N, 2)

    # Generate a csv files used for the grouping script
    df_merge2.to_csv(Freqs2N,index=False)

    # Group two gapped data - generates new gap information
    os.system("python3 DI_boundary_grouping_2N_kj_v4.py --input_file {0} --output_file {1}".format(Freqs2N, Group2N))

    # Merging read information with NewGap information
    d2g_m = MergeReadsGroups(Reads2N, Group2N, 2)

    # Count the New Gap information for each sample and segment
    d2c = CountGroupedDVGs(d2g_m, 2)

    # Determine whether forward/reverse reads have gap - uses new gap information
    b2 = binocheck(d2c, 2)

    # pass df counts, and frequency/rpkm etc.
    # Calculate rpkm and percent read data using new gap information
    f2 = rpkm_and_percent(b2)

    # write all information to csv files
    f2.to_csv(All2N,index=False)

    #write reduced csv files to be used with figures
    f2r = ReduceDF(f2, 2)
    f2r.to_csv(Final2N, index=False)

else:
    print("No reads with two gaps")
