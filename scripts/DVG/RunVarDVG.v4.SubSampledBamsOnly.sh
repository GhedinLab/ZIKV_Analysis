#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=2GB

###############################################################################${ref} ${RUNDIR} ${STRAIN} ${REFSEQ}
########################## DON'T CHANGE ANY CODE BELOW ########################${REF} ${RUNDIR} ${STRAIN} ${REFSEQ};
ref=$1
RUNDIR=$2
STRAIN=$3
REFSEQ=$4
GTF=$5


echo $ref
echo $RUNDIR
echo $STRAIN
echo $REFSEQ

cd ${RUNDIR}

module purge
module load samtools/intel/1.6
module load pysam/intel/python3.6/0.14.1


mkdir ./N_files
mkdir ./IdxStats
#mkdir ./coverage
#mkdir ./consensus

# Do not need coverage/ consensus information from these samples
#generate consensus and coverage information
#could have a cov cutoff for the consensus seq, the default is 50, parameter --cov
#python3 ConsensusFasta.Coverage.v3.py --ref $REFSEQ --var ./FILES/fullvarlist --strain ${STRAIN}

#mv *.coverage.* ./coverage
#mv *.fasta ./consensus

#pull out read name, flags, ref/seg name, 1-based leftmost mapping position, mapq, cigar
for i in $(ls ../PCR_SubSortedBam);do if [[ $i == *.${STRAIN}.sorted.bam ]]; then samtools view ../PCR_SubSortedBam/$i | awk '$6 ~ /N/' | cut -f 1,2,3,4,5,6 | uniq  > $i.$STRAIN.split.txt;fi;done
mv *.split.* N_files/

for i in $(ls ./N_files);do if [[ $i == *.${STRAIN}.split.txt ]]; then sed -i "s/$/\t${i%".$STRAIN.rmd"*}/" ./N_files/$i;fi;done

cat ./N_files/*${STRAIN}.split.txt > ./${STRAIN}.SplitReads.csv

#change from tab forumalted to csv file
#add header to the single idxstat file:
sed -i 's/\t/,/g' ./${STRAIN}.SplitReads.csv

sed -i 1i"readname,flags,segment,left_pos,mapq,cigar,name" ./${STRAIN}.SplitReads.csv

# using feature counts instead of IdxStats
module load subread/intel/1.5.1

featureCounts -T 4 -B -p -t exon -a ${GTF} -o ./${STRAIN}.FeaturesOutput.txt ../PCR_SubSortedBam/*.${STRAIN}.sorted.bam

sed -i '1d' ./${STRAIN}.FeaturesOutput.txt # remove top line

# change to csv file
sed 's/\t/,/g' ./${STRAIN}.FeaturesOutput.txt > ./${STRAIN}.FeaturesOutput.csv

