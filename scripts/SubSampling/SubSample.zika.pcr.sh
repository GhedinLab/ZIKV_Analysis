#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:40:00
#SBATCH --mem=8GB


# subsampling deduplicated zika bam files to look at shannon diversity
# INPUTS (in order):
# FRAC = the fraction of reads to pull from bam bamfiles
# DATA_DIR = the directory path (with / at the end) of where bamfiles are located
# RUNDIR = the directory where to run and put output files

# OUTPUT:
# Sub-sampled alignment files - sorted - and indexed
# Used for variant calling

# to run: sbatch -a 0-241 SubSample.zika.pcr.sh "1" "/scratch/kej310/ZIKA/ZikaDVG/PCR/" "/scratch/kej310/ZIKA/ZikaDVG/" "/scratch/kej310/ZIKA/ZikaDVG/reference/zika/MR766.ZIKA.CDS.fas"
module purge
module load samtools/intel/1.6

# name test
FRAC=$1
#FRAC="0.75"
DATA_DIR=$2
#DATA_DIR='/scratch/kej310/ZIKA/ZikaDVG/test/'
RUNDIR=$3

REFSEQ=$4

F="${DATA_DIR}/*.bam" # pull out names of alignment files

cd ${RUNDIR} # change to run directory

Fs=($(ls $F))  # list bam files

F1=${Fs[$SLURM_ARRAY_TASK_ID]}  # bam name associated with slurm task id

N1=${F1#*"$DATA_DIR"}  #  name of sample used to generate other files downstream

N1=${N1%".bam"}  #  Grab name before the .bam  of file

NAME=$N1

echo $N1

echo $NAME
# count the number of reads in the bamfile
ctit=$(samtools view -c $F1)

echo $ctit

selec=$( bc -l <<< "$FRAC * $ctit") # pull fraction of reads from alignment

echo $selec

rdit=$(printf "%.0f\n" $selec) # round that number

echo $rdit

#  Sub sample the alignment file and generate a sam file
cat <(samtools view -H $F1) <(samtools view -q 255 $F1 | shuf -n $rdit) > ./${NAME}.${FRAC}.sam

# change to bamfile
samtools view -bSq 20 ./${NAME}.${FRAC}.sam > ./${NAME}.${FRAC}.bam

# sort bam file
samtools sort -T ./${NAME}.${FRAC}.sorted \
              -o ./${NAME}.${FRAC}.sorted.bam \
              ./${NAME}.${FRAC}.bam

samtools index ./${NAME}.${FRAC}.sorted.bam

module load pysam/intel/0.10.0

python readreport_v5_1.py --strain "ZIKA_${FRAC}" --infile ./${NAME}.${FRAC}.sorted.bam --ref ${REFSEQ}
