#!/bin/bash

module purge
module load pysam/intel/0.10.0
module load star/intel/2.5.3a
module load samtools/intel/1.6
module load picard/2.8.2

python STAR_align_v3.rename.merge.py --ref zika --ref_fasta MR766.ZIKA.CDS.fas -m1 rawfiles/mate_1 -m2 rawfiles/mate_2 --rename RenameFile.csv
python coverage_readreport.v2.py --reflist zika
python ConsensusFasta.py --ref reference/zika/MR766.ZIKA.CDS.fas --var FILES/fullvarlist/ --strain zika
