#!/bin/bash

module purge
module load pysam/intel/0.10.0
module load samtools/intel/1.6


mkdir N_files
mkdir IdxStats
for i in $(ls ./bamfiles/rmdup_bams); do if [[ $i == *rmd.star.bam ]]; then samtools view ./bamfiles/rmdup_bams/$i | awk '$6 ~ /N/' | cut -f 1,3,4,6 |uniq > $i.split.txt;fi;done
mv *.split.* N_files/
for i in $(ls ./bamfiles/rmdup_bams); do if [[ $i == *rmd.star.bam ]]; then samtools idxstats ./bamfiles/rmdup_bams/$i | cut -f 1,3 | uniq > $i.idxstats.txt;fi;done
mv *.idxstats.* IdxStats/
python grab_idxstat_info.py
python3 FindDI_star_7.py --strain ZIKA --file_direct N_files --ref ZIKA
