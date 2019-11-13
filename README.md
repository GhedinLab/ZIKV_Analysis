# Mapping the evolutionary landscape of Zika virus infection in immunocompromised mice.

Scripts for processing sequencing data collected from mice innoculated with Zika virus. Experiments were performed by the Stapleford Lab at NYU Langone.

Scripts included: 

### Alignment and variant calling:  

+ `STAR.sh` uses the following four scripts, to trim, align, and remove duplicate reads in the sequencing data, call minority variants, pull coverage information, and generate consensus fasta files from the alignments.

  1.  `STAR_align_v3.rename.merge.py` trims, aligns, deduplicates, merges 3 Zika amplicon alignments, and properly names files    
  2.  `readreport_v4_2.py`  calls minor variants for each alignment   
  3.  `coverage_readreport.v2.py` pulls coverage information from snplist files generated using the `readreport_v4_2.py` code   
  4.  `ConsensusFasta.py`   pulls consensus information from the snplist files generated using the `readreport_v4_2.py` code      
  
+ `linegraphPrep_nonsynhack_8.py`  : pulls all minor variants across snplist files (generated using `readreport_v4_2.py`) using preferred thresholds. 

+ `SNVAnalysis.Zika.All.Files.Rmd` and `coverage.v3.zika.R` generate figures using coverage and variant data

### Diversity and Genetic Distance: 
+ `Shannon_Entropy.zika.py` and `Zika.Shannon.Diversity.v2.Rmd`: uses the snplist files (generated using `readreport_v4_2.py`) to calculate Shannon entropy and generate figures.

+ `euc_dist_files_Zika.v4.py` and `euc.dist.zika.v5.Rmd` : pulls data from snplist files (generated using `readreport_v4_2.py`) and calculates the genetic distance using Euclidean distance and generates figures.

### Defective viral genome identification: 
+ `RunDVG.sh` pulls cigar string and alignment information from the alignment files and pipes it into `FindDI_star_7.py` to identify DVG reads.

+ `AppendCoverage.STAR.DVG.pcr.zika.v2.Rmd` and `CompareCluster.ZIKA.pcr.v3.Rmd` Adds coverage information to DVG output files and generates figures.

# Reference:     
In preparation



# Contact:      
Kate Johnson kej310@nyu.edu
