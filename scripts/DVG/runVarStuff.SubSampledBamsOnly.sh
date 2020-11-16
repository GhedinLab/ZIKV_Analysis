USER='kej310@nyu.edu' #email for running information
RUNDIR="/scratch/kej310/ZIKA/ZikaDVG/SubDVG_PCR" #directory where code should run
REF="/scratch/kej310/ZIKA/ZikaDVG/reference/zika" #where indexed reference located for star

# reference information
# strain information in the case of subsampling is the percentage subsampled
#STRAIN="1" #lowercase strain name - helps if you align one sample to multiple refs MUST be lowercase
REFSEQ="/scratch/kej310/ZIKA/ZikaDVG/reference/zika/MR766.ZIKA.CDS.fas" # reference fasta
GTF="/scratch/kej310/ZIKA/ZikaDVG/reference/zika/MR766.ZIKA.CDS.gtf"

################################################

module purge

sbatch --mail-type=END --mail-user=$USER  ./RunVarDVG.v4.SubSampledBamsOnly.sh ${REF} ${RUNDIR} "1" ${REFSEQ} ${GTF}
sbatch --mail-type=END --mail-user=$USER  ./RunVarDVG.v4.SubSampledBamsOnly.sh ${REF} ${RUNDIR} "0.75" ${REFSEQ} ${GTF}
sbatch --mail-type=END --mail-user=$USER  ./RunVarDVG.v4.SubSampledBamsOnly.sh ${REF} ${RUNDIR} "0.50" ${REFSEQ} ${GTF}
sbatch --mail-type=END --mail-user=$USER  ./RunVarDVG.v4.SubSampledBamsOnly.sh ${REF} ${RUNDIR} "0.25" ${REFSEQ} ${GTF}
sbatch --mail-type=END --mail-user=$USER  ./RunVarDVG.v4.SubSampledBamsOnly.sh ${REF} ${RUNDIR} "0.10" ${REFSEQ} ${GTF}
sbatch --mail-type=END --mail-user=$USER  ./RunVarDVG.v4.SubSampledBamsOnly.sh ${REF} ${RUNDIR} "0.05" ${REFSEQ} ${GTF}
