"""
For alignment with STAR aligner
DIRECTORIES MUST HAVE:

A rawfiles directory with the mates divided into two separate directories such as:
rawfiles/mate1
rawfiles/mate2

A reference directory containing the reference and the indexed reference:
reference/ref.fasta

For start to generate a reference idx do the following:
Run: GenerateIndex.sh OR

STAR --runThreadN 12 \
	--runMode genomeGenerate \
	--genomeDir ./ \
	--genomeFastaFiles ./cal_07_09_seg_cds_MP_NS.fasta \
	--genomeSAindexNbases 6

To run:
python STAR_align_2ref.py --ref1 reference_dir1,shortname1 --ref2 reference_dir2,shortname2
			-m1 rawfiles/mate_1 -m2 rawfiles/mate_2 -rename renamefile.csv

Rename file must have the following columns:
OldName NewName

If you have multiple strains you want to align to you should indicate the reference
folder as the strain name shorthand.

reference/h1n1
reference/yam
etc.
"""

import os
import subprocess
from subprocess import call
import glob
import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--ref','-r',required=True,help='Give directory where first reference is located and shorthand name for ref1') #args.ref
parser.add_argument('--ref_fasta','-r_fa',required=True,help='Give directory where first reference is located and shorthand name for ref1') #args.ref
parser.add_argument('--raw_m1','-m1',required=True,help='Give directory where read 1 is located')
parser.add_argument('--raw_m2','-m2',required =True,help='Give directory where read 2 is located')
parser.add_argument('--rename','-name',required =True,help='Do reads need to be renamed?')
args = parser.parse_args()

def GetNames(direc_list):
    """
    INPUT: list of files in directory
    OUTPUT: list of names, and dictionary of names and filename
    """
    mate_dict = {}
    name_list = set()
    count = 0
    for filename in direc_list:
	#print(filename)
	name = filename.split('.')[1] #2,3,4,5
	#print(name)
        if name in name_list:
            count += 1
            name = '{0}_{1}'.format(name,str(count))
            name_list.add(name)
            mate_dict[name] = filename

        elif name not in name_list:
            name_list.add(name)
            mate_dict[name] = filename
    return name_list,mate_dict

def trim(pair1_dir,pair2_dir,mate1_filename,mate2_filename,name):
	"""
	INPUT: read 1 directory and dictionary that contains the filename corresponding
	to the short name given to each file. The name of each file.
	OUTPUT: trimmed reads.
	"""
	files = ["./{0}/{1}".format(pair1_dir,mate1_filename),"./{0}/{1}".format(pair2_dir,mate2_filename)]
	exists = [f for f in files if os.path.isfile(f)]
	if len(exists) == 2:
		trimmit = "java -jar /share/apps/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 -threads 20 ./{0}/{1} ./{2}/{3} ./rawfiles/trimmed_mate1/{4}_trimmed_1.fq ./rawfiles/unpair/{4}.unpair_trimmed_1.fq ./rawfiles/trimmed_mate2/{4}_trimmed_2.fq ./rawfiles/unpair/{4}.unpair_trimmed_2.fq ILLUMINACLIP:/share/apps/trimmomatic/0.36/adapters/NexteraPE-PE.fa:2:30:10:8:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20".format(pair1_dir, mate1_filename,pair2_dir, mate2_filename,name)
		#print(trimmit)
		os.system(trimmit)
	else: #if the files don't exist - then stop running immediately.
		print("Check read 1 and read 2 files - files don't appear to exist")
		sys.exit()


def align(reference_dir,name,strain_name):
	"""
	INPUT: Trimmed reads, and reference file found in the reference directory [make sure to index ref]
	OUTPUT: Alignment files and log files.
	"""
	strain = strain_name.lower()
	files = ["./rawfiles/trimmed_mate1/{0}_trimmed_1.fq".format(name),
			"./rawfiles/trimmed_mate2/{0}_trimmed_2.fq".format(name),
			"./reference/{0}/chrName.txt".format(strain)]
	#print(files)
	exists = [f for f in files if os.path.isfile(f)]
	#print(exists)
	if len(exists) == 3:
		alignit = "STAR --runThreadN 12 --genomeDir ./reference/{0} --readFilesIn ./rawfiles/trimmed_mate1/{1}_trimmed_1.fq ./rawfiles/trimmed_mate2/{1}_trimmed_2.fq --outReadsUnmapped Fastx --outFileNamePrefix {1}.{2}.".format(reference_dir,name,strain)
		#print(alignit)
		os.system(alignit)
	else:
		print("Check your trimmed reads and your reference directory")
		sys.exit()


def sort_rmdup(name,strain_name):
	"""
	INPUT: Alignment file from the STAR aligner.
	OUTPUT: sorted, and duplicate removed alignment file.
	"""
	strain = strain_name.lower()
	exists = os.path.isfile("{0}.{1}.Aligned.out.sam".format(name,strain))
	if exists:
		sort_and_rmdup = ("samtools view -bSq 20 {0}.{1}.Aligned.out.sam > bamfiles/{0}.{1}.star.bam".format(name,strain),
		"samtools sort -T bamfiles/{0}.{1}.sorted -o bamfiles/sorted_bams/{0}.{1}.sorted.star.bam bamfiles/{0}.{1}.star.bam".format(name,strain),
		"samtools index bamfiles/sorted_bams/{0}.{1}.sorted.star.bam".format(name,strain),
		"java -jar $PICARD_JAR MarkDuplicates I=bamfiles/sorted_bams/{0}.{1}.sorted.star.bam O=bamfiles/rmdup_bams/{0}.{1}.rmd.star.bam M=MetricFiles/{0}.{1}.met.star.txt REMOVE_DUPLICATES=true".format(name,strain),
		"samtools index bamfiles/rmdup_bams/{0}.{1}.rmd.star.bam".format(name,strain))
		#print(sort_and_rmdup)
		runit = [os.system(script) for script in sort_and_rmdup]
	else:
		print("Check to see if alignment worked")
		sys.exit()

def snplist_files(name,strain_name,ref_fasta):
	"""
	INPUT: sorted, deduplicated, alignment file.
	OUTPUT: Tim's snplist files
	"""
	strain = strain_name.lower()
	snplistit = "python readreport_v4_2.py --strain {0} --infile bamfiles/rmdup_bams/{1}.{0}.merged.sorted.rmd.star.bam --ref reference/{0}/{2}".format(strain,name,ref_fasta)
	#print(snplistit)
	os.system(snplistit)

def MergeFiles(name,strain_name,merge_list):
	"""
	INPUT: filenames to MergeFiles
	OUTPUT: merged bamfile
	"""
	strain = strain_name.lower()

	path= 'bamfiles/rmdup_bams/'
	ending = '.{0}.rmd.star.bam'.format(strain)

	NewFileList = ['{0}{1}{2}'.format(path,filename,ending) for filename in merge_list]
	print(NewFileList)
	JoinedFiles = ' '.join(NewFileList)
	CommaFiles = ','.join(NewFileList)
	mergeit = "samtools merge -f bamfiles/rmdup_bams/{0}.{2}.merged.sorted.rmd.star.bam {1}".format(name,JoinedFiles,strain)
	#print(mergeit)
	os.system(mergeit)
	return CommaFiles

def GrabMergeList(namelist,MergeSet):
	"""
	INPUT: name list from directories
	OUTPUT: List of merged names to iterate through
	"""
	OutputList = set()
	for name in namelist:
		OutputList.add(MergeSet[name])
	OutputList = list(OutputList)
	return OutputList

##################################################
directory_list = ['bamfiles','bamfiles/sorted_bams','bamfiles/rmdup_bams','MetricFiles',
                    'rawfiles/trimmed_mate1','rawfiles/trimmed_mate2','rawfiles/unpair',
                    'LogFiles','unmapped']

for direc in directory_list: #make the directories necessary to run script
    if not os.path.exists(direc):
        os.makedirs(direc)


r = [str(item) for item in args.ref.split(',')]
ref_fa = [str(item) for item in args.ref_fasta.split(',')]
if len(r) == len(ref_fa):
	ref_dict = dict(zip(r,ref_fa))

else:
	print("You need to have a directory for each reference used")
	sys.out()

pair_1_dir = args.raw_m1
pair_2_dir= args.raw_m2
mate_1 = os.listdir(pair_1_dir) #list the files in the directory
mate_2 = os.listdir(pair_2_dir)

if len(mate_1) == len(mate_2):
	mate_1 = sorted(mate_1)
	mate_2 = sorted(mate_2)
	name_list, mate_1_dict = GetNames(mate_1)
	name_list2, mate_2_dict = GetNames(mate_2)

elif len(mate_1) != len(mate_2):
	print('There is an error with the number of files- check before continuing')
	sys.exit()

rename_dict={}
merged_dict= {}
MergeSetDict = {}
renamedf = pd.read_csv(args.rename,keep_default_na=False) #input the df with new and oldnames
for index,row in renamedf.iterrows():
	rename_dict[row['OldName']] = row['NewName'] #requires a file that inputs the oldname and recognizes new name
	MergeSetDict[row['OldName']] = row['MergedName']
	if row['MergedName'] not in merged_dict:
		merged_dict[row["MergedName"]] = []
		merged_dict[row['MergedName']].append(row['NewName'])
	else:
		merged_dict[row['MergedName']].append(row['NewName'])

#print(merged_dict)
merge_list = GrabMergeList(name_list,MergeSetDict)
print(merge_list)

with open('CheckFiles.star.csv','w') as f, open('MergedFiles.star.csv','w') as m:
	f.write('name,read1,read2,ref\n')
	m.write('name,mergelist,ShouldMerge\n')
	for OldName in name_list:
		name = rename_dict[OldName]
		read_trimming = trim(pair_1_dir,pair_2_dir,mate_1_dict[OldName],mate_2_dict[OldName],name)
		for reference_name in r:
			reference = reference_name.lower()
			f.write(OldName + ',' + name + ',' + mate_1_dict[OldName] + ',' + mate_2_dict[OldName] + ',' + reference +'\n')
			read_align = align(reference,name,reference)
			sortit = sort_rmdup(name,reference)

	for name in merge_list:
		for reference_name in r:
			reference = reference_name.lower()
			MergeList = merged_dict[name]
			mergeit = MergeFiles(name,reference,MergeList)
			#print(mergeit)
			m.write(name+','+mergeit+'\n')
			sortit = "samtools sort -T bamfiles/rmdup_bams/{0}.{1}.merged.sorted.rmd -o bamfiles/rmdup_bams/{0}.{1}.merged.sorted.rmd.star.bam bamfiles/rmdup_bams/{0}.{1}.merged.rmd.star.bam".format(name,reference)
			os.system(sortit)
			indexit= "samtools index bamfiles/rmdup_bams/{0}.{1}.merged.sorted.rmd.star.bam".format(name,reference)
			os.system(indexit)
			variants = snplist_files(name,reference,ref_dict[reference])

f.close()
m.close()


os.system('mv *.Unmapped.* unmapped/')
os.system('mv *Log.final.out* LogFiles/')
os.system('mv *Log.out* LogFiles/')
os.system('rm *.progress.out')
os.system('rm *.tab')
os.system('rm *.sam')
os.system('rm bamfiles/*.star.bam')
os.system('rm -r *tmp')
