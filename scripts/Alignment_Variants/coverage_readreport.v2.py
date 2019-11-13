#make coverage files from the readreport files
#python coverage_readreport.v2.py --reflist CAL09,PR8
import os
import glob
import pandas as pd
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--reflist','-rl',required=True,help='Give reference, or list of references used to align delimted by comma') #args.ref
args = parser.parse_args()

def printer(outfile,name,segment,ntpos,totalcount,ref):
    print>>outfile, '%s,%s,%s,%s,%s' % (name,segment,ntpos,totalcount,ref)

HEADER = "name,segment,ntpos,totalcount,ref"
path = 'FILES/fullvarlist' #this could be an input

if not os.path.exists('coverage'):
    os.makedirs('coverage')

refs = [str(item) for item in args.reflist.split(',')]
for ref in refs:
    #path = 'FILES/fullvarlist/{0}'.format(ref)	
    strain = ref.upper()
    print(strain)
    filenames = glob.glob('{0}/*.{1}.*'.format(path,strain))
    outfile = open("{0}.coverage.csv".format(strain),'w')
    print>>outfile,HEADER
    for snp_file in filenames:
        print(snp_file)
        df = pd.read_csv(snp_file,keep_default_na=False)
        for index, row in df.iterrows():
            printer(outfile,row['name'],row['segment'],row['ntpos'],row['totalcount'],strain)

os.system('mv *.coverage.* coverage/')
