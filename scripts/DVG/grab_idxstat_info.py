import pandas as pd
import numpy as np
import glob
import os

def printer(outfile,name,segment,totalreads):
    print>>outfile,'%s,%s,%s' %(name,segment,totalreads)

path= 'IdxStats'

HEADER = 'name,segment,totalreads'

outfile = open("IdxStats.csv",'w')
print>>outfile,HEADER

for infile in glob.glob( os.path.join(path, '*.idxstats.txt' ) ):#will go through each specified csv file
    print(infile)
    name = infile.split('/')[-1].split('.')[0]
    print(name)
    df = pd.read_csv(infile, sep = '\t', keep_default_na = False,header=None)
    df.columns= ['segment','totalreads']
    df['name'] = name
    print(df)

    for index, row in df.iterrows():
        printer(outfile,row['name'],row['segment'],row['totalreads'])
