import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
import random

def arguments():
	parser  = argparse.ArgumentParser(description="Filters VCF file for Fijiensis (options are species-specific")
	parser.add_argument("-i", "--input", help="path for input", required=True)
	parser.add_argument("-o","--output", help="path for output", required=True)
	parser.add_argument("-m","--mq", help="the threshold MQ to keep. Default=MQ>=1", default=1)
	parser.add_argument("-d","--lowdepth", help="the lower bound depth to drop. Default=10", default=10)
	parser.add_argument("-D","--highdepth", help="the upper bound depth to drop. Default=50", default=50)
	parser.add_argument("-t","--te", help="path for file containing sites within called TEs")
	parser.add_argument("-r","--nearhet", help="frequency near 0.5 that will be dropped as near-het, e.g. 0.1 means ratio 0.4-0.6 will be dropped", default=0.1)
	parser.add_argument("-f","--readfrequency",help="minimum frequency of reads required to match the genotype call at a site per sample", default=0.9)
	parser.add_argument("-p","--processes",help="number of processes (CPU) to run job on",default=10)
	args = parser.parse_args()
	return(args)


def testfunc(df):
	pos = df[[1]]
	return(pos*pos)

####################
# Begin Script
####################

args = arguments()

#get infile info
input = args.input

#get outpath info
outpath = args.output

#reader = pd.read_table(input, sep='\t', chunksize=2, iterator = True, comment='#', header=None)
pool = mp.Pool(args.processes)

iter = 1

for df in reader:
	result = pool.apply_async(testfunc, [df]).get()
	if iter==1:
		result.to_csv(outpath,header=False,index=False)
	else:
		result.to_csv(outpath,header=False,index=False,mode='a')
	iter+=1

