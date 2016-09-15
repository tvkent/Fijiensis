import argparse

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
	args = parser.parse_args()
	return(args)

def emcuezero(line,mq):
	'''
	this line will be tagged for not passing the filter if
	MQ tag is less than the threshold value
	'''

	info = line.split()[7]
	
	if "MQ=" in info:
		sinfo = info.split(';')
		metric = [n for n in sinfo if "MQ=" in n]
		metric = metric[0][3:]
		if float(metric) >= mq:
			filter = False
		else:
			filter = True
	else:
		filter = True

	return(filter)

def depth(line,d,D):

	'''
	this will change the read info for a sample to missing if the depth
	is too low or too high
	'''

	filter = False
	sline = line.split()
	try:
		index = [i for i, s in enumerate(sline[8].split(':')) if 'DP' in s][0]
	except:
		filter = True
		print('DP field not found, dropping site')
	
	if filter==False:
		#loop over the sample info and record index of filtered sample
		samples = sline[9:]
		i=0
		drops=[]
		for sample in samples:
			if '.' not in sample and '/' not in sample:
				dp = sample.split(':')[index]
			
				if float(dp) > d and float(dp) < D:
					sample=sample
				else:
					sample='.'
					drops.append(i)
			i+=1

		#loop over filter indices and change sample info to missing for these

		for drop in drops:
			samples[drop] = '.'

		nonsamples = "\t".join(sline[0:9])
		samples = "\t".join(samples)
		line = nonsamples + '\t' + samples
		filter = False

	return(line, filter)
	
def parse_te(te):
	'''
	open the TE file from jean, parse, and return a usable dictionary of sites per scaffold
	'''

	lines = open(te)
	lines = lines.readlines()
	sites = {}

	#loop over the lines in the file and make dictionary entry of each
	for line in lines:
		line = line.replace("\n","").replace("'","").replace("]","").replace("[","")
		linesstr = str(line).split(',')
		scaf = linesstr[0].split()[0] #get scaffold info (one per line)
		
		pos = [info.lstrip().split()[1] for info in linesstr] #keep only pos info
		
		sites.update({scaf: pos}) #add list to dict entry scaf
	
	return(sites)

def transposno(line, tedict):
	
	'''
	take the scaffold and position of line, if in the tedict, drop the site
	'''

	sline = line.split()
	scaf = sline[0]
	pos = sline[1]
	
	if pos in tedict[scaf]:
		filter = True
	else:
		filter = False
	return(filter)

def readfrequencies(line, readfreq, nearhet):
	
	'''
	drop a site if any one of the samples has a read frequency below the threshold for its called genotype
	also drop if any one of the samples has a read frequency too near 0.5 (probably a bad het call)
	will also drop if any of them have a het call, i.e. GT=./.
	'''
	
	highratio = nearhet + 0.5
	lowratio = 0.5 - nearhet
	sline = line.split()
	info = sline[8].split(':')
	
	hasAD = [field for field in info if field=='AD']
	GTindex = info.index('GT')				
	samples = sline[9:]
	
	if hasAD:
		ADindex = info.index('AD')

		for sample in samples:
			
			if sample!='.' and sample!='./.':	
				ADinfo = sample.split(':')[ADindex]
				if ADinfo!= '.' and ADinfo!='./.' and ',' in ADinfo:
				
					#filter for read frequency
					GTinfo = int(sample.split(':')[GTindex])
					callreads = ADinfo.split(',')[GTinfo]
					denominator = [int(value) for value in ADinfo.split(',')]
					readratio = int(callreads)/sum(denominator)
				
					if readratio < readfreq:
						filter = True
						break
					else:
						filter = False

					#filter for read ratio
					if readratio >= lowratio and readratio <= 0.5:
						filter=True
						break
					elif readratio <= highratio and readratio >= 0.5:
						filter = True
						break
					else:
						filter=False

				else:
					filter = False

			else:
				filter=False
			
	else:
		for sample in samples:
			GTinfo = sample.split(':')[GTindex]
			if '/' in GTinfo:
				filter=True
			else:
				filter=False

	return(filter)

def dropindels(line):

	'''
	automatically drop indels--as defined by a ref or alt call that is longer than one base pair
	'''

	sline = line.split()
	ref = sline[3]
	alt = sline[4]

	filter = False

	#check if multiallelic
	if ',' in ref:
		sref = ref.split(',')
		for call in sref:
			if len(call)>1:
				filter = True
				break
			else:
				filter = False
	else:
		if len(ref)>1:
			filter = True

	if ',' in alt and filter==False:
		salt = alt.split(',')
		for call in salt:
			if len(call)>1:
				filter = True
				break
			else:
				filter = False
	else:
		if len(alt)>1:
			filter = True


	return(filter)

def drophets(line):

	'''
	automatically drop lines with het calls 
	'''

	sline = line.split()
	samples = sline[9:]
	filter = False

	filtervalues = [True if '/' in sample else False for sample in samples]

	if True in filtervalues:
		filter = True
	
	return(filter)

#############
#Begin
#############

args = arguments()

#get infile info
input = args.input
infile = open(input, 'r')
lines = infile.readlines()

#get outpath info
outpath = args.output
outfile = open(outpath,'w')

#if te file provided, parse and return dictionary
if args.te:
	te = args.te
	tedict = parse_te(te)

for line in lines:
	
	if "#" not in line:
		filter = False
		#run filter on TEs
		if args.te:
			filter = transposno(line, tedict)

		#run filter on indels
		if filter == False:
			filter = dropindels(line)

		#run filter on hets
		if filter == False:
			filter = drophets(line)

		#run filter on MQ
		if filter == False:
			mq = args.mq
			filter = emcuezero(line,mq)

		#run depth filter
		if filter == False:
			d = args.lowdepth
			D = args.highdepth
			line, filter = depth(line,d,D)
		
		#run filter on readfrequencies
		if filter == False:
			nearhet = args.nearhet
			readfreq = args.readfrequency
			filter = readfrequencies(line, readfreq, nearhet)

		#either print line or skip if filtered out
		if filter==False:
			outfile.write(line + '\n')
	else:
		outfile.write(line)

outfile.close()
