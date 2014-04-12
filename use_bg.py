"""
Usage example for BedGraphReader
"""
__author__ = 'eranroz'

import BedGraphReader
import numpy as np
import time


BedGraphReader.VERBOSE = False

def readChromeSize():
	from urllib.request import urlopen
	chrom_sizes = dict()
	with urlopen('http://genome.ucsc.edu/goldenPath/help/hg19.chrom.sizes') as chromsizes_f:
		for line in chromsizes_f:
			splitted=line.decode().split('\t')
			if len(splitted)!=2:
				continue
			chrom, size = splitted[0], splitted[1]
			chrom_sizes[chrom]=int(size)/20
	return chrom_sizes

def readWithBedGraphReader(bg_file):
	chrom_sizes=readChromeSize()
	curr=time.time()
	genome = BedGraphReader.load_bedgraph(bg_file, chrom_sizes)
	return time.time()-curr

def simpleRead(bg_file):
	chrom_sizes = readChromeSize()
	curr=time.time()
	genome = dict()
	last_chrome = None
	with open(bg_file) as f:
		for r in f:
			chrom,start,end,score = r.split('\t')
			if last_chrome!=chrom:
				if last_chrome is None:
					chrom_data = np.zeros(chrom_sizes[chrom])
					last_chrome = chrom
				else:
					genome[chrom] = chrom_data

			start = int(start)//20
			end = int(end)//20
			score = int(score)
			
			for i in range(start, end):
				chrom_data[i] = score

	return time.time()-curr

if __name__ == '__main__':
	bg_file = input("Enter bed graph file name:")
	print('Comparing bedgraph reader to naive python parsing')
	print('Starting bed graph reader')
	time_bg = readWithBedGraphReader(bg_file)
	print('finished... Now starting naive implementation')
	time_naive = simpleRead(bg_file)
	print("Bg: %f\t Naive: %f"%(time_bg, time_naive))
	


