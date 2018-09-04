#!/home/semplicio/virtual_environments/ipython/bin/python -u

import Levenshtein
import sys
import re
import gzip

print "absolute path of the R2 read to identify the lollipop barcodes from"

infile = raw_input(">")

ref_barcodes = ["GCTACA", "CAAGTG", "GTCATG", "TTAGCC", "TCCAAC", 
"GTGGAA", "TTGTGG", "TTGGCC", "TCGCTT", "AACTCG", "TGTGAG", "AAGCAC", 
"TGGAAC", "GACACA", "GGTTAC", "CATCAC", "AGAAGG", "TCTACG", "AGATGC", 
"CGCATA", "AGTCGA", "CGTACA"]
ref_barcode_names = ["r1", "r2", "r3", "r4", "r5", 
"r6", "r7", "r8", "r9", "r10", "r11", "r12", 
"r13", "r14", "r15", "r17", "r18", "r19", "r20", 
"r21", "r22", "r23"]

if re.search('fastq', infile, re.IGNORECASE):
	if re.search('fastq.gz', infile, re.IGNORECASE):
		filetype = 'fastq.gz'
	else:
		filetype = 'fastq'
else:
	print("I don't understand your file type. Your files should end either with .fastq or .fastq.gz")
	sys.exit(1)


if filetype == 'fastq.gz':
	fastq = gzip.open(infile, 'rb')
		
elif filetype == 'fastq':
	fastq = open(infile, 'r')
	
bcode_counts = [0 for x in range(len(ref_barcode_names))]

while True:
	try:
		name 	= fastq.next().strip()
		seq		= fastq.next().strip()
		plus 	= fastq.next().strip()
		qual 	= fastq.next().strip()

		bcode = seq[5:11]

		scores = []

		for i in ref_barcodes:
			scores.append(Levenshtein.hamming(i, bcode))


		if min(scores) <= 2 and scores.count(min(scores)) == 1:
			bcode_counts[scores.index(min(scores))] += 1

	except:
		StopIteration

		print "Number of barcodes assigned to which reference:"

		for i in range(len(ref_barcode_names)):
			print "%s\t%s" %(ref_barcode_names[i], bcode_counts[i])

		break


