#!/home/semplicio/virtualenv/bin/python -u

	# what this does:
	# from loci provided as a bed file the script extracts the length of the reads mapping to these loci and outputs a csv or tsv file 
	# meanReadLength outputs the mean read length per loci
	# readLengthDist output the number of reads of specific sizes per loci
	# the bed file should be structure as follows chrom, start, end, name, value, strand
import argparse
import os



def main():

	parser = argparse.ArgumentParser(description = "This tool exracts reads from regions provided as bed file and puts them in individual bam file or as one merged bam file.")
	parser.add_argument('-R', '--region', metavar='path/to/regions.bed', help='The path to the bed file containing the regions you would like to have the reads from.', required=True)
	parser.add_argument('-o', '--out', metavar='output file name', help='The name of the output bam file (no file format extension).', required=True)
	parser.add_argument('-i', '--inBAM', metavar='path/to/input.bam', help='The path to the BAM file to extract reads from.', required=True)
	parser.add_argument('-t', '--type', choices=['mean', 'distribution'], help = 'Defines whether you want the "mean" read length per loci or the "distribution" of read lengths', required=True)
	parser.add_argument('--min', type=int, metavar='min-insert-size', help='The minimum read length you want to keep. Needs --max to work.')
	parser.add_argument('--max', type=int, metavar='max-insert-size', help='The maximum read length you want to keep. Needs --min to work.')
	args = parser.parse_args()

	current_dir = os.getcwd()+"/"

	regions = os.path.abspath(args.region)
	outfile = current_dir+args.out+".tsv"
	bamfile = os.path.abspath(args.inBAM)

	if args.type == "mean":
		meanReadLength(regions, outfile, bamfile)
	elif args.type == "distribution":
		readLengthDist(regions, outfile, bamfile, args.min, args.max)

def meanReadLength(regions, outfile, bamfile):
	import subprocess
	import sys
	import pysam
	import numpy as np

	samfile = pysam.AlignmentFile(bamfile, 'rb')
	regionfile = open(regions, 'r')

	length = np.array([])

	out_collector = []

	while True:
		try:

			line = regionfile.next().strip()
			chrom=line.split("\t")[0]
			start=int(line.split("\t")[1])
			end=int(line.split("\t")[2])

			for read in samfile.fetch(chrom, start, end):
				length = np.append(length, read.infer_read_length())

			if np.size(length) == 0:
				mean = 0.0
			else:
				mean = round(length.mean(),1)

			if len(line.split("\t")) >= 4:
				name=line.split("\t")[3]
			else:
				name=False

			out_collector.append([chrom, start, end, name, mean])

			length = np.array([])

		except StopIteration:
			break
			length = np.array([])

	if name:
		with open(outfile, 'a') as ff:
			for i in out_collector:
				ff.write("%s\t%s\t%s\t%s\t%s\n" %(i[0],i[1],i[2],i[3],i[4]))

	else:
		with open(outfile, 'a') as ff:
			for i in out_collector:
				ff.write("%s\t%s\t%s\t%s\n" %(i[0],i[1],i[2],i[4]))

	out_collector = []
	
	regionfile.close()
	samfile.close()

def readLengthDist(regions, outfile, bamfile, min, max):

	if min and not max:
		print "min and max option required. quitting."
		sys.exit(1)
	if max and not min:
		print "min and max option required. quitting."
		sys.exit(1)

	import subprocess
	import sys
	import pysam
	import numpy as np

	samfile = pysam.AlignmentFile(bamfile, 'rb')
	regionfile = open(regions, 'r')

	length = []
	counts = []
	collector = []

	while True:
		try:

			line = regionfile.next().strip()
			chrom=line.split("\t")[0]
			start=int(line.split("\t")[1])
			end=int(line.split("\t")[2])
			name = line.split("\t")[3]
			strand = line.split("\t")[5]

			for read in samfile.fetch(chrom, start, end):
				length.append(read.infer_read_length())

			for i in range(min, max):
				counts.append(length.count(i))

			
			collector.append("%s\t%s\t%s\t%s\t%s\t%s" %(chrom, start, end, name, strand, "\t".join([str(x) for x in counts])))
			
			length = []
			counts = []

		except StopIteration:
			break

	with open(outfile, 'a') as ff:
		for i in collector:
			ff.write("%s\n" %i)


if __name__ == "__main__":
	main()
