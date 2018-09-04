#!/home/semplicio/virtualenv/bin/python -u

	# what this does
	# takes bam file and extracts specific loci provided in a bed file
	# output as individual bams or as one merged file

import argparse
import os



def main():

	parser = argparse.ArgumentParser(description = "This tool exracts reads from regions provided as bed file and puts them in individual bam file or as one merged bam file.")
	parser.add_argument('-R', '--region', metavar='path/to/regions.bed', help='The path to the bed file containing the regions you would like to have the reads from.', required=True)
	parser.add_argument('-o', '--out', metavar='output file name', help='The name of the output bam file (no file format extension).', required=True)
	parser.add_argument('-b', '--inBAM', metavar='path/to/input.bam', help='The path to the BAM file to extract reads from.', required=True)
	parser.add_argument('-t', '--type', metavar='type of output BAM', choices=['individual','merged','both'], help='Whether you want INDIVIDUAL BAM files of each region, a single BAM file with all reads MERGED, or BOTH', required=True)
	args = parser.parse_args()

	current_dir = os.getcwd()+"/"

	inputfile = os.path.abspath(args.region)
	outfile = current_dir+args.out
	bamfile = os.path.abspath(args.inBAM)

	bamLociMerger(inputfile, outfile, bamfile, args.type)

def bamLociMerger(inputfile, outfile, bamfile, typeofbam):
	import subprocess
	import sys

	if typeofbam == "individual":
		print ""
	elif typeofbam == "merged":
		print ""
	elif typeofbam == "both":
		print ""
	else:
		print "type 'individual', 'merged' or 'both'. No other option allowed here."
		sys.exit(1)

	infile = open(inputfile, 'r')

	i=1

	while True:
		try:

			line = infile.next().strip()
			chrom=line.split("\t")[0]
			start=line.split("\t")[1]
			end=line.split("\t")[2]

			if (typeofbam == "merged") or (typeofbam == "both"):
				if i == 1:
					cmd = "samtools view -h %s %s:%s-%s >> %s" %(bamfile, chrom, start,end, outfile+"_merged.sam")
					subprocess.check_call(cmd, shell=True)
				else:
					cmd = "samtools view %s %s:%s-%s >> %s" %(bamfile, chrom, start,end, outfile+"_merged.sam")
					subprocess.check_call(cmd, shell=True)

			if (typeofbam == "individual") or (typeofbam == "both"):
				cmd = "samtools view -hb %s %s:%s-%s > %s" %(bamfile, chrom, start,end, outfile+str(i)+".bam")
				subprocess.check_call(cmd, shell=True)

			i+=1

		except StopIteration:


			if (typeofbam == "merged") or (typeofbam == "both"):
				cmd = "samtools view %s %s:%s-%s >> %s" %(bamfile, chrom, start,end, outfile+"_merged.sam")
				subprocess.check_call(cmd, shell=True)

			if (typeofbam == "individual") or (typeofbam == "both"):
				cmd = "samtools view -hb %s %s:%s-%s > %s" %(bamfile, chrom, start,end, outfile+i+".bam")
				subprocess.check_call(cmd, shell=True)
			break

	if (typeofbam == "merged") or (typeofbam == "both"):
		cmd = "samtools view -bS %s > %s" %(outfile+"_merged.sam", outfile+"_merged.bam")
		subprocess.check_call(cmd, shell=True)
		cmd2 = "samtools sort %s %s" %(outfile+"_merged.bam", outfile+"_sorted")
		subprocess.check_call(cmd2, shell=True)
		cmd3 = "samtools index %s" %(outfile+"_sorted.bam")
		subprocess.check_call(cmd3, shell=True)
		cmd4 = "bamCoverage -b %s -o %s" %(outfile+"_sorted.bam", outfile+".bw")
		subprocess.check_call(cmd4, shell=True)
		cmd5 = "rm %s %s" %(outfile+"_merged.sam", outfile+"_merged.bam")
		subprocess.check_call(cmd5, shell=True)

	print "done...good bye"

if __name__ == "__main__":
	main()


