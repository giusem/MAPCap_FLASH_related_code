#!/home/semplicio/virtualenv/bin/python -u

# This script uses the regions given in a bed file to extract all reads from a bam file that are overlapping with those regions
# afterwards the readlength is calculated and plotted

import os
import argparse
import subprocess
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def bamLociMerger(inputfile, outfile, bamfile, paired_end):


	if paired_end:

		genome = "/data/repository/organisms/dm6_ensembl/genome_fasta/genome.chrom.sizes"
		
		temp1 = bamfile.split(".bam")[0]+"_mappedPairsOnly.bam"
		temp2 = bamfile.split(".bam")[0]+"_mappedPairsOnly_NameSorted"
		temp3 = bamfile.split(".bam")[0]+".bed"
		temp4 = bamfile.split(".bam")[0]+"_se.bed"
		temp5 = bamfile.split(".bam")[0]+"_se.bam"
		temp6 = bamfile.split(".bam")[0]+"_se_sorted"

		print "Keeping only mapped reads of proper pairs"

		if os.path.exists(temp6+".bam"):
			print "Paired end mode \n Using the bam file from a previous run. \n If you don't want that then (re)move the file."
			pass
		else:
			cmd = "samtools view -b -F 4 -f 0x2 %s > %s" %(bamfile, temp1) #keep only mapped reads and proper pairs
			subprocess.check_call(cmd, shell=True)
			
			print "Sorting new bam file by name"
			cmd = "samtools sort -n %s %s" %(temp1, temp2) #sort by name
			subprocess.check_call(cmd, shell=True)
			
			print "bamToBed"
			cmd = "bamToBed -i %s -bedpe > %s" %(temp2+".bam", temp3) 
			subprocess.check_call(cmd, shell=True)
			
			print "Making single-end bed file"
			cmd = "awk '{print$1,$2,$6,$7,$8,\".\"}' %s | tr [:blank:] \\\\t > %s" %(temp3, temp4)
			subprocess.check_call(cmd, shell=True)

			print "Creating single-end bam file"
			cmd = "bedToBam -i %s -g %s > %s" %(temp4, genome, temp5)
			subprocess.check_call(cmd, shell=True)

			print "Sorting and indexing single-end bam file"
			cmd = "samtools sort %s %s" %(temp5, temp6)
			subprocess.check_call(cmd, shell=True)
			cmd = "samtools index %s" %temp6+".bam"
			subprocess.check_call(cmd, shell=True)

			garbagelist = [temp1, temp2+".bam", temp3, temp4, temp5]
			for i in garbagelist:
				os.remove(i)

		bamfile = temp6+".bam"

	
	print "Extracting the reads of the specific loci"

	infile = open(inputfile, 'r')
	
	i=1

	cmd = "samtools index %s" %bamfile
	subprocess.check_call(cmd, shell=True)

	while True:
		try:

			line = infile.next().strip()
			chrom=line.split("\t")[0]
			start=line.split("\t")[1]
			end=line.split("\t")[2]

			if i == 1:
				cmd = "samtools view -h %s %s:%s-%s >> %s" %(bamfile, chrom, start, end, outfile+"_temp.sam")
				subprocess.check_call(cmd, shell=True)
			else:
				cmd = "samtools view %s %s:%s-%s >> %s" %(bamfile, chrom, start, end, outfile+"_temp.sam")
				subprocess.check_call(cmd, shell=True)

			i+=1

		except StopIteration:


			cmd = "samtools view %s %s:%s-%s >> %s" %(bamfile, chrom, start,end, outfile+"_temp.sam")
			subprocess.check_call(cmd, shell=True)

			
			break

	cmd = "samtools view -bS %s > %s" %(outfile+"_temp.sam", outfile+"_temp.bam")
	subprocess.check_call(cmd, shell=True)
	cmd2 = "samtools sort %s %s" %(outfile+"_temp.bam", outfile)
	subprocess.check_call(cmd2, shell=True)
	cmd3 = "samtools index %s" %(outfile+".bam")
	subprocess.check_call(cmd3, shell=True)
	cmd5 = "rm %s %s" %(outfile+"_temp.sam", outfile+"_temp.bam")
	subprocess.check_call(cmd5, shell=True)

	lociBam = outfile+".bam"

	return lociBam

def bamToReadLength(inputBam):
	
	
	bamtobed_out = inputBam.split(".bam")[0]+".bed"
	cmd = "bamToBed -i %s > %s" %(inputBam, bamtobed_out)
	subprocess.check_call(cmd, shell=True)

	print "done making BED file"

	inputBed = open(bamtobed_out, 'rb')

	
	seq_length_list = []

	while True:
		try:
			line = inputBed.next().strip()
			start = int(line.split("\t")[1])
			end = int(line.split("\t")[2])
			length = end - start

			seq_length_list.append(length)

		except StopIteration:
			break

	inputBed.close()

	print "Plotting the distribution of the reads.\n"

	plot_name = bamtobed_out.split(".bed")[0]+"_lengthDistribution.pdf"
	a = []
	b = []
	for i in range(max(seq_length_list)):
		a.append(i)
		b.append(seq_length_list.count(i))
		with open(plot_name, 'a') as ff:
			ff.write(str(i)+"\t"+str(seq_length_list.count(i))+"\n")

#	plt.figure(figsize=(15,15))
	titlename = bamtobed_out.split("/")[-1].split(".bed")[0]
	plt.plot(a,b)
	plt.title("Length distribution of %s" %titlename)
	plt.ylabel("Counts")
	plt.xlabel("Length")
	plt.savefig(plot_name)
	plt.close



def main():
	parser = argparse.ArgumentParser(description="This script uses the regions given in a bed file to extract all reads from a bam file that are overlappint with those regions. Afterwards the readlength is calculated and plotted")
	parser.add_argument('-R', '--region', metavar='path/to/regions.bed', help='The path to the bed file containing the regions you would like to have the reads from.', required=True)
	parser.add_argument('-o', '--out', metavar='output file name', help='The name of the output bam file (no file format extension).', required=True)
	parser.add_argument('-b', '--inBAM', metavar='path/to/input.bam', help='The path to the BAM file to extract reads from.', required=True)
	parser.add_argument('-pe', '--paired_end', action='store_true', help='In case of paired-end data.')
	args = parser.parse_args()

	current_dir = os.getcwd()+"/"

	inputfile = os.path.abspath(args.region)
	outfile = current_dir+args.out
	bamfile = os.path.abspath(args.inBAM)


	print "Extracting the reads from the regions defined by the bed file."
	lociBam = bamLociMerger(inputfile, outfile, bamfile, args.paired_end)
	print "Calculating the readlength."
	bamToReadLength(lociBam)

	print "Done. Have a lovely day."


if __name__ == "__main__":
	main()

