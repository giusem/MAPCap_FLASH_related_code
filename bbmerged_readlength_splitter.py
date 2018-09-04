#!/home/semplicio/virtualenv/bin/python -u

from __future__ import division
import time
import sys
import os
import argparse
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import gzip
import subprocess
import datetime
import itertools
import string
import textwrap



def main():
	parser = argparse.ArgumentParser(description = 'This tools is thought to be used to:\n'
'1) plot the length distribution of the reads\n'
'2) extract the fragments of the desired lenth interval.\n'
'It takes fastq(.gz) file as input. It was designed for the bbmerged fastq files. But also works on single end reads.\n'
'When --min and --max are given, there will be a fastq.gz output of the desired\n' 
'fragments within these boundaries. \n'
'If --plot is given a plot of the fragment distribution as well as a .txt file\n' 
'with the underlying data is created.\n'
'A good practice is to let the script run one only with the --plot options to\n' 
'get an idea of the distribution. Then you can exract the desired fragments\n' 
'and map them.')
	parser.add_argument('-i', '--input', metavar='path/to/bbmerged/reads', required=True, help='The fastq file of the bbmerged reads.')
	parser.add_argument('--min', type=int, metavar='min-insert-size', help='The minimum read length you want to keep. Needs --max to work.')
	parser.add_argument('--max', type=int, metavar='max-insert-size', help='The maximum read length you want to keep. Needs --min to work.')
	parser.add_argument('--plot', action='store_true', help='Plot the length distribution of all the reads.')
	parser.add_argument('-m', '--mapping', action='store_true', help='Map the extracted reads with HISAT2')
	parser.add_argument('-s', '--species', choices=['hg19', 'dm3', 'dm6', 'mm10', 'mm9', 'dm6_simple'], help='Map to which organism.')
	args = parser.parse_args()

	if args.min and not args.max:
		print "You only defined a minimum insert size but no maximum. Need both. Quitting..."
		sys.exit(1)
	elif args.max and not args.min:
		print "You only defined a maximum insert size but no minimum. Need both. Quitting..."
		sys.exit(1)


	current_dir = os.getcwd()+'/'

	input_file = os.path.realpath(args.input)
	output_name = input_file.split("/")[-1].split(".fastq")[0]

	if re.search('fastq', input_file, re.IGNORECASE):
		if re.search('fastq.gz', input_file, re.IGNORECASE):
			filetype = 'fastq.gz'
		else:
			filetype = 'fastq'
	else:
		print("I don't understand the file type of %s. Your files should end either with .fastq or .fastq.gz" %args.forward)
		sys.exit(1)

	#linecount = 0

	#if filetype == 'fastq.gz':
	#	with gzip.open(args.input, 'r') as ff:
	#		for f in ff:
	#			linecount += 1
	#elif filetype == 'fastq':
	#	with open(args.input, 'r') as ff:
	#		for f in ff:
	#			linecount += 1

	#readcount = int(linecount/4)

	if filetype == 'fastq.gz':
		input_reads = gzip.open(input_file, 'rb')
	elif filetype == 'fastq':
		input_reads = open(input_file, 'rb')

	count = 0
	targets_dict = {}		
	multiplier = 1
	chunk = 100000
	x = 1
	status = 10
	target_count = 0
	seq_length_list = []

	print "Analyzing the reads.\n"
	if os.path.isfile(current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".fastq.gz"):
		print "You ran this script with the same parameters already once here. I am not creating the same fastq.gz file again unless you rename/move/delete the other file."
		pass
	else:
		while True:
			try:
				name = input_reads.next().strip()
				seq = input_reads.next().strip()
				plus = input_reads.next().strip()
				qual = input_reads.next().strip()

				seq_length = len(seq)
				
				seq_length_list.append(seq_length)

				count += 1
				
				#if count >= ((readcount/10)*x):
				#	print str(status)+"%"+" of the reads processed..."
				#	x += 1
				#	status += 10

				
				if args.min and args.max:	
					if args.min <= seq_length <= args.max:
						target_count += 1
						try:
							targets_dict[current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".fastq"].append(["%s\n%s\n+\n%s\n" %(name, seq, qual)])
						except KeyError:
							targets_dict[current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".fastq"] = [["%s\n%s\n+\n%s\n" %(name, seq, qual)]]

					if count == multiplier*chunk:
						for key, value in targets_dict.iteritems():
							with open('%s' %key, 'a+') as ff:
								for i in range(len(value)):
									ff.write('%s' %value[i][0])

						multiplier += 1
						targets_dict = {}

			except StopIteration:
				for key, value in targets_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

					targets_dict = {}

				break

		input_reads.close()

	if os.path.isfile(current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".fastq"):
		cmd = 'gzip %s' %current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".fastq"
		subprocess.check_call(cmd, shell=True)

		print "\nOut of %s sequences, %s were between %s and %s nucleotides.\n" %(count, target_count, args.min, args.max)
	
	if args.plot:
		print "Plotting the distribution of the reads.\n"

		if os.path.isfile(current_dir+output_name+"_length_distribution.txt"):
			print "Generated the count file before. Rename/move/delete old file or I won't generate a new one."
			pass
		else:
			a = []
			b = []
			for i in range(max(seq_length_list)):
				a.append(i)
				b.append(seq_length_list.count(i))
				with open(output_name+"_length_distribution.txt", 'a') as ff:
					ff.write(str(i)+"\t"+str(seq_length_list.count(i))+"\n")

#			plt.figure(figsize=(15,15))
			titlename = output_name
			plt.plot(a,b)
			plt.title("Length distribution of %s" %titlename)
			plt.ylabel("Counts")
			plt.xlabel("Length")
			plt.savefig(output_name+"_length_distribution.pdf")
			plt.close

	if args.mapping:
		print "Mapping the extracted reads.\n"

		if os.path.isfile(current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".bam"):
			print "Mapped this fastq.gz file before. Rename/move/delete the bam file to map this file again."
			pass
		else:
			hisat_input = current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".fastq.gz"

			indices = {'genome':['path/to/hisat_index'], 'hg19':['/data/repository/organisms/GRCh37_ensembl/HISAT2Index/genome'], 'dm3':['/data/akhtar/group/ibonomics/hisat_index/hisat2/dm3/dm3'], 
			'dm6':['/data/repository/organisms/dm6_ensembl/HISAT2Index/genome'], 'mm9':['/data/repository/organisms/GRCm37_ensembl/HISAT2Index/genome'], 'mm10':['/data/repository/organisms/GRCm38_ensembl/HISAT2Index/genome'],
			'dm6_simple':['/data/akhtar/group/Giuseppe/INDEXES/HISAT2index/dm6_simple/dm6_simple']}

			hisat2_path = "/package/hisat2-2.0.0-beta/hisat2"
			samtools_path = "/package/samtools/bin/samtools"
			bamCoverage_path = "/package/deeptools-2.2.2/bin/bamCoverage"	
			
			hisat_options = "--rna-strandness F --max-intronlen 100000"
					
			cmd = "%s %s -p 18 -x %s -U %s 2> %s | %s view -b -F 256 - | %s sort - %s" %(hisat2_path, hisat_options, indices[args.species][0], 
					hisat_input, current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".hisat.log", samtools_path, samtools_path, current_dir+output_name+"_"+str(args.min)+"-"+str(args.max))
			
			subprocess.check_call(cmd, shell=True)

			cmd3 = "%s index %s" %(samtools_path, current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".bam")
			subprocess.check_call(cmd3, shell=True)
			
			cmd4 = "%s -p 18 -b %s -o %s -bs 1 " %(bamCoverage_path, current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".bam", current_dir+output_name+"_"+str(args.min)+"-"+str(args.max)+".bw")
			subprocess.check_call(cmd4, shell=True)


if __name__ == "__main__":
	start_time = time.time()
	print "\n%s Initializing the script.\n" %(time.strftime("%X %x"))
	main()
	end_time = time.time() - start_time
	print "\n%s Done. Elapsed time: %s seconds." %(time.strftime("%X %x"), round(end_time))
