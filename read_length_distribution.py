

import gzip
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import re

def main():

	parser = argparse.ArgumentParser(description = "Calculate the read length of fastq entries and plot the distribution.")
	parser.add_argument('-i', '--input', metavar = 'fastq file name', required = True, help = "Name of the fastq input file.")
	args = parser.parse_args()

	readLength(args.input)


def readLength(input):

	current_dir = os.path.abspath(".")

	input_file = os.path.abspath(input)
	input_name = input.split("/")[-1].split(".fastq.gz")[0]

	if re.search('fastq', input_file, re.IGNORECASE):
		if re.search('fastq.gz', input_file, re.IGNORECASE):
			filetype = 'fastq.gz'
		else:
			filetype = 'fastq'
	else:
		print("I don't understand the file type of %s. Your files should end either with .fastq or .fastq.gz" %input)
		sys.exit(1)


	linecount = 0

	if filetype == 'fastq.gz':
		with gzip.open(input_file, 'r') as ff:
			for f in ff:
				linecount += 1
	elif filetype == 'fastq':
		with open(input_file, 'r') as ff:
			for f in ff:
				linecount += 1

	readcount = int(linecount/4)

	if filetype == 'fastq.gz':
		input_reads = gzip.open(input_file, 'rb')
	elif filetype == 'fastq':
		input_reads = open(input_file, 'rb')

	count = 0
	factor = 1
	status = 10
	seq_length_list = []

	print "Analyzing the reads length distribution...\n"

	while True:
		try:
			name = input_reads.next().strip()
			seq = input_reads.next().strip()
			plus = input_reads.next().strip()
			qual = input_reads.next().strip()

			seq_length = len(seq)

			seq_length_list.append(seq_length)

			count += 1

			if count >= ((readcount/10)*factor):
				print str(status)+"%"+" of the reads processed..."
				factor += 1
				status += 10

		except StopIteration:
			break

	input_reads.close()


	x = []
	y = []

	for i in range(max(seq_length_list)+10):
		x.append(i)
		y.append(seq_length_list.count(i))
		with open(current_dir+"/"+input_name+"_read_length_distribution.txt", 'a') as ff:
			ff.write(str(i)+"\t"+str(seq_length_list.count(i))+"\n")

	plt.plot(x,y)
	plt.title("Length distribution ")
	plt.ylabel("Counts")
	plt.xlabel("Read Length")
	plt.savefig(current_dir+"/"+input_name+"_read_length_distribution.png")
	plt.close

	seq_length_list = []

if __name__ == "__main__":
	main()
