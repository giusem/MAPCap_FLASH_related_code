
import os
import re
import gzip
import subprocess
import argparse
import sys

def main():

	parser = argparse.ArgumentParser(description = "Extract reads of the desired length from fastq files.")

	parser.add_argument('-i', '--input', metavar = "fastq file name", required = True, help = "Name of the fastq input file.")
	parser.add_argument('--min', type=int, metavar='min-insert-size', help='The minimum read length you want to keep. Needs --max to work.')
	parser.add_argument('--max', type=int, metavar='max-insert-size', help='The maximum read length you want to keep. Needs --min to work.')

	args = parser.parse_args()

	readExtractor(args.input, args.min, args.max)


def readExtractor(inputname, min, max):


	if min and not max:
		print "min and max option required. quitting."
		sys.exit(1)
	if max and not min:
		print "min and max option required. quitting."
		sys.exit(1)

	inputname = inputname.split("/")[-1]

	input_file = os.path.realpath(inputname)

	save_dir = input_file.split(inputname)[0]

	if re.search('fastq', input_file, re.IGNORECASE):
		if re.search('fastq.gz', input_file, re.IGNORECASE):
			filetype = 'fastq.gz'
		else:
			filetype = 'fastq'
	else:
		print("I don't understand the file type of %s. Your files should end either with .fastq or .fastq.gz" %input)
		sys.exit(1)

# this part is only for counting the lines

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


# now we open the file for the actual work

	if filetype == 'fastq.gz':
		input_reads = gzip.open(input_file, 'rb')

	elif filetype == 'fastq':
		input_reads = open(input_file, 'rb')

	count = 0
	targets_dict = {}
	multiplier = 1
	chunk = 100000
	factor = 1
	status = 10
	target_count = 0
	seq_length_list = []

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

			if min <= seq_length <= max:
				target_count += 1
				try:
					targets_dict[save_dir+inputname.split(".fastq")[0]+"_"+str(min)+"-"+str(max)+".fastq"].append(["%s\n%s\n+\n%s\n%s\n" %(name, seq, plus, qual)])
				except KeyError:
					targets_dict[save_dir+inputname.split(".fastq")[0]+"_"+str(min)+"-"+str(max)+".fastq"] = [["%s\n%s\n+\n%s\n%s\n" %(name, seq, plus, qual)]]

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
			break

	targets_dict = {}

	input_reads.close()

	if os.path.isfile(save_dir+inputname.split(".fastq")[0]+"_"+str(min)+"-"+str(max)+".fastq"):
		cmd = "gzip %s" %save_dir+inputname.split(".fastq")[0]+"_"+str(min)+"-"+str(max)+".fastq"
		subprocess.check_call(cmd, shell=True)

	print "\nOut of %s sequences, %s were between %s and %s nucleotides.\n" %(readcount, target_count, min, max)

if __name__ == "__main__":
	main()
