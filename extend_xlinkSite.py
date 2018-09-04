

import argparse
import sys
import os
import re

def main():

	parser = argparse.ArgumentParser(description = "Extend the xlinking site by symmetrically to the desired length")

	parser.add_argument('-i', '--input', metavar = 'bedGraph file name', required = True, help = "Name of the bedGraph input file.")
	parser.add_argument('-l', '--length', metavar = 'xlinking site length', required = True, help = "The length the xlinking site should be extended to.")

	args = parser.parse_args()

	extendXsite(args.input, args.length)


def extendXsite(input, length):

	input = input.split("/")[-1]

	input_file = os.path.realpath(input)

	save_dir = input_file.split(input)[0]


	if re.search('bedgraph', input_file, re.IGNORECASE):
		pass
	elif re.search('bg', input_file, re.IGNORECASE):
		pass
	else:
		print("I don't understand the file type of %s. Your files should end either with .bedGraph or .bg" %input)
		sys.exit(1)


	input_reads = open(input_file, 'rb')

	length = float(length)

	add_nt = int(round(length / 2))

	output_dict = {}


	while True:
		try:
			line = input_reads.next().strip()

			chromosome = line.split("\t")[0]
			start = int(line.split("\t")[1])
			end = int(line.split("\t")[2])
			score = line.split("\t")[3]

			newstart = start - add_nt
			newend = end + add_nt - 1

			try:
				output_dict[save_dir+input.split(".")[0]+"_extended.bedGraph"].append(["%s\t%s\t%s\t%s\n" %(chromosome, newstart, newend, score)])
			except KeyError:
				output_dict[save_dir+input.split(".")[0]+"_extended.bedGraph"] = [["%s\t%s\t%s\t%s\n" %(chromosome, newstart, newend, score)]]


			for key, value in output_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			output_dict = {}

		except StopIteration:

			for key, value in output_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			print "at the end"

			output_dict = {}

			break

	input_reads.close()


if __name__ == "__main__":
	main()
