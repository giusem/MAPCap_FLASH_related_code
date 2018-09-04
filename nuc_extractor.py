#!/home/semplicio/virtualenv/bin/python -u

__author__		= "Giuseppe Semplicio"
__copyright__	= "Copyright 2016"
__version__		= "1.0"
__credits__		= "Giuesppe Semplicio"
__maintainer__	= "None"
__email__		= "semplicio@ie-freiburg.mpg.de"
__status__		= "Production"

# This was generated to obtain the first nucleotide of each annotated transcript.
# In general this script uses a bed format file (col1 chr, col2 start is essential) and then uses these
# coordinations to find the corresponding nucleotide in the fasta file.

import sys
import argparse
from collections import Counter
import os

def main():
	parser = argparse.ArgumentParser(description='Extracting the desired nucleotide from a given list of positions and a reference fasta file')
	parser.add_argument('-r', '--reference', metavar='abs_path_reference_fasta', help = 'The absolute path to the input file with the sequences. Must be in fasta format.', required=True)
	parser.add_argument('-i', '--input_bed', metavar='abs_path_bed', help='The absolute patch to the file containing the coordinates (chr, start). Column 1 must be chromosome, column 2 must be position.', required=True)
	parser.add_argument('-t', '--tmp_fa', action='store_true', help='Use a tmp.fa file that was created in a previous run.')
	args = parser.parse_args()

	current_dir = os.getcwd()+"/"
	nuc_extractor(args.reference, args.input_bed, current_dir, args.tmp_fa)
	

def nuc_extractor(reference, in_bed, current_dir, tmp_fa):

	if tmp_fa == False:

		fasta = open(reference, 'r')
		line = fasta.next().strip()

		if line.startswith(">"):
			with open(current_dir+'tmp.fa', 'a') as ff:
				ff.write('%s\n' %line)

		while True:
			try:
				line = fasta.next().strip()

				if line.startswith(">"):
					with open(current_dir+'tmp.fa', 'a') as ff:
						ff.write('\n%s\n' %line)
				else:
					with open(current_dir+'tmp.fa', 'a') as ff:
						ff.write('%s' %line)
			except StopIteration:
				break

	tmp = open(current_dir+'tmp.fa', 'r')
	
	fasta_dict = {}
	count = 1
	while True:
		try:
			line1 = tmp.next().strip()
			line2 = tmp.next().strip()

			fasta_dict[line1.split(">")[1].split(" ")[0]] = line2
		except StopIteration:
			break


	bed = open(in_bed, 'r')

	linecount = 0

	for f in bed:
		linecount += 1

	nuc_list = []
	line_check = bed.next().strip()
	chrom_check = line_check.split("\t")[0]
	chr_chrom_check = "chr"+chrom_check
	
	bed.close

	bed = open(in_bed, 'r')

	if fasta_dict.has_key(chrom_check):
		
		for i in bed:
			count += 1
			if count == (linecount*0.1):
				print "10%"
			if count == (linecount*0.2):
				print "20%"
			if count == (linecount*0.3):
				print "30%" 
			if count == (linecount*0.4):
				print "40%"
			if count == (linecount*0.5):
				print "50%"
			if count == (linecount*0.6):
				print "60%"
			if count == (linecount/0.7):
				print "70%"
			if count == (linecount/0.8):
				print "80%"
			if count == (linecount/0.9):
				print "90%"
			if count == (linecount):
				print "100%"

			chrom = i.split("\t")[0]
			pos = int((i.split("\t")[1]))-1

			nuc = fasta_dict[chrom][pos]
			nuc_list.append(nuc)
	
	elif fasta_dict.has_key(chr_chrom_check):
		
		for i in bed:
			count += 1
			if count == (linecount*0.1):
				print "10%"
			if count == (linecount*0.2):
				print "20%"
			if count == (linecount*0.3):
				print "30%" 
			if count == (linecount*0.4):
				print "40%"
			if count == (linecount*0.5):
				print "50%"
			if count == (linecount*0.6):
				print "60%"
			if count == (linecount/0.7):
				print "70%"
			if count == (linecount/0.8):
				print "80%"
			if count == (linecount/0.9):
				print "90%"
			if count == (linecount):
				print "100%"
			chrom = "chr"+i.split("\t")[0]
			pos = int((i.split("\t")[1]))-1

			nuc = fasta_dict[chrom][pos]
			nuc_list.append(nuc)

	else:
		chrom_check_chr = line_check.split("\t")[0].split("chr")[1]
		if fasta_dict.has_key(chrom_check_chr):
		
			for i in bed:
				count += 1
				if count == (linecount*0.1):
					print "10%"
				if count == (linecount*0.2):
					print "20%"
				if count == (linecount*0.3):
					print "30%" 
				if count == (linecount*0.4):
					print "40%"
				if count == (linecount*0.5):
					print "50%"
				if count == (linecount*0.6):
					print "60%"
				if count == (linecount/0.7):
					print "70%"
				if count == (linecount/0.8):
					print "80%"
				if count == (linecount/0.9):
					print "90%"
				if count == (linecount):
					print "100%"
				chrom = i.split("\t")[0].split("chr")[1]
				pos = int((i.split("\t")[1]))-1
				
				nuc = fasta_dict[chrom][pos]
				nuc_list.append(nuc)

	nuc_counted = Counter(nuc_list)
	with open(current_dir+"nuc_extractor_stats.log", 'a') as ff:
		ff.write("There were %s entries.\nFor the given positions:\n" %count)

	print "\nThere were %s entries in total." %count
	print "For the given positions:"
	for i in 'acgt':
		msg = "\t%s was found %d times" %(i, nuc_counted[i])
		print msg
		with open(current_dir+"nuc_extractor_stats.log", 'a') as ff:
			ff.write(msg)

	

if __name__ == "__main__":
	main()



for i in fastain:
	if i.startswith(">"):
		fasta_dict[i.split(">")[1].split(" ")[0]] = None
		tmp_chr_i = i.split(">")[1].split(" ")[0]
	else:
		fasta_list.append(i)
	