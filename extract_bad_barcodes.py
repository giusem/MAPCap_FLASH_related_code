#!/home/semplicio/virtualenv/bin/python -u


# extracts sequence of desired lentgh and position from fastq file, counts the sequences 
# allows to save reads with desired sequence in fastq format



import re
import subprocess
import gzip
import sys
import time
import argparse

def main():
	parser = argparse.ArgumentParser(description='Extract sequence from fastq file, save underlying read of desired sequence.')
	parser.add_argument('-i', '--input', metavar='path_to_input_fastq', help = 'Path to input file. Can be fastq or fastq.gz', required=True)
	parser.add_argument('-s', '--start', type=int, help='Start position of the sequence', required=True)
	parser.add_argument('-e', '--end', type=int, help='End position of the sequence', required=True)
	parser.add_argument('-o', '--orientation', choices=['forward', 'reverse'], help='Whether the read is from R1 or R2. Reads in reverse orientation (R2) will be reverse complement to match forward reads.')
	parser.add_argument('-bc', '--barcodes', type=str, nargs='*', help='Reads with given barcode sequence, in CAPITALS, will be extracted and save as separate fastq file. Multiple sequences have to be space separated. Example: AGATGC NNNNNN')
	args = parser.parse_args()

	#current_dir = os.getcwd()+"/"
	#nuc_extractor(args.reference, args.input_bed, current_dir, args.tmp_fa)

	print args.barcodes

	barcode_extractor(args.input, args.start, args.end, args.orientation, args.barcodes)

def barcode_extractor(fastq_input, start, end, orientation, barcodes):

	start_pos = start - 1
	end_pos = end - 1

	#if barcodes:
	#	ind_barcode = barcode.split()


	if re.search('fastq.gz', fastq_filename, re.IGNORECASE):
		filetype = 'fastq.gz'
	else:
		filetype = 'fastq'

	if filetype == 'fastq.gz':
		fastq_input = gzip.open(fastq_filename, 'rb')

	elif filetype == 'fastq':
		fastq_input = open(fastq_filename, 'rb')
	else:
		print "the filetype is wrong...quitting"
		sys.exit(1)

	#extract_targets = extract.split()

	extract_dict = {}
	bad_set = []

	while True:
		try:
			name = fastq_input.next().strip()
			seq = fastq_input.next().strip()
			plus = fastq_input.next().strip()
			qual = fastq_input.next().strip()

			if orientation == 'forward':
				index_bcode = seq[start_pos:end_pos]
			elif orientation == 'reverse':
				complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N' : 'N'}
				index_bcode = "".join(complement.get(base, base) for base in reversed(seq))[start_pos:end_pos]
			else:
				print "don't know the orientation of the read...quitting"
				sys.exit(1)

			bad_set.append(index_bcode)


			if index_bcode in barcodes:
				try:
					extract_dict["bad_barcodes_"+index_bcode+"_reads.txt"].append(["%s\n%s\n+\n%s\n"%(name, seq, qual)])
				except KeyError:	
					extract_dict["bad_barcodes_"+index_bcode+"_reads.txt"] = ["%s\n" %seq]


		except StopIteration:
			break



	bad_file = "bad_barcodes.txt"					
	with open(bad_file, 'a') as ff:
		for i in bad_set:
			ff.write(i+'\n')


	cmd = "sort bad_barcodes.txt | uniq -c | sort -k1,1n > bad_barcodes_with_counts.txt"
	subprocess.check_call(cmd, shell=True)


	for key, value in extract_dict.iteritems():
		with open('%s' %key, 'a+') as gg:
			for i in range(len(value)):
				gg.write('%s' %value[i][0])

def printStatus(msg):
	print "%s: %s" % (time.strftime("%X %x"), msg)

if __name__ == "__main__":
	start_time = time.time()
	printStatus("Initializing the script.")
	main()
	end_time = time.time() - start_time
	printStatus("Done. Elapsed time: %s seconds." %round(end_time))