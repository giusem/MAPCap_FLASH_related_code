#!/home/semplicio/virtual_environments/ipython/bin/python -u

import argparse
import sys
import os
import Levenshtein
import time
import re
import gzip
import subprocess
import glob

def processFastqFile(forward, reverse, barcodes, chunk, file_type, gzipp, hamming_distance, current_dir, keep_bad_barcodes, filetype):
	
	sample_names, sample_index_barcodes = extract_barcodes(barcodes)

	index_uniq = {x for x in sample_index_barcodes}

	if len(index_uniq) != len(sample_index_barcodes):
		printStatus("Something is wrong with your barcodes! Are they all unique?")
		sys.exit(1)

	stats = fastqSplitter(forward, reverse, sample_names, sample_index_barcodes, chunk, hamming_distance, current_dir, keep_bad_barcodes, filetype, gzipp)


	print "\n%s\n" %forward
	print "%s" %reverse

	if stats[2] > 0:
		print "\n\n############################################################\n"
		print "Here's some barcode splitting statistics:\n"
		print "Total number of reads:\t\t %s" %(stats[0])
		print "Number of reads that contained \nan unidentifiable barcode:\t %s (%s%% of all reads)" %(stats[2], round(float(stats[2])/float(stats[0])*100))
		for i in range(len(sample_names)):
			print "Reads from sample %s: \t %s" %(sample_names[i], stats[1][i])
		print "############################################################\n\n"

		cmd_bad_barcodes = "sort %sbad_barcodes.txt | uniq -c | sort -k1,1nr > %sbad_barcodes_with_counts.txt" %(current_dir, current_dir)
		subprocess.check_call(cmd_bad_barcodes, shell=True)
		
		rm = "rm %sbad_barcodes.txt" %current_dir
		subprocess.check_call(rm, shell=True)

	else:
		print "\n\n############################################################\n"
		print "Here's some barcode splitting statistics:\n"
		print "Total number of reads:\t %s" %(stats[0])
		for i in range(len(sample_names)):
			print "Reads from sample %s: %s" %(sample_names[i], stats[1][i])
		print "############################################################\n\n"

	if gzipp:
		fastq_files = glob.glob(current_dir+"*.fastq")
		if len(fastq_files) > 0:
			cmd = "gzip %s*.fastq" %current_dir
			subprocess.check_call(cmd, shell=True)

def main():
	parser = argparse.ArgumentParser(description="Demultiplex samples with a lollipop adapter.")
	parser.add_argument('-f', '--forward', metavar='path/to/left_reads_R1.fastq.gz', help='The fastq file with the forward reads', required=True)
	parser.add_argument('-r', '--reverse', metavar='path/to/right_reads_R2.fastq.gz', help='The fastq file with the reverse reads', required=True)
	parser.add_argument('-b', '--barcodes', metavar='path/to/barcodes.fa', help='The absolute path to the fasta file that contains the barcodes.' ,required=True)
	parser.add_argument('-c', '--chunk', type=int, help='Determine the chunk size. The formula is chunk * 100000 per file. Default is 2.', default=2)
	parser.add_argument('-gz', '--gzip', action='store_true', help='Gzip the fastq files.')
	parser.add_argument('-hm', '--hamming_distance', type=int, help="Maximum allowed mutations in the internal index. Default is 1 mutation allowed. For shorter barcodes, can 0 (i.e. perfect match). For longer barcodes can use larger numbers. Experiment with this number for best results", default=1)
	parser.add_argument('-k', '--keep', action='store_true', help='Enabling this option will keep the sequences that could not be assigned to a given reference index barcode')
	args = parser.parse_args()
	
	current_dir = os.path.abspath(args.barcodes).split("barcodes")[0]

	if re.search('fastq', args.forward, re.IGNORECASE):
		if re.search('fastq.gz', args.forward, re.IGNORECASE):
			filetype = 'fastq.gz'
		else:
			filetype = 'fastq'
	else:
		print("I don't understand your file type. Your files should end either with .fastq or .fastq.gz")
		sys.exit(1)

	processFastqFile(args.forward, args.reverse, args.barcodes, args.chunk*100000, filetype, args.gzip, args.hamming_distance, current_dir, args.keep, filetype)

def extract_barcodes(barcodes):
	barcode_names = []
	barcode_seqs =[]
	with open(barcodes, 'r') as ff:
		while True:
			try:
				barcode_names.append(ff.next().strip().split('>')[1])
				barcode_seqs.append(ff.next().strip())
			except StopIteration:
				break
		return barcode_names, barcode_seqs



def fastqSplitter(forward, reverse, sample_names, sample_index_barcodes, chunk, hamming_distance, current_dir, keep_bad_barcodes, filetype, gzipp):
	
	if filetype == 'fastq.gz':
		left_reads = gzip.open(forward, 'rb')
		right_reads = gzip.open(reverse, 'rb')

	elif filetype == 'fastq':
		left_reads = open(forward, 'r')
		right_reads = open(reverse, 'r')

	multiplier = 1
	count = 0 #total number of reads
	alien_count = 0 # not A or B.

	f_reads_dict = {}
	r_reads_dict = {}
	if keep_bad_barcodes:
		alien_reads_f_dict = {}
		alien_reads_r_dict = {}
		bad_set = []
	sample_count = [0 for x in range(len(sample_names))]

	while True:
		try:
			name_f, name_r	= left_reads.next().strip(), right_reads.next().strip()
			seq_f,	seq_r	= left_reads.next().strip(), right_reads.next().strip()
			plus_f, plus_r	= left_reads.next().strip(), right_reads.next().strip()
			qual_f, qual_r	= left_reads.next().strip(), right_reads.next().strip()
			count += 1
			# collect the random barcode from both ends, wherever it is defined

			index_bcode = seq_r[5:11]

			scores = []

			for i in sample_index_barcodes:
				scores.append(Levenshtein.hamming(i, index_bcode))

			# If binary barcode scheme has been used, split accordingly, assuming that A or B refers to biological replicates or just replicates of some sort

			if min(scores) <= hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 and unique
				sample_count[scores.index(min(scores))] += 1
				try:
					f_reads_dict[current_dir + sample_names[scores.index(min(scores))] + "_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)])
					r_reads_dict[current_dir + sample_names[scores.index(min(scores))] + "_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_r, seq_r, qual_r)])
				except KeyError:
					f_reads_dict[current_dir + sample_names[scores.index(min(scores))] + "_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)]]
					r_reads_dict[current_dir + sample_names[scores.index(min(scores))] + "_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_r, seq_r, qual_r)]]

			elif keep_bad_barcodes:
				try:
					alien_reads_f_dict[current_dir+"unknownBarcode_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)])
					alien_reads_r_dict[current_dir+"unknownBarcode_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_r, seq_r, qual_r)])
				except KeyError:
					alien_reads_f_dict[current_dir+"unknownBarcode_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)]]
					alien_reads_r_dict[current_dir+"unknownBarcode_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_r, seq_r, qual_r)]]
				bad_set.append(index_bcode)
				alien_count += 1


			if count == multiplier*chunk:
				for key, value in f_reads_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])
				
				for key, value in r_reads_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				if keep_bad_barcodes:
					for key, value in alien_reads_f_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in alien_reads_r_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

				multiplier += 1
				f_reads_dict = {}
				r_reads_dict = {}
				if keep_bad_barcodes:
					alien_reads_f_dict = {}
					alien_reads_r_dict = {}

		except StopIteration:
			for key, value in f_reads_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			for key, value in r_reads_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			if keep_bad_barcodes:
				for key, value in alien_reads_f_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				for key, value in alien_reads_r_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])


				bad_file = current_dir+"bad_barcodes.txt"
				with open(bad_file, 'a') as ff:
					for i in bad_set:
						ff.write(i+'\n')

			f_reads_dict = {}
			r_reads_dict = {}
			if keep_bad_barcodes:
				alien_reads_f_dict = {}
				alien_reads_r_dict = {}
			break

	left_reads.close()
	right_reads.close()

	if keep_bad_barcodes:
		return count, sample_count, alien_count
	else:
		return count, sample_count, 0

def printStatus(msg):
	print "%s: %s" % (time.strftime("%X %x"), msg)


if __name__ == "__main__":
	start_time = time.time()
	printStatus("Initializing the script.")
	main()
	end_time = time.time() - start_time
	printStatus("Done. Elapsed time: %s seconds." %round(end_time))


