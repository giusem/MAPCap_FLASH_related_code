#!/home/semplicio/virtualenv/bin/python -u

# starting from the fastq read identifier of R1, find the the corresponding R2 read in fastq file and put it in a new fastq file

from __future__ import division
import time
import sys
import os
import re
import gzip
import argparse



def main():
	parser = argparse.ArgumentParser(description="This tool takes the uniquely identified reads from FLASH/MAPCAP and transforms themback to FASTQ files\n")
	parser.add_argument('-b', '--BED', metavar='path/to/BEDfile', help='The BED file containing in the 4th column the read identifier. If not available this can be done by using bedtools bamtobed.' ,required=True)
	parser.add_argument('-f', '--forward', metavar='path/to/forwardReads', help='The R1 file from which the reads should be extracted.')
	parser.add_argument('-r', '--reverse', metavar='path/to/reverseReads', help='The R2 file from which the reads should be extracted.')
	parser.add_argument('-m', '--merged', metavar='path/to/mergedFASTQfile', help='The merged FASTQ file to extract the reads from.')
	parser.add_argument('-o', '--out', metavar='name of output fastq file', help='The name of the output file to which the extracted reads will be added.' ,required=True)
	args = parser.parse_args()

	current_dir = os.getcwd()+"/"
	inputpath = os.path.realpath(args.BED)

	if args.forward and args.reverse:
		forwardpath = os.path.realpath(args.forward)
		reversepath = os.path.realpath(args.reverse)
		fwd_out = current_dir+args.out.split(".fastq")[0]+"_UNIQUE_R1.fastq"
		rvs_out = current_dir+args.out.split(".fastq")[0]+"_UNIQUE_R2.fastq"

		if re.search('bed', args.BED, re.IGNORECASE):
			print
		else:
			print("%s is not a BED file. Your files should end with '.bed'." %args.BED)
			sys.exit(1)

		if re.search('fastq', args.forward, re.IGNORECASE):
			if re.search('fastq.gz', args.forward, re.IGNORECASE):
				filetype1 = 'fastq.gz'
			else:
				filetype1 = 'fastq'
		else:
			print("I don't understand the file type of %s. Your files should end either with .fastq or .fastq.gz" %args.fastq)
			sys.exit(1)

		if re.search('fastq', args.reverse, re.IGNORECASE):
			if re.search('fastq.gz', args.reverse, re.IGNORECASE):
				filetype2 = 'fastq.gz'
			else:
				filetype2 = 'fastq'
		else:
			print("I don't understand the file type of %s. Your files should end either with .fastq or .fastq.gz" %args.fastq)
			sys.exit(1)

	if args.merged:
		mergedpath = os.path.realpath(args.merged)
		merged_out = current_dir+args.out.split(".fastq")[0]+"_UNIQUE_merged.fastq"

		if re.search('fastq', args.merged, re.IGNORECASE):
			if re.search('fastq.gz', args.merged, re.IGNORECASE):
				filetype3 = 'fastq.gz'
			else:
				filetype3 = 'fastq'
		else:
			print("I don't understand the file type of %s. Your files should end either with .fastq or .fastq.gz" %args.mergedFASTQ)
			sys.exit(1)	

	if args.merged and args.forward and args.reverse:		
		uniqeFASTQ(inputpath, forwardpath, reversepath, mergedpath, fwd_out, rvs_out, merged_out, current_dir, filetype1, filetype2, filetype3)

	elif args.merged and not args.forward and not args.reverse:
		uniqeFASTQ(inputpath, None, None, mergedpath, None, None, merged_out, current_dir, None, None, filetype3)
	
	elif args.forward and args.reverse and not args.merged:	
		uniqeFASTQ(inputpath, forwardpath, reversepath, None, fwd_out, rvs_out, None, current_dir, filetype1, filetype2, None)

def uniqeFASTQ(inputpath, forwardpath, reversepath, mergedpath, fwd_out, rvs_out, merged_out, current_dir, file_type1, file_type2, file_type3):

	linecount1 = 0
	linecount2 = 0

	inputbed = open(inputpath, 'r')

	if file_type1 == 'fastq.gz':
		with gzip.open(forwardpath, 'r') as ff:
			for f in ff:
				linecount1 += 1
		readcount1 = int(linecount1/4)
		forward = gzip.open(forwardpath, 'rb')
	elif file_type1 == 'fastq':
		with open(forwardpath, 'r') as ff:
			for f in ff:
				linecount1 += 1
		readcount1 = int(linecount1/4)
		forward = open(forwardpath, 'rb')

	if file_type2 == 'fastq.gz':
		reverse = gzip.open(reversepath, 'rb')
	elif file_type2 == 'fastq':
		reverse = open(reversepath, 'rb')

	if file_type3 == 'fastq.gz':
		with gzip.open(mergedpath, 'r') as ff:
			for f in ff:
				linecount2 += 1
		readcount2 = int(linecount2/4)
		merged = gzip.open(mergedpath, 'rb')
	elif file_type3 == 'fastq':
		with open(mergedpath, 'r') as ff:
			for f in ff:
				linecount2 += 1
		readcount2 = int(linecount2/4)
		merged = open(mergedpath, 'rb')

	# prepare the bed file by taking the IDs and throw them in a list
	ID_set = set()

	while True:
		try:
			line = inputbed.next().strip()
			ID = line.split("\t")[3]
			ID_set.add(ID)

		except StopIteration:
			break

	# compare the IDs from the list with the fastq files
	forward_set = set()
	reverse_set = set()
	merged_set = set()

	if forwardpath and reversepath:
		print "processing unmerged reads now"

		count = 0
		x = 1
		status = 10

		while True:
			try:
				name_f, name_r	= forward.next().strip(), reverse.next().strip()
				seq_f,	seq_r	= forward.next().strip(), reverse.next().strip()
				plus_f, plus_r	= forward.next().strip(), reverse.next().strip()
				qual_f, qual_r	= forward.next().strip(), reverse.next().strip()

				count += 1

				if count >= ((readcount1/10)*x):
					print str(status)+"%"+" of unmerged reads processed..."
					x += 1
					status += 10

				fastqID_f = name_f.split(" ")[0].split("@")[1]
				fastqID_r = name_r.split(" ")[0].split("@")[1]

				if fastqID_f in ID_set:
					forward_set.add("%s\n%s\n+\n%s\n" %(name_f, seq_f, qual_f))
					
				if fastqID_r in ID_set:
					reverse_set.add("%s\n%s\n+\n%s\n" %(name_r, seq_r, qual_r))

				
			except StopIteration:
				
				out_name = open(fwd_out, 'a')
				for i in forward_set:
					out_name.write(i)
				out_name.close()

				out_name = open(rvs_out, 'a')
				for i in reverse_set:
					out_name.write(i)
				out_name.close()

				break
	if mergedpath:
		print "processing the merged reads now"

		count = 0
		x = 1
		status = 10

		while True:
			try:
				name_m = merged.next().strip()
				seq_m = merged.next().strip()
				plus_m = merged.next().strip()
				qual_m = merged.next().strip()

				count += 1

				if count >= ((readcount2/10)*x):
					print str(status)+"%"+" of merged reads processed..."
					x += 1
					status += 10

				fastqID_m = name_m.split(" ")[0].split("@")[1]

				if fastqID_m in ID_set:
					merged_set.add("%s\n%s\n+\n%s\n" %(name_m, seq_m, qual_m))
			
			except StopIteration:

				out_name = open(merged_out, 'a')
				for i in merged_set:
					out_name.write(i)
				out_name.close()

				break


if __name__ == "__main__":
 	print "\nInitialized the script at: %s" %(time.strftime("%X"))
 	main()
 	print "\nDone. Now it's: %s" %(time.strftime("%X"))
	
