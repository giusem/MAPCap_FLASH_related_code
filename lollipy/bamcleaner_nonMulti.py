#!/home/semplicio/virtualenv/bin/python -u
#todo: work with paired-end data
from __future__ import division

__author__		= "Ibrahim Ilik"
__copyright__	= "Copyleft 2015"
__version__		= "0.0.2"
__credits__		= ["Ibrahim Ilik"]
__maintainer__	= "None"
__email__		= "ibrahimus@gmail.com"

import Levenshtein
import time
import argparse
import hashlib
import subprocess
import os
import random
import multiprocessing as mp
import glob
import sys
import pysam

def DupRemover(inputbam, output, bedout, chunk, quality, multi, current_dir, usage):
	bedtools_path = "/package/bedtools2/bin/bedtools"
	samtools_path = "/package/samtools/bin/samtools"
	start = time.time()
	#pre-processing. Generating unique names for temporary files, in case multiple instances of the script is run at the same time.

	mm = hashlib.sha1()
	mm.update(inputbam.split("/")[-1]+str(random.random()))
	tempbed = current_dir+".temp_bed-"+mm.hexdigest()
	tempbed0 = current_dir+".0temp_bed-"+mm.hexdigest()
	tempsam = current_dir+".temp_sam-"+mm.hexdigest()
	ticker = False
	ticker2 = True
	logfile = output.split(".bam")[0]+".log"

	with open(logfile, 'w') as ff:
		msg = "Working on file %s\n" %inputbam.split("/")[-1]
		ff.write("%s: %s" % (time.strftime("%X %x"), msg))
		
	output = output.split(".bam")[0]+".sam"

	# check if the input file is from bbmap or bowtie

	cmd0 = "%s bamtobed -i %s | head > %s" %(bedtools_path, inputbam, tempbed0)
	subprocess.check_call(cmd0, shell=True)
	x = os.path.abspath(tempbed0)
	xx = open(tempbed0, 'r')
	# This routine tries to understand if the input file came from bowtie2 or bbmap (because bbmap generate some funny .sam files). 
	# It may not work for other aligners. Although I guess it might also work.

	
	if os.path.getsize(x) > 0:
		check_file = inputbam.split("/")[-1].split(".")[0]+".lock"
		with open(check_file, 'a') as ff:
			pass

		if len(xx.next().strip().split(" ")) == 1:
			cmd = "%s bamtobed -i %s | sed 's/:/ /g' | sort -k1,1 -k2,2n -k3,3n -k11,11 -k12,12n -k13,13 > %s" %(bedtools_path, inputbam, tempbed)
		elif len(xx.next().strip().split(" ")) == 2:
			ticker2 = False # marks  as "bbmap" output
			cmd = "%s bamtobed -i %s | awk '{print $1,$2,$3,$4,$6,$7}' | tr [:blank:] \\\\t | sed 's/:/ /g' | sort -k1,1 -k2,2n -k3,3n -k11,11 -k12,12n -k13,13 > %s" %(bedtools_path, inputbam, tempbed)
		elif len(xx.next().strip().split(" ")) == 3:
			ticker2 = False#  marks  as "bbmap" output
			cmd = "%s bamtobed -i %s | awk '{print $1,$3,$4,$5,$7,$8}' | tr [:blank:] \\\\t | sed 's/:/ /g' | sort -k1,1 -k2,2n -k3,3n -k11,11 -k12,12n -k13,13 > %s" %(bedtools_path, inputbam, tempbed)
		elif len(xx.next().strip().split(" ")) == 4:
			ticker2 = False#  marks  as "bbmap" output
			cmd = "%s bamtobed -i %s | awk '{print $1,$4,$5,$6,$8,$9}' | tr [:blank:] \\\\t | sed 's/:/ /g' | sort -k1,1 -k2,2n -k3,3n -k11,11 -k12,12n -k13,13 > %s" %(bedtools_path, inputbam, tempbed)
		else:
			print "Unrecognized .bam formatting."
			sys.exit(1)

		xx.close()
		os.remove(tempbed0)


		with open(logfile, 'a') as ff:
			msg = "Generating some temp files.\n"
			ff.write("%s: %s" % (time.strftime("%X %x"), msg))
		
		cmd2 = "%s view -h %s > %s" %(samtools_path, inputbam, tempsam)
		
		with open(logfile, 'a') as ff:
			msg = "Running command: \n\t%s\n"%cmd
			ff.write("%s: %s" % (time.strftime("%X %x"), msg))
		
		subprocess.check_call(cmd, shell=True)
		
		with open(logfile, 'a') as ff:
			msg = "Running command: \n\t%s"%cmd2
			ff.write("%s: %s\n" % (time.strftime("%X %x"), msg))	
		 
		subprocess.check_call(cmd2, shell=True)

		with open(logfile, 'a') as ff:
			msg = "Done. Starting to work on the temp files.\n"
			ff.write("%s: %s" % (time.strftime("%X %x"), msg))	
		
		aa = open(tempbed, 'r')

		seq_set = set()
		name_set = set()
		sam_set = set()
		counter = 0 # reads that are kept
		counter2 = 2 # total number of reads
		counter3 = 0 # reads that are below the quality threshold. 
		counter4 = 0 # PCR dups
		counter5 = 0 # PCRdups and below quality threshold
		multiplier = 1

		# Initiate the whole thing with the first two lines.
		line1 = aa.next()
		line2 = aa.next()

		line1_bcode = line1.split(" ")[-1].split("\t")[0]
		line2_bcode = line2.split(" ")[-1].split("\t")[0]
		line1_qual = int(line1.split("\t")[-2])
		line1_strand = line1.strip().split("\t")[-1]
		line2_strand = line2.strip().split("\t")[-1]
		line1_begin = line1.split("\t")[0] + "," + line1.split("\t")[1]
		line1_end = line1.split("\t")[0] + "," + line1.split("\t")[2]
		line2_begin = line2.split("\t")[0] + "," + line2.split("\t")[1]
		line2_end = line2.split("\t")[0] + "," + line2.split("\t")[2]
		levent = Levenshtein.hamming(line1_bcode, line2_bcode)
		
		#comparing consecutive lines with each other. 

		if line1_strand == line2_strand:
			if line1_strand == "+":
				if line1_begin != line2_begin:
					if line1_qual >= quality:
						if bedout:
							seq_set.add(":".join(line1.split(" ")))
						name_set.add(":".join(line1.split(" ")).split("\t")[3])
						counter += 1
					else:
						counter3 +=1
				else:	
					if levent > 1:
						if line1_qual >= quality:
							if bedout:
								seq_set.add(":".join(line1.split(" ")))
							name_set.add(":".join(line1.split(" ")).split("\t")[3])
							counter += 1
						else:
							counter3 += 1
					else: #if barcodes not different
						if line1_qual >= quality: #if quality ok, then pcr duplicate
							counter4 += 1
						else: #if quality bad, then pcr duplicate of multimapper
							counter5 += 1

			elif line1_strand =="-":
				if line1_end != line2_end:
					if line1_qual >= quality:
						if bedout:
							seq_set.add(":".join(line1.split(" ")))
						name_set.add(":".join(line1.split(" ")).split("\t")[3])
						counter += 1
					else:
						counter3 += 1
				else:	
					if levent > 1:
						if line1_qual >= quality:
							if bedout:
								seq_set.add(":".join(line1.split(" ")))
							name_set.add(":".join(line1.split(" ")).split("\t")[3])
							counter += 1
						else:
							counter3 += 1
					else: #if barcodes not different
						if line1_qual >= quality: #if quality ok, then pcr duplicate
							counter4 += 1
						else: #if quality bad, then pcr duplicate of multimapper
							counter5 += 1
			else:
				printStatus("The strand is neither + or - Something must be wrong, exiting")
				aa.close()
				os.remove(tempbed)
				os.remove(tempsam)
				sys.exit(1)
		else:		
			if line1_qual >= quality:
				if bedout:
					seq_set.add(":".join(line1.split(" ")))
			if line1_qual >= quality:
				name_set.add(":".join(line1.split(" ")).split("\t")[3])
				counter += 1
			else:
				counter3 += 1

		#switching lines, so we can iteratively work with consecutive lines
			
		line1 = line2

		while True:
			try:
				line1_bcode = line1.split(" ")[-1].split("\t")[0]
				line1_qual = int(line1.split("\t")[-2])
				line2 = aa.next()
				counter2 += 1 #this counter holds the total number of lines in the file. 
				line2_bcode = line2.split(" ")[-1].split("\t")[0]
				line1_strand = line1.strip().split("\t")[-1]
				line2_strand = line2.strip().split("\t")[-1]
				line1_begin = line1.split("\t")[0] + "," + line1.split("\t")[1]
				line1_end = line1.split("\t")[0] + "," + line1.split("\t")[2]
				line2_begin = line2.split("\t")[0] + "," + line2.split("\t")[1]
				line2_end = line2.split("\t")[0] + "," + line2.split("\t")[2]			
				levent = Levenshtein.hamming(line1_bcode, line2_bcode)
				
				if line1_strand == line2_strand:
					if line1_strand == "+":
						if line1_begin != line2_begin:
							if line1_qual >= quality:
								if bedout:
									seq_set.add(":".join(line1.split(" ")))
								name_set.add(":".join(line1.split(" ")).split("\t")[3])
								counter += 1
							else:
								counter3 +=1
						else:	
							if levent > 1:
								if line1_qual >= quality:
									if bedout:
										seq_set.add(":".join(line1.split(" ")))
									name_set.add(":".join(line1.split(" ")).split("\t")[3])
									counter += 1
								else:
									counter3 += 1
							else: #if barcodes not different
								if line1_qual >= quality: #if quality ok, then pcr duplicate
									counter4 += 1
								else: #if quality bad, then pcr duplicate of multimapper
									counter5 += 1

					elif line1_strand =="-":
						if line1_end != line2_end:
							if line1_qual >= quality:
								if bedout:
									seq_set.add(":".join(line1.split(" ")))
								name_set.add(":".join(line1.split(" ")).split("\t")[3])
								counter += 1
							else:
								counter3 += 1
						else:	
							if levent > 1:
								if line1_qual >= quality:
									if bedout:
										seq_set.add(":".join(line1.split(" ")))
									name_set.add(":".join(line1.split(" ")).split("\t")[3])
									counter += 1
								else:
									counter3 += 1
							else: #if barcodes not different
								if line1_qual >= quality: #if quality ok, then pcr duplicate
									counter4 += 1
								else: #if quality bad, then pcr duplicate of multimapper
									counter5 += 1
									
					else:
						printStatus("The strand is neither + or - Something must be wrong, exiting.")
						aa.close()
						os.remove(tempbed)
						os.remove(tempsam)
						sys.exit(1)
				else:		
					if line1_qual >= quality:
						if bedout:
							seq_set.add(":".join(line1.split(" ")))
					if line1_qual >= quality:
						name_set.add(":".join(line1.split(" ")).split("\t")[3])
						counter += 1
					else:
						counter3 += 1

				line1 = line2

				# Every 100000 unique reads (or multiples of 100000 as determined by the chunk variable), the set that holds unique reads are flushed. 

				if bedout:
					if counter == multiplier * chunk:
						for i in seq_set:
							with open(bedout, 'a') as ff:
								ff.write(i)
						multiplier += 1
						seq_set = set()

			except StopIteration:
				if bedout:
					for i in seq_set:
						with open(bedout, 'a') as ff:
							ff.write(i)

				allremovedreads = int(counter3+counter4+counter5)

				with open(logfile, 'a') as ff:
					msg ="\nTotal number of reads:\t\t%s\t(%%%s)\n\nNumber of removed reads:\t%s\t(%%%s)\n\tPCR duplicates:\t\t%s\t(%%%s)\n\tBelow quality of %s:\t%s\t(%%%s)\n\tPCR duplicates and\n\tbelow quality of %s:\t%s\t(%%%s)\n\nKept reads:\t\t\t%s\t(%%%s)\n" %(counter2, 
						round((counter2/counter2)*100), allremovedreads, round((1-((counter2-allremovedreads)/counter2))*100), counter4, round((1-((counter2-counter4)/counter2))*100), quality, counter3, round((1-((counter2-counter3)/counter2))*100),
						quality, counter5, round((1-((counter2-counter5)/counter2))*100), counter, round((1-((counter2-counter)/counter2))*100))
					ff.write("%s: %s" % (time.strftime("%X %x"), msg))
				break
		aa.close()

		bb = open(tempsam, 'r')
		list_of_stuff = ['@SQ', '@PG', '@HD']
		counter3 = 0
		multiplier2 = 1
		

		# Capturing the header from the .bam file.

		while True:
			line = bb.next()
			if line.split("\t")[0] in list_of_stuff:
				with open(output, 'a') as ff:
					ff.write(line)
			else:
				if line.split("\t")[0] in name_set:
					sam_set.add(line)
					counter3 += 1
				break


		while True:
			try:
				sam_line = bb.next()
				if ticker2: # If input was from bowtie
					if sam_line.split("\t")[0] in name_set:
						sam_set.add(sam_line)
						counter3 += 1 #this is a useless counter. Should be the same as counter
				else: #if input was from BBMap
					if sam_line.split(" ")[0] in name_set:
						sam_set.add(sam_line.split(" ")[0]+"\t"+"\t".join(sam_line.split("\t")[1:]))
						counter3 += 1 #this is a useless counter. Should be the same as counter				

				if counter3 == multiplier2*chunk:
					output_file = open(output, 'a')
					for i in sam_set:
						output_file.write(i)
					output_file.close()
					
					multiplier2 += 1
					sam_set = set()

			except StopIteration:
				output_file = open(output, 'a')
				for i in sam_set:
					output_file.write(i)
				output_file.close()

				sam_set = set()
				break
		bb.close()

	#cleaning up

		os.remove(tempbed)
		os.remove(tempsam)

		
		cmd_bam = "%s view -Sb %s | samtools sort - %s" %(samtools_path, output, output.split(".sam")[0])
		subprocess.check_call(cmd_bam, shell=True)
		
		cmd_sam = "%s index %s" %(samtools_path, output.split(".sam")[0]+".bam")
		subprocess.check_call(cmd_sam, shell=True)
		os.remove(output)
		os.remove(check_file)

		if usage == "count_spikes":
			with open(logfile, 'a') as ff:
				msg = "Here are some counts on the spikes.\nThe total number of reads mapping to spikes is: %s.\nThey are distributed as follows:\n" %(counter)
				ff.write(msg)
				
			cmd_linecount = "%s view -c %s > linecount.tmp" %(samtools_path, output.split(".sam")[0]+".bam")	
			subprocess.check_call(cmd_linecount, shell=True)
			linecount_file = open("linecount.tmp", 'r')
			linecount = int(linecount_file.next())
			spike_list = []
			bam_file = pysam.AlignmentFile(output.split(".sam")[0]+".bam")

			SQ = bam_file.header['SQ']

			for i in SQ:
				spike_list.append(i['SN'])

			count = 0
			
			for i in spike_list:
				for j in bam_file.fetch(i):
					count+=1
				with open(logfile, 'a') as ff:
					msg = "%s \t %s \t %s \n" %(i, count, round((count/linecount)*100))
					ff.write(msg)
				count = 0


		
	else:				
		
		with open(logfile, 'a') as ff:
			msg = "%s is empty. Continue...\n" %inputbam.split("/")[-1]
			ff.write("%s: %s" % (time.strftime("%X %x"), msg))

		xx.close()
		os.remove(tempbed0)

	
	end = time.time()
	elapsed = round(end-start)
	with open(logfile, 'a') as ff:
			msg ="Done! Elapsed time: %s seconds.\n" %elapsed
			ff.write("%s: %s" % (time.strftime("%X %x"), msg))
			
def main():
	parser = argparse.ArgumentParser(description="Remove PCR-duplicates by making use of random barcodes in uvCLAP data. \
		The header of your fastq file should look like this: @39V34V1:117:H9FK3ADXX:1:2205:20357:12179:AACCCGCCAA 1:N:0:GGTAGC where AACCCGCCAA is the random barcode (of any length)")
	parser.add_argument('-i', '--input', help='The input .bam file. Doesn\'t have to be ordered.', required=True)
	parser.add_argument('-r', '--reverse-reads', help='For paired-end data. (not done yet')
	parser.add_argument('-k', '--keep-bed', help='Keep the cleaned .bed file? If so where?')
	parser.add_argument('-o', '--output', help='Output bam file.', required=True)
	parser.add_argument('-c', '--chunk', type=int, help='Chunk size. Number of reads kept in memory before written to disk. This value will be multiplied by 100000. Default = 1', default = 1)
	parser.add_argument('-m', '--min-quality', type=int, help='Minimum mapping quality. Default is 10', default=10)
	args = parser.parse_args()
	current_dir = "/".join(os.path.abspath(args.input).split("/")[0:-1])+"/"
	
	if args.input and args.output == 'all':
		print "Working on all bam files in this folder. You can find the the cleaned .bam, the index and the log file under the 'cleaned' folder."
		global multi
		multi = True
		if not os.path.exists("cleaned"):
				os.makedirs("cleaned")
		files = glob.glob("*.bam")
		jobs = []
		for i in files:
			p = mp.Process(target=DupRemover, args=(i, "cleaned/"+i.split(".bam")[0]+"_cleaned.bam", args.keep_bed, args.chunk*100000, args.min_quality, multi, current_dir))
			jobs.append(p)
			p.start()
	elif args.input and args.output == 'all':
		print "Both -i and -o should be 'all' for this to work."
		sys.exit(1)
	else:
		multi = False
		if args.output.split(".")[-1] != 'bam':
			print "The output has to be .bam!"
			print args.output
			print args.output.split(".")[-1]
			sys.exit(1)
		DupRemover(inputbam=args.input,  output=args.output, bedout=args.keep_bed, chunk=args.chunk, quality=args.min_quality, multi=multi, current_dir=current_dir, usage=None)

def printStatus(msg):
	print "%s: %s" % (time.strftime("%X %x"), msg)

if __name__ == "__main__":
	multi = False
	start1 = time.time()
	main()
	if not multi:
		end = time.time()
		elapsed = round(end-start1)
		printStatus("Done! Elapsed time: %s seconds." %elapsed)
