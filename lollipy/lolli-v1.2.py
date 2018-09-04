#!/home/semplicio/virtual_environments/ipython/bin/python -u
# change the path to your python interpreter

#changelog:

#lolli_v1.2
#changed p5end to p5end (=5 prime end)
#added p3end (=3 prime end) coverage track generation



#v.3.5 aka lolli_v1.1
#adjusted the hisat settings (--pen-cansplice 3)
#implemented cpm normalization for FLASH and MAPCap bigwigs on top of the raw bigwigs (replaces RPKM normalization for MAPCap)
#also all coverage files are now generated with bamCoverage (no bedgraphs anymore)
#removed paired-end mode which wasn't ready yet and fetish mode which is now obsolete




#v3.4 aka lollipy
# added inner_distance.py to calculate fragment length
# fixed hisat2 outputting multiple alignments by adding -F 256 option to samtools view
# added bbmap dm6_simple index (switches to Giuseppes bbmap for that)
#a bit of restructuring
#folders are only created if needed, bbmerge output name changed and is now in a folder
#some more stuff not worth mentioning
# added an option to map the unmerged reads as paired end (only with hisat2 so far)

#v3.3
#introduced "FETISH-mode", which means Hisat mapping with standard --pen-noncansplice option and CPM normalized coverage output and extended xlinking sites
#giving now the possibilty to keep the multimappers and PCR duplicates (option -k)

#v3.2
#improvements in merged and umerged splitter,  output file not opened every time but only once per chunk dump

#v3.1
#deeptools2 plotEnrichment implemented (count reads over features (std gtf file so far))


#v3
#SGE compatible
#uknown binary reads output
#bad barcode reads output
#minor changes in log output file
#serial run of dupremover

#v2
#spike mapper
#bad_barcodes now in bad_barcodes folder
#dm6_simple index (chr2,3,4,X,Y,rDNA)

#v1
#added bowtie2 only mapping


from __future__ import division
__author__		= ["Ibrahim Ilik, Giuseppe Semplicio"]
__copyright__	= "Copyright 2018"
__version__		= "lollipy v1.2"
__credits__		= ["Ibrahim Ilik, Giuseppe Semplicio"]
__maintainer__	= "None"
__email__		= "ibrahimus@gmail.com"
__status__		= "Production"


import re
import string
import time
import sys
import time
import gzip
import Levenshtein
import argparse
import os
import subprocess
import datetime
import glob
import multiprocessing as mp
import itertools
import pysam
#import numpy as np
#import matplotlib.pyplot as plt
try:
	from bamcleaner_nonMulti import DupRemover as duprem
except ImportError:
	print "bamcleaner.py not found. Stopping."
	sys.exit(1)

try:
	from bamcleaner_nonMulti_splitter import DupRemover as duprem_splitter
except ImportError:
	print "bamcleaner_splitter.py not found. Stopping."
	sys.exit(1)

script_path = os.path.dirname(os.path.abspath(__file__)) + "/"

#### PATHS ####################################################################
## Path defaults to all needed scripts and programs
## updated 13.3.2018

bowtie2_path = "/package/bowtie2-2.2.3/bin/bowtie2"; bowtie2_path_version = "2.2.3"
BBMAP_PATH = "/data/akhtar/group/Ibrahim_data_2014/__BBMap__/bbmap"
bbmap_path_sh = "/data/akhtar/group/Ibrahim_data_2014/__BBMap__/bbmap/bbmap.sh"; BBMAP_PATH_version = "37.56"
bbmerge_path = "/data/akhtar/group/Ibrahim_data_2014/__BBMap__/bbmap/bbmerge.sh"
hisat2_path = "/package/hisat2-2.0.4/hisat2"; hisat2_path_version = "2.0.4"
samtools_path = "/package/samtools-1.1/bin/samtools"; samtools_path_version = "1.1"
#genomeCoverageBed_path = "/package/bedtools2-2.19.0/bin/genomeCoverageBed"; genomeCoverageBed_path_version = "2.19.0"
#bedGraphToBigWig_path = "/home/semplicio/programs/bedGraphToBigWig/bedGraphToBigWig"; bedGraphToBigWig_path_version = "4"
#bedtools_path = "/package/bedtools2-2.19.0/bin/bedtools"; bedtools_path_version = "2.19.0"
plotEnrichment_path = "/package/anaconda3/envs/deeptools-3.0.1/bin/plotEnrichment"; plotEnrichment_path_version = "3.0.1"
bamCoverage_path = "/package/anaconda3/envs/deeptools-3.0.1/bin/bamCoverage"; bamCoverage_path_version = "3.0.1"
inner_distance_path = "/package/RSeQC-2.4/bin/inner_distance.py"; inner_distance_path_version = "2.4"

path_list = [bowtie2_path, bbmap_path_sh, bbmerge_path, hisat2_path, samtools_path, plotEnrichment_path, bamCoverage_path, inner_distance_path]

version_list = [bowtie2_path_version, BBMAP_PATH_version, BBMAP_PATH_version, hisat2_path_version, samtools_path_version, plotEnrichment_path_version, bamCoverage_path_version, bamCoverage_path_version, inner_distance_path_version]

for i in path_list:
	if not os.path.exists(i):
		printStatus("Error! Could not find %s. Maybe the path has changed. Exiting now..." %i)
		sys.exit(1)

#emergency break
#print "this version is not ready for usage yet, but improvements are on the way, have patience, use the older version"
#sys.exit(1)


def processFastqFiles(forward, reverse, merged, barcodes, file_type, bbmerge, ignore, chunk, mapping, gzipp, reps_merged, script_process, gzipafter, user_barcode, hamming_distance, processor, aligner, spikes_mapper, logfile, current_dir, keep_BAM, raw_fwd, raw_rvs, dist):
	# creating folders, extracting and defining barcodes ##

	index_barcode_names, index_barcodes = extract_barcodes(barcodes)

	index_uniq = {x for x in index_barcodes}

	if len(index_uniq) != len(index_barcodes):
		printStatus("Something is wrong with your barcodes! Are they all unique?")
		sys.exit(1)


	# Fill these in as reqired. Should also change the available options down at the main() function, --mapping option.

	indices = {'genome':['path/to/bowtie2_index','bbmap_index_designator_as_build_number','path/to/chrom.sizes', 'path/to/hisat_index', 'path/to/spike_index', 'path/to/gtf', 'path/to/bed'],
				'hg19':['/data/repository/organisms/GRCh37_ensembl/BowtieIndex/genome', 1, '/data/repository/organisms/GRCh37_ensembl/genome_fasta/genome.chrom.sizes',
				'/data/repository/organisms/GRCh37_ensembl/HISAT2Index/genome', None,
				'/data/repository/organisms/GRCh37_ensembl/gencode/release_19/genes.gtf', '/data/repository/organisms/GRCh37_ensembl/gencode/release_19/genes.bed'],
				'dm3':['/data/repository/organisms/dm3_ensembl/BowtieIndex/genome', 2, '/data/repository/organisms/dm3_ensembl/genome_fasta/genome.chrom.sizes',
				'/data/repository/organisms/dm3_ensembl/HISAT2Index/genome', '/data/akhtar/group/Giuseppe/INDEXES/bowtie2indexes/',
				'/data/repository/organisms/dm3_ensembl/Ensembl/release-78/genes.gtf', '/data/repository/organisms/dm3_ensembl/Ensembl/release-78/genes.bed'],
				'dm6':['/data/repository/organisms/dm6_ensembl/BowtieIndex/genome', 6, '/data/repository/organisms/dm6_ensembl/genome_fasta/genome.chrom.sizes',
				'/data/repository/organisms/dm6_ensembl/HISAT2Index/genome', '/data/akhtar/group/Giuseppe/INDEXES/bowtie2indexes/',
				'/data/repository/organisms/dm6_ensembl/Ensembl/release-79/genes.gtf', '/data/repository/organisms/dm6_ensembl/Ensembl/release-79/genes.bed'],
				'mm9':['/data/repository/organisms/GRCm37_ensembl/BowtieIndex/genome', 3, '/data/repository/organisms/GRCm37_ensembl/genome_fasta/genome.chrom.sizes',
				'/data/repository/organisms/GRCm37_ensembl/HISAT2Index/genome', None,
				'/data/repository/organisms/GRCm37_ensembl/gencode/m1/genes.gtf', '/data/repository/organisms/GRCm37_ensembl/gencode/m1/genes.bed'],
				'mm10':['/data/repository/organisms/GRCm38_ensembl/BowtieIndex/genome', 4, '/data/repository/organisms/GRCm38_ensembl/genome_fasta/genome.chrom.sizes',
				'/data/repository/organisms/GRCm38_ensembl/HISAT2Index/genome', None,
				'/data/repository/organisms/GRCm38_ensembl/gencode/m4/genes.gtf', '/data/repository/organisms/GRCm38_ensembl/gencode/m4/genes.bed'],
				'hg19_simple':['/data/akhtar/group/ibonomics/hg19_simple/hg19_simple', 5, '/data/akhtar/group/ibonomics/hg19_simple/hg19_simple.chrom.sizes',
				None, None, '/data/repository/organisms/GRCh37_ensembl/gencode/release_19/genes.gtf',
				'/data/repository/organisms/GRCh37_ensembl/gencode/release_19/genes.bed'],
				'hg38':['/data/repository/organisms/GRCh38_ensembl/BowtieIndex/genome', 7, '/data/repository/organisms/GRCh38_ensembl/genome_fasta/genome.chrom.sizes',
				None, None, '/data/repository/organisms/GRCh38_ensembl/gencode/release_21/genes.gtf',
				'/data/repository/organisms/GRCh38_ensembl/gencode/release_21/genes.bed'],
				'dm6_simple':['/data/akhtar/group/Giuseppe/INDEXES/bowtie2indexes/dm6_simple/dm6_simple', 5, '/data/akhtar/group/Giuseppe/INDEXES/HISAT2index/dm6_simple/dm6_simple.chr.sizes',
				'/data/akhtar/group/Giuseppe/INDEXES/HISAT2index/dm6_simple/dm6_simple', '/data/akhtar/group/Giuseppe/INDEXES/bowtie2indexes/',
				'/data/repository/organisms/dm6_ensembl/Ensembl/release-79/genes.gtf', '/data/repository/organisms/dm6_ensembl/Ensembl/release-79/genes.bed' ],
				'dvir':['/data/repository/organisms/dvir1.3_ensembl/BowtieIndex/genome ', 6, '/data/repository/organisms/dvir1.3_ensembl/genome_fasta/genome.chrom.sizes',
				'/data/repository/organisms/dvir1.3_ensembl/HISAT2Index/genome', None, '/data/repository/organisms/dvir1.3_ensembl/Ensembl/release-30/genes.gtf',
				'/data/repository/organisms/dvir1.3_ensembl/Ensembl/release-30/genes.bed']}

	### iterating over the input files and splitting reads according to their barcodes and read lengths ###

	merged_statistics = []
	unmerged_statistics = []

	#Figure out structure of the internal barcode, and start the splitter accordingly

	if user_barcode.lower() == 'iclip':
		left_user_barcode, right_user_barcode = 'NNNXXXXNN', None
	elif user_barcode.lower() == 'uvclap':
		left_user_barcode, right_user_barcode = 'NNNXXXXXNN', 'YRRYN'
	elif user_barcode.lower() == 'flash':
		left_user_barcode, right_user_barcode = None, 'NNRRNXXXXXXNN'
	else:
		left_user_barcode = None
		right_user_barcode = None
		user_barcode_list = user_barcode.lower().split(',')
		for i in user_barcode_list:
			if re.search('l=', i):
				left_user_barcode = i.split('l=')[1]
			elif re.search('r=', i):
				right_user_barcode = i.split('r=')[1]

	left_params = ([0,0], [0,0], [0,0], [0,0], [0,0])
	right_params = ([0,0], [0,0], [0,0], [0,0], [0,0])

	if left_user_barcode and right_user_barcode:
		left_params = GenerateBinary(left_user_barcode)
		right_params = GenerateBinary(right_user_barcode)
	elif left_user_barcode:
		left_params = GenerateBinary(left_user_barcode)
	elif right_user_barcode:
		right_params = GenerateBinary(right_user_barcode)
	else:
		printStatus("No internal barcodes/indices, or your -ub input is incorrect, or perhaps this is not the right tool for you?")
		sys.exit(1)

	# this is only needed to have the information whether there is going to be a binary barcode or not in advance before starting the unmerged/merged splitter
	# just about what folders/files can be deleted in case of a rerun
	if sum(left_params[2]):
		check_binary_bcode = True
	elif sum(right_params[2]):
		check_binary_bcode = True
	else:
		check_binary_bcode = False


	print "****************************************"
	print "Barcoding scheme:", user_barcode
	print "The barcode structure of forward reads: ", left_user_barcode
	print "The barcode structure of reverse reads: ", right_user_barcode
	print "The 'left' parameters", left_params
	print "The 'right' parameters", right_params
	print "****************************************\n"

	if not ignore: #ignore means "ignore unmerged reads."
		if script_process['BarcodeSplittingUnmerged'] == 'Finished':
			printStatus("Splitting of unmerged R1 and R2 reads was completed in a previous run. Skipping.")
			pass
		elif script_process['BarcodeSplittingUnmerged'] == 'Started':
			printStatus("Barcode splitting has been initiated before but not completed. Deleting previous files and restarting.")
			files_from_before = glob.glob(current_dir+"barcode_sorted/*fastq*")
			for i in files_from_before:
				os.remove(i)
			os.removedirs(current_dir+"bad_barcodes")
			with open(logfile, 'a') as ff:
				ff.write("** BarcodeSplittingUnmerged:Started:Again@ %s\n"%time.strftime("%X %x"))

			unmerged_statistics = UnmergedSplitter(forward, reverse, check_binary_bcode, gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, file_type, left_user_barcode, right_user_barcode, hamming_distance, logfile, current_dir)

		else:
			script_process['BarcodeSplittingUnmerged'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** BarcodeSplittingUnmerged:Started: %s\n"%time.strftime("%X %x"))

			unmerged_statistics = UnmergedSplitter(forward, reverse, check_binary_bcode, gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, file_type, left_user_barcode, right_user_barcode, hamming_distance, logfile, current_dir)

		printStatus("Done with the R1 and R2 reads.")


	if bbmerge:
		if script_process['BarcodeSplittingMerged'] == 'Finished':
			printStatus("Merged reads were split successfully in a previous run. Skipping.")
			pass

		elif script_process['BarcodeSplittingMerged'] == 'Started':
			printStatus("Barcode splitting has been initiated before but not completed. Deleting previous files and restarting.")
			files_from_before = glob.glob(current_dir+"barcode_sorted/*_merged_*")
			for i in files_from_before:
				os.remove(i)
			if ignore:
				os.removedirs(current_dir+"bad_barcodes")
			else:
				files_from_before = glob.glob(current_dir+"bad_barcodes/*merged*")
				for i in files_from_before:
					os.remove(i)
				if check_binary_bcode == True:
					files_from_before = glob.glob(current_dir+"bad_barcodes/*/*merged*")
					for i in files_from_before:
						os.remove(i)

			with open(logfile, 'a') as ff:
				ff.write("** BarcodeSplittingMerged:Started:Again@ %s\n"%time.strftime("%X %x"))

			merged_statistics = MergedSplitter(merged, check_binary_bcode, gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, left_user_barcode, right_user_barcode, hamming_distance, logfile, current_dir)

		else:
			script_process['BarcodeSplittingMerged'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** BarcodeSplittingMerged:Started: %s\n"%time.strftime("%X %x"))
			merged_statistics = MergedSplitter(merged, check_binary_bcode, gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, left_user_barcode, right_user_barcode, hamming_distance, logfile, current_dir)

		printStatus("Done with merged reads.")

	if merged_statistics and unmerged_statistics:
		total_reads = merged_statistics[0] + unmerged_statistics[0]
		total_alien = merged_statistics[1] + unmerged_statistics[1]
		total_bad = merged_statistics[2] + unmerged_statistics[2]
		if total_alien > 0:
			print "\n\n############################################################\n"
			print "Here's some barcode splitting statistics:\n"
			print "Total number of reads:\t %s" %(total_reads)
			print "Total merged reads:\t %s" %merged_statistics[0]
			print "Total unmerged reads:\t %s" %unmerged_statistics[0]
			print "\nNumber of reads that did not conform to the A/B binary barcode scheme:\t %s (%%%s of all reads)" %(total_bad, round((total_bad/total_reads)*100))
			print "\t%s of these came from merged reads\n\t%s of these came from R1 reads" %(merged_statistics[2], unmerged_statistics[2])
			print "Number of reads that contained an unidentifiable barcode:\t %s (%%%s of all reads)" %(total_alien, round((total_alien/total_reads)*100))
			print "\t%s of these came from merged reads\n\t%s of these came from R1 reads\n" %(merged_statistics[1], unmerged_statistics[1])
			for i in range(len(index_barcode_names)):
				print "Reads from replicate A of %s: %s" %(index_barcode_names[i], merged_statistics[3][i] + unmerged_statistics[3][i])
				print "Reads from replicate B of %s: %s\n" %(index_barcode_names[i], merged_statistics[4][i] + unmerged_statistics[4][i])
			print "############################################################\n\n"
		else:
			print "\n\n############################################################\n"
			print "Here's some barcode splitting statistics:\n"
			print "Total number of reads:\t %s" %(total_reads)
			print "Total merged reads:\t %s" %merged_statistics[0]
			print "Total unmerged reads:\t %s" %unmerged_statistics[0]
			print "Number of reads that contained an unidentifiable barcode:\t %s (%%%s of all reads)" %(total_bad, round((total_bad/total_reads)*100))
			print "\t%s of these came from merged reads\n\t%s of these came from R2 reads\n" %(merged_statistics[2], unmerged_statistics[2])
			for i in range(len(index_barcode_names)):
				print "Reads from sample %s: %s" %(index_barcode_names[i], merged_statistics[3][i] + unmerged_statistics[3][i])
			print "############################################################\n\n"

		cmd_bad_barcodes = "sort %sbad_barcodes/bad_barcodes_merged.txt | uniq -c | sort -k1,1nr > %sbad_barcodes/bad_barcodes_merged_with_counts.txt" %(current_dir, current_dir)
		subprocess.check_call(cmd_bad_barcodes, shell=True)
		cmd_bad_barcodes = "sort %sbad_barcodes/bad_barcodes.txt | uniq -c | sort -k1,1nr > %sbad_barcodes/bad_barcodes_with_counts.txt" %(current_dir, current_dir)
		subprocess.check_call(cmd_bad_barcodes, shell=True)

	if bbmerge:
		if dist:
			if script_process['InnerDistance'] == 'Finished':
				printStatus("Inner distance was successfully calculated in a previous run. Skipping.")
				pass
			elif script_process['InnerDistance'] == 'Started':
				aa = glob.glob(current_dir+"inner_distance/*/*")
				bb = glob.glob(current_dir+"inner_distance/*.gz")
				cc = glob.glob(current_dir+"inner_distance/*.bam")
				dd = glob.glob(current_dir+"inner_distance/*.log")
				files_from_before = aa + bb + cc + dd
				for i in files_from_before:
					os.remove(i)
				dirs_from_before = glob.glob(current_dir+"inner_distance/*")
				for i in dirs_from_before:
					os.removedirs(i)
				with open(logfile, 'a') as ff:
					ff.write("** InnerDistance:Started:Again@ %s\n"%time.strftime("%X %x"))
				printStatus("Calculating the fragment size.")
				innerdistance(forward, reverse, raw_fwd, raw_rvs, script_process, indices, mapping, logfile, current_dir, processor)
			else:
				script_process['InnerDistance'] = 'Started'
				with open(logfile, 'a') as ff:
					ff.write("** InnerDistance:Started %s\n"%time.strftime("%X %x"))
				printStatus("Calculating the fragment size.")
				innerdistance(forward, reverse, raw_fwd, raw_rvs, script_process, indices, mapping, logfile, current_dir, processor)
	else:
		if dist:
			if script_process['InnerDistance'] == 'Finished':
				printStatus("Inner distance was successfully calculated in a previous run. Skipping.")
				pass
			elif script_process['InnerDistance'] == 'Started':
				files_from_before = glob.glob(current_dir+"inner_distance/*/*")
				for i in files_from_before:
					os.remove(i)
				dirs_from_before = glob.glob(current_dir+"inner_distance/*")
				for i in dirs_from_before:
					os.removedirs(i)
				with open(logfile, 'a') as ff:
					ff.write("** InnerDistance:Started:Again@ %s\n"%time.strftime("%X %x"))
				printStatus("Calculating the fragment size.")
				innerdistance(None, None, raw_fwd, raw_rvs, script_process, indices, mapping, logfile, current_dir, processor)
			else:
				script_process['InnerDistance'] = 'Started'
				with open(logfile, 'a') as ff:
					ff.write("** InnerDistance:Started %s\n"%time.strftime("%X %x"))
				printStatus("Calculating the fragment size.")
				innerdistance(None, None, raw_fwd, raw_rvs, script_process, indices, mapping, logfile, current_dir, processor)

	if mapping and aligner == 'bt2-bbmap':
		merged_files = glob.glob(current_dir+"barcode_sorted/fastq_files/*merged*")

		if script_process['bowtie2'] == 'Finished':
			printStatus("Mapping with bowtie2 was successfully completed in a previous run. Skipping.")
			pass
		elif script_process['bowtie2'] == 'Started':
			aa = glob.glob(current_dir+"barcode_sorted/bowtie2/bam_files/*.bam")
			bb = glob.glob(current_dir+"barcode_sorted/bowtie2/unmapped/*fastq*")
			cc = glob.glob(current_dir+"barcode_sorted/bowtie2/log/*.log")
			files_from_before = aa + bb + cc
			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** bowtie2:Started:Again@ %s\n"%time.strftime("%X %x"))
			runBowtie2(merged_files, script_process, indices, mapping, processor, ignore, logfile, current_dir)

		else:
			script_process['bowtie2'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** bowtie2:Started: %s\n"%time.strftime("%X %x"))
			runBowtie2(merged_files, script_process, indices, mapping, processor, ignore, logfile, current_dir)

		unmapped_reads = glob.glob(current_dir+"barcode_sorted/bowtie2/unmapped/*_unmapped.fastq.gz")

		if script_process['bbmap'] == 'Finished':
			printStatus("BBMap has finished running before. Skipping.")
			pass
		elif script_process['bbmap'] == 'Started':
			aa = glob.glob(current_dir+"barcode_sorted/bbmap/bam_files/*.bam")
			bb = glob.glob(current_dir+"barcode_sorted/bbmap/log/*.log")
			files_from_before = aa + bb
			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** bbmap:Started:Again@ %s\n"%time.strftime("%X %x"))

			runBBMap(unmapped_reads, script_process, indices, mapping, merged_files, processor, logfile, current_dir)

		else:
			script_process['bbmap'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** bbmap:Started: %s\n"%time.strftime("%X %x"))

			runBBMap(unmapped_reads, script_process, indices, mapping, merged_files, processor, logfile, current_dir)

		bowtie_bams = glob.glob(current_dir+"barcode_sorted/bowtie2/bam_files/*.bam")
		bbmap_bams = glob.glob(current_dir+"barcode_sorted/bbmap/bam_files/*.bam")
		hisat_bams = None

		if script_process['sambedcleaner-bt2'] == 'Finished':
			printStatus("Bowtie2-BAM files were de-duplicated and cleaned up with sambedcleaner. Skipping.")
			pass
		elif script_process['sambedcleaner-bt2'] == 'Started':
			aa = glob.glob(current_dir+"barcode_sorted/bowtie2/cleaned/*/*")
			bb = glob.glob(current_dir+"barcode_sorted/bowtie2/bam_files/*.bw")
			cc = glob.glob(current_dir+"*temp*")
			dd = glob.glob(current_dir+"barcode_sorted/bowtie2/bam_files/bamCoverage.log")
			files_from_before = aa + bb + cc + dd
			for i in files_from_before:
				os.remove(i)
			dirs_from_before = glob.glob(current_dir+"barcode_sorted/bowtie2/cleaned/*")
			for i in dirs_from_before:
				os.removedirs(i)
			with open(logfile, 'a') as ff:
				ff.write("** sambedcleaner-bt2:Started:Again@ %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, hisat_bams, script_process, chunk, logfile, current_dir, keep_BAM)

		else:
			script_process['sambedcleaner-bt2'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** sambedcleaner-bt2:Started: %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, hisat_bams, script_process, chunk, logfile, current_dir, keep_BAM)




		if script_process['sambedcleaner-bb'] == 'Finished':
			printStatus("BBmap-BAM files were de-duplicated and cleaned up with sambedcleaner. Skipping.")
			pass
		elif script_process['sambedcleaner-bb'] == 'Started':
			aa = glob.glob(current_dir+"barcode_sorted/bbmap/cleaned/*/*")
			bb = glob.glob(current_dir+"barcode_sorted/bbmap/bam_files/*.bw")
			cc = glob.glob(current_dir+"*temp*")
			dd = glob.glob(current_dir+"barcode_sorted/bbmap/bam_files/bamCoverage.log")
			files_from_before = aa + bb + cc + dd
			for i in files_from_before:
				os.remove(i)
			dirs_from_before = glob.glob(current_dir+"barcode_sorted/bbmap/cleaned/*")
			for i in dirs_from_before:
				os.removedirs(i)
			with open(logfile, 'a') as ff:
				ff.write("** sambedcleaner-bb:Started:Again@ %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, hisat_bams, script_process, chunk, logfile, current_dir, keep_BAM)

		else:
			script_process['sambedcleaner-bb'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** sambedcleaner-bb:Started: %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, hisat_bams, script_process, chunk, logfile, current_dir, keep_BAM)

		# Merge the bam files from bowtie2 and bbmap
		cleaned_files = glob.glob(current_dir+"barcode_sorted/bowtie2/cleaned/*/*cleaned.bam")

		if script_process['MergeBowtie2BBMap'] == 'Finished':
			printStatus("BBMap and bowtie2 outputs were merged before. Skipping.")
			pass
		elif script_process['MergeBowtie2BBMap'] == 'Started':
			files_from_before = glob.glob(current_dir+"final_files/bams/*bam*")
			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** MergeBowtie2BBMap:Started:Again@ %s\n"%time.strftime("%X %x"))
			BamMerger(cleaned_files, script_process, logfile, current_dir)
		else:
			script_process['MergeBowtie2BBMap'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** MergeBowtie2BBMap:Started: %s\n"%time.strftime("%X %x"))
			BamMerger(cleaned_files, script_process, logfile, current_dir)

		# Generate coverage files for individual replicates
		if script_process['GenerateCov'] == 'Finished':
			printStatus("Coverages were generated before. Skipping")
			pass
		elif script_process['GenerateCov'] == 'Started':
			aa = glob.glob(current_dir+"final_files/bigwigs/*/*.bw")
			bb = glob.glob(current_dir+"final_files/p5end/*/*bw")
			cc = glob.glob(current_dir+"final_files/p3end/*/*bw")
			files_from_before = aa + bb + cc

			for i in files_from_before:
				os.remove(i)

			with open(logfile, 'a') as ff:
				ff.write("** GenerateCov:Started:Again@ %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping, logfile, current_dir, processor)
		else:
			script_process['GenerateCov'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** GenerateCov:Started: %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping, logfile, current_dir, processor)

		if script_process['EnrichmentPlot'] == 'Finished':
			printStatus("Enrichments of features were calculated before. Skipping")
			pass
		elif script_process['EnrichmentPlot'] == 'Started':
			files_from_before = glob.glob(current_dir+"final_files/EnrichmentPlots/*.png")

			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** EnrichmentPlot:Started:Again@ %s\n"%time.strftime("%X %x"))
			plotEnrichment(script_process, indices, mapping, logfile, current_dir, processor, aligner)
		else:
			script_process['EnrichmentPlot'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** EnrichmentPlot:Started: %s\n"%time.strftime("%X %x"))
			plotEnrichment(script_process, indices, mapping, logfile, current_dir, processor, aligner)

	if mapping and aligner == 'hisat':
		merged_files = glob.glob(current_dir+"barcode_sorted/fastq_files/*merged*")

		if script_process['hisat'] == 'Finished':
			printStatus("Mapping with HISAT2 was successfully completed in a previous run. Skipping.")
			pass
		elif script_process['hisat'] == 'Started':
			aa = glob.glob(current_dir+"barcode_sorted/hisat/bam_files/*.bam")
			bb = glob.glob(current_dir+"barcode_sorted/hisat/unmapped/*fastq*")
			cc = glob.glob(current_dir+"barcode_sorted/hisat/log/*.log")
			files_from_before = aa + bb + cc
			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** hisat:Started:Again@ %s\n"%time.strftime("%X %x"))
			runHiSat(merged_files, script_process, indices, mapping, processor, logfile, current_dir)

		else:
			script_process['hisat'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** hisat:Started: %s\n"%time.strftime("%X %x"))
			runHiSat(merged_files, script_process, indices, mapping, processor, logfile, current_dir)

		hisat_bams = glob.glob(current_dir+"barcode_sorted/hisat/bam_files/*.bam")
		bbmap_bams = None
		bowtie_bams = None
		unmapped_reads = glob.glob(current_dir+"barcode_sorted/hisat/unmapped/*_unmapped.fastq.gz")

		if script_process['sambedcleaner'] == 'Finished':
			printStatus("BAM files were de-duplicated and cleaned up with sambedcleaner. Skipping.")
			pass
		elif script_process['sambedcleaner'] == 'Started':
			aa = glob.glob(current_dir+"barcode_sorted/hisat/cleaned/*/*")
			bb = glob.glob(current_dir+"barcode_sorted/hisat/bam_files/*.bw")
			cc = glob.glob(current_dir+"*temp*")
			dd = glob.glob(current_dir+"barcode_sorted/hisat/bam_files/bamCoverage.log")
			files_from_before = aa + bb + cc + dd
			for i in files_from_before:
				os.remove(i)
			dirs_from_before = glob.glob(current_dir+"barcode_sorted/hisat/cleaned/*")
			for i in dirs_from_before:
				os.removedirs(i)
			with open(logfile, 'a') as ff:
				ff.write("** sambedcleaner:Started:Again@ %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, hisat_bams, script_process, chunk, logfile, current_dir, keep_BAM)

		else:
			script_process['sambedcleaner'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** sambedcleaner:Started: %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, hisat_bams, script_process, chunk, logfile, current_dir, keep_BAM)

		# Copy bams to final folder
		cleaned_files = glob.glob(current_dir+"barcode_sorted/hisat/cleaned/*/*cleaned.bam*")

		if script_process['MoveBams'] == 'Finished':
			printStatus("HISAT .bams moved and indexed. Skipping.")
			pass
		elif script_process['MoveBams'] == 'Started':
			files_from_before = glob.glob(current_dir+"final_files/bams/*.bam*")
			for i in files_from_before:
				os.remove(i)

			with open(logfile, 'a') as ff:
				ff.write("** MoveBams:Started:Again@ %s\n"%time.strftime("%X %x"))

			for i in cleaned_files:
				cmd = "cp %s %sfinal_files/bams/%s" %(i, current_dir, i.split('/')[-1])
				subprocess.check_call(cmd, shell=True)


			#moved_bams = glob.glob(current_dir+"final_files/bams/*.bam")

			#for i in moved_bams:
			#	cmd = "samtools index %s" %i
			#	subprocess.check_call(cmd, shell=True)
			script_process['MoveBams'] = 'Finished'
			with open(logfile, 'a') as ff:
				ff.write("** MoveBams:Finished: %s\n"%time.strftime("%X %x"))

		else:
			script_process['MoveBams'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** MoveBams:Started: %s\n"%time.strftime("%X %x"))

			for i in cleaned_files:
				cmd = "cp %s %sfinal_files/bams/%s" %(i, current_dir, i.split('/')[-1])
				subprocess.check_call(cmd, shell=True)

			#moved_bams = glob.glob(current_dir+"final_files/bams/*.bam")

			#for i in moved_bams:
			#	cmd = "samtools index %s" %i
			#	subprocess.check_call(cmd, shell=True)
			script_process['MoveBams'] = 'Finished'
			with open(logfile, 'a') as ff:
				ff.write("** MoveBams:Finished: %s\n"%time.strftime("%X %x"))

		# Generate coverage files for individual replicates
		if script_process['GenerateCov'] == 'Finished':
			printStatus("Coverages were generated before. Skipping")
			pass

		elif script_process['GenerateCov'] == 'Started':
			aa = glob.glob(current_dir+"final_files/bigwigs/*/*.bw")
			bb = glob.glob(current_dir+"final_files/p5end/*/*bw")
			cc = glob.glob(current_dir+"final_files/p3end/*/*bw")
			files_from_before = aa + bb + cc

			for i in files_from_before:
				os.remove(i)

			with open(logfile, 'a') as ff:
				ff.write("** GenerateCov:Started:Again@ %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping, logfile, current_dir, processor)
		else:
			script_process['GenerateCov'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** GenerateCov:Started: %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping, logfile, current_dir, processor)

		if script_process['EnrichmentPlot'] == 'Finished':
			printStatus("Enrichments of features were calculated before. Skipping")
			pass
		elif script_process['EnrichmentPlot'] == 'Started':
			files_from_before = glob.glob(current_dir+"final_files/EnrichmentPlots/*.png")

			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** EnrichmentPlot:Started:Again@ %s\n"%time.strftime("%X %x"))
			plotEnrichment(script_process, indices, mapping, logfile, current_dir, processor, aligner)
		else:
			script_process['EnrichmentPlot'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** EnrichmentPlot:Started: %s\n"%time.strftime("%X %x"))
			plotEnrichment(script_process, indices, mapping, logfile, current_dir, processor, aligner)

	if mapping and aligner == 'bt2':
		merged_files = glob.glob(current_dir+"barcode_sorted/fastq_files/*merged*")

		if script_process['bowtie2'] == 'Finished':
			printStatus("Mapping with Bowtie2 was successfully completed in a previous run. Skipping")
			pass
		elif script_process['bowtie2'] == 'Started':
			aa = glob.glob(current_dir+"barcode_sorted/bowtie2/bam_files/*.bam")
			bb = glob.glob(current_dir+"barcode_sorted/bowtie2/unmapped/*fastq*")
			cc = glob.glob(current_dir+"barcode_sorted/bowtie2/log/*.log")
			files_from_before = aa + bb + cc
			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** bowtie2:Started:Again@ %s\n"%time.strftime("%X %x"))
			runBowtie2(merged_files, script_process, indices, mapping, processor, ignore, logfile, current_dir)

		else:
			script_process['bowtie2'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** bowtie2:Started: %s\n"%time.strftime("%X %x"))
			runBowtie2(merged_files, script_process, indices, mapping, processor, ignore, logfile, current_dir)

		bowtie_bams = glob.glob(current_dir+"barcode_sorted/bowtie2/bam_files/*.bam")
		bbmap_bams = None
		hisat_bams = None
		unmapped_reads = glob.glob(current_dir+"barcode_sorted/bowtie2/unmapped/*_unmapped.fastq.gz")

		if script_process['sambedcleaner'] == 'Finished':
			printStatus("BAM files were de-duplicated and cleaned up with sambedcleaner. Skipping.")
			pass
		elif script_process['sambedcleaner'] == 'Started':
			aa = glob.glob(current_dir+"barcode_sorted/bowtie2/cleaned/*/*")
			bb = glob.glob(current_dir+"barcode_sorted/bowtie2/bam_files/*.bw")
			cc = glob.glob(current_dir+"*temp*")
			dd = glob.glob(current_dir+"barcode_sorted/bowtie2/bam_files/bamCoverage.log")
			files_from_before = aa + bb + cc + dd
			for i in files_from_before:
				os.remove(i)
			dirs_from_before = glob.glob(current_dir+"barcode_sorted/bowtie2/cleaned/*")
			for i in dirs_from_before:
				os.removedirs(i)
			with open(logfile, 'a') as ff:
				ff.write("** sambedcleaner:Started:Again@ %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, hisat_bams, script_process, chunk, logfile, current_dir, keep_BAM)

		else:
			script_process['sambedcleaner'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** sambedcleaner:Started: %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, hisat_bams, script_process, chunk, logfile, current_dir, keep_BAM)

		# Copy bams to final folder
		cleaned_files = glob.glob(current_dir+"barcode_sorted/bowtie2/cleaned/*/*cleaned.bam*")
		#raw_bams = glob.glob(current_dir+"")

		if script_process['MoveBams'] == 'Finished':
			printStatus("bowtie2.bams moved and indexed. Skipping.")
			pass
		elif script_process['MoveBams'] == 'Started':
			files_from_before = glob.glob(current_dir+"final_files/bams/*bam*")
			for i in files_from_before:
				os.remove(i)

			with open(logfile, 'a') as ff:
				ff.write("** MoveBams:Started:Again@ %s\n"%time.strftime("%X %x"))

			for i in cleaned_files:
				cmd = "cp %s %sfinal_files/bams/%s" %(i, current_dir, i.split('/')[-1])
				subprocess.check_call(cmd, shell=True)

			#moved_bams = glob.glob(current_dir+"final_files/bams/*.bam")

			#for i in moved_bams:
			#	cmd = "samtools index %s" %i
			#	subprocess.check_call(cmd, shell=True)
			script_process['MoveBams'] = 'Finished'
			with open(logfile, 'a') as ff:
				ff.write("** MoveBams:Finished: %s\n"%time.strftime("%X %x"))

		else:
			script_process['MoveBams'] == 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** MoveBams:Started: %s\n"%time.strftime("%X %x"))

			for i in cleaned_files:
				cmd = "cp %s %sfinal_files/bams/%s" %(i, current_dir, i.split('/')[-1])
				subprocess.check_call(cmd, shell=True)

			#moved_bams = glob.glob(current_dir+"final_files/bams/*.bam")

			#for i in moved_bams:
			#	cmd = "samtools index %s" %i
			#	subprocess.check_call(cmd, shell=True)
			script_process['MoveBams'] = 'Finished'
			with open(logfile, 'a') as ff:
				ff.write("** MoveBams:Finished: %s\n"%time.strftime("%X %x"))

		# Generate coverage files for individual replicates
		if script_process['GenerateCov'] == 'Finished':
			printStatus("Coverages were generated before. Skipping")
			pass
		elif script_process['GenerateCov'] == 'Started':
			aa = glob.glob(current_dir+"final_files/bigwigs/*/*.bw")
			bb = glob.glob(current_dir+"final_files/p5end/*/*bw")
			cc = glob.glob(current_dir+"final_files/p3end/*/*bw")
			files_from_before = aa + bb + cc

			for i in files_from_before:
				os.remove(i)


			with open(logfile, 'a') as ff:
				ff.write("** GenerateCov:Started:Again@ %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping, logfile, current_dir, processor)
		else:
			script_process['GenerateCov'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** GenerateCov:Started: %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping, logfile, current_dir, processor)

		if script_process['EnrichmentPlot'] == 'Finished':
			printStatus("Enrichments of features were calculated before. Skipping")
			pass
		elif script_process['EnrichmentPlot'] == 'Started':
			files_from_before = glob.glob(current_dir+"final_files/EnrichmentPlots/*.png")

			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** EnrichmentPlot:Started:Again@ %s\n"%time.strftime("%X %x"))
			plotEnrichment(script_process, indices, mapping, logfile, current_dir, processor, aligner)
		else:
			script_process['EnrichmentPlot']  = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** EnrichmentPlot:Started: %s\n"%time.strftime("%X %x"))
			plotEnrichment(script_process, indices, mapping, logfile, current_dir, processor, aligner)

	if spikes_mapper:
		if script_process['spikesMapper'] == 'Finished':
			printStatus("Spikes were already mapped successfully in a previous run. Skipping.")
			pass

		elif script_process['spikesMapper'] == 'Started':
			printStatus("Spikes mapping has been initiated before but not completed. Deleting previous files and restarting.")
			files_from_before = glob.glob(current_dir+"spikes/*.bam")

			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** spikesMapper:Started:Again@ %s\n"%time.strftime("%X %x"))

			spikeMapper(unmapped_reads, script_process, indices, mapping, processor, merged_files, logfile, current_dir, spikes_mapper)

		else:
			script_process['spikesMapper'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** spikesMapper:Started: %s\n"%time.strftime("%X %x"))

			spikeMapper(unmapped_reads, script_process, indices, mapping, processor, merged_files, logfile, current_dir, spikes_mapper)

		spike_bams = glob.glob(current_dir+"spikes/*.bam")

		if script_process['spikesBamCleaner'] == 'Finished':
			printStatus("Spikes were de-duplicated and cleaned up with BamCleaner. Skipping.")
			pass

		elif script_process['spikesBamCleaner'] == 'Started':
			printStatus("BamCleaner on Spikes has been initiated before but not completed. Deleting previous files and restarting.")
			files_from_before = glob.glob(current_dir+"spikes/cleaned/*.bam")
			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** spikesBamCleaner:Started:Again@ %s\n"%time.strftime("%X %x"))

			runSpikeCleaner(spike_bams, script_process, chunk, logfile, current_dir)

		else:
			script_process['spikesBamCleaner'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** spikesBamCleaner:Started: %s\n"%time.strftime("%X %x"))

			runSpikeCleaner(spike_bams, script_process, chunk, logfile, current_dir)

	# Merge replicates
	if reps_merged:
		if script_process['MergeReps'] == 'Finished':
			printStatus("Replicates were merged before. Skipping.")
			pass

		elif script_process['MergeReps'] == 'Started':
			aa = glob.glob("final_files/replicates_merged/bams/*bam*")
			bb = glob.glob("final_files/replicates_merged/bigwigs/raw/*.bw")
			cc = glob.glob("final_files/replicates_merged/bigwigs/CPM_normalized/*.bw")
			dd = glob.glob("final_files/replicates_merged/p5end/raw/*.bw")
			ee = glob.glob("final_files/replicates_merged/p5end/CPM_normalized/*.bw")
			ff = glob.glob("final_files/replicates_merged/p3end/raw/*.bw")
			gg = glob.glob("final_files/replicates_merged/p3end/CPM_normalized/*.bw")
			files_from_before = aa + bb + cc + dd + ee + ff + gg

			for i in files_from_before:
				os.remove(i)
			with open(logfile, 'a') as ff:
				ff.write("** MergeReps:Started:Again@ %s\n"%time.strftime("%X %x"))

			RepMerger(index_barcode_names, script_process, indices, mapping, aligner, logfile, current_dir, processor)

		else:
			script_process['MergeReps'] = 'Started'
			with open(logfile, 'a') as ff:
				ff.write("** MergeReps:Started: %s\n"%time.strftime("%X %x"))

			RepMerger(index_barcode_names, script_process, indices, mapping, aligner, logfile, current_dir, processor)

	if gzipafter:
		if script_process['gzipping'] == 'Finished':
			printStatus("Gzipping of fastq files was done before. Skipping")
			pass

		elif script_process['gzipping'] == 'Started':

			with open(logfile, 'a') as ff:
				ff.write("** gzipping:Started:Again@ %s\n"%time.strftime("%X %x"))

			fastq_files = glob.glob(current_dir+"barcode_sorted/fastq_files/*.fastq")
			if len(fastq_files) > 0:
				cmd = "gzip %sbarcode_sorted/fastq_files/*.fastq" %current_dir
				printStatus("Gzipping all the fastq files generated after barcode splitting.")
				subprocess.check_call(cmd, shell=True)
				printStatus("Gzipping done.")
			fastq_files2 = glob.glob(current_dir+"bad_barcodes/*.fastq")
			if len(fastq_files2) > 0:
				cmd2 = "gzip %sbad_barcodes/*.fastq" %current_dir
				printStatus("Gzipping all the fastq files with unknown binary barcodes.")
				subprocess.check_call(cmd2, shell=True)
				printStatus("Gzipping done.")
			fastq_files3 = glob.glob(current_dir+"bad_barcodes/*/*.fastq")
			if len(fastq_files3) > 0:
				cmd3 = "gzip %sbad_barcodes/*/*.fastq" %current_dir
				printStatus("Gzipping all the fastq files with bad barcodes.")
				subprocess.check_call(cmd3, shell=True)
				printStatus("Gzipping done.")

			script_process['gzipping'] = 'Finished'
			with open(logfile, 'a') as ff:
				ff.write("** gzipping:Finished: %s\n"%time.strftime("%X %x"))

		else:
			script_process['gzipping'] = 'Started'

			with open(logfile, 'a') as ff:
				ff.write("** gzipping:Started: %s\n"%time.strftime("%X %x"))

			fastq_files = glob.glob(current_dir+"barcode_sorted/fastq_files/*.fastq")
			if len(fastq_files) > 0:
				cmd = "gzip %sbarcode_sorted/fastq_files/*.fastq" %current_dir
				printStatus("Gzipping all the fastq files generated after barcode splitting.")
				subprocess.check_call(cmd, shell=True)
				printStatus("Gzipping done.")
			fastq_files2 = glob.glob(current_dir+"bad_barcodes/*.fastq")
			if len(fastq_files2) > 0:
				cmd2 = "gzip %sbad_barcodes/*.fastq" %current_dir
				printStatus("Gzipping all the fastq files with unknown binary barcodes.")
				subprocess.check_call(cmd2, shell=True)
				printStatus("Gzipping done.")
			fastq_files3 = glob.glob(current_dir+"bad_barcodes/*/*.fastq")
			if len(fastq_files3) > 0:
				cmd3 = "gzip %sbad_barcodes/*/*.fastq" %current_dir
				printStatus("Gzipping all the fastq files with bad barcodes.")
				subprocess.check_call(cmd3, shell=True)
				printStatus("Gzipping done.")

			script_process['gzipping'] = 'Finished'
			with open(logfile, 'a') as ff:
				ff.write("** gzipping:Finished: %s\n"%time.strftime("%X %x"))

def main():
	parser = argparse.ArgumentParser(description='Demultiplexer for FLASH data. Takes in fastq files, list of barcodes, throws out the reads separated for their barcode. Can be used also for iCLIP, uvCLAP or PAR-iCLIP. Just set your internal barcoding scheme with the -ub option.')
	parser.add_argument('-f', '--forward', metavar='path/to/left_reads_R1.fastq.gz', help='The fastq file with the forward reads', required=True)
	parser.add_argument('-r', '--reverse', metavar='path/to/right_reads_R2.fastq.gz', help='The fastq file with the reverse reads')
	parser.add_argument('-b', '--barcodes', metavar='path/to/barcodes.fa', help='The absolute path to the fasta file that contains the barcodes.' ,required=True)
	parser.add_argument('-bm', '--bbmerge', action='store_true', help='Use BBmerge (you need to have it in your PATH) to preprocess the reads. Highly recommended. Using this option will create and use the merged data from BBmerge.')
	parser.add_argument('-i', '--ignore', action='store_true', help='Only use the merged file generated by BBMerge and ignore the non-mergable ones, assumption being that there\'s something wrong with them. Use only with the -bm option, otherwise the program won\'t do anything.')
	parser.add_argument('-c', '--chunk', type=int, help='Determine the chunk size. The formula is chunk * 100000 per file. Default is 2.', default=2)
	parser.add_argument('-gz', '--gzip', action='store_true', help='Gzip the fastq files on the fly to save space. Slower. I just use -gza instead of this, but this is also an option.')
	parser.add_argument('-mr', '--reps_merged', action='store_true', help='Merge the biological replicates A and B in addition to the normal stuff. Slower')
	parser.add_argument('-m', '--mapping', choices=['hg19','hg38', 'dm3', 'dm6', 'mm10', 'mm9', 'hg19_simple', 'dm6_simple', 'dvir'])
	parser.add_argument('-gza', '--gzipafter', action='store_true', help='Runs gzip on the split files, AFTER everything else is done. Much faster than -gz option.')
	parser.add_argument('-ub', '--user_barcode', type=str, help='User defined internal barcode. For iCLIP it would be L=NNNXXXXNN; for uvCLAP it is L=NNNXXXXXNN,R=YRRYN; for FLASH it is R=NNYYNXXXXXXNN. Default is FLASH', default='FLASH')
	parser.add_argument('-hm', '--hamming_distance', type=int, help="Maximum allowed mutations in the internal index. Default is 1 mutation allowed. For shorter barcodes, can 0 (i.e. perfect match). For longer barcodes can use larger numbers. Experiment with this number for best results", default=1)
	parser.add_argument('-p', '--processor', type=int, help="Number of processors to use. This is only passed to bowtie2 BBMap, and bbduk. Default = 8 ", default=8)
	parser.add_argument('-mi', '--minimum_insert', type=int, help='Minimun insert size for BBMerge. Default is 30 (with all the internal barcodes, so subtract that to get the actualy minimum insert size you are setting.)', default=30)
	parser.add_argument('-al', '--aligner', type=str, choices=['bt2', 'bbmap', 'bt2-bbmap', 'hisat'], default='bt2-bbmap', help='Aligner to use. Default is bowtie2 first, and BBMap on the reads that couldn\'t be aligned with bowtie2')
	parser.add_argument('-tr', '--trim', type=str, default='TACACTCTTTCCCTACACGACGCTCTTCCGATCT,AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', help='Remove adapters with BBDuk. Default is Illumina adapters')
	parser.add_argument('-s', '--spikes_mapper', choices=['fetish', 'mapcap', 'srmc'], help='Use Bowtie2 to map the spike-ins. Choices are FETISH, MAPCap or smallRNA-MAPCap (srmc).')
	parser.add_argument('-k', '--keep_BAM', action='store_true', help='Wheter you want to save the reads identified as PCR duplicates and/or multimappers in separate BAM files.')
	parser.add_argument('-d', '--dist', action='store_true', help='Calculate the fragment size of the sequenced reads.')
	args = parser.parse_args()

	current_dir = os.path.abspath(args.barcodes).split("barcodes")[0]

	script_process = {"BBMerge":None, "BarcodeSplittingUnmerged":None, "BarcodeSplittingMerged":None, "InnerDistance":None,
	"bowtie2":None, "bbmap":None, "sambedcleaner":None, "sambedcleaner-bb":None, "sambedcleaner-bt2":None, "MergeBowtie2BBMap":None, "MergedSort":None,
	"MergeReps":None, "GenerateCov":None, "hisat":None, "MoveBams":None, "bbmap-only":None,
	"spikesMapper":None, "spikesBamCleaner":None, "EnrichmentPlot":None, "gzipping":None}

	# write out the current packages used and their version
	#if os.path.exists(current_dir+"packages.txt"):
	#	print
	#else:
	#	i=0
	#	for i in range(len(path_list)):
	#		with open(current_dir+"packages.txt", 'a') as ff:
	#			ff.write("%s\t\t\t%s\n" %(path_list[i], version_list[i]))

	logfile = current_dir+"logfile.log"

	try:
		with open(logfile, 'r') as ff:
			printStatus("Reading the logfile.")
			while True:
				try:
					line = ff.next().strip()
					if line.split(" ")[0] == "**":
						process, result = line.split(" ")[1].split(":")[0], line.split(" ")[1].split(":")[1]
						script_process[process] = result
				except StopIteration:
					break
	except IOError:
		printStatus("Fresh start.")
		with open(logfile, 'a') as ff:
			date_time = time.strftime("%X %x")
			ff.write("Starting at: %s\n\nExecuting %s with the command:\n\n%s\n\n----------\nThe following programs are used in this script:\n" %(date_time, __version__, " ".join(sys.argv[:])))
			for i in range(len(path_list)):
				ff.write("%s\t\t\t version: %s\n" %(path_list[i], version_list[i]))
			ff.write("------\n")


	if re.search('fastq', args.forward, re.IGNORECASE):
		if re.search('fastq.gz', args.forward, re.IGNORECASE):
			filetype = 'fastq.gz'
		else:
			filetype = 'fastq'
	else:
		print("I don't understand your file type. Your files should end either with .fastq or .fastq.gz")
		sys.exit(1)

	name_forward = os.path.splitext(args.forward)[0].split('/')[-1].split("R1")[0] + 'bbmergeOUT_R1.fastq.gz'
	name_reverse = os.path.splitext(args.reverse)[0].split('/')[-1].split("R2")[0] + 'bbmergeOUT_R2.fastq.gz'
	name_merged = os.path.splitext(args.forward)[0].split('/')[-1].split("R1.fastq")[0] + 'merged_bbmergeOUT.fastq.gz'

	bbmerge_forward = current_dir+"bbmap_merged/"+name_forward
	bbmerge_reverse = current_dir+"bbmap_merged/"+name_reverse
	bbmerged_outfile = current_dir+"bbmap_merged/"+name_merged

	if args.bbmerge:
		if script_process['BBMerge'] == 'Finished':
			printStatus("Using the output of a previous BBmerge run. If you don't want this, stop the script, delete bbmap_merged.fastq.gz and other BBmerge files and start again.")
			processFastqFiles(bbmerge_forward, bbmerge_reverse, bbmerged_outfile, args.barcodes, filetype, args.bbmerge, args.ignore, args.chunk*100000,
				args.mapping, args.gzip, args.reps_merged, script_process, args.gzipafter, args.user_barcode, args.hamming_distance,
				args.processor, args.aligner, args.spikes_mapper, logfile, current_dir, args.keep_BAM, args.forward, args.reverse, args.dist)
		elif script_process['BBMerge'] == 'Started':
			printStatus("BBMerge was initiated but not completed, deleting partially processed files and restarting.")
			os.remove(bbmerge_forward)
			os.remove(bbmerge_reverse)
			os.remove(bbmerged_outfile)

			RunBBMerge(script_process, args, current_dir, bbmerge_forward, bbmerge_reverse, bbmerged_outfile, args.minimum_insert, logfile)

			print "\n"
			printStatus("Now running the splitter.\n")
			processFastqFiles(bbmerge_forward, bbmerge_reverse, bbmerged_outfile, args.barcodes, filetype, args.bbmerge, args.ignore, args.chunk*100000,
				args.mapping, args.gzip, args.reps_merged, script_process, args.gzipafter, args.user_barcode, args.hamming_distance,
				args.processor, args.aligner, args.spikes_mapper, logfile, current_dir, args.keep_BAM, args.forward, args.reverse, args.dist)
		else:
			script_process['BBMerge'] == 'Started'
			with open(logfile, 'a') as ff:
				ff.write("\n** BBMerge:Started: %s\n"%time.strftime("%X %x"))

			RunBBMerge(script_process, args, current_dir, bbmerge_forward, bbmerge_reverse, bbmerged_outfile, args.minimum_insert, logfile)

			printStatus("Now running the splitter.\n")
			processFastqFiles(bbmerge_forward, bbmerge_reverse, bbmerged_outfile, args.barcodes, filetype, args.bbmerge, args.ignore, args.chunk*100000,
				args.mapping, args.gzip, args.reps_merged, script_process, args.gzipafter, args.user_barcode, args.hamming_distance,
				args.processor, args.aligner, args.spikes_mapper, logfile, current_dir, args.keep_BAM, args.forward, args.reverse, args.dist)
	else:
		processFastqFiles(args.forward, args.reverse, False, args.barcodes, filetype, False, False, args.chunk*100000,
			False, args.gzip, False, script_process, args.gzipafter, args.user_barcode, args.hamming_distance,
			args.processor, args.aligner, args.spikes_mapper, logfile, current_dir, args.keep_BAM, args.forward, args.reverse, args.dist)

####################################################
# Some useful functions here
####################################################


def extract_barcodes(barcodes):
	"""This script reads the barcodes fasta file and returns two lists,
	one with the names of the samples, the other with the barcodes """
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

def printStatus(msg):
	print "%s: %s" % (time.strftime("%X %x"), msg)

def RunBBMerge(script_process, args, current_dir, bbmerge_forward, bbmerge_reverse, bbmerged_outfile, minimum_insert, logfile):
	""" This just runs BBMerge on the input files. Doesn't return anything (probably returns 1 or 0?)
	The output is named bbmap_merged.fastq.gz. This cannot be changed right now, but I guess should be changed at some point."""

	folder = current_dir+"bbmap_merged"
	if not os.path.exists(folder):
		os.makedirs(folder)

	printStatus("Running BBMerge now.\n")
	BBMERGE_COMMAND = "mininsert=%s in1=%s in2=%s out=%s outu1=%s outu2=%s 2>&1"  %(minimum_insert, os.path.abspath(args.forward), os.path.abspath(args.reverse),
		bbmerged_outfile, bbmerge_forward, bbmerge_reverse)
	cmd = "%s %s" %(bbmerge_path, BBMERGE_COMMAND)
	subprocess.check_call(cmd, shell=True)

	printStatus("\nFinished BBmerge.")
	script_process['BBMerge'] = 'Finished'
	with open(logfile, 'a') as ff:
		ff.write("** BBMerge:Finished: %s\n"%time.strftime("%X %x"))

def GenerateBinary(user_barcode):
	"""Takes in a string that is composed of n, x, y or r letters
	such as 'nnnxxxyrrynn' (lower case), returns two sets that contain all possible combinations
	of the binary part (y/r) in DNA language, index of the binary barcode as a list, where the first element
	is the beginning position and the second the end,
	index of the internal index as a list (same as above), and finally a list that contains
	the positions of random nucleotides as its elemenents."""

	barcode = user_barcode.lower()
	code_a = ''.join(''.join(barcode.split('n')).split('x'))
	code_a_index = [ barcode.index(code_a), barcode.index(code_a) + len(code_a)]
	internal_index = ''.join(''.join(''.join(barcode.split('r')).split('y')).split('n'))
	internal_index_index = [barcode.index(internal_index), barcode.index(internal_index) + len(internal_index)]
	all_indices = [x for x in range(len(barcode))]

	#remove binary and internal index from the user barcode, what is left is the positions of the random barcode
	#e.g. read = 'ABNHJKIOLKJHDS'
	#rand = ''.join([read[x] for x in all_indices]) will give the random part as a string

	for i in range(code_a_index[0],code_a_index[1]):
		all_indices.remove(i)

	for i in range(internal_index_index[0],internal_index_index[1]):
		all_indices.remove(i)

	convert = {'y':'r', 'r':'y'}
	code_b = "".join([convert[x] for x in code_a])

	letters = {'r':['A', 'G'], 'y':['C', 'T']}

	indices = list(itertools.product([0, 1], repeat=len(code_a)))
	clist_a = [x for x in code_a]
	clist_b = [x for x in code_b]
	bcode_a = set()
	bcode_b = set()

	for i in indices:
		new_list = []
		for x in range(len(code_a)):
			a = i[x]
			new_list.append(letters[clist_a[x]][a])
		new_code = "".join(new_list)
		bcode_a.add(new_code)

	for i in indices:
		new_list = []
		for x in range(len(code_b)):
			b = i[x]
			new_list.append(letters[clist_b[x]][b])
		new_code = "".join(new_list)
		bcode_b.add(new_code)

	return bcode_a, bcode_b, code_a_index, internal_index_index, all_indices

def UnmergedSplitter(forward, reverse, check_binary_bcode, gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, file_type, left_user_barcode, right_user_barcode, hamming_distance, logfile, current_dir):
	""" This is the first splitter function. This one works on unmerged R1 and R2 files.
	It basically does 3 things at the same time:
	1) Finds to what sample the reads belongs to. Looks for an index that mathces with at most 1 mismatch
	2) See if the reads has an A or B binary barcode, and split accordingly
	3) Collect splitting statistics to quickly check how many reads comes from which sample in the end """

	folders = [current_dir+"bad_barcodes", current_dir+"barcode_sorted", current_dir+"barcode_sorted/fastq_files"]

	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)

	if check_binary_bcode == True:
		bad_barcodes_folders =[]

		for i in index_barcode_names:
			bad_barcodes_folders.append(current_dir+"bad_barcodes/"+i)

		for i in bad_barcodes_folders:
			if not os.path.exists(i):
				os.makedirs(i)

	printStatus("Will process R1 and R2 files now.")

	left_reads = gzip.open(forward, 'rb')
	right_reads = gzip.open(reverse, 'rb')

	multiplier = 1
	count1 = 0 #total number of reads
	alien_count = 0 # not A or B.
	bad_barcode = 0 # A or B but unknown barcode.

	single_reads_f_dict = {}
	single_reads_r_dict = {}
	alien_reads_f_dict = {}
	alien_reads_r_dict = {}
	bad_barcode_reads_f_dict = {}
	bad_barcode_reads_r_dict = {}
	bad_set = []
	a_barcodes_unmerged = [0 for x in range(len(index_barcode_names))]
	b_barcodes_unmerged = [0 for x in range(len(index_barcode_names))]

	while True:
		try:
			name_f, name_r	= left_reads.next().strip(), right_reads.next().strip()
			seq_f,	seq_r	= left_reads.next().strip(), right_reads.next().strip()
			plus_f, plus_r	= left_reads.next().strip(), right_reads.next().strip()
			qual_f, qual_r	= left_reads.next().strip(), right_reads.next().strip()
			count1 += 1
			# collect the random barcode from both ends, wherever it is defined

			if sum(left_params[4]) != 0 and sum(right_params[4]) != 0:
				ran_bcode = ''.join([seq_f[x] for x in left_params[4]]) + ''.join([seq_r[x] for x in right_params[4]])
				left_trim_l = len(left_user_barcode)
				right_trim_l = -len(right_user_barcode)
				left_trim_r = -len(left_user_barcode)
				right_trim_r = len(right_user_barcode)

			elif sum(left_params[4]):
				ran_bcode = ''.join([seq_f[x] for x in left_params[4]])
				left_trim_l = len(left_user_barcode)
				right_trim_l = None
				left_trim_r = -len(left_user_barcode)
				right_trim_r = None

			elif sum(right_params[4]):
				ran_bcode = ''.join([seq_r[x] for x in right_params[4]])
				left_trim_l = None
				right_trim_l = -len(right_user_barcode)
				left_trim_r = None
				right_trim_r = len(right_user_barcode)

			# also add the randomness from the binary part, if it exists
			if sum(left_params[2]):
				ran_bcode += seq_f[left_params[2][0]:left_params[2][1]]
			elif sum(right_params[2]):
				ran_bcode += seq_r[right_params[2][0]:right_params[2][1]]

			#Get the internal index from one side
			if sum(left_params[3]):
				index_bcode = seq_f[left_params[3][0]:left_params[3][1]]
			elif sum(right_params[3]):
				index_bcode = seq_r[right_params[3][0]:right_params[3][1]]

			#Get the binary barcode from one side, if it exists
			if sum(left_params[2]):
				binary_bcode = seq_f[left_params[2][0]:left_params[2][1]]
				list_of_A_barcodes = left_params[0]
				list_of_B_barcodes = left_params[1]
			elif sum(right_params[2]):
				binary_bcode = seq_r[right_params[2][0]:right_params[2][1]]
				list_of_A_barcodes = right_params[0]
				list_of_B_barcodes = right_params[1]
			else:
				binary_bcode = None


			name_f1 = name_f.split(' ')
			name_f1.insert(1, ":" + ran_bcode + " ")
			name_f = ''.join(name_f1)

			name_r1 = name_r.split(' ')
			name_r1.insert(1, ":" + ran_bcode + " ")
			name_r = ''.join(name_r1)

			scores = []

			for i in index_barcodes:
				scores.append(Levenshtein.hamming(i, index_bcode))

			# If binary barcode scheme has been used, split accordingly, assuming that A or B refers to biological replicates or just replicates of some sort
			if binary_bcode:
				if binary_bcode in list_of_A_barcodes: # For YY sequences
					if min(scores) <= hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 and unique
						a_barcodes_unmerged[scores.index(min(scores))] += 1
						try:
							single_reads_f_dict[current_dir+"barcode_sorted/fastq_files/A_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])])
							single_reads_r_dict[current_dir+"barcode_sorted/fastq_files/A_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])])
						except KeyError:
							single_reads_f_dict[current_dir+"barcode_sorted/fastq_files/A_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])]]
							single_reads_r_dict[current_dir+"barcode_sorted/fastq_files/A_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])]]
					else:
						try:
							bad_barcode_reads_f_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/A_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)])
							bad_barcode_reads_f_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/A_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)])
						except KeyError:
							bad_barcode_reads_f_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/A_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)]]
							bad_barcode_reads_f_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/A_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)]]
						bad_barcode += 1
						bad_set.append(index_bcode)

				elif binary_bcode in list_of_B_barcodes: # For RR sequences
					if min(scores) <= hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 and unique
						b_barcodes_unmerged[scores.index(min(scores))] +=1
						try:
							single_reads_f_dict[current_dir+"barcode_sorted/fastq_files/B_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])])
							single_reads_r_dict[current_dir+"barcode_sorted/fastq_files/B_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])])
						except KeyError:
							single_reads_f_dict[current_dir+"barcode_sorted/fastq_files/B_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])]]
							single_reads_r_dict[current_dir+"barcode_sorted/fastq_files/B_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])]]
					else:
						try:
							bad_barcode_reads_f_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/B_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)])
							bad_barcode_reads_f_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/B_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)])
						except KeyError:
							bad_barcode_reads_f_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/B_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)]]
							bad_barcode_reads_f_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/B_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)]]
						bad_set.append(index_bcode)
						bad_barcode += 1

				else:
					try:
						alien_reads_f_dict[current_dir+"bad_barcodes/unknownBinary_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)])
						alien_reads_r_dict[current_dir+"bad_barcodes/unknownBinary_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_r, seq_r, qual_r)])
					except KeyError:
						alien_reads_f_dict[current_dir+"bad_barcodes/unknownBinary_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)]]
						alien_reads_r_dict[current_dir+"bad_barcodes/unknownBinary_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_r, seq_r, qual_r)]]
					alien_count += 1

			else:
				if min(scores) <= hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 and unique
					a_barcodes_unmerged[scores.index(min(scores))] += 1
					try:
						single_reads_f_dict[current_dir+"barcode_sorted/fastq_files/" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])])
						single_reads_r_dict[current_dir+"barcode_sorted/fastq_files/" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])])
					except KeyError:
						single_reads_f_dict[current_dir+"barcode_sorted/fastq_files/" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])]]
						single_reads_r_dict[current_dir+"barcode_sorted/fastq_files/" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])]]
				else:
					try:
						alien_reads_f_dict[current_dir+"bad_barcodes/unknownBarcode_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)])
						alien_reads_r_dict[current_dir+"bad_barcodes/unknownBarcode_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_r, seq_r, qual_r)])
					except KeyError:
						alien_reads_f_dict[current_dir+"bad_barcodes/unknownBarcode_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f, qual_f)]]
						alien_reads_r_dict[current_dir+"bad_barcodes/unknownBarcode_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_r, seq_r, qual_r)]]
					bad_set.append(index_bcode)
					alien_count += 1


			# This is the routine that makes sure that the dictionary is emptied every "chunk" times
			# I guess one could write a function to find out the optimal "chunk" size by checking system recources
			# and figuring out how much one should use. I cannot do this at the moment unfortunately.

			if count1 == multiplier*chunk:
				if gzipp:
					for key, value in single_reads_f_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in single_reads_r_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in bad_barcode_reads_f_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in bad_barcode_reads_r_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in alien_reads_f_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in alien_reads_r_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

				else:
					for key, value in single_reads_f_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])


					for key, value in single_reads_r_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in bad_barcode_reads_f_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])


					for key, value in bad_barcode_reads_r_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in alien_reads_f_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in alien_reads_r_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

				multiplier += 1
				single_reads_f_dict = {}
				single_reads_r_dict = {}
				bad_barcode_reads_f_dict = {}
				bad_barcode_reads_r_dict = {}
				alien_reads_f_dict = {}
				alien_reads_r_dict = {}

		except StopIteration:
			if gzipp:
					for key, value in single_reads_f_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in single_reads_r_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in bad_barcode_reads_f_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in bad_barcode_reads_r_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in alien_reads_f_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in alien_reads_r_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

			else:
				for key, value in single_reads_f_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])


				for key, value in single_reads_r_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				for key, value in bad_barcode_reads_f_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])


				for key, value in bad_barcode_reads_r_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				for key, value in alien_reads_f_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				for key, value in alien_reads_r_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])


			bad_file = current_dir+"bad_barcodes/bad_barcodes.txt"
			with open(bad_file, 'a') as ff:
				for i in bad_set:
					ff.write(i+'\n')

			single_reads_f_dict = {}
			single_reads_r_dict = {}
			bad_barcode_reads_f_dict = {}
			bad_barcode_reads_r_dict = {}
			alien_reads_f_dict = {}
			alien_reads_r_dict = {}

			script_process['BarcodeSplittingUnmerged'] = 'Finished'
			with open(logfile, 'a') as ff:
				ff.write("** BarcodeSplittingUnmerged:Finished: %s\n"%time.strftime("%X %x"))
			break

	left_reads.close()
	right_reads.close()
	return count1, alien_count, bad_barcode, a_barcodes_unmerged, b_barcodes_unmerged, binary_bcode

def MergedSplitter(merged, check_binary_bcode, gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, left_user_barcode, right_user_barcode, hamming_distance, logfile, current_dir):
	""" Same as above, just for the merged reads.
	The barcodes are reverse complemented to make sure that A for unmerged means A for the merged reads. """


	folders = [current_dir+"bad_barcodes", current_dir+"barcode_sorted", current_dir+"barcode_sorted/fastq_files"]

	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)

	if check_binary_bcode == True:
		bad_barcodes_folders =[]

		for i in index_barcode_names:
			bad_barcodes_folders.append(current_dir+"bad_barcodes/"+i)

		for i in bad_barcodes_folders:
			if not os.path.exists(i):
				os.makedirs(i)


	printStatus("Will process merged reads now.")

	merged_reads_dict = {}
	merged_alien_reads_dict = {}
	merged_bad_barcodes_dict = {}
	alien_count = 0
	bad_barcode = 0
	a_barcodes_merged = [0 for x in range(len(index_barcode_names))]
	b_barcodes_merged = [0 for x in range(len(index_barcode_names))]
	count2 = 0
	multiplier = 1
	bad_set = []

	seqprep_merged = gzip.open(merged, 'rb')

	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N' : 'N'}

	while True:
		try:
			merged_name = seqprep_merged.next().strip()
			merged_seq = seqprep_merged.next().strip()
			merged_plus = seqprep_merged.next().strip()
			merged_qual = seqprep_merged.next().strip()

			count2 += 1

			# collect the random barcode from both ends, wherever it is defined
			left_trim = None
			right_trim = None
			if sum(left_params[4]) != 0 and sum(right_params[4]) != 0:
				ran_bcode = ''.join([merged_seq[x] for x in left_params[4]]) + ''.join(["".join(complement.get(base, base) for base in reversed(merged_seq))[x] for x in right_params[4]])
				left_trim = len(left_user_barcode)
				right_trim = -len(right_user_barcode)
			elif sum(left_params[4]):
				ran_bcode = ''.join([merged_seq[x] for x in left_params[4]])
				left_trim = len(left_user_barcode)
				right_trim = None
			elif sum(right_params[4]):
				ran_bcode = ''.join(["".join(complement.get(base, base) for base in reversed(merged_seq))[x] for x in right_params[4]])
				left_trim = None
				right_trim = -len(right_user_barcode)

			# also add the randomness from the binary part, if it exists
			if sum(left_params[2]):
				ran_bcode += merged_seq[left_params[2][0]:left_params[2][1]]
			elif sum(right_params[2]):
				ran_bcode += merged_seq[right_params[2][0]:right_params[2][1]]

			#Get the internal index from one side
			if sum(left_params[3]):
				index_bcode = merged_seq[left_params[3][0]:left_params[3][1]]
			elif sum(right_params[3]):
				index_bcode = "".join(complement.get(base, base) for base in reversed(merged_seq))[right_params[3][0]:right_params[3][1]]

			#Get the binary barcode from one side, if it exists
			if sum(left_params[2]):
				binary_bcode = merged_seq[left_params[2][0]:left_params[2][1]]
				list_of_A_barcodes = left_params[0]
				list_of_B_barcodes = left_params[1]
			elif sum(right_params[2]):
				binary_bcode = "".join(complement.get(base, base) for base in reversed(merged_seq))[right_params[2][0]:right_params[2][1]]
				list_of_A_barcodes = right_params[0]
				list_of_B_barcodes = right_params[1]
			else:
				binary_bcode = None

			merged_name_1 = merged_name.split(' ')
			merged_name_1.insert(1, ":" + ran_bcode + " ")
			merged_name = ''.join(merged_name_1)

			scores = []

			for i in index_barcodes:
				scores.append(Levenshtein.hamming(i, index_bcode))
			if binary_bcode:
				if binary_bcode in list_of_A_barcodes: # For RYYR sequences
					if min(scores) <= hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 (0 or 1) and unique
						a_barcodes_merged[scores.index(min(scores))] += 1
						try:
							merged_reads_dict[current_dir+"barcode_sorted/fastq_files/A_merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"].append(['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])])
						except KeyError:
							merged_reads_dict[current_dir+"barcode_sorted/fastq_files/A_merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"] = [['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])]]
					else:
						try:
							merged_bad_barcodes_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/A_merged_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + ".fastq"].append(["%s\n%s\n+\n%s\n" %(merged_name, merged_seq, merged_qual)])
						except KeyError:
							merged_bad_barcodes_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/A_merged_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + ".fastq"] = [["%s\n%s\n+\n%s\n" %(merged_name, merged_seq, merged_qual)]]

						bad_barcode += 1
						bad_set.append(index_bcode)

				elif binary_bcode in list_of_B_barcodes: # For YRRY sequences
					if min(scores) <= hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 (0 or 1) and unique
						b_barcodes_merged[scores.index(min(scores))] += 1
						try:
							merged_reads_dict[current_dir+"barcode_sorted/fastq_files/B_merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"].append(['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])])
						except KeyError:
							merged_reads_dict[current_dir+"barcode_sorted/fastq_files/B_merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"] = [['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])]]
					else:
						try:
							merged_bad_barcodes_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/B_merged_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + ".fastq"].append(["%s\n%s\n+\n%s\n" %(merged_name, merged_seq, merged_qual)])
						except KeyError:
							merged_bad_barcodes_dict[current_dir+"bad_barcodes/"+index_barcode_names[scores.index(min(scores))]+"/B_merged_bad_barcode_" + index_barcode_names[scores.index(min(scores))] + ".fastq"] = [["%s\n%s\n+\n%s\n" %(merged_name, merged_seq, merged_qual)]]

						bad_barcode += 1
						bad_set.append(index_bcode)

				else:
					try:
						merged_alien_reads_dict[current_dir+"bad_barcodes/unknownBinary_merged.fastq"].append(["%s\n%s\n+\n%s\n" %(merged_name, merged_seq, merged_qual)])
					except KeyError:
						merged_alien_reads_dict[current_dir+"bad_barcodes/unknownBinary_merged.fastq"] = [["%s\n%s\n+\n%s\n" %(merged_name, merged_seq, merged_qual)]]

					alien_count += 1
			else:
				if min(scores) <= hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 (0 or 1) and unique
					a_barcodes_merged[scores.index(min(scores))] += 1
					try:
						merged_reads_dict[current_dir+"barcode_sorted/fastq_files/merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"].append(['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])])
					except KeyError:
						merged_reads_dict[current_dir+"barcode_sorted/fastq_files/merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"] = [['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])]]
				else:
					try:
						merged_alien_reads_dict[current_dir+"bad_barcodes/unknownBarcode_merged.fastq"].append(["%s\n%s\n+\n%s\n" %(merged_name, merged_seq, merged_qual)])
					except KeyError:
						merged_alien_reads_dict[current_dir+"bad_barcodes/unknownBarcode_merged.fastq"] = [["%s\n%s\n+\n%s\n" %(merged_name, merged_seq, merged_qual)]]

					alien_count += 1
					bad_set.append(index_bcode)

			if count2 == multiplier*chunk:
				if gzipp:
					for key, value in merged_reads_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in merged_bad_barcodes_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in merged_alien_reads_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])
				else:
					for key, value in merged_reads_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in merged_bad_barcodes_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in merged_alien_reads_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

				multiplier += 1
				merged_reads_dict = {}
				merged_bad_barcodes_dict = {}
				merged_alien_reads_dict = {}

		except StopIteration:
				if gzipp:
					for key, value in merged_reads_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in merged_bad_barcodes_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in merged_alien_reads_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])
				else:
					for key, value in merged_reads_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in merged_bad_barcodes_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

					for key, value in merged_alien_reads_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])

				bad_file_merged = current_dir+"bad_barcodes/bad_barcodes_merged.txt"
				with open(bad_file_merged, 'a') as ff:
					for i in bad_set:
						ff.write(i+'\n')

				merged_reads_dict = {}
				merged_alien_reads_dict = {}
				merged_bad_barcodes_dict = {}

				script_process['BarcodeSplittingMerged'] = 'Finished'
				with open(logfile, 'a') as ff:
					ff.write("** BarcodeSplittingMerged:Finished: %s\n"%time.strftime("%X %x"))
				break

	seqprep_merged.close()
	return count2, alien_count, bad_barcode, a_barcodes_merged, b_barcodes_merged, binary_bcode

def runBowtie2(merged_files, script_process, indices, mapping, processor, ignore, logfile, current_dir):

	folders = [current_dir+"barcode_sorted/bowtie2",
				current_dir+"barcode_sorted/bowtie2/bam_files",
				current_dir+"barcode_sorted/bowtie2/log",
				current_dir+"barcode_sorted/bowtie2/cleaned",
				current_dir+"barcode_sorted/bowtie2/unmapped"]

	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)

	counting_files = 1
	total_files = len(merged_files)

	for i in merged_files:
		if merged_files[0].split('/')[-1].split('merged')[0] in ['A_', 'B_']:
			if merged_files[0].split(".fastq")[1] == ".gz": #If the files are gzipped to save space.
				both_files = i + "," + i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]+"_R1.fastq.gz"
			else:
				both_files = i + "," + i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]+"_R1.fastq"

			new_name = i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]
		else:
			if merged_files[0].split(".fastq")[1] == ".gz": #If the files are gzipped to save space.
				both_files = i + "," + i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq.gz"
				fwd_fastq = i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq.gz"
				rvs_fastq = i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R2.fastq.gz"
			else:
				both_files = i + "," + i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq"
				fwd_fastq = i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq"
				rvs_fastq = i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R2.fastq"
			new_name = i.split("merged_")[1].split(".fastq")[0]


		cmd = "%s -p %s -x %s --un-gz %s -U %s 2> %s | %s view -b - | %s sort - %s" %(bowtie2_path, processor, indices[mapping][0],
			current_dir+"barcode_sorted/bowtie2/unmapped/"+ new_name.split("/")[-1] +"_unmapped.fastq.gz",
			both_files, current_dir+"barcode_sorted/bowtie2/log/"+ new_name.split("/")[-1] + ".bowtie2.log", samtools_path, samtools_path, current_dir+"barcode_sorted/bowtie2/bam_files/"+ new_name.split("/")[-1])

		cmd2 = "bowtie2 -p %s -x %s --un-gz %s -U %s 2> %s | samtools view -b - | samtools sort - %s" %(processor, indices[mapping][0],
			current_dir+"barcode_sorted/bowtie2/unmapped/"+ new_name.split("/")[-1] +"_unmapped.fastq.gz",
			both_files, current_dir+"barcode_sorted/bowtie2/log/"+ new_name.split("/")[-1] + ".bowtie2.log", current_dir+"barcode_sorted/bowtie2/bam_files/"+ new_name.split("/")[-1])

		printStatus("Running bowtie2 for file (%s of %s): %s " %(counting_files, total_files, i.split("merged_")[-1].split(".fastq.gz")[0]))

		counting_files += 1

		bowtie_logfile = current_dir+"barcode_sorted/bowtie2/log/"+ new_name.split("/")[-1] + ".bowtie2.log"

		with open(bowtie_logfile, 'a') as ff:
			ff.write("%s\n\n" %cmd2)

		subprocess.check_call(cmd, shell=True)

	mapped_bams = glob.glob(current_dir+"barcode_sorted/bowtie2/bam_files/*.bam")

	for i in mapped_bams:
		cmd = "samtools index %s" %i
		subprocess.check_call(cmd, shell=True)


	script_process['bowtie2'] = 'Finished'
	with open(logfile, 'a') as ff:
		ff.write("** bowtie2:Finished: %s\n"%time.strftime("%X %x"))


def runBBMap(unmapped_reads, script_process, indices, mapping, merged_files, processor, logfile, current_dir):

	folders = [current_dir+"barcode_sorted/bbmap",
				current_dir+"barcode_sorted/bbmap/bam_files",
				current_dir+"barcode_sorted/bbmap/log",
				current_dir+"barcode_sorted/bbmap/cleaned"]

	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)

	counting_files2 = 1
	total_files = len(merged_files)

	# Some parameters for BBMap and BBMerge. usejni=t only after compiling the C code!
	BBMAP_OPTIONS = "usejni=f vslow=t maxindel=100k sam=1.3"

	if mapping == "dm6_simple" or mapping == "dvir":
		BBMAP_PATH = "/data/akhtar/group/Giuseppe/bbmap"
		bbmap_path_sh = "/data/akhtar/group/Giuseppe/bbmap/bbmap.sh"
	else:
		BBMAP_PATH = "/data/akhtar/group/Ibrahim_data_2014/__BBMap__/bbmap"
		bbmap_path_sh = "/data/akhtar/group/Ibrahim_data_2014/__BBMap__/bbmap/bbmap.sh"


	for i in unmapped_reads:
		cmd = "cd %s && %s %s threads=%s build=%s in=%s out=%s 2> %s" %(BBMAP_PATH, bbmap_path_sh,
																	BBMAP_OPTIONS, processor,
																	indices[mapping][1],
																	i,
																	current_dir+"barcode_sorted/bbmap/bam_files/"+i.split("/")[-1].split("_unmapped.fastq")[0]+".bam",
																	current_dir+"barcode_sorted/bbmap/log/"+i.split("/")[-1].split("_unmapped.fastq")[0]+".bbmap.log")

		printStatus("Running BBMap on file (%s of %s): %s" %(counting_files2, total_files, i.split("/")[-1]))
		#print "\n", cmd, "\n"
		counting_files2 += 1
		subprocess.check_call(cmd, shell=True)

	script_process['bbmap'] = 'Finished'
	with open(logfile, 'a') as ff:
		ff.write("** bbmap:Finished: %s\n"%time.strftime("%X %x"))

def runHiSat(merged_files, script_process, indices, mapping, processor, logfile, current_dir):

	folders = [current_dir+"barcode_sorted/hisat",
				current_dir+"barcode_sorted/hisat/bam_files",
				current_dir+"barcode_sorted/hisat/log",
				current_dir+"barcode_sorted/hisat/unmapped",
				current_dir+"barcode_sorted/hisat/cleaned"]

	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)



	counting_files = 1
	total_files = len(merged_files)

	for i in merged_files:
		if merged_files[0].split('/')[-1].split('merged')[0] in ['A_', 'B_']:
			if merged_files[0].split(".fastq")[1] == ".gz": #If the files are gzipped to save space.
				both_files = i + "," + i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]+"_R1.fastq.gz"
			else:
				both_files = i + "," + i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]+"_R1.fastq"

			new_name = i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]
		else:
			if merged_files[0].split(".fastq")[1] == ".gz": #If the files are gzipped to save space.
				both_files = i + "," + i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq.gz"
				fwd_fastq = i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq.gz"
				rvs_fastq = i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R2.fastq.gz"
			else:
				both_files = i + "," + i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq"
				fwd_fastq = i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq"
				rvs_fastq = i.split("merged_")[0] + i.split("merged_")[1].split(".fastq")[0]+"_R2.fastq"
			new_name = i.split("merged_")[1].split(".fastq")[0]




		# Options to add/change:
		# --pen-cansplice <int>              penalty for a canonical splice site (0)
		# --pen-noncansplice <int>           penalty for a non-canonical splice site (12)
		# --min-intronlen <int>              minimum intron length (20)
		# --max-intronlen <int>              maximum intron length (500000)
		# --rna-strandness <string>          Specify strand-specific information (unstranded)
		# --ma <int>         match bonus (0 for --end-to-end, 2 for --local)
		# --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)
		# --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
		# --rdg <int>,<int>  read gap open, extend penalties (5,3)

		mammals = ["hg19", "hg19_simple", "hg38", "mm9", "mm10"]


		if mapping in mammals:
			hisat_options = "--rna-strandness F --max-intronlen 100000 --pen-cansplice 3"
		else:
			hisat_options = "--rna-strandness F --max-intronlen 10000 --pen-cansplice 3"


		cmd = "%s %s -p %s -x %s -U %s --un-gz %s 2>> %s | %s view -b -F 256 - | %s sort - %s" %(hisat2_path, hisat_options, processor, indices[mapping][3],
				both_files, current_dir+"barcode_sorted/hisat/unmapped/"+ new_name.split("/")[-1] +"_unmapped.fastq.gz",
				current_dir+"barcode_sorted/hisat/log/"+ new_name.split("/")[-1] + ".hisat.log", samtools_path, samtools_path, current_dir+"barcode_sorted/hisat/bam_files/"+ new_name)

		cmd2 = "hisat2 %s -p %s -x %s -U %s --un-gz %s 2>> %s | samtools view -b -F 256 - | samtools sort - %s" %(hisat_options, processor, indices[mapping][3],
				both_files, current_dir+"barcode_sorted/hisat/unmapped/"+ new_name.split("/")[-1] +"_unmapped.fastq.gz",
				current_dir+"barcode_sorted/hisat/log/"+ new_name.split("/")[-1] + ".hisat.log", current_dir+"barcode_sorted/hisat/bam_files/"+ new_name)

		printStatus("Running HISAT2 for file (%s of %s): %s " %(counting_files, total_files, new_name))

		counting_files += 1

		hisat_logfile = current_dir+"barcode_sorted/hisat/log/"+ new_name.split("/")[-1] + ".hisat.log"

		with open(hisat_logfile, 'a') as ff:
			ff.write("%s\n\n" %cmd2)

		subprocess.check_call(cmd, shell=True)

	mapped_bams = glob.glob(current_dir+"barcode_sorted/hisat/bam_files/*.bam")

	for i in mapped_bams:
		cmd = "samtools index %s" %i
		subprocess.check_call(cmd, shell=True)

	script_process['hisat'] = 'Finished'
	with open(logfile, 'a') as ff:
		ff.write("** hisat:Finished: %s\n"%time.strftime("%X %x"))

def spikeMapper(unmapped_reads, script_process, indices, mapping, processor, merged_files, logfile, current_dir, spike_type):


	if spike_type == "fetish":
		spikepath = indices[mapping][4]+"FETISH_spikes/drosophila/fetish_spikes"
	elif spike_type == "mapcap":
		spikepath = indices[mapping][4]+"MAPCAP_spikes/drosophila/mapcap_spikes"
	elif spike_type == "srmc":
		spikepath = indices[mapping][4]+"smallRNA-MAPCAP_spikes/drosophila/srmc_spikes"


	folders = [current_dir+"spikes",
				current_dir+"spikes/cleaned"]

	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)

	counting_files = 1
	total_files = len(merged_files)

	for i in unmapped_reads:

		cmd = "%s -p %s -x %s -U %s 2> %s | %s view -b - | %s sort - %s" %(bowtie2_path, processor, spikepath,
			i, current_dir+"spikes/"+i.split("/")[-1].split("_unmapped.fastq")[0]+"_spikes.log",
			samtools_path, samtools_path, current_dir+"spikes/"+i.split("/")[-1].split("_unmapped.fastq")[0]+"_spikes")

		cmd2 = "bowtie2 -p %s -x %s -U %s 2>> %s | samtools view -b - | samtools sort - %s" %(processor, spikepath,
			i, current_dir+"spikes/"+i.split("/")[-1].split("_unmapped.fastq")[0]+"_spikes.log",
			current_dir+"spikes/"+i.split("/")[-1].split("_unmapped.fastq")[0]+"_spikes")

		printStatus("Running Bowtie2 on Spikes for file (%s of %s): %s " %(counting_files, total_files, i.split("/")[-1]))

		counting_files += 1

		bowtie_logfile = current_dir+"spikes/"+i.split("/")[-1].split("_unmapped.fastq")[0]+"_spikes.log"

		with open(bowtie_logfile, 'a') as ff:
			ff.write("%s\n\n" %cmd2)

		subprocess.check_call(cmd, shell=True)

	script_process['spikesMapper'] = 'Finished'
	with open(logfile, 'a') as ff:
		ff.write("** spikesMapper:Finished: %s\n"%time.strftime("%X %x"))


def runSamBedCleaner(bowtie_bams, bbmap_bams, hisat_bams, script_process, chunk, logfile, current_dir, keep_BAM):

	folders = [current_dir+"final_files",
				current_dir+"final_files/bams"]

	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)

	if keep_BAM:
		from bamcleaner_nonMulti_splitter import DupRemover as duprem
	else:
		from bamcleaner_nonMulti import DupRemover as duprem

	counting_files = 1
	counting_files2 = 1
	#### sambedcleaner options:
	keep_bed = False
	min_quality = 10 # Reads with quality less than this will be discarded.
	multi = False # This is marker that tells the program that all the .bam files will be processed at the same time.

	if bowtie_bams and bbmap_bams:
		total_files = len(bowtie_bams)
		total_files2 = len(bbmap_bams)

		for i in bowtie_bams:
			newfolder = current_dir+"barcode_sorted/bowtie2/cleaned/"+i.split("/")[-1].split(".bam")[0]
			if not os.path.exists(newfolder):
				os.makedirs(newfolder)

			printStatus("Running BamCleaner for Bowtie2 .bam file (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
			duprem(i, newfolder+"/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi, current_dir, None)

			counting_files+=1

		script_process['sambedcleaner-bt2'] = 'Finished'
		with open(logfile, 'a') as ff:
			ff.write("** sambedcleaner-bt2:Finished: %s\n"%time.strftime("%X %x"))

		for i in bbmap_bams:
			newfolder = current_dir+"barcode_sorted/bbmap/cleaned/"+i.split("/")[-1].split(".bam")[0]
			if not os.path.exists(newfolder):
				os.makedirs(newfolder)

			printStatus("Running BamCleaner for BBMap .bam file (%s of %s): %s" %(counting_files2, total_files2, i.split("/")[-1]))
			duprem(i, newfolder+"/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi, current_dir, None)

			counting_files2+=1

		script_process['sambedcleaner-bb'] = 'Finished'
		with open(logfile, 'a') as ff:
			ff.write("** sambedcleaner-bb:Finished: %s\n"%time.strftime("%X %x"))

	elif bowtie_bams:
		total_files = len(bowtie_bams)

		for i in bowtie_bams:
			newfolder = current_dir+"barcode_sorted/bowtie2/cleaned/"+i.split("/")[-1].split(".bam")[0]
			if not os.path.exists(newfolder):
				os.makedirs(newfolder)

			printStatus("Running BamCleaner for Bowtie2 .bam file (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
			duprem(i, newfolder+"/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi, current_dir, None)

			counting_files+=1

		script_process['sambedcleaner'] = 'Finished'
		with open(logfile, 'a') as ff:
			ff.write("** sambedcleaner:Finished: %s\n"%time.strftime("%X %x"))

	elif bbmap_bams:
		total_files2 = len(bbmap_bams)

		for i in bbmap_bams:
			newfolder = current_dir+"barcode_sorted/bbmap/cleaned/"+i.split("/")[-1].split(".bam")[0]
			if not os.path.exists(newfolder):
				os.makedirs(newfolder)

			printStatus("Running BamCleaner for BBMap .bam file (%s of %s): %s" %(counting_files2, total_files2, i.split("/")[-1]))
			duprem(i, newfolder+"/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi, current_dir, None)

			counting_files2+=1

		script_process['sambedcleaner'] = 'Finished'
		with open(logfile, 'a') as ff:
			ff.write("** sambedcleaner:Finished: %s\n"%time.strftime("%X %x"))

	elif hisat_bams:
		total_files = len(hisat_bams)

		for i in hisat_bams:
			newfolder = current_dir+"barcode_sorted/hisat/cleaned/"+i.split("/")[-1].split(".bam")[0]
			if not os.path.exists(newfolder):
				os.makedirs(newfolder)

			printStatus("Running BamCleaner for HISAT2 .bam file (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
			duprem(i, newfolder+"/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi, current_dir, None)

			counting_files+=1

		script_process['sambedcleaner'] = 'Finished'
		with open(logfile, 'a') as ff:
			ff.write("** sambedcleaner:Finished: %s\n"%time.strftime("%X %x"))

def runSpikeCleaner(spike_bams, script_process, chunk, logfile, current_dir):
	keep_bed = False
	min_quality = 10 # Reads with quality less than this will be discarded.
	multi = False # This is marker that tells the program that all the .bam files will be processed at the same time.

	counting_files = 1

	if spike_bams:
		total_files = len(spike_bams)

		for i in spike_bams:
			duprem(i, current_dir+"spikes/cleaned/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi, current_dir, "count_spikes")
			printStatus("Running SpikeCleaner for HISAT2 .bam file (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))

			counting_files+=1

	script_process['spikesBamCleaner'] = 'Finished'
	with open(logfile, 'a') as ff:
		ff.write("** spikesBamCleaner:Finished: %s\n"%time.strftime("%X %x"))

def BamMerger(cleaned_files, script_process, logfile, current_dir):
	""" This just merges the bowtie2 and bbmap generated .bams to get the best of both worlds. """

	printStatus("Merging and sorting and indexing bowtie2 and BBMap outputs.")

	for i in cleaned_files:
		cleaned_files_prefix = i.split("/")[-1].split("_cleaned")[0]
		cmd = "%s merge -h %s %s %s %s" %(samtools_path, i, current_dir+"final_files/bams/temp-"+cleaned_files_prefix+"_cleaned+merged.bam", i, current_dir+"barcode_sorted/bbmap/cleaned/"+cleaned_files_prefix+"/"+i.split("/")[-1])
		subprocess.check_call(cmd, shell=True)

	merged_bams = glob.glob(current_dir+"final_files/bams/*.bam")

	for i in merged_bams:
		cmd = "%s sort %s %s" %(samtools_path, i, current_dir+"final_files/bams/"+i.split("temp-")[-1].split(".bam")[0])
		printStatus("Sorting .bam file: %s"%i.split("/")[-1])
		subprocess.check_call(cmd, shell=True)

	for i in merged_bams:
		os.remove(i)

	merged_bams2 = glob.glob(current_dir+"final_files/bams/*.bam")

	for i in merged_bams2:
		cmd = "%s index %s" %(samtools_path, i)
		subprocess.check_call(cmd, shell=True)

	script_process['MergeBowtie2BBMap'] = 'Finished'

	with open(logfile, 'a') as ff:
		ff.write("** MergeBowtie2BBMap:Finished: %s\n"%time.strftime("%X %x"))

def RepMerger(index_barcode_names, script_process, indices, mapping, aligner, logfile, current_dir, processor):

	""" With too many profiles at hand, it doesn't make sense to look at replicates all the time.
	This function merges A and B replicates into one and calculates a bunch of coverage files using
	genomeCoverageBed and converts them into bigwigs using bedGraphToBigWig."""

	folders = [current_dir+"final_files/replicates_merged/bams",
	current_dir+"final_files/replicates_merged/bigwigs",
	current_dir+"final_files/replicates_merged/bigwigs/raw",
	current_dir+"final_files/replicates_merged/bigwigs/CPM_normalized",
	current_dir+"final_files/replicates_merged/p5end",
	current_dir+"final_files/replicates_merged/p5end/raw",
	current_dir+"final_files/replicates_merged/p5end/CPM_normalized",
	current_dir+"final_files/replicates_merged/p3end",
	current_dir+"final_files/replicates_merged/p3end/raw",
	current_dir+"final_files/replicates_merged/p3end/CPM_normalized"]


	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)


	printStatus("Merging biological replicates.")
	if aligner == 'bt2-bbmap':
		for i in index_barcode_names:
			name1 = "%sfinal_files/bams/A_%s_cleaned+merged.bam" %(current_dir, i)
			name2 = "%sfinal_files/bams/B_%s_cleaned+merged.bam" %(current_dir, i)
			out = "%sfinal_files/replicates_merged/bams/%s_cleaned_mergedReps.bam" %(current_dir, i)
			cmd = "%s merge -h %s %s %s %s" %(samtools_path, name1, out, name1, name2)
			printStatus("Merging replicates of sample:%s"%i)
			subprocess.check_call(cmd, shell=True)
	else:
		for i in index_barcode_names:
			name1 = "%sfinal_files/bams/A_%s_cleaned.bam" %(current_dir, i)
			name2 = "%sfinal_files/bams/B_%s_cleaned.bam" %(current_dir, i)
			out = "%sfinal_files/replicates_merged/bams/%s_cleaned_mergedReps.bam" %(current_dir, i)
			cmd = "%s merge -h %s %s %s %s" %(samtools_path, name1, out, name1, name2)
			printStatus("Merging replicates of sample:%s"%i)
			subprocess.check_call(cmd, shell=True)

	merged_reps_bams = glob.glob(current_dir+"final_files/replicates_merged/bams/*.bam")

	printStatus("Generating indices for those merged bams.")

	for i in merged_reps_bams:
		cmd = "%s index %s" %(samtools_path, i)
		subprocess.check_call(cmd, shell=True)

	printStatus("Generating coverage files for those merged bams.")

	counting_files = 1
	for i in merged_reps_bams:
		total_files = len(merged_reps_bams)
		cmd1 = "%s -b %s -o %s -bs 1 -p %s --skipNAs  2>> %s" %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/bigwigs/raw/"+i.split("/")[-1].split(".bam")[0]+".bw", processor, current_dir+"final_files/replicates_merged/bigwigs/raw/bamCoverage.log")
		cmd2 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --filterRNAstrand reverse 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/bigwigs/raw/"+i.split("/")[-1].split(".bam")[0]+"_plus.bw", processor, current_dir+"final_files/replicates_merged/bigwigs/raw/bamCoverage.log", current_dir+"final_files/replicates_merged/bigwigs/raw/bamCoverage.log")
		cmd3 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/bigwigs/raw/"+i.split("/")[-1].split(".bam")[0]+"_minus.bw", processor, current_dir+"final_files/bigwigs/raw/bamCoverage.log", current_dir+"final_files/bigwigs/raw/bamCoverage.log")

		cmd4 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM 2>> %s " %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/bigwigs/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_CPM.bw", processor, current_dir+"final_files/bigwigs/CPM_normalized/bamCoverage.log")
		cmd5 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --filterRNAstrand reverse 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/bigwigs/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_plus_CPM.bw", processor, current_dir+"final_files/replicates_merged/bigwigs/CPM_normalized/bamCoverage.log", current_dir+"final_files/replicates_merged/bigwigs/CPM_normalized/bamCoverage.log")
		cmd6 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/bigwigs/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_minus_CPM.bw", processor, current_dir+"final_files/replicates_merged/bigwigs/CPM_normalized/bamCoverage.log", current_dir+"final_files/replicates_merged/bigwigs/CPM_normalized/bamCoverage.log")

		printStatus("Calculating raw and CPM normalized coverage for (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd1, shell=True)
		subprocess.check_call(cmd4, shell=True)
		printStatus("Calculating raw and CPM normalized coverage for the plus strand")
		subprocess.check_call(cmd2, shell=True)
		subprocess.check_call(cmd5, shell=True)
		printStatus("Calculating raw and CPM normalized coverage for the minus strand")
		subprocess.check_call(cmd3, shell=True)
		subprocess.check_call(cmd6, shell=True)

		counting_files+=1


	# 5'end coverage
	# bamCoverage uses std-illumina library prep parameters as default (fr-secondstrand) since our data is fr-firststrand to get the correct stand and filename I had to swap the --filterRNAstrand option
	counting_files = 1
	for i in merged_reps_bams:
		total_files = len(merged_reps_bams)
		cmd1 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset 1 2>> %s " %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p5end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p5end.bw", processor, current_dir+"final_files/replicates_merged/p5end/raw/bamCoverage.log")
		cmd2 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset 1 --filterRNAstrand reverse 2>> %s 1>> %s " %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p5end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p5end_plus_CPM.bw", processor, current_dir+"final_files/replicates_merged/p5end/raw/bamCoverage.log", current_dir+"final_files/replicates_merged/p5end/raw/bamCoverage.log")
		cmd3 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset 1 --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p5end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p5end_minus_CPM.bw", processor, current_dir+"final_files/replicates_merged/p5end/raw/bamCoverage.log", current_dir+"final_files/replicates_merged/p5end/raw/bamCoverage.log")

		cmd4 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset 1 2>> %s " %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p5end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p5end_CPM.bw", processor, current_dir+"final_files/replicates_merged/p5end/CPM_normalized/bamCoverage.log")
		cmd5 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset 1 --filterRNAstrand reverse 2>> %s 1>> %s " %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p5end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p5end_plus_CPM.bw", processor, current_dir+"final_files/replicates_merged/p5end/CPM_normalized/bamCoverage.log", current_dir+"final_files/replicates_merged/p5end/CPM_normalized/bamCoverage.log")
		cmd6 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset 1 --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p5end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p5end_minus_CPM.bw", processor, current_dir+"final_files/replicates_merged/p5end/CPM_normalized/bamCoverage.log", current_dir+"final_files/replicates_merged/p5end/CPM_normalized/bamCoverage.log")

		printStatus("Calculating raw and CPM normalized xlink coverage for (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd1, shell=True)
		subprocess.check_call(cmd4, shell=True)
		printStatus("Calculating raw and CPM normalized xlink coverage for the plus strand (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd2, shell=True)
		subprocess.check_call(cmd5, shell=True)
		printStatus("Calculating raw and CPM normalized xlink coverage for the minus strand (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd3, shell=True)
		subprocess.check_call(cmd6, shell=True)

		counting_files+=1


	# 3'end coverage
	# bamCoverage uses std-illumina library prep parameters as default (fr-secondstrand) since our data is fr-firststrand to get the correct stand and filename I had to swap the --filterRNAstrand option
	counting_files = 1
	for i in merged_reps_bams:
		total_files = len(merged_reps_bams)
		cmd1 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset -1 2>> %s " %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p3end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p3end.bw", processor, current_dir+"final_files/replicates_merged/p3end/raw/bamCoverage.log")
		cmd2 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset -1 --filterRNAstrand reverse 2>> %s 1>> %s " %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p3end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p3end_plus_CPM.bw", processor, current_dir+"final_files/replicates_merged/p3end/raw/bamCoverage.log", current_dir+"final_files/replicates_merged/p3end/raw/bamCoverage.log")
		cmd3 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset -1 --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p3end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p3end_minus_CPM.bw", processor, current_dir+"final_files/replicates_merged/p3end/raw/bamCoverage.log", current_dir+"final_files/replicates_merged/p3end/raw/bamCoverage.log")

		cmd4 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset -1 2>> %s " %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p3end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p3end_CPM.bw", processor, current_dir+"final_files/replicates_merged/p3end/CPM_normalized/bamCoverage.log")
		cmd5 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset -1 --filterRNAstrand reverse 2>> %s 1>> %s " %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p3end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p3end_plus_CPM.bw", processor, current_dir+"final_files/replicates_merged/p3end/CPM_normalized/bamCoverage.log", current_dir+"final_files/replicates_merged/p3end/CPM_normalized/bamCoverage.log")
		cmd6 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset -1 --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/replicates_merged/p3end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p3end_minus_CPM.bw", processor, current_dir+"final_files/replicates_merged/p3end/CPM_normalized/bamCoverage.log", current_dir+"final_files/replicates_merged/p3end/CPM_normalized/bamCoverage.log")

		printStatus("Calculating raw and CPM normalized xlink coverage for (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd1, shell=True)
		subprocess.check_call(cmd4, shell=True)
		printStatus("Calculating raw and CPM normalized xlink coverage for the plus strand (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd2, shell=True)
		subprocess.check_call(cmd5, shell=True)
		printStatus("Calculating raw and CPM normalized xlink coverage for the minus strand (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd3, shell=True)
		subprocess.check_call(cmd6, shell=True)

		counting_files+=1



	script_process['MergeReps'] = 'Finished'
	with open(logfile, 'a') as ff:
		ff.write("** MergeReps:Finished: %s\n"%time.strftime("%X %x"))

def GenerateCoverage(script_process, indices, mapping, logfile, current_dir, processor):

	folders = [current_dir+"final_files/bigwigs",
				current_dir+"final_files/bigwigs/raw",
				current_dir+"final_files/bigwigs/CPM_normalized",
				current_dir+"final_files/p5end",
				current_dir+"final_files/p5end/raw",
				current_dir+"final_files/p5end/CPM_normalized",
				current_dir+"final_files/p3end",
				current_dir+"final_files/p3end/raw",
				current_dir+"final_files/p3end/CPM_normalized"]


	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)


	merged_bams2 = glob.glob(current_dir+"final_files/bams/*.bam")

	counting_files = 1
	for i in merged_bams2:
		total_files = len(merged_bams2)
		cmd1 = "%s -b %s -o %s -bs 1 -p %s --skipNAs  2>> %s" %(bamCoverage_path, i, current_dir+"final_files/bigwigs/raw/"+i.split("/")[-1].split(".bam")[0]+".bw", processor, current_dir+"final_files/bigwigs/raw/bamCoverage.log")
		cmd2 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --filterRNAstrand reverse 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/bigwigs/raw/"+i.split("/")[-1].split(".bam")[0]+"_plus.bw", processor, current_dir+"final_files/bigwigs/raw/bamCoverage.log", current_dir+"final_files/bigwigs/raw/bamCoverage.log")
		cmd3 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/bigwigs/raw/"+i.split("/")[-1].split(".bam")[0]+"_minus.bw", processor, current_dir+"final_files/bigwigs/raw/bamCoverage.log", current_dir+"final_files/bigwigs/raw/bamCoverage.log")

		cmd4 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM 2>> %s " %(bamCoverage_path, i, current_dir+"final_files/bigwigs/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_CPM.bw", processor, current_dir+"final_files/bigwigs/CPM_normalized/bamCoverage.log")
		cmd5 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --filterRNAstrand reverse 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/bigwigs/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_plus_CPM.bw", processor, current_dir+"final_files/bigwigs/CPM_normalized/bamCoverage.log", current_dir+"final_files/bigwigs/CPM_normalized/bamCoverage.log")
		cmd6 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/bigwigs/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_minus_CPM.bw", processor, current_dir+"final_files/bigwigs/CPM_normalized/bamCoverage.log", current_dir+"final_files/bigwigs/CPM_normalized/bamCoverage.log")

		printStatus("Calculating raw and CPM normalized coverage for (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd1, shell=True)
		subprocess.check_call(cmd4, shell=True)
		printStatus("Calculating raw and CPM normalized coverage for the plus strand")
		subprocess.check_call(cmd2, shell=True)
		subprocess.check_call(cmd5, shell=True)
		printStatus("Calculating raw and CPM normalized coverage for the minus strand")
		subprocess.check_call(cmd3, shell=True)
		subprocess.check_call(cmd6, shell=True)

		counting_files+=1


	#5'end coverage
	# bamCoverage uses std-illumina library prep parameters as default (fr-secondstrand) since our data is fr-firststrand to get the correct stand and filename I had to swap the --filterRNAstrand option
	counting_files = 1
	for i in merged_bams2:
		total_files = len(merged_bams2)
		cmd1 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset 1 2>> %s " %(bamCoverage_path, i, current_dir+"final_files/p5end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p5end.bw", processor, current_dir+"final_files/p5end/raw/bamCoverage.log")
		cmd2 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset 1 --filterRNAstrand reverse 2>> %s 1>> %s " %(bamCoverage_path, i, current_dir+"final_files/p5end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p5end_plus_CPM.bw", processor, current_dir+"final_files/p5end/raw/bamCoverage.log", current_dir+"final_files/p5end/raw/bamCoverage.log")
		cmd3 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset 1 --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/p5end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p5end_minus_CPM.bw", processor, current_dir+"final_files/p5end/raw/bamCoverage.log", current_dir+"final_files/p5end/raw/bamCoverage.log")

		cmd4 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset 1 2>> %s " %(bamCoverage_path, i, current_dir+"final_files/p5end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p5end_CPM.bw", processor, current_dir+"final_files/p5end/CPM_normalized/bamCoverage.log")
		cmd5 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset 1 --filterRNAstrand reverse 2>> %s 1>> %s " %(bamCoverage_path, i, current_dir+"final_files/p5end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p5end_plus_CPM.bw", processor, current_dir+"final_files/p5end/CPM_normalized/bamCoverage.log", current_dir+"final_files/p5end/CPM_normalized/bamCoverage.log")
		cmd6 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset 1 --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/p5end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p5end_minus_CPM.bw", processor, current_dir+"final_files/p5end/CPM_normalized/bamCoverage.log", current_dir+"final_files/p5end/CPM_normalized/bamCoverage.log")

		printStatus("Calculating raw and CPM normalized 5' end coverage for (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd1, shell=True)
		subprocess.check_call(cmd4, shell=True)
		printStatus("Calculating raw and CPM normalized 5' end coverage for the plus strand (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd2, shell=True)
		subprocess.check_call(cmd5, shell=True)
		printStatus("Calculating raw and CPM normalized 5' end coverage for the minus strand (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd3, shell=True)
		subprocess.check_call(cmd6, shell=True)

		counting_files+=1

	#3'end coverage
	# bamCoverage uses std-illumina library prep parameters as default (fr-secondstrand) since our data is fr-firststrand to get the correct stand and filename I had to swap the --filterRNAstrand option
	counting_files = 1
	for i in merged_bams2:
		total_files = len(merged_bams2)
		cmd1 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset -1 2>> %s " %(bamCoverage_path, i, current_dir+"final_files/p3end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p3end.bw", processor, current_dir+"final_files/p3end/raw/bamCoverage.log")
		cmd2 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset -1 --filterRNAstrand reverse 2>> %s 1>> %s " %(bamCoverage_path, i, current_dir+"final_files/p3end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p3end_plus_CPM.bw", processor, current_dir+"final_files/p3end/raw/bamCoverage.log", current_dir+"final_files/p3end/raw/bamCoverage.log")
		cmd3 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --Offset -1 --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/p3end/raw/"+i.split("/")[-1].split(".bam")[0]+"_p3end_minus_CPM.bw", processor, current_dir+"final_files/p3end/raw/bamCoverage.log", current_dir+"final_files/p3end/raw/bamCoverage.log")

		cmd4 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset -1 2>> %s " %(bamCoverage_path, i, current_dir+"final_files/p3end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p3end_CPM.bw", processor, current_dir+"final_files/p3end/CPM_normalized/bamCoverage.log")
		cmd5 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset -1 --filterRNAstrand reverse 2>> %s 1>> %s " %(bamCoverage_path, i, current_dir+"final_files/p3end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p3end_plus_CPM.bw", processor, current_dir+"final_files/p3end/CPM_normalized/bamCoverage.log", current_dir+"final_files/p3end/CPM_normalized/bamCoverage.log")
		cmd6 = "%s -b %s -o %s -bs 1 -p %s --skipNAs --normalizeUsing CPM --Offset -1 --filterRNAstrand forward 2>> %s 1>> %s" %(bamCoverage_path, i, current_dir+"final_files/p3end/CPM_normalized/"+i.split("/")[-1].split(".bam")[0]+"_p3end_minus_CPM.bw", processor, current_dir+"final_files/p3end/CPM_normalized/bamCoverage.log", current_dir+"final_files/p3end/CPM_normalized/bamCoverage.log")

		printStatus("Calculating raw and CPM normalized 3' end coverage for (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd1, shell=True)
		subprocess.check_call(cmd4, shell=True)
		printStatus("Calculating raw and CPM normalized 3' end coverage for the plus strand (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd2, shell=True)
		subprocess.check_call(cmd5, shell=True)
		printStatus("Calculating raw and CPM normalized 3' end coverage for the minus strand (%s of %s): %s" %(counting_files, total_files, i.split("/")[-1]))
		subprocess.check_call(cmd3, shell=True)
		subprocess.check_call(cmd6, shell=True)

		counting_files+=1

	script_process['GenerateCov'] = 'Finished'
	with open(logfile, 'a') as ff:
		ff.write("** GenerateCov:Finished: %s\n"%time.strftime("%X %x"))

def plotEnrichment(script_process, indices, mapping, logfile, current_dir, processor, aligner):


	folders = current_dir+"final_files/EnrichmentPlots"

	if not os.path.exists(folders):
		os.makedirs(folders)


	drosophila = ["dm6", "dm6_simple", "dm3"]

	if mapping in drosophila:
		beds = "/data/akhtar/group/Giuseppe/INDEXES/annotations_for_MAPCap/TSS_0-100bp_dedup_reallyNOabundantRNAs.bed /data/akhtar/group/Giuseppe/INDEXES/annotations_for_MAPCap/dm6_CRs.bed /data/akhtar/group/Giuseppe/INDEXES/annotations_for_MAPCap/dm6_intergenics.bed /data/akhtar/group/Giuseppe/INDEXES/annotations_for_MAPCap/dm6_rRNAs.bed /data/akhtar/group/Giuseppe/INDEXES/annotations_for_MAPCap/dm6_snRNAs.bed /data/akhtar/group/Giuseppe/INDEXES/annotations_for_MAPCap/dm6_snoRNAs.bed /data/akhtar/group/Giuseppe/INDEXES/annotations_for_MAPCap/genes.bed"
		abundant_beds = "/data/akhtar/group/Giuseppe/INDEXES/annotations_for_MAPCap/individual_abundantRNAs/*.bed"
		rRNA_beds = "/data/akhtar/group/Giuseppe/INDEXES/annotations_for_MAPCap/individual_rRNAs/*.bed"

		input_bams = glob.glob(current_dir+"final_files/bams/*.bam")
		counting_files = 1
		for i in input_bams:
			if aligner == "bt2":
				both_files = i + " " + current_dir+"barcode_sorted/bowtie2/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"
				raw_bam = current_dir+"barcode_sorted/bowtie2/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"
			elif aligner == "hisat":
				both_files = i + " " + current_dir+"barcode_sorted/hisat/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"
				raw_bam = current_dir+"barcode_sorted/hisat/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"
			elif aligner == "bbmap":
				both_files = i + " " + current_dir+"barcode_sorted/bbmap/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"
				raw_bam = current_dir+"barcode_sorted/bbmap/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"
			elif aligner == "bt2-bbmap":
				both_files = i + " " + current_dir+"barcode_sorted/bowtie2/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"
				raw_bam = current_dir+"barcode_sorted/bowtie2/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"

			total_files = len(input_bams)
			cmd = "%s -b %s -o %s --BED %s --labels cleaned raw --regionLabels TSSs CRs intergenics rRNAs snRNAs snoRNAs genes -p %s --perSample --outRawCounts %s" %(plotEnrichment_path, both_files, current_dir+"final_files/EnrichmentPlots/"+i.split("/")[-1].split("_cleaned.bam")[0]+"_enrichementPlot.png",
			beds, processor, current_dir+"final_files/EnrichmentPlots/"+i.split("/")[-1].split("_cleaned.bam")[0]+"_counts.tsv")

			cmd2 = "%s -b %s -o %s --BED %s --labels raw -p %s --perSample --outRawCounts %s" %(plotEnrichment_path, raw_bam, current_dir+"final_files/EnrichmentPlots/"+i.split("/")[-1].split("_cleaned.bam")[0]+"_individual_abundantRNAs_enrichmentPlot.png",
			abundant_beds, processor, current_dir+"final_files/EnrichmentPlots/"+i.split("/")[-1].split("_cleaned.bam")[0]+"_individual_abundantRNAs_counts.tsv")

			cmd3 = "%s -b %s -o %s --BED %s --labels raw -p %s --perSample --outRawCounts %s" %(plotEnrichment_path, raw_bam, current_dir+"final_files/EnrichmentPlots/"+i.split("/")[-1].split("_cleaned.bam")[0]+"_individual_rRNAs_enrichmentPlot.png",
			rRNA_beds, processor, current_dir+"final_files/EnrichmentPlots/"+i.split("/")[-1].split("_cleaned.bam")[0]+"_individual_rRNAs_counts.tsv")

			printStatus("Calculating enrichment over features for %s (%s of %s)." %(i.split("/")[-1], counting_files, total_files))
			subprocess.check_call(cmd, shell=True)
			subprocess.check_call(cmd2, shell=True)
			subprocess.check_call(cmd3, shell=True)

			counting_files += 1

		script_process['EnrichmentPlot'] = 'Finished'
		with open(logfile, 'a') as ff:
			ff.write("** EnrichmentPlot:Finished: %s\n"%time.strftime("%X %x"))


	else:
		input_bams = glob.glob(current_dir+"final_files/bams/*.bam")
		counting_files = 1
		for i in input_bams:
			if aligner == "bt2":
				both_files = i + " " + current_dir+"barcode_sorted/bowtie2/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"
			elif aligner == "hisat":
				both_files = i + " " + current_dir+"barcode_sorted/hisat/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"
			elif aligner == "bbmap":
				both_files = i + " " + current_dir+"barcode_sorted/bbmap/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"
			elif aligner == "bt2-bbmap":
				both_files = i + " " + current_dir+"barcode_sorted/bowtie2/bam_files/" + i.split("/")[-1].split("_cleaned.bam")[0]+".bam"

			total_files = len(input_bams)
			cmd = "%s -b %s -o %s --BED %s --labels cleaned raw -p %s --perSample --outRawCounts %s" %(plotEnrichment_path, both_files, current_dir+"final_files/EnrichmentPlots/"+i.split("/")[-1].split("_cleaned.bam")[0]+"_enrichmentPlot.png", indices[mapping][5], processor,
			current_dir+"final_files/EnrichmentPlots/"+i.split("/")[-1].split("_cleaned.bam")[0]+"_counts.tsv")

			printStatus("Calculating enrichment over features for %s (%s of %s)." %(i.split("/")[-1], counting_files, total_files))
			subprocess.check_call(cmd, shell=True)

			counting_files += 1

		script_process['EnrichmentPlot'] = 'Finished'
		with open(logfile, 'a') as ff:
			ff.write("** EnrichmentPlot:Finished: %s\n"%time.strftime("%X %x"))

def innerdistance(fwd_reads, rvs_reads, raw_fwd, raw_rvs, script_process, indices, mapping, logfile, current_dir, processor):


	folders = [current_dir+"inner_distance", current_dir+"inner_distance/"+fwd_reads.split("/")[-1].split("_R1")[0], current_dir+"inner_distance/"+raw_fwd.split("/")[-1].split("_R1")[0]]

	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)

	input_reads = [fwd_reads, rvs_reads, raw_fwd, raw_rvs]

	for i in input_reads:
		cmd = "zcat %s | head -n 4000000 > %s" %(i, current_dir+"inner_distance/"+i.split("/")[-1].split(".fastq.gz")[0]+"_subsample.fastq")
		subprocess.check_call(cmd, shell=True)

	sub_fastq = glob.glob(current_dir+"inner_distance/*.fastq")
	for i in sub_fastq:
		cmd = "gzip %s" %(i)
		subprocess.check_call(cmd, shell=True)

	sub_gz = glob.glob(current_dir+"inner_distance/*.gz")

	cmd = "%s -x %s -1 %s -2 %s -p %s --end-to-end --fast 2> %s| %s view -Sb - | %s sort - %s"%(bowtie2_path, indices[mapping][0], sub_gz[0], sub_gz[1], processor, current_dir+"inner_distance/"+sub_gz[0].split("/")[-1].split("_R1")[0]+".transcriptome_mapped.bowtie2.log" ,samtools_path, samtools_path, current_dir+"inner_distance/"+sub_gz[0].split("/")[-1].split("_R1")[0]+".transcriptome_mapped")
	subprocess.check_call(cmd, shell=True)

	cmd = "%s -x %s -1 %s -2 %s -p %s --end-to-end --fast 2> %s| %s view -Sb - | %s sort - %s"%(bowtie2_path, indices[mapping][0], sub_gz[2], sub_gz[3], processor, current_dir+"inner_distance/"+sub_gz[2].split("/")[-1].split("_R1")[0]+".transcriptome_mapped.bowtie2.log" ,samtools_path, samtools_path, current_dir+"inner_distance/"+sub_gz[2].split("/")[-1].split("_R1")[0]+".transcriptome_mapped")
	subprocess.check_call(cmd, shell=True)

	sub_mapped = glob.glob(current_dir+"inner_distance/*.bam")

	for i in sub_mapped:
		cmd = "%s -i %s -o %s -l -500 -u 500 -r %s 2> %s 1> %s" %(inner_distance_path, i, current_dir+"inner_distance/"+i.split("/")[-1].split(".transcriptome_mapped.bam")[0], indices[mapping][6], current_dir+"inner_distance/"+i.split("/")[-1].split(".bam")[0]+".inner_distance.log", current_dir+"inner_distance/"+i.split("/")[-1].split(".bam")[0]+".inner_distance.log")
		subprocess.check_call(cmd, shell=True)

		innerdist_files = glob.glob(current_dir+"inner_distance/"+"*"+i.split("/")[-1].split("transcriptome_mapped.bam")[0]+"*")
		for j in innerdist_files:
			cmd = "mv %s %s" %(j, current_dir+"inner_distance/"+i.split("/")[-1].split(".transcriptome_mapped.bam")[0]+"/")
			subprocess.check_call(cmd, shell=True)

	for i in sub_gz:
		cmd = "mv %s %s" %(i, current_dir+"inner_distance/"+i.split("/")[-1].split("_R")[0])
		subprocess.check_call(cmd, shell=True)


	script_process['InnerDistance'] = 'Finished'

	with open(logfile, 'a') as ff:
		ff.write("** InnerDistance:Finished: %s\n"%time.strftime("%X %x"))

if __name__ == "__main__":
	start_time = time.time()
	printStatus("Initializing the script.")
	main()
	end_time = time.time() - start_time
	printStatus("lollipy successfully completed. Enjoy your results.\nElapsed time: %s seconds." %round(end_time))
