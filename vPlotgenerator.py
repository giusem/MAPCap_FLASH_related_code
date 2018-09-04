#!/home/semplicio/virtualenv/bin/python -u

# extract genomic distance between center of read and desired locus
# useful for V plots

import os
import argparse
import subprocess
import sys
import glob


def getLociSams(inputfile, outfile, bamfile, paired_end):


	if paired_end:
		print "running in paired_end mode"

		genome = "/data/repository/organisms/dm6_ensembl/genome_fasta/genome.chrom.sizes"
		
		temp1 = bamfile.split(".bam")[0]+"_mappedPairsOnly.bam"
		temp2 = bamfile.split(".bam")[0]+"_mappedPairsOnly_NameSorted"
		temp3 = bamfile.split(".bam")[0]+".bed"
		temp4 = bamfile.split(".bam")[0]+"_se.bed"
		temp5 = bamfile.split(".bam")[0]+"_se.bam"
		temp6 = bamfile.split(".bam")[0]+"_se_sorted"

		
		if os.path.exists(temp6+".bam"):
			print "\n\tUsing the bam file from a previous run. \n\tIf you don't want that, then (re)move the bam file."
			if not os.path.exists (temp6+".bam.bai"):
				cmd = "samtools index %s" %temp6+".bam"
				subprocess.check_call(cmd, shell=True)

		else:
			#print "\tKeeping only mapped reads of proper pairs"
			#cmd = "samtools view -b -F 4 -f 0x2 %s > %s" %(bamfile, temp1) #keep only mapped reads and proper pairs
			#subprocess.check_call(cmd, shell=True)
			
			print "\tSorting new bam file by name"
			cmd = "samtools sort -n %s %s" %(bamfile, temp2) #sort by name
			subprocess.check_call(cmd, shell=True)
			
			print "\tbamToBed"
			cmd = "bamToBed -i %s -bedpe > %s" %(temp2+".bam", temp3) 
			subprocess.check_call(cmd, shell=True)
			
			print "\tMaking single-end bed file"
			cmd = "awk '{print$1,$2,$6,$7,$8,\".\"}' %s | tr [:blank:] \\\\t > %s" %(temp3, temp4)
			subprocess.check_call(cmd, shell=True)

			print "\tCreating single-end bam file"
			cmd = "bedToBam -i %s -g %s > %s" %(temp4, genome, temp5)
			subprocess.check_call(cmd, shell=True)

			print "\tSorting and indexing single-end bam file"
			cmd = "samtools sort %s %s" %(temp5, temp6)
			subprocess.check_call(cmd, shell=True)
			cmd = "samtools index %s" %temp6+".bam"
			subprocess.check_call(cmd, shell=True)

			garbagelist = [temp1, temp2+".bam", temp3, temp4, temp5]
			for i in garbagelist:
				os.remove(i)

		bamfile = temp6+".bam"

	else:
		print "running in single end mode"

	if not os.path.exists(bamfile+".bai"):
		print "vPlotgenerator only works on sorted and indexed bam files."
		sys.exit(1)

	if not os.path.exists(outfile+"_lociSams"):
		os.makedirs(outfile+"_lociSams")


	infile = open(inputfile, 'r')
	while True:
		try:

			line = infile.next().strip()
			chrom = line.split("\t")[0]
			start = line.split("\t")[1]
			end = line.split("\t")[2]

			cmd = "samtools view %s %s:%s-%s >> %s" %(bamfile, chrom, start, end, outfile+"_lociSams/"+chrom+"x"+start+"x"+end+".sam")
			subprocess.check_call(cmd, shell=True)

		except StopIteration:
		
			break

	output_folder = os.path.abspath(outfile)+"_lociSams/"
	return output_folder


def getGenomicDistance(sam_folder, dist_min, dist_max, read_min, read_max, outfile, paired_end):

	print "\n\tGenerating matrix for V-plot."
	samfiles = glob.glob(sam_folder+"/*")
	output_list = []

	for i in samfiles:
		name = i.split("/")[-1]
		locus_chr = name.split("x")[0]
		locus_start = int(name.split("x")[1])
		locus_end = int(name.split("x")[2].split(".sam")[0])

		locus_center = locus_start + (locus_end - locus_start)/2

		sam_open = open(i, 'r')

		while True:
			try:
				line = sam_open.next().strip()
				
				if paired_end:
					read_start = int(line.split("\t")[3])
					read_length = int(line.split("\t")[5].split("M")[0])
				else:
					read_start = int(line.split("\t")[3])
					read_length = int(len(line.split("\t")[9]))

				read_end = int(read_start+read_length)
				read_center = read_start+(read_end - read_start)/2
				
				delta = read_center - locus_center

				output_list.append([read_length,delta])

			except StopIteration:

				break


	fileopen = open(output+"_scatter.tsv",'a')
	for i in output_list:
		fileopen.write("%s\t%s\n"%(i[0], i[1]))
	fileopen.close()


	# this converts the distances into frequency of distances for each list in the dictionary
	#for key, value in delta_dict.iteritems():
	#    for i in range(dist_min, dist_max, 1):
	#        try:
	#        	frequency_dict[key].append(value.count(i))
	#        except KeyError:
	#        	frequency_dict[key] = [value.count(i)]

	# create a list of 0s with the same length to fill the empty spots
	#blank_list = []
	#for i in range(dist_min, dist_max,1):
	#	blank_list.append(0)

	# this prints out in column 1 the read length followed by the frequencies, sorted by read length (y axis) and distance (x axis)
	#if read_max == 0:
	#	read_max = max(frequency_dict.viewkeys())
	#	print "\tLongest fragment: %s" %read_max

	#for fragment_length in range(read_min, read_max,1):
	#	if frequency_dict.has_key(fragment_length):
	#		with open('%s' %outfile+".tsv", 'a+') as ff:
	#			ff.write('%s\t' %fragment_length)
	#			for i in frequency_dict[fragment_length]:
	#				ff.write('%s\t' %i)
	#			ff.write('\n')
	#	else:
	#		with open('%s' %outfile+".tsv", 'a+') as ff:
	#			ff.write('%s\t' %fragment_length)
	#			for i in blank_list:
	#				ff.write('%s\t' %i)
	#			ff.write('\n')



def main():
	parser = argparse.ArgumentParser(description="Creates the count matrix for a V plot (Hennikoff et al, 2011, PNAS). Calculates frquencies of distances between fragment center and locus center for each fragment length.")
	parser.add_argument('-R', '--region', metavar='path/to/regions.bed', help='The path to the bed file containing the reference regions.', required=True)
	parser.add_argument('-o', '--out', metavar='output file name', help='The prefix of the output.', required=True)
	parser.add_argument('-b', '--inBAM', metavar='path/to/input.bam', help='The path to the BAM file to extract reads from.', required=True)
	parser.add_argument('-minDist', '--dist_minimum', type=int, help='The lower distance limit of the read center and the locus center.', required=True)
	parser.add_argument('-maxDist', '--dist_maximum', type=int, help='The upper distance limit of the read center and the locus center.', required=True)
	parser.add_argument('-minRead', '--read_minimum', type=int, help='The minimum fragment length to consider. (Default = 0)', default=0)
	parser.add_argument('-maxRead', '--read_maximum', type=int, help='The maximum fragment length to consider. (Default = maximum fragment length)', default=0)	
	parser.add_argument('-pe', '--paired_end', action='store_true', help='Use this option in case of paired-end data.')
	args = parser.parse_args()

	current_dir = os.getcwd()+"/"
	bedfile = os.path.abspath(args.region)
	outfile = current_dir+args.out
	bamfile = os.path.abspath(args.inBAM)

	print "\nBED file: %s" %bedfile
	print "BAM file: %s" %bamfile
	print "Output prefix: %s" %outfile
	print "Min and max distance from target: %s\t%s" %(args.dist_minimum, args.dist_maximum)
	if args.read_maximum == 0:
		print "Min and max fragment length: %s\tmaximum fragment length\n" %(args.read_minimum)
	else:
		print "Min and max fragment length: %s\t%s\n" %(args.read_minimum, args.read_maximum)

	if os.path.exists(outfile+".tsv"):
		print "vPlotgenerator.py was already run with these settings. (Re)move the files or change the output file prefix."
		sys.exit(1)

	sam_output_folder = getLociSams(bedfile, outfile, bamfile, args.paired_end)

	getGenomicDistance(sam_output_folder, args.dist_minimum, args.dist_maximum, args.read_minimum, args.read_maximum, outfile, args.paired_end)

	print "\nDone. Have a lovely day."


if __name__ == "__main__":
	main()