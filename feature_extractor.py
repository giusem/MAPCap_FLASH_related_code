#!/home/semplicio/virtualenv/bin/python -u

#the start is from a GTF file to extract the transcript name of the longest isoform
#then this list is taken to extract all targets form bed12 file

#GTF file has to be sorted by chromosome, start in ascending order and end in descending order and should contain only transcsript annotations
#use:
#sort -k1,1 -k4,4n -k5,5nr dm6_simple_genes.gtf | awk '{if($3 == "transcript") print $0}' > dm6_simple_only_transcripts.gtf

import subprocess 
import os


print "this little tool extracts start- and end-sites of each exon and intron form a GTF format file."
print "now also intergenic regions"
print "and UTR\n"

print "note that exons still contain UTRs"

print "GTF file please..."
input_gtf = raw_input(">")

print "BED file please..."
input_bed = raw_input(">")

print "output prefix please..."
output_prefix = raw_input(">")

print "stupid mouse?"
mouse = raw_input(">")

print "\ngetting transcripts from GTF file"
if mouse:
	cmd = "sort -k1,1 -k4,4n -k5,5n %s > temp.gtf" %input_gtf
	subprocess.check_call(cmd, shell=True)
else:
	cmd = "sort -k1,1 -k4,4n -k5,5n %s | awk '{if($3 == \"transcript\") print $0}' > temp.gtf" %input_gtf
	subprocess.check_call(cmd, shell=True)

print "extracting transcript names of the longest isoform"
gtf = open("temp.gtf", 'r')

if mouse:
	transcript_set = set()

	line1 = gtf.next()
	line2 = gtf.next()

	fbgn1 = line1.split(";")[0].split("\"")[1]
	fbgn2 = line2.split(";")[0].split("\"")[1]

	fbtr1 = line1.split(";")[2].split("\"")[1]
	fbtr2 = line2.split(";")[2].split("\"")[1]

	start1 = int(line1.split("\t")[3])
	start2 = int(line2.split("\t")[3])

	end1 = int(line1.split("\t")[4])
	end2 = int(line2.split("\t")[4])

	ENSMUST1 = line1.split(";")[1].split("\"")[1]
	ENSMUST2 = line2.split(";")[1].split("\"")[1]

	if fbgn1 == fbgn2:
	    if ENSMUST1 == ENSMUST2:
	        line1 = line1
	    else:
	        transcript_set.add("%s\t%s" % (fbgn1, ENSMUST1))
	        line1 = line2
	else:
	    transcript_set.add("%s\t%s" % (fbgn1, ENSMUST1))
	    line1 = line2

	while True:
	    try:
	        line2 = gtf.next()

	        fbgn1 = line1.split(";")[0].split("\"")[1]
	        fbgn2 = line2.split(";")[0].split("\"")[1]

	        fbtr1 = line1.split(";")[2].split("\"")[1]
	        fbtr2 = line2.split(";")[2].split("\"")[1]

	        start1 = int(line1.split("\t")[3])
	        start2 = int(line2.split("\t")[3])

	        end1 = int(line1.split("\t")[4])
	        end2 = int(line2.split("\t")[4])

	        ENSMUST1 = line1.split(";")[1].split("\"")[1]
	        ENSMUST2 = line2.split(";")[1].split("\"")[1]

	        if fbgn1 == fbgn2:
	            if ENSMUST1 == ENSMUST2:
	                line1 = line1
	            else:
	                transcript_set.add("%s\t%s" % (fbgn1, ENSMUST1))
	                line1 = line2
	        else:
	            transcript_set.add("%s\t%s" % (fbgn1, ENSMUST1))
	            line1 = line2
	    except StopIteration:
	        ENSMUST1 = line1.split(";")[1].split("\"")[1]
	        transcript_set.add("%s\t%s" % (fbgn1, ENSMUST1))

	        with open("temp.txt", 'a') as ff:
	            for i in transcript_set:
	                ff.write("%s\n" % i)
	        transcript_set = set()

	        break

	gtf.close()

	cmd = "sort -k1,1 temp.txt > temp_sorted.txt"
	subprocess.check_call(cmd, shell=True)

	cmd = "awk '{print$2}' temp_sorted.txt > col2"
	subprocess.check_call(cmd, shell=True)

	cmd = "while read f; do grep $f %s >> temp.bed; done < col2" %input_bed
	subprocess.check_call(cmd, shell=True)

	cmd = "awk '{print $3-$2}' temp.bed > temp_length"
	subprocess.check_call(cmd, shell=True)

	cmd = "paste temp_sorted.txt temp_length > temp_final.txt"
	subprocess.check_call(cmd, shell=True)

	os.remove("temp.txt")
	os.remove("temp_sorted.txt")
	os.remove("col2")
	os.remove("temp.bed")
	os.remove("temp_length")

	lengthfile = open("temp_final.txt", 'r')


	transcript_set = set()

	line1 = lengthfile.next()
	line2 = lengthfile.next()

	fbgn1 = line1.split("\t")[0]
	fbgn2 = line1.split("\t")[0]

	fbtr1 = line1.split("\t")[1]
	fbtr2 = line2.split("\t")[1]

	start1 = int(line1.split("\t")[2])
	start2 = int(line2.split("\t")[2])


	if fbgn1 == fbgn2:
	    if start1 >= start2:
	        line1 = line1
	    else:
	        line1 = line2
	else:
	    transcript_set.add(fbtr1)
	    line1 = line2

	while True:
	    try:
	        line2 = lengthfile.next()

	        fbgn1 = line1.split("\t")[0]
	        fbgn2 = line2.split("\t")[0]

	        fbtr1 = line1.split("\t")[1]
	        fbtr2 = line2.split("\t")[1]

	        start1 = int(line1.split("\t")[2])
	        start2 = int(line2.split("\t")[2])


	        if fbgn1 == fbgn2:
	            if start1 >= start2:
	                line1 = line1
	            else:
	                line1 = line2
	        else:
	            transcript_set.add(fbtr1)
	            line1 = line2

	    except StopIteration:
	        transcript_set.add(fbtr1)

	        with open("temp.txt", 'a') as ff:
	            for i in transcript_set:
	                ff.write("%s\n" % i)
	        transcript_set = set()

	        break

	os.remove("temp_final.txt")

else:
	transcript_set = set()

	line1 = gtf.next()
	line2 = gtf.next()

	fbgn1 = line1.split(";")[0].split("\"")[1]
	fbgn2 = line2.split(";")[0].split("\"")[1]

	fbtr1 = line1.split(";")[2].split("\"")[1]
	fbtr2 = line2.split(";")[2].split("\"")[1]

	start1 = int(line1.split("\t")[3]) 
	start2 = int(line2.split("\t")[3]) 

	end1 = int(line1.split("\t")[4])
	end2 = int(line2.split("\t")[4])

	
	if fbgn1 == fbgn2:
		#print "same gene"
		if end1 - start1 >= end2 - start2:
			line1 = line1
		#	print "1 bigger than 2"
		else:
			line1 = line2 
		#	print "2 bigger than 1"
	else:
		#print "not same gene"
		transcript_set.add(fbtr1)
		line1 = line2


	while True:
		try:
			line2 = gtf.next()

			fbgn1 = line1.split(";")[0].split("\"")[1]
			fbgn2 = line2.split(";")[0].split("\"")[1]

			fbtr1 = line1.split(";")[2].split("\"")[1]
			fbtr2 = line2.split(";")[2].split("\"")[1]

			start1 = int(line1.split("\t")[3]) 
			start2 = int(line2.split("\t")[3]) 

			end1 = int(line1.split("\t")[4])
			end2 = int(line2.split("\t")[4])

			#print line1, line2
		
			if fbgn1 == fbgn2:
				#print "same gene"
				if end1 - start1 >= end2 - start2:
					line1 = line1
				#	print "1 bigger than 2"
				else:
					line1 = line2 
				#	print "2 bigger than 1"
			else:
				#print "not same gene"
				transcript_set.add(fbtr1)
				line1 = line2

		except StopIteration:
			fbtr1 = line1.split(";")[2].split("\"")[1]
			transcript_set.add(fbtr1)

			with open("temp.txt", 'a') as ff:
				for i in transcript_set:
					ff.write("%s\n" %i)
			transcript_set = set()

			break

	gtf.close()

print "getting the longest isoform from Bed12"
cmd = "while read f; do grep $f %s >> longest_isoform_temp.bed ; done < temp.txt" %(input_bed)
subprocess.check_call(cmd, shell=True)

cmd = "sort -k1,1 -k2,2n -k3,3n longest_isoform_temp.bed > longest_isoform_temp_sorted.bed"
subprocess.check_call(cmd, shell=True)

bed = open("longest_isoform_temp_sorted.bed", 'r')

exon_dict = {}
intron_dict = {}
utr5_dict = {}
utr3_dict = {}

count = 0
multiplier = 1
chunk = 100000

#for exons and introns
print "extract introns, exons and UTRs ..."
while True:
	try:
		line = bed.next().strip()
		#count += 1

		chrom = line.split("\t")[0]
		start = line.split("\t")[1]
		end = line.split("\t")[2]
		utr5_end = line.split("\t")[6]
		utr3_start = line.split("\t")[7]
		direction = line.split("\t")[5]

		#write the UTRs out
				
		if direction == "+":
			if start != utr5_end and end != utr3_start:
				try:
					utr5_dict[output_prefix+"_5utr.bed"].append(["%s\t%s\t%s\n" %(chrom, start, utr5_end)])
				except KeyError:
					utr5_dict[output_prefix+"_5utr.bed"] = [["%s\t%s\t%s\n" %(chrom, start, utr5_end)]]
			if start != utr5_end and end != utr3_start:
				try:
					utr3_dict[output_prefix+"_3utr.bed"].append(["%s\t%s\t%s\n" %(chrom, utr3_start, end)])
				except KeyError:
					utr3_dict[output_prefix+"_3utr.bed"] = [["%s\t%s\t%s\n" %(chrom, utr3_start, end)]]


		elif direction == "-":
			if start != utr5_end and end != utr3_start:
				try:
					utr3_dict[output_prefix+"_3utr.bed"].append(["%s\t%s\t%s\n" %(chrom, start, utr5_end)])
				except KeyError:
					utr3_dict[output_prefix+"_3utr.bed"] = [["%s\t%s\t%s\n" %(chrom, start, utr5_end)]]
			if start != utr5_end and end != utr3_start:
				try:
					utr5_dict[output_prefix+"_5utr.bed"].append(["%s\t%s\t%s\n" %(chrom, utr3_start, end)])
				except KeyError:
					utr5_dict[output_prefix+"_5utr.bed"] = [["%s\t%s\t%s\n" %(chrom, utr3_start, end)]]

		#because these stupid bed12 files do not stick to one type of annotation we need to check each line whether start position
		#and length are two separate columns or not
		#it might have been an artifact from copying one line from bed file. 


		try:
			start_position = line.split("\t")[11].split(",")[0:-1]
			length = line.split("\t")[10].split(",")[0:-1]

		except IndexError:
			coordinates = line.split("\t")[10].split(",0,")
			temp = "0," + coordinates[1]
			start_position = temp.split(",")[0:-1]
			length = coordinates[0].split(",")
						

		end_position = [int(start_position[i]) + int(length[i]) for i in range(len(start_position))]

		if len(start_position) == 1:
			utr5_length = int(utr5_end)-int(start)
			exon_start = [int(utr5_end)]
			exon_end = [utr3_start]
			try:
				exon_dict[output_prefix+"_exons.bed"].append(["%s\t%s\t%s\n" %(chrom, exon_start[0], exon_end[0])])
			except KeyError:
				exon_dict[output_prefix+"_exons.bed"] = [["%s\t%s\t%s\n" %(chrom, exon_start[0], exon_end[0])]]
		else:
			for i in range(len(start_position)):
				if i == 0:
					utr5_length = int(utr5_end)-int(start)
					exon_start = [int(utr5_end)]
					exon_end = [int(exon_start[i]) + int(length[i]) - utr5_length]
					try:
						exon_dict[output_prefix+"_exons.bed"].append(["%s\t%s\t%s\n" %(chrom, exon_start[i], exon_end[i])])
					except KeyError:
						exon_dict[output_prefix+"_exons.bed"] = [["%s\t%s\t%s\n" %(chrom, exon_start[i], exon_end[i])]]

				elif i == int(len(start_position) - 1):
					exon_start = [int(start) + int(start_position[-1])]
					exon_end = [int(utr3_start)]
					try:
						exon_dict[output_prefix+"_exons.bed"].append(["%s\t%s\t%s\n" %(chrom, exon_start[0], exon_end[0])])
					except KeyError:
						exon_dict[output_prefix+"_exons.bed"] = [["%s\t%s\t%s\n" %(chrom, exon_start[0], exon_end[0])]]

				else:
					exon_start = [int(start) + int(start_position[i])]
					exon_end = [int(exon_start[0]) + int(length[i])]
					try:
						exon_dict[output_prefix+"_exons.bed"].append(["%s\t%s\t%s\n" %(chrom, exon_start[0], exon_end[0])])
					except KeyError:
						exon_dict[output_prefix+"_exons.bed"] = [["%s\t%s\t%s\n" %(chrom, exon_start[0], exon_end[0])]]
		
		if len(start_position) > 1:
			exon_start = [int(start) + int(start_position[i]) for i in xrange(len(start_position))]
			exon_end = [int(exon_start[i]) + int(length[i]) for i in xrange(len(exon_start))]
			intron_start = [int(exon_end[i]) for i in xrange(len(exon_end) - 1)]
			intron_end = [int(exon_start[i+1]) for i in xrange(len(exon_start) - 1)]

			for i in xrange(len(intron_start)):
				try:
					intron_dict[output_prefix+"_introns.bed"].append(["%s\t%s\t%s\n" %(chrom, intron_start[i], intron_end[i])])
				except KeyError:
					intron_dict[output_prefix+"_introns.bed"] = [["%s\t%s\t%s\n" %(chrom, intron_start[i], intron_end[i])]]


		if count == multiplier*chunk:
			print"."
			for key, value in exon_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			for key, value in intron_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			for key, value in utr5_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			for key, value in utr3_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			multiplier += 1
			exon_dict = {}
			intron_dict = {}

	except StopIteration:

		for key, value in exon_dict.iteritems():
			with open('%s' %key, 'a+') as ff:
				for i in range(len(value)):
					ff.write('%s' %value[i][0])

		for key, value in intron_dict.iteritems():
			with open('%s' %key, 'a+') as ff:
				for i in range(len(value)):
					ff.write('%s' %value[i][0])

		for key, value in utr5_dict.iteritems():
			with open('%s' %key, 'a+') as ff:
				for i in range(len(value)):
					ff.write('%s' %value[i][0])

		for key, value in utr3_dict.iteritems():
			with open('%s' %key, 'a+') as ff:
				for i in range(len(value)):
					ff.write('%s' %value[i][0])

		multiplier += 1
		exon_dict = {}
		intron_dict = {}

		break

bed.close()

#for intergenic regions
bed = open("longest_isoform_temp_sorted.bed", 'r')

intergenic_dict = {}

count = 0
multiplier = 1
chunk = 100000

print "extract intergenic regions..."
line1 = bed.next()
line2 = bed.next()

count += 1

chrom1 = line1.split("\t")[0]
chrom2 = line2.split("\t")[0]
start1 = line1.split("\t")[1]
start2 = line2.split("\t")[1]
end1 = line1.split("\t")[2]
end2 = line2.split("\t")[2]

if chrom1 == chrom2:
	if start2 > end1:
		try:
			intergenic_dict[output_prefix+"_intergenics.bed"].append(["%s\t%s\t%s\n" %(chrom1, end1, start2)])
		except KeyError:
			intergenic_dict[output_prefix+"_intergenics.bed"] = [["%s\t%s\t%s\n" %(chrom1, end1, start2)]]


line1 = line2

newcount = 0

while True:
	try:
		line2 = bed.next()
	
		count += 1

		chrom1 = line1.split("\t")[0]
		chrom2 = line2.split("\t")[0]
		start1 = line1.split("\t")[1]
		start2 = line2.split("\t")[1]
		end1 = line1.split("\t")[2]
		end2 = line2.split("\t")[2]
		
		if chrom1 == chrom2:
			if start2 > end1:
				try:
					intergenic_dict[output_prefix+"_intergenics.bed"].append(["%s\t%s\t%s\n" %(chrom1, end1, start2)])
				except KeyError:
					intergenic_dict[output_prefix+"_intergenics.bed"] = [["%s\t%s\t%s\n" %(chrom1, end1, start2)]]

		line1 = line2

	

		if count == multiplier*chunk:
			print"."
			for key, value in intergenic_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			multiplier += 1
			intergenic_dict = {}


	except StopIteration:

		for key, value in intergenic_dict.iteritems():
			with open('%s' %key, 'a+') as ff:
				for i in range(len(value)):
					ff.write('%s' %value[i][0])

		multiplier += 1
		intergenic_dict = {}

		break

os.remove("temp.txt")
os.remove("longest_isoform_temp.bed")
os.remove("longest_isoform_temp_sorted.bed")
os.remove("temp.gtf")
bed.close()

print "Done. Enjoy your day."