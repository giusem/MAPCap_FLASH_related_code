#!/home/semplicio/virtualenv/bin/python -u

print "this little tool extracts start- and end-sites of each exon and intron form a bed12 format file."
print "now also intergenic regions"

print "input file please..."
input_bed = raw_input(">")

print "output prefix please..."
output_prefix = raw_input(">")

bed = open(input_bed, 'r')

exon_dict = {}
intron_dict = {}
utr5_dict = {}
utr3_dict = {}

count = 0
multiplier = 1
chunk = 100000

#for exons and introns
print "working on introns, exons and UTRs ..."
while True:
	try:
		line = bed.next().strip()

		count += 1

		chrom = line.split("\t")[0]
		start = line.split("\t")[1]
		end = line.split("\t")[2]
		utr5_end = line.split("\t")[6]
		utr3_start = line.split("\t")[7]

		#write the UTRs out
		try:
			utr5_dict[output_prefix+"_5utr.bed"].append(["%s\t%s\t%s\n" %(chrom, start, utr5_end)])
			utr3_dict[output_prefix+"_3utr.bed"].append(["%s\t%s\t%s\n" %(chrom, utr3_start, end)])
		except KeyError:
			utr5_dict[output_prefix+"_5utr.bed"] = [["%s\t%s\t%s\n" %(chrom, start, utr5_end)]]
			utr3_dict[output_prefix+"_3utr.bed"] = [["%s\t%s\t%s\n" %(chrom, utr3_start, end)]]


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
			


		end_position = [int(start_position[i]) + int(length[i]) for i in xrange(len(start_position))]

		exon_start = [int(start) + int(start_position[i]) for i in xrange(len(start_position))]
		exon_end = [int(start) + int(end_position[i]) for i in xrange(len(end_position))]

		intron_start = [int(exon_end[i]) for i in xrange(len(exon_end) - 1)]
		intron_end = [int(exon_start[i+1]) for i in xrange(len(exon_start) - 1)]

		for i in xrange(len(exon_start)):
			try:
				exon_dict[output_prefix+"_exons.bed"].append(["%s\t%s\t%s\n" %(chrom, exon_start[i], exon_end[i])])
			except KeyError:
				exon_dict[output_prefix+"_exons.bed"] = [["%s\t%s\t%s\n" %(chrom, exon_start[i], exon_end[i])]]

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
bed = open(input_bed, 'r')

intergenic_dict = {}

count = 0
multiplier = 1
chunk = 100000

print "working on the intergenic regions..."
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

bed.close()