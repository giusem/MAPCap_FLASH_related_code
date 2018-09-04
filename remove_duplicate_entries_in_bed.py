# this script is to remove duplicate entries from a bed or bedGraph
# duplicate entries are decided by the 1st and 2nd column (chromosome and start position)
# file must be sorted by chr, start e.g. sort -k1,1 -k2,2n
# if bedgraph than sort the score in reverse as well if you want to keep the highest score e.g. -k4,4nr

print "filename"
infile = raw_input("> ")

outfile = infile.split(".")[0]+"_unique."+infile.split(".")[-1]

output_set = set()
count = 0

infile = open(infile, 'r')

line1 = infile.next()
line2 = infile.next()

chr1 = line1.split("\t")[0]
start1 = line1.split("\t")[1]
chr2 = line2.split("\t")[0]
start2 = line2.split("\t")[1]

if chr1 == chr2:
	if start1 == start2:
		count += 1
	else:
		output_set.add(line1)
else:
	output_set.add(line1)

line1 = line2

while True:
	try:
		chr1 = line1.split("\t")[0]
		start1 = line1.split("\t")[1]

		line2 = infile.next()

		chr2 = line2.split("\t")[0]
		start2 = line2.split("\t")[1]

		if chr1 == chr2:
			if start1 == start2:
				count += 1
			else:
				output_set.add(line1)
		else:
			output_set.add(line1)

		line1 = line2
	except StopIteration:
		with open(outfile, 'a') as ff:
			for i in output_set:
				ff.write(i)

		output_set = set()

		break

print "number of duplicate entries removed: %s" %count

infile.close()
