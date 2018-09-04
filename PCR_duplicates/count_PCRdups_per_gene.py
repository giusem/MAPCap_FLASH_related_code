#!/home/semplicio/virtualenv/bin/python -u

from sys import argv

script, inputfile, outputfile = argv

closestBed = open(inputfile, 'r')

line1 = closestBed.next().strip()
count1=line1.split("\t")[3]
gene1 = line1.split("\t")[7]
line2 = closestBed.next().strip()
count2=line2.split("\t")[3]
gene2 = line2.split("\t")[7]

genecount = int(count1)

output_list = []

if gene1 == gene2:
	genecount = genecount + int(count2)
else:
	output_list.append([gene1, genecount])
	genecount = int(count1)

line1 = line2

while True:
	try:
		count1=line1.split("\t")[3]
		gene1 = line1.split("\t")[7]
		line2 = closestBed.next().strip()
		count2=line2.split("\t")[3]
		gene2 = line2.split("\t")[7]

		if gene1 == gene2:
			genecount = genecount + int(count2)
		else:
			output_list.append([gene1, genecount])
			genecount = int(count1)

		line1 = line2

	except StopIteration:
		break


fileopen = open(outputfile,'a')
for i in output_list:
	fileopen.write("%s\t%s\n"%(i[0], i[1]))
fileopen.close()

closestBed.close()
