#!/home/semplicio/virtualenv/bin/python -u

from sys import argv

script, infile, outfile = argv

aa = open(infile, 'r')

line1 = aa.next()
line2 = aa.next()

chr1 = line1.split("\t")[0]
start1 = line1.split("\t")[1]
count1 = int(line1.split("\t")[-1])
chr2 = line2.split("\t")[0]
start2 = line2.split("\t")[1]
count2 = int(line2.split("\t")[-1])

final_count = count1

output_set = set()

if chr1 == chr2:
	if start1 == start2:
		final_count += count2
	elif start1 != start2:
		output_set.add("%s\t%s\t%s\n" %(chr1, start1, final_count))
		final_count = count2 
elif chr1 != chr2:
	output_set.add("%s\t%s\t%s\n" %(chr1, start1, final_count))
	final_count = count2 

line1 = line2

count = 0
chunk = 100000
multiplier = 1

while True:
	try:
		chr1 = line1.split("\t")[0]
		start1 = line1.split("\t")[1]
		count1 = int(line1.split("\t")[-1])
		line2 = aa.next()
		chr2 = line2.split("\t")[0]
		start2 = line2.split("\t")[1]
		count2 = int(line2.split("\t")[-1])
		
		count += 1
		if chr1 == chr2:
			if start1 == start2:
				final_count += count2
			elif start1 != start2:
				output_set.add("%s\t%s\t%s\n" %(chr1, start1, final_count))
				final_count = count2 
		elif chr1 != chr2:
			output_set.add("%s\t%s\t%s\n" %(chr1, start1, final_count))
			final_count = count2 
		line1 = line2


		if count == multiplier * chunk:
			outfile_open = open(outfile, 'a')
			for i in output_set:
				outfile_open.write(i)
			outfile_open.close()

			multiplier += 1
			output_set = set()

	except StopIteration:
		output_set.add("%s\t%s\t%s\n" %(chr1, start1, final_count))
		outfile_open = open(outfile, 'a')
		for i in output_set:
			outfile_open.write(i)
		outfile_open.close()

		break


