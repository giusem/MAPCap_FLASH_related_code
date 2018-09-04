


# remove duplicate entries from SAF file, based on gene name (last column) and start and end site (column 3 and 4)
# File does not have to be sorted beforehand


import subprocess

print "name of the input file:"
inputfile = raw_input("> ")



print "name of output file:"
output = raw_input("> ")

in_file = open(inputfile, 'r')

print "sorting the input file\n"
cmd = "sort -k 2,2 -k 3,3n %s > %s" %(inputfile, ".tmp_infile_sorted.SAF")
subprocess.check_call(cmd, shell=True)

in_file = open(".tmp_infile_sorted.SAF", 'r')

print "removing duplicate entries from SAF file\n"

line1 = in_file.next()
line2 = in_file.next()

#gene1 = line1.split("\t")[-1]
#gene2 = line2.split("\t")[-1]

chr1 = line1.split("\t")[1]
chr2 = line2.split("\t")[1]

start1 = line1.split("\t")[2]
start2 = line2.split("\t")[2]

#end1 = line1.split("\t")[3]   end position not needed because all transcripts have same length
#end2 = line2.split("\t")[3]

output_set = set()

if chr1 == chr2:
	if start1 != start2:
		output_set.add(line1)

line1 = line2

while True:
	try:
		line2 = in_file.next()
		
		chr1 = line1.split("\t")[1]
		chr2 = line2.split("\t")[1]

		start1 = line1.split("\t")[2]
		start2 = line2.split("\t")[2]

		if chr1 == chr2:
			if start1 != start2:
				output_set.add(line1)

		else:
			output_set.add(line1)
			
		line1 = line2
		

	except StopIteration:
		print "writing to file...\n"
		for i in output_set:
			with open(output, 'a') as ff:
				ff.write(i)
		break
in_file.close()

cmd = "mv %s %s" %(output, ".tmp_output.SAF")
subprocess.check_call(cmd, shell=True)

cmd = "sort -k 2,2 -k 3,3n %s > %s" %(".tmp_output.SAF", output)
subprocess.check_call(cmd, shell=True)

cmd = "rm .tmp_infile_sorted.SAF .tmp_output.SAF"
subprocess.check_call(cmd, shell=True)

print "all done, have a nice day"