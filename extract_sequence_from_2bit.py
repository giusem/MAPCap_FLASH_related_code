


# extract sequence of interest from fasta files
# takes regions from bed file and finds the sequence in the fasta file
# if desired counts the number of unique sequences


import subprocess
import py2bit # pip install git+https://github.com/dpryan79/py2bit
import os




print "name of the 2bit file:"
tb_input = raw_input("> ")
tb_input = os.path.abspath(tb_input)
tb = py2bit.open(tb_input)

print "name of bed file:"
bed_input = raw_input("> ")
bed_input = os.path.abspath(bed_input)
bed = open(bed_input, 'r')

print "name of the output file:"
output_name =raw_input("> ")
output_name = output_name+".fa"

sequence_list = []
count = 0

while True:
	try:
		line = bed.next()
		chrom = line.split("\t")[0]
		start = int(line.split("\t")[1])
		end = int(line.split("\t")[2])

		count += 1
		sequence_list.append([chrom, start, end, tb.sequence(chrom, start, end)])
	except StopIteration:
		break

bed.close()
tb.close()

outfile = open(output_name, 'a')
for i in sequence_list:
	outfile.write(">%s:%s:%s\n%s\n"%(i[0],i[1],i[2],i[3]))

outfile.close()


#with open(output_name, 'a') as ff:
#	ff.write("2bit file: %s\nBed file: %s\nTotal number of sequences: %s\n" %(tb_input, bed_input, count))

#nuc_list = ['A', 'T', 'G', 'C']

#for i in nuc_list:
#	seq_count = sequence_list.count(i)
#	with open(output_name, 'a') as ff:
#		ff.write("%s\t%s\n" %(i, seq_count))


