

# 1. takes bed file and 2bit file as input
# 2. based on bed file locations extract the underlying sequence and write into a fasta format
# 3. fasta output with name from bed file (column 4) if not provided use location as name



import subprocess
import py2bit # pip install git+https://github.com/dpryan79/py2bit
import os


tb_input = "/data/repository/organisms/dm6_ensembl/genome_fasta/genome.2bit"
tb = py2bit.open(tb_input)

print "2bit file: %s" %tb_input

print "name of bed file:"
bed_input = raw_input("> ")
bed_input = os.path.abspath(bed_input)
bed = open(bed_input, 'r')

print "name of the output file:"
output_name =raw_input("> ")
output_name = output_name+".fa"


#check whether there is a column 4 in bed file for names
line = bed.next().strip()
colnum = len(line.split("\t"))

if colnum >= 4:
	namecol = True
else:
	namecol = False
bed.close()

bed = open(bed_input, 'r')
sequence_list = []

while True:
	try:
		line = bed.next().strip()
		chrom = line.split("\t")[0]
		start = int(line.split("\t")[1])
		end = int(line.split("\t")[2])
		if namecol:
			genename = line.split("\t")[3]
		else:	
			genename = "".join(line.split("\t")[0:3])


		sequence_list.append(">"+genename+"\n"+tb.sequence(chrom, start, end))

	except StopIteration:
		break

output = open(output_name, 'a')
for i in sequence_list:
	output.write("%s\n" %i)

output.close()
bed.close()
tb.close()