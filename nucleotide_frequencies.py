
#this takes a bunch of sequences and outputs the nucleotide frequency for each position
#sequences need same length and should only contain "a", "C", "G", "T"
import sys

filename = sys.argv[1]

f = open(filename, 'r')

outname = filename.split(".txt")[0]+"_nucleotide_frequency.txt"

nLines = 0
for i in f:
	seqlength = len(i)-1
	nLines += 1

f.close()

f = open(filename, 'r')

hl = []
for i in range(seqlength):
	hl.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})

for line in f:
	for idx, c in enumerate(line.strip()):
		hl[idx][c] += 1
f.close()

nLines = float(nLines)

with open(outname, 'a') as ff:
	for char in ['A', 'C', 'G', 'T']:
		ff.write("{}\t{}\n".format(char, "\t".join(["{:0.2f}".format(x[char]/nLines) for x in hl])))

hl = []


