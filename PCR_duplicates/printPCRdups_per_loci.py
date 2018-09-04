#!/home/semplicio/virtualenv/bin/python -u


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from sys import argv

script, inputfile, outputfile = argv

infile = open(inputfile,'r')

countlist = []
while True:
	try:
	    line = infile.next().strip()
	    count = line.split("\t")[-1]
	    countlist.append(int(count)) 
	except StopIteration:
		break

x = []
y = []

#maxcount = int(max(countlist))+1

for i in range(100):
	x.append(i)
	y.append(countlist.count(i))

	with open(outputfile.split(".png")[0]+"counts.txt", 'a') as ff:
		ff.write("%s\t%s\n" %(i, countlist.count(i)))

plt.plot(x,y)
plt.title("Frequency of PCR duplicates per locus")
plt.xlabel("Number of PCR duplicates")
plt.ylabel("Number of genomic loci")
plt.savefig(outputfile)
plt.close



