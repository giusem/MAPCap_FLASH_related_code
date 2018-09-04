#!/home/semplicio/virtualenv/bin/python -u


from sys import argv

script, inputfile, outputfile = argv

bed = open(inputfile,'r')

alldict = {}
while True:
	try:
		line = bed.next().strip()
		chr = line.split("\t")[0]
		strand = line.split("\t")[-1]
		if strand == "+":
			start = line.split("\t")[1]
		elif strand == "-":
			start =line.split("\t")[2]
		else:
			print "found an entry where the strand is not specified: \n%s\n" %line

		code = line.split("\t")[-3].split(":")[-1]

		try:
			alldict[outputfile].append(["%s\t%s\t%s\n" %(chr, start, code)])
		except KeyError:
			alldict[outputfile] = [["%s\t%s\t%s\n" %(chr, start, code)]]

	except StopIteration:
		for key, value in alldict.iteritems():
			with open('%s' %key, 'a+') as ff:
				for i in range(len(value)):
					ff.write('%s' %value[i][0])

		break
