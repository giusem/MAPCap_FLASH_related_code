
import gzip
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import re
import Levenshtein

def main():

	parser = argparse.ArgumentParser(description = "Checks if the reads have been wrongly reverse transcribed\n"
												   "by checking for the first 7nt of the reverse complement of\n"
												   "the illumina adaptor up to a hamming distance of 1.")

	parser.add_argument('-i', '--input', metavar = 'fastq file name', required = True, help = "Name of the fastq input file.")

	args = parser.parse_args()

	badRTprimerSite(args.input)


def badRTprimerSite(reverse):

	reverse = reverse.split("/")[-1]

	reverse_file = os.path.realpath(reverse)

	save_dir = reverse_file.split(reverse)[0]

	if re.search('fastq', reverse_file, re.IGNORECASE):
		if re.search('fastq.gz', reverse_file, re.IGNORECASE):
			file_type = 'fastq.gz'
		else:
			file_type = 'fastq'
	else:
		print("I don't understand the file type of %s. Your files should end either with .fastq or .fastq.gz" %input)
		sys.exit(1)


	linecount = 0

	if file_type == 'fastq.gz':
		with gzip.open(reverse_file, 'r') as ff:
			for f in ff:
				linecount += 1
	elif file_type == 'fastq':
		with open(reverse_file, 'r') as ff:
			for f in ff:
				linecount += 1

	readcount = int(linecount/4)

	if file_type == 'fastq.gz':
		fastq_r = gzip.open(reverse_file, 'r')
	elif file_type == 'fastq':
		fastq_r = open(reverse_file, 'r')

	count = 0
	badRTcount = 0
	badRTcount_hamdist = 0
	factor = 1
	status = 10

	allRTprimer_position = []
	RTprimer_position = []
	hamRTprimer_position = []

	while True:
		try:
			name_r	= fastq_r.next().strip()
			seq_r	= fastq_r.next().strip()
			plus_r	= fastq_r.next().strip()
			qual_r	= fastq_r.next().strip()
			count += 1

			if count >= ((readcount/10)*factor):
				print str(status)+"%"+" of the reads processed..."
				factor += 1
				status += 10



			if seq_r.find('CCGATCT') != -1:
				badRTcount += 1
				position = seq_r.find('CCGATCT') + 1
				allRTprimer_position.append(position)
				RTprimer_position.append(position)
				with open("badRTprimer.txt", 'a') as ff:
					ff.write(seq_r+"\n")

			else:
				a = 0
				b = 7
				ref_sequence = 'CCGATCT'


				for i in range(int(len(seq_r))-7):
					seq_r_short = seq_r[a:b]
					if Levenshtein.hamming(seq_r_short, ref_sequence) == 1:
						badRTcount_hamdist += 1
						allRTprimer_position.append(a+1)
						hamRTprimer_position.append(a+1)
						with open("badRTprimer_hamdist1.txt",'a') as ff:
							ff.write(seq_r+"\n")
						break

					else:
						a += 1
						b += 1

		except StopIteration:
			break

	all_bad_count = badRTcount + badRTcount_hamdist

	x = []
	a = []
	for i in range(50):
		x.append(i)
		a.append(allRTprimer_position.count(i))


	b = []
	for i in range(50):
		b.append(RTprimer_position.count(i))


	c = []
	for i in range(50):
		c.append(hamRTprimer_position.count(i))


	plt.subplot(311)
	plt.plot(x,a)
	plt.title("All RT primer")
	plt.ylabel("Counts")
	plt.subplot( 312)
	plt.plot(x,b)
	plt.title("RT primer w/o mutations")
	plt.ylabel("Counts")
	plt.subplot(313)
	plt.plot(x,c)
	plt.title("RT primer with a hamming distance of 1")
	plt.ylabel("Counts")
	plt.xlabel("Position")
	plt.savefig("RTprimer_position.png")
	plt.close


	with open("RTprimerposition_report.txt", 'a') as ff:
		msg1= "Here are some statistics:\n"
		msg2= "Number of reads: %s" %count
		msg3= "Number of reads where nothing wrong was found: %s (%%%s of all reads)" %(count-all_bad_count, round((count-all_bad_count)/count*100))
		msg4= "Number of reads where the primer site was shifted: %s (%%%s of all reads)" %(all_bad_count, round((all_bad_count/count)*100))
		msg5= "Of those %s reads, %s (%%%s) where found without any mutations and %s (%%%s) where found with one mutation." %(all_bad_count, badRTcount, round((badRTcount/all_bad_count)*100), badRTcount_hamdist, round((badRTcount_hamdist/all_bad_count)*100))

		ff.write("%s\n%s\n%s\n%s\n%s" %(msg1, msg2, msg3, msg4, msg5))

	print "Here are some statistics:\n"
	print "Number of reads: %s" %count
	print "Number of reads where nothing wrong was found: %s (%%%s of all reads)" %(count-all_bad_count, round((count-all_bad_count)/count*100))
	print "Number of reads where the primer site was shifted: %s (%%%s of all reads)" %(all_bad_count, round((all_bad_count/count)*100))
	print "Of those %s reads, %s (%%%s) where found without any mutations and %s (%%%s) where found with one mutation." %(all_bad_count, badRTcount, round((badRTcount/all_bad_count)*100), badRTcount_hamdist, round((badRTcount_hamdist/all_bad_count)*100))
	print "----------------------------\n"

	fastq_r.close()

if __name__ == "__main__":
	main()
