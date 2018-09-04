#!/home/semplicio/virtualenv/bin/python -u

import os
import re
import gzip
import subprocess
import argparse
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():

	parser = argparse.ArgumentParser(description = "Look for the barcode\n"
"report the position in a plot as well as in an output file\n"
"to run with single end reads, just load the same file for forward and reverse\n")

	parser.add_argument('-f', '--forward', metavar="forward read file", help='The fastq file that contains the forward reads.' ,required=True)
	parser.add_argument('-r', '--reverse', metavar="reverse read file ", help='The fastq file that contains the reverse reads.' ,required=True)

	args = parser.parse_args()

	adapterPositions(args.forward, args.reverse)

def adapterPositions(forward, reverse):

	forward = forward.split("/")[-1]
	reverse = reverse.split("/")[-1]


	forward_file = os.path.realpath(forward)
	reverse_file = os.path.realpath(reverse)


	save_dir = forward_file.split(forward)[0]

	if re.search('fastq', forward, re.IGNORECASE):
		if re.search('fastq.gz', forward, re.IGNORECASE):
			filetype_f = 'fastq.gz'
		else:
			filetype_f = 'fastq'
	else:
		print("I don't understand the file type of %s. Your files should end either with .fastq or .fastq.gz" %args.forward)
		sys.exit(1)

	if re.search('fastq', reverse, re.IGNORECASE):
		if re.search('fastq.gz', reverse, re.IGNORECASE):
			filetype_r = 'fastq.gz'
		else:
			filetype_r = 'fastq'
	else:
		print("I don't understand the file type of %s. Your files should end either with .fastq or .fastq.gz" %args.reverse)
		sys.exit(1)

	if filetype_f == 'fastq.gz':
		fastq_f = gzip.open(forward, 'rb')
	elif filetype_f == 'fastq':
		fastq_f = open(forward, 'rb')

	if filetype_r == 'fastq.gz':
		fastq_r = gzip.open(reverse, 'rb')
	elif filetype_r == 'fastq':
		fastq_r = open(reverse, 'rb')



	count = 0
	multiplier = 1
	chunk = 100000
	diff_count = 0
	same_count = 0
	noadapt_count = 0

	f_reads_dict = {}
	r_reads_dict = {}
	same_pos_dict = {}
	diff_pos_dict = {}
	no_adapt_dict_f = {}
	no_adapt_dict_r = {}

	adapt_pos = []
	x = []
	y = []


	while True:
		try:
			name_f, name_r	= fastq_f.next().strip(), fastq_r.next().strip()
			seq_f,	seq_r	= fastq_f.next().strip(), fastq_r.next().strip()
			plus_f, plus_r	= fastq_f.next().strip(), fastq_r.next().strip()
			qual_f, qual_r	= fastq_f.next().strip(), fastq_r.next().strip()

			count +=1

			pos_f = seq_f.find('TATGCG')
			pos_r = seq_r.find('TATGCG')

			newseq_f = seq_f.split("TATGCG")[0]+"--TATGCG--"+seq_f.split("TATGCG")[-1]
			newseq_r = seq_r.split("TATGCG")[0]+"--TATGCG--"+seq_r.split("TATGCG")[-1]

			if pos_f != pos_r:
				diff_count += 1
				try:
					f_reads_dict[save_dir+"adapter_diff_positions_R1.fastq"].append(["%s\n%s\n+\n%s\n" %(name_f, seq_f, qual_f)])
					r_reads_dict[save_dir+"adapter_diff_positions_R2.fastq"].append(["%s\n%s\n+\n%s\n" %(name_r, seq_r, qual_r)])
				except KeyError:
					f_reads_dict[save_dir+"adapter_diff_positions_R1.fastq"] = [["%s\n%s\n+\n%s\n" %(name_f, seq_f, qual_f)]]
					r_reads_dict[save_dir+"adapter_diff_positions_R2.fastq"] = [["%s\n%s\n+\n%s\n" %(name_r, seq_r, qual_r)]]

				if pos_f  == -1:
					try:
						diff_pos_dict[save_dir+"adapter_different_positions.txt"].append(["%s\n%s\t%s\t%s\n%s\t%s\t%s\n" %(name_f, count, seq_f, pos_f, count, newseq_r, pos_r)])
					except KeyError:
						diff_pos_dict[save_dir+"adapter_different_positions.txt"] = [["%s\n%s\t%s\t%s\n%s\t%s\t%s\n" %(name_f, count, seq_f, pos_f, count, newseq_r, pos_r)]]

				elif pos_r  == -1:
					try:
						diff_pos_dict[save_dir+"adapter_different_positions.txt"].append(["%s\n%s\t%s\t%s\n%s\t%s\t%s\n" %(name_f, count, newseq_f, pos_f, count, seq_r, pos_r)])
					except KeyError:
						diff_pos_dict[save_dir+"adapter_different_positions.txt"] = [["%s\n%s\t%s\t%s\n%s\t%s\t%s\n" %(name_f, count, newseq_f, pos_f, count, seq_r, pos_r)]]

				else:
					try:
						diff_pos_dict[save_dir+"adapter_different_positions.txt"].append(["%s\n%s\t%s\t%s\n%s\t%s\t%s\n" %(name_f, count, newseq_f, pos_f, count, newseq_r, pos_r)])
					except KeyError:
						diff_pos_dict[save_dir+"adapter_different_positions.txt"] = [["%s\n%s\t%s\t%s\n%s\t%s\t%s\n" %(name_f, count, newseq_f, pos_f, count, newseq_r, pos_r)]]

			elif pos_f == -1 and pos_r == -1:
				noadapt_count += 1
				try:
					no_adapt_dict_f[save_dir+"no_adapter_R1.fastq"].append(["%s\n%s\n%s\n%s\n" %(name_f, seq_f, plus_f, qual_f)])
					no_adapt_dict_r[save_dir+"no_adapter_R2.fastq"].append(["%s\n%s\n%s\n%s\n" %(name_r, seq_r, plus_r, qual_r)])
				except KeyError:
					no_adapt_dict_f[save_dir+"no_adapter_R1.fastq"] = [["%s\n%s\n%s\n%s\n" %(name_f, seq_f, plus_f, qual_f)]]
					no_adapt_dict_r[save_dir+"no_adapter_R2.fastq"] = [["%s\n%s\n%s\n%s\n" %(name_r, seq_r, plus_r, qual_r)]]

			elif pos_f == pos_r:
				same_count += 1
				adapt_pos.append(pos_f)
				try:
					f_reads_dict[save_dir+"adapter_same_positions_R1.fastq"].append(["%s\n%s\n+\n%s\n" %(name_f, seq_f, qual_f)])
					r_reads_dict[save_dir+"adapter_same_positions_R2.fastq"].append(["%s\n%s\n+\n%s\n" %(name_r, seq_r, qual_r)])
				except KeyError:
					f_reads_dict[save_dir+"adapter_same_positions_R1.fastq"] = [["%s\n%s\n+\n%s\n" %(name_f, seq_f, qual_f)]]
					r_reads_dict[save_dir+"adapter_same_positions_R2.fastq"] = [["%s\n%s\n+\n%s\n" %(name_r, seq_r, qual_r)]]

				try:
					same_pos_dict[save_dir+"adapter_same_positions.txt"].append(["%s\n%s\t%s\t%s\n%s\t%s\n" %(name_f, count, newseq_f, pos_f, count, newseq_r)])
				except KeyError:
					same_pos_dict[save_dir+"adapter_same_positions.txt"] = [["%s\n%s\t%s\t%s\n%s\t%s\n" %(name_f, count, newseq_f, pos_f, count, newseq_r)]]

			if count == multiplier*chunk:
				for key, value in f_reads_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				for key, value in r_reads_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				for key, value in diff_pos_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				for key, value in same_pos_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				for key, value in no_adapt_dict_f.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				for key, value in no_adapt_dict_r.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

				multiplier += 1
				f_reads_dict = {}
				r_reads_dict = {}
				diff_pos_dict = {}
				same_pos_dict = {}
				no_adapt_dict_f = {}
				no_adapt_dict_r = {}

		except StopIteration:
			for key, value in f_reads_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			for key, value in r_reads_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			for key, value in diff_pos_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			for key, value in same_pos_dict.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			for key, value in no_adapt_dict_f.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			for key, value in no_adapt_dict_r.iteritems():
				with open('%s' %key, 'a+') as ff:
					for i in range(len(value)):
						ff.write('%s' %value[i][0])

			f_reads_dict = {}
			r_reads_dict = {}
			diff_pos_dict = {}
			same_pos_dict = {}
			no_adapt_dict_f = {}
			no_adapt_dict_r = {}

			break

	fastq_f.close()
	fastq_r.close()

	if os.path.isfile(save_dir+"diff_positions_R1.fastq"):
		cmd = "gzip %s" %save_dir+"diff_positions_R1.fastq"
		cmd2 = "gzip %s" %save_dir+"diff_positions_R2.fastq"
		subprocess.check_call(cmd, shell=True)
		subprocess.check_call(cmd2, shell=True)

	if os.path.isfile(save_dir+"no_adapter_R1.fastq"):
		cmd = "gzip %s" %save_dir+"no_adapter_R1.fastq"
		cmd2 = "gzip %s" %save_dir+"no_adapter_R2.fastq"
		subprocess.check_call(cmd, shell=True)
		subprocess.check_call(cmd2, shell=True)

	if os.path.isfile(save_dir+"same_positions_R1.fastq"):
		cmd = "gzip %s" %save_dir+"same_positions_R1.fastq"
		cmd2 = "gzip %s" %save_dir+"same_positions_R2.fastq"
		subprocess.check_call(cmd, shell=True)
		subprocess.check_call(cmd2, shell=True)

	for i in range(75):
		x.append(i)
		y.append(adapt_pos.count(i))
		with open(save_dir+"adapter_position.txt", 'a') as ff:
			ff.write(str(i)+"\t"+str(adapt_pos.count(i))+"\n")

	plt.plot(x,y)
	plt.title("Adapter positions")
	plt.ylabel("counts")
	plt.ylim((0,max(y)+(max(y)/10)))
	plt.xlabel("position")
	plt.savefig(save_dir+"adapter_position.png")
	plt.close

	adapt_count = count - noadapt_count

	with open("adapterposition_stats.txt", 'a') as ff:
		msg1 = "Here are some statistics:\n"
		msg2 = "Number of reads: %s" %count
		msg3 = "Number of reads where no adapter sequence was found: %s (%%%s of all reads)" %(noadapt_count, round((noadapt_count/count)*100))
		msg4 = "Number of reads with adapter sequences found: %s (%%%s of all reads)" %(adapt_count, round((adapt_count/count)*100))
		msg5 = "Number of reads where the positions of the adapters in the forward and reverse read were the same: %s (%%%s of all reads)" %(same_count, round((same_count/count)*100))
		msg6 = "Number of reads where the positions of the adapters in the forward and reverse reads were NOT the same: %s (%%%s of all reads)" %(diff_count, round((diff_count/count)*100))

		ff.write("%s\n%s\n%s\n%s\n%s\n%s" %(msg1, msg2, msg3, msg4, msg5, msg6))

	print "Here are some statistics:\n"
	print "Number of reads: %s" %count
	print "Number of reads where no adapter sequence was found: %s (%%%s of all reads)" %(noadapt_count, round((noadapt_count/count)*100))
	print "Number of reads with adapter sequences found: %s (%%%s of all reads)" %(adapt_count, round((adapt_count/count)*100))
	print "Number of reads where the positions of the adapters in the forward and reverse read were the same: %s (%%%s of all reads)" %(same_count, round((same_count/count)*100))
	print "Number of reads where the positions of the adapters in the forward and reverse reads were NOT the same: %s (%%%s of all reads)" %(diff_count, round((diff_count/count)*100))

if __name__ == "__main__":
	main()