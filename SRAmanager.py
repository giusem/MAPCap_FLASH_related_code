import subprocess
import glob

dump_path = "/package/sratoolkit/bin/fastq-dump"


with open ("SRAin.txt", 'r') as ff:
	for line in ff:
		print line
		#name = line.split("/")[-1]

		print "getting SRA file\n"
		cmd = "%s" %(line) 
		subprocess.check_call(cmd, shell=True)

print "transforming file into fastq.gz"
sra = glob.glob("*.sra")
for i in sra:
	name = i.split(".")[0]
	print name

	cmd2 = "%s --gzip -A %s %s " %(dump_path, name, i)
	subprocess.check_call(cmd2, shell=True)

		#print "making fastq.gz file\n"
		#cmd2 = "%s -A /data/akhtar/group/Giuseppe/05-20160229_Msl3_patients_polyA_RNAseq/06-public_data_controls/%s /data/akhtar/group/Giuseppe/05-20160229_Msl3_patients_polyA_RNAseq/06-public_data_controls/%s -O /data/akhtar/group/Giuseppe/05-20160229_Msl3_patients_polyA_RNAseq/06-public_data_controls/ | gzip /data/akhtar/group/Giuseppe/05-20160229_Msl3_patients_polyA_RNAseq/06-public_data_controls/%s " %(dump_path, name.split(".")[0], name, name.split(".")[0]+".fastq" )
		#subprocess.check_call(cmd2, shell=True)
		
