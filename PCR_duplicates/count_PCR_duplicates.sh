#!/bin/bash 

virtual #load virtual environment for matplotlib etc.

echo "welcome to Giuseppes PCR duplicate toolkit"
echo "loading scripts"

preparePCRdups=/data/akhtar/group/Giuseppe/Scripts/python_scripts/PCR_duplicates/preparePCRdups.py
printPCRdups=/data/akhtar/group/Giuseppe/Scripts/python_scripts/PCR_duplicates/printPCRdups.py
countPCRdup_per_loci=/data/akhtar/group/Giuseppe/Scripts/python_scripts/PCR_duplicates/countPCRdup_per_loci.py
printPCRdups_per_loci=/data/akhtar/group/Giuseppe/Scripts/python_scripts/PCR_duplicates/printPCRdups_per_loci.py
count_PCRdups_per_gene=/data/akhtar/group/Giuseppe/Scripts/python_scripts/PCR_duplicates/count_PCRdups_per_gene.py

for i in $preparePCRdups $printPCRdups $countPCRdup_per_loci $printPCRdups_per_loci $count_PCRdups_per_gene; do
	echo $i 
	if [ -e $i ]; then
		echo "found"
	else echo "not found"
		exit 1
	fi
done

echo "relative path of the BAM file"
read BAMname
bampath=$(realpath $BAMname)

echo "name of output (no extensions)"
read out_name

echo "was bamToBed already performed? answer with yes or no"
read status

if [ $status == "no" ]; then
	echo "preparing for bamToBed"
	if [ -e $bampath ]; then
		echo "BAM file found: $bampath"
	else echo "BAMfile not found at $bampath ... quitting..."
		exit 1
	fi
	current_dir=$(pwd)
	bedname=$current_dir/${out_name}.bed
	echo "bamToBed"
	bamToBed -i $bampath > $bedname
elif [ $status == "yes" ]; then
	echo "skipping bamToBed. looking for bed file now."
	current_dir=$(pwd)
	bedname=$current_dir/${out_name}.bed
	if [ -e $bedname ]; then
		echo "bed file found: $bedname"
	else echo "bed file not found at: $bedname. either bamToBed was not performed or the ouput name you provided is not the same as the name of the bed file...quitting..."
		exit 1
	fi
fi
pyoutname=$current_dir/${out_name}_temp.txt
echo "preparing data for counting"
python $preparePCRdups $bedname $pyoutname

echo "counting"

sort $pyoutname | uniq -c | sort -k1,1nr | awk '{print$2,$3,$4,$1-1}' | tr [:blank:] \\\t > ${out_name}_PCRdup_counts_per_barcode.txt
#the -1 is to remove one read as the original and keeping the rest as pcr duplicates

#echo "preparing for plotting"
#awk '{if($1 == "2L") print $0}' temp_counted.txt > temp_counted_2L.txt
#awk '{if($1 == "2R") print $0}' temp_counted.txt > temp_counted_2R.txt
#awk '{if($1 == "3L") print $0}' temp_counted.txt > temp_counted_3L.txt
#awk '{if($1 == "3R") print $0}' temp_counted.txt > temp_counted_3R.txt
#awk '{if($1 == "X") print $0}' temp_counted.txt > temp_counted_X.txt

#cat temp_counted_2L.txt  >> ${out_name}_counts.txt
#cat temp_counted_2R.txt  >> ${out_name}_counts.txt
#cat temp_counted_3L.txt  >> ${out_name}_counts.txt
#cat temp_counted_3R.txt  >> ${out_name}_counts.txt
#cat temp_counted_X.txt  >> ${out_name}_counts.txt

#awk '{print$1,$2,$3,$4}' temp_counted_all.txt | tr [:blank:] \\\t > ${out_name}_counts.txt

echo "plotting"
python $printPCRdups ${out_name}_PCRdup_counts_per_barcode.txt $current_dir/${out_name}_PCRdups_counts_per_barcode_PLOT.png

#rm temp_counted_2L.txt 
#rm temp_counted_2R.txt
#rm temp_counted_3L.txt
#rm temp_counted_3R.txt
#rm temp_counted_X.txt
#rm temp_counted.txt
rm $pyoutname

echo "preparing data for counting PCRdups per loci"
sort -k1,1 -k2,2n ${out_name}_PCRdup_counts_per_barcode.txt > temp.txt

python $countPCRdup_per_loci temp.txt ${out_name}_PCRdup_counts_per_loci.txt

rm temp.txt

echo "plotting PCRdups per loci"
python $printPCRdups_per_loci ${out_name}_PCRdup_counts_per_loci.txt ${out_name}_PCRdup_counts_per_loci_PLOT.png

echo "making bedgraph of the PCR duplicates"
sort -k1,1 -k2,2n ${out_name}_PCRdup_counts_per_loci.txt | awk '{if($4 != 0) print$0}' | tr [:blank:] \\\t > ${out_name}_PCRdups.bedGraph

echo "count PCRdups per gene"
sort -k1,1 -k2,2n ${out_name}_PCRdup_counts_per_loci.txt | awk '{print$1,$2,$2+1,$3}' | tr [:blank:] \\\t > newtemp.bed

echo "reference file used:"
echo "/data/akhtar/group/Giuseppe/supplementary/dm6/BED/genes_flat_longest_isoform/dm6_genes_all.bed"

closestBed -a newtemp.bed -b /data/akhtar/group/Giuseppe/supplementary/dm6/BED/genes_flat_longest_isoform/dm6_genes_all.bed > outtemp.bed

python $count_PCRdups_per_gene outtemp.bed pyout.txt

sort -k2,2nr pyout.txt > ${out_name}_PCRdup_counts_per_gene.txt

rm newtemp.bed
rm outtemp.bed
rm pyout.txt

echo "done, have a nice day"
