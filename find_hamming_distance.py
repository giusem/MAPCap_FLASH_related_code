
import Levenshtein

print "input the sequence you want to compare to the reference barcodes"

query = raw_input(">")

ref_barcodes = ["CAAGTG", "GTCATG", "TTAGCC", "GTGGAA", "TGTGAG", "TGGAAC", "GGTTAC", "CATCAC", "AGAAGG", "CGCATA", "AGTCGA", "AACTCG"]

scores = []

for i in ref_barcodes:
	scores.append(Levenshtein.hamming(i, query))

for i in range(len(ref_barcodes)):
	print "%s\t%s" %(ref_barcodes[i], scores[i])