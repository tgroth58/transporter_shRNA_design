#!/bin/bash
#shRNA tool
shRNATool=/home/ted/Documents/Programs/si_shRNA_design/si_shRNA_selector-x64
cat random_dna* > total_random_seqs.fa
$shRNATool -i total_random_seqs.fa -o total_random_seqs_designs.txt
sed -i '1s/.*/&\tefficiency/'  total_random_seqs_designs.txt
awk ' 
BEGIN{count=0}
$0~/efficient/{
	if ($4<=(-28) && $4>=(-33))
		{count=count+1}
}
END{print "There are "count" efficient shRNAs in the random sequences"}
' < total_random_seqs_designs.txt

