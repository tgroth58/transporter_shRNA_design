#!/usr/bin/awk -f


BEGIN{FS="\n";OFS="\n"}
{ 
	if (index($0,">")!=0){
		split($0,geneName,">")
		print geneName[2]
	}
	print $0 >> geneName[2]"_seq.fa"
}
