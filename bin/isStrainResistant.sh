vcf=$1
ref_mutation=$2

cut -f3 $ref_mutation -d ',' | sort -u > mutation_list

while read mutation
do
	gene=`echo $mutation | cut -f1 -d ','`
	mut=`echo $mutation | cut -f2 -d ','`

	if (grep "$gene" $vcf | grep -q $mut)
	then
		echo $mutation | cut -f3 -d ','
	fi
done < $ref_mutation | sort -u > sample_mutation

comm sample_mutation mutation_list -12 | awk '{if(NF){print $0 "\tTRUE"}}' > resistance
comm sample_mutation mutation_list -13 | awk '{if(NF){print $0 "\tFALSE"}}' > sensitivity

cat resistance sensitivity | sort -u > resistance_list

echo -en "sample\t"
cat resistance_list | cut -f1 | tr '\n' '\t' | sed 's/.$//'
printf "\n$(basename $vcf .vcf)\t"
cat resistance_list | cut -f2 | tr '\n' '\t' | sed 's/.$//'
echo ""
rm -rf mutation_list sample_mutation resistance sensitivity resistance_list
