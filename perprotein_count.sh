#DONE--pull ncbiprot, kegggene(eg. tca) and kegggene(eg.tca), pathway
#DONE--merge those two table in R (use all to get all XP lines) to get kegg_xp_to_path_tca.csv

#DONE--protein fasta -> XP list -> array

#$1 is gggsss
#$2 is kegg organism code
#$3 is ncbi accession list

echo "NCBI protein,KOBAS annotations,Kofamscan annotations(*),KEGG species annotations" > $1/per_protein_count_"$2".csv
readarray -t XParray < $3
for each in "${XParray[@]}"
	do
	#grep -c XP kobas
	kobascnt=$(grep -c $each ../../clusterProfiler/$1/kobas_kegg_paths.txt)
	#grep -c XP clusterprofiler
	kofamcnt=$(grep -c $each ../../clusterProfiler/$1/asterisk/merged_"$1"_asterisk.csv)
	#grep -c XP from kegg_xp_to_path_tca.csv
	IFS=. read prot ver <<<"${each}"
	echo $prot
	keggcnt=$(grep -c $prot $1/kegg_xp_to_path_"$2".csv)
	#echo numbers to output file for species
	echo "$each,$kobascnt,$kofamcnt, $keggcnt" >>$1/per_protein_count_"$2".csv
done
