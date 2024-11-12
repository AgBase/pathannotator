#! /bin/bash

#$1 is gggsss
#$2 is tsv output from kofamscan
#$3 is eval cutoff you want to use e.g. 0.05

#USE FOR CMD LINE $3 SPECIFIED CUTOFF
#> $1/filtered"$3"_"$2".txt

#USE FOR PREPROGGED CUTOFFS
> $1/"$2".txt
> $1/filtered0.05_"$2".txt
> $1/filtered0.01_"$2".txt
> $1/filtered_asterisk_"$2".txt
> $1/filtered_asterisk0.05_"$2".txt


sed 's/\t/!/g' $1/$2 > $1/"$2".tmp
sed -i 's/$!\+//g' $1/"$2".tmp

cat $1/"$2".tmp | while IFS='!' read -r asterisk gene_name KO thrshld score E_value KO_definition
do
	decEval=$(echo "$E_value" | awk -F"e" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}') #makes notation into decimal
#=================OPTION 1 DOES WHATEVER YOU SPECIFY FOR EVAL CUTOFF===========
#BECAUSE BASH CANT DO DECIMALS, USE bc. BECAUSE bc CANT DO SCIENTIFIC NOTATION USE awk
#USE TO SPECIFIY EVAL AS $3 ON CMD LINE
#	if (( $( bc -l<<<"$decEval <= $3" ) )); #compare decimal to decimal
#	then
#		echo $decEval
#		echo -e "$asterisk $gene_name $KO $thrshld $score $E_value $KO_definition" >> $1/filtered"$3"_"$2".txt
#	fi
#============END OPTION 1 BEGIN OPTION 2=========================================
#USE TO PROCESS 0.05, 0.01, *, *+0.05 AND FULL (WITH SPACES) AT ONCE
	#FULL RESULTS TABS TO SPACES
	echo $decEval
	if (( $( bc -l<<<"$decEval > 0.05" ) )); #compare decimal to decimal
	then
	        echo -e "$asterisk $gene_name $KO $thrshld $score $E_value $KO_definition" >> $1/"$2".txt
	# EVALUE 0.01 AND TABS TO SPACES
	elif (( $( bc -l<<<"$decEval <= 0.01" ) )); #compare decimal to decimal
	then
		echo -e "$asterisk $gene_name $KO $thrshld $score $E_value $KO_definition" >> $1/filtered0.01_"$2".txt
		echo -e "$asterisk $gene_name $KO $thrshld $score $E_value $KO_definition" >> $1/filtered0.05_"$2".txt
		echo -e "$asterisk $gene_name $KO $thrshld $score $E_value $KO_definition" >> $1/"$2".txt
	# EVALUE 0.05 AND TABS TO SPACES
	else
		echo -e "$asterisk $gene_name $KO $thrshld $score $E_value $KO_definition" >> $1/filtered0.05_"$2".txt
		echo -e "$asterisk $gene_name $KO $thrshld $score $E_value $KO_definition" >> $1/"$2".txt
	fi
done

# EVALUE 0.05 AND ASTERISKED LINE
grep -P "^\*" $1/"$2".txt >> $1/filtered_asterisk_"$2".txt

# ALL ASTERISKED LINES ONLY
grep -P "^\*" $1/filtered0.05_"$2".txt >> $1/filtered_asterisk0.05_"$2".txt

#=====================END OPTION 2===========================

rm $1/"$2".tmp
