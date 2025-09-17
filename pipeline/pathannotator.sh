#! /bin/bash

#CHECK FOR OUTDIR. IF IT DOESN'T EXIST CREATE IT
if [ -f "$outdir"/link_ko_pathway.tsv ]; then rm "$outdir"/link_ko_pathway.tsv; fi
if [ -f "$outdir"/list_pathway.tsv ]; then rm "$outdir"/list_pathway.tsv; fi
if [ -f "$outdir"/conv_ncbi-proteinid_"$keggcode".tsv ]; then rm "$outdir"/conv_ncbi-proteinid_"$keggcode".tsv; fi
if [ -f "$outdir"/link_"$keggcode"_ko.tsv ]; then rm "$outdir"/link_"$keggcode"_ko.tsv; fi
if [ -f "$outdir"/link_pathway_"$keggcode".tsv ]; then rm "$outdir"/link_pathway_"$keggcode".tsv; fi
if [ -f "$outdir"/list_pathway_"$keggcode".tsv ]; then rm "$outdir"/list_pathway_"$keggcode".tsv; fi
if [ -f "$outdir"/deflines.tmp ]; then rm "$outdir"/deflines.tmp; fi
if [ -f "$outdir"/ko_ncbi.tsv ]; then rm "$outdir"/ko_ncbi.tsv; fi
if [ -f "$outdir"/Fbgn_CG.tsv ]; then rm "$outdir"/Fbgn_CG.tsv; fi
if [ -f "$outdir"/Fbgn_groupid.tsv ]; then rm "$outdir"/Fbgn_groupid.tsv; fi
if [ -f "$outdir"/pathway_group_data_latest.tsv ]; then rm $outdir/pathway_group_data_latest.tsv; fi
if [ -f "$outdir"/kofam_filtered_asterisk.txt ]; then rm "$outdir"/kofam_filtered_asterisk.txt; fi
if [ -f "$outdir"/kegg_organisms.txt ]; then rm "$outdir"/kegg_organisms.txt; fi
if [ -f "$outdir"/kegg_org_codes.txt ]; then rm "$outdir"/kegg_org_codes.txt; fi
if [ -f "$outdir"/kegg_orgs_with_codes.txt ]; then rm "$outdir"/kegg_orgs_with_codes.txt; fi
if [ -n "$(ls $outdir/*pathway_group_data_fb* 2>/dev/null)" ]; then rm $outdir/*pathway_group_data_fb*; fi
if [ -n "$(ls $outdir/fbgn_annotation_ID_fb* 2>/dev/null)" ]; then rm $outdir/fbgn_annotation_ID_fb*; fi
if [ -n "$(ls $outdir/dmel-all-translation*.fasta* 2>/dev/null)" ]; then rm $outdir/dmel-all-translation*.fasta*; fi
if [ -n "$(ls $outdir/fbgn_fbtr_fbpp_fb* 2>/dev/null)" ]; then rm $outdir/fbgn_fbtr_fbpp_fb*; fi
if [ -f "$outdir"/Fbgn_fbpp.tsv ]; then rm "$outdir"/Fbgn_fbpp.tsv; fi
if [ -d "$outdir"/tmp ]; then rm -r "$outdir"/tmp; fi
if [ -f "$outdir"/tmp.txt ]; then rm  "$outdir"/tmp.txt; fi
if [ -f "$outdir"/ncbiversion.tmp ]; then rm "$outdir"/ncbiversion.tmp; fi
if [ -f "$outdir"/ncbiver.tsv ]; then rm "$outdir"/ncbiver.tsv; fi
if [ -d "$outdir"/orthofinder ]; then rm -r "$outdir"/orthofinder; fi
if [ -f "$outdir"/UniProt2Reactome_DME.txt ]; then rm "$outdir"/UniProt2Reactome_DME.txt; fi
if [ -n "$(ls $outdir/gp_information.* 2>/dev/null)" ]; then rm "$outdir"/gp_information.*; fi

starttime=$(date +%s)

############################################################################################################################
#SETUP ARGS

while getopts 'd:f:i:k:o:h' option
do
  case "${option}" in
    k) keggcode=${OPTARG};;
    i) input=${OPTARG};;
    d) outdir=${OPTARG};;
    f) flybase=${OPTARG};;
    o) outbase=${OPTARG};;
    h) help=true;;
    \?) echo "No legal parameters were passed. Please run with -h parameter to see help"; exit 1;;
esac
done

##############################################################################################################################

if [[ "$help" = "true" ]]
then
	echo "Help and Usage:
	-h to see this help and usage statement
	-k KEGG species code (NA or related species code if species not in KEGG; default is 'NA')
	   KEGG species codes can be found here: https://www.genome.jp/brite/br08611
	-i input file (protein FASTA without header lines)
	-d (optional: default is '.') output directory
	-f (optional: default is 'NA') Must be either: 'FB' for flybase and DME Reactome annotations or 'NA' for none
	-o outbase (file basename to use for output files)

	KofamScan is used under an MIT License:

	Copyright (c) 2019 Takuya Aramaki

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE."

	exit 0
fi
#######################################################################################################
#SET DEFAULTS IF OPTIONS NOT PROVIDED
if [ -z "${flybase}" ]; then $flybase == 'NA'; fi
if [ -z "${outdir}" ]; then $outdir == '.'; fi
if [ -z "${keggcode}" ]; then $keggcode == 'NA'; fi

if [ ! -d "$outdir" ]; then mkdir -p "$outdir"; fi

#GETTING NUMBER OF AVAILABLE PROCESSORS FOR USE IN THREADING
avail=$(getconf _NPROCESSORS_ONLN)
cpus=$(( $avail - 1 ))


#TESTS WHETHER ACCESSIONS ARE NCBI PROTEIN IDS
acc1=$(head -n1 $input | sed 's/>//g' | sed 's/\s.*$//')
if  [[ $acc1 == NP_* ]] || [[ $acc1 == XP_* ]] || [[ $acc1 == YP_* ]];
then
	ncbi=true
	echo "$acc1 These are NCBI protein IDs."
else
	ncbi=false
fi


#PULLS THE KEGG ORG CODES FILE (NEEDS TO BE IN HERE, NOT PULL_DATA.SH BECAUSE PULL DATA ONLY RUNS IN THE IF STATEMENTS BELOW)
wget https://rest.kegg.jp/list/genome -O $outdir/kegg_organisms.txt
grep ';' $outdir/kegg_organisms.txt > $outdir/kegg_orgs_with_codes.txt
cut -f 2 $outdir/kegg_orgs_with_codes.txt > $outdir/kegg_org_codes.txt
sed -i 's/;.*$//g' $outdir/kegg_org_codes.txt


if [ "$ncbi" == true ] ;
then
	#TAKES FASTA AND CREATES ACCESSION LIST. ACCESSION IS EVERYTHING BEFORE THE FIRST SPACE
	grep ">" $input > $outdir/deflines.tmp
	sed -i 's/>//g' $outdir/deflines.tmp
	sed -i 's/\s.*$//' $outdir/deflines.tmp
	# ADD A TWO COLUMN FILE OR ASSOC ARRAY OF WITH AND WITHOUT VERSION HERE THAT CAN BE USED TO MERGE LATER
	awk 'BEGIN {OFS="\t"} {print $1, $1}' $outdir/deflines.tmp > $outdir/ncbiversion.tmp
	awk 'BEGIN {OFS="\t"} { sub(/\.[0-9]+/, "", $2) }1' $outdir/ncbiversion.tmp > $outdir/ncbiver.tsv
	sed -i 's/.[0-9]$//' $outdir/deflines.tmp
	readarray -t defarray < $outdir/deflines.tmp

	if grep -q $keggcode $outdir/kegg_org_codes.txt; #IF THIS IS A KEGG SPECIES
	then
		#PULL DATA
		echo "This is a KEGG species code. Pulling KEGG API data now."
		cp /FB/* $outdir/
		bash /usr/bin/pull_data.sh $keggcode no $outdir ncbi $flybase

		#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
		if [[ -s "$outdir/link_ko_pathway.tsv" && -s "$outdir/conv_ncbi-proteinid_"$keggcode".tsv" && -s "$outdir/list_pathway.tsv" && -s "$outdir/link_pathway_"$keggcode".tsv" && -s "$outdir/list_pathway_"$keggcode".tsv" && -s "$outdir/link_"$keggcode"_ko.tsv" ]]
		then
    			echo "All KEGG files exist and are not empty."
		else
    			echo "One or more of the specified files are empty or do not exist."
			exit
		fi

		#NEED TO COMPARE DEFLINES.TMP TO SPECIFIED SPECIES CODE AND DECIDE IF THEY ARE THE SAME SPECIES
		echo "${defarray[0]}"
		if grep -q "${defarray[0]}" $outdir/conv_ncbi-proteinid_"$keggcode".tsv; #TESTING IF INPUT IDS ARE THE SAME SPECIES AS THE KEGG CODE
		then
			#IF YES, MERGE FROM API DATA
			echo "IDs are $keggcode species IDs"

			#IF FB AND NOT 'DME' RUN ORTHOFINDER AND PROCEED TO MERGE (INCLUDING FLYBASE)
			if [ "$keggcode" != "dme" ] && [ "$flybase" == "FB" ];
			then
				echo "Performing Flybase annotation".
				mkdir $outdir/orthofinder

				#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
				if [[ -s "$outdir/Fbgn_groupid.tsv" && -s "$outdir/Fbgn_CG.tsv" && -s "$outdir/Fbgn_fbpp.tsv" && -s "$outdir/gp_information.fb" && -s "$outdir/UniProt2Reactome_DME.txt" ]]
				then
    					echo "All FlyBase and Reactome files exist and are not empty."
				else
    					echo "One or more of the specified files are empty or do not exist."
					exit
				fi

				ext="*.faa"
				if [[ $input == $ext ]];
				then
					echo FAA
					noext=$(basename "$input" .faa)
				else
					ext="*.fasta"
					if [[ $input == $ext ]];
					then
						echo FASTA
						noext=$(basename "$input" .fasta)
					else
						ext="*.fa"
						if [[ $input == $ext ]];
						then
							echo FA
                                       			noext=$(basename "$input" .fa)
						else
							echo -e "This FASTA input file does not have an appropriate extension (.fa, .faa, .fasta)"
						fi
					fi
				fi
				echo -e "noext is: $noext"

				cp $input $outdir/orthofinder/"$noext"_cluster.fa

				#RUN ORTHOFINDER WITH FASTAS FROM INPUT SPECIES AND REFERENCE SET
				tar -xvzf /OF/ref_set.tgz -C $outdir/orthofinder
				orthofinder -f $outdir/orthofinder -t $cpus -b $outdir/orthofinder/ref_set

				#MOVE THE Orthologues_dromel_cluster DIR UP TO orthofinder
				mv $outdir/orthofinder/ref_set/OrthoFinder/Results_*/Orthologues/Orthologues_"$noext"_cluster/ $outdir/orthofinder/

			fi

			#MERGE DATA HERE
			echo "Creating annotations output."
			python /usr/bin/merge_data.py $keggcode no $outdir $outdir $flybase $outdir/orthofinder/Orthologues_"$noext"_cluster/"$noext"_cluster__v__dromel_cluster.tsv $outbase

			#CREATE GMT FILE
			python /usr/bin/pathannot_to_gmt.py $outdir/ $outdir/ $outbase
		else
			#IF NO, THEN RUN KOFAM, FILTER, FB, MERGE FROM KOFAM DATA
			echo "IDs are NOT $keggcode species IDs"

			#PULL ADDITIONAL DATA FOR KOFAMSCAN
			echo "Pulling more KEGG API data now."
			cp /FB/* $outdir/
			bash /usr/bin/pull_data.sh $keggcode yes $outdir ncbi $flybase

			#RUN KOFAMSCAN
			echo "Running KofamScan now."
			/usr/bin/kofam_scan/exec_annotation -o $outdir/kofam_result_full.txt -f detail --tmp-dir $outdir/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $input

			#FILTER KOFAM HERE
			echo "Filtering KofamScan results"
			grep -P "^\*" $outdir/kofam_result_full.txt >> $outdir/kofam_filtered_asterisk.txt
	        	awk '{ print $3"\t"$2 }' $outdir/kofam_filtered_asterisk.txt > $outdir/ko_ncbi.tsv
	        	sed -i 's/\..*$//' $outdir/ko_ncbi.tsv

			#IF FB AND NOT 'DME' RUN ORTHOFINDER AND PROCEED TO MERGE (INCLUDING FLYBASE)
			if [ "$keggcode" != "dme" ] && [ "$flybase" == "FB" ];
			then
				echo "Performing Flybase annotation".
				mkdir $outdir/orthofinder

				#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
				if [[ -s "$outdir/Fbgn_groupid.tsv" && -s "$outdir/Fbgn_CG.tsv" && -s "$outdir/Fbgn_fbpp.tsv" && -s "$outdir/gp_information.fb" && -s "$outdir/UniProt2Reactome_DME.txt" ]]
				then
    					echo "All FlyBase and Reactome files exist and are not empty."
				else
    					echo "One or more of the specified files are empty or do not exist."
					exit
				fi

				ext="*.faa"
				if [[ $input == $ext ]];
				then
					echo FAA
					noext=$(basename "$input" .faa)
				else
					ext="*.fasta"
					if [[ $input == $ext ]];
					then
						echo FASTA
						noext=$(basename "$input" .fasta)
					else
						ext="*.fa"
						if [[ $input == $ext ]];
						then
							echo FA
                                        		noext=$(basename "$input" .fa)
						else
							echo -e "This FASTA input file does not have an appropriate extension (.fa, .faa, .fasta)"
						fi
					fi
				fi
				echo -e "noext is: $noext"

				cp $input $outdir/orthofinder/"$noext"_cluster.fa

				#RUN ORTHOFINDER WITH FASTA FROM INPUT SPECIES AND REFERENCE SET
				tar -xvzf /OF/ref_set.tgz -C $outdir/orthofinder
				orthofinder -f $outdir/orthofinder -t $cpus -b $outdir/orthofinder/ref_set

				#MOVE THE Orthologues_dromel_cluster DIR UP TO orthofinder
				mv $outdir/orthofinder/ref_set/OrthoFinder/Results_*/Orthologues/Orthologues_"$noext"_cluster/ $outdir/orthofinder/

			fi

			#MERGE DATA HERE
			echo "Creating annotations output."
			python /usr/bin/merge_data.py $keggcode yes $outdir $outdir $flybase $outdir/orthofinder/Orthologues_"$noext"_cluster/"$noext"_cluster__v__dromel_cluster.tsv $outbase

			#CREATE GMT FILE
			python /usr/bin/pathannot_to_gmt.py $outdir/ $outdir/ $outbase

		fi

	else # ELSE MEANS THE THE CODE IS NOT A KEGG SPECIES CODE

		#PULL DATA
		echo "Pulling KEGG API data."
		cp /FB/* $outdir/
		bash /usr/bin/pull_data.sh $keggcode yes $outdir ncbi $flybase

		#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
		if [[ -s "$outdir/link_ko_pathway.tsv" && -s "$outdir/list_pathway.tsv" ]]
		then
    			echo "All KEGG files exist and are not empty."
		else
    			echo "One or more of the specified files are empty or do not exist."
			exit
		fi

		#RUN KOFAMSCAN
		echo "This is not a KEGG species code. Running KofamScan now."
		/usr/bin/kofam_scan/exec_annotation -o $outdir/kofam_result_full.txt -f detail --tmp-dir $outdir/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $input

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $outdir/kofam_result_full.txt >> $outdir/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $outdir/kofam_filtered_asterisk.txt > $outdir/ko_ncbi.tsv
	        sed -i 's/\..*$//' $outdir/ko_ncbi.tsv

		#IF FB AND NOT 'DME' RUN ORTHOFINDER AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$keggcode" != "dme" ] && [ "$flybase" == "FB" ];
		then
			echo "Performing Flybase annotation".
			mkdir $outdir/orthofinder

			#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
			if [[ -s "$outdir/Fbgn_groupid.tsv" && -s "$outdir/Fbgn_CG.tsv" && -s "$outdir/Fbgn_fbpp.tsv" && -s "$outdir/gp_information.fb" && -s "$outdir/UniProt2Reactome_DME.txt" ]]
			then
				echo "All FlyBase and Reactome files exist and are not empty."
			else
				echo "One or more of the specified files are empty or do not exist."
				exit
			fi


			ext="*.faa"
			if [[ $input == $ext ]];
			then
				echo FAA
				noext=$(basename "$input" .faa)
			else
				ext="*.fasta"
				if [[ $input == $ext ]];
				then
					echo FASTA
					noext=$(basename "$input" .fasta)
				else
					ext="*.fa"
					if [[ $input == $ext ]];
					then
						echo FA
                                       		noext=$(basename "$input" .fa)
					else
						echo -e "This FASTA input file does not have an appropriate extension (.fa, .faa, .fasta)"
					fi
				fi
			fi
			echo -e "noext is: $noext"

			cp $input $outdir/orthofinder/"$noext"_cluster.fa

			#RUN ORTHOFINDER WITH FASTA FROM INPUT SPECIES AND REFERENCE SET
			tar -xvzf /OF/ref_set.tgz -C $outdir/orthofinder
			orthofinder -f $outdir/orthofinder -t $cpus -b $outdir/orthofinder/ref_set

			#MOVE THE Orthologues_dromel_cluster DIR UP TO orthofinder
			mv $outdir/orthofinder/ref_set/OrthoFinder/Results_*/Orthologues/Orthologues_"$noext"_cluster/ $outdir/orthofinder/
		fi

		#MERGE DATA
		echo "Creating annotations output."
		python /usr/bin/merge_data.py $keggcode yes $outdir $outdir $flybase $outdir/orthofinder/Orthologues_"$noext"_cluster/"$noext"_cluster__v__dromel_cluster.tsv $outbase

		#CREATE GMT FILE
		python /usr/bin/pathannot_to_gmt.py $outdir/ $outdir/ $outbase

	fi

else #ELSE MEANS THESE ARE NOT NCBI PROTEIN IDS.

	echo "These are NOT NCBI protein IDs. Proceeding with KofamScan."
	#TAKES FASTA AND CREATES ACCESSION LIST. ACCESSION IS EVERYTHING BEFORE THE FIRST SPACE
	grep ">" $input > $outdir/deflines.tmp
	sed -i 's/>//g' $outdir/deflines.tmp
	sed -i 's/\s.*$//' $outdir/deflines.tmp
	#ADD A TWO COLUMN FILE OF WITH AND WITHOUT VERSION THAT CAN BE USED TO MERGE LATER
	awk 'BEGIN {OFS="\t"} {print $1, $1}' $outdir/deflines.tmp > $outdir/ncbiversion.tmp
	awk 'BEGIN {OFS="\t"} { sub(/\.[0-9]+/, "", $2) }1' $outdir/ncbiversion.tmp > $outdir/ncbiver.tsv

	if grep -q $keggcode $outdir/kegg_org_codes.txt;
	then
		echo "This is a KEGG species".

		#PULL DATA
		echo "Pulling KEGG API data."
		cp /FB/* $outdir/
		bash /usr/bin/pull_data.sh $keggcode yes $outdir non-ncbi $flybase

		#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
		if [[ -s "$outdir/link_ko_pathway.tsv" && -s "$outdir/list_pathway.tsv" && -s "$outdir/link_pathway_"$keggcode".tsv" && -s "$outdir/list_pathway_"$keggcode".tsv" && -s "$outdir/link_"$keggcode"_ko.tsv" ]]
		then
    			echo "All KEGG files exist and are not empty."
		else
    			echo "One or more of the specified files are empty or do not exist."
			exit
		fi


		#RUN KOFAM HERE
		/usr/bin/kofam_scan/exec_annotation -o $outdir/kofam_result_full.txt -f detail --tmp-dir $outdir/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $input

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $outdir/kofam_result_full.txt >> $outdir/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $outdir/kofam_filtered_asterisk.txt > $outdir/ko_ncbi.tsv
	        sed -i 's/\..*$//' $outdir/ko_ncbi.tsv

		#IF FB RUN ORTHOFINDER AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$flybase" == FB ];
		then
			echo "Performing Flybase annotation".
			mkdir $outdir/orthofinder

			#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
			if [[ -s "$outdir/Fbgn_groupid.tsv" && -s "$outdir/Fbgn_CG.tsv" && -s "$outdir/Fbgn_fbpp.tsv" && -s "$outdir/gp_information.fb" && -s "$outdir/UniProt2Reactome_DME.txt" ]]
			then
				echo "All FlyBase and Reactome files exist and are not empty."
			else
				echo "One or more of the specified files are empty or do not exist."
				exit
			fi

			ext="*.faa"
			if [[ $input == $ext ]];
			then
				echo FAA
				noext=$(basename "$input" .faa)
			else
				ext="*.fasta"
				if [[ $input == $ext ]];
				then
					echo FASTA
					noext=$(basename "$input" .fasta)
				else
					ext="*.fa"
					if [[ $input == $ext ]];
					then
						echo FA
                                 		noext=$(basename "$input" .fa)
					else
						echo -e "This FASTA input file does not have an appropriate extension (.fa, .faa, .fasta)"
					fi
				fi
			fi
			echo -e "noext is: $noext"

			cp $input $outdir/orthofinder/"$noext"_cluster.fa

			#RUN ORTHOFINDER WITH FASTA FROM INPUT SPECIES AND REFSET
			tar -xvzf /OF/ref_set.tgz -C $outdir/orthofinder
			orthofinder -f $outdir/orthofinder -t $cpus -b $outdir/orthofinder/ref_set

			#MOVE THE Orthologues_dromel_cluster DIR UP TO orthofinder
			mv $outdir/orthofinder/ref_set/OrthoFinder/Results_*/Orthologues/Orthologues_"$noext"_cluster/ $outdir/orthofinder/
		fi

		#MERGE DATA
		echo "Creating annotation outputs."
		python /usr/bin/merge_data.py $keggcode yes $outdir $outdir $flybase $outdir/orthofinder/Orthologues_"$noext"_cluster/"$noext"_cluster__v__dromel_cluster.tsv $outbase

		#CREATE GMT FILE
		python /usr/bin/pathannot_to_gmt.py $outdir/ $outdir/ $outbase

	else #ELSE MEANS THIS IS NOT A KEGG SPECIES

		echo "This is not a KEGG species".

		#PULL DATA
		echo "Pulling KEGG API data."
		cp /FB/* $outdir/
		bash /usr/bin/pull_data.sh $keggcode yes $outdir non-ncbi $flybase

		#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
		if [[ -s "$outdir/link_ko_pathway.tsv" && -s "$outdir/list_pathway.tsv" ]]
		then
    			echo "All KEGG files exist and are not empty."
		else
    			echo "One or more of the specified files are empty or do not exist."
			exit
		fi

		#RUN KOFAM HERE
		/usr/bin/kofam_scan/exec_annotation -o $outdir/kofam_result_full.txt -f detail --tmp-dir $outdir/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $input

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $outdir/kofam_result_full.txt >> $outdir/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $outdir/kofam_filtered_asterisk.txt > $outdir/ko_ncbi.tsv
	        sed -i 's/\..*$//' $outdir/ko_ncbi.tsv

		#IF FB RUN ORTHOFINDER AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$flybase" == FB ];
		then
			echo "Performing Flybase annotation".
			mkdir $outdir/orthofinder

			#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
			if [[ -s "$outdir/Fbgn_groupid.tsv" && -s "$outdir/Fbgn_CG.tsv" && -s "$outdir/Fbgn_fbpp.tsv" && -s "$outdir/gp_information.fb" && -s "$outdir/UniProt2Reactome_DME.txt" ]]
			then
    				echo "All FlyBase and Reactome files exist and are not empty."
			else
    				echo "One or more of the specified files are empty or do not exist."
				exit
			fi

			ext="*.faa"
			if [[ $input == $ext ]];
			then
				echo FAA
				noext=$(basename "$input" .faa)
			else
				ext="*.fasta"
				if [[ $input == $ext ]];
				then
					echo FASTA
					noext=$(basename "$input" .fasta)
				else
					ext="*.fa"
					if [[ $input == $ext ]];
					then
						echo FA
                                 		noext=$(basename "$input" .fa)
					else
						echo -e "This FASTA input file does not have an appropriate extension (.fa, .faa, .fasta)"
					fi
				fi
			fi
			echo -e "noext is: $noext"

			cp $input $outdir/orthofinder/"$noext"_cluster.fa

			#RUN ORTHOFINDER WITH FASTAS FROM INPUT SPECIES AND REFERENCE SET
			tar -xvzf /OF/ref_set.tgz -C $outdir/orthofinder
			orthofinder -f $outdir/orthofinder -t $cpus -b $outdir/orthofinder/ref_set

			#MOVE THE Orthologues_dromel_cluster DIR UP TO orthofinder
			mv $outdir/orthofinder/ref_set/OrthoFinder/Results_*/Orthologues/Orthologues_"$noext"_cluster/ $outdir/orthofinder/
		fi

		#MERGE DATA
		echo "Creating annotation outputs."
		python /usr/bin/merge_data.py $keggcode yes $outdir $outdir $flybase $outdir/orthofinder/Orthologues_"$noext"_cluster/"$noext"_cluster__v__dromel_cluster.tsv $outbase

		#CREATE GMT FILE
		python /usr/bin/pathannot_to_gmt.py $outdir/ $outdir/ $outbase
	fi
fi

if [ -f "$outdir"/link_ko_pathway.tsv ]; then rm "$outdir"/link_ko_pathway.tsv; fi
if [ -f "$outdir"/list_pathway.tsv ]; then rm "$outdir"/list_pathway.tsv; fi
if [ -f "$outdir"/conv_ncbi-proteinid_"$keggcode".tsv ]; then rm "$outdir"/conv_ncbi-proteinid_"$keggcode".tsv; fi
if [ -f "$outdir"/link_"$keggcode"_ko.tsv ]; then rm "$outdir"/link_"$keggcode"_ko.tsv; fi
if [ -f "$outdir"/link_pathway_"$keggcode".tsv ]; then rm "$outdir"/link_pathway_"$keggcode".tsv; fi
if [ -f "$outdir"/list_pathway_"$keggcode".tsv ]; then rm "$outdir"/list_pathway_"$keggcode".tsv; fi
if [ -f "$outdir"/deflines.tmp ]; then rm "$outdir"/deflines.tmp; fi
if [ -f "$outdir"/ko_ncbi.tsv ]; then rm "$outdir"/ko_ncbi.tsv; fi
if [ -f "$outdir"/Fbgn_CG.tsv ]; then rm "$outdir"/Fbgn_CG.tsv; fi
if [ -f "$outdir"/Fbgn_groupid.tsv ]; then rm "$outdir"/Fbgn_groupid.tsv; fi
if [ -f "$outdir"/pathway_group_data_latest.tsv ]; then rm $outdir/pathway_group_data_latest.tsv; fi
if [ -f "$outdir"/kofam_filtered_asterisk.txt ]; then rm "$outdir"/kofam_filtered_asterisk.txt; fi
if [ -f "$outdir"/kofam_result_full.txt ]; then rm "$outdir"/kofam_result_full.txt; fi
if [ -f "$outdir"/kegg_organisms.txt ]; then rm "$outdir"/kegg_organisms.txt; fi
if [ -f "$outdir"/kegg_org_codes.txt ]; then rm "$outdir"/kegg_org_codes.txt; fi
if [ -f "$outdir"/kegg_orgs_with_codes.txt ]; then rm "$outdir"/kegg_orgs_with_codes.txt; fi
if [ -n "$(ls $outdir/*pathway_group_data_fb* 2>/dev/null)" ]; then rm $outdir/*pathway_group_data_fb*; fi
if [ -n "$(ls $outdir/fbgn_annotation_ID_fb* 2>/dev/null)" ]; then rm $outdir/fbgn_annotation_ID_fb*; fi
if [ -n "$(ls $outdir/dmel-all-translation*.fasta* 2>/dev/null)" ]; then rm $outdir/dmel-all-translation*.fasta*; fi
if [ -n "$(ls $outdir/fbgn_fbtr_fbpp_fb* 2>/dev/null)" ]; then rm $outdir/fbgn_fbtr_fbpp_fb*; fi
if [ -f "$outdir"/Fbgn_fbpp.tsv ]; then rm "$outdir"/Fbgn_fbpp.tsv; fi
if [ -d "$outdir"/tmp ]; then rm -r "$outdir"/tmp; fi
if [ -f "$outdir"/tmp.txt ]; then rm  "$outdir"/tmp.txt; fi
if [ -d "$outdir"/orthofinder/ ]; then rm -r "$outdir"/orthofinder/; fi
if [ -f "$outdir"/ncbiversion.tmp ]; then rm "$outdir"/ncbiversion.tmp; fi
if [ -f "$outdir"/ncbiver.tsv ]; then rm "$outdir"/ncbiver.tsv; fi
if [ -f "$outdir"/UniProt2Reactome_DME.txt ]; then rm "$outdir"/UniProt2Reactome_DME.txt; fi
if [ -n "$(ls $outdir/gp_information.* 2>/dev/null)" ]; then rm "$outdir"/gp_information.*; fi

endtime=$(date +%s)
seconds=$(($endtime - $starttime))
runtime=$(($seconds / 60))
echo "Run time: $runtime minutes"
