#! /bin/bash

#CHECK FOR OUTDIR. IF IT DOESN'T EXIST CREATE IT
if [ ! -d "$3" ]; then mkdir -p "$3"; fi
if [ -f "$3"/link_ko_pathway.tsv ]; then rm "$3"/link_ko_pathway.tsv; fi
if [ -f "$3"/list_pathway.tsv ]; then rm "$3"/list_pathway.tsv; fi
if [ -f "$3"/conv_ncbi-proteinid_"$1".tsv ]; then rm "$3"/conv_ncbi-proteinid_"$1".tsv; fi
if [ -f "$3"/link_"$1"_ko.tsv ]; then rm "$3"/link_"$1"_ko.tsv; fi
if [ -f "$3"/link_pathway_"$1".tsv ]; then rm "$3"/link_pathway_"$1".tsv; fi
if [ -f "$3"/list_pathway_"$1".tsv ]; then rm "$3"/list_pathway_"$1".tsv; fi
if [ -f "$3"/deflines.tmp ]; then rm "$3"/deflines.tmp; fi
if [ -f "$3"/ko_ncbi.tsv ]; then rm "$3"/ko_ncbi.tsv; fi
if [ -f "$3"/Fbgn_CG.tsv ]; then rm "$3"/Fbgn_CG.tsv; fi
if [ -f "$3"/Fbgn_groupid.tsv ]; then rm "$3"/Fbgn_groupid.tsv; fi
if [ -f "$3"/pathway_group_data_latest.tsv ]; then rm $3/pathway_group_data_latest.tsv; fi
if [ -f "$3"/kofam_filtered_asterisk.txt ]; then rm "$3"/kofam_filtered_asterisk.txt; fi
if [ -f "$3"/kegg_organisms.txt ]; then rm "$3"/kegg_organisms.txt; fi
if [ -f "$3"/kegg_org_codes.txt ]; then rm "$3"/kegg_org_codes.txt; fi
if [ -f "$3"/kegg_orgs_with_codes.txt ]; then rm "$3"/kegg_orgs_with_codes.txt; fi
if [ -n "$(ls $3/*pathway_group_data_fb* 2>/dev/null)" ]; then rm $3/*pathway_group_data_fb*; fi
if [ -n "$(ls $3/fbgn_annotation_ID_fb* 2>/dev/null)" ]; then rm $3/fbgn_annotation_ID_fb*; fi
if [ -n "$(ls $3/dmel-all-translation*.fasta* 2>/dev/null)" ]; then rm $3/dmel-all-translation*.fasta*; fi
if [ -n "$(ls $3/fbgn_fbtr_fbpp_fb* 2>/dev/null)" ]; then rm $3/fbgn_fbtr_fbpp_fb*; fi
if [ -f "$3"/Fbgn_fbpp.tsv ]; then rm "$3"/Fbgn_fbpp.tsv; fi
if [ -d "$3"/tmp ]; then rm -r "$3"/tmp; fi
if [ -f "$3"/tmp.txt ]; then rm  "$3"/tmp.txt; fi
if [ -d "$3"/orthofinder/ref_set ]; then rm -r "$3"/orthofinder/ref_set; fi
if [ -f "$3"/orthofinder/*_cluster.fa* ]; then rm -r "$3"/orthofinder/*_cluster.fa*; fi

starttime=$(date +%s)

if [ $1 == "help" ];
then
	echo "Help and Usage:
	There are 4 positional arguments.
	1: KEGG species code (NA or related species code if species not in KEGG; 'help' to see this help and usage statement)
	   KEGG species codes can be found here: https://www.genome.jp/brite/br08611
	2: input file (protein FASTA without header lines)
	3: output directory (must be an existing directory)
	4: 'FB' for flybase annotations, 'NA' for none
	5: outbase (file basename to use for output files)

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
#GETTING NUMBER OF AVAILABLE PROCESSORS FOR USE IN THREADING
avail=$(getconf _NPROCESSORS_ONLN)
cpus=$(( $avail - 1 ))


#TESTS WHETHER ACCESSIONS ARE NCBI PROTEIN IDS
acc1=$(head -n1 $2 | sed 's/>//g' | sed 's/\s.*$//')
if  [[ $acc1 == NP_* ]] || [[ $acc1 == XP_* ]] || [[ $acc1 == YP_* ]];
then
	ncbi=true
	echo "$acc1 These are NCBI protein IDs."
else
	ncbi=false
fi


#PULLS THE KEGG ORG CODES FILE (NEEDS TO BE IN HERE, NOT PULL_DATA.SH BECAUSE PULL DATA ONLY RUNS IN THE IF STATEMENTS BELOW)
wget https://rest.kegg.jp/list/genome -O $3/kegg_organisms.txt
grep ';' $3/kegg_organisms.txt > $3/kegg_orgs_with_codes.txt
cut -f 2 $3/kegg_orgs_with_codes.txt > $3/kegg_org_codes.txt
sed -i 's/;.*$//g' $3/kegg_org_codes.txt

if [ "$ncbi" == true ] ;
then
	#TAKES FASTA AND CREATES ACCESSION LIST. ACCESSION IS EVERYTHING BEFORE THE FIRST SPACE
	grep ">" $2 > $3/deflines.tmp
	sed -i 's/>//g' $3/deflines.tmp
	sed -i 's/\s.*$//' $3/deflines.tmp
	# ADD A TWO COLUMN FILE OR ASSOC ARRAY OF WITH AND WITHOUT VERSION HERE THAT CAN BE USED TO MERGE LATER
	awk 'BEGIN {OFS="\t"} {print $1, $1}' $3/deflines.tmp > $3/ncbiversion.tmp
	awk 'BEGIN {OFS="\t"} { sub(/\.[0-9]+/, "", $2) }1' $3/ncbiversion.tmp > $3/ncbiver.tsv
	sed -i 's/.[0-9]$//' $3/deflines.tmp
	readarray -t defarray < $3/deflines.tmp

	if grep -q $1 $3/kegg_org_codes.txt; #IF THIS IS A KEGG SPECIES
	then
		#PULL DATA
		echo "This is a KEGG species code. Pulling KEGG API data now."
		bash /usr/bin/pull_data.sh $1 no $3 ncbi $4

		#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
		if [[ -s "$3/link_ko_pathway.tsv" && -s "$3/conv_ncbi-proteinid_"$1".tsv" && -s "$3/list_pathway.tsv" && -s "$3/link_pathway_"$1".tsv" && -s "$3/list_pathway_"$1".tsv" && -s "$3/link_"$1"_ko.tsv" ]]
		then
    			echo "All KEGG files exist and are not empty."
		else
    			echo "One or more of the specified files are empty or do not exist."
			exit
		fi

		#NEED TO COMPARE DEFLINES.TMP TO SPECIFIED SPECIES CODE AND DECIDE IF THEY ARE THE SAME SPECIES
		echo "${defarray[0]}"
		if grep -q "${defarray[0]}" $3/conv_ncbi-proteinid_"$1".tsv; #TESTING IF INPUT IDS ARE THE SAME SPECIES AS THE KEGG CODE
		then
			#IF YES, MERGE FROM API DATA
			echo "IDs are $1 species IDs"

			#IF FB AND NOT 'DME' RUN ORTHOFINDER AND PROCEED TO MERGE (INCLUDING FLYBASE)
			if [ "$1" != "dme" ] && [ "$4" == "FB" ];
			then
				echo "Performing Flybase annotation".
				mkdir $3/orthofinder

				#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
				if [[ -s "$3/Fbgn_groupid.tsv" && -s "$3/Fbgn_CG.tsv" && -s "$3/Fbgn_fbpp.tsv" ]]
				then
    					echo "All FlyBase files exist and are not empty."
				else
    					echo "One or more of the specified files are empty or do not exist."
					exit
				fi

				ext="*.faa"
				if [[ $2 == $ext ]];
				then
					echo FAA
					noext=$(basename "$2" .faa)
				else
					ext="*.fasta"
					if [[ $2 == $ext ]];
					then
						echo FASTA
						noext=$(basename "$2" .fasta)
					else
						ext="*.fa"
						if [[ $2 == $ext ]];
						then
							echo FA
                                       			noext=$(basename "$2" .fa)
						else
							echo -e "This FASTA input file does not have an appropriate extension (.fa, .faa, .fasta)"
						fi
					fi
				fi
				echo -e "noext is: $noext"

				cp $2 $3/orthofinder/"$noext"_cluster.fa

				#RUN ORTHOFINDER WITH FASTAS FROM INPUT SPECIES AND REFERENCE SET
				tar -xvzf /OF/ref_set.tgz -C $3/orthofinder
				orthofinder -f $3/orthofinder -t $cpus -b $3/orthofinder/ref_set

				#MOVE THE Orthologues_dromel_cluster DIR UP TO orthofinder
				mv $3/orthofinder/ref_set/OrthoFinder/Results_*/Orthologues/Orthologues_"$noext"_cluster/ $3/orthofinder/

			fi

			#MERGE DATA HERE
			echo "Creating annotations output."
			python /usr/bin/merge_data.py $1 no $3 $3 $4 $3/orthofinder/Orthologues_"$noext"_cluster/"$noext"_cluster__v__dromel_cluster.tsv $5

		else
			#IF NO, THEN RUN KOFAM, FILTER, FB, MERGE FROM KOFAM DATA
			echo "IDs are NOT $1 species IDs"

			#PULL ADDITIONAL DATA FOR KOFAMSCAN
			echo "Pulling more KEGG API data now."
			bash /usr/bin/pull_data.sh $1 yes $3 ncbi $4

			#RUN KOFAMSCAN
			echo "Running KofamScan now."
			/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

			#FILTER KOFAM HERE
			echo "Filtering KofamScan results"
			grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        	awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        	sed -i 's/\..*$//' $3/ko_ncbi.tsv

			#IF FB AND NOT 'DME' RUN ORTHOFINDER AND PROCEED TO MERGE (INCLUDING FLYBASE)
			if [ "$1" != "dme" ] && [ "$4" == "FB" ];
			then
				echo "Performing Flybase annotation".
				mkdir $3/orthofinder

				#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
				if [[ -s "$3/Fbgn_groupid.tsv" && -s "$3/Fbgn_CG.tsv" && -s "$3/Fbgn_fbpp.tsv" ]]
				then
    					echo "All FlyBase files exist and are not empty."
				else
    					echo "One or more of the specified files are empty or do not exist."
					exit
				fi

				ext="*.faa"
				if [[ $2 == $ext ]];
				then
					echo FAA
					noext=$(basename "$2" .faa)
				else
					ext="*.fasta"
					if [[ $2 == $ext ]];
					then
						echo FASTA
						noext=$(basename "$2" .fasta)
					else
						ext="*.fa"
						if [[ $2 == $ext ]];
						then
							echo FA
                                        		noext=$(basename "$2" .fa)
						else
							echo -e "This FASTA input file does not have an appropriate extension (.fa, .faa, .fasta)"
						fi
					fi
				fi
				echo -e "noext is: $noext"

				cp $2 $3/orthofinder/"$noext"_cluster.fa

				#RUN ORTHOFINDER WITH FASTA FROM INPUT SPECIES AND REFERENCE SET
				tar -xvzf /OF/ref_set.tgz -C $3/orthofinder
				orthofinder -f $3/orthofinder -t $cpus -b $3/orthofinder/ref_set

				#MOVE THE Orthologues_dromel_cluster DIR UP TO orthofinder
				mv $3/orthofinder/ref_set/OrthoFinder/Results_*/Orthologues/Orthologues_"$noext"_cluster/ $3/orthofinder/

			fi

			#MERGE DATA HERE
			echo "Creating annotations output."
			python /usr/bin/merge_data.py $1 yes $3 $3 $4 $3/orthofinder/Orthologues_"$noext"_cluster/"$noext"_cluster__v__dromel_cluster.tsv $5

		fi

	else # ELSE MEANS THE THE CODE IS NOT A KEGG SPECIES CODE

		#PULL DATA
		echo "Pulling KEGG API data."
		bash /usr/bin/pull_data.sh $1 yes $3 ncbi $4

		#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
		if [[ -s "$3/link_ko_pathway.tsv" && -s "$3/list_pathway.tsv" ]]
		then
    			echo "All KEGG files exist and are not empty."
		else
    			echo "One or more of the specified files are empty or do not exist."
			exit
		fi

		#RUN KOFAMSCAN
		echo "This is not a KEGG species code. Running KofamScan now."
		/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        sed -i 's/\..*$//' $3/ko_ncbi.tsv

		#IF FB AND NOT 'DME' RUN ORTHOFINDER AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$1" != "dme" ] && [ "$4" == "FB" ];
		then
			echo "Performing Flybase annotation".
			mkdir $3/orthofinder

			#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
			if [[ -s "$3/Fbgn_groupid.tsv" && -s "$3/Fbgn_CG.tsv" && -s "$3/Fbgn_fbpp.tsv" ]]
			then
				echo "All FlyBase files exist and are not empty."
			else
				echo "One or more of the specified files are empty or do not exist."
				exit
			fi


			ext="*.faa"
			if [[ $2 == $ext ]];
			then
				echo FAA
				noext=$(basename "$2" .faa)
			else
				ext="*.fasta"
				if [[ $2 == $ext ]];
				then
					echo FASTA
					noext=$(basename "$2" .fasta)
				else
					ext="*.fa"
					if [[ $2 == $ext ]];
					then
						echo FA
                                       		noext=$(basename "$2" .fa)
					else
						echo -e "This FASTA input file does not have an appropriate extension (.fa, .faa, .fasta)"
					fi
				fi
			fi
			echo -e "noext is: $noext"

			cp $2 $3/orthofinder/"$noext"_cluster.fa

			#RUN ORTHOFINDER WITH FASTA FROM INPUT SPECIES AND REFERENCE SET
			tar -xvzf /OF/ref_set.tgz -C $3/orthofinder
			orthofinder -f $3/orthofinder -t $cpus -b $3/orthofinder/ref_set

			#MOVE THE Orthologues_dromel_cluster DIR UP TO orthofinder
			mv $3/orthofinder/ref_set/OrthoFinder/Results_*/Orthologues/Orthologues_"$noext"_cluster/ $3/orthofinder/
		fi

		#MERGE DATA
		echo "Creating annotations output."
		python /usr/bin/merge_data.py $1 yes $3 $3 $4 $3/orthofinder/Orthologues_"$noext"_cluster/"$noext"_cluster__v__dromel_cluster.tsv $5

	fi

else #ELSE MEANS THESE ARE NOT NCBI PROTEIN IDS.

	echo "These are NOT NCBI protein IDs. Proceeding with KofamScan."
	#TAKES FASTA AND CREATES ACCESSION LIST. ACCESSION IS EVERYTHING BEFORE THE FIRST SPACE
	grep ">" $2 > $3/deflines.tmp
	sed -i 's/>//g' $3/deflines.tmp
	sed -i 's/\s.*$//' $3/deflines.tmp
	#ADD A TWO COLUMN FILE OF WITH AND WITHOUT VERSION THAT CAN BE USED TO MERGE LATER
	awk 'BEGIN {OFS="\t"} {print $1, $1}' $3/deflines.tmp > $3/ncbiversion.tmp
	awk 'BEGIN {OFS="\t"} { sub(/\.[0-9]+/, "", $2) }1' $3/ncbiversion.tmp > $3/ncbiver.tsv

	if grep -q $1 $3/kegg_org_codes.txt;
	then
		echo "This is a KEGG species".

		#PULL DATA
		echo "Pulling KEGG API data."
		bash /usr/bin/pull_data.sh $1 yes $3 non-ncbi $4

		#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
		if [[ -s "$3/link_ko_pathway.tsv" && -s "$3/list_pathway.tsv" && -s "$3/link_pathway_"$1".tsv" && -s "$3/list_pathway_"$1".tsv" && -s "$3/link_"$1"_ko.tsv" ]]
		then
    			echo "All KEGG files exist and are not empty."
		else
    			echo "One or more of the specified files are empty or do not exist."
			exit
		fi


		#RUN KOFAM HERE
		/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        sed -i 's/\..*$//' $3/ko_ncbi.tsv

		#IF FB RUN ORTHOFINDER AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$4" == FB ];
		then
			echo "Performing Flybase annotation".
			mkdir $3/orthofinder

			#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
			if [[ -s "$3/Fbgn_groupid.tsv" && -s "$3/Fbgn_CG.tsv" && -s "$3/Fbgn_fbpp.tsv" ]]
			then
				echo "All FlyBase files exist and are not empty."
			else
				echo "One or more of the specified files are empty or do not exist."
				exit
			fi

			ext="*.faa"
			if [[ $2 == $ext ]];
			then
				echo FAA
				noext=$(basename "$2" .faa)
			else
				ext="*.fasta"
				if [[ $2 == $ext ]];
				then
					echo FASTA
					noext=$(basename "$2" .fasta)
				else
					ext="*.fa"
					if [[ $2 == $ext ]];
					then
						echo FA
                                 		noext=$(basename "$2" .fa)
					else
						echo -e "This FASTA input file does not have an appropriate extension (.fa, .faa, .fasta)"
					fi
				fi
			fi
			echo -e "noext is: $noext"

			cp $2 $3/orthofinder/"$noext"_cluster.fa

			#RUN ORTHOFINDER WITH FASTA FROM INPUT SPECIES AND REFSET
			tar -xvzf /OF/ref_set.tgz -C $3/orthofinder
			orthofinder -f $3/orthofinder -t $cpus -b $3/orthofinder/ref_set

			#MOVE THE Orthologues_dromel_cluster DIR UP TO orthofinder
			mv $3/orthofinder/ref_set/OrthoFinder/Results_*/Orthologues/Orthologues_"$noext"_cluster/ $3/orthofinder/
		fi

		#MERGE DATA
		echo "Creating annotation outputs."
		python /usr/bin/merge_data.py $1 yes $3 $3 $4 $3/orthofinder/Orthologues_"$noext"_cluster/"$noext"_cluster__v__dromel_cluster.tsv $5

	else #ELSE MEANS THIS IS NOT A KEGG SPECIES

		echo "This is not a KEGG species".

		#PULL DATA
		echo "Pulling KEGG API data."
		bash /usr/bin/pull_data.sh $1 yes $3 non-ncbi $4

		#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
		if [[ -s "$3/link_ko_pathway.tsv" && -s "$3/list_pathway.tsv" ]]
		then
    			echo "All KEGG files exist and are not empty."
		else
    			echo "One or more of the specified files are empty or do not exist."
			exit
		fi

		#RUN KOFAM HERE
		/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        sed -i 's/\..*$//' $3/ko_ncbi.tsv

		#IF FB RUN ORTHOFINDER AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$4" == FB ];
		then
			echo "Performing Flybase annotation".
			mkdir $3/orthofinder

			#CHECK IF PULLED DATA FILES ARE PRESENT AND HAVE CONTENT BEFORE CONINUING
			if [[ -s "$3/Fbgn_groupid.tsv" && -s "$3/Fbgn_CG.tsv" && -s "$3/Fbgn_fbpp.tsv" ]]
			then
    				echo "All FlyBase files exist and are not empty."
			else
    				echo "One or more of the specified files are empty or do not exist."
				exit
			fi

			ext="*.faa"
			if [[ $2 == $ext ]];
			then
				echo FAA
				noext=$(basename "$2" .faa)
			else
				ext="*.fasta"
				if [[ $2 == $ext ]];
				then
					echo FASTA
					noext=$(basename "$2" .fasta)
				else
					ext="*.fa"
					if [[ $2 == $ext ]];
					then
						echo FA
                                 		noext=$(basename "$2" .fa)
					else
						echo -e "This FASTA input file does not have an appropriate extension (.fa, .faa, .fasta)"
					fi
				fi
			fi
			echo -e "noext is: $noext"

			cp $2 $3/orthofinder/"$noext"_cluster.fa

			#RUN ORTHOFINDER WITH FASTAS FROM INPUT SPECIES AND REFERENCE SET
			tar -xvzf /OF/ref_set.tgz -C $3/orthofinder
			orthofinder -f $3/orthofinder -t $cpus -b $3/orthofinder/ref_set

			#MOVE THE Orthologues_dromel_cluster DIR UP TO orthofinder
			mv $3/orthofinder/ref_set/OrthoFinder/Results_*/Orthologues/Orthologues_"$noext"_cluster/ $3/orthofinder/
		fi

		#MERGE DATA
		echo "Creating annotation outputs."
		python /usr/bin/merge_data.py $1 yes $3 $3 $4 $3/orthofinder/Orthologues_"$noext"_cluster/"$noext"_cluster__v__dromel_cluster.tsv $5
	fi
fi

if [ -f "$3"/link_ko_pathway.tsv ]; then rm "$3"/link_ko_pathway.tsv; fi
if [ -f "$3"/list_pathway.tsv ]; then rm "$3"/list_pathway.tsv; fi
if [ -f "$3"/conv_ncbi-proteinid_"$1".tsv ]; then rm "$3"/conv_ncbi-proteinid_"$1".tsv; fi
if [ -f "$3"/link_"$1"_ko.tsv ]; then rm "$3"/link_"$1"_ko.tsv; fi
if [ -f "$3"/link_pathway_"$1".tsv ]; then rm "$3"/link_pathway_"$1".tsv; fi
if [ -f "$3"/list_pathway_"$1".tsv ]; then rm "$3"/list_pathway_"$1".tsv; fi
if [ -f "$3"/deflines.tmp ]; then rm "$3"/deflines.tmp; fi
if [ -f "$3"/ko_ncbi.tsv ]; then rm "$3"/ko_ncbi.tsv; fi
if [ -f "$3"/Fbgn_CG.tsv ]; then rm "$3"/Fbgn_CG.tsv; fi
if [ -f "$3"/Fbgn_groupid.tsv ]; then rm "$3"/Fbgn_groupid.tsv; fi
if [ -f "$3"/pathway_group_data_latest.tsv ]; then rm $3/pathway_group_data_latest.tsv; fi
if [ -f "$3"/kofam_filtered_asterisk.txt ]; then rm "$3"/kofam_filtered_asterisk.txt; fi
if [ -f "$3"/kofam_result_full.txt ]; then rm "$3"/kofam_result_full.txt; fi
if [ -f "$3"/kegg_organisms.txt ]; then rm "$3"/kegg_organisms.txt; fi
if [ -f "$3"/kegg_org_codes.txt ]; then rm "$3"/kegg_org_codes.txt; fi
if [ -f "$3"/kegg_orgs_with_codes.txt ]; then rm "$3"/kegg_orgs_with_codes.txt; fi
if [ -n "$(ls $3/*pathway_group_data_fb* 2>/dev/null)" ]; then rm $3/*pathway_group_data_fb*; fi
if [ -n "$(ls $3/fbgn_annotation_ID_fb* 2>/dev/null)" ]; then rm $3/fbgn_annotation_ID_fb*; fi
if [ -n "$(ls $3/dmel-all-translation*.fasta* 2>/dev/null)" ]; then rm $3/dmel-all-translation*.fasta*; fi
if [ -n "$(ls $3/fbgn_fbtr_fbpp_fb* 2>/dev/null)" ]; then rm $3/fbgn_fbtr_fbpp_fb*; fi
if [ -f "$3"/Fbgn_fbpp.tsv ]; then rm "$3"/Fbgn_fbpp.tsv; fi
if [ -d "$3"/tmp ]; then rm -r "$3"/tmp; fi
if [ -f "$3"/tmp.txt ]; then rm  "$3"/tmp.txt; fi
if [ -d "$3"/orthofinder/ref_set ]; then rm -r "$3"/orthofinder/ref_set; fi
if [ -f "$3"/orthofinder/*_cluster.fa* ]; then rm -r "$3"/orthofinder/*_cluster.fa*; fi
if [ -d "$3"/orthofinder/ ]; then rm -r "$3"/orthofinder/; fi
if [ -f "$3"/ncbiversion.tmp ]; then rm "$3"/ncbiversion.tmp; fi
if [ -f "$3"/ncbiver.tsv ]; then rm "$3"/ncbiver.tsv; fi

endtime=$(date +%s)
seconds=$(($endtime - $starttime))
runtime=$(($seconds / 60))
echo "Run time: $runtime minutes"
