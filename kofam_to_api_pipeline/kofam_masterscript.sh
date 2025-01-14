#! /bin/bash

if [ $1 == "help" ];
then
	echo "Help and Usage:
	There are 4 positional arguments.
	1: KEGG species code (NA or related species code if species not in KEGG; 'help' to see this help and usage statement)
	2: input file (protein FASTA without header lines)
	3: output directory (must be an existing directory)
	4: 'FB' for flybase annotations, 'NA' for none"
	exit 0
fi

# WORKS-TESTS WHETHER ACCESSIONS ARE NCBI PROTEIN IDS
acc1=$(head -n1 $2 | sed 's/>//g' | sed 's/\s.*$//')
if  [[ $acc1 == NP_* ]] || [[ $acc1 == XP_* ]] ;
then
	ncbi=true
	echo "$acc1 These are NCBI protein IDs."
else
	ncbi=false
fi


if [ "$ncbi" == true ] ;
then
	# WORKS-TAKES FASTA AND CREATES ACCESSION LIST. ACCESSION IS EVERYTHING BEFORE THE FIRST SPACE
	grep ">" $2 > $3/deflines.tmp
	sed -i 's/>//g' $3/deflines.tmp
	sed -i 's/\s.*$//' $3/deflines.tmp
	sed -i 's/.[0-9]$//' $3/deflines.tmp
	readarray -t defarray < $3/deflines.tmp

	if grep -q $1 /usr/bin/kegg_org_codes.txt; #IF THIS IS A KEGG SPECIES
	then
		#PULL DATA
		echo "This is a KEGG species code. Pulling KEGG API data now."
		bash /usr/bin/pull_data.sh $1 yes $3 ncbi $4

		#NEED TO COMPARE DEFLINES.TMP TO SPECIFIED SPECIES CODE AND DECIDE IF THEY ARE THE SAME SPECIES
		#if i were to pull data then grep the first defline.tmp entry from ncbi_ko to see if it exists...
		echo "${defarray[0]}"
		if grep -q "${defarray[0]}" $3/conv_ncbi-proteinid_"$1".tsv; #TESTING IF INPUT IDS ARE THE SAME SPECIES AS THE KEGG CODE
		then
			#IF YES, THEN PROCEED WITH FB, MERGE FROM API DATA
			echo "IDs are $1 species IDs"

			#IF FB AND NOT DME RUN HMMER AND PROCEED TO MERGE (INCLUDING FLYBASE)
			if [ "$1" != "dme" ] && [ "$4" == "FB" ];
			then
				echo "Performing Flybase annotation".
				avail=$(nproc)
				cpus=$(( $avail - 1 ))
				phmmer --cpu $cpus --tblout $3/FB_phmmer.txt -o /dev/null -E 0.05 $2 $3/dmel-all-translation-*.fasta
 				#PULL MATCHES FROM OUTPUT
				grep -v ^\# $3/FB_phmmer.txt | awk -F " +" '{print $3}' | sort | uniq > $3/phmmacc.txt
				readarray -t phmmarray < $3/phmmacc.txt
				for each in "${phmmarray[@]}"
        			do
					grep -m 1 $each $3/FB_phmmer.txt > $3/phmm_tophits.txt
					awk -F ' +' '{ OFS="\t"; print $1, $3 }'  $3/phmm_tophits.txt >> $3/phmm_matches.txt
				done
			fi

			#MERGE DATA HERE
			echo "Creating annotations output."
			python /usr/bin/merge_data.py $1 no $3 $3 $4

		else
			#IF NO, THEN RUN KOFAM, FILTER, FB, MERGE FROM KOFAM DATA
			echo "IDs are NOT $1 species IDs"

			#PULL ADDITIONAL DATA FOR KOFAMSCAN
			echo "Pulling more KEGG API data now."
			bash /usr/bin/pull_data.sh $1 yes $3 ncbi $4

			#RUN KOFAMSCAN
			echo "This is not a KEGG species code. Running KofamScan now."
			avail=$(nproc)
			cpus=$(( $avail - 1 ))
			/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

			#FILTER KOFAM HERE
			echo "Filtering KofamScan results"
			grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        	awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        	sed -i 's/.[0-9]$//' $3/ko_ncbi.tsv

			#IF FB AND NOT DME RUN HMMER AND PROCEED TO MERGE (INCLUDING FLYBASE)
			if [ "$1" != "dme" ] && [ "$4" == "FB" ];
			then
				echo "Performing Flybase annotation".
				avail=$(nproc)
				cpus=$(( $avail - 1 ))
				phmmer --cpu $cpus --tblout $3/FB_phmmer.txt -o /dev/null -E 0.05 $2 $3/dmel-all-translation-*.fasta
 				#PULL MATCHES FROM OUTPUT
				grep -v ^\# $3/FB_phmmer.txt | awk -F " +" '{print $3}' | sort | uniq > $3/phmmacc.txt
				readarray -t phmmarray < $3/phmmacc.txt
				for each in "${phmmarray[@]}"
        			do
					grep -m 1 $each $3/FB_phmmer.txt > $3/phmm_tophits.txt
					awk -F ' +' '{ OFS="\t"; print $1, $3 }'  $3/phmm_tophits.txt >> $3/phmm_matches.txt
				done
			fi

			#MERGE DATA HERE
			echo "Creating annotations output."
			python /usr/bin/merge_data.py $1 yes $3 $3 $4


		fi

	else # ELSE MEANS THE THE CODE IS NOT A KEGG SPECIES CODE

		#PULL DATA
		echo "Pulling KEGG API data."
		bash /usr/bin/pull_data.sh $1 yes $3 ncbi $4

		#RUN KOFAMSCAN
		echo "This is not a KEGG species code. Running KofamScan now."
		avail=$(nproc)
		cpus=$(( $avail - 1 ))
		/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        sed -i 's/.[0-9]$//' $3/ko_ncbi.tsv

		#IF FB AND NOT DME RUN HMMER AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$1" != "dme" ] && [ "$4" == "FB" ];
		then
			echo "Performing Flybase annotation".
			avail=$(nproc)
			cpus=$(( $avail - 1 ))
			phmmer --cpu $cpus --tblout $3/FB_phmmer.txt -o /dev/null -E 0.05 $2 $3/dmel-all-translation-*.fasta
 			#PULL MATCHES FROM OUTPUT
			grep -v ^\# $3/FB_phmmer.txt | awk -F " +" '{print $3}' | sort | uniq > $3/phmmacc.txt
			readarray -t phmmarray < $3/phmmacc.txt
			for each in "${phmmarray[@]}"
        		do
				grep -m 1 $each $3/FB_phmmer.txt > $3/phmm_tophits.txt
				awk -F ' +' '{ OFS="\t"; print $1, $3 }'  $3/phmm_tophits.txt >> $3/phmm_matches.txt
			done
		fi

		#MERGE DATA
		echo "Creating annotations output."
		python /usr/bin/merge_data.py $1 yes $3 $3 $4
	fi

else #ELSE MEANS THESE ARE NOT NCBI PROTEIN IDS.

	echo "These are NOT NCBI protein IDs. Proceeding with KofamScan."

	if grep -q $1 /usr/bin/kegg_org_codes.txt;
	then
		echo "This is a KEGG species".

		#PULL DATA
		echo "Pulling KEGG API data."
		bash /usr/bin/pull_data.sh $1 yes $3 non-ncbi $4

		#RUN KOFAM HERE
		avail=$(nproc)
		cpus=$(( $avail - 1 ))
		#NEED TO MAKE THIS WORK WITH HMM FILES INSTEAD OF EUKARYOTE.HAL??
		/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        sed -i 's/.[0-9]$//' $3/ko_ncbi.tsv

		#IF FB RUN HMMER AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$4" == FB ];
		then
			echo "Performing Flybase annotation".
			phmmer --cpu $cpus --tblout $3/FB_phmmer.txt -o /dev/null -E 0.05 $2 $3/dmel-all-translation-*.fasta
 			#PULL MATCHES FROM OUTPUT
			grep -v ^\# $3/FB_phmmer.txt | awk -F " +" '{print $3}' | sort | uniq > $3/phmmacc.txt
			readarray -t phmmarray < $3/phmmacc.txt
			for each in "${phmmarray[@]}"
        		do
				grep -m 1 $each $3/FB_phmmer.txt > $3/phmm_tophits.txt
				awk -F ' +' '{ OFS="\t"; print $1, $3 }'  $3/phmm_tophits.txt >> $3/phmm_matches.txt
			done
		fi

		#MERGE DATA
		echo "Creating annotation outputs."
		python /usr/bin/merge_data.py $1 yes $3 $3 $4

	else #ELSE MEANS THIS IS NOT A KEGG SPECIES

		echo "This is not a KEGG species".

		#PULL DATA
		echo "Pulling KEGG API data."
		bash /usr/bin/pull_data.sh $1 yes $3 non-ncbi $4

		#RUN KOFAM HERE
		avail=$(nproc)
		cpus=$(( $avail - 1 ))
		#NEED TO MAKE THIS WORK WITH HMM FILES INSTEAD OF EUKARYOTE.HAL??
		/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        sed -i 's/.[0-9]$//' $3/ko_ncbi.tsv

		#IF FB RUN HMMER AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$4" == FB ];
		then
			echo "Performing Flybase annotation".
			phmmer --cpu $cpus --tblout $3/FB_phmmer.txt -o /dev/null -E 0.05 $2 $3/dmel-all-translation-*.fasta
 			#PULL MATCHES FROM OUTPUT
			grep -h -v ^\# $3/FB_phmmer.txt | awk -F " +" '{print $3}' | sort | uniq > $3/phmmacc.txt
			readarray -t phmmarray < $3/phmmacc.txt
			for each in "${phmmarray[@]}"
        		do
				grep -h -m 1 $each $3/FB_phmmer.txt > $3/phmm_tophits.txt
				awk -F ' +' '{ OFS="\t"; print $1, $3 }'  $3/phmm_tophits.txt >> $3/phmm_matches.txt
			done
		fi

		#MERGE DATA
		echo "Creating annotation outputs."
		python /usr/bin/merge_data.py $1 yes $3 $3 $4
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
if [ -n "$(ls $3/*pathway_group_data_fb* 2>/dev/null)" ]; then rm $3/*pathway_group_data_fb*; fi
if [ -n "$(ls $3/fbgn_annotation_ID_fb* 2>/dev/null)" ]; then rm $3/fbgn_annotation_ID_fb*; fi
if [ -n "$(ls $3/dmel-all-translation*.fasta* 2>/dev/null)" ]; then rm $3/dmel-all-translation*.fasta*; fi
if [ -n "$(ls $3/fbgn_fbtr_fbpp_fb* 2>/dev/null)" ]; then rm $3/fbgn_fbtr_fbpp_fb*; fi
if [ -f "$3"/phmmacc.txt ]; then  rm "$3"/phmmacc.txt; fi
if [ -f "$3"/phmm_tophits.txt ]; then  rm "$3"/phmm_tophits.txt; fi
if [ -f "$3"/phmm_matches.txt ]; then  rm "$3"/phmm_matches.txt; fi
if [ -f "$3"/Fbgn_fbpp.tsv ]; then rm "$3"/Fbgn_fbpp.tsv; fi
if [ -f "$3"/FB_phmmer.txt ]; then rm "$3"/FB_phmmer.txt; fi
if [ -d "$3"/tmp ]; then rm -r "$3"/tmp; fi
if [ -d "$3"/tmp.txt ]; then rm  "$3"/tmp.txt; fi
