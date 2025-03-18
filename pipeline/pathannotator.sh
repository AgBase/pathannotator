#! /bin/bash

#CHECK FOR OUTDIR. IF IT DOESN'T EXIST CREATE IT
if [ ! -d "$3" ]; then mkdir "$3"; fi


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
	5: GFF corresponding to input FASTA
	6: genomics FASTA for input species

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

# WORKS-TESTS WHETHER ACCESSIONS ARE NCBI PROTEIN IDS
acc1=$(head -n1 $2 | sed 's/>//g' | sed 's/\s.*$//')
if  [[ $acc1 == NP_* ]] || [[ $acc1 == XP_* ]] ;
then
	ncbi=true
	echo "$acc1 These are NCBI protein IDs."
else
	ncbi=false
fi


#PULLS THE KEGG ORG CODES FILE (NEEDS TO BE IN HERE, NOT PULL_DATA.SH BECAUSE PULL DATA ONLY RUNS IN THE IF STATEMENTS BELOW
wget https://rest.kegg.jp/list/genome -O $3/kegg_organisms.txt
grep ';' $3/kegg_organisms.txt > $3/kegg_orgs_with_codes.txt
cut -f 2 $3/kegg_orgs_with_codes.txt > $3/kegg_org_codes.txt
sed -i 's/;.*$//g' $3/kegg_org_codes.txt

if [ "$ncbi" == true ] ;
then
	# WORKS-TAKES FASTA AND CREATES ACCESSION LIST. ACCESSION IS EVERYTHING BEFORE THE FIRST SPACE
	grep ">" $2 > $3/deflines.tmp
	sed -i 's/>//g' $3/deflines.tmp
	sed -i 's/\s.*$//' $3/deflines.tmp
	sed -i 's/.[0-9]$//' $3/deflines.tmp
	readarray -t defarray < $3/deflines.tmp

	if grep -q $1 $3/kegg_org_codes.txt; #IF THIS IS A KEGG SPECIES
	then
		#PULL DATA
		echo "This is a KEGG species code. Pulling KEGG API data now."
		bash /usr/bin/pull_data.sh $1 no $3 ncbi $4

		#NEED TO COMPARE DEFLINES.TMP TO SPECIFIED SPECIES CODE AND DECIDE IF THEY ARE THE SAME SPECIES
		echo "${defarray[0]}"
		if grep -q "${defarray[0]}" $3/conv_ncbi-proteinid_"$1".tsv; #TESTING IF INPUT IDS ARE THE SAME SPECIES AS THE KEGG CODE
		then
			#IF YES, MERGE FROM API DATA
			echo "IDs are $1 species IDs"

			#IF FB AND NOT DME RUN DIAMOND AND PROCEED TO MERGE (INCLUDING FLYBASE)
			if [ "$1" != "dme" ] && [ "$4" == "FB" ];
			then
				echo "Performing Flybase annotation".
#RUN agat ON FB AND REPRESENTATIVE SPECIES GFFS--SAVE NEW SINGLE-TRANSCRIPT FASTA
#FLYBASE GFF FILE IS FUCKED UP. NEED TO REMOVE EVERYTHING THAT ISN'T "\tFlyBase\t".
#AGAT ONLY PULLS WHICHEVER FEATURES ARE PARENTS OF THE CDS. GFFREAD ALWAYS PULL TRANSCRIPT IDS EVEN WHEN IT IS AMINO ACID SEQUENCE.
#				tail -n +2 $3/dmel-all*.gff | grep -P "\tFlyBase\t" > $3/dmel.gff.tmp
#				mv $3/dmel.gff.tmp $3/dmel.gff
#				perl /opt/conda/bin/agat_sp_keep_longest_isoform.pl -gff $3/dmel.gff   -o $3/dromel_longest_isoform.gff -c /usr/bin/agat_config.yaml
#				gffread -y $3/dromel_longest_isoform.fa -g $3/dmel-all-chromosome*.fasta $3/dromel_longest_isoform.gff

#				perl /opt/conda/bin/agat_sp_keep_longest_isoform.pl -gff $3/GCF_031307605.1_icTriCast1.1_genomic.gff   -o $3/tricas_longest_isoform.gff -c /usr/bin/agat_config.yaml
#				gffread -y $3/tricas_longest_isoform.fa -g $3/GCF_031307605.1_icTriCast1.1_genomic.fna $3/tricas_longest_isoform.gff

#				perl /opt/conda/bin/agat_sp_keep_longest_isoform.pl -gff $3/GCF_003254395.2_Amel_HAv3.1_genomic.gff   -o $3/apimel_longest_isoform.gff -c /usr/bin/agat_config.yaml
#				gffread -y $3/apimel_longest_isoform.fa -g $3/GCF_003254395.2_Amel_HAv3.1_genomic.fna $3/apimel_longest_isoform.gff

#				perl /opt/conda/bin/agat_sp_keep_longest_isoform.pl -gff $3/GCF_014839805.1_JHU_Msex_v1.0_genomic.gff   -o $3/mansex_longest_isoform.gff -c /usr/bin/agat_config.yaml
#				gffread -y $3/mansex_longest_isoform.fa -g $3/GCF_014839805.1_JHU_Msex_v1.0_genomic.fna $3/mansex_longest_isoform.gff

#				perl /opt/conda/bin/agat_sp_keep_longest_isoform.pl -gff $3/GCF_020184175.1_ASM2018417v2_genomic.gff   -o $3/aphgos_longest_isoform.gff -c /usr/bin/agat_config.yaml
#				gffread -y $3/aphgos_longest_isoform.fa -g $3/GCF_020184175.1_ASM2018417v2_genomic.fna $3/aphgos_longest_isoform.gff

#				perl /opt/conda/bin/agat_sp_keep_longest_isoform.pl -gff $3/GCF_023897955.1_iqSchGreg1.2_genomic.gff   -o $3/schgre_longest_isoform.gff -c /usr/bin/agat_config.yaml
#				gffread -y $3/schgre_longest_isoform.fa -g $3/GCF_023897955.1_iqSchGreg1.2_genomic.fna $3/schgre_longest_isoform.gff

#				sed -i 's/\./\-/g' $3/*longest_isoform.fa

#RUN agat ON INPUT GFF--SAVE NEW SINGLE-TRANSCRIPT FASTA
#AGAT ONLY PULLS WHICHEVER FEATURES ARE PARENTS OF THE CDS. GFFREAD ALWAYS PULL TRANSCRIPT IDS EVEN WHEN IT IS AMINO ACID SEQUENCE.

#				noext=$(basename "$5" .gff)
#				perl /opt/conda/bin/agat_sp_keep_longest_isoform.pl -gff $5   -o $3/"$noext"_longest_isoform.gff -c /usr/bin/agat_config.yaml
#				gffread -y $3/"$noext"_longest_isoform.fa -g $6  $3/"$noext"_longest_isoform.gff
#				sed -i 's/\./\-/g' $3/*longest_isoform.fa

#RUN ORTHOFINDER WITH SINGLE-TRANCRIPT FASTAS FROM INPUT SPECIES AND DROMEL
				mkdir $3/orthofinder
				mv $3/*longest_isoform.fa $3/orthofinder
				orthofinder -f $3/orthofinder -t 12
#PULL MATCHES FROM OUTPUT
#				sed -i '1i Query_ID\tQuery_length\tQuery_start\tQuery_end\tSubject_ID\tSubject_length\tSubject_start\tSubject_end\tE_value\tPercent_ID\tPercent_positive_ID\tGap_openings\tTotal_gaps\tBitscore\tRaw_score\tAlignment_length' $3/diamond_out.tsv
#				awk '{ if(($10 > 70) && ($16/$2 > 0.7) && ($12 < 9) && ($2/$6 <= 1.2)) { print }}' $3/diamond_out.tsv > $3/dia_matches.tsv
#				cut -f 1,5 $3/dia_matches.tsv > $3/FB_diamond.tsv
			fi

			#MERGE DATA HERE
#			echo "Creating annotations output."
#			python /usr/bin/merge_data.py $1 no $3 $3 $4
		else
			#IF NO, THEN RUN KOFAM, FILTER, FB, MERGE FROM KOFAM DATA
			echo "IDs are NOT $1 species IDs"

			#PULL ADDITIONAL DATA FOR KOFAMSCAN
			echo "Pulling more KEGG API data now."
			bash /usr/bin/pull_data.sh $1 yes $3 ncbi $4

			#RUN KOFAMSCAN
			echo "Running KofamScan now."
			avail=$(getconf _NPROCESSORS_ONLN)
			cpus=$(( $avail - 1 ))

			/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

			#FILTER KOFAM HERE
			echo "Filtering KofamScan results"
			grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        	awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        	sed -i 's/.[0-9]$//' $3/ko_ncbi.tsv

			#IF FB AND NOT DME RUN DIAMOND AND PROCEED TO MERGE (INCLUDING FLYBASE)
			if [ "$1" != "dme" ] && [ "$4" == "FB" ];
			then
				echo "Performing Flybase annotation".
				diamond version
				diamond makedb --in $3/dmel-all-translation-*.fasta --db $3/dmel_db
				diamond blastp -q $2 -d $3/dmel_db -o $3/diamond_out.tsv --max-target-seqs 3 --outfmt 6 qseqid qlen qstart qend sseqid slen sstart send evalue pident ppos gapopen gaps bitscore score length
				sed -i '1i Query_ID\tQuery_length\tQuery_start\tQuery_end\tSubject_ID\tSubject_length\tSubject_start\tSubject_end\tE_value\tPercent_ID\tPercent_positive_ID\tGap_openings\tTotal_gaps\tBitscore\tRaw_score\tAlignment_length' $3/diamond_out.tsv

 				#PULL MATCHES FROM OUTPUT
				awk '{ if(($10 > 70) && ($16/$2 > 0.7) && ($12 < 9) && ($2/$6 <= 1.2)) { print }}' $3/diamond_out.tsv > $3/dia_matches.tsv
				cut -f 1,5 $3/dia_matches.tsv > $3/FB_diamond.tsv
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
		avail=$(getconf _NPROCESSORS_ONLN)
		cpus=$(( $avail - 1 ))

		/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        sed -i 's/.[0-9]$//' $3/ko_ncbi.tsv

		#IF FB AND NOT DME RUN DIAMOND AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$1" != "dme" ] && [ "$4" == "FB" ];
		then
			echo "Performing Flybase annotation".
			diamond version
			diamond makedb --in $3/dmel-all-translation-*.fasta --db $3/dmel_db
			diamond blastp -q $2 -d $3/dmel_db -o $3/diamond_out.tsv --max-target-seqs 3 --outfmt 6 qseqid qlen qstart qend sseqid slen sstart send evalue pident ppos gapopen gaps bitscore score length
			sed -i '1i Query_ID\tQuery_length\tQuery_start\tQuery_end\tSubject_ID\tSubject_length\tSubject_start\tSubject_end\tE_value\tPercent_ID\tPercent_positive_ID\tGap_openings\tTotal_gaps\tBitscore\tRaw_score\tAlignment_length' $3/diamond_out.tsv

 			#PULL MATCHES FROM OUTPUT
			awk '{ if(($10 > 70) && ($16/$2 > 0.7) && ($12 < 9) && ($2/$6 <= 1.2)) { print }}' $3/diamond_out.tsv > $3/dia_matches.tsv
			cut -f 1,5 $3/dia_matches.tsv > $3/FB_diamond.tsv
		fi

		#MERGE DATA
		echo "Creating annotations output."
		python /usr/bin/merge_data.py $1 yes $3 $3 $4
	fi

else #ELSE MEANS THESE ARE NOT NCBI PROTEIN IDS.

	echo "These are NOT NCBI protein IDs. Proceeding with KofamScan."

	if grep -q $1 $3/kegg_org_codes.txt;
	then
		echo "This is a KEGG species".

		#PULL DATA
		echo "Pulling KEGG API data."
		bash /usr/bin/pull_data.sh $1 yes $3 non-ncbi $4

		#RUN KOFAM HERE
		avail=$(getconf _NPROCESSORS_ONLN)
		cpus=$(( $avail - 1 ))

		/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        sed -i 's/.[0-9]$//' $3/ko_ncbi.tsv

		#IF FB RUN DIAMOND AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$4" == FB ];
		then
			echo "Performing Flybase annotation".
			diamond version
			diamond makedb --in $3/dmel-all-translation-*.fasta --db $3/dmel_db
			diamond blastp -q $2 -d $3/dmel_db -o $3/diamond_out.tsv --max-target-seqs 3 --outfmt 6 qseqid qlen qstart qend sseqid slen sstart send evalue pident ppos gapopen gaps bitscore score length
			sed -i '1i Query_ID\tQuery_length\tQuery_start\tQuery_end\tSubject_ID\tSubject_length\tSubject_start\tSubject_end\tE_value\tPercent_ID\tPercent_positive_ID\tGap_openings\tTotal_gaps\tBitscore\tRaw_score\tAlignment_length' $3/diamond_out.tsv

 			#PULL MATCHES FROM OUTPUT
			awk '{ if(($10 > 70) && ($16/$2 > 0.7) && ($12 < 9) && ($2/$6 <= 1.2)) { print }}' $3/diamond_out.tsv > $3/dia_matches.tsv
			cut -f 1,5 $3/dia_matches.tsv > $3/FB_diamond.tsv
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
		avail=$(getconf _NPROCESSORS_ONLN)
		cpus=$(( $avail - 1 ))

		/usr/bin/kofam_scan/exec_annotation -o $3/kofam_result_full.txt -f detail --tmp-dir $3/tmp --cpu $cpus -k /data/ko_list -p /data/profiles/eukaryote.hal $2

		#FILTER KOFAM HERE
		echo "Filtering KofamScan results"
		grep -P "^\*" $3/kofam_result_full.txt >> $3/kofam_filtered_asterisk.txt
	        awk '{ print $3"\t"$2 }' $3/kofam_filtered_asterisk.txt > $3/ko_ncbi.tsv
	        sed -i 's/.[0-9]$//' $3/ko_ncbi.tsv

		#IF FB RUN DIAMOND AND PROCEED TO MERGE (INCLUDING FLYBASE)
		if [ "$4" == FB ];
		then
			echo "Performing Flybase annotation".
			diamond version
			diamond makedb --in $3/dmel-all-translation-*.fasta --db $3/dmel_db
			diamond blastp -q $2 -d $3/dmel_db -o $3/diamond_out.tsv --max-target-seqs 3 --outfmt 6 qseqid qlen qstart qend sseqid slen sstart send evalue pident ppos gapopen gaps bitscore score length
			sed -i '1i Query_ID\tQuery_length\tQuery_start\tQuery_end\tSubject_ID\tSubject_length\tSubject_start\tSubject_end\tE_value\tPercent_ID\tPercent_positive_ID\tGap_openings\tTotal_gaps\tBitscore\tRaw_score\tAlignment_length' $3/diamond_out.tsv

 			#PULL MATCHES FROM OUTPUT
			awk '{ if(($10 > 70) && ($16/$2 > 0.7) && ($12 < 9) && ($2/$6 <= 1.2)) { print }}' $3/diamond_out.tsv > $3/dia_matches.tsv
			cut -f 1,5 $3/dia_matches.tsv > $3/FB_diamond.tsv
		fi

		#MERGE DATA
		echo "Creating annotation outputs."
		python /usr/bin/merge_data.py $1 yes $3 $3 $4
	fi
fi

#if [ -f "$3"/link_ko_pathway.tsv ]; then rm "$3"/link_ko_pathway.tsv; fi
#if [ -f "$3"/list_pathway.tsv ]; then rm "$3"/list_pathway.tsv; fi
#if [ -f "$3"/conv_ncbi-proteinid_"$1".tsv ]; then rm "$3"/conv_ncbi-proteinid_"$1".tsv; fi
#if [ -f "$3"/link_"$1"_ko.tsv ]; then rm "$3"/link_"$1"_ko.tsv; fi
#if [ -f "$3"/link_pathway_"$1".tsv ]; then rm "$3"/link_pathway_"$1".tsv; fi
#if [ -f "$3"/list_pathway_"$1".tsv ]; then rm "$3"/list_pathway_"$1".tsv; fi
#if [ -f "$3"/deflines.tmp ]; then rm "$3"/deflines.tmp; fi
#if [ -f "$3"/ko_ncbi.tsv ]; then rm "$3"/ko_ncbi.tsv; fi
#if [ -f "$3"/Fbgn_CG.tsv ]; then rm "$3"/Fbgn_CG.tsv; fi
#if [ -f "$3"/Fbgn_groupid.tsv ]; then rm "$3"/Fbgn_groupid.tsv; fi
#if [ -f "$3"/pathway_group_data_latest.tsv ]; then rm $3/pathway_group_data_latest.tsv; fi
#if [ -f "$3"/kofam_filtered_asterisk.txt ]; then rm "$3"/kofam_filtered_asterisk.txt; fi
#if [ -f "$3"/kegg_organisms.txt ]; then rm "$3"/kegg_organisms.txt; fi
#if [ -f "$3"/kegg_org_codes.txt ]; then rm "$3"/kegg_org_codes.txt; fi
#if [ -f "$3"/kegg_orgs_with_codes.txt ]; then rm "$3"/kegg_orgs_with_codes.txt; fi
#if [ -n "$(ls $3/*pathway_group_data_fb* 2>/dev/null)" ]; then rm $3/*pathway_group_data_fb*; fi
#if [ -n "$(ls $3/fbgn_annotation_ID_fb* 2>/dev/null)" ]; then rm $3/fbgn_annotation_ID_fb*; fi
#if [ -n "$(ls $3/dmel-all-translation*.fasta* 2>/dev/null)" ]; then rm $3/dmel-all-translation*.fasta*; fi
#if [ -n "$(ls $3/fbgn_fbtr_fbpp_fb* 2>/dev/null)" ]; then rm $3/fbgn_fbtr_fbpp_fb*; fi
#if [ -f "$3"/Fbgn_fbpp.tsv ]; then rm "$3"/Fbgn_fbpp.tsv; fi
#if [ -d "$3"/tmp ]; then rm -r "$3"/tmp; fi
#if [ -f "$3"/tmp.txt ]; then rm  "$3"/tmp.txt; fi

endtime=$(date +%s)
seconds=$(($endtime - $starttime))
runtime=$(($seconds / 60))
echo "Run time: $runtime minutes"
