#! /bin/bash

#$1 is the KEGG species code (use NA or related species if species not in KEGG)
#$2 whether kofamscan is necessary ('yes' or 'no')
#$3 output directory
#$4 ncbi status of input FASTA accessions ('ncbi' or 'non-ncbi')
#$5 FB for flybase annotations, NA for none

kegg=$(grep $1 $3/kegg_org_codes.txt)

if [ ! -z "${kegg}" ] && [ $4 == "ncbi" ];
then
	#CAN'T PULL THIS FILE FOR SPECIES THAT AREN'T IN KEGG
	wget https://rest.kegg.jp/conv/ncbi-proteinid/"$1" -O $3/conv_ncbi-proteinid_"$1".tsv
	sed -i 's/ncbi-proteinid\://g' $3/conv_ncbi-proteinid_"$1".tsv
	sed -i "s/$1\://g" $3/conv_ncbi-proteinid_"$1".tsv
fi

if [ $2 == "yes" ];
then
	#THIS PULLS THE DATABASE FILES FOR KOFAMSCAN
	if [ ! -f /data/ko_list ] && [ ! -f /data/ko_list.gz ];
	then
		echo "Getting Kofam ko_list."
		wget https://www.genome.jp/ftp/db/kofam/ko_list.gz -O /data/ko_list.gz
		gunzip /data/ko_list.gz
	elif [ -f /data/ko_list.gz ];
	then
		echo "Decompressing ko_list."
		gunzip /data/ko_list.gz
	else
		echo "ko_list is already present."
	fi


	if [ ! -d /data/profiles ] && [ ! -f /data/profiles.tar.gz ];
	then
		echo "Getting Kofam profiles."
		wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz -O /data/profiles.tar.gz
		tar -xzf /data/profiles.tar.gz -C /data/
		rm /data/profiles.tar.gz
	elif [ -f /data/profiles.tar.gz ] && [ ! -d /data/profiles ];
	then
		echo "profiles.tar.gz present; decompressing."
		tar -xzf /data/profiles.tar.gz -C /data/
		rm /data/profiles.tar.gz
	else
		echo "profiles are already present."
	fi
fi



#THESE HAVE TO BE PULLED FOR ALL SPECIES
wget https://rest.kegg.jp/link/ko/pathway -O $3/link_ko_pathway.tsv
sed -i 's/ko\://g' $3/link_ko_pathway.tsv
sed -i 's/path\://g' $3/link_ko_pathway.tsv
grep -v ko $3/link_ko_pathway.tsv > $3/tmp.txt
mv $3/tmp.txt $3/link_ko_pathway.tsv

wget https://rest.kegg.jp/list/pathway -O $3/list_pathway.tsv

if [ "$1" != "NA" ];
then
	#THESE CAN BE PULLED FOR A RELATED KEGG SPECIES IF DESIRED
	wget https://rest.kegg.jp/link/pathway/"$1" -O  $3/link_pathway_"$1".tsv
	sed -i "s/$1\://g" $3/link_pathway_"$1".tsv
	sed -i 's/path\://g' $3/link_pathway_"$1".tsv

	wget https://rest.kegg.jp/list/pathway/"$1" -O $3/list_pathway_"$1".tsv

	wget https://rest.kegg.jp/link/"$1"/ko -O $3/link_"$1"_ko.tsv
	sed -i 's/ko\://g' $3/link_"$1"_ko.tsv
	sed -i "s/$1\://g" $3/link_"$1"_ko.tsv
fi

if [ "$5" == "FB" ];
then
	echo "Pulling Flybase data now."
	#THIS IS THE HARDCODED VERSION. IDEALLY WE COULD USE ftp:// AND A WILDCARD BUT THAT DOESN'T WORK WITH CURRENT FB SETUP (HTTPS)
	wget https://s3ftp.flybase.org/releases/FB2025_02/precomputed_files/genes/metabolic_pathway_group_data_fb_2025_02.tsv.gz -P $3/
	wget https://s3ftp.flybase.org/releases/FB2025_02/precomputed_files/genes/signaling_pathway_group_data_fb_2025_02.tsv.gz -P $3/
	gunzip -f $3/signaling_pathway_group_data_fb_2025_02.tsv.gz
	gunzip -f $3/metabolic_pathway_group_data_fb_2025_02.tsv.gz

	grep -h -v ^\# $3/signaling_pathway_group_data_fb_2025_02.tsv > $3/pathway_group_data_latest.tsv
	grep -h -v ^\# $3/metabolic_pathway_group_data_fb_2025_02.tsv >> $3/pathway_group_data_latest.tsv
	cut -f 1,3,6 $3/pathway_group_data_latest.tsv > $3/Fbgn_groupid.tsv

	wget https://s3ftp.flybase.org/releases/FB2025_02/precomputed_files/genes/fbgn_annotation_ID_fb_2025_02.tsv.gz -P $3
	gunzip -f $3/fbgn_annotation_ID_fb_2025_02.tsv.gz
	grep -h -v ^\# $3/fbgn_annotation_ID_fb_2025_02.tsv | cut -f 3,5 > $3/Fbgn_CG.tsv
	sed -i 's/Dmel_//g' $3/Fbgn_CG.tsv

	#PULL FBGN TO FBPP FILES
	wget https://s3ftp.flybase.org/releases/FB2025_02/precomputed_files/genes/fbgn_fbtr_fbpp_fb_2025_02.tsv.gz -P $3/
	gunzip -f $3/fbgn_fbtr_fbpp_fb_2025_02.tsv.gz
	grep -v ^\# $3/fbgn_fbtr_fbpp_fb_2025_02.tsv | cut -f 1,3 > $3/Fbgn_fbpp.tsv


	#THIS SYSTEM TRIES TO PULL THE FILES BY GUESSING THE CURRENT VERSION INSTEAD OF USING A WILDCARD.
	#ONE FILE WORKS ONE TIME...THEN IT GETS BLOCKED BY THE SERVER I THINK
	#styr=2025
	#endyr=2025

	#for year in $(seq $styr $endyr)
	#do
    		#for month in $(seq -w 1 12)
		#do
		        #echo "Processing: $year-$month"
			#THIS PULLS THE FLYBASE ANNOTATIONS
			#wget -P "$3" https://s3ftp.flybase.org/releases/current/precomputed_files/genes/metabolic_pathway_group_data_fb_"$year"_"$month".tsv.gz
			#wget -P "$3" https://s3ftp.flybase.org/releases/current/precomputed_files/genes/signaling_pathway_group_data_fb_"$year"_"$month".tsv.gz
                        #wget -P "$3" https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_fb_"$year"_"$month".tsv.gz
                        #wget -P "$3" https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_fbtr_fbpp_fb_"$year"_"$month".tsv.gz
    		#done
	#done
	#THIS REMOVES THE EMPTY FILES FROM THE GENERIC DOWNLOADS AND PROCESSES THE REMAINING FILES
	#find "$3"/ -size 0 -delete
	#gunzip -f "$3"/signaling_pathway_group_data_fb_"$year"_"$month".tsv.gz
	#gunzip -f "$3"/metabolic_pathway_group_data_fb_"$year"_"$month".tsv.gz
	#grep -h -v ^\# $3/signaling_pathway_group_data_fb_"$year"_"$month".tsv > $3/pathway_group_data_latest.tsv
	#grep -h -v ^\# $3/metabolic_pathway_group_data_fb_"$year"_"$month".tsv >> $3/pathway_group_data_latest.tsv
	#cut -f 1,3,6 $3/pathway_group_data_latest.tsv > $3/Fbgn_groupid.tsv

	#gunzip -f $3/fbgn_annotation_ID_fb_"$year"_"$month".tsv.gz
	#grep -h -v ^\# $3/fbgn_annotation_ID_fb_"$year"_"$month".tsv | cut -f 3,5 > $3/Fbgn_CG.tsv
	#sed -i 's/Dmel_//g' $3/Fbgn_CG.tsv

	#gunzip -f $3/fbgn_fbtr_fbpp_fb_"$year"_"$month".tsv.gz
	#grep -v ^\# $3/fbgn_fbtr_fbpp_fb_"$year"_"$month".tsv | cut -f 1,3 > $3/Fbgn_fbpp.tsv
fi


