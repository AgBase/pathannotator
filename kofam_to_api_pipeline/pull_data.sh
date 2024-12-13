#! /bin/bash

#MUST SUPPLY $1 OR $2 OR BOTH
#$1 is the KEGG species code (use NA or related species if species not in KEGG)
#$2 is the filtered * input file ( from kofam, filter kofam)
#$3 output directory
#$4 ncbi status of input FASTA accessions ('ncbi' or 'non-ncbi')


kegg=$(grep $1 /usr/bin/kegg_org_codes.txt)

if [ ! -z "${kegg}" ] && [ $4 == "ncbi" ];
then
	#CAN'T PULL THIS FILE FOR SPECIES THAT AREN'T IN KEGG
	wget https://rest.kegg.jp/conv/ncbi-proteinid/"$1" -O $3/conv_ncbi-proteinid_"$1".tsv
	sed -i 's/ncbi-proteinid\://g' $3/conv_ncbi-proteinid_"$1".tsv
	sed -i "s/$1\://g" $3/conv_ncbi-proteinid_"$1".tsv
else

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
		tar -xzf /data/profiles.tar.gz
		rm /data/profiles.tar.gz
	elif [ -f /data/profiles.tar.gz ] && [ ! -f /data/profiles ];
	then
		echo "profiles.tar.gz present; decompressing."
		tar -xzf /data/profiles.tar.gz
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

if [ "$1" == "dme" ];
then
	echo "Species code is 'dme'; pulling Flybase data now."
	#THIS PULLS THE FLYBASE ANNOTATIONS
	wget -r -nd -A gz --accept-regex "signaling_pathway_group_data*" -P $3/ 'https://ftp.flybase.net/releases/current/precomputed_files/genes/'
	wget -r -nd -A gz --accept-regex "metabolic_pathway_group_data*" -P $3/ 'https://ftp.flybase.net/releases/current/precomputed_files/genes/'
	gunzip $3/signaling_pathway_group_data_*
	gunzip $3/metabolic_pathway_group_data_*
	grep -v ^\# $3/signaling_pathway_group_data_* > $3/pathway_group_data_latest.tsv
	grep -v ^\# $3/metabolic_pathway_group_data_* >> $3/pathway_group_data_latest.tsv
	cut -f 1,3,6 $3/pathway_group_data_latest.tsv > $3/Fbgn_groupid.tsv

	wget -r -nd -A gz --accept-regex "fbgn_annotation_ID_fb*" -P $3/ 'https://ftp.flybase.org/releases/current/precomputed_files/genes/'
	gunzip $3/fbgn_annotation_ID_fb*
	grep -v ^\# $3/fbgn_annotation_ID_fb* | cut -f 3,5 > $3/Fbgn_CG.tsv
	sed -i 's/Dmel_//g' $3/Fbgn_CG.tsv

	#PULL FB PROTEIN FASTA
	wget -r -nd -A gz --accept-regex "dmel-all-translation*" -P $3/ 'https://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/'
        gunzip $3/dmel-all-translation-*

	#PULL FBGN TO FBPP FILES
	wget -r -nd -A gz  --accept-regex "fbgn_fbtr_fbpp_fb*" -P $3/ 'https://ftp.flybase.net/releases/current/precomputed_files/genes/'
	gunzip $3/fbgn_fbtr_fbpp_fb_*
	grep -v ^\# $3/fbgn_fbtr_fbpp_fb_* | cut -f 1,3 > $3/Fbgn_fbpp.tsv
fi
