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
#	#INSTEAD FOR NON-KEGG SPECIES USE THE KOFAM OUTPUT FILTERED FOR ASTERISK AND LIMITED TO TWO COLUMNS (WITHOUT NCBI VERSION)
#	awk '{ print $3"\t"$2 }' $2 > $3/ko_ncbi.tsv
#	sed -i 's/.[0-9]$//' $3/ko_ncbi.tsv

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
	wget http://ftp.flybase.org/releases/FB2024_05/precomputed_files/genes/pathway_group_data_fb_2024_05.tsv.gz -O $3/pathway_group_data_fb_2024_05.tsv.gz
	gunzip $3/pathway_group_data_fb_2024_05.tsv.gz
	grep -v ^\## $3/pathway_group_data_fb_2024_05.tsv | cut -f 1,3,6 > $3/Fbgn_groupid.tsv

	wget http://ftp.flybase.org/releases/FB2024_05/precomputed_files/genes/fbgn_annotation_ID_fb_2024_05.tsv.gz -O $3/fbgn_annotation_ID_fb_2024_05.tsv.gz
	gunzip $3/fbgn_annotation_ID_fb_2024_05.tsv.gz
	grep -v ^\## $3/fbgn_annotation_ID_fb_2024_05.tsv | cut -f 3,5 > $3/Fbgn_CG.tsv
	sed -i 's/Dmel_//g' $3/Fbgn_CG.tsv

	#PULL FB PROTEIN FASTA
        wget https://ftp.flybase.net/releases/current/dmel_r6.60/fasta/dmel-all-translation-r6.60.fasta.gz -O $3/dmel-all-translation-r6.60.fasta.gz
        gunzip $3/dmel-all-translation-r6.60.fasta.gz

	#PULL FBGN TO FBPP FILES
	wget http://ftp.flybase.org/releases/FB2024_05/precomputed_files/genes/fbgn_fbtr_fbpp_fb_2024_05.tsv.gz -O $3/fbgn_fbtr_fbpp_fb_2024_05.tsv.gz
	gunzip $3/fbgn_fbtr_fbpp_fb_2024_05.tsv.gz
	grep -v ^\## $3/fbgn_fbtr_fbpp_fb_2024_05.tsv | cut -f 1,3 > $3/Fbgn_fbpp.tsv
fi
