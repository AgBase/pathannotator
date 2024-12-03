#! /bin/bash

#MUST SUPPLY $1 OR $2 OR BOTH
#$1 is the KEGG species code (use NA or related species if species not in KEGG)
#$2 is the filtered * input file ( from kofam, filter kofam)
#$3 output directory
#$4 ncbi status of input FASTA accessions ('ncbi' or 'non-ncbi')

#THIS WILL HAVE TO CHANGE BECAUSE SOME INPUTS WILL NOT BE NCBI ORIGIN
#RIGHT NOW I HAVE TO USE 'NA' HERE BUT 'DME' IN MERGE SCRIPT AND I NEED DME HERE TO PULL FLYBASE

kegg=$(grep $1 kegg_org_codes.txt)

if [ ! -z "${kegg}" ] && [ $4 == "ncbi" ];
then
	#CAN'T PULL THIS FILE FOR SPECIES THAT AREN'T IN KEGG
	wget https://rest.kegg.jp/conv/ncbi-proteinid/"$1" -O $3/conv_ncbi-proteinid_"$1".tsv
	sed -i 's/ncbi-proteinid\://g' $3/conv_ncbi-proteinid_"$1".tsv
	sed -i "s/$1\://g" $3/conv_ncbi-proteinid_"$1".tsv
else
	#INSTEAD FOR NON-KEGG SPECIES USE THE KOFAM OUTPUT FILTERED FOR ASTERISK AND LIMITED TO TWO COLUMNS (WITHOUT NCBI VERSION)
	awk '{ print $3"\t"$2 }' $2 > $3/ko_ncbi.tsv
	sed -i 's/.[0-9]$//' $3/ko_ncbi.tsv
fi



#THESE HAVE TO BE PULLED FOR ALL SPECIES
wget https://rest.kegg.jp/link/ko/pathway -O $3/link_ko_pathway.tsv
sed -i 's/ko\://g' $3/link_ko_pathway.tsv
sed -i 's/path\://g' $3/link_ko_pathway.tsv
grep -v ko $3/link_ko_pathway.tsv > $3/tmp
mv $3/tmp $3/link_ko_pathway.tsv

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
fi

#THIS PULLS THE DATABASE FILES FOR KOFAMSCAN
#WILL NEED TO UPDATE PATH TO MATCH CONTAINER KOFAM INSTALLATION PATH
if ( ! -f /workdir/ko_list );
then
	wget https://www.genome.jp/ftp/db/kofam/ko_list.gz -O ko_list.gz
	gunzip /workdir/ko_list.gz
fi

if ( ! -d /workdir/profiles );
then
	wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz -O profiles.tar.gz
	tar -xzf /workdir/profiles.tar.gz
fi
