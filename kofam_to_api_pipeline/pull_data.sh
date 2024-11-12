#! /bin/bash

#MUST SUPPLY $1 OR $2 OR BOTH
#$1 is the KEGG species code (use NA if non-KEGG species)
#$2 is the filtered * input file ( from kofam, filter kofam)

#CAN'T PULL THIS FILE FOR SPECIES THAT AREN'T IN KEGG
if grep -q $1 kegg_org_codes.txt;
then
	wget https://rest.kegg.jp/conv/ncbi-proteinid/"$1" -O conv_ncbi-proteinid_"$1".tsv
	sed -i 's/ncbi-proteinid\://g' conv_ncbi-proteinid_"$1".tsv
	sed -i "s/$1\://g" conv_ncbi-proteinid_"$1".tsv
else
	#INSTEAD FOR NON-KEGG SPECIES USE THE KOFAM OUTPUT FILTERED FOR ASTERISK AND LIMITED TO TWO COLUMNS (WITHOUT NCBI VERSION)
	awk '{ print $3"\t"$2 }' $2 > ko_ncbi.tsv
	sed -i 's/.[0-9]$//' ko_ncbi.tsv
fi


#THESE HAVE TO BE PULLED FOR ALL SPECIES
wget https://rest.kegg.jp/link/ko/pathway -O link_ko_pathway.tsv
sed -i 's/ko\://g' link_ko_pathway.tsv
sed -i 's/path\://g' link_ko_pathway.tsv
grep -v ko link_ko_pathway.tsv > tmp
mv tmp link_ko_pathway.tsv

wget https://rest.kegg.jp/list/pathway -O list_pathway.tsv

if [ "$1" != "NA" ];
then
	#THESE CAN BE PULLED FOR A RELATED KEGG SPECIES IF DESIRED
	wget https://rest.kegg.jp/link/pathway/"$1" -O  link_pathway_"$1".tsv
	sed -i "s/$1\://g" link_pathway_"$1".tsv
	sed -i 's/path\://g' link_pathway_"$1".tsv

	wget https://rest.kegg.jp/list/pathway/"$1" -O list_pathway_"$1".tsv

	wget https://rest.kegg.jp/link/"$1"/ko -O link_"$1"_ko.tsv
	sed -i 's/ko\://g' link_"$1"_ko.tsv
	sed -i "s/$1\://g" link_"$1"_ko.tsv
fi
