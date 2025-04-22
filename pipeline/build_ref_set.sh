#! /bin/bash

#mkdir new_ref_set
cd new_ref_set
#cp /AGAT/agat_config.yaml .

#PULL FB GENOME FASTA AND GFF
#wget http://flybase-ftp.s3-website-us-east-1.amazonaws.com/genomes/Drosophila_melanogaster/dmel_r6.62_FB2025_01/gff/dmel-all-r6.62.gff.gz -O dmel-all-r6.62.gff.gz
#wget http://flybase-ftp.s3-website-us-east-1.amazonaws.com/genomes/Drosophila_melanogaster/dmel_r6.62_FB2025_01/fasta/dmel-all-chromosome-r6.62.fasta.gz -O dmel-all-chromosome-r6.62.fasta.gz
#wget http://flybase-ftp.s3-website-us-east-1.amazonaws.com/releases/FB2025_01/precomputed_files/genes/fbgn_fbtr_fbpp_fb_2025_01.tsv.gz
#gunzip -f dmel-all-r6.62.gff.gz
#gunzip -f dmel-all-chromosome-r6.62.fasta.gz
#gunzip -f fbgn_fbtr_fbpp_fb_2025_01.tsv.gz

#PROCESS FB GFF TO GET "\tFlyBase\t" LINES ONLY
#tail -n+2 dmel-all-r6.62.gff > dmel-all-r6.62.gff.tmp
#grep -P "\tFlyBase\t" dmel-all-r6.62.gff.tmp > dmel-all-r6.62_flybase_only.gff

#PULL REFSEQ GENOME FASTAS AND GFFS
#TRIBOLIUM CASTENEUM
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/031/307/605/GCF_031307605.1_icTriCast1.1/GCF_031307605.1_icTriCast1.1_genomic.fna.gz -O GCF_031307605.1_icTriCast1.1_genomic.fna.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/031/307/605/GCF_031307605.1_icTriCast1.1/GCF_031307605.1_icTriCast1.1_genomic.gff.gz -O GCF_031307605.1_icTriCast1.1_genomic.gff.gz
#gunzip -f GCF_031307605.1_icTriCast1.1_genomic.fna.gz
#gunzip -f GCF_031307605.1_icTriCast1.1_genomic.gff.gz

#APIS MELLIFERA
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz -O GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz -O GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz
#gunzip -f GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
#gunzip -f GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz

#PLODIA INTERPUNCTELLA
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/563/975/GCF_027563975.2_ilPloInte3.2/GCF_027563975.2_ilPloInte3.2_genomic.fna.gz -O GCF_027563975.2_ilPloInte3.2_genomic.fna.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/563/975/GCF_027563975.2_ilPloInte3.2/GCF_027563975.2_ilPloInte3.2_genomic.gff.gz -O GCF_027563975.2_ilPloInte3.2_genomic.gff.gz
#gunzip -f GCF_027563975.2_ilPloInte3.2_genomic.fna.gz
#gunzip -f GCF_027563975.2_ilPloInte3.2_genomic.gff.gz

#RHOPALOSIPHUM MAIDIS
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/676/215/GCF_003676215.2_ASM367621v3/GCF_003676215.2_ASM367621v3_genomic.fna.gz -O GCF_003676215.2_ASM367621v3_genomic.fna.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/676/215/GCF_003676215.2_ASM367621v3/GCF_003676215.2_ASM367621v3_genomic.gff.gz -O GCF_003676215.2_ASM367621v3_genomic.gff.gz
#gunzip -f GCF_003676215.2_ASM367621v3_genomic.fna.gz
#gunzip -f GCF_003676215.2_ASM367621v3_genomic.gff.gz

#SCHISTOCERCA SERIALIS CUBENSE
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/864/345/GCF_023864345.2_iqSchSeri2.2/GCF_023864345.2_iqSchSeri2.2_genomic.fna.gz -O GCF_023864345.2_iqSchSeri2.2_genomic.fna.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/864/345/GCF_023864345.2_iqSchSeri2.2/GCF_023864345.2_iqSchSeri2.2_genomic.gff.gz -O GCF_023864345.2_iqSchSeri2.2_genomic.gff.gz
#gunzip -f GCF_023864345.2_iqSchSeri2.2_genomic.fna.gz
#gunzip -f GCF_023864345.2_iqSchSeri2.2_genomic.gff.gz

#ISCHNURA ELEGANS
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/921/293/095/GCF_921293095.1_ioIscEleg1.1/GCF_921293095.1_ioIscEleg1.1_genomic.fna.gz -O GCF_921293095.1_ioIscEleg1.1_genomic.fna.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/921/293/095/GCF_921293095.1_ioIscEleg1.1/GCF_921293095.1_ioIscEleg1.1_genomic.gff.gz -O GCF_921293095.1_ioIscEleg1.1_genomic.gff.gz
#gunzip -f GCF_921293095.1_ioIscEleg1.1_genomic.fna.gz
#gunzip -f GCF_921293095.1_ioIscEleg1.1_genomic.gff.gz


#RUN AGAT TO GET LONGEST ISOFORM
#agat_sp_keep_longest_isoform.pl -gff GCF_921293095.1_ioIscEleg1.1_genomic.gff -o iscele_cluster.gff -c agat_config.yaml
#agat_sp_keep_longest_isoform.pl -gff GCF_023864345.2_iqSchSeri2.2_genomic.gff -o schser_cluster.gff -c agat_config.yaml
#agat_sp_keep_longest_isoform.pl -gff GCF_003676215.2_ASM367621v3_genomic.gff -o rhomai_cluster.gff -c agat_config.yaml
#agat_sp_keep_longest_isoform.pl -gff GCF_027563975.2_ilPloInte3.2_genomic.gff -o ploint_cluster.gff -c agat_config.yaml
#agat_sp_keep_longest_isoform.pl -gff GCF_003254395.2_Amel_HAv3.1_genomic.gff -o apimel_cluster.gff -c agat_config.yaml
#agat_sp_keep_longest_isoform.pl -gff GCF_031307605.1_icTriCast1.1_genomic.gff -o tricas_cluster.gff -c agat_config.yaml
#agat_sp_keep_longest_isoform.pl -gff dmel-all-r6.62_flybase_only.gff -o dromel_cluster.gff -c agat_config.yaml


#USE AGAT EXTRACT SEQUENCES TO PULL FASTA FROM REFEQ SPECIES
#agat_sp_extract_sequences.pl -g iscele_cluster.gff -f GCF_921293095.1_ioIscEleg1.1_genomic.fna -t cds -p --keep_attributes -o iscele_cluster.fa -c agat_config.yaml
#agat_sp_extract_sequences.pl -g schser_cluster.gff -f GCF_023864345.2_iqSchSeri2.2_genomic.fna -t cds -p --keep_attributes -o schser_cluster.fa -c agat_config.yaml
#agat_sp_extract_sequences.pl -g rhomai_cluster.gff -f GCF_003676215.2_ASM367621v3_genomic.fna -t cds -p --keep_attributes -o rhomai_cluster.fa -c agat_config.yaml
#agat_sp_extract_sequences.pl -g ploint_cluster.gff -f GCF_027563975.2_ilPloInte3.2_genomic.fna -t cds -p --keep_attributes -o ploint_cluster.fa -c agat_config.yaml
#agat_sp_extract_sequences.pl -g apimel_cluster.gff -f GCF_003254395.2_Amel_HAv3.1_genomic.fna -t cds -p --keep_attributes -o apimel_cluster.fa -c agat_config.yaml
#agat_sp_extract_sequences.pl -g tricas_cluster.gff -f GCF_031307605.1_icTriCast1.1_genomic.fna -t cds -p --keep_attributes -o tricas_cluster.fa -c agat_config.yaml

#SAME BUT DIFFERENT FOR DROMEL
#agat_sp_extract_sequences.pl -g dromel_cluster.gff -f dmel-all-chromosome-r6.62.fasta -t cds -p -o dromel_cluster.fa --merge -c agat_config.yaml

#REPLACE fbtr IDS WITH fbpp IDS BASED ON MAPPING FILE
grep -v '#' fbgn_fbtr_fbpp_fb_2025_01.tsv > mapping.tmp
cut -f 2,3 mapping.tmp > mapping.tsv

#REPLACING TRASNCRIPT ID WITH PROTEIN ID WITH SEQKIT
#cat mapping.tsv | while IFS="!" read -r transcript protein
#do
	sed -i 's/\s.*$//' dromel_cluster.fa
	seqkit replace -k mapping.tsv  -p '(.+)$' -r '{kv}' dromel_cluster.fa > dromel_cluster_renamed.fa
#done
#RUN ORTHOFINDER WITH SINGLE-TRANSCRIPT (OR NOT) FASTAS FROM INPUT SPECIES AND DROMEL
#tar -xvzf /OF/ref_set.tgz -C $3/orthofinder
#orthofinder -f $3/orthofinder -t $cpus -b $3/orthofinder/ref_set


#rm GCF* dmel*
