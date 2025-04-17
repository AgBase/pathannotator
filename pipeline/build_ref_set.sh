#! /bin/bash

#PULL FB GENOME FASTA AND GFF
#wget http://flybase-ftp.s3-website-us-east-1.amazonaws.com/genomes/Drosophila_melanogaster/dmel_r6.62_FB2025_01/gff/dmel-all-r6.62.gff.gz -O dmel-all-r6.62.gff.gz
#wget http://flybase-ftp.s3-website-us-east-1.amazonaws.com/genomes/Drosophila_melanogaster/dmel_r6.62_FB2025_01/fasta/dmel-all-chromosome-r6.62.fasta.gz -O dmel-all-chromosome-r6.62.fasta.gz
#gunzip -f dmel-all-r6.62.gff.gz
#gunzip -f dmel-all-chromosome-r6.62.fasta.gz

#PROCESS FB GFF TO GET "\tFlyBase\t" LINES ONLY
#tail -n+2 dmel-all-r6.62.gff > dmel-all-r6.62.gff.tmp
#grep -P "\tFlyBase\t" dmel-all-r6.62.gff.tmp > dmel-all-r6.62_flybase_only.gff

#PULL REFSEQ GENOME FASTAS AND GFFS
#TRIBOLIUM CASTENEUM
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/031/307/605/GCF_031307605.1_icTriCast1.1/GCF_031307605.1_icTriCast1.1_genomic.fna.gz -O GCF_031307605.1_icTriCast1.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/031/307/605/GCF_031307605.1_icTriCast1.1/GCF_031307605.1_icTriCast1.1_genomic.gff.gz -O GCF_031307605.1_icTriCast1.1_genomic.gff.gz
gunzip -f GCF_031307605.1_icTriCast1.1_genomic.fna.gz
gunzip -f GCF_031307605.1_icTriCast1.1_genomic.gff.gz

#APIS MELLIFERA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz -O GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz -O GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz
gunzip -f GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
gunzip -f GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz

#PLODIA INTERPUNCTELLA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/563/975/GCF_027563975.2_ilPloInte3.2/GCF_027563975.2_ilPloInte3.2_genomic.fna.gz -O GCF_027563975.2_ilPloInte3.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/563/975/GCF_027563975.2_ilPloInte3.2/GCF_027563975.2_ilPloInte3.2_genomic.gff.gz -O GCF_027563975.2_ilPloInte3.2_genomic.gff.gz
gunzip -f GCF_027563975.2_ilPloInte3.2_genomic.fna.gz
gunzip -f GCF_027563975.2_ilPloInte3.2_genomic.gff.gz

#RHOPALOSIPHUM MAIDIS
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/676/215/GCF_003676215.2_ASM367621v3/GCF_003676215.2_ASM367621v3_genomic.fna.gz -O GCF_003676215.2_ASM367621v3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/676/215/GCF_003676215.2_ASM367621v3/GCF_003676215.2_ASM367621v3_genomic.gff.gz -O GCF_003676215.2_ASM367621v3_genomic.gff.gz
gunzip -f GCF_003676215.2_ASM367621v3_genomic.fna.gz
gunzip -f GCF_003676215.2_ASM367621v3_genomic.gff.gz

#SCHISTOCERCA SERIALIS CUBENSE
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/864/345/GCF_023864345.2_iqSchSeri2.2/GCF_023864345.2_iqSchSeri2.2_genomic.fna.gz -O GCF_023864345.2_iqSchSeri2.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/864/345/GCF_023864345.2_iqSchSeri2.2/GCF_023864345.2_iqSchSeri2.2_genomic.gff.gz -O GCF_023864345.2_iqSchSeri2.2_genomic.gff.gz
gunzip -f GCF_023864345.2_iqSchSeri2.2_genomic.fna.gz
gunzip -f GCF_023864345.2_iqSchSeri2.2_genomic.gff.gz

#ISCHNURA ELEGANS
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/921/293/095/GCF_921293095.1_ioIscEleg1.1/GCF_921293095.1_ioIscEleg1.1_genomic.fna.gz -O GCF_921293095.1_ioIscEleg1.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/921/293/095/GCF_921293095.1_ioIscEleg1.1/GCF_921293095.1_ioIscEleg1.1_genomic.gff.gz -O GCF_921293095.1_ioIscEleg1.1_genomic.gff.gz
gunzip -f GCF_921293095.1_ioIscEleg1.1_genomic.fna.gz
gunzip -f GCF_921293095.1_ioIscEleg1.1_genomic.gff.gz


#RUN AGAT TO GET LONGEST ISOFORM

agat_sp_keep_longest_isoform.pl -gff GCF_921293095.1_ioIscEleg1.1_genomic.gff -o iscele_cluster.gff
agat_sp_keep_longest_isoform.pl -gff GCF_023864345.2_iqSchSeri2.2_genomic.gff -o schser_cluster.gff
agat_sp_keep_longest_isoform.pl -gff GCF_003676215.2_ASM367621v3_genomic.gff -o rhomai_cluster.gff
agat_sp_keep_longest_isoform.pl -gff GCF_027563975.2_ilPloInte3.2_genomic.gff -o ploint_cluster.gff
agat_sp_keep_longest_isoform.pl -gff GCF_003254395.2_Amel_HAv3.1_genomic.gff -o apimel_cluster.gff
agat_sp_keep_longest_isoform.pl -gff GCF_031307605.1_icTriCast1.1_genomic.gff -o tricas_cluster.gff
#agat_sp_keep_longest_isoform.pl -gff dmel-all-r6.62_flybase_only.gff -o dromel_cluster.gff


#USE AGAT EXTRACT SEQUENCES TO PULL FASTA FROM REFEQ SPECIES
agat_sp_extract_sequences.pl -g iscele_cluster.gff -f GCF_921293095.1_ioIscEleg1.1_genomic.fna -t cds -p --keep_attributes -o iscele_cluster.fa
agat_sp_extract_sequences.pl -g schser_cluster.gff -f GCF_023864345.2_iqSchSeri2.2_genomic.fna -t cds -p --keep_attributes -o schser_cluster.fa
agat_sp_extract_sequences.pl -g rhomai_cluster.gff -f GCF_003676215.2_ASM367621v3_genomic.fna -t cds -p --keep_attributes -o rhomai_cluster.fa
agat_sp_extract_sequences.pl -g ploint_cluster.gff -f GCF_027563975.2_ilPloInte3.2_genomic.fna -t cds -p --keep_attributes -o ploint_cluster.fa
agat_sp_extract_sequences.pl -g apimel_cluster.gff -f GCF_003254395.2_Amel_HAv3.1_genomic.fna -t cds -p --keep_attributes -o apimel_cluster.fa
agat_sp_extract_sequences.pl -g tricas_cluster.gff -f GCF_031307605.1_icTriCast1.1_genomic.fna -t cds -p --keep_attributes -o tricas_cluster.fa

#SAME BUT DIFFERENT FOR DROMEL
#agat_sp_extract_sequences.pl -g tricas_cluster.gff -f GCF_031307605.1_icTriCast1.1_genomic.fna -t cds -p --keep_attributes -o tricas_cluster.fa


#USE GFFREAD TO GET PROTEIN FASTA FOR LONGEST ISOFORM--MAY NEED TO TRANSLATE DMEL USING MAPPING FILE

#gffread -C -g GCF_921293095.1_ioIscEleg1.1_genomic.fna -y iscele_cluster.fa iscele_cluster.gff
#gffread -C -g GCF_023864345.2_iqSchSeri2.2_genomic.fna -y schser_cluster.fa schser_cluster.gff
#gffread -C -g GCF_003676215.2_ASM367621v3_genomic.fna -y rhomai_cluster.fa rhomai_cluster.gff
#gffread -C -g GCF_027563975.2_ilPloInte3.2_genomic.fna -y ploint_cluster.fa ploint_cluster.gff
#gffread -C -g GCF_003254395.2_Amel_HAv3.1_genomic.fna -y apimel_cluster.fa apimel_cluster.gff
#gffread -C -g GCF_031307605.1_icTriCast1.1_genomic.fna -y tricas_cluster.fa tricas_cluster.gff
#gffread -C -g dmel-all-chromosome-r6.62.fasta -y dromel_cluster.fa dromel_cluster.gff

#RUN ORTHOFINDER WITH SINGLE-TRANCRIPT (OR NOT) FASTAS FROM INPUT SPECIES AND DROMEL
#tar -xvzf /OF/ref_set.tgz -C $3/orthofinder
#orthofinder -f $3/orthofinder -t $cpus -b $3/orthofinder/ref_set
