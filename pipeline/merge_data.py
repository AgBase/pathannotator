#! /usr/bin/env/python

import pandas as pd
import argparse
import os
from sys import exit


parser = argparse.ArgumentParser()
parser.add_argument('species')
parser.add_argument('kofam')
parser.add_argument('indir')
parser.add_argument('outdir')
parser.add_argument('flybase')
parser.add_argument('orthologs')
parser.add_argument('outbase')
args = parser.parse_args()
species = args.species # kegg speicies code, NA or related species if species not in KEGG
kofam = args.kofam # yes or no
indir = args.indir # directory with outputs from pull_data.sh
outdir = args.outdir # (default is '.') directory where outputs from this will go
flybase = args.flybase # (default if 'NA') FB for Flybase and DME Reactome annotations, NA for none
orthologs = args.orthologs # $3/orthofinder/Orthologues_"$noext"-cluster/"$noext"-cluster__v__dromel-cluster.tsv from pathannotator.sh script
outbase = args.outbase # file basename for output files suppliec to the pathannotator.sh wrapper script

pd.set_option('display.max_columns', None)

#READ API TABLES INTO PANDAS DATAFRAMES
if kofam == "no" and species != "NA":
    ncbi_ver = pd.read_table(f"{indir}/ncbiver.tsv", dtype=str)
    ncbi_spec = pd.read_table(f"{indir}/conv_ncbi-proteinid_{species}.tsv", dtype=str)
    spec_ko = pd.read_table(f"{indir}/link_{species}_ko.tsv", dtype=str)
    spec_pathway = pd.read_table(f"{indir}/link_pathway_{species}.tsv", dtype=str)
    list_pathway_spec = pd.read_table(f"{indir}/list_pathway_{species}.tsv", dtype=str)
    ko_pathway = pd.read_table(f"{indir}/link_ko_pathway.tsv", dtype=str)
    pathway = pd.read_table(f"{indir}/list_pathway.tsv", dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_ver.columns = ['Input_protein_ID_version', 'Input_protein_ID']
    ncbi_spec.columns = ['KEGG_genes_ID', 'Input_protein_ID']
    spec_ko.columns = ['KEGG_KO', 'KEGG_genes_ID']
    ko_pathway.columns = ['KEGG_ref_pathway', 'KEGG_KO']
    pathway.columns = ['KEGG_ref_pathway', 'KEGG_ref_pathway_name']
    spec_pathway.columns = ['KEGG_genes_ID', f"KEGG_{species}_pathway"]
    list_pathway_spec.columns = [f"KEGG_{species}_pathway", f"KEGG_{species}_pathway_name"]
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_ver_spec = pd.merge(ncbi_ver, ncbi_spec, on='Input_protein_ID', how='inner')
    ncbi_ver_spec_ko = pd.merge(ncbi_ver_spec, spec_ko, on='KEGG_genes_ID', how='inner')
    ncbi_ver_spec_ko_pathway = pd.merge(ncbi_ver_spec_ko, ko_pathway, on='KEGG_KO', how='inner')
    ncbi_ver_spec_ko_pathway_pathname = pd.merge(ncbi_ver_spec_ko_pathway, pathway, on='KEGG_ref_pathway', how='left')
    ncbi_ver_spec_ko_pathway_pathname.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
    ncbi_ver_spec_ko_pathway_pathname.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
    ncbi_ver_spec_ko_pathway_pathname = ncbi_ver_spec_ko_pathway_pathname.drop_duplicates()
    ncbi_ver_spec_ko_pathway_pathname.to_csv(f"{outdir}/{outbase}_KEGG_ref.tsv", sep='\t', index=False)
#MERGE DATAFRAMES INTO ONE FOR { species } PATHWAYS
    ncbi_ver_spec_ko_specpath = pd.merge(ncbi_ver_spec_ko, spec_pathway, on='KEGG_genes_ID', how='inner')
    ncbi_ver_spec_ko_specpath_specpathname = pd.merge(ncbi_ver_spec_ko_specpath, list_pathway_spec, on=f"KEGG_{species}_pathway", how='left')
    ncbi_ver_spec_ko_specpath_specpathname.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
    ncbi_ver_spec_ko_specpath_specpathname.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
    ncbi_ver_spec_ko_specpath_specpathname = ncbi_ver_spec_ko_specpath_specpathname.drop_duplicates()
    ncbi_ver_spec_ko_specpath_specpathname.to_csv(f"{outdir}/{outbase}_KEGG_species.tsv", sep='\t', index=False)
#ADD FLYBASE AND REACTOME ANNOTATIONS WHEN DME IS THE SPECIFIED SPECIES
    if flybase == "FB" and species == "dme":
    #READ INTO DATAFRAMES
        fbgn_CG = pd.read_table(f"{indir}/Fbgn_CG.tsv", dtype=str)
        fbgn_path = pd.read_table(f"{indir}/Fbgn_groupid.tsv", dtype=str)
    #ADD HEADERS
        fbgn_CG.columns = ['Flybase_gene', 'KEGG_genes_ID']
        fbgn_path.columns = ['Flybase_pathway_ID', 'Flybase_pathway_name', 'Flybase_gene']
    #MERGE AND OUTPUT TO FILE
        fbgn_CG_path = pd.merge(fbgn_CG, fbgn_path, on='Flybase_gene', how='inner')
        ncbi_ver_spec_ko['KEGG_genes_ID'] = ncbi_ver_spec_ko['KEGG_genes_ID'].str.replace('Dmel_', '')
        fbgn_CG_path_ncbi_ver_spec_ko = pd.merge(ncbi_ver_spec_ko, fbgn_CG_path, on='KEGG_genes_ID', how='inner')
        fbforreact = fbgn_CG_path_ncbi_ver_spec_ko.copy()
        fbgn_CG_path_ncbi_ver_spec_ko.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko = fbgn_CG_path_ncbi_ver_spec_ko[["KEGG_genes_ID","Input_protein_ID_version","Input_protein_ID","KEGG_KO","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_CG_path_ncbi_ver_spec_ko.drop(['Input_protein_ID','KEGG_genes_ID'], axis=1, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko = fbgn_CG_path_ncbi_ver_spec_ko.drop_duplicates()
        fbgn_CG_path_ncbi_ver_spec_ko.to_csv(f"{outdir}/{outbase}_flybase.tsv", sep='\t', index=False)
    #READ REACTOME TABLES INTO DATAFRAMES
        fbuni = pd.read_table(f"{indir}/gp_information.fb", dtype=str, on_bad_lines='warn')
        unireact = pd.read_table(f"{indir}/UniProt2Reactome_DME.txt", dtype=str, header=None, on_bad_lines='warn')
    #ADD HEADERS TO DATARAME COLUMNS
        fbuni.columns = ['Database', 'Flybase_gene', 'Gene_symbol', 'Gene_name', 'KEGG_genes_ID', 'Gene_product_type', 'Taxon_ID', 'Empty', 'UniProt_ID', 'Empty2' ]
        unireact.columns = ['UniProt_ID', 'Reactome_pathway_ID', 'Pathway_URL', 'Reactome_pathway_name', 'Evidence code', 'Species']
    #DROP UNWANTED COLUMNS AND ROWS
        fbuni = fbuni[fbuni['Gene_product_type'] == 'protein']
        fbuni.drop(columns=['Database', 'Gene_symbol', 'Gene_name', 'Gene_product_type', 'Taxon_ID', 'Empty', 'Empty2'], axis=1, inplace=True)
        unireact.drop(columns=['Pathway_URL', 'Evidence code', 'Species'], axis=1, inplace=True)
    #EXPLODE UNIPROT COLUMN
        fbuni['UniProt_ID'] = fbuni['UniProt_ID'].str.split('|')
        fbuni = fbuni.explode('UniProt_ID')
        fbuni[['DB', 'UniProt_ID']] = fbuni["UniProt_ID"].str.split(":", expand=True)
        fbuni.drop('DB', axis=1, inplace=True)
    #MERGE DATAFRAMES FOR REACTOME
        fbunireact = pd.merge(fbuni, unireact, on='UniProt_ID', how='inner')
        fbuniforreact = pd.merge(fbunireact, fbforreact, on='Flybase_gene', how='inner')
        fbuniforreact = fbuniforreact[["Input_protein_ID_version","UniProt_ID","Reactome_pathway_ID","Reactome_pathway_name"]]
        fbuniforreact = fbuniforreact.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'})
        fbuniforreact = fbuniforreact.drop_duplicates()
    #WRITE TABULAR OUTPUT FOR REACTOME PATHWAYS
        fbuniforreact.to_csv(f"{outdir}/{outbase}_reactome.tsv", sep='\t', index=False)
    elif flybase =="FB" and species != "dme":
    #READ INTO DATAFRAMES
        fbgn_CG = pd.read_table(f"{indir}/Fbgn_CG.tsv", dtype=str)
        fbgn_path = pd.read_table(f"{indir}/Fbgn_groupid.tsv", dtype=str)
        fbpp_ortho = pd.read_table(f"{orthologs}", dtype=str)
        fbgn_fbpp = pd.read_table(f"{indir}/Fbgn_fbpp.tsv", dtype=str)
    #REMOVE ORTHOGROUP COLUMN
        fbpp_ortho.drop('Orthogroup', axis=1, inplace=True)
    #ADD HEADERS
        fbgn_path.columns = ['Flybase_pathway_ID', 'Flybase_pathway_name', 'Flybase_gene']
        fbgn_CG.columns = ['Flybase_gene', 'KEGG_genes_ID']
        fbpp_ortho.columns = ['Input_protein_ID', 'Flybase_protein_ID']
        fbgn_fbpp.columns = ['Flybase_gene', 'Flybase_protein_ID']
    #SPLIT ID COLUMN IN FBPP_ORTHO INTO ID AND NAME, KEEP ID
        fbpp_ortho['Input_protein_ID'] = fbpp_ortho['Input_protein_ID'].str.split(' ', n=1).str[0]
    #SPLIT AND EXPLODE TO EXPAND LISTS IN BOTH COLUMNS
        fbpp_ortho["Flybase_protein_ID"] = fbpp_ortho["Flybase_protein_ID"].str.split(", ")
        fbpp_ortho["Input_protein_ID"] = fbpp_ortho["Input_protein_ID"].str.split(", ")
        fbpp_ortho = fbpp_ortho.explode('Flybase_protein_ID').explode('Input_protein_ID')
        fbpp_ortho = fbpp_ortho.sort_values(by=['Flybase_protein_ID', 'Input_protein_ID'])
    #MERGE AND OUTPUT TO FILE
        fbgn_fbpp_ortho = pd.merge(fbgn_fbpp, fbpp_ortho, on='Flybase_protein_ID', how='inner')
        fbgn_fbpp_ortho_path = pd.merge(fbgn_fbpp_ortho, fbgn_path, on='Flybase_gene', how='inner')
        fbgn_fbpp_ortho_path_CG = pd.merge(fbgn_fbpp_ortho_path, fbgn_CG, on='Flybase_gene', how='inner')
        fbforreact = fbgn_fbpp_ortho_path_CG.copy()
        fbgn_fbpp_ortho_path_CG.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG = fbgn_fbpp_ortho_path_CG[["KEGG_genes_ID","Input_protein_ID","Flybase_protein_ID","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_fbpp_ortho_path_CG.drop('KEGG_genes_ID', axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG =fbgn_fbpp_ortho_path_CG.drop_duplicates()
        fbgn_fbpp_ortho_path_CG.to_csv(f"{outdir}/{outbase}_flybase.tsv", sep='\t', index=False)
    #READ REACTOME TABLES INTO DATAFRAMES
        fbuni = pd.read_table(f"{indir}/gp_information.fb", dtype=str, on_bad_lines='warn')
        unireact = pd.read_table(f"{indir}/UniProt2Reactome_DME.txt", dtype=str, header=None, on_bad_lines='warn')
    #ADD HEADERS TO DATARAME COLUMNS
        fbuni.columns = ['Database', 'Flybase_gene', 'Gene_symbol', 'Gene_name', 'KEGG_genes_ID', 'Gene_product_type', 'Taxon_ID', 'Empty', 'UniProt_ID', 'Empty2' ]
        unireact.columns = ['UniProt_ID', 'Reactome_pathway_ID', 'Pathway_URL', 'Reactome_pathway_name', 'Evidence code', 'Species']
    #DROP UNWANTED COLUMNS AND ROWS
        fbuni = fbuni[fbuni['Gene_product_type'] == 'protein']
        fbuni.drop(columns=['Database', 'Gene_symbol', 'Gene_name', 'Gene_product_type', 'Taxon_ID', 'Empty', 'Empty2'], axis=1, inplace=True)
        unireact.drop(columns=['Pathway_URL', 'Evidence code', 'Species'], axis=1, inplace=True)
    #EXPLODE UNIPROT COLUMN
        fbuni['UniProt_ID'] = fbuni['UniProt_ID'].str.split('|')
        fbuni = fbuni.explode('UniProt_ID')
        fbuni[['DB', 'UniProt_ID']] = fbuni["UniProt_ID"].str.split(":", expand=True)
        fbuni.drop('DB', axis=1, inplace=True)
    #MERGE DATAFRAMES FOR REACTOME
        fbunireact = pd.merge(fbuni, unireact, on='UniProt_ID', how='inner')
        fbuniforreact = pd.merge(fbunireact, fbforreact, on='Flybase_gene', how='inner')
        fbuniforreact = fbuniforreact[["Input_protein_ID","UniProt_ID","Reactome_pathway_ID","Reactome_pathway_name"]]
        fbuniforreact = fbuniforreact.drop_duplicates()
    #WRITE TABULAR OUTPUT FOR REACTOME PATHWAYS
        fbuniforreact.to_csv(f"{outdir}/{outbase}_reactome.tsv", sep='\t', index=False)
    else:
        print("You have not requested Flybase annotations.")
elif kofam == "yes" and species == "NA":
#READ API TABLES INTO PANDAS DATAFRAMES
    ncbi_ver = pd.read_table(f"{indir}/ncbiver.tsv", dtype=str)
    ncbi_ko = pd.read_table(f"{indir}/ko_ncbi.tsv", dtype=str)
    ko_pathway = pd.read_table(f"{indir}/link_ko_pathway.tsv", dtype=str)
    pathway = pd.read_table(f"{indir}/list_pathway.tsv", dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_ver.columns = ['Input_protein_ID_version', 'Input_protein_ID']
    ncbi_ko.columns = ['KEGG_KO', 'Input_protein_ID']
    ko_pathway.columns = ['KEGG_ref_pathway', 'KEGG_KO']
    pathway.columns = ['KEGG_ref_pathway', 'KEGG_ref_pathway_name']
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_ver_ko = pd.merge(ncbi_ver, ncbi_ko, on='Input_protein_ID', how='inner')
    ncbi_ver_ko_pathway = pd.merge(ncbi_ver_ko, ko_pathway, on='KEGG_KO', how='inner')
    ncbi_ver_ko_pathway_pathname = pd.merge(ncbi_ver_ko_pathway, pathway, on='KEGG_ref_pathway', how='left')
    ncbi_ver_ko_pathway_pathname = ncbi_ver_ko_pathway_pathname[["Input_protein_ID_version","Input_protein_ID","KEGG_KO","KEGG_ref_pathway","KEGG_ref_pathway_name"]]
    ncbi_ver_ko_pathway_pathname.drop('Input_protein_ID', axis=1, inplace=True)
    ncbi_ver_ko_pathway_pathname.rename(columns={"Input_protein_ID_version": "Input_protein_ID"}, inplace=True)
    ncbi_ver_ko_pathway_pathname = ncbi_ver_ko_pathway_pathname.drop_duplicates()
    ncbi_ver_ko_pathway_pathname.to_csv(f"{outdir}/{outbase}_KEGG_ref.tsv", sep='\t', index=False)
    if flybase == "FB":
    #READ INTO DATAFRAMES
        fbgn_CG = pd.read_table(f"{indir}/Fbgn_CG.tsv", dtype=str)
        fbgn_path = pd.read_table(f"{indir}/Fbgn_groupid.tsv", dtype=str)
        fbpp_ortho = pd.read_table(f"{orthologs}", dtype=str)
        fbgn_fbpp = pd.read_table(f"{indir}/Fbgn_fbpp.tsv", dtype=str)
    #REMOVE ORTHOGROUP COLUMN
        fbpp_ortho.drop('Orthogroup', axis=1, inplace=True)
    #ADD HEADERS
        fbgn_path.columns = ['Flybase_pathway_ID', 'Flybase_pathway_name', 'Flybase_gene']
        fbgn_CG.columns = ['Flybase_gene', 'KEGG_genes_ID']
        fbpp_ortho.columns = ['Input_protein_ID', 'Flybase_protein_ID']
        fbgn_fbpp.columns = ['Flybase_gene', 'Flybase_protein_ID']
    #SPLIT ID COLUMN IN FBPP_ORTHO INTO ID AND NAME, KEEP NAME
        fbpp_ortho['Input_protein_ID'] = fbpp_ortho['Input_protein_ID'].str.split(' ', n=1).str[0]
    #SPLIT AND EXPLODE TO EXPAND LISTS IN BOTH COLUMNS
        fbpp_ortho["Flybase_protein_ID"] = fbpp_ortho["Flybase_protein_ID"].str.split(", ")
        fbpp_ortho["Input_protein_ID"] = fbpp_ortho["Input_protein_ID"].str.split(", ")
        fbpp_ortho = fbpp_ortho.explode('Flybase_protein_ID').explode('Input_protein_ID')
        fbpp_ortho = fbpp_ortho.sort_values(by=['Flybase_protein_ID', 'Input_protein_ID'])
    #MERGE AND OUTPUT TO FILE
        fbgn_fbpp_ortho = pd.merge(fbgn_fbpp, fbpp_ortho, on='Flybase_protein_ID', how='inner')
        fbgn_fbpp_ortho_path = pd.merge(fbgn_fbpp_ortho, fbgn_path, on='Flybase_gene', how='inner')
        fbgn_fbpp_ortho_path_CG = pd.merge(fbgn_fbpp_ortho_path, fbgn_CG, on='Flybase_gene', how='inner')
        fbforreact = fbgn_fbpp_ortho_path_CG.copy()
        fbgn_fbpp_ortho_path_CG.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG = fbgn_fbpp_ortho_path_CG[["KEGG_genes_ID","Input_protein_ID","Flybase_protein_ID","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_fbpp_ortho_path_CG.drop('KEGG_genes_ID', axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG = fbgn_fbpp_ortho_path_CG.drop_duplicates()
        fbgn_fbpp_ortho_path_CG.to_csv(f"{outdir}/{outbase}_flybase.tsv", sep='\t', index=False)
    #READ REACTOME TABLES INTO DATAFRAMES
        fbuni = pd.read_table(f"{indir}/gp_information.fb", dtype=str, on_bad_lines='warn')
        unireact = pd.read_table(f"{indir}/UniProt2Reactome_DME.txt", dtype=str, header=None, on_bad_lines='warn')
    #ADD HEADERS TO DATARAME COLUMNS
        fbuni.columns = ['Database', 'Flybase_gene', 'Gene_symbol', 'Gene_name', 'KEGG_genes_ID', 'Gene_product_type', 'Taxon_ID', 'Empty', 'UniProt_ID', 'Empty2' ]
        unireact.columns = ['UniProt_ID', 'Reactome_pathway_ID', 'Pathway_URL', 'Reactome_pathway_name', 'Evidence code', 'Species']
    #DROP UNWANTED COLUMNS AND ROWS
        fbuni = fbuni[fbuni['Gene_product_type'] == 'protein']
        fbuni.drop(columns=['Database', 'Gene_symbol', 'Gene_name', 'Gene_product_type', 'Taxon_ID', 'Empty', 'Empty2'], axis=1, inplace=True)
        unireact.drop(columns=['Pathway_URL', 'Evidence code', 'Species'], axis=1, inplace=True)
    #EXPLODE UNIPROT COLUMN
        fbuni['UniProt_ID'] = fbuni['UniProt_ID'].str.split('|')
        fbuni = fbuni.explode('UniProt_ID')
        fbuni[['DB', 'UniProt_ID']] = fbuni["UniProt_ID"].str.split(":", expand=True)
        fbuni.drop('DB', axis=1, inplace=True)
    #MERGE DATAFRAMES FOR REACTOME
        fbunireact = pd.merge(fbuni, unireact, on='UniProt_ID', how='inner')
        fbuniforreact = pd.merge(fbunireact, fbforreact, on='Flybase_gene', how='inner')
        fbuniforreact = fbuniforreact[["Input_protein_ID","UniProt_ID","Reactome_pathway_ID","Reactome_pathway_name"]]
        fbuniforreact = fbuniforreact.drop_duplicates()
    #WRITE TABULAR OUTPUT FOR REACTOME PATHWAYS
        fbuniforreact.to_csv(f"{outdir}/{outbase}_reactome.tsv", sep='\t', index=False)
    else:
        print("You have not requested Flybase annotations.")
elif kofam == "yes" and species != "NA":
#READ API TABLES INTO PANDAS DATAFRAMES
    ncbi_ver = pd.read_table(f"{indir}/ncbiver.tsv", dtype=str)
    ncbi_ko = pd.read_table(f"{indir}/ko_ncbi.tsv", dtype=str)
    ko_pathway = pd.read_table(f"{indir}/link_ko_pathway.tsv", dtype=str)
    pathway = pd.read_table(f"{indir}/list_pathway.tsv", dtype=str)
    spec_ko = pd.read_table(f"{indir}/link_{species}_ko.tsv", dtype=str)
    spec_pathway = pd.read_table(f"{indir}/link_pathway_{species}.tsv", dtype=str)
    list_pathway_spec = pd.read_table(f"{indir}/list_pathway_{species}.tsv", dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_ver.columns = ['Input_protein_ID_version', 'Input_protein_ID']
    ncbi_ko.columns = ['KEGG_KO', 'Input_protein_ID']
    spec_ko.columns = ['KEGG_KO', 'KEGG_genes_ID']
    ko_pathway.columns = ['KEGG_ref_pathway', 'KEGG_KO']
    pathway.columns = ['KEGG_ref_pathway', 'KEGG_ref_pathway_name']
    spec_pathway.columns = ['KEGG_genes_ID', f"KEGG_{species}_pathway"]
    list_pathway_spec.columns = [f"KEGG_{species}_pathway", f"KEGG_{species}_pathway_name"]
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_ver_ko = pd.merge(ncbi_ver, ncbi_ko, on='Input_protein_ID', how='inner')
    ncbi_ver_ko_pathway = pd.merge(ncbi_ver_ko, ko_pathway, on='KEGG_KO', how='inner')
    ncbi_ver_ko_pathway_pathname = pd.merge(ncbi_ver_ko_pathway, pathway, on='KEGG_ref_pathway', how='left')
    ncbi_ver_ko_pathway_pathname = ncbi_ver_ko_pathway_pathname[["Input_protein_ID_version","Input_protein_ID","KEGG_KO","KEGG_ref_pathway","KEGG_ref_pathway_name"]]
    ncbi_ver_ko_pathway_pathname.drop('Input_protein_ID', axis=1, inplace=True)
    ncbi_ver_ko_pathway_pathname.rename(columns={"Input_protein_ID_version": "Input_protein_ID"}, inplace=True)
    ncbi_ver_ko_pathway_pathname = ncbi_ver_ko_pathway_pathname.drop_duplicates()
    ncbi_ver_ko_pathway_pathname.to_csv(f"{outdir}/{outbase}_KEGG_ref.tsv", sep='\t', index=False)
#MERGE DATAFRAMES INTO ONE FOR { species } PATHWAYS
    ncbi_ver_spec_ko = pd.merge(ncbi_ver_ko, spec_ko, on='KEGG_KO', how='inner')
    ncbi_ver_spec_ko_specpath = pd.merge(ncbi_ver_spec_ko, spec_pathway, on='KEGG_genes_ID', how='inner')
    ncbi_ver_spec_ko_specpath_specpathname = pd.merge(ncbi_ver_spec_ko_specpath, list_pathway_spec, on=f"KEGG_{species}_pathway", how='left')
    ncbi_ver_spec_ko_specpath_specpathname = ncbi_ver_spec_ko_specpath_specpathname[["KEGG_genes_ID","Input_protein_ID_version","Input_protein_ID","KEGG_KO",f"KEGG_{species}_pathway",f"KEGG_{species}_pathway_name"]]
    ncbi_ver_spec_ko_specpath_specpathname.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
    ncbi_ver_spec_ko_specpath_specpathname.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
    ncbi_ver_spec_ko_specpath_specpathname = ncbi_ver_spec_ko_specpath_specpathname.drop_duplicates()
    ncbi_ver_spec_ko_specpath_specpathname.to_csv(f"{outdir}/{outbase}_KEGG_species.tsv", sep='\t', index=False)
#ADD FLYBASE AND REACTOME ANNOTATIONS
    if flybase == "FB" and species == "dme":
    #READ INTO DATAFRAMES
        fbgn_CG = pd.read_table(f"{indir}/Fbgn_CG.tsv", dtype=str)
        fbgn_path = pd.read_table(f"{indir}/Fbgn_groupid.tsv", dtype=str)
    #ADD HEADERS
        fbgn_CG.columns = ['Flybase_gene', 'KEGG_genes_ID']
        fbgn_path.columns = ['Flybase_pathway_ID', 'Flybase_pathway_name', 'Flybase_gene']
    #MERGE AND OUTPUT TO FILE
        fbgn_CG_path = pd.merge(fbgn_CG, fbgn_path, on='Flybase_gene', how='inner')
        ncbi_ver_spec_ko['KEGG_genes_ID'] = ncbi_ver_spec_ko['KEGG_genes_ID'].str.replace('Dmel_', '')
        fbgn_CG_path_ncbi_ver_spec_ko = pd.merge(ncbi_ver_spec_ko, fbgn_CG_path, on='KEGG_genes_ID', how='inner')
        fbforreact = fbgn_CG_path_ncbi_ver_spec_ko.copy()
        fbgn_CG_path_ncbi_ver_spec_ko.drop(['Input_protein_ID', 'KEGG_genes_ID', 'Flybase_gene'], axis=1, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko = fbgn_CG_path_ncbi_ver_spec_ko.drop_duplicates()
        fbgn_CG_path_ncbi_ver_spec_ko.to_csv(f"{outdir}/{outbase}_flybase.tsv", sep='\t', index=False)
    #READ REACTOME TABLES INTO DATAFRAMES
        fbuni = pd.read_table(f"{indir}/gp_information.fb", dtype=str, on_bad_lines='warn')
        unireact = pd.read_table(f"{indir}/UniProt2Reactome_DME.txt", dtype=str, header=None, on_bad_lines='warn')
    #ADD HEADERS TO DATARAME COLUMNS
        fbuni.columns = ['Database', 'Flybase_gene', 'Gene_symbol', 'Gene_name', 'KEGG_genes_ID', 'Gene_product_type', 'Taxon_ID', 'Empty', 'UniProt_ID', 'Empty2' ]
        unireact.columns = ['UniProt_ID', 'Reactome_pathway_ID', 'Pathway_URL', 'Reactome_pathway_name', 'Evidence code', 'Species']
    #DROP UNWANTED COLUMNS AND ROWS
        fbuni = fbuni[fbuni['Gene_product_type'] == 'protein']
        fbuni.drop(columns=['Database', 'Gene_symbol', 'Gene_name', 'Gene_product_type', 'Taxon_ID', 'Empty', 'Empty2'], axis=1, inplace=True)
        unireact.drop(columns=['Pathway_URL', 'Evidence code', 'Species'], axis=1, inplace=True)
    #EXPLODE UNIPROT COLUMN
        fbuni['UniProt_ID'] = fbuni['UniProt_ID'].str.split('|')
        fbuni = fbuni.explode('UniProt_ID')
        fbuni[['DB', 'UniProt_ID']] = fbuni["UniProt_ID"].str.split(":", expand=True)
        fbuni.drop('DB', axis=1, inplace=True)
    #MERGE DATAFRAMES FOR REACTOME
        fbunireact = pd.merge(fbuni, unireact, on='UniProt_ID', how='inner')
        fbuniforreact = pd.merge(fbunireact, fbforreact, on='Flybase_gene', how='inner')
        fbuniforreact = fbuniforreact[["Input_protein_ID_version","UniProt_ID","Reactome_pathway_ID","Reactome_pathway_name"]]
        fbuniforreact = fbuniforreact.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'})
        fbuniforreact = fbuniforreact.drop_duplicates()
    #WRITE TABULAR OUTPUT FOR REACTOME PATHWAYS
        fbuniforreact.to_csv(f"{outdir}/{outbase}_reactome.tsv", sep='\t', index=False)
    elif flybase =="FB" and species != "dme":
    #READ INTO DATAFRAMES
        fbgn_CG = pd.read_table(f"{indir}/Fbgn_CG.tsv", dtype=str)
        fbgn_path = pd.read_table(f"{indir}/Fbgn_groupid.tsv", dtype=str)
        fbpp_ortho = pd.read_table(f"{orthologs}", dtype=str)
        fbgn_fbpp = pd.read_table(f"{indir}/Fbgn_fbpp.tsv", dtype=str)
    #REMOVE ORTHOGROUP COLUMN
        fbpp_ortho.drop('Orthogroup', axis=1, inplace=True)
    #ADD HEADERS
        fbgn_path.columns = ['Flybase_pathway_ID', 'Flybase_pathway_name', 'Flybase_gene']
        fbgn_CG.columns = ['Flybase_gene', 'KEGG_genes_ID']
        fbpp_ortho.columns = ['Input_protein_ID', 'Flybase_protein_ID']
        fbgn_fbpp.columns = ['Flybase_gene', 'Flybase_protein_ID']
    #SPLIT ID COLUMN IN FBPP_ORTHO INTO ID AND NAME, KEEP ID
        fbpp_ortho['Input_protein_ID'] = fbpp_ortho['Input_protein_ID'].str.split(' ', n=1).str[0]
    #SPLIT AND EXPLODE TO EXPAND LISTS IN BOTH COLUMNS
        fbpp_ortho["Flybase_protein_ID"] = fbpp_ortho["Flybase_protein_ID"].str.split(", ")
        fbpp_ortho["Input_protein_ID"] = fbpp_ortho["Input_protein_ID"].str.split(", ")
        fbpp_ortho = fbpp_ortho.explode('Flybase_protein_ID').explode('Input_protein_ID')
        fbpp_ortho = fbpp_ortho.sort_values(by=['Flybase_protein_ID', 'Input_protein_ID'])
    #MERGE AND OUTPUT TO FILE
        fbgn_fbpp_ortho = pd.merge(fbgn_fbpp, fbpp_ortho, on='Flybase_protein_ID', how='inner')
        fbgn_fbpp_ortho_path = pd.merge(fbgn_fbpp_ortho, fbgn_path, on='Flybase_gene', how='inner')
        fbgn_fbpp_ortho_path_CG = pd.merge(fbgn_fbpp_ortho_path, fbgn_CG, on='Flybase_gene', how='inner')
        fbforreact = fbgn_fbpp_ortho_path_CG.copy()
        fbgn_fbpp_ortho_path_CG.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG = fbgn_fbpp_ortho_path_CG[["KEGG_genes_ID","Input_protein_ID","Flybase_protein_ID","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_fbpp_ortho_path_CG.drop('KEGG_genes_ID', axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG = fbgn_fbpp_ortho_path_CG.drop_duplicates()
        fbgn_fbpp_ortho_path_CG.to_csv(f"{outdir}/{outbase}_flybase.tsv", sep='\t', index=False)
    #READ REACTOME TABLES INTO DATAFRAMES
        fbuni = pd.read_table(f"{indir}/gp_information.fb", dtype=str, on_bad_lines='warn')
        unireact = pd.read_table(f"{indir}/UniProt2Reactome_DME.txt", dtype=str, header=None, on_bad_lines='warn')
    #ADD HEADERS TO DATARAME COLUMNS
        fbuni.columns = ['Database', 'Flybase_gene', 'Gene_symbol', 'Gene_name', 'KEGG_genes_ID', 'Gene_product_type', 'Taxon_ID', 'Empty', 'UniProt_ID', 'Empty2' ]
        unireact.columns = ['UniProt_ID', 'Reactome_pathway_ID', 'Pathway_URL', 'Reactome_pathway_name', 'Evidence code', 'Species']
    #DROP UNWANTED COLUMNS AND ROWS
        fbuni = fbuni[fbuni['Gene_product_type'] == 'protein']
        fbuni.drop(columns=['Database', 'Gene_symbol', 'Gene_name', 'Gene_product_type', 'Taxon_ID', 'Empty', 'Empty2'], axis=1, inplace=True)
        unireact.drop(columns=['Pathway_URL', 'Evidence code', 'Species'], axis=1, inplace=True)
    #EXPLODE UNIPROT COLUMN
        fbuni['UniProt_ID'] = fbuni['UniProt_ID'].str.split('|')
        fbuni = fbuni.explode('UniProt_ID')
        fbuni[['DB', 'UniProt_ID']] = fbuni["UniProt_ID"].str.split(":", expand=True)
        fbuni.drop('DB', axis=1, inplace=True)
    #MERGE DATAFRAMES FOR REACTOME
        fbunireact = pd.merge(fbuni, unireact, on='UniProt_ID', how='inner')
        fbuniforreact = pd.merge(fbunireact, fbforreact, on='Flybase_gene', how='inner')
        fbuniforreact = fbuniforreact[["Input_protein_ID","UniProt_ID","Reactome_pathway_ID","Reactome_pathway_name"]]
        fbuniforreact = fbuniforreact.drop_duplicates()
    #WRITE TABULAR OUTPUT FOR REACTOME PATHWAYS
        fbuniforreact.to_csv(f"{outdir}/{outbase}_reactome.tsv", sep='\t', index=False)
    else:
        print("You have not requested Flybase annotations.")
else:
    print("Not an acceptable combination of arguments.")

