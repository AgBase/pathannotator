
#! /usr/bin/env/python

#module load python miniconda3

#conda create -c conda-forge -n venv-pandas python pandas
#conda init
#conda activate venv-pandas

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('species')
parser.add_argument('kofam')
parser.add_argument('indir')
parser.add_argument('outdir')
parser.add_argument('flybase')
args = parser.parse_args()
species = args.species # kegg speicies code, NA or related species if species not in KEGG
kofam = args.kofam # yes or no
indir = args.indir # directory with outputs from pull_data.sh
outdir = args.outdir # directory where outputs this will go
flybase = args.flybase #FB for Flybase annotations, NA for none

#READ API TABLES INTO PANDAS DATAFRAMES
if kofam == "no" and species != "NA":
    ncbi_spec = pd.read_table(f"{indir}/conv_ncbi-proteinid_{species}.tsv", dtype=str)
    spec_ko = pd.read_table(f"{indir}/link_{species}_ko.tsv", dtype=str)
    spec_pathway = pd.read_table(f"{indir}/link_pathway_{species}.tsv", dtype=str)
    list_pathway_spec = pd.read_table(f"{indir}/list_pathway_{species}.tsv", dtype=str)
    ko_pathway = pd.read_table(f"{indir}/link_ko_pathway.tsv", dtype=str)
    pathway = pd.read_table(f"{indir}/list_pathway.tsv", dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_spec.columns = ['KEGG_genes_ID', 'Input_protein_ID']
    spec_ko.columns = ['KEGG_KO', 'KEGG_genes_ID']
    ko_pathway.columns = ['KEGG_ref_pathway', 'KEGG_KO']
    pathway.columns = ['KEGG_ref_pathway', 'KEGG_ref_pathway_name']
    spec_pathway.columns = ['KEGG_genes_ID', f"KEGG_{species}_pathway"]
    list_pathway_spec.columns = [f"KEGG_{species}_pathway", f"KEGG_{species}_pathway_name"]
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_spec_ko = pd.merge(ncbi_spec, spec_ko, on='KEGG_genes_ID', how='inner')
    ncbi_spec_ko_pathway = pd.merge(ncbi_spec_ko, ko_pathway, on='KEGG_KO', how='inner')
    ncbi_spec_ko_pathway_pathname = pd.merge(ncbi_spec_ko_pathway, pathway, on='KEGG_ref_pathway', how='left')
    ncbi_spec_ko_pathway_pathname.to_csv(f"{outdir}/{species}_KEGG_ref.tsv", sep='\t', index=False)
#MERGE DATAFRAMES INTO ONE FOR { species } PATHWAYS
    ncbi_spec_ko_specpath = pd.merge(ncbi_spec_ko, spec_pathway, on='KEGG_genes_ID', how='inner')
    ncbi_spec_ko_specpath_specpathname = pd.merge(ncbi_spec_ko_specpath, list_pathway_spec, on=f"KEGG_{species}_pathway", how='left')
    ncbi_spec_ko_specpath_specpathname.to_csv(f"{outdir}/{species}_KEGG_species.tsv", sep='\t', index=False)
    #ADD FLYBASE ANNOTATIONS WHEN DME IS THE SPECIFIED SPECIES
    if flybase == "FB" and species == "dme":
        #READ INTO DATAFRAMES
        fbgn_CG = pd.read_table(f"{indir}/Fbgn_CG.tsv", dtype=str)
        fbgn_path = pd.read_table(f"{indir}/Fbgn_groupid.tsv", dtype=str)
        #ADD HEADERS
        fbgn_CG.columns = ['Flybase_gene', 'KEGG_genes_ID']
        fbgn_path.columns = ['Flybase_pathway_ID', 'Flybase_pathway_name', 'Flybase_gene']
        #MERGE AND OUTPUT TO FILE
        fbgn_CG_path = pd.merge(fbgn_CG, fbgn_path, on='Flybase_gene', how='inner')
        ncbi_spec_ko['KEGG_genes_ID'] = ncbi_spec_ko['KEGG_genes_ID'].str.replace('Dmel_', '')
        fbgn_CG_path_ncbi_spec_ko = pd.merge(ncbi_spec_ko, fbgn_CG_path, on='KEGG_genes_ID', how='inner')
        fbgn_CG_path_ncbi_spec_ko.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_CG_path_ncbi_spec_ko = fbgn_CG_path_ncbi_spec_ko[["KEGG_genes_ID","Input_protein_ID","KEGG_KO","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_CG_path_ncbi_spec_ko.to_csv(f"{outdir}/{species}_flybase.tsv", sep='\t', index=False)
    elif flybase =="FB" and species != "dme":
        #READ INTO DATAFRAMES
        fbgn_CG = pd.read_table(f"{indir}/Fbgn_CG.tsv", dtype=str)
        fbgn_path = pd.read_table(f"{indir}/Fbgn_groupid.tsv", dtype=str)
        fbgn_phmm = pd.read_table(f"{indir}/phmm_matches.txt", dtype=str)
        fbgn_fbpp = pd.read_table(f"{indir}/Fbgn_fbpp.tsv", dtype=str)
        #ADD HEADERS
        fbgn_path.columns = ['Flybase_pathway_ID', 'Flybase_pathway_name', 'Flybase_gene']
        fbgn_CG.columns = ['Flybase_gene', 'KEGG_genes_ID']
        fbgn_phmm.columns = ['Flybase_protein_ID', 'Input_protein_ID']
        fbgn_fbpp.columns = ['Flybase_gene', 'Flybase_protein_ID']
        #MERGE AND OUTPUT TO FILE
        fbgn_fbpp_phmm = pd.merge(fbgn_fbpp, fbgn_phmm, on='Flybase_protein_ID', how='inner')
        fbgn_fbpp_phmm_path = pd.merge(fbgn_fbpp_phmm, fbgn_path, on='Flybase_gene', how='inner')
        fbgn_fbpp_phmm_path_CG = pd.merge(fbgn_fbpp_phmm_path, fbgn_CG, on='Flybase_gene', how='inner')
        fbgn_fbpp_phmm_path_CG.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_fbpp_phmm_path_CG = fbgn_fbpp_phmm_path_CG[["KEGG_genes_ID","Input_protein_ID","Flybase_protein_ID","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_fbpp_phmm_path_CG.to_csv(f"{outdir}/HMM_flybase.tsv", sep='\t', index=False)
    else:
        print("You have not requested Flybase annotations.")
elif kofam == "yes" and species == "NA":
#READ API TABLES INTO PANDAS DATAFRAMES
    ncbi_ko = pd.read_table(f"{indir}/ko_ncbi.tsv", dtype=str)
    ko_pathway = pd.read_table(f"{indir}/link_ko_pathway.tsv", dtype=str)
    pathway = pd.read_table(f"{indir}/list_pathway.tsv", dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_ko.columns = ['KEGG_KO', 'Input_protein_ID']
    ko_pathway.columns = ['KEGG_ref_pathway', 'KEGG_KO']
    pathway.columns = ['KEGG_ref_pathway', 'KEGG_ref_pathway_name']
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_ko_pathway = pd.merge(ncbi_ko, ko_pathway, on='KEGG_KO', how='inner')
    ncbi_ko_pathway_pathname = pd.merge(ncbi_ko_pathway, pathway, on='KEGG_ref_pathway', how='left')
    ncbi_ko_pathway_pathname.insert(0, 'KEGG_genes_ID','NA',allow_duplicates=True)
    ncbi_ko_pathway_pathname = ncbi_ko_pathway_pathname[['KEGG_genes_ID', 'Input_protein_ID','KEGG_KO', 'KEGG_ref_pathway', 'KEGG_ref_pathway_name']]
    ncbi_ko_pathway_pathname.to_csv(f"{outdir}/NA_KEGG_ref.tsv", sep='\t', index=False)
    if flybase == "FB":
        #READ INTO DATAFRAMES
        fbgn_CG = pd.read_table(f"{indir}/Fbgn_CG.tsv", dtype=str)
        fbgn_path = pd.read_table(f"{indir}/Fbgn_groupid.tsv", dtype=str)
        fbgn_phmm = pd.read_table(f"{indir}/phmm_matches.txt", dtype=str)
        fbgn_fbpp = pd.read_table(f"{indir}/Fbgn_fbpp.tsv", dtype=str)
        #ADD HEADERS
        fbgn_path.columns = ['Flybase_pathway_ID', 'Flybase_pathway_name', 'Flybase_gene']
        fbgn_CG.columns = ['Flybase_gene', 'KEGG_genes_ID']
        fbgn_phmm.columns = ['Flybase_protein_ID', 'Input_protein_ID']
        fbgn_fbpp.columns = ['Flybase_gene', 'Flybase_protein_ID']
        #MERGE AND OUTPUT TO FILE
        fbgn_fbpp_phmm = pd.merge(fbgn_fbpp, fbgn_phmm, on='Flybase_protein_ID', how='inner')
        fbgn_fbpp_phmm_path = pd.merge(fbgn_fbpp_phmm, fbgn_path, on='Flybase_gene', how='inner')
        fbgn_fbpp_phmm_path_CG = pd.merge(fbgn_fbpp_phmm_path, fbgn_CG, on='Flybase_gene', how='inner')
        fbgn_fbpp_phmm_path_CG.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_fbpp_phmm_path_CG = fbgn_fbpp_phmm_path_CG[["KEGG_genes_ID","Input_protein_ID","Flybase_protein_ID","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_fbpp_phmm_path_CG.to_csv(f"{outdir}/HMM_flybase.tsv", sep='\t', index=False)
elif kofam == "yes" and species != "NA":
#READ API TABLES INTO PANDAS DATAFRAMES
    ncbi_ko = pd.read_table(f"{indir}/ko_ncbi.tsv", dtype=str)
    ko_pathway = pd.read_table(f"{indir}/link_ko_pathway.tsv", dtype=str)
    pathway = pd.read_table(f"{indir}/list_pathway.tsv", dtype=str)
    spec_ko = pd.read_table(f"{indir}/link_{species}_ko.tsv", dtype=str)
    spec_pathway = pd.read_table(f"{indir}/link_pathway_{species}.tsv", dtype=str)
    list_pathway_spec = pd.read_table(f"{indir}/list_pathway_{species}.tsv", dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_ko.columns = ['KEGG_KO', 'Input_protein_ID']
    spec_ko.columns = ['KEGG_KO', 'KEGG_genes_ID']
    ko_pathway.columns = ['KEGG_ref_pathway', 'KEGG_KO']
    pathway.columns = ['KEGG_ref_pathway', 'KEGG_ref_pathway_name']
    spec_pathway.columns = ['KEGG_genes_ID', f"KEGG_{species}_pathway"]
    list_pathway_spec.columns = [f"KEGG_{species}_pathway", f"KEGG_{species}_pathway_name"]
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_ko_pathway = pd.merge(ncbi_ko, ko_pathway, on='KEGG_KO', how='inner')
    ncbi_ko_pathway_pathname = pd.merge(ncbi_ko_pathway, pathway, on='KEGG_ref_pathway', how='left')
    ncbi_ko_pathway_pathname = ncbi_ko_pathway_pathname[["Input_protein_ID","KEGG_KO","KEGG_ref_pathway","KEGG_ref_pathway_name"]]
    ncbi_ko_pathway_pathname.insert(0, 'KEGG_genes_ID','NA',allow_duplicates=True)
    ncbi_ko_pathway_pathname.to_csv(f"{outdir}/{species}_KEGG_ref.tsv", sep='\t', index=False)
#MERGE DATAFRAMES INTO ONE FOR { species } PATHWAYS
    ncbi_spec_ko = pd.merge(ncbi_ko, spec_ko, on='KEGG_KO', how='inner')
    ncbi_spec_ko_specpath = pd.merge(ncbi_spec_ko, spec_pathway, on='KEGG_genes_ID', how='inner')
    ncbi_spec_ko_specpath_specpathname = pd.merge(ncbi_spec_ko_specpath, list_pathway_spec, on=f"KEGG_{species}_pathway", how='left')
    ncbi_spec_ko_specpath_specpathname = ncbi_spec_ko_specpath_specpathname[["KEGG_genes_ID","Input_protein_ID","KEGG_KO",f"KEGG_{species}_pathway",f"KEGG_{species}_pathway_name"]]
    ncbi_spec_ko_specpath_specpathname.to_csv(f"{outdir}/{species}_KEGG_species.tsv", sep='\t', index=False)
    #ADD FLYBASE ANNOTATIONS
    if flybase == "FB" and species == "dme":
        #READ INTO DATAFRAMES
        fbgn_CG = pd.read_table(f"{indir}/Fbgn_CG.tsv", dtype=str)
        fbgn_path = pd.read_table(f"{indir}/Fbgn_groupid.tsv", dtype=str)
        #ADD HEADERS
        fbgn_CG.columns = ['Flybase_gene', 'KEGG_genes_ID']
        fbgn_path.columns = ['Flybase_pathway_ID', 'Flybase_pathway_name', 'Flybase_gene']
        #MERGE AND OUTPUT TO FILE
        fbgn_CG_path = pd.merge(fbgn_CG, fbgn_path, on='Flybase_gene', how='inner')
        ncbi_spec_ko['KEGG_genes_ID'] = ncbi_spec_ko['KEGG_genes_ID'].str.replace('Dmel_', '')
        fbgn_CG_path_ncbi_spec_ko = pd.merge(ncbi_spec_ko, fbgn_CG_path, on='KEGG_genes_ID', how='inner')
        fbgn_CG_path_ncbi_spec_ko.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_CG_path_ncbi_spec_ko = fbgn_CG_path_ncbi_spec_ko[["KEGG_genes_ID","Input_protein_ID","KEGG_KO","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_CG_path_ncbi_spec_ko.to_csv(f"{outdir}/{species}_flybase.tsv", sep='\t', index=False)
    elif flybase =="FB" and species != "dme":
        #READ INTO DATAFRAMES
        fbgn_CG = pd.read_table(f"{indir}/Fbgn_CG.tsv", dtype=str)
        fbgn_path = pd.read_table(f"{indir}/Fbgn_groupid.tsv", dtype=str)
        fbgn_phmm = pd.read_table(f"{indir}/phmm_matches.txt", dtype=str)
        fbgn_fbpp = pd.read_table(f"{indir}/Fbgn_fbpp.tsv", dtype=str)
        #ADD HEADERS
        fbgn_path.columns = ['Flybase_pathway_ID', 'Flybase_pathway_name', 'Flybase_gene']
        fbgn_CG.columns = ['Flybase_gene', 'KEGG_genes_ID']
        fbgn_phmm.columns = ['Flybase_protein_ID', 'Input_protein_ID']
        fbgn_fbpp.columns = ['Flybase_gene', 'Flybase_protein_ID']
        #MERGE AND OUTPUT TO FILE
        fbgn_fbpp_phmm = pd.merge(fbgn_fbpp, fbgn_phmm, on='Flybase_protein_ID', how='inner')
        fbgn_fbpp_phmm_path = pd.merge(fbgn_fbpp_phmm, fbgn_path, on='Flybase_gene', how='inner')
        fbgn_fbpp_phmm_path_CG = pd.merge(fbgn_fbpp_phmm_path, fbgn_CG, on='Flybase_gene', how='inner')
        fbgn_fbpp_phmm_path_CG.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_fbpp_phmm_path_CG = fbgn_fbpp_phmm_path_CG[["KEGG_genes_ID","Input_protein_ID","Flybase_protein_ID","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_fbpp_phmm_path_CG.to_csv(f"{outdir}/HMM_flybase.tsv", sep='\t', index=False)
    else:
        print("You have not requested Flybase annotations.")
else:
    print("Not an acceptable combination of arguments.")
