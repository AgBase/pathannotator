
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
args = parser.parse_args()
species = args.species

#species is the KEGG species code

#READ API TABLES INTO PANDAS DATAFRAMES
if os.path.exists(f"conv_ncbi-proteinid_{species}.tsv"):
    ncbi_spec = pd.read_table(f"conv_ncbi-proteinid_{species}.tsv", dtype=str)
    spec_ko = pd.read_table(f"link_{species}_ko.tsv", dtype=str)
    spec_pathway = pd.read_table(f"link_pathway_{species}.tsv", dtype=str)
    list_pathway_spec = pd.read_table(f"list_pathway_{species}.tsv", dtype=str)
    ko_pathway = pd.read_table('link_ko_pathway.tsv', dtype=str)
    pathway = pd.read_table('list_pathway.tsv', dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_spec.columns = [species, 'ncbi']
    spec_ko.columns = ['ko', species]
    ko_pathway.columns = ['pathway', 'ko']
    pathway.columns = ['pathway', 'pathname']
    spec_pathway.columns = [species, f"{species}pathway"]
    list_pathway_spec.columns = [f"{species}pathway", f"{species}pathname"]
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_spec_ko = pd.merge(ncbi_spec, spec_ko, on=species, how='left')
    ncbi_spec_ko_pathway = pd.merge(ncbi_spec_ko, ko_pathway, on='ko', how='left')
    ncbi_spec_ko_pathway_pathname = pd.merge(ncbi_spec_ko_pathway, pathway, on='pathway', how='left')
    ncbi_spec_ko_pathway_pathname.to_csv(f"{species}_direct_KEGG_ref.tsv", sep='\t', index=False)
#MERGE DATAFRAMES INTO ONE FOR { species } PATHWAYS
    ncbi_spec_ko_specpath = pd.merge(ncbi_spec_ko, spec_pathway, on=species, how='left')
    ncbi_spec_ko_specpath_specpathname = pd.merge(ncbi_spec_ko_specpath, list_pathway_spec, on=f"{species}pathway", how='left')
    ncbi_spec_ko_specpath_specpathname.to_csv(f"{species}_direct_KEGG_{species}.tsv", sep='\t', index=False)
else:
    ncbi_ko = pd.read_table("ko_ncbi.tsv", dtype=str)
    ko_pathway = pd.read_table('link_ko_pathway.tsv', dtype=str)
    pathway = pd.read_table('list_pathway.tsv', dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_ko.columns = ['ko', 'ncbi']
    ko_pathway.columns = ['pathway', 'ko']
    pathway.columns = ['pathway', 'pathname']
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_ko_pathway = pd.merge(ncbi_ko, ko_pathway, on='ko', how='left')
    ncbi_ko_pathway_pathname = pd.merge(ncbi_ko_pathway, pathway, on='pathway', how='left')
    ncbi_ko_pathway_pathname.to_csv(f"direct_KEGG_ref.tsv", sep='\t', index=False)

#==========================
#THINGS AS THEY WERE
#==============================
#ADD HEADERS TO DATAFRAME COLUMNS
#spec_ko.columns = ['ko', species]
#ko_pathway.columns = ['pathway', 'ko']
#pathway.columns = ['pathway', 'pathname']
#spec_pathway.columns = [species, f"{species}pathway"]
#list_pathway_spec.columns = [f"{species}pathway", f"{species}pathname"]

#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
#ncbi_spec_ko = pd.merge(ncbi_spec, spec_ko, on=species, how='left')
#ncbi_spec_ko_pathway = pd.merge(ncbi_spec_ko, ko_pathway, on='ko', how='left')
#ncbi_spec_ko_pathway_pathname = pd.merge(ncbi_spec_ko_pathway, pathway, on='pathway', how='left')
#ncbi_spec_ko_pathway_pathname.to_csv(f"{species}_direct_KEGG_ref.tsv", sep='\t', index=False)

#MERGE DATAFRAMES INTO ONE FOR { species } PATHWAYS
#ncbi_spec_ko_specpath = pd.merge(ncbi_spec_ko, spec_pathway, on=species, how='left')
#ncbi_spec_ko_specpath_specpathname = pd.merge(ncbi_spec_ko_specpath, list_pathway_spec, on=f"{species}pathway", how='left')
#ncbi_spec_ko_specpath_specpathname.to_csv(f"{species}_direct_KEGG_{species}.tsv", sep='\t', index=False)
