
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
args = parser.parse_args()
species = args.species # kegg speicies code, NA or related species if species not in KEGG
kofam = args.kofam # yes or no
indir = args.indir # directory with outputs from pull_data.sh
outdir = args.outdir # directory where outputs this will go

#READ API TABLES INTO PANDAS DATAFRAMES
if kofam == "no" and species != "NA":
    ncbi_spec = pd.read_table(f"{indir}/conv_ncbi-proteinid_{species}.tsv", dtype=str)
    spec_ko = pd.read_table(f"{indir}/link_{species}_ko.tsv", dtype=str)
    spec_pathway = pd.read_table(f"{indir}/link_pathway_{species}.tsv", dtype=str)
    list_pathway_spec = pd.read_table(f"{indir}/list_pathway_{species}.tsv", dtype=str)
    ko_pathway = pd.read_table(f"{indir}/link_ko_pathway.tsv", dtype=str)
    pathway = pd.read_table(f"{indir}/list_pathway.tsv", dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_spec.columns = ['species', 'ncbi']
    spec_ko.columns = ['ko', 'species']
    ko_pathway.columns = ['pathway', 'ko']
    pathway.columns = ['pathway', 'pathname']
    spec_pathway.columns = ['species', f"{species}pathway"]
    list_pathway_spec.columns = [f"{species}pathway", f"{species}pathname"]
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_spec_ko = pd.merge(ncbi_spec, spec_ko, on='species', how='inner')
    ncbi_spec_ko_pathway = pd.merge(ncbi_spec_ko, ko_pathway, on='ko', how='inner')
    ncbi_spec_ko_pathway_pathname = pd.merge(ncbi_spec_ko_pathway, pathway, on='pathway', how='left')
    ncbi_spec_ko_pathway_pathname.to_csv(f"{outdir}/{species}_direct_KEGG_ref.tsv", sep='\t', index=False)
#MERGE DATAFRAMES INTO ONE FOR { species } PATHWAYS
    ncbi_spec_ko_specpath = pd.merge(ncbi_spec_ko, spec_pathway, on='species', how='inner')
    ncbi_spec_ko_specpath_specpathname = pd.merge(ncbi_spec_ko_specpath, list_pathway_spec, on=f"{species}pathway", how='left')
    ncbi_spec_ko_specpath_specpathname.to_csv(f"{outdir}/{species}_direct_KEGG_{species}.tsv", sep='\t', index=False)
elif kofam == "yes" and species == "NA":
#READ API TABLES INTO PANDAS DATAFRAMES
    ncbi_ko = pd.read_table(f"{indir}/ko_ncbi.tsv", dtype=str)
    ko_pathway = pd.read_table(f"{indir}/link_ko_pathway.tsv", dtype=str)
    pathway = pd.read_table(f"{indir}/list_pathway.tsv", dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_ko.columns = ['ko', 'ncbi']
    ko_pathway.columns = ['pathway', 'ko']
    pathway.columns = ['pathway', 'pathname']
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_ko_pathway = pd.merge(ncbi_ko, ko_pathway, on='ko', how='inner')
    ncbi_ko_pathway_pathname = pd.merge(ncbi_ko_pathway, pathway, on='pathway', how='left')
    ncbi_ko_pathway_pathname.to_csv(f"{outdir}/direct_KEGG_ref.tsv", sep='\t', index=False)
elif kofam == "yes" and species != "NA":
#READ API TABLES INTO PANDAS DATAFRAMES
    ncbi_ko = pd.read_table(f"{indir}/ko_ncbi.tsv", dtype=str)
    ko_pathway = pd.read_table(f"{indir}/link_ko_pathway.tsv", dtype=str)
    pathway = pd.read_table(f"{indir}/list_pathway.tsv", dtype=str)
    spec_ko = pd.read_table(f"{indir}/link_{species}_ko.tsv", dtype=str)
    spec_pathway = pd.read_table(f"{indir}/link_pathway_{species}.tsv", dtype=str)
    list_pathway_spec = pd.read_table(f"{indir}/list_pathway_{species}.tsv", dtype=str)
#ADD HEADERS TO DATAFRAME COLUMNS
    ncbi_ko.columns = ['ko', 'ncbi']
    spec_ko.columns = ['ko', 'species']
    ko_pathway.columns = ['pathway', 'ko']
    pathway.columns = ['pathway', 'pathname']
    spec_pathway.columns = ['species', f"{species}pathway"]
    list_pathway_spec.columns = [f"{species}pathway", f"{species}pathname"]
#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
    ncbi_ko_pathway = pd.merge(ncbi_ko, ko_pathway, on='ko', how='inner')
    ncbi_ko_pathway_pathname = pd.merge(ncbi_ko_pathway, pathway, on='pathway', how='left')
    ncbi_ko_pathway_pathname = ncbi_ko_pathway_pathname[["ncbi","ko","pathway","pathname"]]
    ncbi_ko_pathway_pathname.to_csv(f"{outdir}/{species}_direct_KEGG_ref.tsv", sep='\t', index=False)
#MERGE DATAFRAMES INTO ONE FOR { species } PATHWAYS
    ncbi_spec_ko = pd.merge(ncbi_ko, spec_ko, on='ko', how='inner')
    ncbi_spec_ko_specpath = pd.merge(ncbi_spec_ko, spec_pathway, on='species', how='inner')
    ncbi_spec_ko_specpath_specpathname = pd.merge(ncbi_spec_ko_specpath, list_pathway_spec, on=f"{species}pathway", how='left')
    ncbi_spec_ko_specpath_specpathname = ncbi_spec_ko_specpath_specpathname[["species","ncbi","ko",f"{species}pathway",f"{species}pathname"]]
    ncbi_spec_ko_specpath_specpathname.to_csv(f"{outdir}/{species}_direct_KEGG_{species}.tsv", sep='\t', index=False)
else:
    "Not an acceptable combination of arguments."
