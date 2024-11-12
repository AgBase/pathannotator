#! /usr/bin/env/python

#module load python miniconda3

#conda create -c conda-forge -n venv-pandas python pandas
#conda init
#conda activate venv-pandas

import pandas as pd

#READ API TABLES INTO PANDAS DATAFRAMES
ncbi_ame = pd.read_table('conv_ncbi-proteinid_ame.tsv', dtype=str)
ame_ko = pd.read_table('link_ame_ko.tsv', dtype=str)
ko_pathway = pd.read_table('link_ko_pathway.tsv', dtype=str)
pathway = pd.read_table('list_pathway.tsv', dtype=str)
ame_pathway = pd.read_table('link_pathway_ame.tsv', dtype=str)
list_pathway_ame = pd.read_table('list_pathway_ame.tsv', dtype=str)

#ADD HEADERS TO DATAFRAME COLUMNS
ncbi_ame.columns = ['ame', 'ncbi']
ame_ko.columns = ['ko', 'ame']
ko_pathway.columns = ['pathway', 'ko']
pathway.columns = ['pathway', 'pathname']
ame_pathway.columns = ['ame', 'amepathway']
list_pathway_ame.columns = ['amepathway', 'amepathname']

#MERGE DATAFRAMES INTO ONE FOR REFERENCE PATHWAYS
ncbi_ame_ko = pd.merge(ncbi_ame, ame_ko, on='ame', how='left')

ncbi_ame_ko_pathway = pd.merge(ncbi_ame_ko, ko_pathway, on='ko', how='left')

ncbi_ame_ko_pathway_pathname = pd.merge(ncbi_ame_ko_pathway, pathway, on='pathway', how='left')

ncbi_ame_ko_pathway_pathname.to_csv('ame_direct_KEGG_ref.tsv', sep='\t', index=False)

#MERGE DATAFRAMES INTO ONE FOR ame PATHWAYS
ncbi_ame_ko_amepath = pd.merge(ncbi_ame_ko, ame_pathway, on='ame', how='left')

ncbi_ame_ko_amepath_amepathname = pd.merge(ncbi_ame_ko_amepath, list_pathway_ame, on='amepathway', how='left')

ncbi_ame_ko_amepath_amepathname.to_csv('ame_direct_KEGG_ame.tsv', sep='\t', index=False)
