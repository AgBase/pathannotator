#! /usr/bin/python3

import pandas as pd
import argparse
import os
import re
import glob
from sys import exit


parser = argparse.ArgumentParser()
parser.add_argument('pathannotator')
parser.add_argument('outdir')
parser.add_argument('outbase')

args = parser.parse_args()
pathannotator = args.pathannotator #directory containing pathannotator output
outdir = args.outdir #directory to contain results file (dir must exist)
outbase = args.outbase #basename for output file provided to pathannotator.sh

pd.set_option('display.max_columns', None)


#FIND THE OUTPUT FILES TO COMBINE
keggref = f"{pathannotator}/*KEGG_ref.tsv"
keggreffile = glob.glob(keggref)
keggreffile = str(keggreffile[0])

keggspec = f"{pathannotator}/*KEGG_species.tsv"
keggspecfile = glob.glob(keggspec)
keggspecfile = str(keggspecfile[0])

flybasetsv = f"{pathannotator}/*flybase.tsv"
flybasefile = glob.glob(flybasetsv)
flybasefile = str(flybasefile[0])

reactometsv = f"{pathannotator}/*reactome.tsv"
reactomefile = glob.glob(reactometsv)
reactomefile = str(reactomefile[0])

#READ TABLES INTO PANDAS DATAFRAMES AND ADD HEADERS TO IPRS
kr = pd.read_table(f"{keggreffile}", dtype=str)
kr.columns = ['Input_protein_ID', 'KEGG_KO', 'KEGG_ref_pathway', 'KEGG_ref_pathway_name']
ks = pd.read_table(f"{keggspecfile}", dtype=str)
ks.columns = ['Input_protein_ID', 'KEGG_KO', 'KEGG_species_pathway', 'KEGG_species_pathway_name']
fb = pd.read_table(f"{flybasefile}", dtype=str)
fb.columns = ['Input_protein_ID', 'KEGG_KO', 'Flybase_pathway_ID', 'Flybase_pathway_name']
rt = pd.read_table(f"{reactomefile}", dtype=str)
rt.columns = ['Input_protein_ID', 'UniProt_ID', 'Reactome_pathway_ID', 'Reactome_pathway_name']

#DROP UNWANTED COLUMNS (KEEP PROTEIN ACCESSIONS (1), INTERPRO ANNOTATIONS (12) AND PATHWAY ANNOTATIONS (15)
kr = kr.drop(columns=['KEGG_KO', 'KEGG_ref_pathway_name'])
ks = ks.drop(columns=['KEGG_KO', 'KEGG_species_pathway_name'])
fb = fb.drop(columns=['KEGG_KO', 'Flybase_pathway_name'])
rt = rt.drop(columns=['UniProt_ID', 'Reactome_pathway_name'])

#MAKE HEADERS MATCH FOR ALL DFS
kr.columns = ['Input_protein_ID', 'Pathway_or_domain']
ks.columns = ['Input_protein_ID', 'Pathway_or_domain']
fb.columns = ['Input_protein_ID', 'Pathway_or_domain']
rt.columns = ['Input_protein_ID', 'Pathway_or_domain']

#ADD SECOND COLUMN TO EACH DF WITH 'OPTIONAL DESCRIPTION' FOR GMT FORMAT
kr.insert(loc=1, column='Description', value='KEGG_reference_pathway')
ks.insert(loc=1, column='Description', value='KEGG_species_pathway')
fb.insert(loc=1, column='Description', value='FlyBase_pathway')
rt.insert(loc=1, column='Description', value='Reactome_pathway')

#REMOVE DUPLICATE ROWS
kr = kr.drop_duplicates()
ks = ks.drop_duplicates()
fb = fb.drop_duplicates()
rt = rt.drop_duplicates()

#JOIN DFS INTO ONE
alltogether = pd.concat([kr, ks, fb, rt])

#RUN GROUPBY
alltogether = alltogether.groupby(['Pathway_or_domain', 'Description'])['Input_protein_ID'].agg(list).reset_index()

#SPLIT LIST INTO SEPARATE COLUMNS
justlist = alltogether['Input_protein_ID'].apply(pd.Series)
alltogether = pd.concat([alltogether.drop('Input_protein_ID', axis=1), justlist], axis=1)

#PRINT TO TSV WITHOUT HEADER
alltogether.to_csv(f"{outdir}/{outbase}_all_pathways.gmt", sep='\t', header=False, index=False)

#REMOVE TRAILING TABS
gmtfile=f"{outdir}/{outbase}_all_and_pathways.gmt"
notab=f"{outdir}/gmt.tmp"


with open(gmtfile, 'r') as infile, open(notab, 'w') as outfile:
    for line in infile:
        cleaned_line = re.sub(r'\t+\n', '\n', line)
        outfile.write(cleaned_line)

os.rename(notab, gmtfile)
