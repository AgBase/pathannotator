#! /usr/bin/python3

import pandas as pd
import argparse
import os
import re
import glob
from sys import exit


parser = argparse.ArgumentParser()
parser.add_argument('pathannotator')
parser.add_argument('interproscan')
parser.add_argument('outdir')

args = parser.parse_args()
pathannotator = args.pathannotator #directory containing pathannotator output
interproscan = args.interproscan #directory containing interproscan output
outdir = args.outdir #directory to contain results file (dir must exist)

pd.set_option('display.max_columns', None)

#FIND THE OUTPUT FILES TO COMBINE
iprstsv = f"{interproscan}/*.tsv"
iprsfile = glob.glob(iprstsv)
iprsfile = str(iprsfile[0])

keggref = f"{pathannotator}/*KEGG_ref.tsv"
keggreffile = glob.glob(keggref)
keggreffile = str(keggreffile[0])

keggspec = f"{pathannotator}/*KEGG_species.tsv"
keggspecfile = glob.glob(keggspec)
keggspecfile = str(keggspecfile[0])

flybasetsv = f"{pathannotator}/*flybase.tsv"
flybasefile = glob.glob(flybasetsv)
flybasefile = str(flybasefile[0])

#GET BASE FILE NAME FROM IPRS OUTPUT TO USE AS NAME OF GMT
basename = os.path.basename(iprsfile)
basename = os.path.splitext(basename)
base = basename[0]

#READ TABLES INTO PANDAS DATAFRAMES AND ADD HEADERS TO IPRS
iprs = pd.read_table(f"{iprsfile}", dtype=str)
iprs.columns = ['Protein_accession', 'Sequence_MD5_digest', 'Sequence_length', 'Analysis', 'Signature_accession', 'Signature_description', 'Start_location', 'Stop_location', 'Score', 'Status', 'Date', 'Interpro_annotation_accession', 'Interpro_annotation_description', 'GO_annotations', 'Pathway_annotations']
kr = pd.read_table(f"{keggreffile}", dtype=str)
ks = pd.read_table(f"{keggspecfile}", dtype=str)
ks.columns = ['Input_protein_ID', 'KEGG_KO', 'KEGG_species_pathway', 'KEGG_species_pathway_name']
fb = pd.read_table(f"{flybasefile}", dtype=str)

#DROP UNWANTED COLUMNS (KEEP PROTEIN ACCESSIONS (1), INTERPRO ANNOTATIONS (12) AND PATHWAY ANNOTATIONS (15)
iprsdom = iprs.drop(columns=['Sequence_MD5_digest', 'Sequence_length', 'Analysis', 'Signature_accession', 'Signature_description', 'Start_location', 'Stop_location', 'Score', 'Status', 'Date', 'Interpro_annotation_description', 'GO_annotations', 'Pathway_annotations'])
iprspath = iprs.drop(columns=['Sequence_MD5_digest', 'Sequence_length', 'Analysis', 'Signature_accession', 'Signature_description', 'Start_location', 'Stop_location', 'Score', 'Status', 'Date', 'Interpro_annotation_accession', 'Interpro_annotation_description', 'GO_annotations'])
kr = kr.drop(columns=['KEGG_KO', 'KEGG_ref_pathway_name'])
ks = ks.drop(columns=['KEGG_KO', 'KEGG_species_pathway_name'])
fb = fb.drop(columns=['Flybase_protein_ID', 'Flybase_pathway_name'])

#MAKE HEADERS MATCH FOR ALL DFS
iprsdom.columns = ['Input_protein_ID', 'Pathway_or_domain']
iprspath.columns = ['Input_protein_ID', 'Pathway_or_domain']
kr.columns = ['Input_protein_ID', 'Pathway_or_domain']
ks.columns = ['Input_protein_ID', 'Pathway_or_domain']
fb.columns = ['Input_protein_ID', 'Pathway_or_domain']

#EXPLODE PATHWAY ANNOTATION IPRS COLUMN
iprspath= iprspath.astype(str)
iprspath = iprspath.groupby('Input_protein_ID')['Pathway_or_domain'].apply(lambda x: '|'.join(x)).reset_index()
iprspath['Pathway_or_domain'] = iprspath['Pathway_or_domain'].str.split('|')
iprspath = iprspath.explode('Pathway_or_domain')

#ADD SECOND COLUMN TO EACH DF WITH 'OPTIONAL DESCRIPTION' FOR GMT FORMAT
iprsdom.insert(loc=1, column='Description', value='Interpro_domain')
kr.insert(loc=1, column='Description', value='KEGG_reference_pathway')
ks.insert(loc=1, column='Description', value='KEGG_species_pathway')
fb.insert(loc=1, column='Description', value='FlyBase_pathway')

iprspath[['Description', 'Pathway_or_domain']] = iprspath["Pathway_or_domain"].str.split(":", expand=True)

#REMOVE ROWS WITH - OR EMPTY COLUMNS
indices_to_drop = iprspath[iprspath['Pathway_or_domain'] == '-'].index
iprspath = iprspath.drop(indices_to_drop)

indices_to_drop = iprsdom[iprsdom['Pathway_or_domain'] == '-'].index
iprsdom = iprsdom.drop(indices_to_drop)

#REMOVE DUPLICATE ROWS
iprspath = iprspath.drop_duplicates()
iprsdom = iprsdom.drop_duplicates()
kr = kr.drop_duplicates()
ks = ks.drop_duplicates()
fb = fb.drop_duplicates()

#JOIN DFS INTO ONE
alltogether = pd.concat([iprsdom, iprspath, kr, ks, fb])

#RUN GROUPBY
alltogether = alltogether.groupby(['Pathway_or_domain', 'Description'])['Input_protein_ID'].agg(list).reset_index()

#SPLIT LIST INTO SEPARATE COLUMNS
justlist = alltogether['Input_protein_ID'].apply(pd.Series)
alltogether = pd.concat([alltogether.drop('Input_protein_ID', axis=1), justlist], axis=1)

#PRINT TO TSV WITHOUT HEADER
alltogether.to_csv(f"{outdir}/{base}_domains_and_pathways.gmt", sep='\t', header=False, index=False)

#REMOVE TRAILING TABS

gmtfile=f"{outdir}/{base}_domains_and_pathways.gmt"
notab=f"{outdir}/gmt.tmp"


with open(gmtfile, 'r') as infile, open(notab, 'w') as outfile:
    for line in infile:
        print (line)
        cleaned_line = re.sub(r'\t+\n', '\n', line)
        print (cleaned_line)
        outfile.write(cleaned_line)

os.rename(notab, gmtfile)
