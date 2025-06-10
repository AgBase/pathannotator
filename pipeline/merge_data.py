
#! /usr/bin/env/python

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('species')
parser.add_argument('kofam')
parser.add_argument('indir')
parser.add_argument('outdir')
parser.add_argument('flybase')
parser.add_argument('orthologs')
args = parser.parse_args()
species = args.species # kegg speicies code, NA or related species if species not in KEGG
kofam = args.kofam # yes or no
indir = args.indir # directory with outputs from pull_data.sh
outdir = args.outdir # directory where outputs this will go
flybase = args.flybase #FB for Flybase annotations, NA for none
orthologs = args.orthologs # $3/orthofinder/Orthologues_"$noext"-cluster/"$noext"-cluster__v__dromel-cluster.tsv from pathannotator.sh script

pd.set_option('display.max_columns', None)
#################
#FUNCTIONS
#################
def combine_series_lists(s1, s2):
    combined = []
    for val1, val2 in zip(s1, s2):
        list1 = val1 if isinstance(val1, list) else []
        list2 = val2 if isinstance(val2, list) else []
        combined.append(list1 + list2)
    return pd.Series(combined)
#####################
####################

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
    ncbi_ver_spec_ko_pathway_pathname.to_csv(f"{outdir}/{species}_KEGG_ref.tsv", sep='\t', index=False)
    keggref = ncbi_ver_spec_ko_pathway_pathname[["Input_protein_ID","KEGG_ref_pathway"]]
    keggref.columns = ['Input_protein_ID','pathway']
    keggref.loc[:, 'pathway'] = 'KEGG:' + keggref['pathway'].astype(str)
    acckeggref = keggref.groupby('Input_protein_ID')['pathway'].agg(list).reset_index()
    annkeggref = keggref.groupby('pathway')['Input_protein_ID'].agg(list).reset_index()
#MERGE DATAFRAMES INTO ONE FOR { species } PATHWAYS
    ncbi_ver_spec_ko_specpath = pd.merge(ncbi_ver_spec_ko, spec_pathway, on='KEGG_genes_ID', how='inner')
    ncbi_ver_spec_ko_specpath_specpathname = pd.merge(ncbi_ver_spec_ko_specpath, list_pathway_spec, on=f"KEGG_{species}_pathway", how='left')
    ncbi_ver_spec_ko_specpath_specpathname.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
    ncbi_ver_spec_ko_specpath_specpathname.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
    ncbi_ver_spec_ko_specpath_specpathname.to_csv(f"{outdir}/{species}_KEGG_species.tsv", sep='\t', index=False)
    keggspec = ncbi_ver_spec_ko_specpath_specpathname[["Input_protein_ID",f"KEGG_{species}_pathway"]]
    keggspec.columns = ['Input_protein_ID','pathway']
    keggspec.loc[:, 'pathway'] = 'KEGG:' + keggspec['pathway'].astype(str)
    acckeggspec = keggspec.groupby('Input_protein_ID')['pathway'].agg(list).reset_index()
    annkeggspec = keggspec.groupby('pathway')['Input_protein_ID'].agg(list).reset_index()
#MERGE FOR AGGREGATED OUTPUTS
    acckegg = pd.merge(acckeggref, acckeggspec, on='Input_protein_ID', how='outer')
    annkegg = pd.merge(annkeggref, annkeggspec, on='pathway', how='outer')
#    acckegg['pathway'] = acckegg['pathway_x'].astype(str) + acckegg['pathway_y'].astype(str)
#    acckegg['pathway'] = acckegg['pathway_x'].fillna('') + acckegg['pathway_y'].fillna('')
#    annkegg['Input_protein_ID'] = annkegg['Input_protein_ID_x'].astype(str) + annkegg['Input_protein_ID_y'].astype(str)
#    annkegg['Input_protein_ID'] = annkegg['Input_protein_ID_x'].fillna('') + annkegg['Input_protein_ID_y'].fillna('')
    acckegg['pathway'] = combine_series_lists(acckegg['pathway_x'], acckegg['pathway_y'])
    annkegg['Input_protein_ID'] = combine_series_lists(annkegg['Input_protein_ID_x'], annkegg['Input_protein_ID_y'])
    acckegg.to_csv(f"{outdir}/{species}_acc_pathways.tsv", sep='\t', index=False)
    annkegg.to_csv(f"{outdir}/{species}_pathways_acc.tsv", sep='\t', index=False)
#    acckegg.drop(['pathway_x', 'pathway_y'], axis=1, inplace=True)
#    annkegg.drop(['Input_protein_ID_x', 'Input_protein_ID_y'], axis=1, inplace=True)
#    acckegg.to_csv(f"{outdir}/{species}_acc_pathways.tsv", sep='\t', index=False)
#    annkegg.to_csv(f"{outdir}/{species}_pathways_acc.tsv", sep='\t', index=False)
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
        ncbi_ver_spec_ko['KEGG_genes_ID'] = ncbi_ver_spec_ko['KEGG_genes_ID'].str.replace('Dmel_', '')
        fbgn_CG_path_ncbi_ver_spec_ko = pd.merge(ncbi_ver_spec_ko, fbgn_CG_path, on='KEGG_genes_ID', how='inner')
        fbgn_CG_path_ncbi_ver_spec_ko.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko = fbgn_CG_path_ncbi_ver_spec_ko[["KEGG_genes_ID","Input_protein_ID","KEGG_KO","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_CG_path_ncbi_ver_spec_ko.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko.to_csv(f"{outdir}/{species}_flybase.tsv", sep='\t', index=False)
    #MERGE FOR AGGREGATED OUTPUTS
        fbpath = fbgn_CG_path_ncbi_ver_spec_ko[["Input_protein_ID","Flybase_pathway_ID"]]
        fbpath.columns = ['Input_protein_ID','pathway']
        fbpath.loc[:, 'pathway'] = 'Flybase:' + fbpath['pathway'].astype(str)
        accfbpath = fbpath.groupby('Input_protein_ID')['pathway'].agg(list).reset_index()
        annfbpath = fbpath.groupby('pathway')['Input_protein_ID'].agg(list).reset_index()
        acckeggfb = pd.merge(acckegg, accfbpath, on='Input_protein_ID', how='outer')
        annkeggfb = pd.merge(annkegg, annfbpath, on='pathway', how='outer')
        acckeggfb['pathway'] = acckeggfb['pathway_x'] + acckeggfb['pathway_y']
        annkeggfb['Input_protein_ID'] = annkeggfb['Input_protein_ID_x'] + annkeggfb['Input_protein_ID_y']
        acckeggfb.drop(['pathway_x', 'pathway_Y'], axis=1, inplace=True)
        annkeggfb.drop(['Input_protein_ID_x', 'Input_protein_ID_Y'], axis=1, inplace=True)
        acckeggfb.to_csv(f"{outdir}/{species}_acc_pathways.tsv", sep='\t', index=False)
        annkeggfb.to_csv(f"{outdir}/{species}_pathways_acc.tsv", sep='\t', index=False)
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
    #SPLIT AND EXPLODE TO GET LISTS IN BOTH COLUMNS
        fbpp_ortho["Flybase_protein_ID"] = fbpp_ortho["Flybase_protein_ID"].str.split(", ")
        fbpp_ortho = fbpp_ortho.explode("Flybase_protein_ID")
        fbpp_ortho["Input_protein_ID"] = fbpp_ortho["Input_protein_ID"].str.split(", ")
        fbpp_ortho = fbpp_ortho.explode("Input_protein_ID")
        fbpp_ortho = fbpp_ortho.sort_values(by=['Flybase_protein_ID', 'Input_protein_ID'])
    #MERGE AND OUTPUT TO FILE
        fbgn_fbpp_ortho = pd.merge(fbgn_fbpp, fbpp_ortho, on='Flybase_protein_ID', how='inner')
        fbgn_fbpp_ortho_path = pd.merge(fbgn_fbpp_ortho, fbgn_path, on='Flybase_gene', how='inner')
        fbgn_fbpp_ortho_path_CG = pd.merge(fbgn_fbpp_ortho_path, fbgn_CG, on='Flybase_gene', how='inner')
        fbgn_fbpp_ortho_path_CG.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG = fbgn_fbpp_ortho_path_CG[["KEGG_genes_ID","Input_protein_ID","Flybase_protein_ID","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_fbpp_ortho_path_CG.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
        fbgn_fbpp_ortho_path_CG.to_csv(f"{outdir}/Orthofinder_flybase.tsv", sep='\t', index=False)
    #MERGE FOR AGGREGATED OUTPUTS
        fbpath = fbgn_fbpp_ortho_path_CG[["Input_protein_ID","Flybase_pathway_ID"]]
        fbpath.columns = ['Input_protein_ID','pathway']
        fbpath.loc[:, 'pathway'] = 'Flybase:' + fbpath['pathway'].astype(str)
        accfbpath = fbpath.groupby('Input_protein_ID')['pathway'].agg(list).reset_index()
        annfbpath = fbpath.groupby('pathway')['Input_protein_ID'].agg(list).reset_index()
        acckeggfb = pd.merge(acckegg, accfbpath, on='Input_protein_ID', how='outer')
        annkeggfb = pd.merge(annkegg, annfbpath, on='pathway', how='outer')
        acckeggfb['pathway'] = acckeggfb['pathway_x'] + acckeggfb['pathway_y']
        annkeggfb['Input_protein_ID'] = annkeggfb['Input_protein_ID_x'] + annkeggfb['Input_protein_ID_y']
        acckeggfb.drop(['pathway_x', 'pathway_Y'], axis=1, inplace=True)
        annkeggfb.drop(['Input_protein_ID_x', 'Input_protein_ID_Y'], axis=1, inplace=True)
        acckeggfb.to_csv(f"{outdir}/{species}_acc_pathways.tsv", sep='\t', index=False)
        annkeggfb.to_csv(f"{outdir}/{species}_pathways_acc.tsv", sep='\t', index=False)
        annkeggfb.to_csv(f"{outdir}/{species}_pathways_acc.tsv", sep='\t', index=False)
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
    ncbi_ver_ko_pathway_pathname.insert(0, 'KEGG_genes_ID','NA',allow_duplicates=True)
    ncbi_ver_ko_pathway_pathname = ncbi_ver_ko_pathway_pathname[['KEGG_genes_ID', 'Input_protein_ID','KEGG_KO', 'KEGG_ref_pathway', 'KEGG_ref_pathway_name']]
    ncbi_ver_ko_pathway_pathname.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
    ncbi_ver_ko_pathway_pathname.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
    ncbi_ver_ko_pathway_pathname.to_csv(f"{outdir}/NA_KEGG_ref.tsv", sep='\t', index=False)
#MERGE FOR AGGREGATED OUTPUTS
    keggref = ncbi_ver_spec_ko_pathway_pathname[["Input_protein_ID","KEGG_ref_pathway"]]
    keggref.columns = ['Input_protein_ID','pathway']
    keggref.loc[:, 'pathway'] = 'KEGG:' + keggref['pathway'].astype(str)
    acckeggref = keggref.groupby('Input_protein_ID')['pathway'].agg(list).reset_index()
    annkeggref = keggref.groupby('pathway')['Input_protein_ID'].agg(list).reset_index()
    acckeggref.to_csv(f"{outdir}/{species}_acc_pathways.tsv", sep='\t', index=False)
    annkeggref.to_csv(f"{outdir}/{species}_pathways_acc.tsv", sep='\t', index=False)
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
    #SPLIT AND EXPLODE TO GET LISTS IN BOTH COLUMNS
        fbpp_ortho["Flybase_protein_ID"] = fbpp_ortho["Flybase_protein_ID"].str.split(", ")
        fbpp_ortho = fbpp_ortho.explode("Flybase_protein_ID")
        fbpp_ortho["Input_protein_ID"] = fbpp_ortho["Input_protein_ID"].str.split(", ")
        fbpp_ortho = fbpp_ortho.explode("Input_protein_ID")
        fbpp_ortho = fbpp_ortho.sort_values(by=['Flybase_protein_ID', 'Input_protein_ID'])
    #MERGE AND OUTPUT TO FILE
        fbgn_fbpp_ortho = pd.merge(fbgn_fbpp, fbpp_ortho, on='Flybase_protein_ID', how='inner')
        fbgn_fbpp_ortho_path = pd.merge(fbgn_fbpp_ortho, fbgn_path, on='Flybase_gene', how='inner')
        fbgn_fbpp_ortho_path_CG = pd.merge(fbgn_fbpp_ortho_path, fbgn_CG, on='Flybase_gene', how='inner')
        fbgn_fbpp_ortho_path_CG.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG = fbgn_fbpp_ortho_path_CG[["KEGG_genes_ID","Input_protein_ID","Flybase_protein_ID","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_fbpp_ortho_path_CG.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
        fbgn_fbpp_ortho_path_CG.to_csv(f"{outdir}/Orthofinder_flybase.tsv", sep='\t', index=False)
    #MERGE FOR AGGREGATED OUTPUTS
        fbpath = fbgn_fbpp_ortho_path_CG[["Input_protein_ID","Flybase_pathway_ID"]]
        fbpath.columns = ['Input_protein_ID','pathway']
        fbpath.loc[:, 'pathway'] = 'Flybase:' + fbpath['pathway'].astype(str)
        accfbpath = fbpath.groupby('Input_protein_ID')['pathway'].agg(list).reset_index()
        annfbpath = fbpath.groupby('pathway')['Input_protein_ID'].agg(list).reset_index()
        acckeggfb = pd.merge(acckeggref, accfbpath, on=['Input_protein_ID', 'pathway'], how='outer')
        annkeggfb = pd.merge(annkeggref, annfbpath, on=['Input_protein_ID', 'pathway'], how='outer')
        acckeggfb.to_csv(f"{outdir}/{species}_acc_pathways.tsv", sep='\t', index=False)
        annkeggfb.to_csv(f"{outdir}/{species}_pathways_acc.tsv", sep='\t', index=False)
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
    print (ncbi_ver_ko)
    ncbi_ver_ko_pathway = pd.merge(ncbi_ver_ko, ko_pathway, on='KEGG_KO', how='inner')
    ncbi_ver_ko_pathway_pathname = pd.merge(ncbi_ver_ko_pathway, pathway, on='KEGG_ref_pathway', how='left')
    ncbi_ver_ko_pathway_pathname = ncbi_ver_ko_pathway_pathname[["Input_protein_ID","KEGG_KO","KEGG_ref_pathway","KEGG_ref_pathway_name"]]
    print (ncbi_ver_ko_pathway_pathname)
    ncbi_ver_ko_pathway_pathname.insert(0, 'KEGG_genes_ID','NA',allow_duplicates=True)
    print (ncbi_ver_ko_pathway_pathname)
    ncbi_ver_ko_pathway_pathname.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
    ncbi_ver_ko_pathway_pathname.rename(columns={"Input_protein_ID_version": "Input_protein_ID"}, inplace=True)
    ncbi_ver_ko_pathway_pathname.to_csv(f"{outdir}/{species}_KEGG_ref.tsv", sep='\t', index=False)
    keggref = ncbi_ver_ko_pathway_pathname[["Input_protein_ID", "KEGG_ref_pathway"]]
    keggref.columns = ['Input_protein_ID','pathway']
    keggref.loc[:, 'pathway'] = 'KEGG:' + keggref['pathway'].astype(str)
    acckeggref = keggref.groupby('Input_protein_ID')['pathway'].agg(list).reset_index()
    annkeggref = keggref.groupby('pathway')['Input_protein_ID'].agg(list).reset_index()
#MERGE DATAFRAMES INTO ONE FOR { species } PATHWAYS
    ncbi_ver_spec_ko = pd.merge(ncbi_ver_ko, spec_ko, on='KEGG_KO', how='inner')
    ncbi_ver_spec_ko_specpath = pd.merge(ncbi_ver_spec_ko, spec_pathway, on='KEGG_genes_ID', how='inner')
    ncbi_ver_spec_ko_specpath_specpathname = pd.merge(ncbi_ver_spec_ko_specpath, list_pathway_spec, on=f"KEGG_{species}_pathway", how='left')
    ncbi_ver_spec_ko_specpath_specpathname = ncbi_ver_spec_ko_specpath_specpathname[["KEGG_genes_ID","Input_protein_ID","KEGG_KO",f"KEGG_{species}_pathway",f"KEGG_{species}_pathway_name"]]
    ncbi_ver_spec_ko_specpath_specpathname.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
    ncbi_ver_spec_ko_specpath_specpathname.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
    ncbi_ver_spec_ko_specpath_specpathname.to_csv(f"{outdir}/{species}_KEGG_species.tsv", sep='\t', index=False)
    keggspec = ncbi_ver_spec_ko_specpath_specpathname[["Input_protein_ID","KEGG_{species}_pathway"]]
    keggspec.columns = ['Input_protein_ID','pathway']
    keggspec.loc[:, 'pathway'] = 'KEGG:' + keggspec['pathway'].astype(str)
    acckeggspec = keggref.groupby('Input_protein_ID')['pathway'].agg(list).reset_index()
    annkeggspec = keggref.groupby('pathway')['Input_protein_ID'].agg(list).reset_index()
#MERGE FOR AGGREGATED OUTPUTS
    acckegg = pd.merge(acckeggref, acckeggspec, on=['Input_protein_ID', 'pathway'], how='outer')
    annkegg = pd.merge(annkeggref, annkeggspec, on=['Input_protein_ID', 'pathway'], how='outer')
    acckegg.to_csv(f"{outdir}/{species}_acc_pathways.tsv", sep='\t', index=False)
    annkegg.to_csv(f"{outdir}/{species}_pathways_acc.tsv", sep='\t', index=False)
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
        ncbi_ver_spec_ko['KEGG_genes_ID'] = ncbi_ver_spec_ko['KEGG_genes_ID'].str.replace('Dmel_', '')
        fbgn_CG_path_ncbi_ver_spec_ko = pd.merge(ncbi_ver_spec_ko, fbgn_CG_path, on='KEGG_genes_ID', how='inner')
        fbgn_CG_path_ncbi_ver_spec_ko.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko = fbgn_CG_path_ncbi_ver_spec_ko[["KEGG_genes_ID","Input_protein_ID","KEGG_KO","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_CG_path_ncbi_ver_spec_ko.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
        fbgn_CG_path_ncbi_ver_spec_ko.to_csv(f"{outdir}/{species}_flybase.tsv", sep='\t', index=False)
    #MERBE FOR AGGREGATED OUTPUTS
        fbpath = fbgn_CG_path_ncbi_ver_spec_ko[["Input_protein_ID","Flybase_pathway_ID"]]
        fbpath.columns = ['Input_protein_ID','pathway']
        fbpath.loc[:, 'pathway'] = 'Flybase:' + fbpath['pathway'].astype(str)
        accfbpath = fbpath.groupby('Input_protein_ID')['pathway'].agg(list).reset_index()
        annfbpath = fbpath.groupby('pathway')['Input_protein_ID'].agg(list).reset_index()
        acckeggfb = pd.merge(acckegg, accfbpath, on=['Input_protein_ID', 'pathway'], how='outer')
        annkeggfb = pd.merge(annkegg, annfbpath, on=['Input_protein_ID', 'pathway'], how='outer')
        acckeggfb.to_csv(f"{outdir}/{species}_acc_pathways.tsv", sep='\t', index=False)
        annkeggfb.to_csv(f"{outdir}/{species}_pathways_acc.tsv", sep='\t', index=False)
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
    #SPLIT AND EXPLODE TO GET LISTS IN BOTH COLUMNS
        fbpp_ortho["Flybase_protein_ID"] = fbpp_ortho["Flybase_protein_ID"].str.split(", ")
        fbpp_ortho = fbpp_ortho.explode("Flybase_protein_ID")
        fbpp_ortho["Input_protein_ID"] = fbpp_ortho["Input_protein_ID"].str.split(", ")
        fbpp_ortho = fbpp_ortho.explode("Input_protein_ID")
        fbpp_ortho = fbpp_ortho.sort_values(by=['Flybase_protein_ID', 'Input_protein_ID'])
    #MERGE AND OUTPUT TO FILE
        fbgn_fbpp_ortho = pd.merge(fbgn_fbpp, fbpp_ortho, on='Flybase_protein_ID', how='inner')
        fbgn_fbpp_ortho_path = pd.merge(fbgn_fbpp_ortho, fbgn_path, on='Flybase_gene', how='inner')
        fbgn_fbpp_ortho_path_CG = pd.merge(fbgn_fbpp_ortho_path, fbgn_CG, on='Flybase_gene', how='inner')
        fbgn_fbpp_ortho_path_CG.drop('Flybase_gene', axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG = fbgn_fbpp_ortho_path_CG[["KEGG_genes_ID","Input_protein_ID","Flybase_protein_ID","Flybase_pathway_ID","Flybase_pathway_name"]]
        fbgn_fbpp_ortho_path_CG.drop(['Input_protein_ID', 'KEGG_genes_ID'], axis=1, inplace=True)
        fbgn_fbpp_ortho_path_CG.rename(columns={'Input_protein_ID_version': 'Input_protein_ID'}, inplace=True)
        fbgn_fbpp_ortho_path_CG.to_csv(f"{outdir}/Orthofinder_flybase.tsv", sep='\t', index=False)
    #MERGE FOR AGGREGATED OUTPUT
        fbpath = fbgn_fbpp_ortho_path_CG[["Input_protein_ID","Flybase_pathway_ID"]]
        fbpath.columns = ['Input_protein_ID','pathway']
        fbpath.loc[:, 'pathway'] = 'Flybase:' + fbpath['pathway'].astype(str)
        accfbpath = fbpath.groupby('Input_protein_ID')['pathway'].agg(list).reset_index()
        annfbpath = fbpath.groupby('pathway')['Input_protein_ID'].agg(list).reset_index()
        acckeggfb = pd.merge(acckegg, accfbpath, on=['Input_protein_ID', 'pathway'], how='outer')
        annkeggfb = pd.merge(annkegg, annfbpath, on=['Input_protein_ID', 'pathway'], how='outer')
        acckeggfb.to_csv(f"{outdir}/{species}_acc_pathways.tsv", sep='\t', index=False)
        annkeggfb.to_csv(f"{outdir}/{species}_pathways_acc.tsv", sep='\t', index=False)
    else:
        print("You have not requested Flybase annotations.")
else:
    print("Not an acceptable combination of arguments.")

