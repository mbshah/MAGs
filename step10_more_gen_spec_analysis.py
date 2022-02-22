import random
#import subprocess
import sys
import statsmodels.stats.multitest
import pandas as pd
from scipy import stats
#import numpy as np
#import math
import json

csv_file_path = "../WMG/Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/tax_annot.csv"
#x=subprocess.check_output(["Rscript","/home/hb0358/PycharmProjects/mbs_general/MAGs/Step009-2_more_statistics.R"])
kraken_profile_folder = "../WMG/Kraken2_NT_Scaffolds_profiles2/"
# kraken_profile_folder = "/home/manan/disk2/UDE/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/"
kegg_db_folder = kraken_profile_folder + "../../ancilary/kegg/"
kegg_db_file = kegg_db_folder + "kegg.db"
kegg_database = json.load(open(kegg_db_file))  # Genes:name,def,pathways; Pathways:name,classes


def kegg_info_accumulate(genes_list,org):
    pathway_count={}
    class_count={}
    for kid in genes_list:
        kid_pathways=kegg_database["Genes"][kid]["pathways"]
        for pathway in kid_pathways:
            if pathway in pathway_count: pathway_count[pathway]=pathway_count[pathway]+1
            else: pathway_count[pathway]=1
            if not kegg_database["Pathways"][pathway]["class"].strip() =="":
                pway_class=kegg_database["Pathways"][pathway]["class"].strip().split(";")
                for kclass in pway_class:
                    kclass="kc_"+kclass.strip().replace(" ","_").replace(",","")
                    if kclass in class_count: class_count[kclass]=class_count[kclass]+1
                    else:class_count[kclass]=1
    return pathway_count, class_count


def remake_gsk(genes_list,gsk_table):
    new_gsk=pd.DataFrame()
    info_cols=gsk_table[[col for col in gsk_table if not col.startswith(("kc","ko","K"))]]
    for org,entry in gsk_table.iterrows():
        pgenes = []
        for gene in genes_list:
            if gsk_table.loc[org,gene]==1:
                pgenes.append(gene)
                new_gsk.loc[org, gene] = 1
            else:
                new_gsk.loc[org, gene] = 0
        pathway_count,class_count=kegg_info_accumulate(pgenes,org)
        for pathway in pathway_count:
            if not pathway in new_gsk.keys():new_gsk.loc[:,pathway]=0
            new_gsk.loc[org,pathway]=pathway_count[pathway]
        if len(class_count)>1:
            for kclass in class_count:
                if not kclass in new_gsk.keys():new_gsk.loc[:,kclass]=0
                new_gsk.loc[org,kclass]=class_count[kclass]
    ret_gsk=pd.concat([new_gsk,info_cols],axis=1).fillna(0)
    return ret_gsk


def table_processing():
    #table read
    gs_table = pd.read_csv(csv_file_path, sep=",", header=0, index_col=1)

    #table without kc and ko cols
    info_table=gs_table[[col for col in gs_table if not col.startswith(("kc", "ko","K"))]]
    #print(info_table.shape)
    gsnk_table = gs_table[[col for col in gs_table if col.startswith("K")]]
    #print(gsnk_table.shape)
    gsnk_table = gsnk_table.loc[gsnk_table.sum(axis=1) > 0, gsnk_table.sum(axis=0) > 0]
    #print(gsnk_table.shape)
    gsnk_table=pd.concat([gsnk_table,info_table],axis=1,join="outer",sort=False).fillna(0)
    #print(gsnk_table.shape)


    # individual tables for genralists and specialists:

    g_table = gsnk_table[gsnk_table["cl"] == 1]
    s_table = gsnk_table[gsnk_table["cl"] == 2]
    g_len_old=len([col for col in g_table if col.startswith('K')])
    s_len_old=len([col for col in s_table if col.startswith('K')])
    print(f"Old no of genes for Generalists: {g_len_old}\tSpecialists: {s_len_old}", file=sys.stderr)

    #removing non all zero genes from g and s tables
    # lengths of g_table and s_table at the end of this section is all genes present only in respective category
    g_table = g_table[[col for col in g_table if col.startswith("K")]]
    g_table = g_table.loc[(g_table.sum(axis=1)) > 0, g_table.sum(axis=0) > 0]
    s_table = s_table[[col for col in s_table if col.startswith("K")]]
    s_table = s_table.loc[(s_table.sum(axis=1)) > 0, s_table.sum(axis=0) > 0]
    g_len_new = len([col for col in g_table if col.startswith('K')])
    s_len_new = len([col for col in s_table if col.startswith('K')])
    print(f"New no of genes for Generalists: {g_len_new}\tSpecialists: {s_len_new}", file=sys.stderr)

    only_in_g=[gene for gene in [col for col in gsnk_table if col.startswith('K')] if gene not in [col for col in s_table if col.startswith("K")]]
    og_table=remake_gsk(only_in_g,gsnk_table)
    og_table.to_csv(
        "/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/tax_annot_og.csv")
    only_in_s=[gene for gene in [col for col in gsnk_table if col.startswith('K')] if gene not in [col for col in g_table if col.startswith("K")]]
    os_table = remake_gsk(only_in_s, gsnk_table)
    os_table.to_csv(
        "/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/tax_annot_os.csv")
    poino=only_in_s+only_in_g
    present_in_both= [gene for gene in [col for col in gsnk_table if col.startswith('K')] if gene not in poino]
    print (f"So {len(only_in_g)} genes present only in Generalist and {len(only_in_s)} genes are only present in Specialists",file=sys.stderr)
    print( f"Genes present in both are = {len(present_in_both)}",file=sys.stderr)

    #select KiDs randomly from g and s tables
    sample_size=4000
    random.seed(250)
    rand_KID1 = random.sample(present_in_both,sample_size)
    rand_KIDs=[col for col in s_table if col in rand_KID1]
    rand_KIDg=[col for col in g_table if col in rand_KID1]
    #rand_KID1= list(set(rand_KIDg+rand_KIDs))
    print(f"Randomly Selected {len(rand_KIDg)} Gen and {len(rand_KIDs)} spec genes and combined to get {len(rand_KID1)} total genes",file=sys.stderr)

    #calculating pathway and class information for randomly selected KiDs
    print(f"Creating annotation table for randomly selected Kegg IDs", file=sys.stderr)
    ngs_table=remake_gsk(rand_KID1,gsnk_table)
    #pd.set_option('display.max_columns', None)
    ngs_table=ngs_table.sort_index(axis=1)
    ngs_table.to_csv("/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/tax_annot2.csv")
    #print(ngs_table.sort_index(axis=1).to_csv())
    print(f"Variance Testing in progress by Fligner Test", file=sys.stderr)

    #creating g_table and s_table with all information
    g_table = ngs_table[gsnk_table["cl"] == 1]
    s_table = ngs_table[gsnk_table["cl"] == 2]

    #CLASS test
    kc_table = {}
    ko_table = {}
    kid_table = {}
    info_table={}
    info_params=("GC", )
    for col in ngs_table:
        #print (col)
        if col.startswith("kc"):
            y = stats.fligner(g_table[col], s_table[col])
            kc_table[col] = [abs(float(y[0])), y[1]]
        if col.startswith("ko"):
            y = stats.fligner(g_table[col], s_table[col])
            ko_table[col] = [abs(float(y[0])), y[1]]
        if col.startswith("K"):
            y = stats.fligner(g_table[col], s_table[col])
            kid_table[col] = [abs(float(y[0])), y[1]]
        #if col.startswith("total_Len_Scaffolds"):
        #    y= stats.fligner(g_table[col], s_table[col])
        #    print(y,file=sys.stderr)
    kct = pd.DataFrame.from_dict(kc_table, orient="index", columns=["s", "p"])
    kct = kct.sort_values("p")
    kct = kct.dropna(thresh=2)
    kot = pd.DataFrame.from_dict(ko_table, orient="index", columns=["s", "p"])
    kot = kot.sort_values("p")
    kot = kot.dropna(thresh=2)
    kidt = pd.DataFrame.from_dict(kid_table, orient="index", columns=["s", "p"])
    kidt = kidt.sort_values("p")
    kidt = kidt.dropna(thresh=2)
    #print(len(kidt))
    print(f"Correcting Multiple Testing Errors", file=sys.stderr)
    new_kcp = list(statsmodels.stats.multitest.multipletests(list(kct['p']),method="fdr_bh",is_sorted=True)[1])
    new_kop = list(statsmodels.stats.multitest.multipletests(list(kot['p']), method="fdr_bh", is_sorted=True)[1])
    i=0
    for index, x in kct.iterrows():
        x['p']=new_kcp[i]
        i=i+1
        if x["p"] < 0.005:
            print(f"{index}\t{x['s']}\t{x['p']}\t{str(index).replace('^kc_','')}")
    i = 0
    for index, x in kot.iterrows():
        x['p'] = new_kop[i]
        i = i + 1
        if x["p"] < 0.005:
            print(f"{index}\t{x['s']}\t{x['p']}\t{str(kegg_database['Pathways'][index]['name']).replace(', ','_')}")
    #for index, x in kidt.iterrows():
    #    if x["p"] < 0.005:
    #        print(f"{index}\t{x['s']}\t{x['p']}\t{str(kegg_database['Genes'][index]['def']).replace(', ','_')}")
    print(f"Done...", file=sys.stderr)


table_processing()
