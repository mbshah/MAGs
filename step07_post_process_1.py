# steps to run:
# 1. genemarkS on cluster genomes
# 2. KEGG KOFamKoala
# 3. run checkm on clusters
# 3. irep
# 4. growthpred
import collections
import csv
import json
import subprocess
import urllib.request

import pandas as pd
from pathlib import Path
import config as config
import os

infolder = config.outfolder + "clusters" + config.postfix + "_" + config.method + "/"
genemark_folder = "../tools/gms2_linux_64/"
growthpred_folder = "../tools/growthpred-v1.07/"
os.environ["GROWTHPRED_SHARE"] = growthpred_folder + "/shared/"
os.environ["GROWTHPRED_LIBEX"] = growthpred_folder + "/Programs/"
kofam_folder="../tools/KOfamKoala/kofam_scan-1.3.0/"
kegg_db_folder="../ancilary/kegg/"
genomes_metadata = {}
sample_metadata = {}
bins_metadata = {}


def load_data():
    global genomes_metadata
    global sample_metadata
    global bins_metadata
    clusters_file = infolder + "clusters_wphyp_sited.tsv"
    sample_metadata_file = config.sample_metadata
    bins_metadata_file = config.outfolder + "/all_sample_summary.tsv"
    with open(clusters_file, "r") as cf:
        header = next(cf)
        header_ar = header.split("\t")
        for line in cf:
            entry = line.split("\t")
            key = entry[0]
            genomes_metadata[key] = {}
            i = 0
            while i < len(header_ar):
                genomes_metadata[key][header_ar[i]] = entry[i]
                i = i + 1
    with open(sample_metadata_file, "r") as sf:
        header = next(sf)
        header_ar = header.split("\t")
        for line in sf:
            entry = line.split("\t")
            key = entry[1]
            sample_metadata[key] = {}
            i = 0
            while i < len(header_ar):
                sample_metadata[key][header_ar[i]] = entry[i]
                i = i + 1
    with open(bins_metadata_file, "r") as bf:
        header = next(bf)
        header_ar = header.split("\t")
        for line in bf:
            entry = line.split("\t")
            key = entry[0]
            bins_metadata[key] = {}
            i = 0
            while i < len(header_ar):
                bins_metadata[key][header_ar[i]] = entry[i]
                i = i + 1


def geneprediction_generate():
    for genome in genomes_metadata:
        genome_dir = infolder + genome + "/"
        genome_fasta = genome_dir + genome + ".fasta"
        gene_pred_fnn = genome_dir + genome + "_gms.fna"
        gene_pred_faa = genome_dir + genome + "_gms.faa"
        gms_out = genome_dir + genome + "_gms.lst"
        command = f"perl {genemark_folder}/gms2.pl --seq {genome_fasta} --output {gms_out} --fnn {gene_pred_fnn} --faa {gene_pred_faa} --genome-type bacteria"
        print(command)
        os.system(command)
        os.system("rm ./GMS2.mod")


def get_ogt(members):
    temp_sum = 0
    if str(members).endswith(";"): members = members[0:-1]
    members = members.replace('"', '')
    members2 = members.split(";")
    for member in members2:
        if member != "":
            sample = "_".join(member.split("_")[0:2])
            temp = sample_metadata[sample]["temperature"]
            temp_sum = temp_sum + float(temp)
    temp_avg = temp_sum / len(members2)
    return temp_avg


def run_growthpred():
    for genome in genomes_metadata:
        genome_dir = infolder + genome + "/"
        gene_pred_fnn = genome_dir + genome + "_gms.fna"
        code = 0
        gp_out = "gp_out"
        gp_out_final = genome_dir + genome + "_growthpred.result"
        ogt = get_ogt(genomes_metadata[genome]["Members"])
        command = f"python2 {growthpred_folder}/growthpred-v1.07.py -b -g {gene_pred_fnn} -o {gp_out} -c {code} -T {ogt} -s -S"
        print (command)
        os.system(command)
        print(command)
        command=f"mv gp_out.results {gp_out_final}"
        os.system(command)
        #command=f'grep "Predicted minimum generation time" {infolder}/*/*'
        out = (subprocess.Popen(["grep", "Predicted minimum generation time", f'{gp_out_final}'],
                                stdout=subprocess.PIPE).communicate()[0]).decode("utf-8")
        out = out.replace("\n", "").replace("Predicted minimum generation time:  ", "")
        print(f"predicted growth rate= {out}")


def run_kofamkoala():
    for genome in genomes_metadata:
        genome_dir=infolder + genome + "/"
        gene_pred_faa = genome_dir + genome + "_gms.faa"
        config_file=kofam_folder+"config.yml"
        kofam_raw_output=genome_dir+genome+"_ko_raw.tsv"
        kofam_filtered_output=genome_dir+genome+"_ko_filtered.tsv"
        kofam_exe=f"{kofam_folder}exec_annotation"
        command=f"{kofam_exe} -c {config_file} -f detail-tsv -o {kofam_raw_output} {gene_pred_faa}"
        print(command)
        os.system(command)
        command=f"head -n1 {kofam_raw_output} > {kofam_filtered_output}"
        print (command)
        os.system(command)
        command=f'grep "^*" {kofam_raw_output} >>{kofam_filtered_output}'
        print (command)
        os.system(command)


def run_checkm():
    checkm_in=infolder+"cluster_reps/"
    checkm_out=infolder+"checkm_stg2_out/"
    checkm_tsv=infolder+"checkm_reps.tsv"
    command=f"checkm lineage_wf -t 8 -x fasta --tab_table -q {checkm_in} {checkm_out} >{checkm_tsv}"
    print(command)
    os.system(command)


def find_genes():
    genes_of_interest_list=[]
    genes_of_interest={}
    genomes_genes={}
    with open(config.scg_list_file,"r") as scg_list_file:
        for line in scg_list_file:
            col_of_interest=3
            gene=line.split("\t")[col_of_interest-1].strip().replace('"',"")
            if not gene=="":
                genes_of_interest_list.append(gene)
                genes_of_interest[gene]=[]
    for genome in genomes_metadata:
        genomes_genes[genome]=[]
        genome_annotation_file=f"{infolder}/{genome}/{genome}_annotation.tsv"
        col_of_gene_name=8
        with open(genome_annotation_file,"r") as annot_file:
            header=next(annot_file)
            for line in annot_file:
                entry=line.strip().split("\t")
                if len(entry)>col_of_gene_name-1:
                    gene=entry[col_of_gene_name-1]
                    for igene in genes_of_interest_list:
                        if gene.__contains__(igene):
                            if not genome in genes_of_interest[igene]:genes_of_interest[igene].append(genome)
                            if not igene in genomes_genes[genome]:genomes_genes[genome].append(igene)
    #summarize_gene_dict(genomes_genes,genes_of_interest_list,0)
    summarize_gene_dict(genes_of_interest,genomes_metadata.keys(),300)


def summarize_gene_dict(dict,master_list,cutoff):
    newdict={}
    for key in dict:
        length=len(dict[key])
        missing_individuals=[x for x in master_list if x not in dict[key]]
        newdict[key]=int(length)
        #print (f"{length}\t{key}") if length > cutoff else print("",end="")
    ordered_dict={k: v for k, v in sorted(newdict.items(), key=lambda item: item[1])}
    for k, v in ordered_dict.items(): print(k, v, dict[k])


kegg_db_file=kegg_db_folder+ "kegg.db"
if Path(kegg_db_file).is_file():
    kegg_database=json.load(open(kegg_db_file))
else:
    kegg_database = {"Genes": {}, "Pathways": {}}


def ko_profiler():
    table= {}
    table_pathway={}
    table_subSystem={}
    header=[]
    uniq_classes=[]
    uniq_pathways=[]
    all_KIDs=[]
    #infolder = config.outfolder + "clusters" + config.postfix + "_" + config.method + "/"
    org_tsv_out = infolder + "clusters_wphyp_sited.tsv"
    with open(org_tsv_out,"r")as org_tsv:
        header=next(org_tsv).strip().split("\t")
        for line in org_tsv:
            entry=line.strip().split("\t")
            print(f"Identifying KOs/Pathways/SubSystem for {entry[0]}")
            table[entry[0]]={}
            table[entry[0]]["Cluster"]=entry[0]
            table_pathway[entry[0]] = {}
            table_pathway[entry[0]]["Cluster"] = entry[0]
            table_subSystem[entry[0]] = {}
            table_subSystem[entry[0]]["Cluster"] = entry[0]
            #for i in range(0,len(header)):
            #    table[entry[0]][header[i]]=entry[i]
            kofam_out=f"{infolder}/{entry[0]}/{entry[0]}_ko_filtered.tsv"
            number_of_KOs=int(os.popen(f"grep \"^*\" {kofam_out} |wc -l").read())
            org_kids=[]
            org_pathways=[]
            org_classes= []
            #print(table[entry[0]]["taxonomy_ID"])
            if os.path.getsize(kofam_out)==0:pass
            else:
                with open(kofam_out) as kofam_tsv:
                    next(kofam_tsv)
                    for kLine in kofam_tsv:
                        kEntry=kLine.split("\t")
                        kID=kEntry[2]
                        org_kids.append(kID)
            kIDs_c=dict(collections.Counter(org_kids))
            #pathways = getpathways(kIDs_c)
            #org_pathways=org_pathways+pathways
            #pathway_count=dict(collections.Counter(org_pathways))
            #org_classes=org_classes+getkclasses(pathway_count)
            #class_count=dict(collections.Counter(org_classes))
            #print(class_count)
            pathway_count,class_count=kegg_info_accumulate(list(kIDs_c.keys()))
            for keggID in kIDs_c:
                if keggID not in all_KIDs: all_KIDs.append(keggID)
                table[entry[0]][keggID]=1
            for kpathway in pathway_count:
                if kpathway not in uniq_pathways:uniq_pathways.append(kpathway)
                table_pathway[entry[0]][kpathway]=pathway_count[kpathway]
            for kclass in class_count:
                if kclass not in uniq_classes:uniq_classes.append(kclass)
                table_subSystem[entry[0]][kclass]=class_count[kclass]
            number_of_KOs=len(kIDs_c)
            table[entry[0]]["no_of_KOs"]=number_of_KOs

    #writing_KOs
    t_o=org_tsv_out.replace("_wphyp_sited.tsv","")+"_ko_summary.tsv"
    header2=header[0:1]
    if "no_of_KOs" not in header2: header2.append("no_of_KOs")
    for keggID in all_KIDs:
        if keggID not in header2:header2.append(keggID)
    x=open(t_o,"w")
    x.write("\t".join(header2))
    x.write("\n")
    for entry in table:
        line=[]
        for value_h in header2:
            if value_h in table[entry]:value=table[entry][value_h]
            else:value=0
            line.append(str(value))
        x.write("\t".join(line)+"\n")
    x.close()
    #writing pathways
    t2_o = org_tsv_out.replace("_wphyp_sited.tsv", "") + "_pathway_summary.tsv"
    header2 = header[0:1]
    for kpathway in uniq_pathways:
        if kpathway not in header2:header2.append(kpathway)
    x = open(t2_o, "w")
    x.write("\t".join(header2))
    x.write("\n")
    for entry in table_pathway:
        line = []
        for value_h in header2:
            if value_h in table_pathway[entry]:
                value = table_pathway[entry][value_h]
            else:
                value = 0
            line.append(str(value))
        x.write("\t".join(line) + "\n")
    x.close()
    #writing subSystems
    t3_o = org_tsv_out.replace("_wphyp_sited.tsv", "") + "_subSystem_summary.tsv"
    header2 = header[0:1]
    for kclass in uniq_classes:
        if kclass not in header2:header2.append(kclass)
    x = open(t3_o, "w")
    x.write("\t".join(header2))
    x.write("\n")
    for entry in table_subSystem:
        line = []
        for value_h in header2:
            if value_h in table_subSystem[entry]:
                value = table_subSystem[entry][value_h]
            else:
                value = 0
            line.append(str(value))
        x.write("\t".join(line) + "\n")
    x.close()

    #writing KEgg_DB for quicker access
    kegg_map = json.dumps(kegg_database)
    f = open(kegg_db_file, "w")
    f.write(kegg_map)
    f.close()


def kegg_info_accumulate(genes_list):
    pathway_count={}
    class_count={}
    for kid in genes_list:
        kid_pathways=getpathways2(kid)
        for pathway in kid_pathways:
            if pathway in pathway_count: pathway_count[pathway]=pathway_count[pathway]+1
            else: pathway_count[pathway]=1
            pway_class=getkclasses2(pathway)
            for kclass in pway_class:
                if kclass in class_count: class_count[kclass]=class_count[kclass]+1
                else:class_count[kclass]=1
    return pathway_count, class_count


def getpathways2(id):
    kFile = kegg_db_folder + id
    if Path(kFile).is_file():
        pass
    else:
        url = "http://rest.kegg.jp/get/" + id
        # print(f"retreiving {url} as {kFile}")
        urllib.request.urlretrieve(url, kFile)
    pathwayflag=0
    pathways=[]
    with open(kFile, "r") as file:
        kegg_database["Genes"][id]={}
        kegg_database["Genes"][id]["pathways"]=[]
        for line in file:
            if line[0].isupper():
                pathwayflag=0
            if pathwayflag==1:
                pathway=line.split()[0].strip().replace(" ","_").replace(",","_")
                pathways.append(pathway)
                kegg_database["Genes"][id]["pathways"].append(pathway)
            if line.startswith("PATHWAY"):
                pathwayflag=1
                pathway=line.split()[1]
                pathways.append(pathway)
                kegg_database["Genes"][id]["pathways"].append(pathway)
            if line.startswith("NAME"):
                kegg_database["Genes"][id]["name"]=" ".join(line.split(" ")[1:]).strip()
            if line.startswith("DEFINITION"):
                kegg_database["Genes"][id]["def"]=" ".join(line.split(" ")[1:]).strip()
    return pathways


def getkclasses2(koID):
    koClass = []
    koFile = kegg_db_folder + koID
    if Path(koFile).is_file():
        pass
    else:
        url = "http://rest.kegg.jp/get/" + koID
        # print(f"retreiving {url} as {kFile}")
        urllib.request.urlretrieve(url, koFile)
    with open(koFile, "r") as file:
        kegg_database["Pathways"][koID]={}
        kegg_database["Pathways"][koID]["class"]=""
        orthoflag=0
        kegg_database["Pathways"][koID]["orthos"]=[]
        for line in file:
            if orthoflag==1:
                if line.startswith(r"[A-Z]"): orthoflag=0
                else:
                    ortho=line.split()[0].strip().replace(" ","_").replace(",","_")
                    kegg_database["Pathways"][koID]["orthos"].append(ortho)
            if line.startswith("ORTHOLOGY"):
                orthoflag=1
                ortho=line.split()[1]
                kegg_database["Pathways"][koID]["orthos"].append(ortho)
            if line.startswith("NAME"):
                kegg_database["Pathways"][koID]["name"]=" ".join(line.split(" ")[1:]).strip()
            if line.startswith("CLASS"):
                kegg_database["Pathways"][koID]["class"]=list(filter(None, line.strip().split('  ')))[1]
                kClass = list(filter(None, line.split('  ')))[1].split(";")
                for i in kClass:
                    i1="kc_"+i.strip().replace(" ","_").replace(",","_")
                    koClass.append(i1)
    return koClass


load_data()
#geneprediction_generate()
#run_growthpred()
#run_kofamkoala()
#run_checkm()
find_genes()
#ko_profiler()
