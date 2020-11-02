# steps to run:
# 1. genemarkS on cluster genomes
# 2. KEGG KOFamKoala
# 3. run checkm on clusters
# 3. irep
# 4. growthpred
import subprocess

import files as config
import os

infolder = config.outfolder + "clusters" + config.postfix + "_" + config.method + "/"
genemark_folder = "../tools/gms2_linux_64/"
growthpred_folder = "../tools/growthpred-v1.07/"
os.environ["GROWTHPRED_SHARE"] = growthpred_folder + "/shared/"
os.environ["GROWTHPRED_LIBEX"] = growthpred_folder + "/Programs/"
kofam_folder="../tools/KOfamKoala/kofam_scan-1.3.0/"
genomes_metadata = {}
sample_metadata = {}
bins_metadata = {}


def load_data():
    global genomes_metadata
    global sample_metadata
    clusters_file = infolder + "clusters.tsv"
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
        #os.system(command)
        command=f"head -n1 {kofam_raw_output} > {kofam_filtered_output}"
        print (command)
        #os.system(command)
        command=f'grep "^*" {kofam_raw_output} >>{kofam_filtered_output}'
        print (command)
        #os.system(command)


def run_checkm():
    checkm_in=infolder+"cluster_reps/"
    checkm_out=infolder+"checkm_stg2_out/"
    checkm_tsv=infolder+"checkm_reps.tsv"
    command=f"checkm lineage_wf -t 8 -x fasta --tab_table -q {checkm_in} {checkm_out} >{checkm_tsv}"
    print(command)
    os.system(command)


def gen_spec_listing():
    profile=""

load_data()
#geneprediction_generate()
#run_growthpred()
run_kofamkoala()
#run_checkm()
