import statistics
import subprocess
import config as config
from Bio import SeqIO
import os
import step03_tax_profiler as tp
import pandas as pd
from pathlib import Path
import shutil
import json
import time
import re
from bs4 import BeautifulSoup

drep_folder = config.outfolder + "dRep__gff/"
genomes_folder = config.outfolder + "good_bins/"

postfix = config.postfix
method = config.method
infolder = drep_folder
threads = 15
taxdumpdir = "../ancilary/new_taxdump"
kraken_path = "~/manan_WMG/tools/kraken2/kraken2"
kraken_db = "~/manan_WMG/ancilary/kraken_bact/"
siteMetadata = {}
siteClAbundance = {}
phyla_colors = {
    "Pseudomonadota": "#88CCEE",
    "Bacteroidota": "#44AA99",
    "Actinomycetota": "#332288",
    "Cyanobacteriota": "#AAAA00",
    "Verrucomicrobiota": "#EE8866",
    "Planctomycetota": "#77AADD",
    "Chloroflexota": "#117733",
    "Bdellovibrionota": "#DDCC77",
    "Deinococcota": "#999933",
    "Armatimonadota": "#CC6677",
    "Bdellovibrionota_C": "#882255",
    "Acidobacteriota": "#AA4499",
    "Firmicutes": "#DDDDDD",
    "Gemmatimonadota": "#BBCC33",
    "Myxococcota": "#D41159"
}


def fasta_info(fastafile):
    mem_fasta = open(fastafile, "r")
    cluster_len = 0
    gc_count = 0
    max_len = 0
    for record in SeqIO.parse(mem_fasta, "fasta"):
        gc_count = gc_count + record.seq.count("G") + record.seq.count("C")
        cluster_len = cluster_len + record.__len__()
        if max_len < record.__len__(): max_len = record.__len__()
    return gc_count, cluster_len, max_len


def summarize_reassembly():
    infile = infolder + "post_cluster_summary.csv"
    tmpfile = infolder + "post_cluster_data.csv"
    segregated_folder = infolder + "Genome_data"
    Path(segregated_folder).mkdir(parents=True, exist_ok=True)
    o = open(tmpfile, "w")
    with open(infile, "r") as file:
        header = next(file)
        header = header.strip() + ",reassembly_gc,reassembly_size_mb,reassembly_longestContig_kb,coding_perc\n"
        o.write(header)
        for line in file:
            entry = line.strip().split(",")
            clustername = entry[0]
            cluster_folder = segregated_folder + "/" + clustername + "/"
            Path(cluster_folder).mkdir(exist_ok=True)
            rep_fasta = entry[1]
            fastafile = infolder + "/dereplicated_genomes/" + rep_fasta
            prodigal_fna = infolder + "/data/prodigal/" + rep_fasta + ".fna"
            prodigal_faa = infolder + "/data/prodigal/" + rep_fasta + ".faa"
            prodigal_gff = infolder + "/data/prodigal/" + rep_fasta + ".gff"
            prodigal_fna_dest = cluster_folder + rep_fasta + ".fna"
            prodigal_faa_dest = cluster_folder + rep_fasta + ".faa"
            prodigal_gff_dest = cluster_folder + rep_fasta + ".gff"
            if not Path(prodigal_fna_dest).is_file(): shutil.copy2(prodigal_fna, prodigal_fna_dest)
            if not Path(prodigal_faa_dest).is_file(): shutil.copy2(prodigal_faa, prodigal_faa_dest)
            if not Path(prodigal_gff_dest).is_file(): shutil.copy2(prodigal_gff, prodigal_gff_dest)
            members = entry[16].strip().split(";")
            for member in members:
                if member == rep_fasta:
                    source_fasta = fastafile
                    dest_fasta = cluster_folder + rep_fasta
                else:
                    source_fasta = genomes_folder + member
                    dest_fasta = cluster_folder + member
                if not Path(dest_fasta).is_file(): shutil.copy2(source_fasta, dest_fasta)
            coding_perc = get_perc_coding(rep_fasta, "total")
            gc_count, cluster_len, max_len = fasta_info(fastafile)
            gc_perc = (gc_count / cluster_len) * 100
            entry.append(str(float("{:.2f}".format(gc_perc))))
            entry.append(str(float("{:.2f}".format(cluster_len / 1000000))))
            entry.append(str(float("{:.2f}".format(max_len / 1000))))
            entry.append(str(coding_perc))
            new_line = ",".join(entry)
            o.write(new_line + "\n")
    o.close()


def retax_kraken(infolder1=infolder):
    clfolder = infolder1
    cfile_in = clfolder + "clusters_tmp.tsv"
    with open(cfile_in, "r+") as c:
        header = next(c)
        for cluster_line in c:
            cluster_info = cluster_line.strip().split("\t")
            cluster = cluster_info[0]
            # if cluster in ("EUL_0","EUL_37"):
            cl_path = clfolder + cluster + "/"
            infasta = cl_path + cluster + ".fasta"
            kraken_out_report = cl_path + cluster + ".kreport.tsv"
            command = kraken_path + " --db " + kraken_db + " --threads " + str(
                threads) + " --report " + kraken_out_report + " " + infasta
            os.system(command)
            genus, perc, genus_id = indentify_genus(kraken_out_report, cluster)


def indentify_genus(report_file, cluster):
    max_perc = 0
    ret_gen = ""
    ret_gid = ""
    with open(report_file, "r") as report:
        for entry_line in report:
            entry = entry_line.split("\t")
            if entry[4] == "0": entry[4] = "1"
            if entry[4] in tp.nodes:
                pass
            else:
                entry[4] = tp.find_missing_tax(entry[4])
            if tp.nodes[entry[4]] == "genus":
                if float(entry[0]) > max_perc:
                    max_perc = float(entry[0])
                    ret_gen = tp.names[entry[4]]
                    ret_gid = entry[4]
    min_assignment = 10
    global count
    if float(max_perc) < min_assignment: count = count + 1;print(
        f"{count}: warning cluster {cluster} has been assigned genus {ret_gen} with less the {min_assignment}% assignment({max_perc})")
    return ret_gen, max_perc, ret_gid
    # print(entry[0]+"\t"+tp.nodes[str(entry[4])]+"\t"+tp.names[entry[4]])


def find_phyla(taxa):
    tax_id = taxa.strip()
    if taxa not in tp.names: tax_id = tp.find_missing_tax(taxa)
    tax_id = tax_id.strip()
    lin = tp.lineage_db[tax_id]
    lineage = lin.split()
    phyla = "Other"
    clas = "Other"
    for level in lineage:
        l_node = tp.nodes[level]
        if l_node == "phylum":
            phyla = tp.names[level]
        if l_node == "class":
            clas = tp.names[level]
    if phyla == "Other": print("----")
    retval = clas if phyla == "Pseudomonadota" else phyla
    return retval


def find_sites(bins_list):
    bin_list = bins_list.replace('"', '')[0:-1]
    bins = bin_list.split(";")
    site_list = []
    for bin in bins:
        site = "_".join(bin.strip().split("_")[0:2])
        site_list.append(site)
    ret_str = ";".join(site_list)
    return ret_str


def getSiteData(sites, parameter):
    if str(sites).endswith(""): sites = "".join(sites.rsplit("\"", 1))
    if str(sites).endswith(";"): sites = "".join(sites.rsplit(";", 1))
    sites_arr = sites.split(";")
    parameter_array = []
    for site in sites_arr:
        site = site.strip()
        parameter_value = siteMetadata[site][parameter]
        if parameter_value not in parameter_array: parameter_array.append(
            float(parameter_value) if not parameter.__contains__("level") else parameter_value)
    # if len(sites_arr) >1:print(f"{sites}\t{parameter_array}")
    # return str(";".join(parameter_array))
    if parameter.__contains__("level"):
        retval = "; ".join(parameter_array)
    else:
        retval = parameter_array[0] if len(parameter_array) < 2 else statistics.mean(parameter_array)
        # retval = "; ".join(parameter_array)
    return retval


def getWtSiteData(sites, parameter, cluster):
    if str(sites).endswith(""): sites = "".join(sites.rsplit("\"", 1))
    if str(sites).endswith(";"): sites = "".join(sites.rsplit(";", 1))
    sites_arr = sites.split(";")
    parameter_array = []
    total_count = 0
    for site in sites_arr:
        site = site.strip()
        parameter_value = float(siteMetadata[site][parameter]) * float(siteClAbundance[cluster][site])
        total_count = total_count + float(siteClAbundance[cluster][site])
        parameter_array.append(float(parameter_value))
    retval = sum(parameter_array) / total_count
    print(cluster,retval, sites_arr)
    return retval


def getSiteDataArray(sites, parameter):
    if str(sites).endswith(""): sites = "".join(sites.rsplit("\"", 1))
    if str(sites).endswith(";"): sites = "".join(sites.rsplit(";", 1))
    sites_arr = sites.split(";")
    parameter_array = []
    for site in sites_arr:
        site = site.strip()
        parameter_value = siteMetadata[site][parameter]
        parameter_array.append(str(parameter_value))
    # if len(sites_arr) >1:print(f"{sites}\t{parameter_array}")
    retval = str("; ".join(parameter_array))
    return retval


def load_site_metadata():
    seperator = "\t"
    with open(config.sample_metadata, "r") as md_file:
        header = next(md_file).strip().replace('"', '').split(seperator)
        for line in md_file:
            entry = line.strip().replace('"', '').split(seperator)
            site = entry[0]
            for i in range(1, len(entry)):
                if site not in siteMetadata: siteMetadata[site] = {}
                siteMetadata[site][header[i].strip()] = entry[i]
    with open(drep_folder + "cluster_abundance_profile_normalized_deseq2.tsv", "r") as cl_abundance:
        header_sites = next(cl_abundance)
        header_sites = header_sites.split()[1:]
        for line in cl_abundance:
            entry = line.split()
            cluster = entry[0]
            counts = entry[1:]
            for i in range(len(header_sites)):
                #print(f"{i}\t{header_sites[i]}\t{counts[i]}")
                if cluster not in siteClAbundance:
                    siteClAbundance[cluster] = {}
                siteClAbundance[cluster][header_sites[i]] = counts[i]


def get_perc_coding(Cluster, ret_type="total"):
    cluster_lstfile = drep_folder + "data/prodigal/" + Cluster + ".gff"  ##change
    combined_arr = []
    forward_arr = []
    reverse_arr = []
    t_sum = 0
    start = 3
    end = 4
    strand_col = 6
    with open(cluster_lstfile, "r") as gff:
        gff_header = next(gff)
        this_combined = []
        this_forward = []
        this_reverse = []
        for line in gff:
            line = line.strip()
            if line.startswith("# Sequence Data:"):
                combined_arr = [*combined_arr, *this_combined]
                forward_arr = [*forward_arr, *this_forward]
                reverse_arr = [*reverse_arr, *this_reverse]
                seq_data = line.split(";")
                seq_data[0] = seq_data[0].split(":")[1].strip()
                seq_len = int(seq_data[1].split("=")[1].strip())
                t_sum = t_sum + seq_len
                this_reverse = [0] * (seq_len)
                this_combined = [0] * (seq_len)
                this_forward = [0] * (seq_len)
            elif line.startswith("# Mod"):
                pass
            else:
                line_entry = line.split()
                ent_strand = line_entry[strand_col]

                if ent_strand.__contains__("+"):
                    ent_start = int(line_entry[start].replace(">", "").replace("<", ""))
                    ent_end = int(line_entry[end].replace(">", "").replace("<", ""))
                    for i in range(ent_start, ent_end):
                        # print(f"{i}\t{forward_arr}\t{ent_cid}\t{cl}")
                        this_forward[i] = 1
                        this_combined[i] = 1
                else:
                    ent_end = int(line_entry[end].replace(">", "").replace("<", ""))
                    ent_start = int(line_entry[start].replace(">", "").replace("<", ""))
                    for i in range(ent_start, ent_end):
                        # print(f"{i}\t{len(reverse_arr)}\t{ent_cid}\t{Cluster}")
                        this_reverse[i] = 1
                        this_combined[i] = 1
        combined_arr = [*combined_arr, *this_combined]
        forward_arr = [*forward_arr, *this_forward]
        reverse_arr = [*reverse_arr, *this_reverse]
    combined_perc = combined_arr.count(1) * 100 / len(combined_arr)
    forward_perc = forward_arr.count(1) * 100 / len(forward_arr)
    reverse_perc = reverse_arr.count(1) * 100 / len(reverse_arr)
    ret_val = forward_perc if ret_type == "forward" else reverse_perc if ret_type == "reverse" else combined_perc
    return ret_val


def addSiteMetaData(working_file=infolder + "post_cluster_data.csv"):
    # working_file = infolder + "post_cluster_data.csv"
    my_table = pd.read_csv(working_file)
    if "site_list" not in my_table.columns: my_table['site_list'] = my_table.members_list.apply(lambda x: find_sites(x))
    if "site_pH" not in my_table.columns: my_table['site_pH'] = my_table.site_list.apply(
        lambda x: getSiteData(x, "pH1"))
    if "site_Elevation" not in my_table.columns: my_table['site_Elevation'] = my_table.site_list.apply(
        lambda x: getSiteData(x, "elevation"))
    if "site_Temperature" not in my_table.columns: my_table['site_Temperature'] = my_table.site_list.apply(
        lambda x: getSiteData(x, "temperature"))
    if "site_TP" not in my_table.columns: my_table['site_TP'] = my_table.site_list.apply(
        lambda x: getSiteData(x, "TP"))
    if "weighted_TP" not in my_table.columns: my_table['weighted_TP'] = my_table.apply(
        lambda x: getWtSiteData(x["site_list"], "TP", x["cluster_ID"]), axis=1)
    if "weighted_TN" not in my_table.columns: my_table['weighted_TN'] = my_table.apply(
        lambda x: getWtSiteData(x["site_list"], "TNs", x["cluster_ID"]), axis=1)
    if "weighted_DOC" not in my_table.columns: my_table['weighted_DOC'] = my_table.apply(
        lambda x: getWtSiteData(x["site_list"], "DOC", x["cluster_ID"]), axis=1)
    # if "site_pH_Array" not in my_table.columns: my_table['site_pH_Array'] = my_table.site_list.apply(
    #    lambda x: getSiteDataArray(x, "pH1"))
    # if "site_Elevation_Array" not in my_table.columns: my_table['site_Elevation_Array'] = my_table.site_list.apply(
    #    lambda x: getSiteDataArray(x, "elevation"))
    # if "site_Temperature_Array" not in my_table.columns: my_table['site_Temperature_Array'] = my_table.site_list.apply(
    #    lambda x: getSiteDataArray(x, "temperature"))
    # if "site_TP_Array" not in my_table.columns: my_table['site_TP_Array'] = my_table.site_list.apply(
    #    lambda x: getSiteDataArray(x, "TP"))
    my_table.to_csv(working_file, index=False)


def addgtdbtkInfo():
    working_file = infolder + "post_cluster_data.csv"
    magsTable = pd.read_csv(working_file)
    if "taxonomy" in magsTable.columns: magsTable = magsTable.drop(["taxonomy", "tax_confidence"], axis=1)
    gtdb_assignment_file = infolder + "gtdbtk_classify_out/gtdbtk.bac120.summary.tsv"
    gtdb_table = pd.read_csv(gtdb_assignment_file, sep="\t")
    if "gtdb_phyla" not in magsTable.columns: magsTable["gtdb_phyla"] = magsTable.Representative_Genome.apply(
        lambda x: getgtdbInfo(1, gtdb_table.loc[
            gtdb_table['user_genome'] == str(x).replace(".fasta", ""), "classification"].tolist()[0]))
    if "gtdb_genus" not in magsTable.columns: magsTable["gtdb_genus"] = magsTable.Representative_Genome.apply(
        lambda x: getgtdbInfo(5, gtdb_table.loc[
            gtdb_table['user_genome'] == str(x).replace(".fasta", ""), "classification"].tolist()[0]))
    if "gtdb_species" not in magsTable.columns: magsTable["gtdb_species"] = magsTable.Representative_Genome.apply(
        lambda x: getgtdbInfo(6, gtdb_table.loc[
            gtdb_table['user_genome'] == str(x).replace(".fasta", ""), "classification"].tolist()[0]))
    if "gtdb_taxonomy" not in magsTable.columns: magsTable["gtdb_taxonomy"] = magsTable.Representative_Genome.apply(
        lambda x: gtdb_table.loc[gtdb_table['user_genome'] == str(x).replace(".fasta", ""), "classification"].tolist()[
            0])
    if "phyla_colors" not in magsTable.columns: magsTable['phyla_colors'] = magsTable.gtdb_phyla.apply(
        lambda x: phyla_colors[x])

    magsTable.to_csv(working_file, index=False)


def getgtdbInfo(lev, tax_str):
    tax = str(tax_str).strip().split(";")
    taxa = tax[lev]
    taxa_refined = taxa.split("__")[-1]
    if taxa_refined == "": taxa_refined = "N/A"
    # print(str(tax_str),taxa_refined)
    return taxa_refined


def create_members_table():
    cluster_table = infolder + "post_cluster_data.csv"
    members_df = pd.DataFrame(
        columns=["cluster", "size_rep", "size_member", "site_list", "gtdb_phylum", "gtdb_genus", "gtdb_species", "abundance_m"])
    cluster_df = pd.read_csv(cluster_table, index_col="cluster_ID")
    for row in cluster_df.iterrows():
        members = row[1]["members_list"].split(";")
        for member in members:
            cluster = row[0]
            rep_size = row[1]["reassembly_size_mb"]
            rep_genus = row[1]["gtdb_genus"]
            rep_phyla = row[1]["gtdb_phyla"]
            rep_speceies = row[1]["gtdb_species"]
            mem_size = fasta_info(f"outfile/dRep__gff/Genome_data/{row[0]}/" + member)[1]
            mem_size = str(float("{:.2f}".format(mem_size / 1000000)))
            site = "_".join(member.strip().split("_")[0:2])
            abundance_mem=siteClAbundance[cluster][site]
            # print(f"{member}\t{rep_size}\t{mem_size}\t{site}")
            members_df.loc[member] = [cluster, rep_size, mem_size, site, rep_phyla, rep_genus, rep_speceies, abundance_mem]
    members_df.index.name = "bin"
    members_df.to_csv(infolder + "members_table.csv")
    addSiteMetaData(infolder + "members_table.csv")


def calculate_codon_usage(rg, cl):
    codonDict = {'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0, 'ACC': 0, 'ACG': 0,
                 'ACT': 0, 'AGA': 0, 'AGC': 0, 'AGG': 0, 'AGT': 0, 'ATA': 0, 'ATC': 0, 'ATG': 0, 'ATT': 0, 'CAA': 0,
                 'CAC': 0,
                 'CAG': 0, 'CAT': 0, 'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, 'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0,
                 'CTA': 0,
                 'CTC': 0, 'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0, 'GAT': 0, 'GCA': 0, 'GCC': 0, 'GCG': 0,
                 'GCT': 0,
                 'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, 'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0, 'TAA': 0, 'TAC': 0,
                 'TAG': 0,
                 'TAT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0, 'TTA': 0,
                 'TTC': 0,
                 'TTG': 0, 'TTT': 0}
    rg_fna_loc = infolder + "Genome_data/" + cl + "/" + rg + ".fna"
    for record in SeqIO.parse(rg_fna_loc, "fasta"):
        for codon_no in range(1, int(len(record) / 3) + 1):
            codon_start = (codon_no - 1) * 3
            codon_seq = record.seq[codon_start:codon_start + 3]
            if not str(codon_seq).__contains__("N"):
                codonDict[codon_seq] = codonDict[codon_seq] + 1
    scale_factor = 1.0 / sum(codonDict.values())
    scaled_codonDict = {k: codonDict[k] * scale_factor for k in codonDict.keys()}
    return scaled_codonDict


def get_sigmafactor(cl):
    annotation_file = infolder + "MicrobeAnnotator_out/annotation_results/" + cl + ".faa.annot"
    command = f"grep sigma {annotation_file}|grep factor |cut -f1 |sort |uniq |wc -l"
    no_of_sigma = subprocess.run([command], stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8").strip()
    return no_of_sigma


def get_CDS_count(cl, rg):
    gff_file = f"{infolder}Genome_data/{cl}/{rg}.gff"
    command = f"grep -e CDS {gff_file}|wc -l"
    cds = subprocess.run([command], stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8").strip()
    return cds


def add_functional_data(working_file=infolder + "post_cluster_data.csv"):
    # working_file = infolder + "post_cluster_data.csv"
    my_table = pd.read_csv(working_file)
    if "sigma_factor" not in my_table.columns: my_table['sigma_factor'] = my_table.Representative_Genome.apply(
        lambda x: get_sigmafactor(x))
    if "no_of_CDS" not in my_table.columns: my_table["no_of_CDS"] = my_table.cluster_ID.apply(
        lambda x: get_CDS_count(x, my_table.loc[my_table.cluster_ID == x, "Representative_Genome"].tolist()[0]))
    my_table.to_csv(working_file, index=False)

    codon_file = working_file.replace(".csv", "_codon_table.csv")
    if not os.path.exists(codon_file):
        codon_table = {}
        for cluster in my_table["cluster_ID"]:
            rg = my_table.loc[my_table.cluster_ID == cluster, "Representative_Genome"].tolist()[0]
            new_row = calculate_codon_usage(rg, cluster)
            codon_table[cluster] = new_row
        nct = pd.DataFrame(codon_table).T
        print(codon_table)
        nct.to_csv(codon_file)


def functioanl_data_organize(working_file=infolder + "post_cluster_data.csv"):
    my_table = pd.read_csv(working_file)
    for cluster in my_table.index:
        clusterId = my_table.loc[my_table.index == cluster, "cluster_ID"].tolist()[0]
        rep_genome = my_table.loc[my_table.index == cluster, "Representative_Genome"].tolist()[0]
        microbe_annot_file = f"{infolder}MicrobeAnnotator_out/annotation_results/{rep_genome}.faa.annot"
        annot_table = pd.read_csv(microbe_annot_file, sep="\t")
        gff_file = f"{infolder}Genome_data/{clusterId}/{rep_genome}.gff"
        outfile = f"{infolder}Genome_data/{clusterId}/{rep_genome}_annotations.csv"
        print(f"processing {clusterId}")
        final_annot_df = pd.DataFrame(columns=["product_id", "product", "ko", "source"])
        with open(gff_file) as gff:
            for line in gff:
                if not line.startswith("#"):
                    entry1 = line.strip().split("\t")
                    entry2 = entry1[-1].strip().split(";")
                    cds_id = entry1[0] + "_" + entry2[0].split("_")[-1]
                    prot_id = ""
                    source = ""
                    crc_count = 0
                    product = ""
                    prot_id_arr = annot_table.loc[annot_table.query_id == cds_id, "protein_id"]
                    ko = ""
                    if len(prot_id_arr) > 1:
                        source = ";".join(set(annot_table.loc[annot_table.query_id == cds_id, "database"].tolist()))
                        for x in prot_id_arr:
                            if str(x).__contains__("nan"):
                                ko_arr = annot_table.loc[annot_table.query_id == cds_id, "ko_number"].tolist()
                                ko_arr = list(filter(lambda j: str(j) != "nan", ko_arr))
                                prot_id = prot_id + ";" + str(ko_arr[0])
                                product = product + ";" + str(
                                    annot_table.loc[annot_table.query_id == cds_id, "ko_product"].tolist()[0])
                                ko = ";".join(set(ko_arr))
                            else:
                                prot_id = prot_id + ";" + str(x)
                                product = product + ";" + str(
                                    annot_table.loc[annot_table.protein_id == x, "product"].tolist()[0])
                        prot_id = prot_id.replace(";", "", 1)
                        product = ";".join(set(product.replace(";", "", 1).split(";")))
                        # print(product)
                    elif len(prot_id_arr) == 1:
                        pid = prot_id_arr.tolist()[0]
                        source = annot_table.loc[annot_table.query_id == cds_id, "database"].tolist()[0]
                        if str(pid).__contains__("nan"):
                            prot_id = annot_table.loc[annot_table.query_id == cds_id, "ko_number"].tolist()[0]
                            product = str(annot_table.loc[annot_table.query_id == cds_id, "ko_product"].tolist()[0])
                            ko = prot_id
                        else:
                            prot_id = str(pid)
                            product = str(
                                annot_table.loc[annot_table.protein_id == pid, "product"].tolist()[0])
                    if prot_id:
                        final_annot_df.loc[cds_id] = [prot_id, product, ko, source]
                    else:
                        final_annot_df.loc[cds_id] = [None, None, None, None]
        final_annot_df.index.name = "cds_id"
        final_annot_df.to_csv(outfile)


def create_ko_table(working_file=infolder + "post_cluster_data.csv"):
    my_table = pd.read_csv(working_file)
    kotable={}
    ko_uniqlist= []
    for cluster in my_table.index:
        #if cluster <5:
            clusterId = my_table.loc[my_table.index == cluster, "cluster_ID"].tolist()[0]
            rep_genome = my_table.loc[my_table.index == cluster, "Representative_Genome"].tolist()[0]
            microbe_annot_kofile = f"{infolder}MicrobeAnnotator_out/annotation_results/{rep_genome}.faa.ko"
            ko_list=os.popen(f"cat {microbe_annot_kofile} |sort |uniq ").read().split("\n")
            ko_list.remove("NA")
            ko_list.remove("ko_number")
            ko_list.remove("")
            kotable[clusterId]= ko_list
            for x in ko_list:
                #kotable[clusterId][x]=1
                if x not in ko_uniqlist:
                    ko_uniqlist.append(x)
    ko_profile_out=infolder+"ko_profile.tsv"
    kout=open(ko_profile_out,"w")
    #ko_uniqlist.sort()
    header=["Data"]+ko_uniqlist
    kout.write("\t".join(header))
    for cluster in my_table.index:
        #if cluster <5:
            clusterID = my_table.loc[my_table.index == cluster, "cluster_ID"].tolist()[0]
            kout.write(f"\n{clusterID}")
            for kid in ko_uniqlist:
                kout.write("\t1") if kid in kotable[clusterID] else kout.write("\t0")


# summarize_reassembly()
# addgtdbtkInfo()
load_site_metadata()
addSiteMetaData()
#create_members_table()
# functioanl_data_organize()
# add_functional_data()
#create_ko_table()


# cluster_reassembly("average")
# tp.tax_initiate(taxdumpdir)
# retax_kraken()
# patric_sandq()

# patric_retrieve_sort()
##command used for phylophlan:
# phylophlan_metagenomic -i ../cluster_reps/ -o output_metagenomic --nproc 15 -n 1 -d SGB.Jul20 -e 'fasta' --verbose 2>&1 | tee phylophlan_metagenomic.log ----database_folder ../../ancilary/phylophlan_databases
# addData()

## Command used for microbe annotator
# microbeannotator -l prodigal_outnames_per_cluster.txt -d ../../MicrobeAnnotatorDbs/ -o MicrobeAnnotator_out -m sword -t 5 -p 4
