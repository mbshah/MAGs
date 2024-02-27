import statistics
import subprocess
import config as config
from Bio import SeqIO
import os
import step03_tax_profiler as tp
import json
import time
import re
from bs4 import BeautifulSoup

postfix = config.postfix
method = config.method
infolder = config.outfolder + "clusters" + postfix + "_" + method + "/"
threads = 15
taxdumpdir = "../ancilary/new_taxdump"
kraken_path = "~/manan_WMG/tools/kraken2/kraken2"
kraken_db = "~/manan_WMG/ancilary/kraken_bact/"
siteMetadata = {}
siteClAbundance={}


def cluster_reassembly(method):
    clfolder = infolder
    cfile_in = clfolder + "clusters.tsv"
    # os.system("bash conda activate amos")
    # os.system("source ~/anaconda3/bin/activate /home/hb0358/anaconda3/envs/amos")
    with open(cfile_in, "r+") as c:
        header = next(c)
        for cluster_line in c:
            cluster_info = cluster_line.strip().split("\t")
            cluster = cluster_info[0]
            cluster_folder = clfolder + cluster + "/"
            if cluster_info[2].endswith(";"): cluster_info[2] = cluster_info[2][0:-1]
            members = cluster_info[2].strip().split(";")
            avglen = cluster_info[5]
            combined_fasta = cluster_folder + cluster + "_raw_combined.fasta"
            command = "cat " + cluster_folder + "*.fasta >" + combined_fasta
            os.system(command)
            ##un-comment the assembler you want to run.
            # print (" running canu on the cluster "+cluster)
            # canudir = "/home/hb0358/PycharmProjects/mbs_general/tools/canu-2.0/Linux-amd64/"
            # canu_out = cluster_folder + "canu_out/"
            # command = canudir + "bin/canu -nanopore-raw " + combined_fasta + " genomeSize=" + avglen.__str__() + "m -p " + cluster + " -d " + canu_out + " -minInputCoverage=0.1 -stopOnLowCoverage=0.1"
            # os.system(command)
            command = "toAmos -s " + combined_fasta + " -o " + cluster_folder + cluster + ".afg"
            os.system(command)
            command = "minimus2 " + cluster_folder + cluster + " -D REFCOUNT=0"
            os.system(command)


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
    infile = infolder + "clusters.tsv"
    tmpfile = infolder + "clusters_tmp.tsv"
    o = open(tmpfile, "w")
    with open(infile, "r") as file:
        header = next(file)
        header = header.strip() + "\treassembly_gc\treassembly_size_mb\treassembly_longestContig_kb\n"
        o.write(header)
        for line in file:
            entry = line.split("\t")
            entry[-1] = entry[-1].strip()
            clustername = entry[0]
            fastafile = infolder + clustername + "/" + clustername + ".fasta"
            if int(entry[1]) == 1:
                command = "rm -r " + infolder + clustername + "/" + clustername + "*"
                os.system(command)
                command = "cat " + infolder + clustername + "/""*.fasta >" + fastafile
                os.system(command)
            gc_count, cluster_len, max_len = fasta_info(fastafile)
            gc_perc = (gc_count / cluster_len) * 100
            entry.append(str(float("{:.2f}".format(gc_perc))))
            entry.append(str(float("{:.2f}".format(cluster_len / 1000000))))
            entry.append(str(float("{:.2f}".format(max_len / 1000))))
            new_line = "\t".join(entry)
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
            patric_JSON_maker(cluster, genus, genus_id)


def patric_JSON_maker(cluster, genus, gen_id):
    contigs_file_local = f"{infolder}/{cluster}/{cluster}.fasta"
    contigs_file_remote = f"/mananbshah@patricbrc.org/unidue/infolder/{cluster}.fasta"
    command = ("p3-cp " + contigs_file_local + " ws:" + contigs_file_remote + " -m fasta=contigs")
    os.system(command)
    json_file_path = f"{infolder}/{cluster}/{cluster}_particCGA.json"
    json_content = {
        "genome_size": "5M",
        "pilon_iter": 2,
        "debug_level": 0,
        "queue_nowait": "0",
        "min_contig_cov": 5,
        "domain": "Bacteria",
        "input_type": "contigs",
        "public": "0",
        "code": "11",
        "min_contig_len": 300,
        "trim": "0",
        "skip_indexing": "0",
        "taxonomy_id": gen_id,
        "scientific_name": cluster + "_" + genus,
        "contigs": contigs_file_remote,
        "output_file": f"{cluster}_{genus}",
        "racon_iter": 2,
        "output_path": "/mananbshah@patricbrc.org/unidue/outfile",
        "recipe": "auto"
    }
    with open(json_file_path, "w") as json_file:
        json.dump(json_content, json_file)


count = 0


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


def patric_sandq():
    import os.path
    patric_submissions = {}
    for dirpath, dirnames, filenames in os.walk(infolder):
        for filename in [f for f in filenames if f.endswith(".json")]:
            patric_submissions[os.path.basename(dirpath)] = os.path.join(dirpath, filename)
    prr = 0
    for cluster in sorted(patric_submissions):
        if prr <= 5:
            prr = prr + 1
            print(f"submit job for {cluster}")
            command = "appserv-start-app ComprehensiveGenomeAnalysis " + patric_submissions[
                cluster] + " \"mananbshah@patricbrc.org/unidue/\""
            # os.system(command)
        if prr == 5:
            prr = 0
            time.sleep(0)
            # print("sleeping for 10 mins to give time for the jobs to finish")


qlty_c = re.compile(
    r"Based on the annotation statistics and a comparison to other genomes in PATRIC within this same species, this genome appears to be of (.*) quality\.")
pat_dict_headers = []


def extract_from_patric(html_file):
    global pat_dict_headers
    ret_dict = {}
    soup = BeautifulSoup(open(html_file), "html.parser")
    for x in soup.find_all('section'):
        s2 = BeautifulSoup(str(x), "html.parser")
        for child in s2.recursiveChildGenerator():
            if (child.name):
                if str(child.text).__contains__("Based"):
                    text = str(child.text).replace("\n", "")
                    text = " ".join(text.split())
                    ret_dict["Patric_Quality"] = qlty_c.search(text).group(1)
    for x in soup.find_all('table'):
        s3 = BeautifulSoup(str(x), 'html.parser')
        for child in s3.recursiveChildGenerator():
            if child.name == 'tr':
                text = str(child.text).replace("\n", "~~").replace("\s", "")
                # print (text)
                t_arr = text.split("~~")
                t2 = []
                for i in t_arr:
                    i = i.lstrip()
                    if i != "": t2.append(i)
                if len(t2) > 1:
                    key = "-".join(t2[0:-1])
                    value = t2[-1]
                    ret_dict[key] = value
                    # if key=="Source":print(html_file+"\t"+value)
                    pat_dict_headers.append(key)
    pat_dict_headers = sorted(list(set(pat_dict_headers)))
    pat_dict_headers.remove("Source")
    pat_dict_headers.remove("AMR Mechanism")
    return ret_dict


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


def patric_retrieve_sort():
    patric_submissions = {}
    for dirpath, dirnames, filenames in os.walk(infolder):
        for filename in [f for f in filenames if f.endswith(".json")]:
            if filename.startswith("EUL"):
                patric_submissions[os.path.basename(dirpath)] = (dirpath, filename)
    genus = ""
    for cluster in patric_submissions:
        # print(cluster)
        cluster_dir = patric_submissions[cluster][0]
        json_filename = "/".join(patric_submissions[cluster])
        with open(json_filename, "r") as json_file:
            data = json.load(json_file)
        genus_id = data["taxonomy_id"]
        genus = f"{tp.names[genus_id]}({genus_id})"
        patric_html_onserver = str(
            "ws:" + data["output_path"] + "/." + data["output_file"] + "/FullGenomeReport.html").replace(" ", "\ ")
        newhtml = cluster_dir + "/" + cluster + "_patricOut.html"
        patric_html_onpc = cluster_dir + "/FullGenomeReport.html"
        if not os.path.isfile(newhtml):
            command = f"p3-cp {patric_html_onserver} {cluster_dir}"
            os.system(command)
            command = f"mv {patric_html_onpc} {newhtml}"
            os.system(command)
        patric_annot_onserver = str(
            "ws:" + data["output_path"] + "/." + data["output_file"] + "/.annotation/annotation.txt").replace(" ", "\ ")
        newannotfile = cluster_dir + "/" + cluster + "_annotation.tsv"
        patric_annot_onpc = cluster_dir + "/annotation.txt"
        patric_submissions[cluster] = extract_from_patric(newhtml)
        # print(patric_submissions[cluster])
        patric_submissions[cluster]["Kraken_Genus"] = genus
        patric_submissions[cluster]["Kraken_Phyla"] = find_phyla(genus.rsplit("(")[1].replace(")", ""))
        if not os.path.isfile(newannotfile):
            command = f"p3-cp {patric_annot_onserver} {cluster_dir}"
            os.system(command)
            command = f"mv {patric_annot_onpc} {newannotfile}"
            os.system(command)
    finalfile = infolder + "clusters.tsv"
    tmpfile = infolder + "clusters_tmp.tsv"
    o = open(finalfile, "w")
    pat_dict_headers.insert(0, "Patric_Quality")
    pat_dict_headers.insert(0, "Kraken_Genus")
    pat_dict_headers.insert(0, "Kraken_Phyla")
    # pat_dict_headers.insert(0,"site_list")
    with open(tmpfile, "r") as clusterfile:
        header = str(next(clusterfile)).strip() + "\t" + "\t".join(pat_dict_headers) + "\n"
        o.write(header)
        for line in clusterfile:
            linearr = line.strip().split("\t")
            cluster = linearr[0]
            for key in pat_dict_headers:
                if key in patric_submissions[cluster]:
                    value = patric_submissions[cluster][key]
                else:
                    value = "-"
                linearr.append(value)
            newline = "\t".join(linearr) + "\n"
            o.write(newline)
    o.close()


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

def getWtSiteData(sites,parameter,cluster):
    if str(sites).endswith(""): sites = "".join(sites.rsplit("\"", 1))
    if str(sites).endswith(";"): sites = "".join(sites.rsplit(";", 1))
    sites_arr = sites.split(";")
    parameter_array = []
    total_count=0
    for site in sites_arr:
        site=site.strip()
        parameter_value = float(siteMetadata[site][parameter]) * float(siteClAbundance[cluster][site])
        total_count=total_count+float(siteClAbundance[cluster][site])
        parameter_array.append(float(parameter_value))
    retval=sum(parameter_array)/total_count
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
        header = next(md_file).strip().split(seperator)
        for line in md_file:
            entry = line.strip().split(seperator)
            site = entry[0]
            for i in range(1, len(entry)):
                if site not in siteMetadata: siteMetadata[site] = {}
                siteMetadata[site][header[i].strip()] = entry[i]
                # print(f"{site}\t{header[i]}\t{entry[i]}")
    with open(config.outfolder +"cluster_abundance_profile_normalized_deseq.tsv","r") as cl_abundance:
        header_sites=next(cl_abundance)
        header_sites=header_sites.split()
        for line in cl_abundance:
            entry=line.split()
            cluster=entry[0]
            counts=entry[1:]
            for i in range(len(header_sites)):
                if cluster in siteClAbundance:pass
                else:siteClAbundance[cluster]={}
                siteClAbundance[cluster][header_sites[i]]=counts[i]


def addData():
    import pandas as pd
    phylophlan_tax_file = infolder + "phylophlan/output_metagenomic.tsv"
    new_taxa = {}
    with open(phylophlan_tax_file, "r") as newt:
        newt_head = next(newt)
        for line in newt:
            line = line.strip().replace("__", "~")
            bin = line.split("\t")[0]
            taxa = line.split("\t")[1].split(":")[2]
            phyla = "Other"
            for lvl in taxa.split("|"):
                if lvl.startswith("p~"):
                    phyla = lvl.replace("p~", "")
            new_taxa[bin] = {}
            new_taxa[bin]["taxa"] = taxa
            new_taxa[bin]["phyla"] = phyla
            # print(f"{bin}\t{taxa}\t{phyla}")
    working_file = infolder + "clusters.tsv"
    my_table = pd.read_csv(working_file, sep="\t")
    if "phylophlan_phyla" not in my_table.columns: my_table["phylophlan_phyla"] = my_table.Cluster.apply(
        lambda x: new_taxa[x]["phyla"])
    my_table["Kraken_Phyla"] = my_table.Kraken_Genus.apply(lambda x: find_phyla(x.rsplit("(")[1].replace(")", "")))
    my_table.to_csv(infolder + "clusters_wphyp.tsv", sep="\t", index=False)


def get_perc_coding(Cluster,retTyp="total"):
    cl=Cluster
    avg_list=[]
    max_sizes = {}
    cluster_lstfile = infolder + cl + "/" + cl + "_annotation.tsv"
    total_len = 0
    coded_len = 0
    combined_total_len=0
    combined_coded_len = 0
    # finding max sizes of contigs
    with open(cluster_lstfile, "r") as clusterfile:
        header = next(clusterfile)
        for line in clusterfile:
            line_entry = line.strip().split("\t")
            ent_cid = line_entry[0]
            ent_strand = line_entry[6]
            if ent_strand.__contains__("+"):
                ent_end = int(line_entry[5].replace(">", "").replace("<", ""))
            else:
                ent_end = int(line_entry[4].replace(">", "").replace("<", ""))
            if not ent_cid in max_sizes:
                max_sizes[ent_cid] = ent_end
            else:
                if max_sizes[ent_cid] < ent_end: max_sizes[ent_cid] = ent_end
    # populating contigs
    with open(cluster_lstfile, "r") as clusterfile:
        header = next(clusterfile)
        cid = ""
        forward_arr = [0]
        reverse_arr = [0]
        combined_arr=[0]
        contig_no = 0
        for line in clusterfile:
            line_entry = line.strip().split("\t")
            ent_cid = line_entry[0]
            ent_strand = line_entry[6]
            if cid != ent_cid:
                forward_arr.pop(0)
                reverse_arr.pop(0)
                combined_arr.pop(0)
                # print(str(forward_arr) + "\n" + str(reverse_arr))
                combined_total_len=combined_total_len+len(combined_arr)
                combined_coded_len=combined_coded_len+combined_arr.count(1)
                total_len = total_len + len(forward_arr) + len(reverse_arr)
                coded_len = coded_len + forward_arr.count(1) + reverse_arr.count(1)
                if len(forward_arr)>0:
                    avg_list.append((forward_arr.count(1) + reverse_arr.count(1))*100/(len(forward_arr)
                                                                                       + len(reverse_arr)))
                contig_no += 1
                cid = ent_cid
                forward_arr = [0] * (max_sizes[ent_cid] + 1)
                reverse_arr = [0] * (max_sizes[ent_cid] + 1)
                combined_arr = [0] * (max_sizes[ent_cid] + 1)
            if ent_strand.__contains__("+"):
                ent_start = int(line_entry[4].replace(">", "").replace("<", ""))
                ent_end = int(line_entry[5].replace(">", "").replace("<", ""))
                for i in range(ent_start, ent_end + 1):
                    # print(f"{i}\t{forward_arr}\t{ent_cid}\t{cl}")
                    forward_arr[i] = 1
                    combined_arr[i]=1
            else:
                ent_end = int(line_entry[4].replace(">", "").replace("<", ""))
                ent_start = int(line_entry[5].replace(">", "").replace("<", ""))
                for i in range(ent_start, ent_end + 1):
                    # print(f"{i}\t{len(forward_arr)}\t{ent_cid}\t{cl}")
                    reverse_arr[i] = 1
                    combined_arr[i] = 1
        forward_arr.pop(0)
        reverse_arr.pop(0)
        combined_arr.pop(0)
        # print(str(forward_arr) + "\n" + str(reverse_arr))
        total_len = total_len + len(forward_arr) + len(reverse_arr)
        coded_len = coded_len + forward_arr.count(1) + reverse_arr.count(1)
        coded_perc = (coded_len / total_len) * 100
        # print(f"{cl}\t{total_len}\t{coded_perc}\t{non_coded_perc}")
        coded_perc_cont=sum(avg_list)/len(avg_list)

        combined_total_len = combined_total_len + len(combined_arr)
        combined_coded_len = combined_coded_len + combined_arr.count(1)
        combined_perc=(combined_coded_len*100)/combined_total_len


        if retTyp=="Cont":return coded_perc_cont
        elif retTyp=="merged_total":return combined_perc
        else:return(coded_perc)

def get_sigmafactor(cl):
    annotation_file=infolder + cl + "/" + cl + "_annotation.tsv"
    command=f"cut {annotation_file} -f 8 |grep \"sigma factor\" |wc -l"
    no_of_sigma=subprocess.run([command], stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8").strip()
    return no_of_sigma


def addSiteMetaData():
    import pandas as pd
    working_file = infolder + "clusters_wphyp_sited.tsv"
    my_table = pd.read_csv(working_file, sep="\t")
    if "site_list" not in my_table.columns: my_table['site_list'] = my_table.Members.apply(lambda x: find_sites(x))
    if "site_pH" not in my_table.columns: my_table['site_pH'] = my_table.site_list.apply(
        lambda x: getSiteData(x, "pH1"))
    if "site_Elevation" not in my_table.columns: my_table['site_Elevation'] = my_table.site_list.apply(
        lambda x: getSiteData(x, "elevation"))
    if "site_Temperature" not in my_table.columns: my_table['site_Temperature'] = my_table.site_list.apply(
        lambda x: getSiteData(x, "temperature"))
    if "site_TP" not in my_table.columns: my_table['site_TP'] = my_table.site_list.apply(
        lambda x: getSiteData(x, "TP"))
    if "weighted_TP" not in my_table.columns: my_table['weighted_TP'] = my_table.apply(
        lambda x: getWtSiteData(x["site_list"], "TP",x["Cluster"]),axis=1)
    if "site_pH_Array" not in my_table.columns: my_table['site_pH_Array'] = my_table.site_list.apply(
        lambda x: getSiteDataArray(x, "pH1"))
    if "site_Elevation_Array" not in my_table.columns: my_table['site_Elevation_Array'] = my_table.site_list.apply(
        lambda x: getSiteDataArray(x, "elevation"))
    if "site_Temperature_Array" not in my_table.columns: my_table['site_Temperature_Array'] = my_table.site_list.apply(
        lambda x: getSiteDataArray(x, "temperature"))
    if "site_TP_Array" not in my_table.columns: my_table['site_TP_Array'] = my_table.site_list.apply(
        lambda x: getSiteDataArray(x, "TP"))
    if "perc_coding" not in my_table.columns: my_table['perc_coding'] = my_table.Cluster.apply(
        lambda x: get_perc_coding(x))
    if "perc_coding_cont" not in my_table.columns: my_table['perc_coding_cont'] = my_table.Cluster.apply(
        lambda x: get_perc_coding(x,"Cont"))
    if "perc_coding_combined" not in my_table.columns: my_table['perc_coding_combined'] = my_table.Cluster.apply(
        lambda x: get_perc_coding(x,"merged_total"))
    if "sigma_factor" not in my_table.columns: my_table['sigma_factor'] = my_table.Cluster.apply(
        lambda x: get_sigmafactor(x))
    my_table.to_csv(infolder + "clusters_wphyp_sited.tsv", sep="\t", index=False)

def addgtdbtkInfo():
    import pandas as pd
    working_file = infolder + "clusters_wphyp_sited.tsv"
    magsTable = pd.read_csv(working_file, sep="\t")
    gtdb_assignment_file=infolder + "gtdbtk/gtdbtk_classify_out/classify/gtdbtk.bac120.summary.tsv"
    gtdb_table=pd.read_csv(gtdb_assignment_file,sep="\t")
    if "gtdb_phyla" not in magsTable.columns: magsTable["gtdb_phyla"] = magsTable.Cluster.apply(
        lambda x:getgtdbInfo(1,gtdb_table.loc[gtdb_table['user_genome']==x,"pplacer_taxonomy"].tolist()[0]))
    if "gtdb_genus" not in magsTable.columns: magsTable["gtdb_genus"] = magsTable.Cluster.apply(
        lambda x:getgtdbInfo(5,gtdb_table.loc[gtdb_table['user_genome']==x,"pplacer_taxonomy"].tolist()[0]))
    magsTable.to_csv(infolder + "clusters_gtdbtk.tsv", sep="\t", index=False)

def getgtdbInfo(lev, tax_str):
    tax=str(tax_str).strip().split(";")
    taxa=tax[lev]
    taxa_refined=taxa.split("__")[-1]
    if taxa_refined=="":taxa_refined="N/A"
    #print(str(tax_str),taxa_refined)
    return taxa_refined


# cluster_reassembly("average")
# summarize_reassembly()
# tp.tax_initiate(taxdumpdir)
# retax_kraken()
# patric_sandq()
#load_site_metadata()
# patric_retrieve_sort()
##command used for phylophlan:
# phylophlan_metagenomic -i /mnt/xio_5t/botany/MAGs/outfile/dRep__gff/Rep_genomes/ -o output_metagenomic_phlophlan --nproc 15 -n 1 -d SGB.Jul20 -e 'fasta' --database_folder /mnt/xio_5t/botany/phylophlan_databases --verbose 2>&1 | tee phylophlan_metagenomic.log
# phylophlan -i ../cluster_reps/ -d phylophlan --diversity low -f supermatrix_aa.cfg --genome_extension fasta --nproc 15
# addData()
#addSiteMetaData()
addgtdbtkInfo()
