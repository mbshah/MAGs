import files as config
from Bio import SeqIO
import os
import tax_profiler as tp
import json
import time
import re
from bs4 import BeautifulSoup

postfix = config.postfix
method = config.method
infolder = config.outfolder + "clusters" + postfix + "_" + method + "/"
threads = 15
taxdumpdir = "/home/hb0358/PycharmProjects/mbs_general/ancilary/new_taxdump"


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


def retax_kraken():
    kraken_path = "~/kraken2_int/kraken2"
    kraken_db = "~/mbs_workfolder/ancilary/kraken_bact/"
    clfolder = infolder
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
                    #if key=="Source":print(html_file+"\t"+value)
                    pat_dict_headers.append(key)
    pat_dict_headers = sorted(list(set(pat_dict_headers)))
    pat_dict_headers.remove("Source")
    pat_dict_headers.remove("AMR Mechanism")
    return ret_dict


def patric_retrieve_sort():
    patric_submissions = {}
    for dirpath, dirnames, filenames in os.walk(infolder):
        for filename in [f for f in filenames if f.endswith(".json")]:
            patric_submissions[os.path.basename(dirpath)] = (dirpath, filename)
    genus=""
    for cluster in patric_submissions:
        cluster_dir = patric_submissions[cluster][0]
        json_filename = "/".join(patric_submissions[cluster])
        with open(json_filename, "r") as json_file:
            data = json.load(json_file)
        genus_id=data["taxonomy_id"]
        genus=f"{tp.names[genus_id]}({genus_id})"
        patric_html_onserver = str(
            "ws:" + data["output_path"] + "/." + data["output_file"] + "/FullGenomeReport.html").replace(" ", "\ ")
        newhtml = cluster_dir + "/" + cluster + "_patricOut.html"
        patric_html_onpc = cluster_dir + "/FullGenomeReport.html"
        command = f"p3-cp {patric_html_onserver} {cluster_dir}"
        os.system(command)
        command = f"mv {patric_html_onpc} {newhtml}"
        os.system(command)
        patric_submissions[cluster] = extract_from_patric(newhtml)
        patric_submissions[cluster]["Kraken_Genus"]=genus
    finalfile = infolder + "clusters.tsv"
    tmpfile = infolder + "clusters_tmp.tsv"
    o=open(finalfile,"w")
    pat_dict_headers.insert(0,"Patric_Quality")
    pat_dict_headers.insert(0, "Kraken_Genus")
    with open(tmpfile,"r") as clusterfile:
        header=str(next(clusterfile)).strip()+"\t"+"\t".join(pat_dict_headers)+"\n"
        o.write(header)
        for line in clusterfile:
            linearr=line.strip().split("\t")
            cluster=linearr[0]
            for key in pat_dict_headers:
                if key in patric_submissions[cluster]:value=patric_submissions[cluster][key]
                else:value="-"
                linearr.append(value)
            newline="\t".join(linearr)+"\n"
            o.write(newline)
    o.close()




#cluster_reassembly("average")
#summarize_reassembly()
#tp.tax_initiate(taxdumpdir)
#retax_kraken()
#patric_sandq()
#patric_retrieve_sort()
