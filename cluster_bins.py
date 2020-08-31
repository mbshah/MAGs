#!/usr/bin/env python
import files as config
import re
import subprocess
import multiprocessing as mp
import sys
import os
from skbio import DistanceMatrix
from skbio.tree import nj
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

samples = list(config.files_dict.keys())
infolder = config.outfolder
postfix = "_dastool_90_10"
infile = infolder + "/all_sample_summary.tsv"
outfile = infolder + "ani_results" + postfix + ".tsv"
tmp_folder = infolder + "/ani_" + postfix + "_tmp/"
matrix_file = infolder + "dist_matrix" + postfix + ".tsv"
tree_file = infolder + "tree" + postfix + ".txt"
inBins = {}
threads = 5
if len(sys.argv) > 1:
    threads = sys.argv[1]
    infile = sys.argv[2]
    outfile = sys.argv[3]
print("threads: " + threads.__str__())
print("infile " + infile)
print("outfile " + outfile)
clusters = {}
dist_matrix = {}
dm = []
bins_list = []
bin_contam_compl = {}


def build_bin_fasta_list():
    with open(infile, "r") as tsv:
        header = next(tsv)
        for line in tsv:
            entry = line.split("\t")
            isfilterok = config.filterok(
                config.filter_check(config.max_compl, config.min_compl, entry[config.compl_col]),
                config.filter_check(config.max_contam, config.min_contam, entry[config.contam_col]))
            if isfilterok:
                bin_name = entry[0]
                sample_name = "_".join(bin_name.split("_")[0:2])
                fasta_file = infolder + "/" + sample_name + "/bins/" + bin_name + ".fasta"
                inBins[bin_name] = fasta_file
                bin_contam_compl[bin_name] = {"compl": entry[config.compl_col], "contam": entry[config.contam_col]}


def ani_calc(ref):
    o = open(tmp_folder + "__tmp__" + ref, "a+")
    for i in inBins:
        if i == ref:
            pass
        else:
            cmd = ["perl", "ani-m.pl", inBins[ref], inBins[i]]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            so = proc.communicate()
            result = re.findall(r'b\'(.*)\\t(.*)\\t(.*)\\t(.*)\\t(.*)\\n\'', str(so[0]))
            result = result[0]
            res_str = "\t".join(result)
            o.write(res_str + "\n")
    o.close()


def merge_ani_results():
    onlyfiles = [f for f in os.listdir(tmp_folder) if os.path.isfile(os.path.join(tmp_folder, f))]
    for file in onlyfiles:
        if file.__contains__("__tmp__"):
            os.system("cat " + tmp_folder + file + " >>" + outfile)
            # subprocess.Popen(["rm", tmp_folder + file])


def distance_matrix():
    global dist_matrix
    global dm
    with open(outfile, "r") as f:
        header = next(f)
        for line in f:
            entry = line.split("\t")
            dist = "1"
            if float(entry[4]) >= 50:
                dist = "{:.3f}".format((1 - (float(entry[2]) / 100)))
            b0 = entry[0].replace(".fasta", "")
            b1 = entry[1].replace(".fasta", "")
            if b0 in dist_matrix:
                dist_matrix[b0][b1] = dist
            else:
                dist_matrix[b0] = {}
                dist_matrix[b0][b1] = dist
    global bins_list
    bins_list = list(sorted(dist_matrix.keys()))
    o = open(matrix_file, "w")
    o.write("bins\t" + "\t".join(bins_list))
    for i in bins_list:
        o.write(i + "\t")
        bin_arr = []
        for x in bins_list:
            if x in dist_matrix[i]:
                if dist_matrix[i][x] == dist_matrix[x][i]:
                    pass
                else:
                    dist_matrix[i][x] = dist_matrix[x][i]
                dist = str(dist_matrix[i][x])
            else:
                dist = 0
            o.write(str(dist) + "\t")
            bin_arr.append(dist)
        o.write("\n")
        dm.append(bin_arr)
    o.close()


def initialize():
    o = open(outfile, "w+")
    o.write("f1\tf2\tANI\tstdev\tmatchperc+\n")
    o.close()


def segregate_cluster(method):
    clfolder = infolder + "clusters" + postfix + "_" + method + "/"
    clusterfile = clfolder + "clusters.tsv"
    file_list = str(os.popen("ls " + infolder).read()).split("\n")
    if "clusters" + postfix + "_" + method in file_list:
        os.system("rm -r " + clfolder + "/*")
    else:
        print("making output folder for clusters")
        os.mkdir(clfolder)
    o = open(clusterfile, "w+")
    o.write(
        "Cluster\tNumberofMembers\tMembers\tAvg_gc\ttotal_len_mb\tavg_len_mb\tmax_len\tavg_completeness\tavg_contam\n")
    print("Number of Clusters found is " + str(len(clusters)) + " by method " + method)
    # more_than_one_mem=0
    for cluster in clusters:
        members = clusters[cluster]
        o.write(cluster + "\t" + str(len(members)) + "\t")
        cluster_folder = clfolder + cluster + "/"
        os.mkdir(cluster_folder)
        GC_count = 0
        gc_len = 0
        max_len = 0
        total_compl = 0
        total_contam = 0
        for x in members:
            member_fasta = infolder + "_".join(str(x).split("_")[0:2]) + "/bins/" + x + ".fasta"
            os.system("cp " + member_fasta + " " + cluster_folder)
            mem_fasts = open(member_fasta, "r")
            mem_len = 0
            # print(bin_contam_compl[x])
            total_compl = total_compl + float(bin_contam_compl[x]["compl"])
            total_contam = total_contam + float(bin_contam_compl[x]["contam"])
            for record in SeqIO.parse(mem_fasts, "fasta"):
                GC_count = GC_count + record.seq.count("G") + record.seq.count("C")
                gc_len = gc_len + record.__len__()
                mem_len = mem_len + record.__len__()
            if mem_len > max_len: max_len = mem_len
            o.write("\"" + x + "\";")
        totlen = float("{:.2f}".format((gc_len / 1000) / 1000))
        avggc = float("{:.2f}".format((GC_count / gc_len) * 100))
        avglen = float("{:.2f}".format(totlen / len(members)))
        avgcompl = float("{:.2f}".format(total_compl / len(members)))
        avgcontam = float("{:.2f}".format(total_contam / len(members)))
        maxlen = float("{:.2f}".format(max_len / 1000000))
        o.write("\t" + str(avggc) + "\t" + str(totlen) + "\t" + str(avglen) + "\t" + str(maxlen) + "\t" + str(
            avgcompl) + "\t" + str(avgcontam) + "\n")
    o.close()


def clustering(method):
    global clusters
    clusters = {}
    # writing tree file
    o = open(tree_file, "w")
    matrix = DistanceMatrix(dm, bins_list)
    tree = nj(matrix)
    o.write(str(tree))
    o.close()

    # Actual clustering
    mat = np.array(dm)
    y = mat[np.triu_indices(len(bins_list),
                            1)]  # converting the matrix to vector (also doable by scikit suqareform but that was not working properly)
    linkage_mtrx = linkage(y, method)
    ct = cut_tree(linkage_mtrx, height=[0.05])
    x = 0
    while x < len(ct):
        cluster_number = ct[x][0]
        clustername = "EUL_" + str(cluster_number)
        member = bins_list[x]
        if clustername in clusters:
            clusters[clustername].append(member)
        else:
            clusters[clustername] = [member]

        x = x + 1

    # PLOTTING dendrogram
    plt.figure(figsize=(40, 10))
    dendrogram(linkage_mtrx, labels=bins_list, color_threshold=0.05)
    plt.title("Clustering Dendrogram by " + method + " Method")
    plt.ylim(0, 0.1)
    plt.axhline(0.05, linestyle="--", c="black")
    plt.savefig(infolder + "Clusters_" + method + postfix + ".svg", dpi=1000, bbox_inches='tight')


def main():
    build_bin_fasta_list()
    print("number of good bins is " + str(len(inBins)))
    initialize()
    pool = mp.Pool(threads)
    pool.map(ani_calc, [ref for ref in inBins])
    merge_ani_results()
    distance_matrix()
    # for method in ("single","average","complete"):
    for method in ['average']:
        clustering(method)
        segregate_cluster(method)


if __name__ == "__main__":
    main()
