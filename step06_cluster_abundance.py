import os
import shutil
from pathlib import Path

import pandas as pd

import config as config
import fileinput
import concurrent.futures

infolder = config.outfolder + "dRep__gff/"
raw_fastq_folder = config.raw_fastq
rep_folder=infolder+"Rep_genomes/"
Path(rep_folder).mkdir(exist_ok=True)
processes=13
cluster_profile={}
myclusters=[]
samples=[]
cluster_table = infolder + "post_cluster_data.csv"
cluster_df = pd.read_csv(cluster_table, index_col="cluster_ID")
dryrun=True


def create_custom_db():
    clusters = cluster_df.index.tolist()
    clusters_folder=infolder+"Genome_data/"
    rep_list=[]
    for cl in clusters:
        cl_f=clusters_folder+cl+"/"
        rep_file=cl_f+cluster_df.loc[cl,'Representative_Genome'].strip()
        dest_file=rep_folder+cl+"_rep.fasta"
        command = f"cp -f {rep_file} {dest_file}"
        os.system(command)
        rep_list.append(dest_file)
        contig_number = 1
        for line in fileinput.input(dest_file, inplace=True):
            line = str(line.rstrip())
            if line.startswith(">"):
                newline = f">ctg{contig_number}|{cl}"
                contig_number = contig_number + 1
                print(newline)
            else:
                print(line)
    files_str = ",".join(rep_list)
    command = f"bowtie2-build {files_str} {rep_folder}reps_bt_index"
    os.system(command)
        #print(f"{rep_file}\t{Path(rep_file).is_file()}")



def map_and_filter(sample,read1,read2,bt2index,file):
    sample_folder=config.outfolder+sample+"/"
    threads=processes
    output_tmp=sample_folder+"cluster_remapping_"+sample+"_tmp_"+file+".sam"
    output = sample_folder + "cluster_remapping_" + sample + ".sam"
    bkp_loc_sam="/mnt/biodiv/Manan/MAGs/outfile/"+sample+"/cluster_remapping_"+sample+".sam"

    if not os.path.exists(bkp_loc_sam):
        command=f"bowtie2 -x {bt2index} -1 {read1} -2 {read2} -p {threads} -S {output_tmp}"
        print(command)
        if not dryrun:os.system(command)
        command = f"samtools view -@ {threads} -G 4 {output_tmp} >> {output}"
        print(command)
        if not dryrun:os.system(command)
        command=f"rm {output_tmp}"
        print(command)
        if not dryrun:os.system(command)
    return output


def read_samfiles(infile):
    clus_tmp=[]
    sample = os.path.basename(os.path.dirname(infile))
    outfile=infile+"_out.list"
    o = open(outfile, "a")
    #print(f"extracting from {infile}")
    tmp_map = {}
    with open(infile, "r") as sam:
        for line in sam:
            line_arr = line.split("\t")
            key = line_arr[0]
            match = line_arr[2]
            if key in tmp_map:
                if tmp_map[key][1] == match:
                    tmp_map[key][0] = tmp_map[key][0] + 1
                    if tmp_map[key][0] == 2:
                        o.write(f"{sample}\t{key}\t{match}\n")
                        cluster=match.split("|")[1]
                        clus_tmp.append(f"{sample};{cluster}")
                else:
                    dict(tmp_map).pop(key)
            else:
                tmp_map[key] = [1, match]
    o.close()
    return clus_tmp


def add_to_cluster_profile(sample,cluster):
    global myclusters
    global cluster_profile
    #print (sample,cluster)
    if cluster in myclusters:pass
    else:myclusters.append(cluster)
    if sample in cluster_profile:
        #print("here")
        if cluster in cluster_profile[sample]:
            cluster_profile[sample][cluster]=cluster_profile[sample][cluster]+1
        else: cluster_profile[sample][cluster]=1
    else:
        cluster_profile[sample]={}
        cluster_profile[sample][cluster]=1


def parallel_read_samfile(samfile):
    #num_lines = sum(1 for line in open(samfile))
    path=samfile.replace(os.path.basename(samfile),"")
    print (f"processeing{samfile}")
    command=f"split -n l/{processes} -d {samfile} {path}tmp"
    os.system(command)
    chunks=[]
    outfile=infolder + "cluster_abundance_list.dump"
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in [f for f in filenames if f.__contains__("tmp")]:
            fullpath=dirpath+filename
            chunks.append(fullpath)
    clus_tmp_combined=[]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for clus_temp in executor.map(read_samfiles, chunks):
            clus_tmp_combined=clus_tmp_combined+clus_temp
    #print(clus_tmp_combined)
    for entry in clus_tmp_combined:
        entry_arr=entry.split(";")
        add_to_cluster_profile(entry_arr[0],entry_arr[1])
    for chunk in chunks:
        chunk_out=chunk+"_out.list"
        command=f"cat {chunk_out} >>{outfile}"
        os.system(command)


def mapping_queue():
    intsv=config.outfolder+"all_sample_summary.tsv"
    global samples
    global myclusters
    global cluster_profile
    with open(intsv,"r")as tsvfile:
        header=next(tsvfile)
        for line in tsvfile:
            entry=line.strip().split("\t")
            sample="_".join(entry[0].split("_")[0:2])
            if sample in samples:pass
            else:samples.append(sample)
    #samples = ["EULs_Z301WO","EULs_N051NO"]
    for sample in sorted(samples):
        raw_sample_folder=raw_fastq_folder+sample+"/"
        file_pairs=[]
        for files in os.listdir(raw_sample_folder):
            file_pair="_".join(files.split("_")[0:-1])
            if file_pair in file_pairs:pass
            else:file_pairs.append(file_pair)
        sam_file=""
        for set in file_pairs:
            r1=raw_sample_folder+set+"_1.fq.gz"
            r2 = raw_sample_folder+set + "_2.fq.gz"
            index = rep_folder + "reps_bt_index"
            #print (f"\n\n{sample}\t{r1}\t{r2}\t{index}\n")
            sam_file=map_and_filter(sample, r1, r2, index,set)
        path_to_samfile_M = "/mnt/biodiv/Manan/" + "/".join(sam_file.split("/")[-4:])
        command=f"mv {sam_file} {path_to_samfile_M}"
        if not dryrun:
            print(command)
            os.system(command)
        if not os.path.exists(sam_file):
            print(f"Copying from {path_to_samfile_M}")
            shutil.copy(path_to_samfile_M,sam_file)
        parallel_read_samfile(sam_file)
        for dirpath, dirnames, filenames in os.walk(os.path.dirname(sam_file)):
            for filename in [f for f in filenames if f.__contains__("tmp")]:
                os.system(f"rm {dirpath}/{filename}")
        os.remove(sam_file) ##also for remote storage of samfiles, hash out if not required
    o1=open(infolder+"cluster_abundance_profile.tsv","w")
    clusters=sorted(myclusters)
    o1.write("abundance\t")
    for cluster in clusters:
        o1.write(cluster+"\t")
    o1.write("\n")
    for sample in sorted(samples):
        o1.write (sample+"\t")
        for cluster in sorted(clusters):
            value= 0
            if cluster in cluster_profile[sample]:
                value=cluster_profile[sample][cluster]
            o1.write(str(value)+"\t")
        o1.write("\n")
    o1.close()


    #o2 = open(infolder + "cluster_abundance_profile_transformed.tsv", "w")
    #clusters = sorted(myclusters)
    #o2.write("abundance\t")
    #for cluster in clusters:
    #    o2.write(cluster + "\t")
    #o2.write("\n")
    #for sample in sorted(samples):
    #    o2.write(sample + "\t")
    #    for cluster in sorted(clusters):
    #        value = 0
    #        if cluster in cluster_profile[sample]:
    #            value = cluster_profile[sample][cluster]
    #        o2.write(str(value) + "\t")
    #    o2.write("\n")
    #o2.close()


#create_custom_db()
o= open(infolder + "cluster_abundance_list.dump", "w")
o.close()
mapping_queue()
