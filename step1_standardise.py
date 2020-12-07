#!/usr/bin/env python
import os
import config as f
from pathlib import Path

files_dict = dict(f.files_dict)
curation_file = f.curation_file
cur_sum = {}


def standardise(sample, file, outfolder):
    print("starting with sample:" + str(sample))
    outfolder = outfolder + str(sample)
    Path(outfolder + "/bins").mkdir(parents=True, exist_ok=True)
    file_id = Path(file).stem
    count_unclassified = 0
    count_classified = 0
    count_binned = 0
    bins = {}
    scaffolds_dict = {}
    # opening of binfile to standardise and summarize
    print("reading sample bins")
    with open(file, "r") as binfile:
        next(binfile)  # header_line
        for entry in binfile:
            entry_arr = entry.split("\t")
            if entry_arr[-1].isspace():
                count_unclassified = count_unclassified + 1
            else:
                count_classified = count_classified + 1
                oldbinid = entry_arr[-1]
                if oldbinid.startswith("maxbin") or oldbinid[0].isdigit():
                    pass
                else:
                    count_binned = count_binned + 1
                    newbinid = oldbinid.replace(file_id, sample)
                    # print(oldbinid+"\t"+ newbinid)
                    newbinid = newbinid.strip()
                    if newbinid in bins:  # modify values of existing bin
                        count = int(bins[newbinid]["count"]) + 1
                        bins[newbinid]["count"] = count
                        sumgc = (float(bins[newbinid]["gcsum"]) + float(entry_arr[1]))
                        bins[newbinid]["gcsum"] = sumgc
                        if float(entry_arr[1]) < float(bins[newbinid]["mingc"]): bins[newbinid]["mingc"] = entry_arr[1]
                        if float(entry_arr[1]) > float(bins[newbinid]["maxgc"]): bins[newbinid]["maxgc"] = entry_arr[1]
                        sumcov = (float(bins[newbinid]["covsum"]) + float(entry_arr[2]))
                        bins[newbinid]["covsum"] = sumcov
                        if float(entry_arr[2]) < float(bins[newbinid]["mincov"]): bins[newbinid]["mincov"] = entry_arr[
                            2]
                        if float(entry_arr[2]) > float(bins[newbinid]["maxcov"]): bins[newbinid]["maxcov"] = entry_arr[
                            2]
                        sumlen = (float(bins[newbinid]["lensum"]) + float(entry_arr[3]))
                        bins[newbinid]["lensum"] = sumlen
                        if float(entry_arr[3]) < float(bins[newbinid]["minlen"]): bins[newbinid]["minlen"] = entry_arr[
                            3]
                        if float(entry_arr[3]) > float(bins[newbinid]["maxlen"]): bins[newbinid]["maxlen"] = entry_arr[
                            3]
                    else:  # initialize values for new bin
                        gc = entry_arr[1]
                        cov = entry_arr[2]
                        length = entry_arr[3]
                        bins[newbinid] = {}
                        bins[newbinid]["count"] = 1
                        bins[newbinid]["gcsum"] = gc
                        bins[newbinid]["mingc"] = gc
                        bins[newbinid]["maxgc"] = gc
                        bins[newbinid]["covsum"] = cov
                        bins[newbinid]["mincov"] = cov
                        bins[newbinid]["maxcov"] = cov
                        bins[newbinid]["lensum"] = length
                        bins[newbinid]["minlen"] = length
                        bins[newbinid]["maxlen"] = length
                    scaffold = ".".join(entry_arr[0].rsplit("_", 1))
                    bin_fasta_list = outfolder + "/bins/" + newbinid + ".fasta"
                    scaffolds_dict[scaffold] = bin_fasta_list
    print("summarizing sample")
    outfile = outfolder + "/summary_all.tsv"
    outfile2 = outfolder + "/../summary_all.tsv"
    o2 = open(outfile2, "a+")
    o = open(outfile, "w+")
    o.write("Bin\tCount\tavg_cg\tmin_gc\t_max_gc\tavg_cov\tmin_cov\tmax_cov\ttotal_len\tmin_len\tmax_len"
            "\tcompleteness\tcontamination\n")
    for i in bins:
        sample_ar = i.split("_")
        samplename = sample_ar[0] + "_" + sample_ar[1]
        avg_gc = round(bins[i]["gcsum"] / bins[i]["count"])
        avg_cov = round(bins[i]["covsum"] / bins[i]["count"])
        len_in_mb = round(bins[i]["lensum"] / 1000000, 2)
        name_gc = sample_ar[-2]
        name_cov = sample_ar[-1]
        sum_key = str(samplename) + "_" + str(name_gc) + "_" + str(name_cov)
        completness = 0
        contamination = 1000
        if sum_key in cur_sum:
            completness = cur_sum[sum_key]["compl"]
        else:
            print(sum_key + " not found in curated list more info:" + i)
        if sum_key in cur_sum:
            contamination = cur_sum[sum_key]["contam"]

        line_str = (str(i) + "\t" + str(bins[i]["count"]) + "\t" +
                    str(avg_gc) + "\t" + str(bins[i]["mingc"]) + "\t" + str(bins[i]["maxgc"]) + "\t" +
                    str(avg_cov) + "\t" + str(bins[i]["mincov"]) + "\t" + str(bins[i]["maxcov"]) + "\t" +
                    str(len_in_mb) + "\t" + str(bins[i]["minlen"]) + "\t" + str(bins[i]["maxlen"]) + "\t" +
                    str(completness) + "\t" + str(contamination) + "\n")
        o.write(line_str)
        o2.write(line_str)
    o.close()
    o2.close()
    source_fasta = f.fasta_folder + sample + "_min1000.fasta"
    extract_from_fasta(scaffolds_dict, source_fasta)  # define here wether to use extraction
    print("oldbins: " + str(count_classified) + "\tnon_binned: " + str(
        count_unclassified) + "\trefined_binned: " + str(count_binned) + "\tNo_of_Bins: " + str(len(bins)))


def extract_from_fasta(scaff_dict, source_fasta):
    print("extracting scaffolds to bins")
    with open(source_fasta) as source:
        for line in source:
            header = line.strip()
            seq = next(source).strip()
            m_temp = (str(header).replace(">", "")).strip()
            if m_temp in scaff_dict:
                x = open(scaff_dict[m_temp], "a+")
                x.write(header + "\n")
                x.write(seq + "\n")
                x.close()
                # command1 = "echo \\>" + m_temp + " >>" + scaff_dict[m_temp]
                # command2 = "echo " + seq + " >>" + scaff_dict[m_temp]
                #os.system(command1)
                #os.system(command2)


def curation_mem(cf):
    summary_dict = {}
    with open(cf, "r") as curfile:
        next(curfile)
        for line in curfile:
            line_arr = line.split("\t")
            sample_name = line_arr[0]
            avg_gc = line_arr[8]
            bname = str(line_arr[7]).split("_")
            avg_cov = bname[-1]
            smu_dict_key = sample_name + "_" + avg_gc + "_" + avg_cov
            compl = line_arr[9]
            contam = line_arr[10]
            summary_dict[smu_dict_key] = dict([("compl", compl), ("contam", contam)])
    return summary_dict


def initiate_outfolders(outfolder):
    # os.system("rm -r "+outfolder)
    Path(outfolder).mkdir(parents=True, exist_ok=True)
    outfile2 = outfolder + "/summary_all.tsv"
    o = open(outfile2, "w+")
    o.write("Bin\tCount\tavg_cg\tmin_gc\t_max_gc\tavg_cov\tmin_cov\tmax_cov\ttotal_len\tmin_len\tmax_len"
            "\tcompleteness\tcontamination\n")
    o.close()


def main():
    global cur_sum
    cur_sum = curation_mem(curation_file)
    outfolder = f.outfolder
    initiate_outfolders(outfolder)
    for samp in files_dict:
        sample = samp
        file = files_dict[samp]
        standardise(sample, file,outfolder)


if __name__ == "__main__":
    main()
