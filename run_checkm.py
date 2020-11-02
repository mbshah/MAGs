#!/usr/bin/env python
import os
import files as f

files_dict = dict(f.files_dict)


def checkm_summary(outfolder, sample):
    print("running checkm for " + sample)
    checkm_out = outfolder + "/" + sample + "/" + "/checkm_out/"
    checkm_in = outfolder + "/" + sample + "/bins"
    checkm_tsv = outfolder + "/" + sample + "/checkm.tsv"
    os.system(f"checkm lineage_wf -t 8 -x fasta --tab_table -q {checkm_in} {checkm_out} >{checkm_tsv}")


def walklevel(some_dir, level=1):
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]


def checkm_stand_merge(outfolder):
    maxlevel = 1
    fullsummary = outfolder + "/all_sample_summary.tsv"
    o2_flag = 1
    for (dirpath, dirnames, filenames) in walklevel(outfolder, maxlevel):
        path_array = str(dirpath).split("/")
        if len(path_array) == 2 and path_array[1].startswith("EUL"):
            sample = path_array[1]
            checkm_dict = {}
            new_outfile = dirpath + "/" + sample + "_summary.tsv"
            with open(dirpath + "/checkm.tsv") as c:
                header = next(c)

                for line in c:
                    line_arr = line.strip().split("\t")
                    bin = line_arr[0]
                    newline = "\t".join(line_arr[1:])
                    checkm_dict[bin] = newline
            with open(dirpath + "/summary_all.tsv") as c:
                temph = next(c)
                header = str(temph).strip() + header.strip().replace("Bin Id", "")
                o = open(new_outfile, "w+")
                o.write(header + "\n")
                o2 = open(fullsummary, "a+")
                if o2_flag == 1: o2.write(header + "\n");o2_flag = 0
                for line in c:
                    line_arr = line.strip().split("\t")
                    bin = line_arr[0]
                    newline = "\t".join(line_arr[1:]) + "\t" + checkm_dict[bin]
                    o.write(bin + "\t" + newline + "\n")
                    o2.write(bin + "\t" + newline + "\n")
                o.close()
                o2.close()


def main():
    outfolder = f.outfolder
    for samp in files_dict:
        checkm_summary(outfolder, samp)
        checkm_stand_merge(outfolder)


if __name__ == "__main__":
    main()
