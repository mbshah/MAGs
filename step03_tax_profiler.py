#!/usr/bin/env python
import os
from pathlib import Path
import re
import config as config
import sys
import subprocess

# input parameters
list_file = config.outfolder+"/dRep_workfolder/post_cluster_data.csv"
taxonomy_col = 27  # zero index # MAINLY if you want to use checkm output, if using dastool+ubin output, sample and taxonomy definintion in function profilize
usecheckm = True
compl_col = config.compl_col  # zero index
min_compl = config.min_compl
max_compl = config.max_compl
contam_col = config.contam_col
min_contam = config.min_contam
max_contam = config.max_contam
taxaofinterest = ()  # comma seperated list of taxa of interest
outdir = config.outfolder + "/profiles_checkm_90_10"  # where output will be saved
if len(sys.argv) > 1: outdir = sys.argv[1]
taxdumpdir = config.taxdumpdir  # must download latest tax dump and extract it and provide path here

# empty variables
lineage_db = {}
names = {}
nodes = {}
names_arr_dict = {}
profiles = {"l1": {}, "l2": {}, "l3": {}, "l4": {}, "l5": {}, "l6": {}, "l7": {}}
samples = []

# default_universal variable {DO NOT CHANGE IF U DONT KNOW WHAT YOU ARE DOING
lin_req_dict = {"Bacteria": ["superkingdom", "phylum", "class", "order", "family", "genus", "species"],
                "Viruses": ["superkingdom", "order", "family", "genus", "species"],
                "Eukaryota": ["superkingdom", "kingdom", "phylum", "class", "order", "family", "species"],
                "Archaea": ["superkingdom", "phylum", "class", "order", "family", "genus", "species"],
                "metagenomes": ["no rank", "no rank", "species"]}


# defining Functions


def tax_initiate(taxdump_dir):  # function to copy NCBI taxonomy to memory
    taxidlineage = taxdump_dir + "/" + "taxidlineage.dmp"
    with open(taxidlineage, "r") as lf:
        for entry in lf:
            entry_arr = entry.split("|")
            global lineage_db
            lineage_db[entry_arr[0].strip()] = entry_arr[1].strip()
    print("loaded " + str(len(lineage_db)) + " lineage entries into memory")

    namesfile = taxdump_dir + "/" + "names.dmp"
    with open(namesfile, "r") as lf:
        for entry in lf:
            entry_arr = entry.split("|")
            if entry_arr[3].strip().startswith("scientific name"):
                global names
                names[entry_arr[0].strip()] = entry_arr[1].strip()
                global names_arr_dict

                if entry_arr[1].strip() in names_arr_dict:
                    temp = [names_arr_dict[entry_arr[1].strip()][0], entry_arr[0].strip()]
                    names_arr_dict[entry_arr[1].strip()] = temp
                else:
                    temp = [entry_arr[0].strip()]
                    names_arr_dict[entry_arr[1].strip()] = temp
    print("loaded " + str(len(names)) + " names entries into memory")

    nodesfile = taxdump_dir + "/" + "nodes.dmp"
    with open(nodesfile, "r") as lf:
        for entry in lf:
            entry_arr = entry.split("|")
            global nodes
            nodes[entry_arr[0].strip()] = entry_arr[2].strip()
    print("loaded " + str(len(names)) + " nodes entries into memory")


def lineage(taxname, node):
    try:
        taxids = names_arr_dict[taxname]
    except KeyError:
        if taxname.endswith("2"): taxname = re.sub('2$', "", taxname)
        if taxname == "unclassified": taxname = "Bacteria"
        try:
            taxids = names_arr_dict[taxname]
        except KeyError:
            return "SKIP_IT"
    taxid = taxids[0]
    if len(taxids) > 1:
        for x in taxids:
            if nodes[x][0] == node:
                taxid = x
    lineage_str = ""
    if taxid in lineage_db:
        lineage_str = lineage_db[taxid]
    else:
        print("error " + str(taxid) + " not found in taxid, trying a fix by finding in delnode (return 2) or merged(return newtaxid)")
        taxid=find_missing_tax(taxid)
    lineage_str = str(lineage_str + " " + str(taxid)).strip()
    if lineage_str.startswith("\s") or lineage_str == "": print("note on " + str(taxid))
    return lineage_str


def filter_check(upper, lower, value):
    if float(lower) <= float(value) <= float(upper):
        return True
    else:
        return False


def filterok(a, b):
    if a == True and b == True:
        return True
    else:
        return False


def taxonomy_refiner(tax, db):
    taxID_lin = str(tax).split(" ")
    try:
        mysk = names[taxID_lin[1]]
    except:
        try:
            mysk=names[taxID_lin[0]]
        except:
            print(str(names[tax]) + " in " + db+ " is not a regular superkingdom")
            exit()
    if mysk in lin_req_dict:pass
    else: mysk=names[taxID_lin[0]]
    named_lineage = {}
    #print(taxID_lin)
    for levelID in taxID_lin:
        if nodes[levelID] in lin_req_dict[mysk]:
            named_lineage[nodes[levelID]] = names[levelID]
    prevlevel = ""
    named_lin_arr = []
    for label in lin_req_dict[mysk]:
        if label in named_lineage:
            lin_level2 = named_lineage[label]
            lin_level = named_lineage[label]
        else:
            lin_level = "other_" + prevlevel
            lin_level2 = "other_" + prevlevel + "_" + label
            parts = str(lin_level2).partition("other")
            lin_level2 = parts[0] + parts[1] + parts[2].replace("other_", "")
        # print(lin_level)
        named_lin_arr.append(lin_level2)
        prevlevel = lin_level
    named_lineage_str = ";".join(named_lin_arr)
    for x in taxaofinterest:
        if named_lineage_str.__contains__(x): print(x + "\t" + db)
    return named_lineage_str


def add_to_profile(sample, lineage_str):
    lin = lineage_str.split(";")
    lin_l1 = lin[0]
    lin_l2 = ";".join(lin[0:2])
    lin_l3 = ";".join(lin[0:3])
    lin_l4 = ";".join(lin[0:4])
    lin_l5 = ";".join(lin[0:5])
    lin_l6 = ";".join(lin[0:6])
    lin_l7 = ";".join(lin[0:7])
    x1 = [lin_l1, lin_l2, lin_l3, lin_l4, lin_l5, lin_l6, lin_l7]
    for x in range(1, 8):
        lin = x1[x - 1]
        lev = "l" + str(x)
        if lin in profiles[lev]:
            if sample in profiles[lev][lin]:
                profiles[lev][lin][sample] = profiles[lev][lin][sample] + 1
            else:
                profiles[lev][lin][sample] = 1
        else:
            profiles[lev][lin] = {}
            profiles[lev][lin][sample] = 1


def profilize(infile, col):
    with open(infile) as tsv:
        header = next(tsv)
        for line in tsv:
            entry = line.split("\t")
            bin_arr = entry[0].split("_")
            sample = "_".join(bin_arr[0:2])
            tax = " ".join(bin_arr[2:-2])  # IF NOT USING CHECKM please change this
            lin_node = "c"
            if sample not in samples: samples.append(sample)
            #isfilterok = config.filterok(config.filter_check(max_compl, min_compl, entry[compl_col]),
            #                             config.filter_check(max_contam, min_contam, entry[contam_col]))
            #print(entry[contam_col])
            isfilterok = True
            if isfilterok:
                if usecheckm:
                    print(entry[col])
                    lin_node = entry[col][0]
                    tax = re.sub('\(.*\)$', "", entry[col]).strip()
                    tax = re.sub('\w__', "", tax)
                entrylineage = lineage(tax, lin_node)
                if entrylineage == "SKIP_IT":
                    pass
                else:
                    refined_lineage = taxonomy_refiner(entrylineage, entry[0])
                    add_to_profile(sample, refined_lineage)
    for z in (range(1, 8)):
        lev = "l" + str(z)
        outfile = outdir + "/" + lev + "_" + str(lin_req_dict["Bacteria"][z - 1]) + ".tsv"
        headers = "taxa/sample\t" + "\t".join(sorted(samples)) + "\n"
        out = open(outfile, "w")
        out.write(headers)
        for y in profiles[lev]:
            out.write(y + "\t")
            for x in sorted(samples):
                count = 0
                if x in profiles[lev][y]: count = profiles[lev][y][x]
                out.write(str(count) + "\t")
            out.write("\n")


def find_missing_tax(tax):
    command=["grep","-r","^"+tax,taxdumpdir]
    #command= f"grep {tax} {taxdumpdir}/*"
    proc=subprocess.run(command,capture_output=True)
    result=proc.stdout.decode("utf-8").strip()

    if result.__contains__("delnode"):
        return str(2)
    elif result.__contains__("merged"):
        newtax=result.split(":")[1].split("|")[1].replace("\s","")
        #print(tax, result)
        return newtax
    else:
        return str(2)


def main():
    print("min_compl " + str(min_compl))
    tax_initiate(taxdumpdir)
    os.system("rm -r " + outdir)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    profilize(list_file, taxonomy_col)


if __name__ == "__main__":
    main()
