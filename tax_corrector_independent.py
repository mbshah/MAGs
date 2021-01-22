import re
import subprocess
import sys

table_file="/home/hb0358/PycharmProjects/mbs_general/some_test_files/filtered_blast_table.csv"
table_seperator=","
tax_col=13
count_col_start=14
count_col_end=-1  #use -1 for last
taxdumpdir = "/home/hb0358/PycharmProjects/mbs_general/ancilary/new_taxdump"  # must download latest tax dump and extract it and provide path here
taxaofinterest = ()  # comma seperated list of taxa of interest
output_type="original"  #for single table use "single", for level_wise summary use "summary"
mode="profile" #if input is a sample vs tax profile use "profile", for assignemnt output that needs counting use "count"

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
                "metagenomes": ["no rank", "no rank", "species"],
                "unclassified entries": ["no rank", "no rank", "species"]}


# defining Functions


def tax_initiate(taxdump_dir):  # function to copy NCBI taxonomy to memory
    taxidlineage = taxdump_dir + "/" + "taxidlineage.dmp"
    with open(taxidlineage, "r") as lf:
        for entry in lf:
            entry_arr = entry.split("|")
            global lineage_db
            lineage_db[entry_arr[0].strip()] = entry_arr[1].strip()
    sys.stderr.write("loaded " + str(len(lineage_db)) + " lineage entries into memory\n")
    sys.stderr.flush()

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
    sys.stderr.write("loaded " + str(len(names)) + " names entries into memory\n")
    sys.stderr.flush()

    nodesfile = taxdump_dir + "/" + "nodes.dmp"
    with open(nodesfile, "r") as lf:
        for entry in lf:
            entry_arr = entry.split("|")
            global nodes
            nodes[entry_arr[0].strip()] = entry_arr[2].strip()
    sys.stderr.write("loaded " + str(len(names)) + " nodes entries into memory\n")
    sys.stderr.flush()


def taxonomy_refiner(tax, db):
    taxID_lin = str(tax).split(" ")
    try:
        mysk = names[taxID_lin[1]]
    except:
        try:
            mysk=names[taxID_lin[0]]
        except:
            sys.stderr.write(str(names[tax]) + " in " + db+ " is not a regular superkingdom")
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


def lineage(taxname):
    try:
        taxids = names_arr_dict[taxname]
    except KeyError:
        if taxname.endswith("2"): taxname = re.sub('2$', "", taxname)
        if taxname.__contains__("unclassified"):taxname = "uncultured bacterium"

        try:
            taxids = names_arr_dict[taxname]
        except KeyError:
            return "SKIP_IT"
    taxid = taxids[0]
    lineage_str = ""
    if taxid in lineage_db:
        lineage_str = lineage_db[taxid]
    else:
        sys.stderr.write("error " + str(taxid) + " not found in taxid, trying a fix by finding in delnode (return 2) or merged(return newtaxid)")
        taxid=find_missing_tax(taxid)
    lineage_str = str(lineage_str + " " + str(taxid)).strip()
    if lineage_str.startswith("\s") or lineage_str == "": print("note on " + str(taxid))
    return lineage_str


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


def read_and_process_table():
    with open(table_file,"r") as table:
        header=next(table)
        global samples
        if mode=="profile":
            samples=header[count_col_start:count_col_end]
        if output_type == "original":
            output_file = table_file + ("_out.tsv")
            print(header)
        for line in table:
            entry=line.strip().split(table_seperator)
            lineage_given=entry[tax_col-1].split(";")
            taxname=lineage_given[-1].strip()
            taxid=names_arr_dict[taxname]
            entrylineage=lineage(taxname)
            if entrylineage == "SKIP_IT":
                pass
            else:
                refined_lineage = taxonomy_refiner(entrylineage, entry[0])
                if mode=="count":
                   add_to_profile(sample, refined_lineage)
                if output_type=="original":
                    output_file=table_file+("_out.tsv")
                    entry[tax_col-1]=refined_lineage
                    new_line=table_seperator.join(entry)
                    print(new_line)

tax_initiate(taxdumpdir)
#print(lineage_db["1673428"])
read_and_process_table()