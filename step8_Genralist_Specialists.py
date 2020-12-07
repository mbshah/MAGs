import step3_tax_profiler as tp
import config as config
from pathlib import Path

kraken_oufolder="../WMG/Kraken2_NT_Scaffolds/"
kraken_profile_folder="../WMG/Kraken2_NT_Scaffolds_profiles/"
tax_col=2 #zero_based
ignoretax=(1,131567,0)
sampleClassificationCounter={}
samples=sorted(config.files_dict.keys())
profiles={}


def kraken_profilize():
    tp.tax_initiate(config.taxdumpdir)
    for sample in samples:
        sampleClassificationCounter[sample]={}
        sampleClassificationCounter[sample]["classified"]=0
        sampleClassificationCounter[sample]["unclassified"]=0
        kraken_out_file=kraken_oufolder+"/"+sample+"_min1000.out"
        print("\rprocessing: "+kraken_out_file,end="")
        with open(kraken_out_file) as raw:
            for line in raw:
                entry=line.split("\t")
                if entry[0]=="U":
                    sampleClassificationCounter[sample]["unclassified"]=sampleClassificationCounter[sample]["unclassified"]+1
                else:
                    tax=entry[tax_col]
                    if int(tax)in ignoretax:
                        sampleClassificationCounter[sample]["unclassified"] = sampleClassificationCounter[sample][
                                                                                  "unclassified"] + 1
                    else:
                        if tax in tp.names:
                            pass
                        else:
                            ntax = tp.find_missing_tax(tax)
                            tax = ntax.strip()
                        fullLineage=str(tp.lineage_db[tax]).split(" ")
                        fullLineage.append(tax)
                        if fullLineage[0]=="":fullLineage.pop(0)
                        root_lineage=fullLineage[0].strip()
                        #print(fullLineage)
                        if int(root_lineage) in (2787854,2787823):
                            sampleClassificationCounter[sample]["unclassified"] = sampleClassificationCounter[sample][
                                                                                      "unclassified"] + 1
                        else:
                            sampleClassificationCounter[sample]["classified"] = sampleClassificationCounter[sample][
                                                                                    "classified"] + 1
                            entryLineage=tp.lineage(tp.names[tax],tp.nodes[tax][0])
                            refined_lineage = tp.taxonomy_refiner(entryLineage, sample)
                            add_to_profiles(sample,refined_lineage)


def add_to_profiles(sample,lineage):
    lineage_list=lineage.split(";")
    superkingdom=lineage_list[0]
    if superkingdom in profiles:pass
    else: profiles[superkingdom]={}
    current_org={}
    i=1
    for level in tp.lin_req_dict[superkingdom]:
        current_org[level]=";".join(lineage_list[0:i])
        i=i+1
        if level in profiles[superkingdom]:pass
        else: profiles[superkingdom][level]={}
        if sample in profiles[superkingdom][level]:pass
        else: profiles[superkingdom][level][sample]={}
        if current_org[level] in profiles[superkingdom][level][sample]:
            profiles[superkingdom][level][sample][current_org[level]]=profiles[superkingdom][level][sample][current_org[level]]+1
        else:
            profiles[superkingdom][level][sample][current_org[level]]=1


def profile_to_files():
    superkingdom_counts={}
    Path(kraken_profile_folder).mkdir(parents=True, exist_ok=True)
    for superkingdom in tp.lin_req_dict:
        supkm_folder=kraken_profile_folder+"/"+superkingdom
        if superkingdom in profiles:
            Path(supkm_folder).mkdir(parents=True, exist_ok=True)
            for lev in tp.lin_req_dict[superkingdom]:
                if lev=="superkingdom":
                    #print(profiles[superkingdom][lev])
                    for samp in profiles[superkingdom][lev]:
                        if samp in superkingdom_counts:pass
                        else:superkingdom_counts[samp]={}
                        for sk in profiles[superkingdom][lev][samp]:
                            superkingdom_counts[samp][sk]=profiles[superkingdom][lev][samp][sk]
                else:
                    outfile = supkm_folder + "/" + lev +".tsv"
                    headers = "taxa/sample\t" + "\t".join(sorted(samples)) + "\n"
                    out = open(outfile, "w")
                    out.write(headers)
                    orgs=set()
                    for samp in profiles[superkingdom][lev]:
                        orgs.update(profiles[superkingdom][lev][samp].keys())
                    for organism in orgs:
                        out.write(organism+"\t")
                        for samp in profiles[superkingdom][lev]:
                            count=0
                            if organism in profiles[superkingdom][lev][samp]:count=profiles[superkingdom][lev][samp][organism]
                            out.write(str(count)+"\t")
                        out.write("\n")
                    out.close()
    summary_file=kraken_profile_folder+"/"+"Counts_overview.tsv"
    sum_out=open(summary_file,"w")
    sum_header="data\t"+"\t".join(profiles)+"\tUnclassified\tclassified\t%classified+total_counts"
    sum_out.write(sum_header+"\n")
    for sample in superkingdom_counts:
        sum_out.write(sample+"\t")
        for supkm in profiles:
            sum_out.write(str(superkingdom_counts[sample][supkm])+"\t")
        total_scaffolds=(sampleClassificationCounter[sample]["classified"]+sampleClassificationCounter[sample]["unclassified"])
        percenclassified=(sampleClassificationCounter[sample]["classified"]/total_scaffolds)*100
        sum_out.write(str(sampleClassificationCounter[sample]["unclassified"])+
                      "\t"+str(sampleClassificationCounter[sample]["classified"])+"\t"+str(percenclassified)+
                      "\t"+str(total_scaffolds)+"\n")


kraken_profilize()
profile_to_files()
#print(profiles["Bacteria"]["phylum"]['EULs_S301BU'].keys())
