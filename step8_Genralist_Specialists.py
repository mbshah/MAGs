import operator
import step3_tax_profiler as tp
import config as config
from pathlib import Path
import collections
import json
import os

kraken_oufolder="../WMG/Kraken2_NT_Scaffolds/"
kraken_profile_folder="../WMG/Kraken2_NT_Scaffolds_profiles2/"
tax_col=2 #zero_based
ignoretax=(1,131567,0)
sampleClassificationCounter={}
samples=sorted(config.files_dict.keys())
profiles={}
scafid_lineage={}
scafLin_jsonfile=kraken_profile_folder+"scafLin.map"
lin_class_jsonfile=kraken_profile_folder+"linclass.map"
scafLinMAG_jasonfile=kraken_profile_folder+"scafLinMAG.map"
mag_lin_jsonfile=kraken_profile_folder+"MAGlin.map"
linScaf_jsonfile=kraken_profile_folder+"linScaf.map"

lineage_scafs={}
gen_spec_dict={}
segregate_scaffolds=True
lineage_mag={}
mag_lineage={}


def kraken_profilize():
    for sample in samples:
        sampleClassificationCounter[sample]={}
        sampleClassificationCounter[sample]["classified"]=0
        sampleClassificationCounter[sample]["unclassified"]=0
        kraken_out_file=kraken_oufolder+"/"+sample+"_min1000.out"
        print("\rprocessing: "+kraken_out_file,end="")
        with open(kraken_out_file) as raw:
            for line in raw:
                entry=line.split("\t")
                scaf_id=entry[1].strip()
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
                            entryLineage=tp.lineage(tp.names[tax],"x")
                            refined_lineage = str(tp.taxonomy_refiner(entryLineage, sample)).strip().replace(" ","_").replace("'","_")
                            add_to_profiles(sample,refined_lineage)
                            scafid_lineage[scaf_id]=(refined_lineage,tax)
                            if refined_lineage in lineage_scafs:
                                lineage_scafs[refined_lineage]=lineage_scafs[refined_lineage]+";"+scaf_id
                            else:
                                lineage_scafs[refined_lineage]=scaf_id
    linScaf_map=json.dumps(lineage_scafs)
    f=open(linScaf_jsonfile,"w")
    f.write(linScaf_map)
    f.close()
    scafLin_map = json.dumps(scafid_lineage)
    f = open(scafLin_jsonfile,"w")
    f.write(scafLin_map)
    f.close()


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
            prefix=1
            for lev in tp.lin_req_dict[superkingdom]:
                if lev=="superkingdom":
                    #print(profiles[superkingdom][lev])
                    for samp in profiles[superkingdom][lev]:
                        if samp in superkingdom_counts:pass
                        else:superkingdom_counts[samp]={}
                        for sk in profiles[superkingdom][lev][samp]:
                            superkingdom_counts[samp][sk]=profiles[superkingdom][lev][samp][sk]
                else:
                    outfile = supkm_folder + "/" + "l"+str(prefix)+"_"+lev + ".tsv"
                    prefix=prefix+1
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


def log_otu_abundance_plot():
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import weibull_min
    normalized_file=kraken_profile_folder+"/Bacteria/l6_species_normalized_deseq2.tsv"
    count_array=[]
    exists_array=[]
    with open(normalized_file,"r") as nfile:
        header=next(nfile)
        for line in nfile:
            entry=line.split("\t")
            otu_name=entry.pop(0)
            total_count=0
            exists_count=0
            for value in entry:
                value=float(value)
                if value>0:
                    total_count=total_count+value
                    exists_count=exists_count+1
            avg_count=total_count/exists_count
            log_total_count=np.log(total_count)
            log_avg_count=np.log(avg_count)
            log_exists_count=np.log(exists_count)
            count_array.append(log_avg_count)
            exists_array.append(log_exists_count)
            #print(otu_name+":\t"+str(log_exists_count)+"\t"+str(log_total_count)+"\t"+str(log_avg_count)+"\n")
    fig,ax=plt.subplots()
    ax.scatter(count_array, exists_array,marker=".",color="red")
    plt.xlabel("abundance")
    plt.ylabel("log of occurrence")
    plt.gca().invert_xaxis()
    #fig.savefig(kraken_profile_folder+"/Bacteria/l6_species_normalized_deseq2_pyplot.png")
    plt.show()

def run_minimap():
    infilelist=config.files_dict.keys()
    for sample in infilelist:
        infile=config.fasta_folder+sample+"_min1000.fasta"
        ref=config.outfolder+"clusters_dastool_90_10_average/cluster_reps/all.fa"
        outfolder=config.outfolder+"scaffold_to_cluster_map"
        Path(outfolder).mkdir(parents=True, exist_ok=True)
        outfile=outfolder+f"/{sample}_scaffold_alignment.paf"
        threads=18
        cmd=f"minimap2 -x asm5 {ref} {infile} -o {outfile} -t {threads} --secondary=no"
        #with secondary have 1019537 scafolds in awk '{print $1}' * |sort |uniq |wc -l
        print(cmd)
        os.system(cmd)


def read_minimap_paf():
    scafid_lineage = json.load(open(scafLin_jsonfile))
    global lineage_mag
    infilelist = config.files_dict.keys()
    scafid_lineage2={}
    scafid_lineage3={}
    unique_magnames=[]
    for sample in infilelist:
        infile = config.outfolder + "scaffold_to_cluster_map" + f"/{sample}_scaffold_alignment.paf"
        with open(infile,"r")as paffile:
            for line in paffile:
                entry=line.split()
                scaffold=entry[0]
                mag=entry[5].split("|")[1]
                scaf_len10perc=float(entry[1])-(float(entry[1])*0.1)
                aln_len=entry[9]
                lineage=[]
                if scaffold in scafid_lineage:
                    lineage=scafid_lineage[scaffold]
                    if len(lineage)<=2:
                        mag_dict = {mag:0}
                        lineage.append(mag_dict)
                        lineage.append(scaf_len10perc)
                    if mag in lineage[2]:
                        lineage[2][mag]=int(lineage[2][mag])+int(aln_len)
                    else:
                        lineage[2][mag]=aln_len
                    scafid_lineage2[scaffold]=lineage
    for scaffold in scafid_lineage2:
        mags=scafid_lineage2[scaffold][2]
        mags2={}
        for mag in mags:
            if float(mags[mag])>float(scafid_lineage2[scaffold][3]):
                mags2[mag]=mags[mag]
                if mag in unique_magnames:pass
                else: unique_magnames.append(mag)
        if len(mags2)>=1:
            lineage=[scafid_lineage2[scaffold][0],scafid_lineage2[scaffold][1],mags2]
            scafid_lineage3[scaffold]=lineage
    for scaffold in scafid_lineage3:
        lineage = scafid_lineage3[scaffold][0]
        tax_no = scafid_lineage3[scaffold][1]
        mags = scafid_lineage3[scaffold][2]
        print(f"{scaffold}has multiple MAGs\t{mags}") if len(mags)>1 else print("",end="")
        mag=list(mags.keys())[0]
        if lineage in lineage_mag:
            lineage_mag[lineage].append(mag)
            lineage_mag[lineage]=list(lineage_mag[lineage])
        else:
            lineage_mag[lineage]=[mag]
        if mag in mag_lineage:
            mag_lineage[mag].append(lineage)
        else:
            mag_lineage[mag]=[lineage]
    #for scaffold in scafid_lineage3:
    #    print (f"{scaffold}\t{scafid_lineage3[scaffold]}")
    #for lineage in lineage_mag:
    #    print(f"{lineage}\t{len(set(lineage_mag[lineage]))}\t{collections.Counter(lineage_mag[lineage])}")
    for mag in mag_lineage:
        #print(mag)
        #final_lineage=""
        if len(mag_lineage[mag])>1:
            final_lineage=lca_finder(mag_lineage[mag])
        else:
            final_lineage=mag_lineage[mag][0]
        mag_lineage[mag]=final_lineage
        #print(f"{mag}\t{final_lineage}\t{len(mag_lineage[mag])}\t{collections.Counter(mag_lineage[mag]).most_common(1)}")
    #print(lineage_mag["Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;Acidovorax;Acidovorax_sp._T1"])
    scafLinMAG_map = json.dumps(scafid_lineage3)
    f = open(scafLinMAG_jasonfile, "w")
    f.write(scafLinMAG_map)
    f.close()
    mag_lineage_map = json.dumps(mag_lineage)
    f = open(mag_lin_jsonfile, "w")
    f.write(mag_lineage_map)
    f.close()


def post_gen_spec_classify():
    classify_file=kraken_profile_folder+"/Bacteria/Genralists_Specialists/gen_spec1000.tsv"
    classify_file_out = kraken_profile_folder + "/Bacteria/Genralists_Specialists/gen_specret.tsv"
    gen_fasta_folder=kraken_profile_folder+"/Bacteria/Genralists_Specialists/GEN/fasta/"
    Path(gen_fasta_folder).mkdir(parents=True, exist_ok=True)
    spec_fasta_folder = kraken_profile_folder + "/Bacteria/Genralists_Specialists/SPEC/fasta/"
    Path(spec_fasta_folder).mkdir(parents=True, exist_ok=True)
    ns_fasta_folder = kraken_profile_folder + "/Bacteria/Genralists_Specialists/NS/fasta/"
    Path(ns_fasta_folder).mkdir(parents=True, exist_ok=True)

    scafid_lineage = json.load(open(scafLin_jsonfile))
    lineage_scafs = json.load(open(linScaf_jsonfile))
    mag_lineage = json.load(open(mag_lin_jsonfile))
    #print(scafid_lineage.keys())
    with open (classify_file) as cf:
        header=next(cf)
        for line in cf:
            entry=line.strip().split("\t")
            mag=""
            lineage=entry[0]
            if lineage in mag_lineage.values():
                entry[5]=entry[5]+"_m"
            gen_spec_dict[lineage]=entry[5]
            scaf=lineage_scafs[lineage]
            scafs=scaf.split(";")
            #print(f"{scafid_lineage}")
            tax = scafid_lineage[scafs[0]][1]
            #print(f"{scafs[0]}\t{tax}")
            out_file=open(classify_file_out,"w")
            out_file.write("\tRes\tColor\n")
            for lineage in gen_spec_dict:
                colour="#44AA99"
                if gen_spec_dict[lineage]=="GENERALIST":
                    colour="#332288"
                if gen_spec_dict[lineage]=="SPECIALIST":
                    colour="#117733"
                if gen_spec_dict[lineage]=="NON SIGNIFICANT":
                    colour="#44AA99"
                if gen_spec_dict[lineage] == "GENERALIST_m":
                    colour = "#88CCEE"
                if gen_spec_dict[lineage] == "SPECIALIST_m":
                    colour = "#DDCC77"
                if gen_spec_dict[lineage] == "NON SIGNIFICANT_m":
                    colour = "#CC6677"
                line=f"{lineage}\t{gen_spec_dict[lineage]}\t{colour}\n"
                out_file.write(line)
            out_file.close()
            if segregate_scaffolds==True:
                if entry[5].startswith("SPEC"):
                    tax_folder=spec_fasta_folder+"\""+tax+"\""+"/"
                    t2=spec_fasta_folder+tax+"/"
                elif entry[5].startswith("GEN"):
                    tax_folder=gen_fasta_folder+"\""+tax+"\""+"/"
                    t2=gen_fasta_folder+tax+"/"
                else:
                    tax_folder=ns_fasta_folder+"\""+tax+"\""+"/"
                    t2=ns_fasta_folder+tax+"/"
                Path(t2).mkdir(parents=True, exist_ok=True)
                for i in scafs:
                    sample="_".join(str(i).split("_")[0:2])
                    sample_fasta=config.fasta_folder+"/"+sample
                    cmd="grep -A1 "+i+" "+sample_fasta+"_min1000.fasta >>"+tax_folder+"org.fasta"
                    #cat bash_commands.sh |parallel -j 10 {}
                    # ##dosent work buffer issue, creates fasta files with incorrect formats.
                    # ###instead split into 15 chunks and ran 15 chunks in parallel
                    os.system(cmd)
                    #print(cmd)
            #print(tax+"\t"+scaf)
    gen_spec_map=json.dumps(gen_spec_dict)
    f = open(lin_class_jsonfile, "w")
    f.write(gen_spec_map)
    f.close()


def lca_finder(lineage_list):
    overdict={}
    parent_dict={}
    spkg=[]
    superkingdom=""
    return_lineage=[]
    for lineage in lineage_list:
        lineagearray=list(map(str.strip,lineage.split(";")))
        level=0
        if lineagearray[level].__contains__("other"):pass
        else:spkg.append(lineagearray[level])
    superkingdom=collections.Counter(spkg).most_common(1)[0][0]
    req_taxa = tp.lin_req_dict[superkingdom]
    for lineage in lineage_list:
        lineagearray = list(map(str.strip, lineage.split(";")))
        i = 0
        if lineagearray[0]==superkingdom:
            for taxa_level in req_taxa:
                #if lineagearray[i].__contains__("other"):continue
                lineagearray[i]=f"{taxa_level[0:2]}_{lineagearray[i]}"
                if taxa_level in overdict:
                    if taxa_level == "superkingdom":
                        parent_dict[lineagearray[i]]="NA"
                    else:
                        parent_dict[lineagearray[i]]=lineagearray[i-1]
                    #print(f"{lineagearray[i-1]}\t{lineagearray[i]}")
                    overdict[taxa_level].append(lineagearray[i])
                else:
                    if taxa_level == "superkingdom":
                        parent_dict[lineagearray[i]]="NA"
                    else:
                        parent_dict[lineagearray[i]]=lineagearray[i-1]
                    overdict[taxa_level]=[lineagearray[i]]
                i=i+1
    overdict2={}
    for req_taxa_lvl in req_taxa:
        overdict2[req_taxa_lvl]=dict(collections.Counter(overdict[req_taxa_lvl]))
        #print(f"{req_taxa_lvl}\t{overdict2[req_taxa_lvl]}")
    weighted_overdict=overdict2
    #print(req_taxa)
    #
    #print(sorted(parent_dict))
    for i in reversed(range(0,len(req_taxa))):
        level=req_taxa[i]
        if level == "superkingdom": continue
        for taxa in overdict2[level]:
            if taxa == superkingdom:pass
            else:
                parent=parent_dict[taxa]
                #print(f"----{parent}--{taxa}----")
                #parent_node=tp.nodes[[-1]]
                weight=overdict2[req_taxa[i-1]][parent]
                old_value=overdict2[req_taxa[i]][taxa]
                new_value=old_value+weight
                weighted_overdict[req_taxa[i]][taxa]=new_value
                #print (f"{req_taxa[i-1]}\t{parent}\t\t{taxa}\t{old_value}\t{weight}\t{new_value}")
    for i in reversed(range(0,len(req_taxa))):
        level=req_taxa[i]
        if level=="superkingdom":continue
        most_abundant = max(weighted_overdict[level].items(), key=operator.itemgetter(1))[0]
        compare_value = float(weighted_overdict[level][most_abundant]) - (
                    float(weighted_overdict[level][most_abundant]) * 0.1)
        sister_orgs = []
        for org in weighted_overdict[level]:
            #print(f"{compare_value}\t{weighted_overdict[level][org]}\t{org}")
            if float(weighted_overdict[level][org]) >= compare_value:
                sister_orgs.append(org)
        if len(sister_orgs)==1:
            x=sister_orgs[0]
            x_str="_".join(x.split("_")[1:])
            #return_lineage.insert(0, x_str)
            while x != "NA":
                x_str = "_".join(x.split("_")[1:])
                parent_x=parent_dict[x]
                parent_x_str="_".join(parent_x.split("_")[1:])
                return_lineage.insert(0,x_str)
                x=parent_x
            break
        else:
            continue
    original_len=len(return_lineage)
    while len(return_lineage)<len(req_taxa):
        x=len(req_taxa)-len(return_lineage)
        return_lineage.append(f"other_{return_lineage[original_len-1]}_{req_taxa[len(req_taxa)-x]}")
    return_lineage=";".join(return_lineage)
    #print(return_lineage)
    return return_lineage

#tp.tax_initiate(config.taxdumpdir)
#kraken_profilize()
#profile_to_files()
#log_otu_abundance_plot()
#print(list(lineage_scafs.keys())[1])
#run_minimap()
#read_minimap_paf()
post_gen_spec_classify()
#mag_lineage=json.load(open(mag_lin_jsonfile))
#lca_finder(mag_lineage["EUL_226"])
