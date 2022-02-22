import collections
import json
import os
import subprocess
import urllib.request
from pathlib import Path
import re


seq_len_re=re.compile(r";seqlen=(\d+);")
gc_re=re.compile(r";gc_cont=(\d+.\d+);")
kraken_profile_folder = "../WMG/Kraken2_NT_Scaffolds_profiles2/"
#kraken_profile_folder = "/home/manan/disk2/UDE/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/"

mag_lin_jsonfile = kraken_profile_folder + "MAGlin.map"
scafLin_jsonfile = kraken_profile_folder + "scafLin.map"
linScaf_jsonfile = kraken_profile_folder + "linScaf.map"
scafLinMAG_jasonfile = kraken_profile_folder + "scafLinMAG.map"
lin_class_jsonfile = kraken_profile_folder + "linclass.map"
org_annotation_jsonfile = kraken_profile_folder + "tax_annot.map"
kegg_db_folder=kraken_profile_folder+"../../ancilary/kegg/"
org_tsv_out=kraken_profile_folder + "tax_annot.tsv"
mag_lineage = json.load(open(mag_lin_jsonfile))
scaf_lin = json.load(open(scafLin_jsonfile))
lin_scaf = json.load(open(linScaf_jsonfile))
lin_class = json.load(open(lin_class_jsonfile))
lin_mag = {}
complete_org_annotation = {}
for mag in mag_lineage:
    lineage = mag_lineage[mag]
    if lineage in lin_mag:
        pass
    else:
        lin_mag[lineage] = []
    lin_mag[lineage].append(mag)

#for lineage in lin_mag:
#    print(f"{lineage}\t{len(lin_mag[lineage])}\t{lin_mag[lineage]}")
#print(f"{len(mag_lineage.keys())}")
kegg_db_file=kegg_db_folder+ "kegg.db"
if Path(kegg_db_file).is_file():
    kegg_database=json.load(open(kegg_db_file))
else:
    kegg_database = {"Genes": {}, "Pathways": {}}


def run_prodigal():
    list_of_orgs = lin_scaf.keys()
    for org in list_of_orgs:
        if org in lin_class:
            classification = str(lin_class[org])
            cl = "NS" if classification.startswith("NON SIGNIFICANT") else "GEN" if classification.startswith(
                "GENERALIST") else "SPEC" if classification.startswith("SPECIALIST") else "XXXXXXXXXXXXXXXXXX"
            tax = scaf_lin[lin_scaf[org].split(";")[0]][1]
            fasta = kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/org.fasta"
            gff = kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/porg.gff"
            faa = kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/porg.faa"
            fna = kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/porg.fna"
            # os.system(f"rm {gff}")
            cmd = f"prodigal -i {fasta} -o {gff} -q -p meta -a {faa} -d {fna}"
            # cmd = f"../tools/gms2_linux_64/gms2.pl --seq {fasta} --output {gff} --genome-type bacteria"
            print(cmd)
            #os.system(cmd)
            ##cat bash_commands.sh| parallel -j 15 -k -u {}


def read_prodigal():
    tsv=open(org_tsv_out,"w")
    tsv.write(
        f"taxonomy_ID\tInferred_Lineage\tNumber_fo_Scaffolds\tScaffolds_with_CDS\tHasMAGS\tClassification\tAvg_Len_of_CDS\tTotal_number_of_cds\ttotal_len_ofCDS\ttotal_Len_Scaffolds\tGC\tGrowth_rate\n")
    list_of_orgs = lin_scaf.keys()
    for org in list_of_orgs:
        if org in lin_class:
            classification = str(lin_class[org])
            cl = "NS" if classification.startswith("NON SIGNIFICANT") else "GEN" if classification.startswith(
                "GENERALIST") else "SPEC" if classification.startswith("SPECIALIST") else "XXXXXXXXXXXXXXXXXX"
            mags = []
            tax = scaf_lin[lin_scaf[org].split(";")[0]][1]
            fasta = kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/org.fasta"
            gff = kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/porg.gff"
            num_of_scafs=0
            sequence_len_scafs=0
            cmd = "awk \'!/^>/{gc+=gsub(/[gGcC]/,\"\"); at+=gsub(/[aAtT]/,\"\");} END{ printf (gc*100)/(gc+at) }\' " + fasta
            # print(cmd)
            gc = float("{:.2f}".format(float(os.popen(cmd).read().replace(",", "."))))
            complete_org_annotation[org] = {}
            complete_org_annotation[org]["name"] = org
            complete_org_annotation[org]["tax"] = tax
            complete_org_annotation[org]["fasta"] = fasta
            complete_org_annotation[org]["gff"] = gff
            complete_org_annotation[org]["class"] = cl
            complete_org_annotation[org]["num_of_scafs"] = num_of_scafs
            complete_org_annotation[org]["total_len_of_scafs"]=sequence_len_scafs
            if org in lin_mag:
                mags = lin_mag[org]
            complete_org_annotation[org]["mags"] = mags
            sequence_len_scafs,cds_len_on_org,num_of_scafs,total_num_cds_on_org,num_of_scafs_with_cds=gff_processor(gff)
            avg_cds_len_org = cds_len_on_org / total_num_cds_on_org if total_num_cds_on_org > 0 else 0
            hasmags = True if len(mags) > 0 else False
            genome_dir=kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/"
            growth_rate="NA"#run_growthrate_pred(tax,genome_dir)
            kofam_raw_output = Path(genome_dir + tax + "_ko_raw.tsv")
            #if kofam_raw_output.is_file(): print("here")
            #run_kofamkoala(tax,genome_dir)
            tsv.write(
                f"{tax}\t{org}\t{num_of_scafs}\t{num_of_scafs_with_cds}\t{hasmags}\t{cl}\t{avg_cds_len_org}\t{total_num_cds_on_org}\t{cds_len_on_org}\t{sequence_len_scafs}\t{gc}\t{growth_rate}\n")
    org_annot_map = json.dumps(complete_org_annotation)
    f = open(org_annotation_jsonfile, "w")
    f.write(org_annot_map)
    f.close()
    tsv.close()


kofam_folder="../tools/KOfamKoala/kofam_scan-1.3.0/"


def run_kofamkoala(tax,genome_dir):
    gene_pred_faa = genome_dir + "/porg.faa"
    config_file = kofam_folder + "config.yml"
    kofam_raw_output = genome_dir + tax + "_ko_raw.tsv"
    kofam_filtered_output = genome_dir + tax + "_ko_filtered.tsv"
    kofam_exe = f"{kofam_folder}exec_annotation"
    command = f"{kofam_exe} -c {config_file} -f detail-tsv -o {kofam_raw_output} {gene_pred_faa}"
    print(command)
    os.system(command)
    command = f"head -n1 {kofam_raw_output} > {kofam_filtered_output}"
    print(command)
    os.system(command)
    command = f'grep "^*" {kofam_raw_output} >>{kofam_filtered_output}'
    print(command)
    os.system(command)


growthpred_folder = "../tools/growthpred-v1.07/"


def run_growthrate_pred(tax,genome_dir):
    gene_pred_fnn = genome_dir + "/porg.fna"
    code = 0
    gp_out = genome_dir + "gp_out"
    gp_out_final = genome_dir + tax + "_growthpred.result"
    #ogt = get_ogt(genomes_metadata[genome]["Members"])
    command = f"python2 {growthpred_folder}/growthpred-v1.07.py -b -g {gene_pred_fnn} -o {gp_out} -c {code} -t -m -s -S"
    print(command)
    os.system(command)
    command = f"mv gp_out.results {gp_out_final}"
    os.system(command)
    # command=f'grep "Predicted minimum generation time" {infolder}/*/*'
    out = (subprocess.Popen(["grep", "Predicted minimum generation time", f'{gp_out_final}'],
                            stdout=subprocess.PIPE).communicate()[0]).decode("utf-8")
    retval=out if out.__contains__("Predicted minimum generation time") else "0"
    retval = retval.replace("\n", "").replace("Predicted minimum generation time:  ", "")
    return retval


def q_kofam():
    list_of_orgs = lin_scaf.keys()
    for org in list_of_orgs:
        if org in lin_class:
            classification = str(lin_class[org])
            cl = "NS" if classification.startswith("NON SIGNIFICANT") else "GEN" if classification.startswith(
                "GENERALIST") else "SPEC" if classification.startswith("SPECIALIST") else "XXXXXXXXXXXXXXXXXX"
            tax = scaf_lin[lin_scaf[org].split(";")[0]][1]
            fasta = kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/org.fasta"
            gff = kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/porg.gff"
            genome_dir = kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/"
            kofam_raw_output = Path(genome_dir + tax + "_ko_raw.tsv")
            if kofam_raw_output.is_file():pass
            #    run_kofamkoala(tax, genome_dir)
            else:
                run_kofamkoala(tax, genome_dir)


def summarize_kofam_results():
    table= {}
    header=[]
    uniq_classes=[]
    uniq_pathways=[]
    all_KIDs=[]
    with open(org_tsv_out,"r")as org_tsv:
        header=next(org_tsv).strip().split("\t")
        for line in org_tsv:
            entry=line.strip().split("\t")
            table[entry[0]] =table[entry[0]]={}
            for i in range(0,len(header)):
                table[entry[0]][header[i]]=entry[i]
            kofam_out=(Path(kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{entry[5]}/fasta/{entry[0]}/{entry[0]}_ko_filtered.tsv"))
            number_of_KOs=int(os.popen(f"grep \"^*\" {kofam_out} |wc -l").read())
            org_kids=[]
            org_pathways=[]
            org_classes= []
            #print(table[entry[0]]["taxonomy_ID"])
            if os.path.getsize(kofam_out)==0:pass
            else:
                with open(kofam_out) as kofam_tsv:
                    next(kofam_tsv)
                    for kLine in kofam_tsv:
                        kEntry=kLine.split("\t")
                        kID=kEntry[2]
                        org_kids.append(kID)
            kIDs_c=dict(collections.Counter(org_kids))
            #pathways = getpathways(kIDs_c)
            #org_pathways=org_pathways+pathways
            #pathway_count=dict(collections.Counter(org_pathways))
            #org_classes=org_classes+getkclasses(pathway_count)
            #class_count=dict(collections.Counter(org_classes))
            #print(class_count)
            pathway_count,class_count=kegg_info_accumulate(list(kIDs_c.keys()))
            for keggID in kIDs_c:
                if keggID not in all_KIDs: all_KIDs.append(keggID)
                table[entry[0]][keggID]=1
            for kpathway in pathway_count:
                if kpathway not in uniq_pathways:uniq_pathways.append(kpathway)
                table[entry[0]][kpathway]=pathway_count[kpathway]
            for kclass in class_count:
                if kclass not in uniq_classes:uniq_classes.append(kclass)
                table[entry[0]][kclass]=class_count[kclass]
            number_of_KOs=len(kIDs_c)
            table[entry[0]]["no_of_KOs"]=number_of_KOs
    t_o=org_tsv_out.replace(".tsv","")+"_ko_summary.tsv"
    if "no_of_KOs" not in header: header.append("no_of_KOs")
    for kclass in uniq_classes:
        if kclass not in header:header.append(kclass)
    for kpathway in uniq_pathways:
        if kpathway not in header:header.append(kpathway)
    for keggID in all_KIDs:
        if keggID not in header:header.append(keggID)
    x=open(t_o,"w")
    x.write("\t".join(header))
    x.write("\n")
    for entry in table:
        line=[]
        for value_h in header:
            if value_h in table[entry]:value=table[entry][value_h]
            else:value=0
            line.append(str(value))
        x.write("\t".join(line)+"\n")
    x.close()
    kegg_map = json.dumps(kegg_database)
    f = open(kegg_db_file, "w")
    f.write(kegg_map)
    f.close()


def kegg_info_accumulate(genes_list):
    pathway_count={}
    class_count={}
    for kid in genes_list:
        kid_pathways=getpathways2(kid)
        for pathway in kid_pathways:
            if pathway in pathway_count: pathway_count[pathway]=pathway_count[pathway]+1
            else: pathway_count[pathway]=1
            pway_class=getkclasses2(pathway)
            for kclass in pway_class:
                if kclass in class_count: class_count[kclass]=class_count[kclass]+1
                else:class_count[kclass]=1
    return pathway_count, class_count


def getpathways2(id):
    kFile = kegg_db_folder + id
    if Path(kFile).is_file():
        pass
    else:
        url = "http://rest.kegg.jp/get/" + id
        # print(f"retreiving {url} as {kFile}")
        urllib.request.urlretrieve(url, kFile)
    pathwayflag=0
    pathways=[]
    with open(kFile, "r") as file:
        kegg_database["Genes"][id]={}
        kegg_database["Genes"][id]["pathways"]=[]
        for line in file:
            if line[0].isupper():
                pathwayflag=0
            if pathwayflag==1:
                pathway=line.split()[0].strip().replace(" ","_").replace(",","_")
                pathways.append(pathway)
                kegg_database["Genes"][id]["pathways"].append(pathway)
            if line.startswith("PATHWAY"):
                pathwayflag=1
                pathway=line.split()[1]
                pathways.append(pathway)
                kegg_database["Genes"][id]["pathways"].append(pathway)
            if line.startswith("NAME"):
                kegg_database["Genes"][id]["name"]=" ".join(line.split(" ")[1:]).strip()
            if line.startswith("DEFINITION"):
                kegg_database["Genes"][id]["def"]=" ".join(line.split(" ")[1:]).strip()
    return pathways


def getkclasses2(koID):
    koClass = []
    koFile = kegg_db_folder + koID
    if Path(koFile).is_file():
        pass
    else:
        url = "http://rest.kegg.jp/get/" + koID
        # print(f"retreiving {url} as {kFile}")
        urllib.request.urlretrieve(url, koFile)
    with open(koFile, "r") as file:
        kegg_database["Pathways"][koID]={}
        kegg_database["Pathways"][koID]["class"]=""
        orthoflag=0
        kegg_database["Pathways"][koID]["orthos"]=[]
        for line in file:
            if orthoflag==1:
                if line.startswith(r"[A-Z]"): orthoflag=0
                else:
                    ortho=line.split()[0].strip().replace(" ","_").replace(",","_")
                    kegg_database["Pathways"][koID]["orthos"].append(ortho)
            if line.startswith("ORTHOLOGY"):
                orthoflag=1
                ortho=line.split()[1]
                kegg_database["Pathways"][koID]["orthos"].append(ortho)
            if line.startswith("NAME"):
                kegg_database["Pathways"][koID]["name"]=" ".join(line.split(" ")[1:]).strip()
            if line.startswith("CLASS"):
                kegg_database["Pathways"][koID]["class"]=list(filter(None, line.strip().split('  ')))[1]
                kClass = list(filter(None, line.split('  ')))[1].split(";")
                for i in kClass:
                    i1="kc_"+i.strip().replace(" ","_").replace(",","_")
                    koClass.append(i1)
    return koClass


def gff_processor(gff_file):
    total_len=0
    coded_len=0
    scaf_flag=0
    scaf_arr=[]
    scaf_count=0
    cds_count=0
    p_cds=0
    scaf_with_cds=0
    with open(gff_file, "r")as gff:
        for line in gff:
            line=line.strip()
            if scaf_flag==1:
                if line.startswith("//"):
                    total_len=total_len+int(scaf_len)
                    coded_len=coded_len+len(scaf_arr)
                    if cds_count>p_cds:scaf_with_cds=scaf_with_cds+1
                    #print(f"{total_len}\t{coded_len}")
                if line.startswith("CDS"):
                    cds_count=cds_count+1
                    coords=line.split()[1].replace("complement(","").replace(")","").replace("<","").replace(">","").split("..")
                    start=coords[0]
                    end=coords[1]
                    if len(scaf_arr)>1:
                        scaf_arr=scaf_arr+list(range(int(start),(int(end))))
                    else:
                        scaf_arr=list(range(int(start),(int(end))))
            if line.startswith("DEFINITION"):
                scaf_flag=1
                scaf_len=seq_len_re.search(line).group(1)
                scaf_gc=float(gc_re.search(line).group(1))
                #print(scaf_len)
                scaf_arr=[]
                scaf_count=scaf_count+1
                p_cds=cds_count
    return total_len, coded_len, scaf_count, cds_count,scaf_with_cds

def svg_edit():
    svg_file="/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/cca_gen_spec.svg"
    svg_out="/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/cca_gen_spec2.svg"
    with open(svg_file,"r") as file:
        for line in file:
            print(line)


#run_prodigal()
#q_kofam()
#read_prodigal()
summarize_kofam_results()
print(kegg_database["Pathways"]["ko05212"]["orthos"])
#readkfiles("/home/manan/disk2/UDE/mbs_general/ancilary/kegg/K00073")
#gff_processor("/home/manan/disk2/UDE/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/GEN/fasta/22/porg.gff")
