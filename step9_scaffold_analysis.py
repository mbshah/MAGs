import json
import os
import subprocess

kraken_profile_folder = "../WMG/Kraken2_NT_Scaffolds_profiles2/"

mag_lin_jsonfile = kraken_profile_folder + "MAGlin.map"
scafLin_jsonfile = kraken_profile_folder + "scafLin.map"
linScaf_jsonfile = kraken_profile_folder + "linScaf.map"
scafLinMAG_jasonfile = kraken_profile_folder + "scafLinMAG.map"
lin_class_jsonfile = kraken_profile_folder + "linclass.map"
org_annotation_jsonfile = kraken_profile_folder + "tax_annot.map"
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
            num_of_scafs,sequence_len_scafs,gc = additional_stats(fasta)
            complete_org_annotation[org] = {}
            complete_org_annotation[org]["name"] = org
            complete_org_annotation[org]["tax"] = tax
            complete_org_annotation[org]["fasta"] = fasta
            complete_org_annotation[org]["gff"] = gff
            complete_org_annotation[org]["class"] = cl
            complete_org_annotation[org]["num_of_scafs"] = num_of_scafs
            complete_org_annotation[org]["total_len_of_scafs"]=sequence_len_scafs
            complete_org_annotation[org]["scaffolds_CDS"] = {}
            complete_org_annotation[org]["scaf_info"] = {}
            num_of_scafs_with_cds = 0
            cds_len_on_org = 0
            total_num_cds_on_org = 0
            if org in lin_mag:
                mags = lin_mag[org]
            complete_org_annotation[org]["mags"] = mags
            with open(gff, "r") as org_file:
                scaf_name = ""
                scaf_flag = 0
                cds_id = 1
                for line in org_file:
                    line = line.strip()
                    if line.startswith("//"):
                        complete_org_annotation[org]["scaf_info"][scaf_name] = {}
                        complete_org_annotation[org]["scaf_info"][scaf_name]["num_of_CDS"] = cds_id - 1
                        if complete_org_annotation[org]["scaf_info"][scaf_name]["num_of_CDS"] > 0:
                            num_of_scafs_with_cds = num_of_scafs_with_cds + 1
                            total_len_cds = 0
                            for cds in complete_org_annotation[org]["scaffolds_CDS"][scaf_name]:
                                total_len_cds = total_len_cds + int(
                                    complete_org_annotation[org]["scaffolds_CDS"][scaf_name][cds]["cds_len"])
                            avg_len = total_len_cds / complete_org_annotation[org]["scaf_info"][scaf_name]["num_of_CDS"]
                            complete_org_annotation[org]["scaf_info"][scaf_name]["avg_len_of_CDS"] = avg_len
                            complete_org_annotation[org]["scaf_info"][scaf_name]["total_len_of_CDS"] = total_len_cds
                            cds_len_on_org = cds_len_on_org + total_len_cds
                        scaf_flag = 0
                        scaf_name = ""
                        cds_id = 1
                    if line.startswith("DEFINITION"):
                        scaf_flag = 1
                        scaf_name = line.split()[1].split(";")[2].split("=")[1].strip("\"").strip()
                        complete_org_annotation[org]["scaffolds_CDS"][scaf_name] = {}
                    if scaf_flag == 1:
                        if line.startswith("CDS"):
                            cds_name = tax + "_" + scaf_name + "_cds" + str(cds_id)
                            cds_id = cds_id + 1
                            if line.__contains__("complement"):
                                line = line.replace("complement(", "").replace(")", "")
                                cds_start = int(line.split()[1].split("..")[1].replace(">", ""))
                                cds_end = int(line.split()[1].split("..")[0].replace("<", ""))
                                cds_len = cds_start - cds_end
                            else:
                                cds_start = int(line.split()[1].split("..")[0].replace("<", ""))
                                cds_end = int(line.split()[1].split("..")[1].replace(">", ""))
                                cds_len = cds_end - cds_start
                            total_num_cds_on_org = total_num_cds_on_org + 1
                            complete_org_annotation[org]["scaffolds_CDS"][scaf_name][cds_name] = {}
                            complete_org_annotation[org]["scaffolds_CDS"][scaf_name][cds_name]["cds_start"] = cds_start
                            complete_org_annotation[org]["scaffolds_CDS"][scaf_name][cds_name]["cds_end"] = cds_end
                            complete_org_annotation[org]["scaffolds_CDS"][scaf_name][cds_name]["cds_len"] = cds_len
            avg_cds_len_org = cds_len_on_org / total_num_cds_on_org if total_num_cds_on_org > 0 else 0
            hasmags = True if len(mags) > 0 else False
            genome_dir=kraken_profile_folder + f"/Bacteria/Genralists_Specialists/{cl}/fasta/{tax}/"
            growth_rate="NA"#run_growthrate_pred(tax,genome_dir)
            #run_kofamkoala(tax,genome_dir)
            tsv.write(
                f"{tax}\t{org}\t{num_of_scafs}\t{num_of_scafs_with_cds}\t{hasmags}\t{cl}\t{avg_cds_len_org}\t{total_num_cds_on_org}\t{cds_len_on_org}\t{sequence_len_scafs}\t{gc}\t{growth_rate}\n")
    org_annot_map = json.dumps(complete_org_annotation)
    f = open(org_annotation_jsonfile, "w")
    f.write(org_annot_map)
    f.close()
    tsv.close()


def additional_stats(fasta):
    num_of_scafs = int(os.popen(f"grep \">\" {fasta} |wc -l").read())
    len_of_scafs= int(os.popen(f"grep -v \">\" {fasta} |wc -m").read())
    actual_len=len_of_scafs-num_of_scafs
    gc=float("{:.2f}".format(float(os.popen("awk \'!/^>/{gc+=gsub(/[gGcC]/,\"\"); at+=gsub(/[aAtT]/,\"\");} END{ printf (gc*100)/(gc+at) }\' "+fasta).read())))
    return num_of_scafs, actual_len, gc


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


#run_prodigal()
read_prodigal()
