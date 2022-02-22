import config as f
import os

postfix="_dastool_90_10"
method="average"
tax_profile_out_DIR = "/profiles_"+postfix


rule all:
    input:
         f.outfolder + "summary_all.tsv",
         f.outfolder + "/all_sample_summary.tsv",
         f.outfolder + tax_profile_out_DIR + "/l1_superkingdom.tsv",
         f.outfolder + tax_profile_out_DIR + "/l2_phylum.tsv",
         f.outfolder + tax_profile_out_DIR + "/l3_class.tsv",
         f.outfolder + tax_profile_out_DIR + "/l4_order.tsv",
         f.outfolder + tax_profile_out_DIR + "/l5_family.tsv",
         f.outfolder + tax_profile_out_DIR + "/l6_genus.tsv",
         f.outfolder + tax_profile_out_DIR + "/l7_species.tsv",
         f.outfolder + "ani_results"+postfix+".tsv",
         f.outfolder + "clusters" + postfix +"_"+ method+ "/clusters.tsv"


rule standardise:
    input:
         cur=f.curation_file
    output:f.outfolder + "summary_all.tsv"
    shell:'echo step1_standardise.py {input.cur}'

rule run_checkm:
    input:f.outfolder + "summary_all.tsv"
    output: f.outfolder + "/all_sample_summary.tsv"
    shell:'echo step2_run_checkm.py'

rule tax_profiler:
    input:f.outfolder + "all_sample_summary.tsv"
    output:
          f.outfolder + tax_profile_out_DIR + "/l1_superkingdom.tsv",
          f.outfolder + tax_profile_out_DIR + "/l2_phylum.tsv",
          f.outfolder + tax_profile_out_DIR + "/l3_class.tsv",
          f.outfolder + tax_profile_out_DIR + "/l4_order.tsv",
          f.outfolder + tax_profile_out_DIR + "/l5_family.tsv",
          f.outfolder + tax_profile_out_DIR + "/l6_genus.tsv",
          f.outfolder + tax_profile_out_DIR + "/l7_species.tsv"
    shell:"echo step3_tax_profiler.py " + f.outfolder + tax_profile_out_DIR

rule cluster_bins:
    input:f.outfolder + "all_sample_summary.tsv"
    output:f.outfolder + "ani_results"+postfix+".tsv", f.outfolder + "clusters" + postfix +"_"+ method+ "/clusters.tsv"
    threads:5
    shell:'echo step4_clusterbins.py {threads} {input} {output}'

rule post_cluster:
    input: f.outfolder + "clusters" + postfix +"_"+ method+ "/clusters.tsv"
    output: f.outfolder + "clusters" + postfix +"_"+ method+ "/clusters.tsv"
    shell: 'echo step5_post_cluster.py'

rule cluster_abundance:
    input:f.raw_fastq, f.outfolder + "clusters" + postfix +"_"+ method+ "/clusters.tsv"
    output: f.outfolder+"cluster_abundance_profile.tsv", f.outfolder+"cluster_abundance_list.dump"
    shell:'echo cluster_abundance.py'

rule kraken_db_maker:
    output: f.anciliaryfolder+"/krakendb_nt/seqid2taxid.map",f.anciliaryfolder+"/krakendb_nt/taxonomy/names.dmp"
    shell:'kraken2_dbmaker.sh'

rule kraken_assign_tax:
    input:fasta=f.fasta_folder+"/"+"{{sample}}_min1000.fasta",database=f.anciliaryfolder+"/krakendb_nt/seqid2taxid.map" #fasta parameter can be expand()value as well, combined output has read ID which can be used to distinguish between samples in later stage
    output: f.outfolder+"/retax_combined.kout"
    params:
        db=f.anciliaryfolder+"/krakendb_nt",
        threads=f.processes
    shell: "../tools/kraken2/kraken2 --output {output} --threads {params.threads} --db {params.db} {input.fasta}"

