import files as f

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
    shell:'echo standardise.py {input.cur}'

rule run_checkm:
    input:f.outfolder + "summary_all.tsv"
    output: f.outfolder + "/all_sample_summary.tsv"
    shell:'echo run_checkm.py'

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
    shell:"echo tax_profiler.py " + f.outfolder + tax_profile_out_DIR

rule cluster_bins:
    input:f.outfolder + "all_sample_summary.tsv"
    output:f.outfolder + "ani_results"+postfix+".tsv", f.outfolder + "clusters" + postfix +"_"+ method+ "/clusters.tsv"
    threads:5
    shell:'echo clusterbins.py {threads} {input} {output}'

rule post_cluster:
    input: f.outfolder + "clusters" + postfix +"_"+ method+ "/clusters.tsv"
    output: f.outfolder + "clusters" + postfix +"_"+ method+ "/clusters.tsv"
    shell: 'echo post_cluster.py'

rule cluster_abundance:
    input:f.raw_fastq, f.outfolder + "clusters" + postfix +"_"+ method+ "/clusters.tsv"
    output: f.outfolder+"cluster_abundance_profile.tsv", f.outfolder+"cluster_abundance_list.dump"
    shell:'echo cluster_abundance.py'

rule post_process_1:
    input: