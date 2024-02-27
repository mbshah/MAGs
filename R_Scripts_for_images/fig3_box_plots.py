import math
import pandas as pd
import config as config #project config script
import seaborn as sns
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

infolder = config.outfolder + "/dRep__gff/"
clusters_file = infolder + "post_cluster_data_Rver.csv"
KO_file = infolder + "ko_profile.tsv"
pathways_file = infolder + "clusters_pathway_summary.tsv"
subSystem_file = infolder + "clusters_subSystem_summary.tsv"
codon_table_file = infolder + "post_cluster_data_codon_table.csv"
members_file = infolder + "members_table_super_populated.csv"
master_table = pd.read_csv(clusters_file, index_col="cluster_ID")
KOtable = pd.read_csv(KO_file, sep="\t", index_col="Data")
codon_table = pd.read_csv(codon_table_file, index_col="Cluster")
members_table = pd.read_csv(members_file, index_col="bin")
streamlined_table = master_table[(master_table["reassembly_size_mb"] <= 2)]
non_streamlined_table = master_table[(master_table["reassembly_size_mb"] >= 4)]
intermediate_streamlined_table = master_table[
    (master_table["reassembly_size_mb"] < 4) & (master_table["reassembly_size_mb"] > 2)]
members_abundance_profile = pd.read_csv(infolder + "cluster_abundance_profile_normalized_deseq2.tsv", sep="\t",
                                        index_col="abundance").transpose()
read_counts = pd.read_csv(config.read_counts_file, sep="\t", index_col="Sample")


def load_from_table():
    global master_table
    global members_table
    global non_streamlined_table
    global streamlined_table
    global intermediate_streamlined_table

    master_table["GC_lev"] = master_table.reassembly_gc.apply(
        lambda x: "low" if x < 45 else "high" if x > 55 else "Balanced")
    master_table["Temp_lev"] = master_table.site_Temperature.apply(lambda x: "low" if x < 20 else "high")
    master_table["pH_lev"] = master_table.site_pH.apply(
        lambda x: "Acidic" if x < 6.8 else "Basic" if x > 7.2 else "Neutral")
    master_table["TP_lev"] = master_table.weighted_TP.apply(lambda x: "low" if x < 10 else "high")
    master_table["Coding_lev"] = master_table.coding_perc.apply(lambda x: "Low" if x < 94 else "High")
    master_table["Genome_Size"] = master_table.reassembly_size_mb.apply(
        lambda x: "Small" if x <= 2 else "Large" if x >= 3.5 else "Medium")
    master_table["TP_log"] = master_table.weighted_TP.apply(lambda x: math.log(x))
    streamlined_table = master_table[(master_table["reassembly_size_mb"] <= 2)]
    non_streamlined_table = master_table[(master_table["reassembly_size_mb"] >= 4)]
    intermediate_streamlined_table = master_table[
        (master_table["reassembly_size_mb"] < 4) & (master_table["reassembly_size_mb"] > 2)]
    if "GC_content" not in members_table.columns:members_table["GC_content"]=members_table.cluster.apply(
        lambda x: master_table.loc[x,"reassembly_gc"])
    if "log_DP" not in members_table.columns: members_table["log_DP"] = members_table["DP.y"].apply(
        lambda x: math.log(x))
    if "log_TP" not in members_table.columns: members_table["log_TP"] = members_table["TP.y"].apply(
        lambda x: math.log(x))
    if "log_DN" not in members_table.columns: members_table["log_DN"] = members_table["DN.y"].apply(
        lambda x: math.log(x))
    if "abundance_val" not in members_table.columns: members_table["abundance_val"]=members_table.apply(
        lambda row:members_abundance_profile.loc[row["site_list"],row["cluster"]],axis=1)
    if "TNs" not in members_table.columns: members_table["TNs"]=members_table.apply(
        lambda row:float(row["DN.y"])+float(row["NH4-N[µg/l]"])+float(row["NO3-N[µg/l]"]),axis=1)
    members_table.to_csv(members_file)

def figure_box_plots():
    fig,axs=plt.subplots(2,3,figsize=(10,6))

    x = "Genome_Size"
    y = "reassembly_gc"
    plot_box = sns.boxplot(data=master_table, y=y, x=x, ax=axs[0,0],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Large"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Large"), ("Small", "Medium"), ("Large", "Medium")],
                        order=["Small", "Medium", "Large"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome size (Mb)', ylabel='GC content (%)')
    axs[0,0].set_title('(A)')

    x = "Genome_Size"
    y = "Total_Redundancy"
    plot_box = sns.boxplot(data=master_table, y=y, x=x,ax=axs[0,1],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Large"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Large"), ("Small", "Medium"), ("Large", "Medium")],
                        order=["Small", "Medium", "Large"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome size (Mb)', ylabel='No. of redundant\ngenes in a genome')
    axs[0,1].set_title('(B)')

    x = "Genome_Size"
    y = "coding_perc"
    plot_box = sns.boxplot(data=master_table, y=y, x=x,ax=axs[0,2],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Large"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Large"), ("Small", "Medium"), ("Large", "Medium")],
                        order=["Small", "Medium", "Large"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome size (Mb)', ylabel='Coding regions (%)')
    axs[0,2].set_title('(C)')

    x = "Genome_Size"
    y = "pc_abc_transporters"
    plot_box = sns.boxplot(data=master_table, y=y, x=x, ax=axs[1,0],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Large"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Large"), ("Small", "Medium"), ("Large", "Medium")],
                        order=["Small", "Medium", "Large"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome size (Mb)', ylabel='Fraction of genes for\nABC transporters (%)')
    axs[1,0].set_title('(D)')

    x = "Genome_Size"
    y = "pc_sigmafactors"
    plot_box = sns.boxplot(data=master_table, y=y, x=x, ax=axs[1,1],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Large"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Large"), ("Small", "Medium"), ("Large", "Medium")],
                        order=["Small", "Medium", "Large"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome size (Mb)', ylabel='Fraction of sigma\nfactor genes (%)')
    axs[1,1].set_title('(E)')

    x = "Genome_Size"
    y = "pc_cell_motility"
    plot_box = sns.boxplot(data=master_table, y=y, x=x, ax=axs[1,2],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Large"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Large"), ("Small", "Medium"), ("Large", "Medium")],
                        order=["Small", "Medium", "Large"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome size (Mb)', ylabel='Fraction of genes for\ncell motility (%)')
    axs[1,2].set_title('(F)')
    plt.tight_layout()
    plt.show()

load_from_table()
figure_box_plots()
