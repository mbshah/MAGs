import math
import os
import subprocess
import pandas as pd
import random
import config as config
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats as sci
import statistics as st
from statannot import add_stat_annotation
import re

###Some import statements are for commented lines LEAVE THEM IN


infolder = config.outfolder + "clusters" + config.postfix + "_" + config.method + "/"
clusters_file = infolder + "clusters_wphyp_sited.tsv"
KO_file = infolder + "clusters_ko_summary.tsv"
pathways_file = infolder + "clusters_pathway_summary.tsv"
subSystem_file = infolder + "clusters_subSystem_summary.tsv"
codon_table_file = infolder + "Cluster_Codon_Table.tsv"
master_table = pd.read_csv(clusters_file, sep="\t", index_col="Cluster")
KOtable = pd.read_csv(KO_file, sep="\t", index_col="Cluster")
pathwaysTable = pd.read_csv(pathways_file, sep="\t", index_col="Cluster")
subSystemTable = pd.read_csv(subSystem_file, sep="\t", index_col="Cluster")
codon_table = pd.read_csv(codon_table_file, sep="\t", index_col="Cluster")
kegg_db_folder = "../ancilary/kegg/"
members_table = pd.DataFrame(
    columns=['Member', "site_pH", "site_Temperature", "site_Elevation", "site_TP", "wt_TP", "log_TP", "Phylum", "Genus",
             "org_size", "abundance", "MAG"])
streamlined_table = master_table[(master_table["reassembly_size_mb"] <= 2)]
non_streamlined_table = master_table[(master_table["reassembly_size_mb"] >= 4)]
intermediate_streamlined_table = master_table[
    (master_table["reassembly_size_mb"] < 4) & (master_table["reassembly_size_mb"] > 2)]
#members_abundance_profile = pd.read_csv(config.outfolder + "cluster_abundance_profile_transformed.tsv", sep="\t",
#                                         index_col="abundance")
members_abundance_profile = pd.read_csv(config.outfolder + "cluster_abundance_profile_normalized_deseq2.tsv", sep="\t",
                                        index_col="Cluster")
read_counts = pd.read_csv(config.read_counts_file, sep="\t", index_col="Sample")


# print(codon_table.index)

def make_ko_category_table():
    small_genomes = master_table[master_table["Genome_Size"] == "Small"].index.values
    kosmall = KOtable[KOtable.index.isin(small_genomes)].drop(columns="no_of_KOs")
    big_genomes = master_table[master_table["Genome_Size"] == "Big"].index.values
    kobig = KOtable[KOtable.index.isin(big_genomes)].drop(columns="no_of_KOs")
    intermediate_genomes = master_table[master_table["Genome_Size"] == "Medium"].index.values
    kointermediate = KOtable[KOtable.index.isin(intermediate_genomes)].drop(columns="no_of_KOs")
    gene_data_table_file = infolder + "gene_differential_presence.tsv"
    out = open(gene_data_table_file, "w")
    out.write("KO\tSmall\tBig\tIntermediate\n")
    for ko in KOtable.drop(columns="no_of_KOs"):
        small_perc = sum(kosmall[ko]) / len(small_genomes)
        big_perc = sum(kobig[ko]) / len(big_genomes)
        intermediate_perc = sum(kointermediate[ko]) / len(intermediate_genomes)
        out.write(f"{ko}\t{small_perc}\t{big_perc}\t{intermediate_perc}\n")
    out.close()


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
    master_table["Coding_lev"] = master_table.perc_coding_combined.apply(lambda x: "Low" if x < 94 else "High")
    master_table["Genome_Size"] = master_table.reassembly_size_mb.apply(
        lambda x: "Small" if x <= 2 else "Big" if x >= 3.5 else "Medium")
    master_table["TP_log"] = master_table.weighted_TP.apply(lambda x: math.log(x))
    master_table["cluster_abundance"] = master_table.apply(
        lambda x: round(sum(transformed_abundace_profile.loc[x.name,]), 3), axis=1)
    streamlined_table = master_table[(master_table["reassembly_size_mb"] <= 2)]
    non_streamlined_table = master_table[(master_table["reassembly_size_mb"] >= 4)]
    intermediate_streamlined_table = master_table[
        (master_table["reassembly_size_mb"] < 4) & (master_table["reassembly_size_mb"] > 2)]
    master_table["CDS"] = master_table.CDS.apply(lambda x: int(x.replace(",", "")))
    # print(f"{max(master_table['cluster_abundance'])}\t{min(master_table['cluster_abundance'])}")
    master_table.to_csv(infolder + "clusters_wphyp_sited_m.tsv", sep="\t")
    x = 0
    # ot=open(f"{infolder}merged_fna_all_clusters.fna","w") ##to run CUB seperately

    for line in master_table.iterrows():
        data = line[1]
        wt_TP = data["weighted_TP"]
        genus = data["Kraken_Genus"]
        phylum = data["Kraken_Phyla"]
        ph = data["site_pH_Array"].split("; ")
        temp = data["site_Temperature_Array"].split("; ")
        elev = data["site_Elevation_Array"].split("; ")
        tp = data["site_TP_Array"].split("; ")
        members = data["Members"].replace('"', '').split(";")
        no_of_members = data["NumberofMembers"]
        size = data["reassembly_size_mb"]
        log_TP = data["TP_log"]
        for i in range(0, no_of_members):
            x = x + 1
            abund = members_abundance_profile.loc[data.name, "_".join(members[i].split("_")[0:2])]
            genus2 = re.sub("(\d+)", "", genus).replace("()", "")
            row = (
                [members[i], float(ph[i]), float(temp[i]), float(elev[i]), (float(tp[i])), wt_TP, log_TP,
                 phylum, genus2,
                 size, abund, data.name])
            members_table.loc[x] = row
        # fna_file=infolder+data.name+"/"+data.name+"_gms.fna"
        # command=f"grep -v '^>' {fna_file}"
        # sequence_combined=subprocess.run([command], stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8").strip().replace("\n","")
        # sequence_combined_entry=f">{data.name}\n{sequence_combined}\n"
        # ot.write(sequence_combined_entry)
    # ot.close()


def abundance_transform(abundance_profile, method="hellinger", optional=None):
    abundance_profile = pd.DataFrame(abundance_profile)
    if optional is None:
        optional = {}
    if method == "hellinger":
        read_counts = pd.DataFrame(optional)
        samples = abundance_profile.columns
        orgs = abundance_profile.index
        for sample in samples:
            total_reads = abundance_profile.sum(axis=0)[sample]
            for org in orgs:
                abundance_value = abundance_profile.loc[org][sample]
                transformed_abundance = math.sqrt(abundance_value)
                abundance_profile.loc[org][sample] = transformed_abundance
                #print(f"{sample}\t{org}\t{abundance_value}\t{total_reads}")
    # print(abundance_profile)
    return abundance_profile


def correlator(data, fact1, fact2,method):
    if method=="anova":
        corr=sci.f_oneway(data[data[fact1]=='Small'][fact2].values,data[data[fact1]=='Medium'][fact2].values,data[data[fact1]=='Big'][fact2].values)
        print(f'anova of {fact1} vs {fact2} in the table is {corr}')
        return corr.pvalue
    elif method == "pearson":
        corr = sci.pearsonr(data[fact1], data[fact2])
        print(f'Pearsons correlation of {fact1} against {fact2}: {corr}')
        return corr.pvalue


def getKOsforPathway(pathway):
    pathwayfile = kegg_db_folder + "/" + pathway
    command = f"grep -l {pathway} {kegg_db_folder}/* |grep \"K\""
    # print(command)
    out = subprocess.run([command], stdout=subprocess.PIPE, shell=True)
    ko_list = out.stdout.decode("utf-8").strip().replace(kegg_db_folder + "/", "").split()
    KO_stats = {}
    for ko in ko_list:
        S_pres = 0
        NS_pres = 0
        IS_pres = 0
        for org in KOtable.index:
            ispresent = True if ko in KOtable and KOtable.loc[org, ko] > 0 else False
            if ispresent == True:
                if org in streamlined_table.index.values:
                    S_pres = 1
                if org in non_streamlined_table.index.values:
                    NS_pres = 2
                if org in intermediate_streamlined_table.index.values:
                    IS_pres = 4
        status = S_pres + NS_pres + IS_pres
        KO_stats[ko] = status
    return KO_stats


def ko_distinguisher(ko):
    S_pres = 0
    NS_pres = 0
    IS_pres = 0
    for org in KOtable.index:
        ispresent = True if ko in KOtable and KOtable.loc[org, ko] > 0 else False
        if ispresent == True:
            if org in streamlined_table.index.values:
                S_pres = 1
            if org in non_streamlined_table.index.values:
                NS_pres = 2
            if org in intermediate_streamlined_table.index.values:
                IS_pres = 4
    status = S_pres + NS_pres + IS_pres
    return status


def rgb_to_hex(rgb):
    rgb = tuple([int(255 * val) for val in rgb])
    return '#' + ''.join([hex(val)[2:] for val in rgb]).upper()


def pathway_plots():  ####NOT WORKING USING PATHVIEW IN R INSTEAD
    color_pellet = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499",
                    "#882255"]  # Tol colorblind_pallate from Paul Tol https://davidmathlogic.com/colorblind/
    from IPython.core.display import Image
    from Bio import SeqIO
    import Bio.KEGG.REST as kg
    from Bio.KEGG.KGML import KGML_parser
    from Bio.Graphics.KGML_vis import KGMLCanvas
    from Bio.Graphics.ColorSpiral import ColorSpiral
    pathway_of_int = "ko00790"
    # KO_hash=getKOsforPathway(pathway_of_int)
    pathway = KGML_parser.read(kg.kegg_get(pathway_of_int, "kgml"))
    for element in pathway.orthologs:
        names = element.name.replace("ko:", "").split()
        for i in range(len(element.graphics)):
            print(element.graphics[i])
            # ortho = names[i]
            # color_code = ko_distinguisher(ortho)
            # color = color_pellet[color_code]
            # x=element.graphics[i]
            # print(f"{ortho}\t{color_code}")
        # print(element)
    canvas = KGMLCanvas(pathway)
    canvas.show_genes = True
    canvas.import_imagemap = True
    canvas.draw("figure.pdf")


def print_full(x):
    pd.set_option('display.max_rows', len(x))
    # pd.set_option('display.max_columns',len(x.columns))
    print(x)
    pd.reset_option('display.max_rows')
    # pd.reset_option('display.max_columns')


def data_plots2():
    color_pellet = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499",
                    "#882255"]  # Tol colorblind_pallate from Paul Tol https://davidmathlogic.com/colorblind/

    # FIG1
    # org="Actinobacteria"
    # org_list=list(master_table.Kraken_Phyla.unique())
    # for org in org_list:
    #    if str(org).__contains__("Ol"):continue
    #    reduced_table=master_table[master_table["Kraken_Phyla"]==org]
    #    #print(reduced_table)
    #    #sns.kdeplot(reduced_table["reassembly_size_mb"],label=org+str(reduced_table.shape))
    # plt.legend()
    # plt.show()

    # Fig 4.1.1 (scatter_plot) phyla_sperated log TP vs reassembly Size
    # print(master_table.columns)
    # x = sns.FacetGrid(data=master_table,
    #                   col="phylophlan_phyla", col_wrap=4, height=4, aspect=.7,
    #                   #hue=master_table.NumberofMembers,
    #                   #palette=color_pellet
    #                   )
    # x.map(sns.scatterplot, "reassembly_size_mb", "TP_log")
    # plt.show()

    # Fig 4.1.2 (scatter_plot)
    # sns.scatterplot(x=master_table.TP_log, y=master_table.reassembly_size_mb,
    #                 size=master_table.cluster_abundance, hue=master_table.NumberofMembers,
    #                 palette=color_pellet, sizes=(20,200)
    #                 )
    # sns.set_palette("dark")
    # plt.show()

    # Fig 4.1.3 (scatter_plot with members_table) ##NEWLY MODIFIED FOR NEW ABUNDANCE
    ax = sns.scatterplot(x=members_table.log_TP, #.apply(lambda x: "low" if x < 10 else "high"),
                         y=members_table.org_size,  # .apply(lambda x: "Small" if x <= 2 else "Big" if x >= 3.5 else "Medium"),
                         size=members_table.abundance,
                         # hue=master_table.NumberofMembers,
                         palette=color_pellet,
                         sizes=(2, 200),
                         # order=["Small","Medium","Big"]
                         )
    sns.set_palette("dark")
    # ax.get_legend().remove()
    plt.legend(title="square root\ntransformed\nabundance")
    plt.show()

    # fig 4.2 (stripplot)
    # sns.stripplot(x=master_table.TP_lev, y=master_table.reassembly_size_mb,
    #               # col=master_table.phylophlan_phyla, col_wrap=4, height=4, aspect=.7,
    #               size=master_table.cluster_abundance,
    #               hue=master_table.NumberofMembers,
    #               palette=color_pellet,
    #               sizes=(10, 100)
    #               )
    # sns.set_palette("dark")
    # plt.show()

    # Fig 5.1 Scatter perc_coding vs genome size with marginal density
    # scatter_plot = sns.JointGrid(y=master_table.perc_coding_combined, x=master_table.reassembly_size_mb,
    #                              # col=master_table.phylophlan_phyla, col_wrap=4, height=4, aspect=.7,
    #                              # size=master_table.cluster_abundance, hue=master_table.NumberofMembers,
    #                              palette=color_pellet,
    #                              xlim=(0, 8)
    #                              )
    # scatter_plot.plot_joint(sns.scatterplot)
    # scatter_plot.ax_joint.axvline(2, ls="--")
    # scatter_plot.plot_marginals(sns.kdeplot)
    # scatter_plot.ax_joint.axvline(4, ls="--")
    # plt.show()

    # Fig 5.1 scatter plot according to phyla, popular phyla
    # req_phyla=["Proteobacteria","Actinobacteria","Bacteroidetes","Cyanobacteria"]
    # x = sns.FacetGrid(data=master_table[master_table['phylophlan_phyla'].apply(lambda x: any([y in x for y in req_phyla]))],
    #                   col="phylophlan_phyla", col_wrap=4, height=4, aspect=.7,
    #                   #hue=master_table.NumberofMembers,
    #                   #palette=color_pellet
    #                   )
    # x.map(sns.scatterplot, "reassembly_size_mb", "perc_coding_combined")
    # plt.show()

    # Fig6.1 boxplots
    # x = "Genome_Size"
    # y = "CDS"
    # plot_box = sns.boxplot(data=master_table, y=y, x=x,
    #                        # col=master_table.phylophlan_phyla, col_wrap=4, height=4, aspect=.7,
    #                        palette="rocket",
    #                        order=["Small", "Medium", "Big"])
    # add_stat_annotation(plot_box, data=master_table, x=x, y=y,
    #                     box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
    #                     order=["Small", "Medium", "Big"],
    #                     test="Mann-Whitney", text_format='star', loc='inside', verbose=2)
    # sns.set_palette("dark")
    # plt.show()

    # Fig 6.2 Scatter perc_coding vs genome size without marginal density but with with HUE and colorbar
    # x="reassembly_size_mb"
    # y="perc_coding_combined"
    # hue="reassembly_gc"
    # norm=plt.Normalize(master_table[hue].min(),master_table[y].max())
    # sm= plt.cm.ScalarMappable(cmap="RdYlBu",norm=norm)
    # ax = sns.scatterplot(data=master_table, x=x,y=y,
    #                     # col=master_table.phylophlan_phyla, col_wrap=4, height=4, aspect=.7,
    #                     # size=master_table.cluster_abundance,
    #                     hue=hue,
    #                     palette="RdYlBu",
    #                     # xlim=(0, 8)
    #                     )
    # ax.axvline(2, ls="--")
    # ax.get_legend().remove()
    # ax.figure.colorbar(sm)
    # ax.axvline(4, ls="--")
    # plt.show()

    # FIg 6.3 CUB Heatmap
    # sorted_master=master_table.sort_values(by=["Genome_Size","Cluster"])
    # universal_codon=pd.read_csv("Assembled_Metagenomes/codon_table.tsv",sep="\t", index_col="Codon").sort_values(by="AA")
    # new=codon_table.reindex(sorted_master.index)
    # new=new[universal_codon.index]
    # new["g_size"]=new.apply(lambda x: sorted_master.loc[x.name]["Genome_Size"],axis=1)
    # new["gs2"]=new.apply(lambda x: sorted_master.loc[x.name]["reassembly_size_mb"],axis=1)
    # #print_full(new.columns)
    # for i in new.columns:
    #     if not i in ["g_size","gs2"]:
    #         #relation=correlator(new,"g_size",i,"anova")
    #         relation2=correlator(new,"gs2",i,"pearson")
    #         #print(f"{i} & genome_size {relation}")
    # color_dict=dict(zip(["Small","Big","Medium"],["#332288", "#117733", "#44AA99"]))
    # row_cols=master_table["Genome_Size"].map(color_dict)
    # new=new.drop(["g_size","gs2"],axis=1)
    # x=sns.clustermap(new,row_colors=row_cols,col_cluster=False,row_cluster=False,xticklabels=True)
    # #plt.legend(x.row_colors)
    # plt.show()

    # Reduced_Table_plots
    # reduced_table = master_table[(master_table["reassembly_size_mb"] <= 2) | (master_table["reassembly_size_mb"] >= 5)]

    # x = sns.FacetGrid(data=members_table,hue="site_Temperature")
    # x.map(sns.scatterplot,"log_TP", "size")
    # x = sns.FacetGrid(data=non_streamlined_table, col="phylophlan_phyla", col_wrap=4, height=4, aspect=.7)
    # x.map(sns.histplot, "reassembly_size_mb", "site_TP")
    # x.set_titles("Non_Strem_lined")
    # print(f"{streamlined_table.shape}\t{non_streamlined_table.shape}")
    # g=sns.JointGrid(data=master_table, x="reassembly_size_mb", y="perc_coding_combined",hue="TP_lev")
    # g.plot_joint(sns.scatterplot)
    # g.plot_marginals(sns.histplot)
    # plt.legend()
    # plt.show()


def data_plots():
    level = "Genus"
    reduced_table = members_table.groupby(level).filter(lambda x: len(x) >= 2)
    # print(reduced_table.shape)
    mag_count = {}
    for i in reduced_table.index:
        mags = reduced_table.loc[i]["MAG"]
        org = reduced_table.loc[i][level]
        if not org in mag_count:
            mag_count[org] = []
            mag_count[org].append(mags)
        else:
            if mags not in mag_count[org]:
                mag_count[org].append(mags)
    for org in mag_count:
        mag_count[org] = len(mag_count[org])
    mag_count = pd.DataFrame.from_dict(mag_count, orient='index', columns=["lev_count"])
    mag_count = mag_count.reindex(reduced_table[level].value_counts().index)

    # sns.kdeplot(data=reduced_table, x="site_TP")
    # sns.scatterplot(x="site_Temperature", y="site_Elevation",hue=reduced_table["site_pH"],size=reduced_table["site_TP"], data=reduced_table)
    # corr = correlator(reduced_table, "site_Temperature", "site_Elevation", "pearson")

    fig, axs = plt.subplots(1, 5, gridspec_kw=dict(width_ratios=(7, 7, 7, 1, 1)), sharey=True)
    fig.subplots_adjust(left=0.2, wspace=0.01)

    sns.boxplot(y=reduced_table[level], x=reduced_table["site_pH"], ax=axs[0],
                order=reduced_table[level].value_counts().index, color="#AA4499")
    axs[0].set_xlabel("Preferred pH")
    # axs[0].set_ylabel("")

    sns.boxplot(y=reduced_table[level], x=reduced_table["site_Temperature"], ax=axs[1],
                order=reduced_table[level].value_counts().index, color="#AA4499")
    axs[1].set_ylabel("")
    axs[1].set_xlabel("Preferred Temperature")

    sns.boxplot(y=reduced_table[level], x=reduced_table["site_TP"], ax=axs[2],
                order=reduced_table[level].value_counts().index, color="#AA4499")
    axs[2].set_ylabel("")
    axs[2].set_xlabel("Preferred\nTP")

    sns.countplot(y=reduced_table[level], ax=axs[4], order=reduced_table[level].value_counts().index, color="#AA4499")
    axs[4].set_xlabel("No of\nLakes")
    axs[4].set_ylabel("")

    sns.barplot(x=mag_count["lev_count"], y=mag_count.index, ax=axs[3], color="#AA4499")
    axs[3].set_xlabel("MAGs")
    axs[3].set_ylabel("")

    # print(len(reduced_table[level].value_counts().index))

    # with mpl.rc_context(rc={'interactive': False}):
    #    plt.show(block=False)
    plt.ioff()
    plt.show()


transformed_abundace_profile = abundance_transform(members_abundance_profile, optional=read_counts)
load_from_table()
# make_ko_category_table()
data_plots2()
# data_plots()
