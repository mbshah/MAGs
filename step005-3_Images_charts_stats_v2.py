import math
import os
import subprocess
import pandas as pd
import numpy as np
import random
import config as config
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats as sci
import statistics as st
from statannot import add_stat_annotation
import re

###Some import statements are for commented lines LEAVE THEM IN


infolder = config.outfolder + "/dRep__gff/"
clusters_file = infolder + "post_cluster_data.csv"
KO_file = infolder + "clusters_ko_summary.tsv"
pathways_file = infolder + "clusters_pathway_summary.tsv"
subSystem_file = infolder + "clusters_subSystem_summary.tsv"
codon_table_file = infolder + "post_cluster_data_codon_table.csv"
members_file = infolder + "members_table.csv"
master_table = pd.read_csv(clusters_file, index_col="cluster_ID")
#KOtable = pd.read_csv(KO_file, sep="\t", index_col="Cluster")
#pathwaysTable = pd.read_csv(pathways_file, sep="\t", index_col="Cluster")
#subSystemTable = pd.read_csv(subSystem_file, sep="\t", index_col="Cluster")
codon_table = pd.read_csv(codon_table_file, index_col="Cluster")
#kegg_db_folder = "../ancilary/kegg/"
members_table = pd.read_csv(members_file, index_col="bin")
streamlined_table = master_table[(master_table["reassembly_size_mb"] <= 2)]
non_streamlined_table = master_table[(master_table["reassembly_size_mb"] >= 4)]
intermediate_streamlined_table = master_table[
    (master_table["reassembly_size_mb"] < 4) & (master_table["reassembly_size_mb"] > 2)]
#members_abundance_profile = pd.read_csv(infolder + "cluster_abundance_profile_transformed_deseq2.tsv", sep="\t",
#                                        index_col="abundance").transpose()
members_abundance_profile = pd.read_csv(infolder + "cluster_abundance_profile_normalized_deseq2.tsv", sep="\t",
                                        index_col="abundance")
read_counts = pd.read_csv(config.read_counts_file, sep="\t", index_col="Sample")


# print(codon_table.index)

class SeabornFig2Grid():

    def __init__(self, seaborngrid, fig,  subplot_spec):
        self.fig = fig
        self.sg = seaborngrid
        self.subplot = subplot_spec
        if isinstance(self.sg, sns.axisgrid.FacetGrid) or \
            isinstance(self.sg, sns.axisgrid.PairGrid):
            self._movegrid()
        elif isinstance(self.sg, sns.axisgrid.JointGrid):
            self._movejointgrid()
        self._finalize()

    def _movegrid(self):
        """ Move PairGrid or Facetgrid """
        self._resize()
        n = self.sg.axes.shape[0]
        m = self.sg.axes.shape[1]
        self.subgrid = gridspec.GridSpecFromSubplotSpec(n,m, subplot_spec=self.subplot)
        for i in range(n):
            for j in range(m):
                self._moveaxes(self.sg.axes[i,j], self.subgrid[i,j])

    def _movejointgrid(self):
        """ Move Jointgrid """
        h= self.sg.ax_joint.get_position().height
        h2= self.sg.ax_marg_x.get_position().height
        r = int(np.round(h/h2))
        self._resize()
        self.subgrid = gridspec.GridSpecFromSubplotSpec(r+1,r+1, subplot_spec=self.subplot)

        self._moveaxes(self.sg.ax_joint, self.subgrid[1:, :-1])
        self._moveaxes(self.sg.ax_marg_x, self.subgrid[0, :-1])
        self._moveaxes(self.sg.ax_marg_y, self.subgrid[1:, -1])

    def _moveaxes(self, ax, gs):
        #https://stackoverflow.com/a/46906599/4124317
        ax.remove()
        ax.figure=self.fig
        self.fig.axes.append(ax)
        self.fig.add_axes(ax)
        ax._subplotspec = gs
        ax.set_position(gs.get_position(self.fig))
        ax.set_subplotspec(gs)

    def _finalize(self):
        plt.close(self.sg.fig)
        self.fig.canvas.mpl_connect("resize_event", self._resize)
        self.fig.canvas.draw()

    def _resize(self, evt=None):
        self.sg.fig.set_size_inches(self.fig.get_size_inches())


def make_ko_category_table():
    master_table2=master_table[master_table["gtdb_phyla"]=="Pseudomonadota"]
    small_genomes = master_table2[master_table2["Genome_Size"] == "Small"].index.values
    kosmall = KOtable[KOtable.index.isin(small_genomes)].drop(columns="no_of_KOs")
    big_genomes = master_table2[master_table2["Genome_Size"] == "Big"].index.values
    kobig = KOtable[KOtable.index.isin(big_genomes)].drop(columns="no_of_KOs")
    intermediate_genomes = master_table2[master_table2["Genome_Size"] == "Medium"].index.values
    kointermediate = KOtable[KOtable.index.isin(intermediate_genomes)].drop(columns="no_of_KOs")
    gene_data_table_file = infolder + "gene_differential_presence_proteo.tsv"
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

    #master_table["GC_lev"] = master_table.reassembly_gc.apply(
    #    lambda x: "low" if x < 45 else "high" if x > 55 else "Balanced")
    #master_table["Temp_lev"] = master_table.site_Temperature.apply(lambda x: "low" if x < 20 else "high")
    #master_table["pH_lev"] = master_table.site_pH.apply(
    #    lambda x: "Acidic" if x < 6.8 else "Basic" if x > 7.2 else "Neutral")
    #master_table["TP_lev"] = master_table.weighted_TP.apply(lambda x: "low" if x < 10 else "high")
    #master_table["Coding_lev"] = master_table.coding_perc.apply(lambda x: "Low" if x < 94 else "High")
    master_table["Genome_Size"] = master_table.reassembly_size_mb.apply(
        lambda x: "Small" if x <= 2 else "Big" if x >= 3.5 else "Medium")
    #master_table["TP_log"] = master_table.weighted_TP.apply(lambda x: math.log(x))
    # master_table["cluster_abundance"] = master_table.apply(
    #    lambda x: transformed_abundace_profile.loc[x.name], axis=1)
    streamlined_table = master_table[(master_table["reassembly_size_mb"] <= 2)]
    non_streamlined_table = master_table[(master_table["reassembly_size_mb"] >= 4)]
    intermediate_streamlined_table = master_table[
        (master_table["reassembly_size_mb"] < 4) & (master_table["reassembly_size_mb"] > 2)]
    #master_table["CDS"] = master_table.CDS.apply(lambda x: int(x.replace(",", "")))
    # print(f"{max(master_table['cluster_abundance'])}\t{min(master_table['cluster_abundance'])}")
    #master_table["gtdb_phyla"]=master_table.Kraken_Phyla.apply(lambda x:"Pseudomonadota" if x in ("Alphaproteobacteri",
    #                                                                                                 "Betaproteobacteria",
    #                                                                                                 "Gammaproteobacteria",
    #                                                                                                 "Deltaproteobacteria",
    #                                                                                                 "Oligoflexia"
    #                                                                                                 ) else x)
    #master_table.to_csv(infolder + "clusters_wphyp_sited_m.tsv", sep="\t")
    #x = 0
    # ot=open(f"{infolder}merged_fna_all_clusters.fna","w") ##to run CUB seperately
        # fna_file=infolder+data.name+"/"+data.name+"_gms.fna"
        # command=f"grep -v '^>' {fna_file}"
        # sequence_combined=subprocess.run([command], stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8").strip().replace("\n","")
        # sequence_combined_entry=f">{data.name}\n{sequence_combined}\n"
        # ot.write(sequence_combined_entry)
    # ot.close()


def abundance_transform(abundance_profile, method="hellinger", optional=None):
    abundance_profile2 = pd.DataFrame(abundance_profile).copy(deep=True).astype(float)
    if optional is None:
        optional = {}
    if method == "hellinger":
        read_counts = pd.DataFrame(optional)
        samples = abundance_profile.columns
        orgs = abundance_profile.index
        for sample in samples:
            total_reads = abundance_profile.sum(axis=0)[sample]
            total_reads = read_counts.loc[sample]["Cleaned"]
            for org in orgs:
                abundance_value = abundance_profile.loc[org][sample]
                transformed_abundance = ((abundance_value)*10000/ total_reads)
                abundance_profile2.loc[org][sample] = float(transformed_abundance)
                # print(f"{sample}\t{org}\t{abundance_value}\t{total_reads}\t{transformed_abundance}")
    # print_full(abundance_profile)
    return abundance_profile2


def correlator(data, fact1, fact2, method):
    if method == "anova":
        corr = sci.f_oneway(data[data[fact1] == 'Small'][fact2].values, data[data[fact1] == 'Medium'][fact2].values,
                            data[data[fact1] == 'Big'][fact2].values)
        #print(f'anova of {fact1} vs {fact2} in the table is {corr}')
        return corr.pvalue
    elif method == "pearson":
        corr = sci.pearsonr(data[fact1], data[fact2])
        print(f'Pearsons correlation of {fact1} against {fact2}: {corr}')
        return corr[1]




def rgb_to_hex(rgb):
    rgb = tuple([int(255 * val) for val in rgb])
    return '#' + ''.join([hex(val)[2:] for val in rgb]).upper()



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
    # org="Actinomycetota"
    # org_list=list(master_table.Kraken_Phyla.unique())
    # for org in org_list:
    #    if str(org).__contains__("Ol"):continue
    #    reduced_table=master_table[master_table["Kraken_Phyla"]==org]
    #    #print(reduced_table)
    #    #sns.kdeplot(reduced_table["reassembly_size_mb"],label=org+str(reduced_table.shape))
    # plt.legend()
    # plt.show()

    # (scatter_plot) phyla_sperated log TP vs reassembly Size
    # print(master_table.columns)
    # x = sns.FacetGrid(data=master_table,
    #                   col="gtdb_phyla", col_wrap=4, height=4, aspect=.7,
    #                   #hue=master_table.NumberofMembers,
    #                   #palette=color_pellet
    #                   )
    # x.map(sns.scatterplot, "reassembly_size_mb", "TP_log")
    # plt.show()

    # (scatter_plot)
    # sns.scatterplot(x=master_table.TP_log, y=master_table.reassembly_size_mb,
    #                 size=master_table.cluster_abundance, hue=master_table.NumberofMembers,
    #                 palette=color_pellet, sizes=(20,200)
    #                 )
    # sns.set_palette("dark")
    # plt.show()

    # (scatter_plot with members_table) ##NEWLY MODIFIED FOR NEW ABUNDANCE BOX plot
    # try_tab = pd.DataFrame(columns=["org", "tp", "org_size"])
    # j = 0
    # for x in members_table.index:
    #     abundance = members_table.abundance[x]
    #     tp = members_table.site_TP[x]
    #     org_size = members_table.org_size[x]
    #     org = members_table.Member[x]
    #     for i in range(1,round(abundance*10000)+1):
    #         j = j + 1
    #         try_tab.loc[j] = [org, tp, org_size]
    # ax = sns.scatterplot(x=members_table.log_TP,#.apply(lambda x: "low" if x < 10 else "high"),
    #                      y=members_table.org_size,# .apply(lambda x: "Small" if x <= 2 else "Big" if x >= 3.5 else "Medium"),
    #                      size=members_table.abundance,
    #                      # hue=master_table.NumberofMembers,
    #                      # palette=color_pellet,
    #                      sizes=(2, 200),
    #                      # order=["Small","Medium","Big"]
    #                      )
    #add_stat_annotation(ax, data=master_table, x=try_tab.tp.apply(lambda x: "low" if x < 10 else "high"), y=try_tab.org_size,
    #                    box_pairs=[("low", "high")],
    #                    #order=["Small", "Medium", "Big"],
    #                    test="Mann-Whitney", text_format='star', loc='inside', verbose=2)
    # sns.set_palette("dark")
    # # ax.get_legend().remove()
    # plt.legend(title="square root\ntransformed\nabundance")
    # plt.show()



    # (stripplot)
    # sns.stripplot(x=master_table.TP_lev, y=master_table.reassembly_size_mb,
    #               # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
    #               size=master_table.cluster_abundance,
    #               hue=master_table.NumberofMembers,
    #               palette=color_pellet,
    #               sizes=(10, 100)
    #               )
    # sns.set_palette("dark")
    # plt.show()

    # Scatter perc_coding vs genome size with marginal density
    # scatter_plot = sns.JointGrid(y=master_table.coding_perc, x=master_table.reassembly_size_mb,
    #                              # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
    #                              # size=master_table.cluster_abundance, hue=master_table.NumberofMembers,
    #                              palette=color_pellet,
    #                              xlim=(0, 8)
    #                              )
    # scatter_plot.plot_joint(sns.scatterplot)
    # scatter_plot.ax_joint.axvline(2, ls="--")
    # scatter_plot.plot_marginals(sns.kdeplot)
    # scatter_plot.ax_joint.axvline(4, ls="--")
    # plt.show()

    # scatter plot according to phyla, popular phyla
    # req_phyla=["Pseudomonadota","Actinomycetota","Bacteroidota","Cyanobacteriota"]
    # x="reassembly_size_mb"
    # y="coding_perc"
    # hue="reassembly_gc"
    # norm = plt.Normalize(master_table[hue].min(), master_table[y].max())
    # sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=norm)
    # ax = sns.FacetGrid(data=master_table[master_table['gtdb_phyla'].apply(lambda x: any([y in x for y in req_phyla]))],
    #                   col="gtdb_phyla", col_wrap=4, height=4, aspect=1,
    #                   hue=hue,
    #                   palette="RdYlBu",
    #                   )
    # ax.map(sns.scatterplot, x, y,edgecolor='#000000')
    # #ax.get_legend().remove()
    # #ax.figure.colorbar(sm)
    # plt.colorbar(sm)
    # plt.show()

    # boxplots
    # x = "Genome_Size"
    # y = "CDS"
    # plot_box = sns.boxplot(data=master_table, y=y, x=x,
    #                        # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
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
    # y="coding_perc"
    # hue="reassembly_gc"
    # norm=plt.Normalize(master_table[hue].min(),master_table[y].max())
    # sm= plt.cm.ScalarMappable(cmap="RdYlBu",norm=norm)
    # ax = sns.scatterplot(data=master_table, x=x,y=y,
    #                     # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
    #                     # size=master_table.cluster_abundance,
    #                     hue=hue,
    #                     palette="RdYlBu",
    #                     # xlim=(0, 8)
    #                      edgecolor="#000000",
    #                     )
    # ax.axvline(2, ls="--")
    # ax.get_legend().remove()
    # ax.set(xlabel='Genome Size (Mb)', ylabel='Percent Coding (%)')
    # ax.figure.colorbar(sm).set_label("GC Content")
    # ax.axvline(4, ls="--")
    # plt.show()




def fig_env_rel():
    color_pellet = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499",
                    "#882255"]
    level = "gtdb_genus"
    reduced_table = members_table.groupby(level).filter(lambda x: len(x) >= 2)
    # print(reduced_table.shape)
    mag_count = {}
    for i in reduced_table.index:
        mags = reduced_table.loc[i]["cluster"]
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
    #corr = correlator(reduced_table, "site_Temperature", "site_Elevation", "pearson")

    fig, axs = plt.subplots(1, 5, gridspec_kw=dict(width_ratios=(1,1,7, 7, 7)), sharey="row")
    fig.subplots_adjust(left=0.2, wspace=0.01)

    #plot1 ax=0
    sns.countplot(y=reduced_table[level], ax=axs[0], order=reduced_table[level].value_counts().index, color="#AA4499")
    axs[0].set_xlabel("No of\nLakes")
    # axs[0].set_ylabel("")

    # plot2 ax=1
    sns.barplot(x=mag_count["lev_count"], y=mag_count.index, ax=axs[1], color="#AA4499")
    axs[1].set_xlabel("MAGs")
    axs[1].set_ylabel("")

    # plot3 ax=2
    sns.boxplot(y=reduced_table[level], x=reduced_table["site_pH"], ax=axs[2],
                order=reduced_table[level].value_counts().index, color="#AA4499")
    axs[2].set_xlabel("pH")
    axs[2].set_ylabel("")

    # plot4 ax=3
    sns.boxplot(y=reduced_table[level], x=reduced_table["site_TP"], ax=axs[3],
                order=reduced_table[level].value_counts().index, color="#AA4499")
    axs[3].set_xscale("log")
    axs[3].set_ylabel("")
    axs[3].set_xlabel("Total Phosphorus")

    # plot5 ax=4
    sns.boxplot(y=reduced_table[level], x=reduced_table["site_Temperature"], ax=axs[4],
                order=reduced_table[level].value_counts().index, color="#AA4499")
    axs[4].set_ylabel("")
    axs[4].set_xlabel("Temperature")

    # plot6 ax=5 if enabled change the above line for defining subplots to: fig, axs = plt.subplots(1, 5, gridspec_kw=dict(width_ratios=(1,1,7, 7, 7)), sharey=True)
    #sns.boxplot(y=reduced_table[level], x=reduced_table["site_Elevation"], ax=axs[5],
    #            order=reduced_table[level].value_counts().index, color="#AA4499")
    #axs[5].set_ylabel("")
    #axs[5].set_xlabel("Elevation")


    # print(len(reduced_table[level].value_counts().index))

    # with mpl.rc_context(rc={'interactive': False}):
    #    plt.show(block=False)
    plt.ioff()
    plt.show()
    #print_full(members_table.groupby(level))


def fig_tpvsSize():
    # Fig 4.1.3 (scatter_plot with members_table) ##NEWLY MODIFIED FOR NEW ABUNDANCE BOX plot WEIGHTED (abundance/10000) takes time
    #the /10000 gives same result as without but is 1000 times faster. for the p value use R's function from library(weights) to perform weighted t.test
    #wtd.t.test(x=subset(data,tp_lev=="low")$size_mem, y=subset(data,tp_lev=="high")$size_mem, weight=subset(data,tp_lev=="low")$abundance_m, weighty = subset(data,tp_lev=="high")$abundance_m, samedata=F)
    try_tab = pd.DataFrame(columns=["org", "tp", "org_size"])
    j = 0
    for x in members_table.index:
        abundance = members_table.abundance_m[x]
        tp = members_table.site_TP[x]
        org_size = members_table.size_rep[x]
        org = x
        for i in range(1, round(abundance/10000) + 1):
            j = j + 1
            try_tab.loc[j] = [org, tp, org_size]
        #print(x)
    ax = sns.boxplot(x=try_tab.tp.apply(lambda x: "Low" if x < 10 else "High"),
                         y=try_tab.org_size,
                         order=["Low","High"]
                         )
    add_stat_annotation(ax, data=try_tab, x=try_tab.tp.apply(lambda x: "Low" if x < 10 else "High"),
                        y=try_tab.org_size,
                        box_pairs=[("Low","High")],
                        order=["Low","High"],
                        text_annot_custom=["pvalue=0.006"],
                        perform_stat_test=False, pvalues=[0.006],
                        text_format='simple', loc='inside', verbose=2)
    sns.set_palette("dark")
    ax.set_xlabel("Phosphorus Level")
    ax.set_ylabel("Genome Size (Mb)")
    plt.show()




def fig_GC_anal():
    # Fig 6.1 boxplots
    # x = "Genome_Size"
    # y = "reassembly_gc"
    # plot_box = sns.boxplot(data=master_table, y=y, x=x,
    #                        # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
    #                        palette="rocket",
    #                        order=["Small", "Medium", "Big"])
    # add_stat_annotation(plot_box, data=master_table, x=x, y=y,
    #                     box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
    #                     order=["Small", "Medium", "Big"],
    #                     test="t-test_welch", loc='inside', verbose=2,text_format='simple')
    # sns.set_palette("dark")
    # plot_box.set(xlabel='Genome Size', ylabel='GC (%)')
    # plt.show()

    # Fig 6.2 Scatter perc_coding vs genome size without marginal density but with with HUE and colorbar
    x="reassembly_size_mb"
    y="coding_perc"
    hue="reassembly_gc"
    norm=plt.Normalize(master_table[hue].min(),master_table[y].max())
    sm= plt.cm.ScalarMappable(cmap="RdYlBu",norm=norm)
    ax = sns.scatterplot(data=master_table, x=x,y=y,
                        # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                        # size=master_table.cluster_abundance,
                        hue=hue,
                        palette="RdYlBu",
                        # xlim=(0, 8)
                         edgecolor="#000000",
                        )
    ax.axvline(2, ls="--")
    ax.get_legend().remove()
    ax.set(xlabel='Genome Size (Mb)', ylabel='Percent Coding (%)')
    ax.figure.colorbar(sm).set_label("GC Content")
    ax.axvline(4, ls="--")
    plt.show()

    #scatter phyla wise
    req_phyla=["Pseudomonadota","Actinomycetota","Bacteroidota","Cyanobacteriota"]
    x="reassembly_size_mb"
    y="coding_perc"
    hue="reassembly_gc"
    norm = plt.Normalize(master_table[hue].min(), master_table[y].max())
    sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=norm)
    ax = sns.FacetGrid(data=master_table[master_table['gtdb_phyla'].apply(lambda x: any([y in x for y in req_phyla]))],
                      col="gtdb_phyla", col_wrap=2, height=4, aspect=1,col_order=req_phyla,
                      hue=hue,
                      palette="RdYlBu",
                      )
    ax.map(sns.scatterplot, x, y,edgecolor='#000000')
    ax.set_axis_labels("Genome Size (Mb)", "Coding Regions (%)")
    ax.set_titles(col_template="{col_name}")
    #ax.get_legend().remove()
    #plt.colorbar(sm).set_label("GC Content")
    plt.show()

def figure_7():
    # Fig 7.1 boxplots
    x = "Genome_Size"
    y = "sigma_factor"
    plot_box = sns.boxplot(data=master_table, y=y, x=x,
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Big"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
                        order=["Small", "Medium", "Big"],
                        test="t-test_welch", loc='inside', verbose=2,text_format='simple')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome Size', ylabel='Sigma Factor Genes')
    plt.show()

    # Fig 7.2 Scatter perc_coding vs genome size without marginal density but with with HUE and colorbar
    x="reassembly_size_mb"
    y="coding_perc"
    hue="sigma_factor"
    norm=plt.Normalize(master_table[hue].min(),master_table[y].max())
    sm= plt.cm.ScalarMappable(cmap="RdYlBu",norm=norm)
    ax = sns.scatterplot(data=master_table, x=x,y=y,
                        # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                        # size=master_table.cluster_abundance,
                        hue=hue,
                        palette="RdYlBu",
                        # xlim=(0, 8)
                         edgecolor="#000000",
                        )
    ax.axvline(2, ls="--")
    ax.get_legend().remove()
    ax.set(xlabel='Genome Size (Mb)', ylabel='Percent Coding (%)')
    ax.figure.colorbar(sm).set_label("No of Sigma Factor Genes")
    ax.axvline(4, ls="--")
    plt.show()

    #scatter phyla wise
    req_phyla=["Pseudomonadota","Actinomycetota","Bacteroidota","Cyanobacteriota"]
    x="reassembly_size_mb"
    y="coding_perc"
    hue="sigma_factor"
    norm = plt.Normalize(master_table[hue].min(), master_table[y].max())
    sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=norm)
    ax2 = sns.FacetGrid(data=master_table[master_table['gtdb_phyla'].apply(lambda x: any([y in x for y in req_phyla]))],
                      col="gtdb_phyla", col_wrap=4, height=4, aspect=1,
                      hue=hue,
                      palette="RdYlBu",
                      )
    ax2.map(sns.scatterplot, x, y,edgecolor='#000000')
    ax2.set_axis_labels("Genome Size (Mb)", "Coding Regions (%)")
    ax2.set_titles(col_template="{col_name}")
    #ax.get_legend().remove()
    #ax.figure.colorbar(sm)
    plt.colorbar(sm).set_label("No of Sigma Factor Genes")
    plt.show()


def figure_scatter_marginal_density():
    color_pellet = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499",
                    "#882255"]
    # Fig 5.1 Scatter perc_coding vs genome size with marginal density
    scatter_plot = sns.JointGrid(y=master_table.coding_perc, x=master_table.reassembly_size_mb,
                                 # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                                 # size=master_table.cluster_abundance, hue=master_table.NumberofMembers,
                                 palette=color_pellet,
                                 xlim=(0, 8)
                                 )
    scatter_plot.plot_joint(sns.scatterplot)
    scatter_plot.ax_joint.axvline(2, ls="--")
    scatter_plot.plot_marginals(sns.kdeplot)
    scatter_plot.ax_joint.axvline(4, ls="--")
    scatter_plot.set_axis_labels("Genome Size (Mb)","Coding Regions (%)")
    plt.show()


def figure_codons():
    from matplotlib.lines import Line2D
    # Codon Preference Graphs
    sorted_master = master_table.copy(deep=True).sort_values(by=["Genome_Size", "cluster_ID"])
    universal_codon = pd.read_csv("Assembled_Metagenomes/codon_table.tsv", sep="\t", index_col="CODONS").sort_values(
        by="AA")
    new = codon_table.reindex(sorted_master.index)
    new = new[universal_codon.index]
    new["g_size"] = new.apply(lambda x: sorted_master.loc[x.name]["Genome_Size"], axis=1)
    new["Genome_Size"] = new.apply(lambda x: sorted_master.loc[x.name]["reassembly_size_mb"], axis=1)

    # print_full(new[new["phyla"]=="Actinomycetota"])
    # Run and Print correlation tests
    # for i in new.columns:
    #     if not i in ["g_size","Genome_Size"]:
    #         #relation=correlator(new,"g_size",i,"anova")
    #         relation2=correlator(new,"Genome_Size",i,"pearson")
    #         #print(f"{i} & genome_size {relation2}")
    # color_dict=dict(zip(["Small","Big","Medium"],["#332288", "#117733", "#44AA99"]))
    # row_cols=master_table["Genome_Size"].map(color_dict)
    # x = sns.regplot(x=new['GCT'],y=new['Genome_Size'])
    # new=new.drop(["g_size","gs2"],axis=1)
    # x=sns.clustermap(new,row_colors=row_cols,col_cluster=False,row_cluster=False,xticklabels=True)
    # plt.legend(x.row_colors)
    # plt.show()
    new["phyla"] = new.apply(lambda x: sorted_master.loc[x.name]["gtdb_phyla"], axis=1)
    new_acti = new[new["phyla"] == "Actinomycetota"]
    new_bacti = new[new["phyla"] == "Bacteroidota"]
    new_proteo = new[new["phyla"] == "Pseudomonadota"]
    new_cyano = new[new["phyla"] == "Cyanobacteriota"]
    fig, axs = plt.subplots(3, 4, sharex="col")
    ftc = ["GC", "GG", "AC"]
    lwb = ["A", "T", "G", "C"]
    for axx in (0, 1, 2):
        for axy in (0, 1, 2, 3):
            codon = str(ftc[axx] + lwb[axy])
            sns.regplot(y=new[codon], x=new['Genome_Size'], color="black", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            sns.regplot(y=new_acti[codon], x=new_acti['Genome_Size'], color="#E69F00", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            sns.regplot(y=new_bacti[codon], x=new_bacti['Genome_Size'], color="#56B4E9", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            sns.regplot(y=new_proteo[codon], x=new_proteo['Genome_Size'], color="#009E73", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            sns.regplot(y=new_cyano[codon], x=new_cyano['Genome_Size'], color="#D55E00", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            axs[axx, axy].set_ylabel("")
            axs[axx, axy].set_xlabel("")
            axs[axx, axy].set_title(codon,fontsize=20)
            axs[axx, axy].set_ylim(0, )
    axs[0, 0].set_ylabel(f"Alanine",fontsize=20)
    axs[1, 0].set_ylabel(f"Glycine",fontsize=20)
    axs[2, 0].set_ylabel(f"Threonine",fontsize=20)
    axs[2, 0].set_xlabel("Genome Size (Mb)",fontsize=20)
    axs[2, 1].set_xlabel("Genome Size (Mb)",fontsize=20)
    axs[2, 2].set_xlabel("Genome Size (Mb)",fontsize=20)
    axs[2, 3].set_xlabel("Genome Size (Mb)",fontsize=20)
    legend_elements=[Line2D([0], [0], marker='o', markerfacecolor='black' , color='black',   label='All', markersize=5, alpha=0.2),
                     Line2D([0], [0], marker='o', markerfacecolor="#E69F00",color="#E69F00", label='Actinomycetota', markersize=5, alpha=0.5),
                     Line2D([0], [0], marker='o', markerfacecolor="#56B4E9",color="#56B4E9", label='Bacteroidota', markersize=5, alpha=0.5),
                     Line2D([0], [0], marker='o', markerfacecolor="#009E73",color="#009E73", label='Pseudomonadota', markersize=5, alpha=0.5),
                     Line2D([0], [0], marker='o', markerfacecolor="#D55E00",color="#D55E00", label='Cyanobacteriota',  markersize=5, alpha=0.5)
                     ]
    fig.legend(legend_elements,
               ["All","Actinomycetota","Bacteroidota","Pseudomonadota","Cyanobacteriota"],
               loc=8,
               ncol=5,
               prop={"size":20})
    plt.show()


def fig_size_box_plots():
    x = "Genome_Size"
    y = "reassembly_gc"
    print(f"{x}\t{y}")
    plot_box = sns.boxplot(data=master_table, y=y, x=x,
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Big"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
                        order=["Small", "Medium", "Big"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='simple')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome Size (Mb)', ylabel='GC (%)')
    plt.show()

    x = "Genome_Size"
    y = "sigma_factor"
    print(f"{x}\t{y}")
    plot_box = sns.boxplot(data=master_table, y=y, x=x,
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Big"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
                        order=["Small", "Medium", "Big"],
                        test="t-test_welch", loc='inside', verbose=2,text_format='simple')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome Size (Mb)', ylabel='No of Sigma Factor Genes')
    plt.show()

    x = "Genome_Size"
    y = "coding_perc"
    print(f"{x}\t{y}")
    plot_box = sns.boxplot(data=master_table, y=y, x=x,
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Big"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
                        order=["Small", "Medium", "Big"],
                        test="t-test_welch", loc='inside', verbose=2,text_format='simple')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome Size (Mb)', ylabel='Coding Regions (%)')
    plt.show()


def figure_phyla_scatter():
    cols_phyla=["Pseudomonadota","Actinomycetota","Bacteroidota","Cyanobacteriota"]
    rows_params=["site_TP"]
    fig, axs = plt.subplots(1, 4, sharex="col",sharey="row")
    for axx in range(len(rows_params)):
        for axy in range(len(cols_phyla)):
            y_param=rows_params[axx]
            x_param="size_rep"
            req_phyla=cols_phyla[axy]
            members_table2=members_table
            members_table2=members_table2.replace(["Alphaproteobacteria","Betaproteobacteria","Gammaroteobacteria","Deltaproteobacteria","Oligoflexia"],"Pseudomonadota")
            new_data=members_table2[members_table2['gtdb_phylum']==req_phyla]
            sns.regplot(data=new_data, y=y_param, x=x_param,
                        color="black", ax=axs[axy], ci=None,
                        scatter_kws={'alpha': 0.5, 's': 15})
            corr=sci.pearsonr(new_data[x_param],new_data[y_param])
            corr_v='%.02f'%corr[0]
            corr_p='%.02e'%corr[1]
            #axs[axx, axy].text(1,y=max(new_data[y_param]), s='%.08f'%corr[1])
            axs[axy].set_title(f"{req_phyla}\ncorr:{corr_v}  p-value:{corr_p}", fontsize=12) if axx==0 else axs[axy].set_title(f"corr:{corr_v}  p-value:{corr_p}", fontsize=12)
            axs[axy].set_ylabel("")
            axs[axy].set(yscale='log')
            axs[axy].set_xlabel("Genome Size (Mb)")
    axs[0].set_ylabel(f"Log of total phosphorus of extraction site")

    plt.show()

def figure_phyla_scatter2():
    cols_phyla=["Pseudomonadota","Actinomycetota","Bacteroidota","Cyanobacteriota"]
    rows_params=["reassembly_gc","sigma_factor"]
    fig, axs = plt.subplots(2, 4, sharex="col",sharey="row")

    for axx in range(len(rows_params)):
        for axy in range(len(cols_phyla)):
            hue_param=rows_params[axx]
            y_param = "coding_perc"
            x_param="reassembly_size_mb"
            req_phyla=cols_phyla[axy]
            new_data=master_table[master_table['gtdb_phyla']==req_phyla]
            norm = plt.Normalize(new_data[hue_param].min(), new_data[y_param].max())
            sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=norm)
            sns.scatterplot(data=new_data, y=y_param, x=x_param,
                        hue=hue_param, ax=axs[axx, axy],palette="RdYlBu",legend=None)
                        #ci=None,
                        #scatter_kws={'alpha': 0.5, 's': 10})
            corr=sci.pearsonr(new_data[x_param],new_data[y_param])
            #axs[axx, axy].text(1,y=max(new_data[y_param]), s='%.08f'%corr[1])
            axs[axx, axy].set_title(f"{req_phyla}\nc:{corr[0]}\np:{corr[1]}", fontsize=8) if axx==0 else axs[axx, axy].set_title(f"c:{corr[0]}\np:{corr[1]}", fontsize=8)
            axs[axx, axy].set_ylabel("Coding Regions (%)")if axy==0 else axs[axx, axy].set_ylabel("")
            axs[axx, axy].set_xlabel("Genome Size (Mb)")if axx==1 else axs[axx, axy].set_xlabel("")
            fig.colorbar(sm, ax=axs[axx, axy]).set_label("GC (%)" if hue_param=="reassembly_gc" else "Number of Sigma Factor Genes")
    plt.show()

def figure_phyla_scatter3():
    cols_phyla=["Pseudomonadota","Actinomycetota","Bacteroidota","Cyanobacteriota"]
    rows_params=["reassembly_gc","coding_perc","sigma_factor"]
    fig, axs = plt.subplots(3, 4, sharex="col",sharey="row")
    for axx in range(len(rows_params)):
        for axy in range(len(cols_phyla)):
            y_param=rows_params[axx]
            x_param="reassembly_size_mb"
            req_phyla=cols_phyla[axy]
            new_data=master_table[master_table['gtdb_phyla']==req_phyla]
            sns.regplot(data=new_data, y=y_param, x=x_param,
                        color="black", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.5, 's': 10})
            corr=sci.pearsonr(new_data[x_param],new_data[y_param])
            #axs[axx, axy].text(1,y=max(new_data[y_param]), s='%.08f'%corr[1])
            axs[axx, axy].set_title(f"{req_phyla}\nc:{corr[0]}\np:{corr[1]}", fontsize=8) if axx==0 else axs[axx, axy].set_title(f"c:{corr[0]}\np:{corr[1]}", fontsize=8)
            axs[axx, axy].set_ylabel("")
            axs[axx, axy].set_xlabel("Genome Size (MB)")if axx==2 else axs[axx, axy].set_xlabel("")
    axs[0, 0].set_ylabel(f"GC (%)")
    axs[1, 0].set_ylabel(f"Coding Regions (%)")
    axs[2, 0].set_ylabel(f"No of Sigma-factor Genes")

    plt.show()

#transformed_abundance_profile = abundance_transform(members_abundance_profile, optional=read_counts)
load_from_table()
#make_ko_category_table()
#data_plots2() #Draft set of plots, workspace to try plots
#fig_env_rel() #genus/phyl members bar plot with pH Temp Elevation with number of lakes and MAGs
#fig_tpvsSize() #box plot TP_Level #takes 5 minutes to make the plot because of the weights
#figure_5() #scatterplot with marginal density
#fig_GC_anal() #GC Analysis
#figure_7() #sigmafactor analysis
figure_codons()
##iteration for draft 3
#figure_codons()
#figure_phyla_scatter()
#figure_scatter_marginal_density()
#figure_phyla_scatter2()
#figure_phyla_scatter3()
#fig_size_box_plots()
