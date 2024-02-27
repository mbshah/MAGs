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
clusters_file = infolder + "post_cluster_data_Rver.csv"
KO_file = infolder + "ko_profile.tsv"
pathways_file = infolder + "clusters_pathway_summary.tsv"
subSystem_file = infolder + "clusters_subSystem_summary.tsv"
codon_table_file = infolder + "post_cluster_data_codon_table.csv"
members_file = infolder + "members_table_super_populated.csv"
master_table = pd.read_csv(clusters_file, index_col="cluster_ID")
KOtable = pd.read_csv(KO_file, sep="\t", index_col="Data")
# pathwaysTable = pd.read_csv(pathways_file, sep="\t", index_col="Cluster")
# subSystemTable = pd.read_csv(subSystem_file, sep="\t", index_col="Cluster")
codon_table = pd.read_csv(codon_table_file, index_col="Cluster")
# kegg_db_folder = "../ancilary/kegg/"
members_table = pd.read_csv(members_file, index_col="bin")
streamlined_table = master_table[(master_table["reassembly_size_mb"] <= 2)]
non_streamlined_table = master_table[(master_table["reassembly_size_mb"] >= 4)]
intermediate_streamlined_table = master_table[
    (master_table["reassembly_size_mb"] < 4) & (master_table["reassembly_size_mb"] > 2)]
members_abundance_profile = pd.read_csv(infolder + "cluster_abundance_profile_normalized_deseq2.tsv", sep="\t",
                                        index_col="abundance").transpose()
# members_abundance_profile = pd.read_csv(config.outfolder + "cluster_abundance_profile_normalized_deseq2.tsv", sep="\t",
#                                        index_col="Cluster")
read_counts = pd.read_csv(config.read_counts_file, sep="\t", index_col="Sample")

class SeabornFig2Grid():
    def __init__(self, seaborngrid, fig, subplot_spec):
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
        self.subgrid = gridspec.GridSpecFromSubplotSpec(n, m, subplot_spec=self.subplot)
        for i in range(n):
            for j in range(m):
                self._moveaxes(self.sg.axes[i, j], self.subgrid[i, j])

    def _movejointgrid(self):
        """ Move Jointgrid """
        h = self.sg.ax_joint.get_position().height
        h2 = self.sg.ax_marg_x.get_position().height
        r = int(np.round(h / h2))
        self._resize()
        self.subgrid = gridspec.GridSpecFromSubplotSpec(r + 1, r + 1, subplot_spec=self.subplot)

        self._moveaxes(self.sg.ax_joint, self.subgrid[1:, :-1])
        self._moveaxes(self.sg.ax_marg_x, self.subgrid[0, :-1])
        self._moveaxes(self.sg.ax_marg_y, self.subgrid[1:, -1])

    def _moveaxes(self, ax, gs):
        # https://stackoverflow.com/a/46906599/4124317
        ax.remove()
        ax.figure = self.fig
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
    master_table2 = master_table  # [master_table["gtdb_phyla"]=="Pseudomonadota"]
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

    master_table["GC_lev"] = master_table.reassembly_gc.apply(
        lambda x: "low" if x < 45 else "high" if x > 55 else "Balanced")
    master_table["Temp_lev"] = master_table.site_Temperature.apply(lambda x: "low" if x < 20 else "high")
    master_table["pH_lev"] = master_table.site_pH.apply(
        lambda x: "Acidic" if x < 6.8 else "Basic" if x > 7.2 else "Neutral")
    master_table["TP_lev"] = master_table.weighted_TP.apply(lambda x: "low" if x < 10 else "high")
    master_table["Coding_lev"] = master_table.coding_perc.apply(lambda x: "Low" if x < 94 else "High")
    master_table["Genome_Size"] = master_table.reassembly_size_mb.apply(
        lambda x: "Small" if x <= 2 else "Big" if x >= 3.5 else "Medium")
    master_table["TP_log"] = master_table.weighted_TP.apply(lambda x: math.log(x))
    # master_table["cluster_abundance"] = master_table.apply(lambda x: transformed_abundace_profile.loc[x.name], axis=1)
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
    # master_table["CDS"] = master_table.CDS.apply(lambda x: int(x.replace(",", "")))
    # print(f"{max(master_table['cluster_abundance'])}\t{min(master_table['cluster_abundance'])}")
    # master_table["gtdb_phyla"]=master_table.Kraken_Phyla.apply(lambda x:"Pseudomonadota" if x in ("Alphaproteobacteri",
    #                                                                                                 "Betaproteobacteria",
    #                                                                                                 "Gammaproteobacteria",
    #                                                                                                 "Deltaproteobacteria",
    #                                                                                                 "Oligoflexia"
    #                                                                                                 ) else x)
    # master_table.to_csv(infolder + "clusters_wphyp_sited_m.tsv", sep="\t")
    # x = 0
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
                transformed_abundance = ((abundance_value) * 10000 / total_reads)
                abundance_profile2.loc[org][sample] = float(transformed_abundance)
                # print(f"{sample}\t{org}\t{abundance_value}\t{total_reads}\t{transformed_abundance}")
    # print_full(abundance_profile)
    return abundance_profile2


def correlator(data, fact1, fact2, method):
    if method == "anova":
        corr = sci.f_oneway(data[data[fact1] == 'Small'][fact2].values, data[data[fact1] == 'Medium'][fact2].values,
                            data[data[fact1] == 'Big'][fact2].values)
        # print(f'anova of {fact1} vs {fact2} in the table is {corr}')
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
    # add_stat_annotation(ax, data=master_table, x=try_tab.tp.apply(lambda x: "low" if x < 10 else "high"), y=try_tab.org_size,
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


def figure_codons():
    from matplotlib.lines import Line2D
    # Codon Preference Graphs
    #sorted_master = master_table.copy(deep=True).sort_values(by=["Genome_Size", "cluster_ID"])
    sorted_members_table = members_table.copy(deep=True).sort_values(by=["size_rep", "bin"])
    universal_codon = pd.read_csv("Assembled_Metagenomes/codon_table.tsv", sep="\t", index_col="CODONS").sort_values(
        by="AA")
    members_codon_table=pd.DataFrame(columns=codon_table.columns,index=members_table.index)
    for member in members_table.index:
        cluster_id=members_table.loc[member,"cluster"]
        for codon in codon_table.columns:
            codon_value=codon_table.loc[cluster_id, codon]
            members_codon_table.loc[member,codon]=codon_value
    new = members_codon_table.reindex(sorted_members_table.index)
    new = new[universal_codon.index]
    #new["g_size"] = new.apply(lambda x: sorted_master.loc[x.name]["Genome_Size"], axis=1)
    new["Measured_Param"] = new.apply(lambda x: sorted_members_table.loc[x.name]["size_rep"], axis=1)
    y_lab="Genome Size (Mb)"
    # print_full(new[new["phyla"]=="Actinomycetota"])
    # Run and Print correlation tests
    # for i in new.columns:
    #     if not i in ["g_size","Measured_Param"]:
    #         #relation=correlator(new,"g_size",i,"anova")
    #         relation2=correlator(new,"Measured_Param",i,"pearson")
    #         #print(f"{i} & Measured_Param {relation2}")
    # color_dict=dict(zip(["Small","Big","Medium"],["#332288", "#117733", "#44AA99"]))
    # row_cols=master_table["Measured_Param"].map(color_dict)
    # x = sns.regplot(x=new['GCT'],y=new['Measured_Param'])
    # new=new.drop(["g_size","gs2"],axis=1)
    # x=sns.clustermap(new,row_colors=row_cols,col_cluster=False,row_cluster=False,xticklabels=True)
    # plt.legend(x.row_colors)
    # plt.show()
    new["phyla"] = new.apply(lambda x: sorted_members_table.loc[x.name]["gtdb_phylum"], axis=1)
    new_acti = new[new["phyla"] == "Actinomycetota"]
    new_bacti = new[new["phyla"] == "Bacteroidota"]
    new_proteo = new[new["phyla"] == "Pseudomonadota"]
    new_cyano = new[new["phyla"] == "Cyanobacteriota"]
    new_verru = new[new["phyla"] == "Verrucomicrobiota"]
    fig, axs = plt.subplots(3, 4, sharex="col")
    ftc = ["GC", "GG", "AC"]
    lwb = ["A", "T", "G", "C"]
    for axx in (0, 1, 2):
        for axy in (0, 1, 2, 3):
            codon = str(ftc[axx] + lwb[axy])
            sns.regplot(y=new[codon], x=new['Measured_Param'], color="black", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            sns.regplot(y=new_acti[codon], x=new_acti['Measured_Param'], color="#E69F00", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            sns.regplot(y=new_bacti[codon], x=new_bacti['Measured_Param'], color="#56B4E9", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            sns.regplot(y=new_proteo[codon], x=new_proteo['Measured_Param'], color="#009E73", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            sns.regplot(y=new_cyano[codon], x=new_cyano['Measured_Param'], color="#D55E00", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            sns.regplot(y=new_verru[codon], x=new_verru['Measured_Param'], color="#AA4499", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.2, 's': 10})
            correlator(new_proteo, codon, 'Measured_Param', "pearson")
            axs[axx, axy].set_ylabel("")
            axs[axx, axy].set_xlabel("")
            axs[axx, axy].set_title(codon, fontsize=20)
            axs[axx, axy].set_ylim(0, )
    axs[0, 0].set_ylabel(f"Alanine", fontsize=20)
    axs[1, 0].set_ylabel(f"Glycine", fontsize=20)
    axs[2, 0].set_ylabel(f"Threonine", fontsize=20)
    axs[2, 0].set_xlabel(y_lab, fontsize=20)
    axs[2, 1].set_xlabel(y_lab, fontsize=20)
    axs[2, 2].set_xlabel(y_lab, fontsize=20)
    axs[2, 3].set_xlabel(y_lab, fontsize=20)
    legend_elements = [
        Line2D([0], [0], marker='o', markerfacecolor='black', color='black', label='All', markersize=5, alpha=0.2),
        Line2D([0], [0], marker='o', markerfacecolor="#E69F00", color="#E69F00", label='Actinomycetota', markersize=5,
               ),
        Line2D([0], [0], marker='o', markerfacecolor="#56B4E9", color="#56B4E9", label='Bacteroidota', markersize=5,
               ),
        Line2D([0], [0], marker='o', markerfacecolor="#009E73", color="#009E73", label='Pseudomonadota', markersize=5,
               ),
        Line2D([0], [0], marker='o', markerfacecolor="#D55E00", color="#D55E00", label='Cyanobacteriota', markersize=5,
               ),
        Line2D([0], [0], marker='o', markerfacecolor="#AA4499", color="#AA4499", label='Verrucomicrobiota', markersize=5,
               )
        ]
    fig.legend(legend_elements,
               ["All", "Actinomycetota", "Bacteroidota", "Pseudomonadota", "Cyanobacteriota", "Verrucomicrobiota"],
               loc=8,
               ncol=6,
               prop={"size": 15})
    plt.show()


def figure_box_plots():
    fig,axs=plt.subplots(2,3)

    x = "Genome_Size"
    y = "reassembly_gc"
    plot_box = sns.boxplot(data=master_table, y=y, x=x, ax=axs[0,0],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Big"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
                        order=["Small", "Medium", "Big"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome Size (Mb)', ylabel='GC (%)')
    axs[0,0].set_title('A:GC Content')

    x = "Genome_Size"
    y = "Total_Redundancy"
    plot_box = sns.boxplot(data=master_table, y=y, x=x,ax=axs[0,1],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Big"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
                        order=["Small", "Medium", "Big"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome Size (Mb)', ylabel='No of redundant\ngenes in a genome')
    axs[0,1].set_title('B:Redundant Genes')

    x = "Genome_Size"
    y = "coding_perc"
    plot_box = sns.boxplot(data=master_table, y=y, x=x,ax=axs[0,2],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Big"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
                        order=["Small", "Medium", "Big"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome Size (Mb)', ylabel='Coding Regions (%)')
    axs[0,2].set_title('C:Percent Coding in the MAGS')

    x = "Genome_Size"
    y = "KOs_for_ABC_Transporter"
    plot_box = sns.boxplot(data=master_table, y=y, x=x, ax=axs[1,0],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Big"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
                        order=["Small", "Medium", "Big"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome Size (Mb)', ylabel='No. of genes for\nABC transporters')
    axs[1,0].set_title('D:ABC Transporters Genes')

    x = "Genome_Size"
    y = "sigma_factor"
    plot_box = sns.boxplot(data=master_table, y=y, x=x, ax=axs[1,1],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Big"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
                        order=["Small", "Medium", "Big"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome Size (Mb)', ylabel='No. of sigma\nfactor genes')
    axs[1,1].set_title('E:Sigma Factor Genes')

    x = "Genome_Size"
    y = "kc_Cell_motility"
    plot_box = sns.boxplot(data=master_table, y=y, x=x, ax=axs[1,2],
                           # col=master_table.gtdb_phyla, col_wrap=4, height=4, aspect=.7,
                           palette="rocket",
                           order=["Small", "Medium", "Big"])
    add_stat_annotation(plot_box, data=master_table, x=x, y=y,
                        box_pairs=[("Small", "Big"), ("Small", "Medium"), ("Big", "Medium")],
                        order=["Small", "Medium", "Big"],
                        test="t-test_welch", loc='inside', verbose=2, text_format='star')
    sns.set_palette("dark")
    plot_box.set(xlabel='Genome Size (Mb)', ylabel='No. of genes for\ncell motility')
    axs[1,2].set_title('F:Cell Motility Genes')




    plt.show()


def figure_phyla_scatter():
    cols_phyla = ["Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota"]
    rows_params = ["site_TP"]
    fig, axs = plt.subplots(1, 4, sharex="col", sharey="row")
    for axx in range(len(rows_params)):
        for axy in range(len(cols_phyla)):
            y_param = rows_params[axx]
            x_param = "size_rep"
            req_phyla = cols_phyla[axy]
            members_table2 = members_table
            members_table2 = members_table2.replace(
                ["Alphaproteobacteria", "Betaproteobacteria", "Gammaroteobacteria", "Deltaproteobacteria",
                 "Oligoflexia"], "Pseudomonadota")
            new_data = members_table2[members_table2['gtdb_phylum'] == req_phyla]
            sns.regplot(data=new_data, y=y_param, x=x_param,
                        color="black", ax=axs[axy], ci=None,
                        scatter_kws={'alpha': 0.5, 's': 15})
            corr = sci.pearsonr(new_data[x_param], new_data[y_param])
            corr_v = '%.02f' % corr[0]
            corr_p = '%.02e' % corr[1]
            # axs[axx, axy].text(1,y=max(new_data[y_param]), s='%.08f'%corr[1])
            axs[axy].set_title(f"{req_phyla}\ncorr:{corr_v}  p-value:{corr_p}", fontsize=12) if axx == 0 else axs[
                axy].set_title(f"corr:{corr_v}  p-value:{corr_p}", fontsize=12)
            axs[axy].set_ylabel("")
            axs[axy].set(yscale='log')
            axs[axy].set_xlabel("Genome Size (Mb)")
    axs[0].set_ylabel(f"Log of total phosphorus of extraction site")

    plt.show()


def figure_phyla_scatter2():
    cols_phyla = ["Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota"]
    rows_params = ["reassembly_gc", "sigma_factor"]
    fig, axs = plt.subplots(2, 5, sharex="col", sharey="row")

    for axx in range(len(rows_params)):
        for axy in range(len(cols_phyla)):
            hue_param = rows_params[axx]
            y_param = "coding_perc"
            x_param = "reassembly_size_mb"
            req_phyla = cols_phyla[axy]
            new_data = master_table[master_table['gtdb_phyla'] == req_phyla]
            norm = plt.Normalize(new_data[hue_param].min(), new_data[y_param].max())
            sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=norm)
            sns.scatterplot(data=new_data, y=y_param, x=x_param,
                            hue=hue_param, ax=axs[axx, axy], palette="RdYlBu", legend=None)
            # ci=None,
            # scatter_kws={'alpha': 0.5, 's': 10})
            corr = sci.pearsonr(new_data[x_param], new_data[y_param])
            # axs[axx, axy].text(1,y=max(new_data[y_param]), s='%.08f'%corr[1])
            axs[axx, axy].set_title(f"{req_phyla}\nc:{corr[0]}\np:{corr[1]}", fontsize=10) if axx == 0 else axs[
                axx, axy].set_title(f"c:{corr[0]}\np:{corr[1]}", fontsize=10)
            axs[axx, axy].set_ylabel("Coding Regions (%)") if axy == 0 else axs[axx, axy].set_ylabel("")
            axs[axx, axy].set_xlabel("Genome Size (Mb)") if axx == 1 else axs[axx, axy].set_xlabel("")
            if hue_param == "reassembly_gc":
                fig.colorbar(sm, ax=axs[axx, axy]).set_label("GC (%)" if axy==4 else"")
            else:
                fig.colorbar(sm, ax=axs[axx, axy]).set_label("Number of Sigma Factor Genes" if axy == 4 else "")
    plt.show()


def figure_phyla_scatter3():
    cols_phyla = ["Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota"]
    rows_params = ["reassembly_gc", "coding_perc", "sigma_factor"]
    fig, axs = plt.subplots(3, 4, sharex="col", sharey="row")
    for axx in range(len(rows_params)):
        for axy in range(len(cols_phyla)):
            y_param = rows_params[axx]
            x_param = "reassembly_size_mb"
            req_phyla = cols_phyla[axy]
            new_data = master_table[master_table['gtdb_phyla'] == req_phyla]
            sns.regplot(data=new_data, y=y_param, x=x_param,
                        color="black", ax=axs[axx, axy], ci=None,
                        scatter_kws={'alpha': 0.5, 's': 10})
            corr = sci.pearsonr(new_data[x_param], new_data[y_param])
            # axs[axx, axy].text(1,y=max(new_data[y_param]), s='%.08f'%corr[1])
            axs[axx, axy].set_title(f"{req_phyla}\nc:{corr[0]}\np:{corr[1]}", fontsize=8) if axx == 0 else axs[
                axx, axy].set_title(f"c:{corr[0]}\np:{corr[1]}", fontsize=8)
            axs[axx, axy].set_ylabel("")
            axs[axx, axy].set_xlabel("Genome Size (MB)") if axx == 2 else axs[axx, axy].set_xlabel("")
    axs[0, 0].set_ylabel(f"GC (%)")
    axs[1, 0].set_ylabel(f"Coding Regions (%)")
    axs[2, 0].set_ylabel(f"No of Sigma-factor Genes")

    plt.show()


def members_scatter():
    cols_phyla = ["Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota"]
    cols_genus=["Limnohabitans","Rhodoferax","Rhodoluna","Planktophila","UBA2093","UBA3006","Sediminibacterium","Flavobacterium","SXYR01","Cyanobium","UBA953"]
    level="gtdb_phylum"
    req_taxa=cols_phyla
    rows_params = ["log_TP", "log_DP", "log_DN"]
    fig, axs = plt.subplots(3, len(req_taxa), sharex="col", sharey="row")
    for axx in range(len(rows_params)):
        for axy in range(len(req_taxa)):
            y_param = rows_params[axx]
            x_param = "size_rep"
            req_phyla = req_taxa[axy]
            new_data = members_table[members_table[level] == req_phyla]
            sns.regplot(data=new_data, y=y_param, x=x_param,
                         ax=axs[axx, axy], ci=None,
                        scatter_kws={'hue':"GC_content",'alpha': 0.5})
            corr = sci.pearsonr(new_data[x_param], new_data[y_param])
            # axs[axx, axy].text(1,y=max(new_data[y_param]), s='%.08f'%corr[1])
            axs[axx, axy].set_title(f"{req_phyla}\nc:{corr[0]}\np:{corr[1]}", fontsize=8) if axx == 0 else axs[
                axx, axy].set_title(f"c:{corr[0]}\np:{corr[1]}", fontsize=8)
            axs[axx, axy].set_ylabel("")
            axs[axx, axy].set_xlabel("GC Percent") if axx == 2 else axs[axx, axy].set_xlabel("")
    axs[0, 0].set_ylabel(f"TP")
    axs[1, 0].set_ylabel(f"DP")
    axs[2, 0].set_ylabel(f"DN")

    plt.show()


def nutrients_kdeplot():
    cols_param = ["ALL"]
    rows_params = ["log_DP", "log_TP", "log_DN"]
    fig, axs = plt.subplots(len(rows_params), len(cols_param))
    for axx in range(len(rows_params)):
        for axy in range(len(cols_param)):
            y_param = rows_params[axx]
            x_param = "GC_content"
            #req_phyla = cols_phyla[axy]
            #new_data = members_table[members_table['gtdb_phylum'] == req_phyla]
            new_data=members_table
            sns.histplot(data=new_data,x=y_param,ax=axs[axx])
            axs[axx].set_xlabel(y_param)
    plt.show()

# transformed_abundance_profile = abundance_transform(members_abundance_profile, optional=read_counts)
load_from_table()


# make_ko_category_table()
# data_plots2() #Draft set of plots, workspace to try plots
# figure_3() #genus/phyl members bar plot with pH Temp Elevation with number of lakes and MAGs
# figure_4() #box plot TP_Level #takes 5 minutes to make the plot because of the weights
# figure_5() #scatterplot with marginal density
# figure_6() #GC Analysis
# figure_7() #sigmafactor analysis
# figure_codons()
##iteration for draft 3
#figure_codons()
# figure_phyla_scatter()
# figure_scatter_marginal_density()
# figure_phyla_scatter2()
# figure_phyla_scatter3()


##iteration for draft 4
figure_codons()
#members_scatter()
#figure_phyla_scatter2()
#nutrients_kdeplot()
#figure_box_plots()
