import pandas as pd
import config
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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

    fig, axs = plt.subplots(1, 5, gridspec_kw=dict(width_ratios=(1,1,7, 7, 7)), sharey="row", figsize=(10,15))
    fig.subplots_adjust(left=0.2, wspace=0.01)

    #plot1 ax=0
    sns.countplot(y=reduced_table[level], ax=axs[0], order=reduced_table[level].value_counts().index, color="#AA4499")
    axs[0].set_xlabel("No of\nLakes")
    axs[0].set_ylabel("Genus")

    # plot2 ax=1
    sns.barplot(x=mag_count["lev_count"], y=mag_count.index, ax=axs[1], color="#AA4499")
    axs[1].set_xlabel("No of\nMAGs")
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
    axs[3].set_xlabel("Natural log of total phosphorus (µg/l)")

    # plot5 ax=4
    sns.boxplot(y=reduced_table[level], x=reduced_table["site_Temperature"], ax=axs[4],
                order=reduced_table[level].value_counts().index, color="#AA4499")
    axs[4].set_ylabel("")
    axs[4].set_xlabel("Recorded Temperature (C)")

    # plot6 ax=5 if enabled change the above line for defining subplots to: fig, axs = plt.subplots(1, 5, gridspec_kw=dict(width_ratios=(1,1,7, 7, 7)), sharey=True)
    #sns.boxplot(y=reduced_table[level], x=reduced_table["site_Elevation"], ax=axs[5],
    #            order=reduced_table[level].value_counts().index, color="#AA4499")
    #axs[5].set_ylabel("")
    #axs[5].set_xlabel("Elevation")


    # print(len(reduced_table[level].value_counts().index))

    # with mpl.rc_context(rc={'interactive': False}):
    #    plt.show(block=False)
    plt.ioff()
    plt.tight_layout()
    plt.show()
    #print_full(members_table.groupby(level))

fig_env_rel()