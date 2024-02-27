import math
import pandas as pd
import config as config #project config script
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as sci

infolder = config.outfolder + "/dRep__gff/"
clusters_file = infolder + "post_cluster_data_Rver.csv"
KO_file = infolder + "ko_profile.tsv"
pathways_file = infolder + "clusters_pathway_summary.tsv"
subSystem_file = infolder + "clusters_subSystem_summary.tsv"
codon_table_file = infolder + "post_cluster_data_codon_table.csv"
members_file = infolder + "members_table_super_populated.csv"
codon_file_list = "../Assembled_Metagenomes/codon_table.tsv"
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
        lambda x: "Small" if x <= 2 else "Big" if x >= 3.5 else "Medium")
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


def figure_codons():
    from matplotlib.lines import Line2D
    # Codon Preference Graphs
    #sorted_master = master_table.copy(deep=True).sort_values(by=["Genome_Size", "cluster_ID"])
    sorted_members_table = members_table.copy(deep=True).sort_values(by=["size_rep", "bin"])
    universal_codon = pd.read_csv(codon_file_list, sep="\t", index_col="CODONS").sort_values(
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
    y_lab="Genome size (Mb)"
    new["phyla"] = new.apply(lambda x: sorted_members_table.loc[x.name]["gtdb_phylum"], axis=1)
    new_acti = new[new["phyla"] == "Actinomycetota"]
    new_bacti = new[new["phyla"] == "Bacteroidota"]
    new_proteo = new[new["phyla"] == "Pseudomonadota"]
    new_cyano = new[new["phyla"] == "Cyanobacteriota"]
    new_verru = new[new["phyla"] == "Verrucomicrobiota"]
    fig, axs = plt.subplots(3, 4, sharex="col",figsize=(10,10))
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
    plt.tight_layout()
    plt.show()

load_from_table()
figure_codons()
