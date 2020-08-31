curation_file = "Assembled_Metagenomes//Curated_MAGs.tsv.txt"
raw_fastq="/mnt/xio/botany/mbs_workfolder/wmg_raw/"
outfolder= "outfile/"
fasta_folder="s_fasta/"
method="average"
postfix="_dastool_90_10"
compl_col = 11  # zero index
min_compl =90
max_compl = 100
contam_col = 12  # zero index
min_contam = 0
max_contam = 10

files_dict={
"EULd_A251SC":"Assembled_Metagenomes/refined_bins/EULd_A251SC.tsv",
"EULd_S102LR":"Assembled_Metagenomes/refined_bins/EULd_S102LR.tsv",
"EULd_S191UT":"Assembled_Metagenomes/refined_bins/EULd_S191UT.tsv",
"EULd_S201PO":"Assembled_Metagenomes/refined_bins/EULd_S201PO.tsv",
"EULd_S271TO":"Assembled_Metagenomes/refined_bins/EULd_S271TO.tsv",
"EULd_Z011OB":"Assembled_Metagenomes/refined_bins/EULd_Z011OB.tsv",
"EULs_A031AN":"Assembled_Metagenomes/refined_bins/EULs_A031AN.tsv",
"EULs_A093GO":"Assembled_Metagenomes/refined_bins/EULs_A093GO.tsv",
"EULs_A111AU":"Assembled_Metagenomes/refined_bins/EULs_A111AU.tsv",
"EULs_A132OS":"Assembled_Metagenomes/refined_bins/EULs_A132OS.tsv",
"EULs_A152WI":"Assembled_Metagenomes/refined_bins/EULs_A152WI.tsv",
"EULs_A191GI":"Assembled_Metagenomes/refined_bins/EULs_A191GI.tsv",
"EULs_N051NO":"Assembled_Metagenomes/refined_bins/EULs_N051NO.tsv",
"EULs_N092FA":"Assembled_Metagenomes/refined_bins/EULs_N092FA.tsv",
"EULs_N141KU":"Assembled_Metagenomes/refined_bins/EULs_N141KU.tsv",
"EULs_N282HE":"Assembled_Metagenomes/refined_bins/EULs_N282HE.tsv",
"EULs_O111BA":"Assembled_Metagenomes/refined_bins/EULs_O111BA.tsv",
"EULs_O121RA":"Assembled_Metagenomes/refined_bins/EULs_O121RA.tsv",
"EULs_O151BU":"Assembled_Metagenomes/refined_bins/EULs_O151BU.tsv",
"EULs_O241PL":"Assembled_Metagenomes/refined_bins/EULs_O24PL.tsv",
"EULs_S022PA":"Assembled_Metagenomes/refined_bins/EULs_S022PA.tsv",
"EULs_S042SE":"Assembled_Metagenomes/refined_bins/EULs_S043SE.tsv",
"EULs_S153TR":"Assembled_Metagenomes/refined_bins/EULs_S153TR.tsv",
"EULs_S171MA":"Assembled_Metagenomes/refined_bins/EULS_S171MA.tsv",
"EULs_S193ME":"Assembled_Metagenomes/refined_bins/EULs_S193ME.tsv",
"EULs_S301MM":"Assembled_Metagenomes/refined_bins/EULs_S301MM.tsv",
"EULs_Z042CP":"Assembled_Metagenomes/refined_bins/EULs_Z042CP.tsv",
"EULs_Z071SI":"Assembled_Metagenomes/refined_bins/EULs_Z071SI.tsv",
"EULs_Z111VV":"Assembled_Metagenomes/refined_bins/EULs_Z111VV.tsv",
"EULs_Z122OU":"Assembled_Metagenomes/refined_bins/EULs_Z122OU.tsv",
"EULs_Z221LA":"Assembled_Metagenomes/refined_bins/EULs_Z221LA.tsv",
"EULs_Z301WO":"Assembled_Metagenomes/refined_bins/EULs_Z301WO.tsv",
"EULs_N031EI":"Assembled_Metagenomes/refined_bins/N031EI.tsv",
"EULs_N121SA":"Assembled_Metagenomes/refined_bins/N121SA.tsv",
"EULs_N163PI":"Assembled_Metagenomes/refined_bins/N163PI.tsv",
"EULs_N261LU":"Assembled_Metagenomes/refined_bins/N261LU.tsv",
"EULs_O022TU":"Assembled_Metagenomes/refined_bins/O022TU.tsv",
"EULs_O072PA":"Assembled_Metagenomes/refined_bins/O072PA.tsv",
"EULs_O102VI":"Assembled_Metagenomes/refined_bins/O102VI.tsv",
"EULs_O251KR":"Assembled_Metagenomes/refined_bins/O251KR.tsv",
"EULs_S222CL":"Assembled_Metagenomes/refined_bins/S222CL.tsv",
"EULs_S281PM":"Assembled_Metagenomes/refined_bins/EULs_S281PM.tsv",
"EULs_S301BU":"Assembled_Metagenomes/refined_bins/S301BU.tsv",
"EULs_Z053CV":"Assembled_Metagenomes/refined_bins/Z053CV.tsv"
}


def filter_check(upper, lower, value):
    if float(lower) <= float(value) <= float(upper):
        return True
    else:
        return False


def filterok(a, b):
    if a == True and b == True:
        return True
    else:
        return False