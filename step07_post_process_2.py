import os
import re

import pandas as pd
from pathlib import Path

import config as config

infolder = config.outfolder + "dRep__gff/"
kegg_db_folder = "/home/hb0358/PycharmProjects/mbs_general/some_test_files/"
kegg_module_defs_file = infolder + "kegg_tmp_modules_db.csv"
kegg_module_reacs_file = infolder + "kegg_tmp_reacs_db.csv"
reac_df = pd.DataFrame(columns=["KOs"])
reac_df.index.name = "reaction_id"
def_df = pd.DataFrame(columns=["definition", "reaction_kos", "def_kos", "reaction_kos_list", "orthology_section_KOs"])
def_df.index.name = "module_id"


def parse_reaction_entry(entry):
    entry_arr = str(entry).split("\n")
    ko_flag = 0
    reaction_id = entry_arr[0].split()[1]
    ko_list = ""
    for subentry in entry_arr:
        subentry_arr = subentry.strip().split()
        if subentry_arr[0].startswith("ORTHOLOG"):
            ko_list = ko_list + subentry_arr[1]
            ko_flag = 1
        if subentry_arr[0].startswith("K"):
            ko_list = ko_list + f",{subentry_arr[0]}"
    return reaction_id.strip(), ko_list.split(",")


def parse_module_entry(entry):
    entry_array = entry.strip().split("\n")
    module_id = entry_array[0].split()[1]
    name_line = " ".join(entry_array[1].split()[1:])
    def_line = ""
    # print(module_id,end="\t")
    reaction_flag = 0
    reaction_list = []
    reaction_kos = []
    orthology_sec = []
    orthology_sec_flag = 0
    def_flag = 0
    for subentry in entry_array:
        if not subentry.strip() == "":
            subentry_arr = subentry.strip().split()
            if orthology_sec_flag == 1:
                if subentry_arr[0].startswith("CLASS"):
                    orthology_sec_flag = 0
                else:
                    orthology_sec.extend([subentry_arr[0]])
            if def_flag == 1:
                if subentry_arr[0].startswith("ORTHO"):
                    def_flag = 0
                else:
                    def_line = def_line + " " + " ".join(subentry_arr)
            if reaction_flag == 1:
                if not subentry_arr[0].startswith("R"):
                    reaction_flag = 0
                else:
                    reaction_list.extend(subentry_arr[0].replace("+", ",").split(","))
            if subentry_arr[0].startswith("REACTION"):
                reaction_list.extend(subentry_arr[1].replace("+", ",").split(","))
                reaction_flag = 1
            if subentry_arr[0].startswith("DEFI"):
                def_line = " ".join(subentry_arr[1:])
                def_flag = 1
            if subentry_arr[0].startswith("ORTHOLOGY"):
                orthology_sec.extend([subentry_arr[1]])
                orthology_sec_flag = 1
    for reaction in reaction_list:
        r_kos = reac_df.loc[reaction, "KOs"]
        reaction_kos.extend(r_kos)
    ret_r_kos = list(sorted(set(reaction_kos)))
    return module_id, def_line, ret_r_kos, orthology_sec


def load_modules_kegg():
    global reac_df
    if Path(kegg_module_reacs_file).is_file():
        reac_df = pd.read_csv(kegg_module_reacs_file, index_col="reaction_id")
        reac_df["KOs"] = reac_df.KOs.apply(lambda x: x.split(";") if not type(x) == float else [])
    else:
        reac_db_file = kegg_db_folder + "module/reaction"
        if not Path(reac_db_file).is_file(): exit(f"reactions_db no fount at {reac_db_file} please extract the "
                                                  f"reaction/reaction file from the ligand subfolder of KEGG database"
                                                  f" and place in the module folder")
        modules_db = os.popen(f"cat {reac_db_file}").read().strip().split("///")
        for entry in modules_db:
            if not entry.strip() == "":
                reaction_id, ko_list_reac = parse_reaction_entry(entry.strip())
                reac_df.loc[reaction_id] = ([";".join(ko_list_reac)])
        reac_df.to_csv(kegg_module_reacs_file)

    global def_df
    if Path(kegg_module_defs_file).is_file():
        def_df = pd.read_csv(kegg_module_defs_file, index_col="module_id")
        def_df["reaction_kos_list"] = def_df.reaction_kos_list.apply(
            lambda x: x.split(";") if not type(x) == float else [])
        def_df["orthology_section_KOs"] = def_df.orthology_section_KOs.apply(
            lambda x: x.split(";") if not type(x) == float else [])
    else:
        module_db_file = kegg_db_folder + "module/module" + '.gz'
        if not Path(module_db_file).is_file(): exit(f"Modules_db no fount at {module_db_file}")
        modules_db = os.popen(f"zcat {module_db_file}").read().strip().split("///")
        for entry in modules_db:
            if not entry.strip() == "":
                module_id, def_line, react_kos, orthology = parse_module_entry(entry.strip())
                def_kos = list(sorted(set(list(filter(None, re.split('\W|\s', def_line))))))
                def_df.loc[module_id] = (
                [def_line, len(react_kos), len(def_kos), ";".join(react_kos), ";".join(orthology)])
        def_df.to_csv(kegg_module_defs_file)
    # print(reac_df.head())
    # print(def_df.head())
    return def_df


def calc_r1(mag_kos, module, module_completeness, mag):
    module_deff = def_df.loc[module, 'definition']
    nos = str(module_deff).count(" ") + 1
    nocb = str(module_deff).count("((") + str(module_deff).count("))")
    nosa = nos - nocb
    sub_factor = module_completeness * nosa / 100
    sub_factor = sub_factor.__round__()

    kos_in_module_deff = list(filter(None, re.split('\W|\s', module_deff)))
    kos_in_module_deff = def_df.loc[module, 'reaction_kos_list']
    ko_counter = 0
    for ko in mag_kos:
        if ko in kos_in_module_deff: ko_counter += 1
    red = ko_counter - sub_factor
    return red


def calc_r2(mag_kos, module):
    module_orthologs_list = def_df.loc[module, 'orthology_section_KOs']
    redunddancy_count = 0
    completeness_count = 0
    no_of_enz_cats = 0
    no_of_total_kos = 0
    for sub_list in module_orthologs_list:
        no_of_enz_cats = no_of_enz_cats + 1
        sub_list_arr = sub_list.split(',')
        no_of_total_kos = no_of_total_kos + len(sub_list_arr)
        intersection_list = list(set(sub_list_arr).intersection(set(mag_kos)))
        if len(intersection_list) > 1: redunddancy_count = redunddancy_count + len(intersection_list) - 1
        if len(intersection_list) > 0: completeness_count = completeness_count + 1
    no_of_possible_redundancies = no_of_total_kos - no_of_enz_cats
    redunddancy_percent = (
                                      redunddancy_count / no_of_possible_redundancies) * 100 if no_of_possible_redundancies > 0 else 0
    # print(f"Total no of Enz_catt: {no_of_enz_cats}\n"
    #      f"total no_of_kos: {no_of_total_kos}\n"
    #      f"Completness:\t{completeness_count}\n"
    #      f"redundancy:\t{redunddancy_count}\n"
    #      f"redundancy_percent:\t{redunddancy_percent}")
    return redunddancy_percent





def remove_minus(mag_kos, eval_str):
    all_minus_terms=[m for m in re.finditer(r"(-K\w+)", eval_str)]
    #print(eval_str)
    #print(all_minus_terms)
    re_redundundancy=0
    for min_term in all_minus_terms:
        min_term_term=min_term.group()
        if min_term_term.replace("-","") in mag_kos: re_redundundancy+=1
        eval_str=eval_str.replace(min_term_term,"")
    eval_str_new=eval_str
    #print(eval_str_new)
    #print(re_redundundancy)
    return eval_str_new,re_redundundancy


def remove_plus(mag_kos, module_def_line):
    plus_blocks=[m for m in re.finditer(r"(\+?K\w+\+[K\w+]+)", module_def_line)]
    complete_pl_blks=[]
    all_blks=[]
    for pl_vlc in plus_blocks:
        pl_vlc_def=pl_vlc.group()
        if pl_vlc_def.startswith("+") or pl_vlc_def.endswith("+"):pass
        else:
            #print(f"solving\t {pl_vlc_def}")
            pl_vlc_list=pl_vlc_def.split("+")
            intersection_list=list(set(pl_vlc_list).intersection(set(mag_kos)))
            if len(intersection_list)==len(pl_vlc_list):
                complete_pl_blks.append(pl_vlc_def)
            all_blks.append(pl_vlc_def)
    return complete_pl_blks, all_blks


def calc_r4(mag_kos,module):
    mag_kos=mag_kos.copy()
    module_def_line = str(def_df.loc[module, "definition"]).replace("-- ","")
    #print(module_def_line)
    redundacy = 0
    temp_kocount = 0
    if module_def_line.__contains__("-K"):module_def_line,redundacy=remove_minus(mag_kos,module_def_line)
    bracket_count=str(module_def_line).count("(")+str(module_def_line).count(")")
    #print(module_def_line)
    while bracket_count>0:
        groups = [m for m in re.finditer(r"(\([^(]*?\))", module_def_line)]
        for group in groups:
            group_str=group.group()
            #print(group_str)
            if module_def_line.__contains__("-K"):
                module_def_line, redundacy_2 = remove_minus(mag_kos, module_def_line)
                redundacy=redundacy+redundacy_2
            if len([m for m in re.finditer(r"(K\w+\+K\w+)", module_def_line)]) > 0:
                complete_pl_blk, all_pl_blks = remove_plus(mag_kos, module_def_line)
                for pl_block in all_pl_blks:
                    tmp_ko_nm = "KTP_" + str(temp_kocount).zfill(2)
                    temp_kocount += 1
                    if pl_block in complete_pl_blk: mag_kos.append(tmp_ko_nm)
                    module_def_line = module_def_line.replace(pl_block, tmp_ko_nm)
            group_match=group_str.replace("(","").replace(")","").replace("+"," ").split(",")
            temp_ko="KTM_"+str(temp_kocount).zfill(2)
            temp_kocount+=1
            for sub_group in group_match:
                if sub_group.__contains__(" "):
                    sub_group_match=sub_group.split(" ")
                    sub_temp_ko="KTX_"+str(temp_kocount).zfill(2)
                    temp_kocount+=1
                    sub_grp_match_intersection_list=list(set(sub_group_match).intersection(set(mag_kos)))
                    is_complete_inner=True if len(sub_grp_match_intersection_list)==len(sub_group_match) else False
                    if is_complete_inner:
                        group_match=[sub_temp_ko if item==sub_group else item for item in group_match]
                        mag_kos.append(sub_temp_ko)
                    else:
                        group_match = [sub_temp_ko if item == sub_group else item for item in group_match]
                    #print(sub_temp_ko,sub_group_match,sub_grp_match_intersection_list,is_complete_inner)
            intersection_list=list(set(group_match).intersection(set(mag_kos)))
            is_complete_OUTER=True if len(intersection_list)>0 else False
            if is_complete_OUTER:
                mag_kos.append(temp_ko)
                module_def_line=module_def_line.replace(group_str,temp_ko)
            else:
                module_def_line = module_def_line.replace(group_str, temp_ko)
            redundacy=redundacy+len(intersection_list)-1 if len(intersection_list)>0 else redundacy+0
            bracket_count = str(module_def_line).count("(") + str(module_def_line).count(")")
            #print(temp_ko,group_match,intersection_list,is_complete_OUTER)
    if len([m for m in re.finditer(r"(K\w+\+K\w+)", module_def_line)])>0:
        complete_pl_blk, all_pl_blks=remove_plus(mag_kos,module_def_line)
        for pl_block in all_pl_blks:
            tmp_ko_nm="KTP_"+str(temp_kocount).zfill(2)
            temp_kocount+=1
            if pl_block in complete_pl_blk: mag_kos.append(tmp_ko_nm)
            module_def_line = module_def_line.replace(pl_block, tmp_ko_nm)
    module_def_final_analysis_match=module_def_line.split(" ")
    completeness_arr=[]
    for block in module_def_final_analysis_match:
        if len(str(block))==6:
            if block in mag_kos:
                completeness_arr.append(True)
            else: completeness_arr.append(False)
        else:
            if block.__contains__(","):
                block_arr = block.split(",")
                sub_block_completness=[]
                for sub_block in block_arr:
                    if len(sub_block)==6:
                        if sub_block in mag_kos:
                            sub_block_completness.append(True)
                        else: sub_block_completness.append(False)
                    else:
                        print (module,module_def_line, sub_block, "ERROR")
    #final_intersection_list=list(set(module_def_final_analysis_match).intersection(set(mag_kos)))
    #is_complete_Final = False if False in completeness_arr or len(completeness_arr)==0 else True
    pc_complete=sum(completeness_arr)*100/len(completeness_arr) if len(completeness_arr)>0 else 0
    #print (module_def_final_analysis_match,final_intersection_list,is_complete_Final, redundacy)
    return redundacy, pc_complete


def calc_redundancy_in_modules():
    load_modules_kegg()
    module_completeness_file = infolder + "MicrobeAnnotator_out/metabolic_summary__module_completeness.tsv"
    mags_file = infolder + "post_cluster_data.csv"
    ko_profile_file = infolder + "ko_profile.tsv"
    ko_profile = pd.read_csv(ko_profile_file, sep="\t", index_col="Data").transpose()
    module_completeness_df = pd.read_csv(module_completeness_file, sep="\t", index_col='module')
    mags_df = pd.read_csv(mags_file, index_col='cluster_ID')
    modules_list = list(module_completeness_df.index)
    MAGs_list = list(mags_df.index)
    red_df = pd.DataFrame(columns=MAGs_list, index=modules_list)
    compl_df = pd.DataFrame(columns=MAGs_list, index=modules_list)
    for mag in MAGs_list:
        mag_kos_te = list(ko_profile[ko_profile[mag] == 1][mag].index)
        magrep = mags_df.loc[mag, 'Representative_Genome'] + ".faa.ko"
        for module in modules_list:
            #if mag == "EUL_095" and module=="M00155:
            r4, c1 = calc_r4(mag_kos_te, module)
            #print(module, magrep, r4, c1) if c1==True else print(end="")
            red_df.loc[module, mag] = r4
            compl_df.loc[module, mag] = c1
    red_df.index.name="module_id"
    red_df.to_csv(infolder + "kegg_module_redundancy.csv")
    compl_df.index.name = "module_id"
    compl_df.to_csv(infolder + "kegg_reCal_module_completeness.csv")


calc_redundancy_in_modules()
