samples=("EULd_A251SC" "EULd_S102LR" "EULd_S191UT" "EULd_S201PO" "EULd_S271TO" "EULd_Z011OB" "EULs_A031AN" "EULs_A093GO" "EULs_A111AU" "EULs_A132OS" "EULs_A152WI" "EULs_A191GI" "EULs_N031EI" "EULs_N051NO" "EULs_N092FA" "EULs_N121SA" "EULs_N141KU" "EULs_N163PI" "EULs_N261LU" "EULs_N263WI" "EULs_N282HE" "EULs_O022TU" "EULs_O072PA" "EULs_O102VI" "EULs_O111BA" "EULs_O121RA" "EULs_O151BU" "EULs_O241PL" "EULs_O251KR" "EULs_S022PA" "EULs_S042SE" "EULs_S153TR" "EULs_S171MA" "EULs_S193ME" "EULs_S222CL" "EULs_S281PM" "EULs_S301BU" "EULs_S301MM" "EULs_Z042CP" "EULs_Z053CV" "EULs_Z071SI" "EULs_Z111VV" "EULs_Z122OU" "EULs_Z221LA" "EULs_Z301WO")
folder="Assembled_Metagenomes/Bins/"

for y in "${samples[@]}"
do
	mkdir "${folder}../sorted/${y}"
done


for i in $(ls $folder)
do
	i_name=${i:0:11}
	for y in "${samples[@]}"
	do
		#echo "$i_name \t $y"
		if [ "$y" = "$i_name" ] ; then
		file1="$folder$i"
		file2="$folder/../sorted/$y/$i"
			#cp $file1 $file2
		fi
	done
done

for y in "${samples[@]}"
do
	python Scripts/merge_bins_to_overview.py "$folder/../sorted/$y/$y""_min1000_underscore.overview.txt" "$folder/../sorted/$y/$y""_metaspades3.12_das_tool_DASTool_scaffolds2bin.txt" 
done
