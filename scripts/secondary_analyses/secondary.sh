#!/bin/bash

working_dir=$1
organism_tag=$2 
cpus=$3
app_dir=$4
srr_list=$5

# printf "\nworking_dir=$1
# organism_tag=$2
# cpus=$3
# app_dir=$4\n"

mkdir -p ${working_dir}/secondary_analyses

#roary analysis
echo -e "\nRoary analysis"
mkdir -p ${working_dir}/secondary_analyses/roary ${working_dir}/secondary_analyses/roary/roary_output
# conda activate "bacterial-genomics-tutorial"
roary -p $cpus -r -s -v -y -f ${working_dir}/secondary_analyses/roary ${working_dir}/pan_files/pseudogenome/gffs/*.gff \
    1> ${working_dir}/secondary_analyses/roary/stdout.txt 2> ${working_dir}/secondary_analyses/roary/stderr.txt
    [ -s ${working_dir}/secondary_analyses/roary/stderr.txt ] && cat ${working_dir}/secondary_analyses/roary/stderr.txt || \
    echo -e "Roary successful\n" | tee -a "${working_dir}/secondary_analyses/roary/roary.log"
mv -t $working_dir/secondary_analyses/roary/roary_output $working_dir/secondary_analyses/roary_**/* 
rm -fd $working_dir/secondary_analyses/roary_*

#post-roary analysis
echo -e "\nPost-roary analysis"
python "$app_dir/scripts/secondary_analyses/roary.py" $working_dir

### COG and KEGG
echo -e "\nPost-COG and -KEGG analysis"
python "$app_dir/scripts/secondary_analyses/cell_function.py" $working_dir $app_dir
python "$app_dir/scripts/secondary_analyses/met_function.py" $working_dir $app_dir

#dnds analysis
echo -e "\ndn/ds analysis"
mkdir -p ${working_dir}/secondary_analyses/dnds
touch ${working_dir}/secondary_analyses/dnds/dn_ds_analysis_summary.csv ${working_dir}/secondary_analyses/dnds/dnds_analysis.log
printf "id,all_total_num,all_dnds_sum,all_average,pseudogenes_total_num,pseudogenes_dnds_sum,pseudogenes_average,genes_total_num,genes_dnds_sum,genes_average\n"> \
    ${working_dir}/secondary_analyses/dnds/dn_ds_analysis_summary.csv

# for fasta in $(find ${working_dir} -maxdepth 1 -type d -name "SRR*"); do id=${fasta##*/}; source ${app_dir}/scripts/dn_ds_analysis.sh $id ${working_dir}; done | \
#     tee -a ${working_dir}/secondary_analyses/dnds/dnds_analysis.log
for fasta in ${srr_list[@]}; do id=${fasta##*/}; echo "$id"; source ${app_dir}/scripts/dn_ds_analysis.sh $id ${working_dir}; done | \
    tee -a ${working_dir}/secondary_analyses/dnds/dnds_analysis.log

# while read -r line; do
#     source ${app_dir}/scripts/dn_ds_analysis.sh $line ${working_dir}; 
# done < $input_file | \
#     tee -a ${working_dir}/secondary_analyses/dnds/dnds_analysis.log

python "${app_dir}/scripts/secondary_analyses/selection.py" $working_dir