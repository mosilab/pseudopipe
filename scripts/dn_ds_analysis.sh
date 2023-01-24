#!/bin/bash

seqid=$1
working_dir=$2
folder_path="${2}/${1}"
id="MU${seqid:3}"


# printf "seqid=$1
# working_dir=$2
# folder_path=${folder_path}
# id=${id}\n\n"

# dnds_working_dir="${folder_path}/dn_ds_analysis"
# mkdir -p $working_dir
# "$root_base_dir/pan_files/$genome/faas/${id}.faa"

# myvar="$folder_path"/dn_ds_analysis/"${id}"_dfast.fna
# echo "$myvar"

grep -A 1 "MG" ${working_dir}/pan_files/pseudogenome/fnas/${id}.fna > ${folder_path}/dn_ds_analysis/${id}_dfast.fna
makeblastdb -in ${folder_path}/rt_output/prokka_annotation/${id}.linear.ffn -title ${folder_path}/rt_output/prokka_annotation/${id}.linear.ffn -dbtype nucl -out ${folder_path}/rt_output/prokka_annotation/${id}.linear.ffn.db
blastn -db ${folder_path}/rt_output/prokka_annotation/${id}.linear.ffn.db -query ${folder_path}/dn_ds_analysis/${id}_dfast.fna -outfmt 6 -out ${folder_path}/dn_ds_analysis/blastn_results.txt -max_target_seqs 1 -subject_besthit -perc_identity 100
awk '{print $2}' ${folder_path}/dn_ds_analysis/blastn_results.txt > ${folder_path}/dn_ds_analysis/found_ref.ids | grep -f - ${folder_path}/rt_output/prokka_annotation/${id}.linear.ffn > ${folder_path}/dn_ds_analysis/found_ref.fna

#all genes
total_num=$(wc -l ${folder_path}/pseudofinder/PF_output_sleuth/sleuth_report.csv | cut -d' ' -f1)
dnds_sum=$(awk -F"," 'FNR>1{x+=$23}END{print x}' ${folder_path}/pseudofinder/PF_output_sleuth/sleuth_report.csv)
average=$(echo "$dnds_sum $total_num" | awk '{print $1/($2-1)}')
echo "All($total_num) found average dn/ds = $average" 

#pseudogenes
cut -d' ' -f1 ${folder_path}/combined_pseudogenome_unique.txt | sed 's/>//' | grep -v "MG" | sed 's/PF_//' > ${folder_path}/pseudofinder/PF_output_sleuth/ids_sleuth_analysis_prokka_pseuds.txt
cat  ${folder_path}/pseudofinder/PF_output_sleuth/ids_sleuth_analysis_prokka_pseuds.txt ${folder_path}/dn_ds_analysis/found_ref.ids | uniq > ${folder_path}/pseudofinder/PF_output_sleuth/ids_sleuth_analysis_all_pseuds.ids

grep -f ${folder_path}/pseudofinder/PF_output_sleuth/ids_sleuth_analysis_all_pseuds.ids ${folder_path}/pseudofinder/PF_output_sleuth/sleuth.finder.blast | awk '{print $1}' | uniq | grep -f - ${folder_path}/pseudofinder/PF_output_sleuth/sleuth_report.csv > ${folder_path}/pseudofinder/PF_output_sleuth/pseuds_sleuth_report.csv

pseudogenes_total_num=$(wc -l ${folder_path}/pseudofinder/PF_output_sleuth/pseuds_sleuth_report.csv | cut -d' ' -f1)
pseudogenes_dnds_sum=$(awk -F"," 'FNR>1{x+=$23}END{print x}' ${folder_path}/pseudofinder/PF_output_sleuth/pseuds_sleuth_report.csv)
pseudogenes_average=$(echo "$pseudogenes_dnds_sum $pseudogenes_total_num" | awk '{print $1/$2}')
echo "Pseudogenes($pseudogenes_total_num) average dn/ds = $pseudogenes_average" 

#genes
grep -vf ${folder_path}/pseudofinder/PF_output_sleuth/ids_sleuth_analysis_all_pseuds.ids ${folder_path}/pseudofinder/PF_output_sleuth/sleuth.finder.blast | awk '{print $1}' | uniq | grep -f - ${folder_path}/pseudofinder/PF_output_sleuth/sleuth_report.csv > ${folder_path}/pseudofinder/PF_output_sleuth/genes_sleuth_report.csv
genes_total_num=$(wc -l ${folder_path}/pseudofinder/PF_output_sleuth/genes_sleuth_report.csv | cut -d' ' -f1)
genes_dnds_sum=$(awk -F"," 'FNR>1{x+=$23}END{print x}' ${folder_path}/pseudofinder/PF_output_sleuth/genes_sleuth_report.csv)
genes_average=$(echo "$genes_dnds_sum $genes_total_num" | awk '{print $1/$2}')
echo "Genes($genes_total_num) average dn/ds = $genes_average" 

#printf "id,all_total_num,all_dnds_sum,all_average,\
#pseudogenes_total_num,pseudogenes_dnds_sum,pseudogenes_average,\
#genes_total_num,genes_dnds_sum,genes_average\n"> dn_ds_analysis_summary.csv

printf "$id,$total_num,$dnds_sum,$average,\
$pseudogenes_total_num,$pseudogenes_dnds_sum,$pseudogenes_average,\
$genes_total_num,$genes_dnds_sum,$genes_average\n" >> "${working_dir}/secondary_analyses/dnds/dn_ds_analysis_summary.csv"
