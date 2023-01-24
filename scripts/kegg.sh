#!/bin/bash

base_dir=$1
seqid=${base_dir##*/}
id=$2
cpus=$3
if ! [[ -z "${5##*/}" ]]; then microbeannotatorDB_path=$5; else microbeannotatorDB_path=''; fi
root_base_dir=${base_dir%/*}

####### KEGG

if [[ "$microbeannotatorDB_path" != "''" ]]; then 
    for element in pseudogenome whole_genome functional_genome; do
        mkdir -p $root_base_dir/pan_files/${element}/KEGG/$id
        microbeannotator -i $root_base_dir/pan_files/${element}/faas/$id.faa -d $microbeannotatorDB_path -o $root_base_dir/pan_files/${element}/KEGG/$id -m diamond \
        -p $cpus --refine \
        1> "$root_base_dir/pan_files/${element}/KEGG/kegg_analysis.log" 2> "$root_base_dir/pan_files/${element}/KEGG/kegg_analysis_error.log"
        [ -s $root_base_dir/pan_files/${element}/KEGG/kegg_analysis_error.log ] && \
            cat $root_base_dir/pan_files/${element}/KEGG/kegg_analysis_error.log || \
        echo -e "KEGG analysis completed\n" | tee -a "$root_base_dir/pan_files/${element}/KEGG/kegg_analysis.log" #-t 3 
    done
fi