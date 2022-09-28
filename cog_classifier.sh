#!/bin/bash

base_dir=$1
seqid=${base_dir##*/}
id=$2
cpus=$3
root_base_dir=${base_dir%/*}
# mkdir -p whole_genome/COG/$id
# COGclassifier -i whole_genome/faas/$id.faa -o whole_genome/COG/$id -t 6 \
#     1> "whole_genome/COG/cog_analysis.log" 2>  >> "whole_genome/COG/cog_analysis_error.log"
#     [ -s whole_genome/COG/cog_analysis_error.log ] && cat whole_genome/COG/cog_analysis_error.log || \
#     echo -e "COG analysis completed\n" | tee -a "whole_genome/COG/cog_analysis.log"

mkdir -p $root_base_dir/pan_files/pseudogenome/COG/$id
COGclassifier -i $root_base_dir/pan_files/pseudogenome/faas/$id.faa -o $root_base_dir/pan_files/pseudogenome/COG/$id -t $cpus \
    1> "$root_base_dir/pan_files/pseudogenome/COG/cog_analysis.log" 2> "$root_base_dir/pan_files/pseudogenome/COG/cog_analysis_error.log"
    [ -s $root_base_dir/pan_files/pseudogenome/COG/cog_analysis_error.log ] && \
        cat $root_base_dir/pan_files/pseudogenome/COG/cog_analysis_error.log || \
    echo -e "COG analysis completed\n" | tee -a "$root_base_dir/pan_files/pseudogenome/COG/cog_analysis.log"