#!/bin/bash

#Need quast and checkm-genome

base_dir=$1
seqid=${1##*/}
assembled_genome="$base_dir/rt_output/corrected"
cpus=$2
[[ $cpus ]] && cpus=$cpus || cpus=4
echo "Using $cpus threads"
assembly_qc_method=$3
[[ $assembly_qc_method ]] && u50=$assembly_qc_method || u50=''
# id="${4}${seqid:3}"
id="${4}"

GENOME_SIZE_OUTPUT="$base_dir/${seqid}-genome-size.txt"
seq_pair_sizes=`head -n1 ${GENOME_SIZE_OUTPUT}`
size_list=(${seq_pair_sizes//\\n/ })

# Assembly QC Parameters

contig_thresholds='0,1000,10000,100000,250000,1000000'
plots_format='pdf'


mkdir -p "$base_dir/assembly_qc"

# Choosing between QUAST ans U50 => Use this for production
if [[ "$assembly_qc_method" = "U50" ]]; then 
    echo "Starting U50 QC"
    u50_opts=" --circos --contig-thresholds ${contig_thresholds} --plots-format ${plots_format} --gene-finding --threads $cpus"
    # U50 bowtie2 mummer bedtools
    ./U50/U50_Terminal_SUBMISSION_SCRIPT.sh \
        "GCA_000013925.2_ASM1392v2_genomic.fna" \
        "${assembled_genome}/${id}.corrected.fasta" \
        "$base_dir/u50" \
        $u50_opts \
        --threads $cpus --seed 20 > "$base_dir/u50.stdout.txt" 2> "$base_dir/u50.stderr.txt"
    mv -f "$base_dir/u50" "$base_dir/assembly_qc/u50"
else
    # QUAST
    echo "Starting QUAST QC"
    GENOME_SIZE=${size_list[0]}
    est_ref_size=""
    if [ "${GENOME_SIZE}" = "0" ]; then
        est_ref_size="--est-ref-size ${GENOME_SIZE}"
    fi
    quast "${assembled_genome}/${id}.corrected.fasta" ${est_ref_size} \
        -o "$base_dir/quast" \
        --threads ${cpus} \
        #--glimmer \
        --circos \
        --contig-thresholds ${contig_thresholds} \
        --plots-format ${plots_format} > "$base_dir/quast.stdout.txt" 2> "$base_dir/quast.stderr.txt"
    mv "$base_dir/quast/quast.log" $base_dir
    mv "$base_dir/quast" "$base_dir/assembly_qc/quast"
fi

if [ -f "$base_dir/assembly_qc/u50/out_sorted_Assembly_Statistics.txt" ]; then
    u50_stats=$(cat "$base_dir/assembly_qc/u50/out_sorted_Assembly_Statistics.txt" | egrep "U*50(%| )")
    echo -e "\nU50 Statistics\n$u50_stats" | tee -a "$base_dir/assembly_qc_check.txt"
    completeness_line=$(tail -n1 $u50_stats | sed 's/: / /g; s/%//g')
    completeness_line_array=($completeness_line)
    completeness=${completeness_line_array[1]}
    if (( $(echo  "$completeness >= 95.0" | bc -l ) )); then
        echo "Assembled genome is reliable for further analysis. Genome completeness is >= 95" | tee -a "$base_dir/assembly_qc_check.txt"

    elif (( $(echo  "$completeness < 95.0" | bc -l ) )) && (( $(echo  "$calcuated_u50 >= 90.0" | bc -l ) )); then
        echo "Assembled genome is not reliable for further analysis. Genome completeness is < 95" | tee -a "$base_dir/assembly_qc_check.txt"

    else
        echo "Assembled genome may suffer from a diverse range of issues eg: substantial assembly or gene calling errors; lineage-specific gene loss" | tee -a "$base_dir/assembly_qc_check.txt"
    fi
else 
    echo "U50 was unsuccessful. Check u50.stderr.txt and u50.stdout.txt file for more detail"
fi

