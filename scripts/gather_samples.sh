#!/bin/bash

#### Need fastq-scan and mash
echo $1 $2 $3
root_dir=${1%/*}
base_dir=$1
runtype=$2
app_dir=$3
single_end=''
[ "$runtype" = "SE" ] && single_end='true' || single_end='false'

skip_fastq_check="false"
min_basepairs=2241820
min_reads=7472
export min_coverage=10
min_proportion=0.5
min_genome_size=100000
max_genome_size=18040666
attempts=3
use_ena="false"
no_cache="false"

OPTS="--sample $root_dir/${1##*/} --min_basepairs ${min_basepairs} --min_reads ${min_reads} --min_proportion ${min_proportion}"
if [ "$single_end" = 'false' ]; then
    mv -t ${1} "${1}_1.fastq.gz" "${1}_2.fastq.gz"
    gzip -cd ${1}/${1##*/}_1.fastq.gz | fastq-scan > $base_dir/r1.json
    gzip -cd $base_dir/${1##*/}_2.fastq.gz | fastq-scan > $base_dir/r2.json
else
    mv "${1}.fastq.gz" ${1} 
    gzip -cd $base_dir/${1##*/}.fastq.gz | fastq-scan > $base_dir/r.json
fi

if ! reformat.sh in1=$base_dir/${1##*/}_1.fastq.gz in2=$base_dir/${1##*/}_2.fastq.gz qin='auto' out=/dev/null 2> $base_dir/${1##*/}-paired-end-error.txt; then
    ERROR=1
    echo "${1##*/} FASTQs contains an error. Please check the input FASTQs.
        Further analysis is discontinued." | \
    sed 's/^\\s*//' >> $base_dir/${1##*/}-paired-end-error.txt
else
    rm -f $base_dir/${1##*/}-paired-end-error.txt
fi
if $app_dir/scripts/check-fastqs.py --fq1 $base_dir/r1.json --fq2 $base_dir/r2.json ${OPTS}; then
    ERROR=1
fi
rm $base_dir/r1.json $base_dir/r2.json

GENOME_SIZE_OUTPUT="$base_dir/${1##*/}-genome-size.txt"
FASTQS="-r $base_dir/${1##*/}_1.fastq.gz $base_dir/${1##*/}_2.fastq.gz"

mash sketch -o test -k 31 -m 3 ${FASTQS} 2>&1 | \
grep "Estimated genome size:" | \
awk '{if($4){printf("%d\\n", $4)}} END {if (!NR) print "0"}' > ${GENOME_SIZE_OUTPUT}
rm -rf $base_dir/test.msh
seq_pair_sizes=`head -n1 ${GENOME_SIZE_OUTPUT}`
size_list=(${seq_pair_sizes//\\n/ })
ESTIMATED_GENOME_SIZE=${size_list[0]}

# Check if second pass is needed
if [ ${ESTIMATED_GENOME_SIZE} -gt "${max_genome_size}" ] || [ ${ESTIMATED_GENOME_SIZE} -lt "${min_genome_size}" ]; then
    # Probably high coverage, try increasing number of kmer copies to 10
    M="-m 10"
    if [ ${ESTIMATED_GENOME_SIZE} -lt "${min_genome_size}" ]; then
        # Probably low coverage, try decreasing the number of kmer copies to 1
        M="-m 1"
    fi
    mash sketch -o test -k 31 ${M} ${FASTQS} 2>&1 | \
        grep "Estimated genome size:" | \
        awk '{if($4){printf("%d\\n", $4)}} END {if (!NR) print "0"}' > ${GENOME_SIZE_OUTPUT}
    rm -rf $base_dir/test.msh
fi
 # Check final estimate
seq_pair_sizes=`head -n1 ${GENOME_SIZE_OUTPUT}`
size_list=(${seq_pair_sizes//\\n/ })
ESTIMATED_GENOME_SIZE=${size_list[0]}
if [ ${ESTIMATED_GENOME_SIZE} -gt "${max_genome_size}" ]; then
    rm ${GENOME_SIZE_OUTPUT}
    echo "${1##*/} estimated genome size (${ESTIMATED_GENOME_SIZE} bp) exceeds the maximum
            allowed genome size (${max_genome_size} bp). If this is unexpected, please
            investigate ${1##*/} to determine a cause (e.g. metagenomic, contaminants, etc...).
            Otherwise, adjust the --max_genome_size parameter to fit your need. Further analysis
            of ${1##*/} will be discontinued." | \
    sed 's/^\\s*//' > $base_dir/${1##*/}-genome-size-error.txt
elif [ ${ESTIMATED_GENOME_SIZE} -lt "${min_genome_size}" ]; then
    rm ${GENOME_SIZE_OUTPUT}
    echo "${1##*/} estimated genome size (${ESTIMATED_GENOME_SIZE} bp) is less than the minimum
            allowed genome size (${min_genome_size} bp). If this is unexpected, please
            investigate ${1##*/} to determine a cause (e.g. metagenomic, contaminants, etc...).
            Otherwise, adjust the --min_genome_size parameter to fit your need. Further analysis
            of ${1##*/} will be discontinued." | \
    sed 's/^\\s*//' > $base_dir/${1##*/}-genome-size-error.txt
else
    echo "${1##*/} sequences gathered"
fi

export genome_size=$ESTIMATED_GENOME_SIZE


