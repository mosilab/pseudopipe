#$/usr/bin/env perl
##$/bain/bash

# Need assembly-scan, flash, makeblastdb, shovill, skesa

base_dir=$1
seqid=${1##*/}
runtype=$2
cpus=$3
assembly_ram=$4
id="${5}${seqid:3}"
[[ $assembly_ram ]] && assembly_ram=$assembly_ram || assembly_ram=4
[ "$runtype" = "SE" ] && single_end='true' || single_end='false'
keep_all_files='false'
[[ $cpus ]] && cpus=$cpus || cpus=4
echo "Using $cpus threads"
min_genome_size=100000
max_genome_size=18040666

GENOME_SIZE_OUTPUT="$base_dir/${1##*/}-genome-size.txt"
seq_pair_sizes=`head -n1 ${GENOME_SIZE_OUTPUT}`
size_list=(${seq_pair_sizes//\\n/ })

results="$base_dir/qc"
fq=($results/${1##*/}_1.fastq.gz $results/${1##*/}_2.fastq.gz)

# Assemble Genome Parameters
shovill_assembler='skesa'
dragonflye_assembler='flye'
use_unicycler='false'
min_contig_len=500
min_contig_cov=2
contig_namefmt=''
shovill_opts=''
shovill_kmers=''
trim=''
no_stitch=''
no_corr=''
unicycler_mode='normal'
min_polish_size=10000
min_component_size=1000
min_dead_end_size=1000
no_miniasm='false'
no_polish='false'
no_rotate='false'
racon_steps=1
medaka_steps=0
medaka_model=''

#Shovill
[[ "$contig_namefmt" ]] && contig_namefmt=$contig_namefmt || contig_namefmt="${id}_%05d"
shovill_ram=$assembly_ram
[[ "$shovill_opts" ]] && opts="--opts '${shovill_opts}'" || opts=""
[[ "$shovill_kmers" ]] && kmers="--kmers '${shovill_kmers}'" || kmers=""
[[ "$no_stitch" ]] && nostitch="--nostitch" || nostitch=""
[[ "$no_corr" ]] && nocorr="--nocorr" || nocorr=""
[[ "$single_end" = 'false' ]] && shovill_mode="shovill --R1 ${fq[0]} --R2 ${fq[1]} ${nostitch}" || shovill_mode="shovill-se --SE ${fq[0]}"
shovill_opts="--assembler ${shovill_assembler} --depth 0 --noreadcorr ${opts} ${kmers} ${nocorr}"
# echo $contig_namefmt $opts $kmers $nostitch $nocorr $shovill_mode $shovill_opts

dragonflye_opts=''

# Merge Shovill/Dragonfly opts
[[ "$runtype" = "ont" ]] && assembler_wf="dragonflye" || assembler_wf="shovill"
[[ "$runtype" = "ont" ]] && assembler_mode="dragonflye --reads ${fq[0]}" || assembler_mode=$shovill_mode
[[ "$runtype" = "ont" ]] && assemnber_opts=$dragonflye_opts || assemnber_opts=$shovill_opts

# echo $assembler_wf $assembler_mode $assemnber_opts

OUTDIR="$base_dir/assembly"
mkdir -p $OUTDIR
GENOME_SIZE=${size_list[0]}

# Shovill or Dragonflye
${assembler_mode} --gsize ${GENOME_SIZE} \
    --outdir ${OUTDIR} \
    --force \
    --minlen ${min_contig_len} \
    --mincov ${min_contig_cov} \
    --namefmt "${contig_namefmt}" \
    --keepfiles \
    --cpus ${cpus} \
    --ram ${shovill_ram} ${assemnber_opts} > ${OUTDIR}/${assembler_wf}.stdout.txt 2> ${OUTDIR}/${assembler_wf}.stderr.txt
mv ${OUTDIR}/contigs.fa "${OUTDIR}/${id}.fna"

# Rename Graphs
if [ -f "${OUTDIR}/contigs.gfa" ]; then
        mv ${OUTDIR}/contigs.gfa ${OUTDIR}/${shovill_assembler}-unpolished.gfa
    elif [ -f "${OUTDIR}/contigs.fastg" ]; then
        mv ${OUTDIR}/contigs.fastg ${OUTDIR}/${shovill_assembler}-unpolished.gfa
    elif [ -f "${OUTDIR}/contigs.LastGraph" ]; then
        mv ${OUTDIR}/contigs.LastGraph ${OUTDIR}/${shovill_assembler}-unpolished.gfa
    fi
if [ -f "${OUTDIR}/flye-info.txt" ]; then
    mv ${OUTDIR}/flye-info.txt ${OUTDIR}/flye.log
fi

# Check quality of assembly
TOTAL_CONTIGS=`grep -c "^>" ${OUTDIR}/${id}.fna || true`
# echo $TOTAL_CONTIGS
touch "total_contigs_${TOTAL_CONTIGS}"
if [ $TOTAL_CONTIGS -gt 0 ]; then
    assembly-scan "${OUTDIR}/${id}.fna" > ${OUTDIR}/${1##*/}.json 2> ${OUTDIR}/assembly-scan.stderr.txt
    TOTAL_CONTIG_SIZE=$(grep "total_contig_length" ${OUTDIR}/${1##*/}.json | tr -d -c 0-9)
    # echo $TOTAL_CONTIG_SIZE
    if [ "${TOTAL_CONTIG_SIZE}" -lt ${min_genome_size} ]; then
        mv ${OUTDIR}/${1##*/}.fna ${OUTDIR}/${1##*/}-error.fna
        mv ${OUTDIR}/${1##*/}.json ${OUTDIR}/${1##*/}-error.json
        echo "${1##*/} assembled size (${TOTAL_CONTIG_SIZE} bp) is less than the minimum allowed genome
                size (${min_genome_size} bp). If this is unexpected, please investigate ${1##*/} to
                determine a cause (e.g. metagenomic, contaminants, etc...) for the poor assembly.
                Otherwise, adjust the --min_genome_size parameter to fit your need. Further assembly
                based analysis of ${1##*/} will be discontinued." | \
        sed 's/^\\s*//' > $base_dir/${1##*/}-assembly-error.txt
    # else
        # Make BLASTDB
        # mkdir -p $base_dir/blastdb
        # cat "${OUTDIR}/${id}.fna" | makeblastdb -dbtype "nucl" -title "Assembled contigs for ${1##*/}" -out $base_dir/blastdb/${1##*/}
    fi
else
    mv "${OUTDIR}/${id}.fna" "${OUTDIR}/${id}-error.fna"
    echo "${1##*/} assembled successfully, but 0 contigs were formed. Please investigate
            ${1##*/} to determine a cause (e.g. metagenomic, contaminants, etc...) for this
            outcome. Further assembly-based analysis of ${1##*/} will be discontinued." | \
    sed 's/^\\s*//' > $base_dir/${1##*/}-assembly-error.txt
fi
# Cleanup and compress
if [ "${keep_all_files}" = "false" ]; then
    # Remove intermediate files
    rm -rfv ${OUTDIR}/shovill.bam* ${OUTDIR}/shovill-se.bam* ${OUTDIR}/flash.extendedFrags*  \
            ${OUTDIR}/flash.notCombined* ${OUTDIR}/*.fq.gz ${OUTDIR}/00*.gfa \
            ${OUTDIR}/flye/ ${OUTDIR}/flye.fasta* ${OUTDIR}/raven/  \
            ${OUTDIR}/raven.fasta* ${OUTDIR}/raven.cereal ${OUTDIR}/miniasm/ ${OUTDIR}/miniasm.fasta* \
            ${OUTDIR}/megahit/ ${OUTDIR}/megahit.fasta* \
            ${OUTDIR}/velvet.fasta* ${OUTDIR}/velvet/
fi
find ${OUTDIR} -maxdepth 1 -name "*.log" | xargs -I {} mv {} "$base_dir"
mv -t $base_dir "*.txt" "*.log" "total_contigs*"