#!/bin/bash
# Need bbmap(bbduk), fastq-scan, lighter, pigz

base_dir=$1
runtype=$2

#From general params
coverage=100
single_end=''
[ "$runtype" = "SE" ] && single_end='true' || single_end='false'
cpus=$3
[[ $cpus ]] && cpus=$cpus || cpus=4
echo "Using $cpus threads"
keep_all_files="false"
min_coverage=10
min_basepairs=2241820
min_reads=7472

GENOME_SIZE_OUTPUT="$base_dir/${1##*/}-genome-size.txt"
seq_pair_sizes=`head -n1 ${GENOME_SIZE_OUTPUT}`
size_list=(${seq_pair_sizes//\\n/ })

#QC Reads Parameters
# Illumina Reads
skip_qc="false"
skip_qc_plots="false"
skip_error_correction="false"
adapters=''
adapter_k=23
phix=''
phix_k=31
ktrim='r'
mink=11
hdist=1
tpe='t'
tbo='t'
qtrim='rl'
trimq=6
maq=10
minlength=35
ftm=5
tossjunk='t'
qout='33'
maxcor=1
sampleseed=42

fq=($base_dir/${1##*/}_1.fastq.gz $base_dir/${1##*/}_2.fastq.gz)

results=$base_dir/qc
mkdir -p $results
ERROR=0
GENOME_SIZE=${size_list[0]}
MIN_COVERAGE=$((min_coverage * GENOME_SIZE))
TOTAL_BP=$((coverage*GENOME_SIZE))

# echo $GENOME_SIZE $MIN_COVERAGE $TOTAL_BP ${fq[0]}

# options.ignore = [ '-genome-size.txt', extra]
if [ -z "${fq[1]}" ]; then single_end='true'; else single_end='false'; fi
if [ "$runtype" = 'assembly' ]; then is_assembly='true'; else is_assembly='false'; fi
if [ "$runtype" = 'assembly' ]; then qin='qin=33'; else qin='qin=auto'; fi
if [ "$adapters" ]; then adapters="path/adapters"; else adapters='adapters'; fi
if [ "$phix" ]; then phix="path/phix"; else phix='phix'; fi
if [ "$single_end" = 'true' ]; then adapter_opts=''; else adapter_opts="in2=${fq[1]} out2=$base_dir/adapter-r2.fq"; fi
if [ "$single_end" = 'true' ]; then phix_opts=''; else phix_opts="in2=$base_dir/adapter-r2.fq out2=$base_dir/phix-r2.fq"; fi
if [ "$single_end" = 'true' ]; then lighter_opts=''; else lighter_opts="-r $base_dir/phix-r2.fq"; fi
if [ "$single_end" = 'true' ]; then reformat_opts=''; else reformat_opts="in2=$base_dir/phix-r2.cor.fq out2=$base_dir/subsample-r2.fq"; fi

# echo $single_end $is_assembly $qin $adapters $phix $adapter_opts $phix_opts $lighter_opts $reformat_opts

# Illumina Reads
# Remove Adapters
bbduk.sh -Xmx1g \
    in=${fq[0]} out=$base_dir/adapter-r1.fq ${adapter_opts} \
    ref=${adapters} \
    k=${adapter_k} \
    ktrim=${ktrim} \
    mink=${mink} \
    hdist=${hdist} \
    tpe=${tpe} \
    tbo=${tbo} \
    threads=${cpus} \
    ftm=${ftm} \
    ${qin} ordered=t \
    overwrite=t \
    stats=$base_dir/bbduk-adapter.stdout.txt 2> $base_dir/bbduk-adapter.stderr.txt
# Remove PhiX
bbduk.sh -Xmx1g \
    in=$base_dir/adapter-r1.fq out=$base_dir/phix-r1.fq ${phix_opts} \
    ref=${phix} \
    k=${phix_k} \
    hdist=${hdist} \
    tpe=${tpe} \
    tbo=${tbo} \
    qtrim=${qtrim} \
    trimq=${trimq} \
    minlength=${minlength} \
    minavgquality=${maq} \
    ${qin} qout=${qout} \
    tossjunk=${tossjunk} \
    threads=${cpus} \
    overwrite=t \
    ordered=t stats=$base_dir/bbduk-phix.stdout.txt 2> $base_dir/bbduk-phix.stderr.txt

# Error Correction
if [[ "${runtype}" = "ont" ]]; then
    echo "Skipping error correction. Have a recommended ONT error corrector? Let me know!"
else
    if [ "${skip_error_correction}" = "false" ]; then
        lighter -od . -r $base_dir/phix-r1.fq ${lighter_opts} -od $base_dir -K 31 ${GENOME_SIZE} -maxcor 1 -zlib 0 -t ${cpus} 1> $base_dir/lighter.stdout.txt 2> $base_dir/lighter.stderr.txt
    else
        echo "Skipping error correction"
        ln -s $base_dir/phix-r1.fq $base_dir/phix-r1.cor.fq
        if [ "${single_end}" = "false" ]; then
            ln -s $base_dir/phix-r2.fq $base_dir/phix-r2.cor.fq
        fi
    fi
fi

# Reduce Coverage
if [ ${TOTAL_BP} -gt 0 ]; then
    if [[ "${runtype}" = "ont" ]]; then
        rasusa -i $base_dir/filt-r1.fq \
                -c ${coverage} \
                -g ${GENOME_SIZE} \
                -s ${sampleseed} 1> $base_dir/subsample-r1.fq 2> $base_dir/rasusa.stderr.txt
    else
        reformat.sh \
            in=$base_dir/phix-r1.cor.fq out=$base_dir/subsample-r1.fq ${reformat_opts} \
            samplebasestarget=${TOTAL_BP} \
            sampleseed=${sampleseed} \
            overwrite=t 1> $base_dir/reformat.stdout.txt 2> $base_dir/reformat.stderr.txt
    fi
else
    echo "Skipping coverage reduction"
    ln -s $base_dir/phix-r1.cor.fq $base_dir/subsample-r1.fq
    if [ "${single_end}" = "false" ]; then
        ln -s $base_dir/phix-r2.cor.fq $base_dir/subsample-r2.fq
    fi
fi

# Compress
if [ "${single_end}" = "false" ]; then
    pigz -p ${cpus} -c -n $base_dir/subsample-r1.fq > $results/${1##*/}_1.fastq.gz
    pigz -p ${cpus} -c -n $base_dir/subsample-r2.fq > $results/${1##*/}_2.fastq.gz
else
    pigz -p ${cpus} -c -n $base_dir/subsample-r1.fq > $results/${1##*/}.fastq.gz
fi
if [ "${keep_all_files}" = "false" ]; then
    # Remove intermediate FASTQ files
    rm $base_dir/*.fq
fi

# Quality stats before and after QC
mkdir -p $results/summary/
# fastq-scan
if [ "${single_end}" = "false" ]; then
    # Paired-End Reads
    gzip -cd ${fq[0]} | fastq-scan -g ${GENOME_SIZE} > $results/summary/${1##*/}_1_original.json
    gzip -cd ${fq[1]} | fastq-scan -g ${GENOME_SIZE} > $results/summary/${1##*/}_2_original.json
    gzip -cd $results/${1##*/}_1.fastq.gz | fastq-scan -g ${GENOME_SIZE} > $results/summary/${1##*/}_1_final.json
    gzip -cd $results/${1##*/}_2.fastq.gz | fastq-scan -g ${GENOME_SIZE} > $results/summary/${1##*/}_2_final.json
else
    # Single-End Reads
    gzip -cd ${fq[0]} | fastq-scan -g ${GENOME_SIZE} > $results/summary/${1##*/}_original.json
    gzip -cd $results/${1##*/}.fastq.gz | fastq-scan -g ${GENOME_SIZE} > $results/summary/${1##*/}_final.json
fi
# FastQC and NanoPlot
if [[ "${skip_qc_plots}" = "false" ]]; then
    if [[ "${runtype}" = "ont" ]]; then
        mkdir $results/summary/${1##*/}-original $results/summary/${1##*/}-final
        NanoPlot ${nanoplot_opts} \
            --threads ${cpus} \
            --fastq ${fq[0]} \
            --outdir $results/summary/${1##*/}-original/ \
            --prefix ${1##*/}-original_
        cp $results/summary/${1##*/}-original/${1##*/}-original_NanoPlot-report.html $results/summary/${1##*/}-original_NanoPlot-report.html
        tar -cvf - $results/summary/${1##*/}-original/ | pigz --best -p ${cpus} > $results/summary/${1##*/}-original_NanoPlot.tar.gz
        NanoPlot ${nanoplot_opts} \
            --threads ${cpus} \
            --fastq $results/${1##*/}.fastq.gz \
            --outdir $results/summary/${1##*/}-final/ \
            --prefix ${1##*/}-final_
        cp $results/summary/${1##*/}-final/${1##*/}-final_NanoPlot-report.html $results/summary/${1##*/}-final_NanoPlot-report.html
        tar -cvf - $results/summary/${1##*/}-final/ | pigz --best -p ${cpus} > $results/summary/${1##*/}-final_NanoPlot.tar.gz
        rm -rf $results/summary/${1##*/}-original/ $results/summary/${1##*/}-final/
    else
        if [ "${single_end}" = "false" ]; then
            # Paired-End Reads
            # ln -s ${fq[0]} $base_dir/${1##*/}_1_original.fastq.gz -f
            # ln -s ${fq[1]} $base_dir/${1##*/}_2_original.fastq.gz -f
            # ln -s $results/${1##*/}_1.fastq.gz $base_dir/${1##*/}_1_final.fastq.gz -f
            # ln -s $results/${1##*/}_2.fastq.gz $base_dir/${1##*/}_2_final.fastq.gz -f
            fastqc --noextract -f fastq -t ${cpus} ${fq[0]} \
                ${fq[1]} $results/${1##*/}_1.fastq.gz $results/${1##*/}_2.fastq.gz \
                -o $base_dir 1> $base_dir/fastqc.stdout.txt 2> $base_dir/fastqc.stderr.txt
        else
            # Single-End Reads
            # ln -s ${fq[0]} $base_dir/${1##*/}_original.fastq.gz -f
            # ln -s $results/${1##*/}.fastq.gz $base_dir/${1##*/}_final.fastq.gz -f
            fastqc --noextract -f fastq -t ${cpus} $base_dir/${1##*/}_original.fastq.gz $base_dir/${1##*/}_final.fastq.gz -o $base_dir \
                1> $base_dir/fastqc.stdout.txt 2> $base_dir/fastqc.stderr.txt
        fi
        mv $base_dir/*_fastqc.html $base_dir/*_fastqc.zip $results/summary/
    fi
fi
# Final QC check
gzip -cd $results/*.fastq.gz | fastq-scan -g ${GENOME_SIZE} > $base_dir/temp.json
FINAL_BP=$(grep "total_bp" $base_dir/temp.json | sed -r "s,[^0-9.]*,,g")
FINAL_READS=$(grep "read_total" $base_dir/temp.json | sed -r "s,[^0-9.]*,,g")
echo "Total basepairs: $FINAL_BP Total reads: $FINAL_READS"
rm $base_dir/temp.json
if [ ${FINAL_BP} -lt ${MIN_COVERAGE} ]; then
    ERROR=1
    echo "After QC, ${1##*/} FASTQ(s) contain ${FINAL_BP} total basepairs. This does
            not exceed the required minimum ${MIN_COVERAGE} bp (${min_coverage}x coverage). Further analysis 
            is discontinued." | \
    sed 's/^\\s*//' > $base_dir/${1##*/}-low-sequence-depth-error.txt
elif [ ${FINAL_BP} -lt ${min_basepairs} ]; then
    ERROR=1
    echo "After QC, ${1##*/} FASTQ(s) contain ${FINAL_BP} total basepairs. This does
            not exceed the required minimum ${min_basepairs} bp. Further analysis
            is discontinued." | \
    sed 's/^\\s*//' > $base_dir/${1##*/}-low-sequence-depth-error.txt
fi
if [ ${FINAL_READS} -lt ${min_reads} ]; then
    # Prevent ONT samples from being caught by this
    if [ ${FINAL_BP} -lt ${MIN_COVERAGE} ]; then
        ERROR=1
        echo "After QC, ${1##*/} FASTQ(s) contain ${FINAL_READS} total reads. This does
                not exceed the required minimum ${min_reads} reads count. Further analysis
                is discontinued." | \
        sed 's/^\\s*//' > $base_dir/${1##*/}-low-read-count-error.txt
    fi
fi
if [ "${is_assembly}" = "true" ]; then
    touch $results/reads-simulated-from-assembly.txt
fi
if [ "${ERROR}" -eq "1" ]; then
    if [ "${single_end}" == "false" ]; then
        mv $results/${1##*/}_1.fastq.gz $results/${1##*/}_1.error-fastq.gz
        mv $results/${1##*/}_2.fastq.gz $results/${1##*/}_2.error-fastq.gz
    else
        mv $results/${1##*/}.fastq.gz $results/${1##*/}.error-fastq.gz
    fi
else
    echo "QC done!!"
fi