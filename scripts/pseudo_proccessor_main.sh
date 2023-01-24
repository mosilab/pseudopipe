#!/bin/bash
# conda activate "bacterial-genomics-tutorial"

cpus=$2
[[ $cpus ]] && cpus=$cpus || cpus=4
# echo "Using $cpus threads"
base_dir=$1
seqid=${1##*/}
id="$5"
dnds=$6
scaffold=$7
predictor_type=$8
stage5=$9
app_dir=${0%/*}

if [[ "$4" != "''" ]]; then ancestor_genome_gbff="$4.gbff"; ancestor_genome_fna="$4.fna"; else ancestor_genome_gbff=""; ancestor_genome_fna=""; fi
if [[ "$3" != "''"  ]]; then reference_genome_gbff="$3.gbff"; reference_genome_fna="$3.fna"; else reference_genome_gbff=""; reference_genome_fna=""; fi
if [[ "${10}" != "''" ]]; then pseudofinder_path=${10}; else pseudofinder_path="${app_dir}/pseudofinder256"; fi

printf "reference_genome_fna="$reference_genome_fna"
reference_genome_gbff="$reference_genome_gbff"
ancestor_genome_gbff=$ancestor_genome_gbff
ancestor_genome_fna=$ancestor_genome_fna
base_dir=$1
seqid=${1##*/}
id=$5
dnds=$6
scaffold=$7
predictor_type=$8
stage5=$stage5
pseudofinder_path=${10}\n\n"

prokka_pseudos_annotator=${app_dir}/scripts/prokka_anno.pl

mkdir -p "$base_dir/rt_output"

if [[ "$scaffold" = "T" ]] || [[ "$stage5" = "a" ]]; then 
    mkdir -p "$base_dir/rt_output/corrected"

    echo -e "Correcting misassemblies with ragtag\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    ragtag.py correct "$reference_genome_fna" "$base_dir/assembly/skesa.fasta" -w --remove-small -o "$base_dir/rt_output/corrected" \
        1> $base_dir/ragtag_cor.stdout.txt 2> $base_dir/ragtag_cor.stderr.txt
        [ -s $base_dir/ragtag_cor.stderr.txt ] && cat $base_dir/ragtag_cor.stderr.txt || \
        echo -e "Ragtag correction successful\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    grep ">" "$base_dir/rt_output/corrected/ragtag.correct.fasta" | wc -l | echo "$(</dev/stdin) contigs after correction" | tee -a "$base_dir/pseudogenome_processing.txt"
    echo -e "################################################################################################\n"

    echo -e "Scaffolding with ragtag\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    ragtag.py scaffold "$reference_genome_fna" "$base_dir/rt_output/corrected/ragtag.correct.fasta" -w -o "$base_dir/rt_output" \
        1> $base_dir/ragtag_scaf.stdout.txt 2> $base_dir/ragtag_scaf.stderr.txt
        [ -s $base_dir/ragtag_scaf.stderr.txt ] && cat $base_dir/ragtag_scaf.stderr.txt || \
        echo -e "Ragtag scaffolding successful\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    grep ">" "$base_dir/rt_output/ragtag.scaffold.fasta" | wc -l | echo "$(</dev/stdin) scaffolds generated" | tee -a "$base_dir/pseudogenome_processing.txt"
    # cp "$base_dir/rt_output/ragtag.scaffold.fasta" "$base_dir/rt_output/corrected/ragtag.scaffold.fasta" 
    #KP work problems with ragtag solution
    sed -E '/^>/!s/[^ATCGN]*(>|[^ATCGN]).?(>|[^ATCGN])//g' "$base_dir/rt_output/ragtag.scaffold.fasta" > "$base_dir/rt_output/ragtag.scaffold_corrected.fasta"
    #sed -E '/^>/!s/[^ATCGN]*(>|[^ATCGN]).?(>|[^ATCGN])//g' "$base_dir/rt_output/ragtag.scaffold_corrected.fasta" > "$base_dir/rt_output/ragtag.scaffold_corrected.fasta"
    # cp "$base_dir/rt_output/ragtag.scaffold.fasta" "$base_dir/rt_output/ragtag.scaffold_corrected.fasta"
    echo -e "################################################################################################\n"
elif [[ "$scaffold" = "F" ]]; then
    if [ -f "$base_dir/${seqid}.fna" ]; then
        cp -f "$base_dir/${seqid}.fna" "$base_dir/rt_output/ragtag.scaffold_corrected.fasta"
    fi
fi

function genome_assembly_qc () {
    echo -e "Genome assembly QC on final assembly\n"
    start_assembly_qc=`date +%s`
    # ./assembly_qc.sh $base_dir $cpus $assembly_qc_method
    ${base_dir}/scripts/assembly_qc.sh $base_dir $cpus $assembly_qc_method $id
    end_assembly_qc=`date +%s`
    runtime=$((end_assembly_qc-start_assembly_qc))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Genome assembly QC runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Genome assembly QC,$hours:$minutes:$seconds" >> $base_dir/timer.txt
    echo -e "###########################################################################\n"
}

function prokka_annotation () {
    echo -e "Annotating with Prokka\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    mkdir -p "$base_dir/rt_output/prokka_annotation"
    prokka --cpus $cpus --kingdom "Bacteria" --locustag $id --outdir "$base_dir/rt_output/prokka_annotation" --prefix $id \
        --strain "$id" \
        --addgenes --force --compliant "$base_dir/rt_output/ragtag.scaffold_corrected.fasta" \
        2> $base_dir/Prokka.stdout.txt 1> $base_dir/Prokka.stderr.txt
    #--genus "Mycobacterium" --species "ulcerans"
    success_str=$(grep -o "Annotation finished successfully" $base_dir/Prokka.stdout.txt)
    [[ "$success_str" = "Annotation finished successfully" ]] && echo "Prokka $success_str\n" || \
        cat $base_dir/Prokka.stdout.txt | tee -a "$base_dir/pseudogenome_processing.txt"
    grep "Walltime used" $base_dir/rt_output/prokka_annotation/$id.log | tee -a "$base_dir/pseudogenome_processing.txt"
    echo -e "################################################################################################\n"
}

#"Pseudogene annotation with Prokka"
function prokka_pseudo_annotation () {
    echo -e "Writing pseudogenes with Prokka\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    $prokka_pseudos_annotator "$base_dir/rt_output/prokka_annotation/$id.faa"  \
        1> $base_dir/Prokka_pseudo.stdout.txt 2> $base_dir/Prokka_pseudo.stderr.txt
    [ -s $base_dir/Prokka_pseudo.stderr.txt ] && cat $base_dir/Prokka_pseudo.stderr.txt || \
        echo -e "Pseudogenes predicted successfully with Prokka\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    echo -e "################################################################################################\n"
}

#Annotating with DFAST
function DFAST_annotation () {
    conda activate "bactopia_manually"
    echo -e "Annotating with DFAST\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    # conda activate 'DFAST'
    if [[ $ancestor_genome_gbff ]]; then
        dfast --genome "$base_dir/rt_output/ragtag.scaffold_corrected.fasta" --strain "$id" \
            --locus_tag_prefix "${id}_DF" --reference $ancestor_genome_gbff \
            --force -o "$base_dir/rt_output/DFAST_output" \
            1> $base_dir/DFAST.stdout.txt 2> $base_dir/DFAST.stderr.txt
    else
        dfast --genome "$base_dir/rt_output/ragtag.scaffold_corrected.fasta" --strain "$id" \
            --locus_tag_prefix "${id}_DF" \
            --force -o "$base_dir/rt_output/DFAST_output" \
            1> $base_dir/DFAST.stdout.txt 2> $base_dir/DFAST.stderr.txt
    fi
    checking_error=$(grep "Aborting" $base_dir/DFAST.stdout.txt)
    [ -s $base_dir/DFAST.stderr.txt ] || [[ $checking_error ]] && cat $base_dir/DFAST.stderr.txt || \
        echo -e "DFAST annotation successful\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    grep "Total running time" $base_dir/rt_output/DFAST_output/application.log | tee -a "$base_dir/pseudogenome_processing.txt"
    echo -e "################################################################################################\n"
}

#Pseudogene prediction with pseudofinder
function pseudofinder_prediction () {
    echo -e "Predicting pseudogenes with Pseudofinder using Prokka annotation\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    conda activate "pseudofinder"
    start_PF=`date +%s`
    if [[ $ancestor_genome_gbff ]] && [[ $dnds ]]; then
        ${pseudofinder_path}/pseudofinder.py annotate --genome "${seqid}/rt_output/prokka_annotation/${id}.gbk" --outprefix "${id}" \
            --database "${pseudofinder_path}/swissprot.dmnd" --threads $cpus --diamond \
            -ref $ancestor_genome_gbff \
            -op "$seqid/PF_output" -dnds $dnds \
            1> $seqid/pseudofinder.stdout.txt 2> $seqid/pseudofinder.stderr.txt
    else
        ${pseudofinder_path}/pseudofinder.py annotate --genome "${seqid}/rt_output/prokka_annotation/${id}.gbk" --outprefix "${id}" \
            --database "${pseudofinder_path}/swissprot.dmnd" --threads $cpus --diamond \
            -op "$seqid/PF_output" \
            1> $seqid/pseudofinder.stdout.txt 2> $seqid/pseudofinder.stderr.txt
    fi
    [ -s $seqid/pseudofinder.stderr.txt ] && cat $seqid/pseudofinder.stderr.txt || \
        echo -e "Pseudogenes predicted successfully with Pseudofinder\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    mkdir -p "$base_dir/pseudofinder"
    mv -f $base_dir/PF_* $base_dir/pseudofinder
    end_PF=`date +%s`
    PF_runtime=$((end_PF-start_PF))
    hours=$((PF_runtime / 3600)); minutes=$(( (PF_runtime % 3600) / 60 )); seconds=$(( (PF_runtime % 3600) % 60 )) 
    echo "Pseudofinder PF_runtime: $hours:$minutes:$seconds (hh:mm:ss)" | tee -a "$base_dir/pseudogenome_processing.txt"
    echo -e "################################################################################################\n"
}

#Processing of predicted pseudogenes
function pseudogenes_extractor () {
    echo -e "Extracting predicted pseudogenes from $1\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    # conda activate "pseudofinder"
    pseudogene_processor="${app_dir}/scripts/pseudo_proccessor_sub.py"
    DFAST_gbk="$base_dir/rt_output/DFAST_output/genome.gbk"
    prokka_gbk="$base_dir/rt_output/prokka_annotation/$id.gbk"
    PF_gff="$base_dir/pseudofinder/PF_output_pseudos.gff"
    prokka_pseudos_tsv="$base_dir/rt_output/prokka_annotation/$id.faa.prokka_pseudos.tsv"
    outdir="$base_dir/pseudogenes"
    mkdir -p $outdir

    python $pseudogene_processor \
        $DFAST_gbk \
        $prokka_gbk \
        $PF_gff \
        $prokka_pseudos_tsv \
        $outdir $predictor_type\
        1> $base_dir/pseudogene_processor.stdout.txt 2> $base_dir/pseudogene_processor.stderr.txt
    [ -s $base_dir/pseudogene_processor.stderr.txt ] && cat $base_dir/pseudogene_processor.stderr.txt || \
        echo -e "Predicted pseudogenes extracted successfully\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    echo -e "################################################################################################\n"
}

#Removing duplicates with CDHIT
function cleaner_cdhit () {
    conda activate "bactopia_manually"
    echo -e "Removing duplicates with CDHIT\n" | tee -a "$base_dir/pseudogenome_processing.txt"
    mkdir -p "$base_dir/cdhit" 
    cd-hit -i "$base_dir/pseudogenes/combined_pseudos.fasta" -o "$base_dir/cdhit/${id}_cdhit" -c 0.95 -n 5 -g 1 -M 6000 -T $cpus -A 0.95 -s 0.9 \
        1> $base_dir/cd_hit.stdout.txt 2> $base_dir/cd_hit.stderr.txt
    [ -s $base_dir/cd_hit.stderr.txt ] && cat $base_dir/cd_hit.stderr.txt || \
        cp "$base_dir/cdhit/${id}_cdhit" "$base_dir/cdhit/${id}_cdhit.fasta"
    grep ">" "$base_dir/cdhit/${id}_cdhit.fasta" | wc -l | echo -e "$(</dev/stdin) unique pseudogenes found\n" | tee -a "$base_dir/pseudogenome_processing.txt"
}

genome_assembly_qc

if [[ "$stage5" = "b" ]]; then
    prokka_annotation
elif [[ "$stage5" = "c" ]]; then
    prokka_pseudo_annotation
elif [[ "$stage5" = "d" ]]; then
    DFAST_annotation
elif [[ "$stage5" = "e" ]]; then
    pseudofinder_prediction
elif [[ "$stage5" = "f" ]]; then
    pseudogenes_extractor
elif [[ "$stage5" = "g" ]]; then
    cleaner_cdhit
elif [[ "$predictor_type" = "prokka" ]]; then
    prokka_annotation
    prokka_pseudo_annotation
    pseudogenes_extractor $predictor_type
    cleaner_cdhit
elif [[ "$predictor_type" = "PF" ]]; then
    prokka_annotation
    pseudofinder_prediction
    pseudogenes_extractor $predictor_type
    cleaner_cdhit
elif [[ "$predictor_type" = "DFAST" ]]; then
    DFAST_annotation
    pseudogenes_extractor $predictor_type
    cleaner_cdhit
elif [[ "$predictor_type" = "all" ]]; then
    prokka_annotation
    prokka_pseudo_annotation
    pseudofinder_prediction
    DFAST_annotation
    pseudogenes_extractor $predictor_type
    cleaner_cdhit
fi
