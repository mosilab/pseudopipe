#!/bin/bash
start=`date +%s`
# echo $@
# app=${BASH_SOURCE[0]}
app=$0
app_dir="${app%/*}"

conda activate "bactopia_manually"

help()
{
    echo "Usage: bash -i main.sh 
               -i | --fastaID
               -o | --organism_tag
               -k | --keep_transposease
               -s | --scaffold_genome
               -t | --pseudogene_predictor
    Running genome assembly
               -r | --runtype
               -c | --cpus
               -m | --ram
               -a | --assembly_qc
    Optional arguments
               -e | --reference_genome
               -n | --ancestor_genome
               -d | --dnds
               -p | --pipeline_stage
               -u | --stage_5_specfic_run
               -v | --verbose
               -h | --help"

    set -e
}

# SHORT=i:,r:,c:,m:,a:,e:,n::,d:,o:,k:,p:,s:,t:,h
# LONG=fastaID:,runtype:,cpus:,ram:,assembly_qc:,reference_genome:,ancestor_genome:,\
#     dnds:,organism_tag:,keep_transposease:,pipeline_stage:,scaffold_genome:,pseudogene_predictor:,help
# OPTS=$(getopt -a -n main.sh --options $SHORT --longoptions $LONG -- "$@")

VALID_ARGUMENTS=$# # Returns the count of arguments that are in short or long options

if [ "$VALID_ARGUMENTS" -eq 0 ]; then
  help
fi

# eval set -- "$OPTS"

###TO-DO add if statements to link inputs and control for wrong combinations of inputs
###TO-DO add overwriting to args parsable
while [[ "$#" -gt 0 ]];
do
  case "$1" in
    -i | --fastaID )
      absolute_path="$2"
      shift 2
      ;;
    -r | --runtype )
      runtype="$2"
      shift 2
      ;;
    -c | --cpus )
        cpus="$2"
        shift 2
        ;;
    -m | --ram )
        ram="$2"
        shift 2
        ;;
    -a | --assembly_qc )
        assembly_qc_method="$2"
        shift 2
        ;;
    -e | --reference_genome )
        reference_genome="$2"
        shift 2
        ;;
    -n | --ancestor_genome )
        ancestor_genome="$2"
        shift 2
        ;;
    -d | --dnds )
        dnds="$2"
        shift 2
        ;;    
    -o | --organism_tag  )
        organism_tag="$2"
        shift 2
        ;;    
    -k | --keep_transposease )
        keep_transposease="$2"
        shift 2
        ;;    
    -p | --pipeline_stage )
        stage="$2"
        shift 2
        ;;  
    -u | --stage_5_specfic_run )
        stage5="$2"
        shift 2
        ;;    
    -s | --scaffold_genome )
        scaffold="$2"
        shift 2
        ;;   
    -t | --pseudogene_predictor )
        predictor_type="$2"
        shift 2
        ;;   
    -v | --verbose )
        verbose="$2"
        shift 2
        ;;
    -h | --help )
        help
        break
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      help
      break
      ;;
  esac
done

# printf "app_dir=${app_dir}
# base_dir=$absolute_path
# runtype=$runtype
# cpus=$cpus
# ram=$ram
# assembly_qc_method=$assembly_qc_method
# reference_genome=$reference_genome
# ancestor_genome=$ancestor_genome
# dnds=$dnds
# organism_tag=$organism_tag
# keep_transposease="${keep_transposease}"
# stage="${stage}"
# stage5="${stage5}"
# scaffold="${scaffold}"
# predictor_type="${predictor_type}"\n\n"

base_dir=$absolute_path
seqid=${absolute_path##*/}
mkdir -p $base_dir

#Controlling for missing args and correct arg combos
if [[ "${ancestor_genome##*/}" = "''" ]]; then ancestor_genome=""; fi
if [[ -z "${stage##*/}" ]]; then stage5=1; else stage=$stage; fi
if [[ -z "${stage5##*/}" ]]; then stage5="''"; fi
if [[ -z "${dnds##*/}" ]]; then dnds="''"; fi
if ! [[ -z "${organism_tag##*/}" ]]; then organism_tag="${organism_tag}${seqid:3}"; else organism_tag=${absolute_path##*/}; fi

echo -e "Starting $seqid\n"
# base_dir=$1
# runtype=$2 #either blank for PE or any letter for SE; preferably SE
[ "$runtype" = "SE" ] && echo "Using single-end reads" || echo "Using paired-end reads"
# cpus=$3
# ram=$4
# assembly_qc_method=$5
! [[ -z "$reference_genome" ]] && reference_genome=$reference_genome || reference_genome="''"
! [[ -z "$ancestor_genome" ]]  &&  ancestor_genome=$ancestor_genome || ancestor_genome="''"
# dnds=$8
# organism_tag=$9
# keep_transposease="${10}"
# stage="${11}"
# scaffold="${12}"

printf "############### Inputs #######################
app_dir=$app_dir
base_dir=$base_dir
runtype=$runtype
cpus=$cpus
ram=$ram
assembly_qc_method=$assembly_qc_method
reference_genome=$reference_genome
ancestor_genome=$ancestor_genome
dnds=$dnds
organism_tag=$organism_tag
keep_transposease="${keep_transposease}"
stage="${stage}"
stage5="${stage5}"
scaffold="${scaffold}"
predictor_type="${predictor_type}"\n\n"

re='^[0-9]+$'
if [ $cpus ] ;then cpus=$cpus; elif [[ "$cpus" =~ "$re" ]] ; then echo -e 'Enter a number as input for the number of cpus to be used'; fi

if [ $stage -le 1 ]; then
    echo -e "###########################################################################\n"
    #echo -e "${base_dir##*/} started at $(( (start / 3600) % 60 )):$(( (start % 3600) / 60 )):$(( (start % 3600) % 60 ))" 
    echo -e "${base_dir##*/} started at $(date | awk '{print $4}')"
    echo -e "Gathering reads\n" 
    start_gather=`date +%s`
    ${app_dir}/gather_samples.sh $base_dir
    end_gather=`date +%s`
    runtime=$((end_gather-start_gather))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Gathering reads runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo -e ">$base_dir\nGathering reads,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 2 ]; then
    echo -e "###########################################################################\n"
    echo -e "Reads QC\n"
    start_qc=`date +%s`
    ${app_dir}/qc_reads.sh $base_dir $runtype $cpus
    end_qc=`date +%s`
    runtime=$((end_qc-start_qc))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Reads QC runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Reads QC,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 3 ]; then
    echo -e "###########################################################################\n"
    echo -e "Assembling genome\n"
    start_assembly=`date +%s`
    ${app_dir}/genome_assembler.sh $base_dir $runtype $cpus $ram $organism_tag
    end_assembly=`date +%s`
    runtime=$((end_assembly-start_assembly))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Genome assembly runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Genome assembly,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 4 ]; then
    echo -e "###########################################################################\n"
    echo -e "Genome assembly QC\n"
    start_assembly_qc=`date +%s`
    # ./assembly_qc_scale.sh $base_dir $cpus $assembly_qc_method $organism_tag
    ${app_dir}/assembly_qc.sh $base_dir $cpus $assembly_qc_method $organism_tag
    end_assembly_qc=`date +%s`
    runtime=$((end_assembly_qc-start_assembly_qc))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Genome assembly QC runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Genome assembly QC,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 5 ]; then
    echo -e "###########################################################################\n"
    echo -e "Predicting pseudogenome\n"
    start_pseudogene_proc=`date +%s`
    source ${app_dir}/pseudo_proccessor_main.sh $base_dir $cpus $reference_genome $ancestor_genome $organism_tag $dnds \
        $scaffold $predictor_type $stage5
    end_pseudogene_proc=`date +%s`
    runtime=$((end_pseudogene_proc-start_pseudogene_proc))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Pseudogenome processing runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Pseudogenome processing,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 6 ]; then
    conda activate 'bactopia_manually'
    echo -e "###########################################################################\n"
    echo -e "Cleaning and writing pseudogenome\n"
    start_pseudogene_cleaning=`date +%s`
    bash -i ${app_dir}/cleaning_pseudogenome.sh $base_dir $organism_tag $keep_transposease $predictor_type ${app_dir}
    end_pseudogene_cleaning=`date +%s`
    runtime=$((end_pseudogene_cleaning-start_pseudogene_cleaning))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Pseudogenome cleaning runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Pseudogenome cleaning,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 7 ]; then
    conda activate 'cog_classifier'
    echo -e "###########################################################################\n"
    echo -e "Running COG analysis\n"
    start_COG_analysis=`date +%s`
    source ${app_dir}/cog_classifier.sh $base_dir $organism_tag $cpus
    end_COG_analysis=`date +%s`
    runtime=$((end_COG_analysis-start_COG_analysis))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "COG analysis runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "COG analysis,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

# find $base_dir -maxdepth 1 -type f -name "*.txt" -o -name "*.log" -o -name "total_contigs*" -o -name "*gz" | xargs -I{} rm "{}" 
# find $base_dir -maxdepth 1 -name "shovill" | xargs -I{} rm -r "{}" 
# rm -rf ./test.msh ./*.bt2 ./total_contigs*

# mv -t $base_dir ./test.msh ./*.bt2 ./total_contigs*
# rm -rf $base_dir/*.bt2 $base_dir/*.fastq.gz
mkdir -p "$base_dir/log_files"
mv -t $base_dir/log_files $base_dir/*out.txt $base_dir/*err.txt $base_dir/*processing.txt $base_dir/timer.txt 

end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
echo "Pipeline runtime: $hours:$minutes:$seconds (hh:mm:ss)"
echo "$base_dir,$hours:$minutes:$seconds" | tee -a $base_dir/timer.txt timer.txt
echo -e "################################################################################################\n"