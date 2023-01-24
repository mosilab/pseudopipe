#!/bin/bash

############# Global Error handling
# set -o errexit
# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

########## Pipeline start point
start=`date +%s`
# echo $@
# app=${BASH_SOURCE[0]}
# app=$0
app_dir=$( cd -- "$( dirname -- "{BASH_RESOURCE[0]}" )" &> /dev/null && pwd)

conda activate "pseudopipe"

help()
{
    echo "Usage: bash -i main.sh 
               -i | --fastaID (Path with SRR ID)
    Running genome assembly
               -r | --runtype (whether it is single or paired end. Only paired end Illumina reads working now); default=PE 
               -c | --cpus (The number of CPUs to be used in GB)
               -m | --ram (The amount of RAM to be used in GB)
               -a | --assembly_qc (U50 or QUAST); default=U50
    Optional arguments
               -o | --organism_tag (String of letters that define an organism. This is incorporateed in to the organism ID); default=BAC
               -k | --keep_transposease (T/F); default=T
               -s | --scaffold_genome (T/F); ; default=F
               -t | --pseudogene_predictor (PF, prokka, DFAST or all: choice determines te tools used); default=all
               -e | --reference_genome (Reference genome for genome scaffolding)
               -n | --ancestor_genome
               -T | --pseudofinder_path
               -M | --microbeannotatorDB (path to Microbeannotator DB) 
               -d | --dnds (dN/dS cutoff for determining a pseudogene with Pseudofinder); default=20
               -p | --stage_start (The stage of starting the pipeline); default=1
               -E | --stage_end (The stage of ending the pipeline); default=8
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
    -I | --input_file )
      input_file="$2"
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
    -p | --stage_start )
        pipeline_stage_start="$2"
        shift 2
        ;;  
    -E | --stage_end )
        pipeline_stage_end="$2"
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
    -T | --pseudofinder_path )
        predictor_type="$2"
        shift 2
        ;;
    -M | --microbeannotatorDB )
        mannotatorDB="$2"
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
      set -e
      ;;
  esac
done


########## Setting defaults
#Controlling for missing args and correct arg combos
if [[ "${ancestor_genome##*/}" = "''" ]]; then ancestor_genome=""; fi
# if [[ -z "${stage##*/}" ]]; then stage5=1; else stage=$stage; fi
if [[ -z "${stage5##*/}" ]]; then stage5="''"; fi
if [[ -z "${dnds##*/}" ]]; then dnds="''"; fi
if [[ -z "${keep_transposease##*/}" ]]; then keep_transposease="T"; fi
if [[ -z "${scaffold##*/}" ]]; then scaffold="F"; fi
if [[ -z "${assembly_qc_method##*/}" ]]; then assembly_qc_method="U50"; fi
if [[ -z "${dnds##*/}" ]]; then dnds=20; fi
if [[ -z "${organism_tag##*/}" ]]; then organism_tag="BAC"; fi
if [[ -z "${predictor_type##*/}" ]]; then predictor_type="all"; fi
[[ ! -z "$pipeline_stage_start" ]] && pipeline_stage_start=$pipeline_stage_start || pipeline_stage_start=1
[[ ! -z "$pipeline_stage_end" ]] && pipeline_stage_end=$pipeline_stage_end || pipeline_stage_end=8
! [[ -z "$reference_genome" ]] && reference_genome=$reference_genome || reference_genome="''"
! [[ -z "$ancestor_genome" ]]  &&  ancestor_genome=$ancestor_genome || ancestor_genome="''"
! [[ -z "$mannotatorDB" ]]  &&  ancestor_genome=$mannotatorDB || mannotatorDB="''"


scripts_dir="${app_dir}/scripts"

printf "############### Inputs #######################
app_dir=$app_dir
scripts_dir="${app_dir}/scripts"
base_dir=$input_file
runtype=$runtype
cpus=$cpus
ram=$ram
assembly_qc_method=$assembly_qc_method
reference_genome=$reference_genome
ancestor_genome=$ancestor_genome
dnds=$dnds
organism_tag=$organism_tag
keep_transposease="${keep_transposease}"
stage_start="${pipeline_stage_start}"
stage_stop=${pipeline_stage_end}
stage5="${stage_5_specfic_run}"
scaffold="${scaffold}"
predictor_type="${predictor_type}"
pseudofinder_path="${pseudofinder_path}"\n\n"

re='^[0-9]+$'
if [ $cpus ] ;then cpus=$cpus; elif [[ "$cpus" =~ "$re" ]] ; then echo -e 'Enter a number as input for the number of cpus to be used'; set -e; fi
    
function stage_1 ()
{
    if [ 1 -ge $3 ] && [ $4 -le 1 ]; then
        echo -e "###########################################################################\n"
        #echo -e "${base_dir##*/} started at $(( (start / 3600) % 60 )):$(( (start % 3600) / 60 )):$(( (start % 3600) % 60 ))" 
        echo -e "${base_dir##*/} started at $(date | awk '{print $4}')"
        echo -e "Gathering reads\n" 
        start_gather=`date +%s`
        ${2}/gather_samples.sh $1
        end_gather=`date +%s`
        runtime=$((end_gather-start_gather))
        hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
        echo "Gathering reads runtime: $hours:$minutes:$seconds (hh:mm:ss)"
        echo -e ">$base_dir\nGathering reads,$hours:$minutes:$seconds" >> $base_dir/timer.txt
    fi
}

function stage_2 ()
{
    if [ 2 -ge $5 ] && [ $6 -le 2 ]; then
        echo -e "###########################################################################\n"
        echo -e "Reads QC\n"
        start_qc=`date +%s`
        ${4}/qc_reads.sh $1 $2 $3
        end_qc=`date +%s`
        runtime=$((end_qc-start_qc))
        hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
        echo "Reads QC runtime: $hours:$minutes:$seconds (hh:mm:ss)"
        echo "Reads QC,$hours:$minutes:$seconds" >> $base_dir/timer.txt
    fi
}

function stage_3 ()
{
    if [ 3 -ge $7 ] && [ $8 -le 4 ]; then
        echo -e "###########################################################################\n"
        echo -e "Assembling genome\n"
        start_assembly=`date +%s`
        ${6}/genome_assembler.sh $1 $2 $3 $4 $5
        end_assembly=`date +%s`
        runtime=$((end_assembly-start_assembly))
        hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
        echo "Genome assembly runtime: $hours:$minutes:$seconds (hh:mm:ss)"
        echo "Genome assembly,$hours:$minutes:$seconds" >> $base_dir/timer.txt
    fi
}

function stage_4 ()
{
    if [ 4 -ge $6 ] && [ $7 -le 4 ]; then
        echo -e "###########################################################################\n"
        echo -e "Genome assembly QC\n"
        start_assembly_qc=`date +%s`
        # ./assembly_qc_scale.sh $base_dir $cpus $assembly_qc_method $organism_tag
        ${5}/assembly_qc.sh $1 $2 $3 $4
        end_assembly_qc=`date +%s`
        runtime=$((end_assembly_qc-start_assembly_qc))
        hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
        echo "Genome assembly QC runtime: $hours:$minutes:$seconds (hh:mm:ss)"
        echo "Genome assembly QC,$hours:$minutes:$seconds" >> $base_dir/timer.txt
    fi
}

function stage_5 ()
{
    if [ 5 -le ${12} ] && [ ${13} -le 4 ]; then
        echo -e "###########################################################################\n"
        echo -e "Predicting pseudogenome\n"
        start_pseudogene_proc=`date +%s`
        source ${11}/pseudo_proccessor_main.sh $1 $2 $3 $4 $5 $6 \
            $7 $8 $9 ${10}
        end_pseudogene_proc=`date +%s`
        runtime=$((end_pseudogene_proc-start_pseudogene_proc))
        hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
        echo "Pseudogenome processing runtime: $hours:$minutes:$seconds (hh:mm:ss)"
        echo "Pseudogenome processing,$hours:$minutes:$seconds" >> $base_dir/timer.txt
    fi
}

function stage_6 ()
{
    if [ 6 -le $7 ] && [ $8 -le 6 ]; then
        conda activate 'pseudopipe'
        echo -e "###########################################################################\n"
        echo -e "Cleaning and writing pseudogenome\n"
        start_pseudogene_cleaning=`date +%s`
        bash -i ${1}/cleaning_pseudogenome.sh $2 $3 $4 $5 ${6}
        end_pseudogene_cleaning=`date +%s`
        runtime=$((end_pseudogene_cleaning-start_pseudogene_cleaning))
        hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
        echo "Pseudogenome cleaning runtime: $hours:$minutes:$seconds (hh:mm:ss)"
        echo "Pseudogenome cleaning,$hours:$minutes:$seconds" >> $base_dir/timer.txt
    fi
}

function stage_7 ()
{
    if [ 7 -le $5 ] && [ $6 -le 7 ]; then
        conda activate 'cog_classifier'
        echo -e "###########################################################################\n"
        echo -e "Running COG analysis\n"
        start_COG_analysis=`date +%s`
        source ${1}/cog.sh $2 $3 $4
        end_COG_analysis=`date +%s`
        runtime=$((end_COG_analysis-start_COG_analysis))
        hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
        echo "COG analysis runtime: $hours:$minutes:$seconds (hh:mm:ss)"
        echo "COG analysis,$hours:$minutes:$seconds" >> $base_dir/timer.txt
    fi
}

function stage_8 ()
{
    if [ 8 -le $6 ] && [ $7 -le 8 ]; then
        conda activate 'pseudopipe'
        echo -e "###########################################################################\n"
        echo -e "Running KEGG analysis\n"
        start_KEGG_analysis=`date +%s`
        source ${1}/kegg.sh $2 $3 $4 $5
        end_KEGG_analysis=`date +%s`
        runtime=$((end_KEGG_analysis-start_KEGG_analysis))
        hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
        echo "KEGG analysis runtime: $hours:$minutes:$seconds (hh:mm:ss)"
        echo "KEG analysis,$hours:$minutes:$seconds" >> $base_dir/timer.txt
    fi
}

function stage_9 ()
{
        # conda activate 'pseudopipe'
        echo -e "###########################################################################\n"
        echo -e "Running secondary analyses\n"
        start_sec_analysis=`date +%s`
        # echo $1 $5
        source ${5}/secondary_analyses/secondary.sh $1 $2 $3 $4 $6
        end_sec_analysis=`date +%s`
        runtime=$((end_sec_analysis-start_sec_analysis))
        hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
        echo "Secondary analyses runtime: $hours:$minutes:$seconds (hh:mm:ss)"
        echo "Secondary analyses,$hours:$minutes:$seconds\n" >> $base_dir/timer.txt
}

if [ $input_file ] ;then
    srr_list=()
    while IFS='\n' read -r line || [ "$line" ]; do
        working_dir=${input_file%/*}
        base_dir="${working_dir}/${line}"
        seqid=$line
        mkdir -p $base_dir
        # echo -e "$line\n${base_dir}\n${working_dir}\n${scripts_dir}"
        if ! [[ -z "${organism_tag##*/}" ]]; then organism_tagged="${organism_tag}${seqid:3}"; touch "$base_dir/orgtag_$organism_tag"; \
            else organism_tagged=${base_dir##*/}; #touch "$base_dir/orgtag_$organism_tagged"; 
            fi

        echo -e "\nStarting $seqid"
        # base_dir=$1
        # runtype=$2 #either blank for PE or any letter for SE; preferably SE
        [ "$runtype" = "SE" ] && echo "Using single-end reads" || echo "Using paired-end reads"
        # echo $pipeline_stage_start 

        stage_1 $base_dir ${scripts_dir} $pipeline_stage_start $pipeline_stage_end
        
        stage_2 $base_dir $runtype $cpus ${scripts_dir} $pipeline_stage_start $pipeline_stage_end
        stage_3 $base_dir $runtype $cpus $ram $organism_tag ${scripts_dir} $pipeline_stage_start $pipeline_stage_end
        stage_4 $base_dir $cpus $assembly_qc_method $organism_tag ${scripts_dir} $pipeline_stage_start $pipeline_stage_end
        stage_5 $base_dir $cpus $reference_genome $ancestor_genome $organism_tag $dnds \
            $scaffold $predictor_type $stage5 $pseudofinder_path ${scripts_dir} $pipeline_stage_start $pipeline_stage_end
        stage_6 ${scripts_dir} $base_dir $organism_tagged $keep_transposease $predictor_type ${app_dir} $pipeline_stage_start $pipeline_stage_end
        stage_7 ${scripts_dir} $base_dir $organism_tagged $cpus $pipeline_stage_start $pipeline_stage_end
        stage_8 ${scripts_dir} $base_dir $organism_tagged $cpus $mannotatorDB $pipeline_stage_start $pipeline_stage_end

        srr_list+=( "$line" )
        # echo -e "${srr_list[@]}"
        organism_tag=$( basename "$base_dir/orgtag_$organism_tag" | cut -d_ -f2 )
    done < $input_file #input should have LF line endings(thereofre input file has to be in unix) or $(perl -pe 's/\r\n|\n|\r/\n/g' "$INPUTFILE")
    # echo -e "\n${srr_list[@]}"

    stage_9 ${working_dir} $organism_tag $cpus $app_dir ${scripts_dir} $srr_list

    rm -f "$base_dir/orgtag_$organism_tag" "$base_dir/orgtag_$organism_tagged"
    mkdir -p "$base_dir/log_files"
    mv -t $base_dir/log_files $base_dir/*out.txt $base_dir/*err.txt $base_dir/*processing.txt $base_dir/timer.txt 2>/dev/null || :
else
    base_dir=$absolute_path
    seqid=${absolute_path##*/}
    mkdir -p $base_dir
    echo -e "Starting $seqid\n"
    # base_dir=$1
    # runtype=$2 #either blank for PE or any letter for SE; preferably SE
    [ "$runtype" = "SE" ] && echo "Using single-end reads" || echo "Using paired-end reads"
    if ! [[ -z "${organism_tag##*/}" ]]; then organism_tagged="${organism_tag}${seqid:3}"; else organism_tagged=${absolute_path##*/}; fi

    stage_1 $base_dir ${scripts_dir}
    stage_2 $base_dir $runtype $cpus ${scripts_dir}
    stage_3 $base_dir $runtype $cpus $ram $organism_tag ${scripts_dir}
    stage_4 $base_dir $cpus $assembly_qc_method $organism_tag ${scripts_dir}
    stage_5 $base_dir $cpus $reference_genome $ancestor_genome $organism_tag $dnds \
        $scaffold $predictor_type $stage5 $pseudofinder_path ${scripts_dir}
    stage_6 ${scripts_dir} $base_dir $organism_tagged $keep_transposease $predictor_type ${app_dir}
    stage_7 ${scripts_dir} $base_dir $organism_tagged $cpus
    stage_8 ${scripts_dir} $base_dir $organism_tagged $cpus $mannotatorDB

    mkdir -p "$base_dir/log_files"
    mv -t $base_dir/log_files $base_dir/*out.txt $base_dir/*err.txt $base_dir/*processing.txt $base_dir/timer.txt 
fi

# find $base_dir -maxdepth 1 -type f -name "*.txt" -o -name "*.log" -o -name "total_contigs*" -o -name "*gz" | xargs -I{} rm "{}" 
# find $base_dir -maxdepth 1 -name "shovill" | xargs -I{} rm -r "{}" 
# rm -rf ./test.msh ./*.bt2 ./total_contigs*

# mv -t $base_dir ./test.msh ./*.bt2 ./total_contigs*
# rm -rf $base_dir/*.bt2 $base_dir/*.fastq.gz


end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
echo "Pipeline runtime: $hours:$minutes:$seconds (hh:mm:ss)" | tee -a $base_dir/timer.txt timer.txt
echo -e "################################################################################################\n"