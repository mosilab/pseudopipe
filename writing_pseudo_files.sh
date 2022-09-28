#!/bin/bash

seqid=$1
id="$2"
predictor_type=$3
#id=$seqid
root_base_dir=${seqid%/*}
# echo "$root_base_dir\n"
mkdir -p "$root_base_dir/pan_files"


echo -e "\nWriting pseudogenome gff, faa and fna"
mkdir -p "$root_base_dir/pan_files/pan-pseudogenome/gffs"

# writing gff header
LC_ALL=C grep "^##" $seqid/rt_output/prokka_annotation/$id.gff | head -n -1 > $root_base_dir/pan_files/pan-pseudogenome/gffs/${id}.gff 

function DFAST_gff_writer () {
    ids=$( grep "DF" ${seqid}/combined_pseudogenome_unique.txt | cut -d' ' -f1 | tr -d '>' )

    # writing records from original annotation files
    #TO-DO Find a better way to eliminate all duplicate ids that affect roary analysis
    LC_ALL=C egrep -h "$ids" $seqid/rt_output/DFAST_output/genome.gff >> $root_base_dir/pan_files/pan-pseudogenome/gffs/${id}.gff

    # writing genome sequence at the end of gff file
    sed -n '/^##FASTA/,$p' $seqid/rt_output/DFAST_output/genome.gff >> $root_base_dir/pan_files/pan-pseudogenome/gffs/${id}.gff
}

function prokka_gff_writer () {
    ids1=$( grep -v "DF" ${seqid}/combined_pseudogenome_unique.txt | cut -d' ' -f1 | sed 's/PF_//' | tr -d '>' | sort )
    #grep -v "DF" ${seqid}/combined_pseudogenome_unique.txt | cut -d' ' -f1 | sed 's/PF_//' | tr -d '>' | sort | sed 's/>//' 
    #grep -v "DF" ${seqid}/combined_pseudogenome_unique.txt | cut -d' ' -f1 | sed 's/PF_//' | tr -d '>' #| xargs printf "%s|" #| head -c -1 

    # writing records from original annotation files
    #TO-DO Find a better way to eliminate all duplicate ids that affect roary analysis
    LC_ALL=C egrep -h "$ids1" $seqid/rt_output/prokka_annotation/$id.gff | sed '/#sequence/d' | sed '/>/d' >> $root_base_dir/pan_files/pan-pseudogenome/gffs/${id}.gff

    # writing genome sequence at the end of gff file
    sed -e '1,/^##FASTA/ d' $seqid/rt_output/prokka_annotation/$id.gff >> $root_base_dir/pan_files/pan-pseudogenome/gffs/${id}.gff
}

# linearize fastas (prokka and DFAST)
[[ -f $seqid/rt_output/prokka_annotation/$id.linear.faa && -f $seqid/rt_output/DFAST_output/$id.linear.cds.fna ]] && echo -e "\nLinear fasta files present\n" || \
    echo -e "\nLinearizing fasta files\n"
    for file_type in ffn faa; do 
        awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' \
        < $seqid/rt_output/prokka_annotation/$id.${file_type} > $seqid/rt_output/prokka_annotation/$id.linear.${file_type}
    done
    for file_type in cds.fna protein.faa; do 
        awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' \
        < $seqid/rt_output/DFAST_output/${file_type} > $seqid/rt_output/DFAST_output/$id.linear.${file_type}
    done

for genome in pseudogenome functional_genome; do
    # writing fastas
    mkdir -p "$root_base_dir/pan_files/$genome/faas"
    mkdir -p "$root_base_dir/pan_files/$genome/fnas"

    if [[ "$predictor_type" = "prokka" ]] || [[ "$predictor_type" = "PF" ]]; then
        ids2=$(grep -v 'DF' "${seqid}/combined_pseudogenome_unique.txt"  | cut -d' ' -f1 | sed 's/PF_//')
    elif [[ "$predictor_type" = "DFAST" ]]; then
        ids3=$(grep 'DF' "${seqid}/combined_pseudogenome_unique.txt"  | cut -d' ' -f1)
    elif [[ "$predictor_type" = "all" ]]; then
        ids2=$(grep -v 'DF' "${seqid}/combined_pseudogenome_unique.txt"  | cut -d' ' -f1 | sed 's/PF_//')
        ids3=$(grep 'DF' "${seqid}/combined_pseudogenome_unique.txt"  | cut -d' ' -f1)
    fi

    if [[ "$genome" = "pseudogenome" ]]; then

        if [[ "$predictor_type" = "prokka" ]] || [[ "$predictor_type" = "PF" ]]; then
            idstwo_count=$(grep -c ">" <<< $ids2)
            pseuds=$idstwo_count
            # faa
            grep -A 1 "$ids2" "$seqid/rt_output/prokka_annotation/$id.linear.faa" --no-group-separator > "$root_base_dir/pan_files/$genome/faas/${id}.faa"
            # fna
            grep -A 1 "$ids2"  "$seqid/rt_output/prokka_annotation/$id.linear.ffn" --no-group-separator > "$root_base_dir/pan_files/$genome/fnas/${id}.fna"
            
        elif [[ "$predictor_type" = "DFAST" ]]; then
            idsthree_count=$(grep -c ">" <<< $ids3)
            pseuds=$idsthree_count
            # faa
            grep -A 1 "$ids3" "$seqid/rt_output/DFAST_output/$id.linear.protein.faa" --no-group-separator >> "$root_base_dir/pan_files/$genome/faas/${id}.faa"
            # fna
            grep -A 1 "$ids3" "$seqid/rt_output/DFAST_output/$id.linear.cds.fna" --no-group-separator >> "$root_base_dir/pan_files/$genome/fnas/${id}.fna"
        
        elif [[ "$predictor_type" = "all" ]]; then
            idstwo_count=$(grep -c ">" <<< $ids2)
            idsthree_count=$(grep -c ">" <<< $ids3)
            pseuds=$(( idstwo_count + idsthree_count ))        
            # faa
            grep -A 1 "$ids2" "$seqid/rt_output/prokka_annotation/$id.linear.faa" --no-group-separator > "$root_base_dir/pan_files/$genome/faas/${id}.faa"
            grep -A 1 "$ids3" "$seqid/rt_output/DFAST_output/$id.linear.protein.faa" --no-group-separator >> "$root_base_dir/pan_files/$genome/faas/${id}.faa"
            # fna
            grep -A 1 "$ids2"  "$seqid/rt_output/prokka_annotation/$id.linear.ffn" --no-group-separator > "$root_base_dir/pan_files/$genome/fnas/${id}.fna"
            grep -A 1 "$ids3" "$seqid/rt_output/DFAST_output/$id.linear.cds.fna" --no-group-separator >> "$root_base_dir/pan_files/$genome/fnas/${id}.fna"
        fi
    
    elif [[ "$genome" = "functional_genome" ]]; then
        echo -e "Writing functional genome files"
        # all=$(grep ">" "$seqid/rt_output/prokka_annotation/$id.linear.faa")
        # func=$(grep -v "$ids2" <<< $all | cut -d' ' -f1 | tr -d '>') 
        all=$(grep -c ">" "$seqid/rt_output/prokka_annotation/$id.linear.faa")
        func=$(( all - pseuds ))
        echo -e "Whole genome => $all, Functional genome => $func"
        # grep -c "KO" <<< $ids2
        # grep -c "KO" <<< $func
        # # faa
        # # split -l 50 $func $genome/faas/pattern-file.split.
        # # while read line ; do grep -A 1 "$line" $seqid/rt_output/prokka_annotation/$id.linear.faa --no-group-separator > "$genome/faas/${id}.faa" ; done <<< $func

        # # for CHUNK in pattern-file.split.* ; do
        # #         grep -A 1 -f "$CHUNK" $seqid/rt_output/prokka_annotation/$id.linear.faa --no-group-separator > "$genome/faas/${id}.faa"
        # # done
        # # rm pattern-file.split.*

        # grep -c -A 1 "$func" "$seqid/rt_output/prokka_annotation/$id.linear.faa" --no-group-separator > "$genome/faas/${id}.faa"
        # # grep -A 1 -v "$ids3" "$seqid/rt_output/DFAST_output/$id.linear.protein.faa" --no-group-separator >> "$genome/faas/${id}.faa"
        # # fna
        # grep -c -A 1 "$func"  "$seqid/rt_output/prokka_annotation/$id.linear.ffn" --no-group-separator > "$genome/fnas/${id}.fna"
        # # grep -A 1 -v "$ids3" "$seqid/rt_output/DFAST_output/$id.linear.cds.fna" --no-group-separator >> "$genome/fnas/${id}.fna"

        # echo $CONDA_DEFAULT_ENV
        conda activate "bactopia_manually"
        cd-hit-2d -i $root_base_dir/pan_files/pseudogenome/faas/${id}.faa -i2 ${seqid}/rt_output/prokka_annotation/${id}.linear.faa \
            -o $root_base_dir/pan_files/$genome/faas/${id}_cdhit -c 0.99 -n 5 -g 1 -M 2000 -T 4 > $root_base_dir/pan_files/$genome/faas/func_genome_analysis.log
        mv $root_base_dir/pan_files/$genome/faas/${id}_cdhit $root_base_dir/pan_files/$genome/faas/${id}.faa

        cd-hit-est-2d -i $root_base_dir/pan_files/pseudogenome/fnas/${id}.fna -i2 ${seqid}/rt_output/prokka_annotation/${id}.linear.ffn \
            -o $root_base_dir/pan_files/$genome/fnas/${id}_cdhit -c 0.99 -n 8 -g 1 -M 2000 -T 4 > $root_base_dir/pan_files/$genome/fnas/func_genome_analysis.log
        mv $root_base_dir/pan_files/$genome/fnas/${id}_cdhit $root_base_dir/pan_files/$genome/fnas/${id}.fna 
        rm $root_base_dir/pan_files/$genome/faas/${id}_cdhit.clstr $root_base_dir/pan_files/$genome/fnas/${id}_cdhit.clstr

        faa=$(grep -c ">" $root_base_dir/pan_files/$genome/faas/${id}.faa)
        fna=$(grep -c ">" $root_base_dir/pan_files/$genome/fnas/${id}.fna)
        echo -e "$faa and $fna prokka-based functional proteins and genes respectively\n"
        if [ $faa -ne $fna ]; then
            if [  $faa -eq $func ]; then
                # takes its ids and get them from fna
                faa_all=$(grep ">" "$root_base_dir/pan_files/$genome/faas/${id}.faa" | cut -d' ' -f1 | tr -d '>')
                grep -A 1 "$faa_all" $root_base_dir/pan_files/$genome/fnas/${id}.fna > $root_base_dir/pan_files/$genome/fnas/${id}_.fna 
                mv $root_base_dir/pan_files/$genome/fnas/${id}_.fna $root_base_dir/pan_files/$genome/fnas/${id}.fna 
                grep -c ">" "$root_base_dir/pan_files/$genome/fnas/${id}.fna"
            elif [ $fna -eq $func ]; then
                # takes its ids and get them from faa
                fnaa_all=$(grep ">" "$root_base_dir/pan_files/$genome/fnas/${id}.fna" | cut -d' ' -f1 | tr -d '>')
                grep -A 1 "$faa_all" $root_base_dir/pan_files/$genome/faas/${id}.faa > $root_base_dir/pan_files/$genome/faas/${id}_.faa
                mv $root_base_dir/pan_files/$genome/faas/${id}_.faa $root_base_dir/pan_files/$genome/faas/${id}.faa
                grep -c ">" "$root_base_dir/pan_files/$genome/faas/${id}.faa"
            fi
        fi
    fi
done

# copying whole genome files
echo -e "Copying whole genome files\n"
mkdir -p $root_base_dir/pan_files/whole_genome/faas $root_base_dir/pan_files/whole_genome/fnas \
    $root_base_dir/pan_files/whole_genome/gffs 
#id="${seqid##*/}"; id="$id" 
cp $seqid/rt_output/prokka_annotation/$id.linear.faa $root_base_dir/pan_files/whole_genome/faas/$id.faa 
cp $seqid/rt_output/prokka_annotation/$id.linear.ffn $root_base_dir/pan_files/whole_genome/fnas/$id.fna 
cp $seqid/rt_output/prokka_annotation/$id.gff $root_base_dir/pan_files/whole_genome/gffs/$id.gff 