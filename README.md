# **Pseudopipe**


## Overview

Pseudopipe is a commandline tool for bacteria pseudogene predictor and pan-pseudogenome analyzer. It wraps three different pseudogene prediction tool including Prokka, DFAST and Pseudofinder, and also provides a combined version. This tool was built to help with understanding the evolution of bacteria from the perspective on the accumulation of pseudogenes. It further provides COG analysis of pseudogenes and other quantitative and qualitative analysis of bacterial pseudogenes to provide insight into the variations in pseudogenes that may indicate evolutionary trends in bacteria. 

Pseudopipe takes Illumina reads or an assemble whole genome and assembles and annotates the bacterial genome from stage 1 to stage 4 based on the input parameters. Pseudogene prediction can be done by one tool or all tools to provide a comprehensive list of predicted pseudogenes in stage 5. Downstream analysis are then performed to provide insights into the pan-pseudogenome of different bacterial strains in stage 7 

![Pipeline flow diagram](https://github.com/eddykay310/demeter/blob/main/img/flow.png?raw=true)

## Installation

## Clone repo
```
git clone https://github.com/mosilab/pseudopipe.git
```
## Create environments
```
conda create --name bactopia_manually --file <path>/pseudopipe/requirements/base_requirement.txt
```
```
conda create --name pseudofinder --file <path>/pseudopipe/requirements/PF_requirement.txt
```
```
conda create --name cog_classifier --file <path>/pseudopipe/requirements/cog_requirement.txt
```

## Usage

**Basic**

Using paired-end Illumina reads

```
bash -i path/pseudopipe/main.sh -i <path>/[SRA gzipped fasta] -r PE -c 4 -m 8 -a U50 -t all
```

Using assembled whole genome

```
bash -i path/pseudopipe/main.sh -i <path>/[SRA gzipped fasta] -r PE -c 4 -m 8 -a U50 -p 5 -t all
```

**Options**

```
[ -fd | --fasta ID ]	Path to input file (sra for Illumina paired-end reads or assembled whole genome fasta file)
[ -rt | --runtype ]	Paired-end (for now with single-end reads in the works)
[ -c | --cpus ]	Number of cpus to be used
[ -rm | --ram ]	Amount of RAM to be usd
[ -aq | --assembly_qc ]	Type of assembly QC to use
[ -ref | --reference_genome ]	Reference genome to be used for correcting missambles and scaffolding
[ -anc | --ancestor_genome ]	Ancestor genome for calculating selection pressure (dn/ds)
[ -d | --dnds ]	dn/ds threshold to be used for determining a pseudogene
[ -ot | --organism_tag ]	Tag used for assigning unique labels during analysis (eg. gene ids during annotation)
[ -kt | --keep_transposease ]	Keep transposeases as pseudogenes 
[ -ps | --pipeline_stage ]	Stage to start pipeline
[ -s | --scaffold_genome ]	Scaffold genome (need a reference genome)
[ -pt | --pseudogene_predictor ]	Pseudogene tool to use [Prokka (prokka), DFAST (DFAST), Psedofinder (PF), All (all) ]
[ -h | --help ]	Help
```

## Feedback

Please file questions, bugs or ideas to the [Issue Tracker](https://github.com/mosilab/pseudopipe/issues)

## License

[BSD3-3-Clause license](https://github.com/eddykay310/demeter/blob/v0.1.0/LICENSE)

## Citation

Not published yet.

## Author

- Edwin Kyei-Baffour
- Web: https://github.com/eddykay310
