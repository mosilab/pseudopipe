# **Pseudopipe**


## Table of Contents

------

- [Table of Contents](#table-of-contents)
- [Overview](#overview)
- [Conda Installation](#conda-installation)
- [Installation](#installation)
  * [microbeannotator evironment](#microbeannotator-evironment)
- [Usage](#usage)
- [Feedback](#feedback)
- [License](#license)
- [Citation](#citation)
- [Author](#author)


## Overview

------

Pseudopipe is a commandline tool for bacteria pseudogene predictor and pan-pseudogenome analyzer. It wraps three different pseudogene prediction tool including Prokka, DFAST and Pseudofinder. This tool was built to help with understanding the evolution of bacteria from the perspective on the accumulation of pseudogenes. It further provides COG analysis of pseudogenes and other quantitative and qualitative analysis of bacterial pseudogenes to provide insight into the variations in pseudogenes that may indicate evolutionary trends in bacteria. 

Pseudopipe takes Illumina reads or an assemble whole genome and assembles and annotates the bacterial genome from stage 1 to stage 4 based on the input parameters. Pseudogene prediction can be done by one tool or all tools to provide a comprehensive list of predicted pseudogenes in stage 5. Downstream analysis are then performed to provide insights into the pan-pseudogenome of different bacterial strains in stages six, seven and eight.

![image-20220712133902287](https://github.com/mosilab/pseudopipe/blob/main/img/flow.png)


## Conda Installation
Pseudopipe uses Conda. You can install if you do not have it.

Download and install Miniconda by following the instructions at https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html.

Add conda channels:
```
conda config --add channels conda-forge
conda config --add channels bioconda
```

## Installation

```
## Clone repo
git clone https://github.com/mosilab/pseudopipe.git

cd pseudopipe

## Create environments
conda create --name pseudopipe --file requirements/base_requirement.txt
pip install -r requirements/pip_requirements.txt
conda create --name pseudofinder --file requirements/PF_requirement.txt
conda create --name cog_classifier --file requirements/cog_requirement.txt
```
### microbeannotator evironment
Check https://github.com/cruizperez/MicrobeAnnotator for microbeannotator installation. Installation should be done in the pseudopipe environment


## Usage

**Basic**

Using paired-end Illumina reads

```
bash -i path/demeter/main.sh -i path/[SRA gzipped fasta (SRR****.fastq.gz); don't add the .fastq.gz to the input name] -r PE -c 4 -m 8 -a U50 -t all
```

Using assembled whole genome

```
bash -i path/demeter/main.sh -i path/[assembled fasta file] -r PE -c 4 -m 8 -p 5 -t all
```
**NB** Pipeline must be started from stage five when using assembled genome. Scaffolding has to be skipped. Also, the genome must be put in a folder as shown below.
```
....Main folder (folder ID name which is passed as input to the pipeline, eg ABN21)
........ABN21.fna (this is the assembled genome file)
```

If scaffolding is used in the pipeline, then the assembled genome has to placed in the directory below.
```
....Main folder (folder ID name which is passed as input to the pipeline)
........assembly
..............skesa.fasta (this is the assembled genome file)
```


**Options**

```
Required argument
            -i | --fastaID (Path with SRR ID) or the input file option (-I)
            -I | --input_file (List of reads or genome Ids
            -c | --cpus (The number of CPUs to be used in GB)
            -m | --ram (The amount of RAM to be used in GB)
Running genome assembly
            -r | --runtype (whether it is single or paired end. Only paired end Illumina reads working now); default=PE 
            -a | --assembly_qc (U50 or QUAST); default=U50
Optional arguments
            -o | --organism_tag (String of letters that define an organism. This is incorporateed in to the organism ID); default=BAC
            -k | --keep_transposease (T/F); default=T
            -s | --scaffold_genome (T/F); ; default=F
            -t | --pseudogene_predictor (PF, prokka, DFAST or all: choice determines te tools used); default=all
            -e | --reference_genome (Reference genome for genome scaffolding)
            -n | --ancestor_genome (Genome needed for dN/dS calculation)
            -T | --pseudofinder_path
            -M | --microbeannotatorDB (path to Microbeannotator DB) 
            -d | --dnds (dN/dS cutoff for determining a pseudogene with Pseudofinder); default=20
            -p | --stage_start (The stage of starting the pipeline); default=1
            -E | --stage_end (The stage of ending the pipeline); default=8
            -u | --stage_5_specfic_run
            -v | --verbose
            -h | --help
```

## Feedback

Please file questions, bugs or ideas to the [Issue Tracker](https://github.com/mosilab/pseudopipe/issues)

## License

[BSD3-3-Clause license](https://github.com/mosilab/pseudopipe/blob/v0.1.0/LICENSE)

## Citation

Not published yet.

## Author

- Edwin Kyei-Baffour
- Web: https://github.com/eddykay310
