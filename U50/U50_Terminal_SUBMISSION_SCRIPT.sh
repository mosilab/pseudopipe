#############################################################################
#Author: Christina Castro
#Date: March 31, 2016
#Affiliation: Centers for Disease Control and Prevention
#			  Division of Viral Diseases
#			  Advanced Molecular Detection
#############################################################################

#USAGE: sh U50_SUBMISSION_SCRIPT.sh <reference.fasta> <contigs.fasta>


Reference=${1}
Contigs=${2}
output_dir=$3
u50_opts=$4

#-------------------------------- [ Calculating N50 from de novo assembly] --------------------------------#

quast ${Contigs} -R ${Reference} --min-contig 10 -o $output_dir	$u50_opts 	#Run quast report on the contigs from de novo assembly 

#-------------------------- [ Running a reference assembly, and calculating U50] --------------------------#

bowtie2-build -f ${Reference} REF						#Build reference index

bowtie2 -x REF -f -U ${Contigs} -S $output_dir/out.sam				#Map contigs to reference output sam file

samtools view -bS -o $output_dir/out.bam $output_dir/out.sam			#Convert sam to bam

bedtools bamtobed -i $output_dir/out.bam > $output_dir/out.bed			#Convert bam to bed

sortBed -i $output_dir/out.bed > $output_dir/out_sorted.bed			#Sort bed file by chromosome then by start position in ascending order

python -W ignore ./U50/U50.py ${Reference} $output_dir/out_sorted.bed			#Run U50
