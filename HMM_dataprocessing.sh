#Ancestor GE A and B R1 and R2 are generated using joelle's R1 and R2 scripts on labpubs 
#USAGE: command line input should look like $bash HMM_processing.sh /work/Ancestors/*R1*

mkdir -p sam 
mkdir -p bam
mkdir -p coverage

#use input argument @ for unspecified number of arguments
for arg in "$@"
do
#strips out sample part of filename, can then use the basename to specify the R1 and R2 files with the same basename
        file=$(basename "$arg" .R1.fastq)
#aligns paired end fastqs to reference
	bwa mem -R "@RG\tID:"${file}"" /storage/reference_genomes/tetrahymena_thermophila/mic/mic.genome.fasta "${file}".R1.fastq "${file}".R2.fastq > sam/"${file}".sam
#converts sam to bam, sorts and indexes bam
        samtools view -S -b sam/"${file}".sam > bam/"${file}".bam
        samtools sort -o bam/"${file}"_sorted.bam bam/"${file}".bam
        samtools rmdup bam/"${file}"_sorted.bam bam/"${file}"_sorted_rmdup.bam
        samtools index bam/"${file}"_sorted_rmdup.bam
done 

#merge bams from every sample 
samtools merge bam/AncestorsGE_tomic_sorted_rmdup.bam bam/*_sorted_rmdup.bam 

#get depth per position
samtools depth -a bam/AncestorsGE_tomic_sorted_rmdup.bam > coverage/AncestorsGE_tomic_coverage.txt
#partition out coverage data per chromosome
cd coverage
awk '{print >> ($1 "_AncestorGE_tomic_coverage.txt")}' AncestorGE_tomic_coverage.txt

