#!/bin/bash

# Set FASTQ, FASTA and output directory
fastqdir="/path/to/my/fastq_directory"
fastadir="/path/to/my/fasta_directory"
workdir="/path/to/my/working_directory"
aligndir="${workdir}/SAM_Files"

mkdir -p ${fastqdir}
mkdir -p ${fastadir}
mkdir -p ${aligndir}

# Provide sample information and other variables
sample1="WT"
sample2="KO"
organism="Homo_sapiens"
threads=16

# Generate Hisat2 index
nice hisat2-build -p ${threads} ${fastadir}/${organism}*.fa ${fastadir}/${organism}_genome

# Align reads to reference genome with Hisat2
for sample in ${sample1} ${sample2}; do
    nice hisat2 --end-to-end --no-discordant --no-mixed --sensitive -p ${threads} -x ${fastadir}/${organism}_genome -1 ${fastqdir}/${sample}_R1.fastq -2 ${fastqdir}/${sample}_R2.fastq >${aligndir}/${sample}.sam
done
