#!/bin/bash

# Set working and output directory
workdir="/path/to/my/working_directory"
aligndir="${workdir}/SAM_Files"

# Provide sample information and other variables
sample1="WT"
sample2="KO"
IntronType="MinorIntrons"
threads=8

for sample in ${sample1} ${sample2}; do
    echo "Isolating all unique, unspliced reads mapping to ${IntronType} in ${sample}... |" date

    # Create BAM files with uniquely aligned, unspliced reads
    awk '{if (substr($1,1,1)=="@" || ($NF=="NH:i:1" && index($10,"NNNNN")==0)) print $0}' ${aligndir}/${sample}.sam >${aligndir}/${sample}_uniqueReads.sam
    awk '{if (substr($1,1,1)=="@" || $6!~"N") print $0}' ${aligndir}/${sample}_uniqueReads.sam >${aligndir}/${sample}_uniqueUnsplicedReads.sam

    samtools view -bS ${aligndir}/${sample}_uniqueUnsplicedReads.sam >${aligndir}/${sample}_uniqueUnsplicedReads.bam
    samtools sort -@ ${threads} ${aligndir}/${sample}_uniqueUnsplicedReads.bam -o ${aligndir}/${sample}_uniqueUnsplicedReads_sorted.bam

    # Create BAM files with uniquely aligned, spliced reads
    awk '{if (substr($1,1,1)=="@" || $6~"N") print $0}' ${aligndir}/${sample}_uniquereads.sam >${aligndir}/${sample}_uniqueSplicedReads.sam
    awk '{if (substr($1,1,1)=="@" || $6!~"N1M") print $0}' ${aligndir}/${sample}_uniqueSplicedReads.sam >${aligndir}/${sample}_uniqueSplicedReads_exclude1M.sam

    samtools view -bS ${aligndir}/${sample}_uniqueSplicedReads_exclude1M.sam >${aligndir}/${sample}_uniqueSplicedReads.bam

    # Remove intermediate files not used in downstream analyses
    rm -f ${aligndir}/${sample}_uniqueReads.sam
    rm -f ${aligndir}/${sample}_uniqueUnsplicedReads.sam
    rm -f ${aligndir}/${sample}_uniqueSplicedReads.sam
    rm -f ${aligndir}/${sample}_uniqueSplicedReads_exclude1M.sam
    rm -f ${aligndir}/${sample}_uniqueUnsplicedReads.bam
done
