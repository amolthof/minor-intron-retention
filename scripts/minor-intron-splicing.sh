#!/bin/bash

# Set working and output directory
workdir="/path/to/my/working_directory"
aligndir="${workdir}/SAM_Files"
outdir="${workdir}/MinorIntronRetention"

mkdir -p ${outdir}

# Provide sample information and other variables
sample1="WT"
sample2="KO"
organism="Homo_sapiens"
genome="hg38Ens95"

# Provide location of BEDfiles
BEDdir="/path/to/my/BEDfile_directory"
IntronType="MinorIntrons"
BedFileROI=${BEDdir}/${organism}.${genome}_${IntronType}_RegionOfInterest.bed
BedFileIntrons=${BEDdir}/${organism}.${genome}_${IntronType}_Introns.bed
BedFile5SSExons=${BEDdir}/${organism}.${genome}_${IntronType}_5SSExons.bed
BedFile3SSExons=${BEDdir}/${organism}.${genome}_${IntronType}_3SSExons.bed

# Isolate reads aligning to exon-exon junctions resulting from proper minor intron splicing
for sample in ${sample1} ${sample2}; do
    echo "Extracting exon-exon junction reads in ${sample}... |" date

    intersectBed -wa -s -abam ${aligndir}/${sample}_uniqueSplicedReads.bam -b ${BedFileROI} >${outdir}/${sample}_splicedReadsSpanningFlankingExons.bam
    intersectBed -wa -s -split -v -abam ${outdir}/${sample}_splicedReadsSpanningFlankingExons.bam -b ${BedFileIntrons} >${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}.bam
    intersectBed -wa -s -split -abam ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}.bam -b ${BedFile5SSExons} >${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons.bam
    intersectBed -wa -s -split -abam ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons.bam -b ${BedFile3SSExons} >${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons_intersect3SSExons.bam

    bedtools bamtobed -i ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons_intersect3SSExons.bam >${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons_intersect3SSExons.bed
    intersectBed -wb -s -a ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons_intersect3SSExons.bed -b ${BedFileIntrons} >${outdir}/${sample}_spliced${IntronType}_CAT1.bed

    awk -F"\t" '{print $10"::"$7":"$8"-"$9"("$12")"}' ${outdir}/${sample}_spliced${IntronType}_CAT1.bed | sort | uniq -c | awk '{print $2,$1}' >${outdir}/${sample}_spliced${IntronType}_CAT1_count.tmp
    awk -F"\t" '{print $4"::"$1":"$2"-"$3"("$6")"}' ${BedFileIntrons} | sort | uniq >${outdir}/${IntronType}_key.txt
    key-merge ${outdir}/${IntronType}_key.txt ${outdir}/${sample}_spliced${IntronType}_CAT1_count.tmp | awk '{if (NF==2) print "CAT1-"$1,$2; else print "CAT1-"$1,"0"}' | sort | uniq >${outdir}/${sample}_spliced${IntronType}_CAT1_counted.txt
done
