#!/bin/bash

# Set working and output directory
workdir="/path/to/my/working_directory"
aligndir="${workdir}/SAM_Files"
fastafile="/path/to/my/fasta_file"
outdir="${workdir}/MinorIntronRetention"

# Provide sample information and other variables
sample1="WT"
sample2="KO"
organism="Homo_sapiens"
genome="hg38Ens95"

# Provide location of BEDfiles
BEDdir="/path/to/my/BEDfile_directory"
IntronType="MinorIntrons"
BedFileIntrons=${BEDdir}/${organism}.${genome}_${IntronType}_Introns.bed
BedFile5SSIntrons=${BEDdir}/${organism}.${genome}_${IntronType}_5SSIntrons.bed
BedFile3SSIntrons=${BEDdir}/${organism}.${genome}_${IntronType}_3SSIntrons.bed
BedFileExons=${BEDdir}/${organism}.${genome}_${IntronType}_FlankingExons.bed

# Calculate MSI values for intron retention
for sample in ${sample1} ${sample2}; do
    # Compute intron coverage
    echo "Computing intron coverage in ${sample}... |" date

    samtools faidx ${fastafile} -o ${outdir}/${organism}_${genome}.fai
    awk -v OFS="\t" '{print $1, $2}' ${outdir}/${organism}_${genome}.fai >${outdir}/${organism}_${genome}.txt
    coverageBed -g ${outdir}/${organism}_${genome}.txt -a ${BedFileIntrons} -b ${aligndir}/${sample}_uniqueUnsplicedReads_sorted.bam -s -split -sorted | awk -F"\t" '{print $4"::"$1":"$2"-"$3"("$6")",$10}' >${outdir}/${sample}_coverageBed_${IntronType}.bed

    # Isolate and quantify reads aligning to exon-intron boundaries resulting from disrupted minor intron splicing
    echo "Getting reads mapping to exon-intron boundaries... |" date

    intersectBed -wa -abam ${aligndir}/${sample}_uniqueUnsplicedReads_sorted.bam -b ${BedFileExons} -s >${outdir}/${sample}_intersectBed_Exons.bam
    intersectBed -wa -abam ${outdir}/${sample}_intersectBed_Exons.bam -b ${BedFile5SSIntrons} -s | intersectBed -wb -abam - -b ${BedFileIntrons} -bed -s | awk '{print $16"::"$13":"$14"-"$15"("$18")"}' | sort | uniq -c | awk '{print $2,$1}' | sort | uniq >${outdir}/${sample}_intersectBed_Exons_5SSBoundaryFilter.tmp
    intersectBed -wa -abam ${outdir}/${sample}_intersectBed_Exons.bam -b ${BedFile3SSIntrons} -s | intersectBed -wb -abam - -b ${BedFileIntrons} -bed -s | awk '{print $16"::"$13":"$14"-"$15"("$18")"}' | sort | uniq -c | awk '{print $2,$1}' | sort | uniq >${outdir}/${sample}_intersectBed_Exons_3SSBoundaryFilter.tmp

    key-merge ${outdir}/${IntronType}_key.txt ${outdir}/${sample}_intersectBed_Exons_5SSBoundaryFilter.tmp | awk '{if (NF==2) print $1,$2; else print $1,"0"}' >${outdir}/${sample}_intersectBed_Exons_5SSBoundaryFilter_allCount.tmp
    key-merge ${outdir}/${IntronType}_key.txt ${outdir}/${sample}_intersectBed_Exons_3SSBoundaryFilter.tmp | awk '{if (NF==2) print $1,$2; else print $1,"0"}' >${outdir}/${sample}_intersectBed_Exons_3SSBoundaryFilter_allCount.tmp

    NumIntrons=$(awk 'END {print NR}' ${outdir}/${IntronType}_key.txt)

    # Calculate mis-splicing index (MSI) for each minor intron
    echo "Getting read counts for ${NumIntrons} introns in ${sample}... |" date

    for a in $(seq 1 "${NumIntrons}"); do

        Intron=$(awk '{print $1}' ${outdir}/${IntronType}_key.txt | head -n "${a}" | tail -n 1)
        ReadA1=$(awk -v intron="${Intron}" '{if ($1==intron) print $2}' ${outdir}/${sample}_intersectBed_Exons_5SSBoundaryFilter_allCount.tmp)
        ReadA2=$(awk -v intron="${Intron}" '{if ($1==intron) print $2}' ${outdir}/${sample}_intersectBed_Exons_3SSBoundaryFilter_allCount.tmp)
        ReadN=$(awk -v intron="${Intron}" '{if (substr($1,6)==intron) print $2}' ${outdir}/${sample}_spliced${IntronType}_CAT1_counted.txt)
        Coverage=$(awk -v intron="${Intron}" '{if ($1==intron) print $2}' ${outdir}/${sample}_coverageBed_${IntronType}.bed)
        MSI=$(awk -v A1="${ReadA1}" -v A2="${ReadA2}" -v N="${ReadN}" 'BEGIN {sum=(A1+A2); if (sum==0) print "0"; else print (sum/(sum+2*N))}')
        Filter=$(awk -v A1="${ReadA1}" -v A2="${ReadA2}" -v c="${Coverage}" 'BEGIN {sum=(A1+A2); if (A1>=1 && A2>=1 && sum>=4 && c>0.95) print "yes"; else print "no"}')

        echo "${Intron} ${sample} ${ReadA1} ${ReadA2} ${ReadN} ${Coverage}" >>${outdir}/${sample}_${IntronType}_counts.tmp
        echo "${Intron} ${sample} ${MSI} ${Filter}" >>${outdir}/${sample}_${IntronType}_MSI.tmp
    done
done

# Combine MSI values of samples in one master output file
key-merge ${outdir}/*_${IntronType}_counts.tmp >${outdir}/${IntronType}_counts_allSamples.txt
key-merge ${outdir}/*_${IntronType}_MSI.tmp >${outdir}/${IntronType}_MSI_allSamples.txt

awk '{pass="no"; for (i=4;i<=NF;i+=3) {if ($i=="yes") pass="yes"}; if (pass=="yes") print $0}' ${outdir}/${IntronType}_MSI_allSamples.txt >${outdir}/${IntronType}_MSI_passInOneSample.txt

rm -f ${outdir}/*.tmp
