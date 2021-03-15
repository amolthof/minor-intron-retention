## Determining the level of minor intron retention.

The purpose of this bash script is to determine the level of minor intron retention. This is established by quantifying reads aligned to exon-minor intron boundaries. This information is then used in combination with the output of the minor-intron-splicing script to generate mis-splicing indices for each minor intron. 

Dependencies:

BEDTools: https://bedtools.readthedocs.io/en/latest/<br>
awk: https://www.gnu.org/software/gawk/gawk.html<br>
key-merge:
___

**Provide information about sample, SAM/BAM file location and output directory**

    sample="WT"
    genome="GRCh38"

    workdir="/path/to/my/working_directory"
    aligndir="${workdir}/SAM_Files"
    outdir="${workdir}/MinorIntronRetention"

**Provide location of BEDfiles**

Several BEDfiles are used to isolate reads that support proper splicing of minor introns. These were downloaded from the Minor Intron Database [(MIDB)](https://midb.pnb.uconn.edu/) and contain minor intron coordinates corresponding to GRCh38 Ensembl v.84. It is ***essential*** that the genome version used to align the data is the same as the one used to obtain minor intron coordinates.

    BEDdir="/path/to/my/BEDfile_directory"
    IntronType="MinorIntrons"
    BedFileROI=${BEDdir}/
    BedFileIntrons=${BEDdir}/
    BedFile5SSExons=${BEDdir}/
    BedFile3SSExons=${BEDdir}/
    ${BedFileExons}=${BEDdir}/

**Compute intron coverage**

Depending on the size of the BAM file, coverageBed may use a lot of memory. This problem can be bypassed by using a sorted BAM, BED and genome file as input. The genome file needs to consist of all chromosomes listed in the same order as the BAM file header, as well as the size of each chromosome. Be wary of scaffold names, as these might be sorted differently in BAM and txt files.

    echo "Computing intron coverage in ${sample}... |" `date`
	
    FASTAindex="/path/to/my/fasta_directory"
    awk -v OFS="\t" '{print $1, $2}' ${FASTAindex} > ${genome}.txt
    coverageBed -g ${genome}.txt -a ${BedFileIntrons} -b ${aligndir}/${sample}_uniqueUnsplicedReads_sorted.bam -s -split -sorted | awk -F"\t" '{print $4"::"$1":"$2"-"$3"("$6")",$10}' > ${outdir}/${sample}_coverageBed_${IntronType}.bed

**Isolate and quantify reads aligning to exon-intron boundaries resulting from disrupted minor intron splicing**

    echo "Getting reads mapping to exon-intron boundaries... |" `date`
	
    intersectBed -wa -abam ${aligndir}/${sample}_uniqueUnsplicedReads_sorted.bam -b ${BedFileExons} -s > ${outdir}/${sample}_intersectBed_Exons.bam
	intersectBed -wa -abam ${outdir}/${sample}_intersectBed_Exons.bam -b ${BedFile5SSIntrons} -s | intersectBed -wb -abam - -b ${BedFileIntrons} -bed -s | awk '{print $16"::"$13":"$14"-"$15"("$18")"}' | sort | uniq -c | awk '{print $2,$1}' | sort | uniq > ${outdir}/${sample}_intersectBed_Exons_5SSBoundaryFilter.tmp
	intersectBed -wa -abam ${outdir}/${sample}_intersectBed_Exons.bam -b ${BedFile3SSIntrons} -s | intersectBed -wb -abam - -b ${BedFileIntrons} -bed -s | awk '{print $16"::"$13":"$14"-"$15"("$18")"}' | sort | uniq -c | awk '{print $2,$1}' | sort | uniq > ${outdir}/${sample}_intersectBed_Exons_3SSBoundaryFilter.tmp

	key-merge ${outdir}/${IntronType}_key.txt ${outdir}/${sample}_intersectBed_Exons_5SSBoundaryFilter.tmp | awk '{if (NF==2) print $1,$2; else print $1,"0"}' > ${outdir}/${sample}_intersectBed_Exons_5SSBoundaryFilter_allCount.tmp
	key-merge ${outdir}/${IntronType}_key.txt ${outdir}/${sample}_intersectBed_Exons_3SSBoundaryFilter.tmp | awk '{if (NF==2) print $1,$2; else print $1,"0"}' > ${outdir}/${sample}_intersectBed_Exons_3SSBoundaryFilter_allCount.tmp

**Calculate mis-splicing index (MSI) for each minor intron**

    NumIntrons=`awk 'END {print NR}' ${outdir}/${IntronType}_key.txt`

	echo "Getting read counts for ${NumIntrons} introns in ${sample}... |" `date`
	
	for a in `seq 1 ${NumIntrons}`; do

	    Intron=`awk '{print $1}' ${outdir}/${IntronType}_key.txt | head -n ${a} | tail -n 1`
		ReadA1=`awk -v intron=${Intron} '{if ($1==intron) print $2}' ${outdir}/${sample}_intersectBed_Exons_5SSBoundaryFilter_allCount.tmp`
		ReadA2=`awk -v intron=${Intron} '{if ($1==intron) print $2}' ${outdir}/${sample}_intersectBed_Exons_3SSBoundaryFilter_allCount.tmp`
		ReadN=`awk -v intron=${Intron} '{if (substr($1,6)==intron) print $2}' ${outdir}/${sample}_spliced${IntronType}_CAT1_counted.txt`
		Coverage=`awk -v intron=${Intron} '{if ($1==intron) print $2}' ${outdir}/${sample}_coverageBed_${IntronType}.bed`
		MSI=`awk -v A1=${ReadA1} -v A2=${ReadA2} -v N=${ReadN} 'BEGIN {sum=(A1+A2); if (sum==0) print "0"; else print (sum/(sum+2*N))}'`
		Filter=`awk -v A1=${ReadA1} -v A2=${ReadA2} -v c=${Coverage} 'BEGIN {sum=(A1+A2); if (A1>=1 && A2>=1 && sum>=4 && c>0.95) print "yes"; else print "no"}'`

		echo ${Intron} ${sample} ${ReadA1} ${ReadA2} ${ReadN} ${Coverage} >> ${outdir}/${sample}_${IntronType}_counts.tmp
		echo ${Intron} ${sample} ${MSI} ${Filter} >> ${outdir}/${sample}_${IntronType}_MSI.tmp
	done
