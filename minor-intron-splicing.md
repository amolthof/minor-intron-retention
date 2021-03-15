## Determining the level of effective minor intron splicing.

The purpose of this bash script is to determine the level of proper minor intron splicing. This is established by quantifying reads aligned to exon-exon junctions resulting from the splicing of a minor intron.

Dependencies:

BEDTools: https://bedtools.readthedocs.io/en/latest/<br>
awk: https://www.gnu.org/software/gawk/gawk.html<br>
key-merge:
___

**Provide information about sample, SAM/BAM file location and output directory**

    sample="WT"
    workdir="/path/to/my/working_directory"
    aligndir="${workdir}/SAM_Files"
    outdir="${workdir}/MinorIntronRetention"

    mkdir -p ${outdir}


**Provide location of BEDfiles**

Several BEDfiles are used to isolate reads that support proper splicing of minor introns. These were downloaded from the Minor Intron Database [(MIDB)](https://midb.pnb.uconn.edu/) and contain minor intron coordinates corresponding to GRCh38 Ensembl v.84. It is ***essential*** that the genome version used to align the data is the same as the one used to obtain minor intron coordinates.

    BEDdir="/path/to/my/BEDfile_directory"
    IntronType="MinorIntrons"
    BedFileROI=${BEDdir}/
    BedFileIntrons=${BEDdir}/
    BedFile5SSExons=${BEDdir}/
    BedFile3SSExons=${BEDdir}/

**Isolate reads aligning to exon-exon junctions resulting from proper minor intron splicing**

    echo "Extracting exon-exon junction reads in ${sample}...|" `date`
    
    intersectBed -wa -s -abam ${aligndir}/${sample}_uniqueSplicedReads.bam -b ${BedFileROI} > ${outdir}/${sample}_splicedReadsSpanningFlankingExons.bam
    intersectBed -wa -s -split -v -abam ${outdir}/${sample}_splicedReadsSpanningFlankingExons.bam -b ${BedFileIntrons} > ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}.bam
    intersectBed -wa -s -split -abam ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}.bam -b ${BedFile5SSExons} > ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons.bam
    intersectBed -wa -s -split -abam ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons.bam -b ${BedFile3SSExons} > ${OutputFiles}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons_intersect3SSExons.bam

    bedtools bamtobed -i ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons_intersect3SSExons.bam > ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons_intersect3SSExons.bed
    intersectBed -wb -s -a ${outdir}/${sample}_splicedReadsSpanningFlankingExons_spliced${IntronType}_intersect5SSExons_intersect3SSExons.bed -b ${BedFileIntrons} > ${outdir}/${sample}_spliced${IntronType}_CAT1.bed

    awk -F"\t" '{print $10"::"$7":"$8"-"$9"("$12")"}' ${outdir}/${sample}_spliced${IntronType}_CAT1.bed | sort | uniq -c | awk '{print $2,$1}' > ${outdir}/${sample}_spliced${IntronType}_CAT1_count.tmp
    awk -F"\t" '{print $4"::"$1":"$2"-"$3"("$6")"}' ${BedFileIntrons} | sort | uniq > ${outdir}/${IntronType}_key.txt
    key-merge ${outdir}/${IntronType}_key.txt ${outdir}/${sample}_spliced${IntronType}_CAT1_count.tmp | awk '{if (NF==2) print "CAT1-"$1,$2; else print "CAT1-"$1,"0"}' | sort | uniq > ${outdir}/${sample}_spliced${IntronType}_CAT1_counted.txt