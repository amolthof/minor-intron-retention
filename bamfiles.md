## Extracting uniquely aligned, spliced and unspliced reads into BAM files.

The purpose of this bash script is to create two BAM files containing reads that align to a single location in the genome. One of the BAM files will only contain paired-end,  spliced reads originating from exon-exon junctions, whereas the other BAM file will only contain paired-end, unspliced reads that align to exons and/or introns. 

Dependencies:

SAMtools: http://samtools.sourceforge.net/<br>
awk: https://www.gnu.org/software/gawk/gawk.html
___

**Provide info on sample and location of SAM files**

    workdir="/path/to/my/working_directory"
    aligndir="${workdir}/SAM_Files"

    sample="WT"
    IntronType="MinorIntrons"

**Set number of threads to be used. This is done to speed things up**

    threads=8

**Create BAM files containing uniquely mapped, spliced or unspliced reads**

To determine whether reads are spliced (i.e. originate from an exon-exon junction), information in the CIGAR of the SAM file is used. Spliced reads that are anchored on one side with a single nucleotide are discarded because BEDTools cannot handle these appropriately. 

    echo "Isolating all unique, unspliced reads mapping to ${IntronType} in ${sample}... |" `date`

	# Create BAM files with uniquely aligned, unspliced reads
    awk '{if (substr($1,1,1)=="@" || ($NF=="NH:i:1" && index($10,"NNNNN")==0)) print $0}' ${aligndir}/${sample}.sam > ${aligndir}/${sample}_uniqueReads.sam
	awk '{if (substr($1,1,1)=="@" || $6!~"N") print $0}' ${aligndir}/${sample}_uniqueReads.sam > ${aligndir}/${sample}_uniqueUnsplicedReads.sam
	
    samtools view -bS ${aligndir}/${sample}_uniqueUnsplicedReads.sam > ${aligndir}/${sample}_uniqueUnsplicedReads.bam
	samtools sort -@ ${threads} ${aligndir}/${sample}_uniqueUnsplicedReads.bam -o ${aligndir}/${sample}_uniqueUnsplicedReads_sorted.bam

    # Create BAM files with uniquely aligned, spliced reads
    awk '{if (substr($1,1,1)=="@" || $6~"N") print $0}' ${aligndir}/${sample}_uniquereads.sam > ${aligndir}/${sample}_uniqueSplicedReads.sam
	awk '{if (substr($1,1,1)=="@" || $6!~"N1M") print $0}' ${aligndir}/${sample}_uniqueSplicedReads.sam > ${aligndir}/${sample}_uniqueSplicedReads_exclude1M.sam
	
    samtools view -bS ${aligndir}/${sample}_uniqueSplicedReads_exclude1M.sam > ${aligndir}/${sample}_uniqueSplicedReads.bam

**Remove intermediate files not used in downstream analyses**

    rm -f ${aligndir}/${sample}_uniqueReads.sam
	rm -f ${aligndir}/${sample}_uniqueUnsplicedReads.sam
	rm -f ${aligndir}/${sample}_uniqueSplicedReads.sam
	rm -f ${aligndir}/${sample}_uniqueSplicedReads_exclude1M.sam
	rm -f ${aligndir}/${sample}_uniqueUnsplicedReads.bam

