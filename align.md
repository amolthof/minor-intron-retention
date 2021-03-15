This bash script provides instructions for aligning paired-end Illumina sequencing reads against a reference genome using Hisat2.
___
The purpose of this script is to obtain alignments for Illumina-generated paired-end sequencing data. This was achieved by aligning the reads to a reference genome using Hisat2. However, other alignment tools such as STAR may also be used to generate SAM files.  

*It is assumed that the FASTQ reads have been checked for quality. Any filtering, adapter trimming etc. that might be required should be done prior to running this script.*

Dependencies:
Hisat2: https://daehwankimlab.github.io/hisat2/
SAMtools: http://samtools.sourceforge.net/
___
**Set FASTQ and result directory**

    fastqdir="/path/to/my/fastq_directory"
    workdir="/path/to/my/working_directory"
    aligndir="${workdir}/SAM_Files"
    mkdir -p ${fastqdir}
    mkdir -p ${aligndir}

**Provide sample information**

    sample="WT"
    organism="Homo_sapiens"

**Set number of threads to be used. This is done to speed things up**

    threads=16

**Generate Hisat2 index**

These were generated using GRCh38 Ensembl v.84 genome sequences downloaded from the Ensembl database. However, genome files obtained from other sources can be used as well.

    hisat2-build 

**Align to the reference genome**

    hisat2-align


