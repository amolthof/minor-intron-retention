## Aligning paired-end Illumina sequencing reads against a reference genome using Hisat2.

The purpose of this bash script is to align Illumina-generated paired-end sequencing data against a reference genome. This is done using Hisat2, but other alignment tools such as STAR may also be used to generate SAM files.  

*It is assumed that the FASTQ reads have been checked for quality. Any filtering, adapter trimming etc. that might be required should be done prior to running this script.*

Dependencies:

Hisat2: https://daehwankimlab.github.io/hisat2/ <br>

___
**Set FASTQ, FASTA and output directory**

    fastqdir="/path/to/my/fastq_directory"
    fastadir="/path/to/my/fasta_directory"
    workdir="/path/to/my/working_directory"
    aligndir="${workdir}/SAM_Files"
    mkdir -p ${fastqdir}
    mkdir -p ${fastadir}
    mkdir -p ${aligndir}

**Provide sample information**

    sample="WT"
    organism="Homo_sapiens"

**Set number of threads to be used. This is done to speed things up**

    threads=16

**Generate Hisat2 index**

These were generated using GRCh38 Ensembl v.95 genome sequences downloaded from the Ensembl database. However, genome files obtained from other sources can be used as well.

    nice hisat2-build -p ${threads} ${fastadir}/${organism}*.fa ${fastadir}/${organism}_genome

**Align to the reference genome**

Several of the options used for alignment are not listed in the Hisat2 manual, but originate from the Bowtie2 implementation. It might be necessary to include additional options (e.g. read trimming) or modify current options to optimize the alignment rates. 

    nice hisat2 --end-to-end --no-discordant --no-mixed --sensitive -p ${threads} -x ${fastadir}/${organism}_genome -1 ${fastqdir}/${sample}_R1.fastq -2 ${fastqdir}/${sample}_R2.fastq > ${aligndir}/${sample}.sam



