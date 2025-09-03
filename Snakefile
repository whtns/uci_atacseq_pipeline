
EMAIL = "kstachel@uci.edu"

onstart:
   shell("mail -s 'STARTED' {EMAIL} < {log}")

onsuccess:
   shell("mail -s 'DONE' {EMAIL} < {log}")

onerror:
   shell("mail -s 'ERROR' {EMAIL} < {log}")

# ATAC-seq Snakemake Workflow
# This workflow automates quality control, trimming, alignment, filtering, deduplication, chromosome filtering, and peak calling for ATAC-seq data.

import os
import glob



# Load configuration file
configfile: "config.yaml"

# Directory and file paths from config
FASTQ_DIR = config["fastq_dir"]           # Directory containing FASTQ files
TRIMMED_DIR = config["trimmed_dir"]       # Directory for trimmed FASTQ files
BOWTIE2_DIR = config["bowtie2_dir"]       # Directory for Bowtie2 alignment outputs
MACS2_DIR = config["macs2_dir"]           # Directory for MACS2 peak calling outputs
FASTQC_DIR = config["fastqc_dir"]         # Directory for FastQC outputs
ADAPTERS = config["adapters"]             # Adapter file for Trimmomatic
BOWTIE2_INDEX = config["bowtie2_index"]   # Bowtie2 index basename
BLACKLIST = config["blacklist"]           # Blacklist regions file
PICARD_JAR = config["picard_jar"]         # Picard tools jar file
TRIMMOMATIC_JAR = config["trimmomatic"]   # Trimmomatic jar file
BEDGRAPH_TO_BIGWIG = config["bedgraph_to_bigwig"]  # Path to bedGraphToBigWig executable


# Find all sample prefixes (assuming paired-end: *_trimmed_1P.fq.gz and *_trimmed_2P.fq.gz)
def get_samples():
    # Detect samples based on trimmed FASTQ files
    r1_files = glob.glob(os.path.join(FASTQ_DIR, "*_1.fq.gz"))
    samples = [os.path.basename(f).replace("_1.fq.gz", "") for f in r1_files]
    return samples

SAMPLES = get_samples()
print(SAMPLES)  # Print detected sample names


# Final target: MACS2 peak files for all samples
rule all:
    input:
        expand(f"{MACS2_DIR}/{{sample}}_macs2_treat_pileup.bw", sample=SAMPLES)

rule fastqc:
    input:
        r1 = f"{FASTQ_DIR}/{{sample}}_1.fq.gz",
        r2 = f"{FASTQ_DIR}/{{sample}}_2.fq.gz"
    output:
        html = f"{FASTQC_DIR}/{{sample}}_fastqc.html"
    params:
        outdir = FASTQC_DIR
    threads: 4
    shell:
        """
        module load fastqc/0.11.9
        fastqc -t {threads} {input.r1} {input.r2} -o {params.outdir}
        mkdir -p {params.outdir}/html
        mv {params.outdir}/*.html {params.outdir}/html
        module unload fastqc/0.11.9
        """

rule trimmomatic:
    input:
        r1 = f"{FASTQ_DIR}/{{sample}}_1.fq.gz",
        r2 = f"{FASTQ_DIR}/{{sample}}_2.fq.gz"
    output:
        r1_paired = f"{TRIMMED_DIR}/{{sample}}_trimmed_1P.fq.gz",
        r1_unpaired = f"{TRIMMED_DIR}/{{sample}}_trimmed_1U.fq.gz",
        r2_paired = f"{TRIMMED_DIR}/{{sample}}_trimmed_2P.fq.gz",
        r2_unpaired = f"{TRIMMED_DIR}/{{sample}}_trimmed_2U.fq.gz"
    params:
        adapters = ADAPTERS,
        outdir = TRIMMED_DIR,
        jar = TRIMMOMATIC_JAR,
        trimmed = f"{TRIMMED_DIR}/{{sample}}_trimmed.fq.gz"
    threads: 4
    resources:
        mem_mb = 24000
    shell:
        """
        java -jar {params.jar} PE -threads {threads} -phred33 \
        -baseout {params.trimmed} \
        {input.r1} \
        {input.r2} \
        ILLUMINACLIP:{params.adapters}:2:30:10 \
        SLIDINGWINDOW:4:20 \
        MINLEN:20
        """


rule bowtie2_align:
    input:
        trimmed_1p = f"{TRIMMED_DIR}/{{sample}}_trimmed_1P.fq.gz",
        trimmed_2p = f"{TRIMMED_DIR}/{{sample}}_trimmed_2P.fq.gz"
    output:
        bam = f"{BOWTIE2_DIR}/{{sample}}_bowtie2.bam",
        log = f"{BOWTIE2_DIR}/{{sample}}_alignment.log"
    params:
        index = BOWTIE2_INDEX
    threads: 8
    resources:
        mem_mb = 48000
    shell:
        """
        module load bowtie2/2.4.1
        module load samtools/1.10
        mkdir -p {BOWTIE2_DIR}
        bowtie2 --phred33 --threads {threads} --dovetail --maxins 2000 \
        -x {params.index} -1 {input.trimmed_1p} -2 {input.trimmed_2p} \
        2> {output.log} | samtools view -@ {threads} -bhF 4 -f 2 -q 30 > {output.bam}
        module unload bowtie2/2.4.1
        module unload samtools/1.10
        """


# Remove reads overlapping blacklist regions and sort BAM
rule samtools_filter_sort:
    input:
        bam = f"{BOWTIE2_DIR}/{{sample}}_bowtie2.bam",
        blacklist = BLACKLIST
    output:
        sorted_bam = temp(f"{BOWTIE2_DIR}/{{sample}}_noexclu.sorted.bam")
    threads: 8
    resources:
        mem_mb = 48000
    shell:
        """
        module load bedtools2/2.30.0
        module load samtools/1.10
        bedtools intersect -v -abam {input.bam} -b {input.blacklist} | \
        samtools sort -@ {threads} -o {output.sorted_bam}
        module unload bedtools2/2.30.0
        module unload samtools/1.10
        """


# Remove duplicate reads and index BAM using Picard
rule picard_dedup_index:
    input:
        sorted_bam = f"{BOWTIE2_DIR}/{{sample}}_noexclu.sorted.bam"
    output:
        dedup_bam = temp(f"{BOWTIE2_DIR}/{{sample}}_dedup_reads.bam"),
        metrics = f"{BOWTIE2_DIR}/{{sample}}_metrics.txt",
        dedup_bam_index = f"{BOWTIE2_DIR}/{{sample}}_dedup_reads.bam.bai"
    params:
        picard_jar = PICARD_JAR
    threads: 8
    resources:
        mem_mb = 48000
    shell:
        """
        module load samtools/1.10
        module load picard-tools/2.27.1
        java -Xmx2g -jar {params.picard_jar} MarkDuplicates \
        INPUT={input.sorted_bam} OUTPUT={output.dedup_bam} METRICS_FILE={output.metrics} REMOVE_DUPLICATES=True
        samtools index {output.dedup_bam}
        module unload samtools/1.10
        module unload picard-tools/2.27.1
        """


# Filter BAM to retain only standard chromosomes and index
rule bam_chr_filter:
    input:
        dedup_bam = f"{BOWTIE2_DIR}/{{sample}}_dedup_reads.bam"
    output:
        cleaned_bam = f"{BOWTIE2_DIR}/{{sample}}_cleaned.bam",
        cleaned_bam_index = f"{BOWTIE2_DIR}/{{sample}}_cleaned.bam.bai"
    threads: 8
    resources:
        mem_mb = 48000
    shell:
        """
        module load samtools/1.10
        samtools view -b {input.dedup_bam} \
        chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 \
        chr17 chr18 chr19 chr20 chr21 chr22 chrX > {output.cleaned_bam}
        samtools index {output.cleaned_bam}
        module unload samtools/1.10
        """


# Call peaks using MACS2
rule macs2:
    input:
        bam = f"{BOWTIE2_DIR}/{{sample}}_dedup_reads.bam"
    output:
        peaks = f"{MACS2_DIR}/{{sample}}_macs2_peaks.narrowPeak",
        bedgraph = f"{MACS2_DIR}/{{sample}}_macs2_treat_pileup.bdg"
    params:
        outdir = MACS2_DIR
    threads: 8
    resources:
        mem_mb = 48000
    shell:
        """
        module load macs/2.2.7.1
        macs2 callpeak -t {input.bam} \
        -p 0.0001 -g mm -f BAMPE --nomodel --nolambda -B --keep-dup all --call-summits \
        -n {wildcards.sample}_macs2 --outdir {params.outdir}
        echo DONE
        module unload macs/2.2.7.1
        """

## Convert MACS2 bedGraph to bigWig
rule bedgraphtobigwig:
    input:
        bedgraph = f"{MACS2_DIR}/{{sample}}_macs2_treat_pileup.bdg"
    output:
        bigwig = f"{MACS2_DIR}/{{sample}}_macs2_treat_pileup.bw"
    params:
        chrom_sizes = config["chrom_sizes"],  # Add chrom_sizes to config.yaml
        bedgraph_to_bigwig = BEDGRAPH_TO_BIGWIG  # Path to bedGraphToBigWig executable
    threads: 2
    shell:
        """
        {params.bedgraph_to_bigwig} {input.bedgraph} {params.chrom_sizes} {output.bigwig}
        """