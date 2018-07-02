#!/usr/bin/env python

localrules:
    bowtie2_build,
    index_bam,

#align to combined genome with Tophat2 (single genome only if no samples have spike-ins), without reference transcriptome
rule bowtie2_build:
    input:
        fasta = config["combinedgenome"]["fasta"] if SISAMPLES else config["genome"]["fasta"]
    output:
        expand(config["tophat2"]["bowtie2-index-path"] + "/{{basename}}.{num}.bt2", num=[1,2,3,4]),
        expand(config["tophat2"]["bowtie2-index-path"] + "/{{basename}}.rev.{num}.bt2", num=[1,2])
    params:
        idx_path = config["tophat2"]["bowtie2-index-path"],
    conda: "../envs/tophat2.yaml"
    log: "logs/bowtie2_build_{basename}.log"
    shell: """
        (bowtie2-build {input.fasta} {params.idx_path}/{wildcards.basename}) &> {log}
        """

rule align:
    input:
        expand(config["tophat2"]["bowtie2-index-path"] + "/" + (config["combinedgenome"]["name"] if SISAMPLES else config["genome"]["name"]) + ".{num}.bt2", num = [1,2,3,4]),
        expand(config["tophat2"]["bowtie2-index-path"] + "/" + (config["combinedgenome"]["name"] if SISAMPLES else config["genome"]["name"]) + ".rev.{num}.bt2", num=[1,2]),
        fastq = "fastq/cleaned/{sample}_net-seq-clean.fastq.gz" if config["random-hexamer"] else "fastq/cleaned/{sample}_net-seq-trimmed.fastq.gz"
    output:
        aligned = "alignment/{sample}/accepted_hits.bam",
        unaligned = "alignment/{sample}/unmapped.bam",
        summary = "alignment/{sample}/align_summary.txt",
    params:
        idx_path = config["tophat2"]["bowtie2-index-path"],
        basename = config["combinedgenome"]["name"] if SISAMPLES else config["genome"]["name"],
        read_mismatches = config["tophat2"]["read-mismatches"],
        read_gap_length = config["tophat2"]["read-gap-length"],
        read_edit_dist = config["tophat2"]["read-edit-dist"],
        min_anchor_length = config["tophat2"]["min-anchor-length"],
        splice_mismatches = config["tophat2"]["splice-mismatches"],
        min_intron_length = config["tophat2"]["min-intron-length"],
        max_intron_length = config["tophat2"]["max-intron-length"],
        max_insertion_length = config["tophat2"]["max-insertion-length"],
        max_deletion_length = config["tophat2"]["max-deletion-length"],
        max_multihits = config["tophat2"]["max-multihits"],
        segment_mismatches = config["tophat2"]["segment-mismatches"],
        segment_length = config["tophat2"]["segment-length"],
        min_coverage_intron = config["tophat2"]["min-coverage-intron"],
        max_coverage_intron = config["tophat2"]["max-coverage-intron"],
        min_segment_intron = config["tophat2"]["min-segment-intron"],
        max_segment_intron = config["tophat2"]["max-segment-intron"],
    conda: "../envs/tophat2.yaml"
    threads : config["threads"]
    log: "logs/align/align-{sample}.log"
    shell:
        """
        (tophat2 --read-mismatches {params.read_mismatches} --read-gap-length {params.read_gap_length} --read-edit-dist {params.read_edit_dist} -o alignment/{wildcards.sample} --min-anchor-length {params.min_anchor_length} --splice-mismatches {params.splice_mismatches} --min-intron-length {params.min_intron_length} --max-intron-length {params.max_intron_length} --max-insertion-length {params.max_insertion_length} --max-deletion-length {params.max_deletion_length} --num-threads {threads} --max-multihits {params.max_multihits} --library-type fr-firststrand --segment-mismatches {params.segment_mismatches} --no-coverage-search --segment-length {params.segment_length} --min-coverage-intron {params.min_coverage_intron} --max-coverage-intron {params.max_coverage_intron} --min-segment-intron {params.min_segment_intron} --max-segment-intron {params.max_segment_intron} --b2-sensitive {params.idx_path}/{params.basename} {input.fastq}) &> {log}
        """

rule select_unique_mappers:
    input:
        "alignment/{sample}/accepted_hits.bam"
    output:
        "alignment/{sample}_net-seq-uniquemappers.bam"
    threads: config["threads"]
    log: "logs/select_unique_mappers/select_unique_mappers-{sample}.log"
    shell: """
        (samtools view -buh -q 50 -@ {threads} {input} | samtools sort -T .{wildcards.sample} -@ {threads} - > {output}) &> {log}
        """

rule remove_PCR_duplicates:
    input:
        "alignment/{sample}_net-seq-uniquemappers.bam"
    output:
        "alignment/{sample}_net-seq-noPCRduplicates.bam"
    log: "logs/remove_PCR_duplicates/remove_PCR_duplicates-{sample}.log"
    shell: """
        (python scripts/removePCRdupsFromBAM.py {input} {output}) &> {log}
        (rm {input}) &>> {log}
        """

#indexing is required for separating species by samtools view
rule index_bam:
    input:
        "alignment/{sample}_net-seq-noPCRduplicates.bam" if config["random-hexamer"] else "alignment/{sample}_net-seq-uniquemappers.bam"
    output:
        "alignment/{sample}_net-seq-noPCRduplicates.bam.bai" if config["random-hexamer"] else "alignment/{sample}_net-seq-uniquemappers.bam.bai"
    log : "logs/index_bam/index_bam-{sample}.log"
    shell: """
        (samtools index {input}) &> {log}
        """

rule bam_separate_species:
    input:
        bam = "alignment/{sample}_net-seq-noPCRduplicates.bam" if config["random-hexamer"] else "alignment/{sample}_net-seq-uniquemappers.bam",
        bai = "alignment/{sample}_net-seq-noPCRduplicates.bam.bai" if config["random-hexamer"] else "alignment/{sample}_net-seq-uniquemappers.bam.bai",
        chrsizes = config["combinedgenome"]["chrsizes"]
    output:
        "alignment/{sample}_net-seq-noPCRduplicates-{species}.bam" if config["random-hexamer"] else "alignment/{sample}_net-seq-uniquemappers-{species}.bam"
    params:
        filterprefix = lambda wc: config["combinedgenome"]["spikein_prefix"] if wc.species=="experimental" else config["combinedgenome"]["experimental_prefix"],
        prefix = lambda wc: config["combinedgenome"]["experimental_prefix"] if wc.species=="experimental" else config["combinedgenome"]["spikein_prefix"]
    threads: config["threads"]
    log: "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    shell: """
        (samtools view -h {input.bam} $(grep {params.prefix} {input.chrsizes} | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') | grep -v -e 'SN:{params.filterprefix}' | sed 's/{params.prefix}//g' | samtools view -bh -@ {threads} -o {output} -) &> {log}
        """

