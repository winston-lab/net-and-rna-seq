#!/usr/bin/env python

rule clean_reads:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"]
    output:
        fastq = "fastq/cleaned/{sample}_net-seq-trimmed.fastq.gz",
        log = "logs/clean_reads/clean_reads-{sample}.log"
    params:
        adapter = config["cutadapt"]["adapter"],
        trim_qual = config["cutadapt"]["trim_qual"]
    threads: config["threads"]
    shell: """
        (cutadapt --cores={threads} --adapter={params.adapter} --cut=-1 --trim-n --nextseq-trim={params.trim_qual} --length-tag='length=' --minimum-length=12 --output={output.fastq} {input}) &> {output.log}
        """

rule extract_molec_barcode:
    input:
        "fastq/cleaned/{sample}_net-seq-trimmed.fastq.gz",
    output:
        fastq = "fastq/cleaned/{sample}_net-seq-clean.fastq.gz",
        barcodes = "qual_ctrl/molec_barcode/barcodes-{sample}.tsv",
        ligation = "qual_ctrl/molec_barcode/ligation-{sample}.tsv"
    threads: config["threads"]
    log: "logs/remove_molec_barcode/remove_molec_barcode-{sample}.log"
    shell: """
        (python scripts/extract_molecular_barcode.py {input} {output.fastq} {output.barcodes} {output.ligation}) &> {log}
        """

