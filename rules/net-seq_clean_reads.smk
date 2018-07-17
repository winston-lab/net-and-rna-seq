#!/usr/bin/env python

# process random hexamer, if present
rule extract_molec_barcode:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"]
    output:
        fastq = f"fastq/cleaned/{{sample}}_{ASSAY}-barcode-processed.fastq.gz",
        barcodes = "qual_ctrl/molec_barcode/barcodes-{sample}.tsv",
        ligation = "qual_ctrl/molec_barcode/ligation-{sample}.tsv"
    threads: config["threads"]
    log: "logs/remove_molec_barcode/remove_molec_barcode-{sample}.log"
    shell: """
        (python scripts/extract_molecular_barcode.py {input} {output.fastq} {output.barcodes} {output.ligation}) &> {log}
        """

# trim adapter sequences from 3' end of read
# if RNA-seq, also trim poly-A sequences
# (the end the poly-A is on may be different depending direction of sequencing)
# also do quality trimming, assuming Nextseq 2-color chemistry
# (sequence loss of poly-G sequences from non-Nextseq platforms should be negligible)
rule clean_reads:
    input:
        lambda wc: f"fastq/cleaned/{wc.sample}_{ASSAY}-barcode-processed.fastq.gz" if config["random-hexamer"] else SAMPLES[wc.sample]["fastq"]
    output:
        fastq = f"fastq/cleaned/{{sample}}_{ASSAY}-clean.fastq.gz",
        log = "logs/clean_reads/clean_reads-{sample}.log"
    params:
        adapter = config["cutadapt"]["adapter"],
        trim_qual = config["cutadapt"]["trim_qual"],
        polya_command = lambda wc: {"netseq": [],
                                    "rnaseq": { True: """-a "A{100}" -n 2""",
                                                False:"""-g "A{100}" -n 2"""}
                                    }.get(ASSAY).get(config["sequence-from-5prime"])
    threads: config["threads"]
    shell: """
        (cutadapt --cores={threads} --adapter={params.adapter} {params.polya_command} --cut=-1 --trim-n --nextseq-trim={params.trim_qual} --length-tag='length=' --minimum-length=12 --output={output.fastq} {input}) &> {output.log}
        """

