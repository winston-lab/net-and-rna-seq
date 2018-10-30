#!/usr/bin/env python

localrules:
    build_combined_genome,
    bowtie2_build,
    index_bam,

basename = "{exp_name}_{exp_fasta}_{si_name}_{si_fasta}".format(exp_name = config["genome"]["name"],
                                                                exp_fasta = os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0],
                                                                si_name = config["spike_in"]["name"],
                                                                si_fasta = os.path.splitext(os.path.basename(config["spike_in"]["fasta"]))[0]) if SISAMPLES else os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0]

rule build_combined_genome:
    input:
        experimental = os.path.abspath(build_annotations(config["genome"]["fasta"])),
        spikein = config["spike_in"]["fasta"] if SISAMPLES else []
    output:
        "{directory}/{bn}.fa".format(directory = os.path.split(config["genome"]["fasta"])[0], bn=basename),
    params:
        exp_name = config["genome"]["name"],
        si_name = config["spike_in"]["name"] if SISAMPLES else []
    log: "logs/build_combined_genome.log"
    shell: """
        (sed 's/>/>{params.exp_name}_/g' {input.experimental} | \
        cat - <(sed 's/>/>{params.si_name}_/g' {input.spikein}) > {output}) &> {log}
        """

#align to combined genome with Tophat2 (single genome only if no samples have spike-ins), without reference transcriptome
rule bowtie2_build:
    input:
        "{directory}/{bn}.fa".format(directory = os.path.split(config["genome"]["fasta"])[0], bn=basename) if SISAMPLES else os.path.abspath(build_annotations(config["genome"]["fasta"])),
    output:
        expand(config["tophat2"]["index-path"] + "/{{basename}}.{num}.bt2", num=[1,2,3,4]),
        expand(config["tophat2"]["index-path"] + "/{{basename}}.rev.{num}.bt2", num=[1,2])
    params:
        idx_path = config["tophat2"]["index-path"],
    conda: "../envs/tophat2.yaml"
    log: "logs/bowtie2_build_{basename}.log"
    shell: """
        (bowtie2-build {input} {params.idx_path}/{wildcards.basename}) &> {log}
        """

rule align:
    input:
        expand("{directory}/{bn}.{num}.bt2", directory = config["tophat2"]["index-path"],
                                             bn = basename,
                                             num = [1,2,3,4]),
        expand("{directory}/{bn}.rev.{num}.bt2", directory = config["tophat2"]["index-path"],
                                             bn = basename,
                                             num = [1,2]),
        fastq = f"fastq/cleaned/{{sample}}_{ASSAY}-clean.fastq.gz"
    output:
        aligned = "alignment/{sample}/accepted_hits.bam",
        unaligned = "alignment/{sample}/unmapped.bam",
        summary = "alignment/{sample}/align_summary.txt",
    params:
        idx_path = config["tophat2"]["index-path"],
        library_type = "fr-secondstrand" if ASSAY=="rnaseq" and config["sequence-from-5prime"] else "fr-firststrand",
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
        (tophat2 --read-mismatches {params.read_mismatches} --read-gap-length {params.read_gap_length} --read-edit-dist {params.read_edit_dist} -o alignment/{wildcards.sample} --min-anchor-length {params.min_anchor_length} --splice-mismatches {params.splice_mismatches} --min-intron-length {params.min_intron_length} --max-intron-length {params.max_intron_length} --max-insertion-length {params.max_insertion_length} --max-deletion-length {params.max_deletion_length} --num-threads {threads} --max-multihits {params.max_multihits} --library-type {params.library_type} --segment-mismatches {params.segment_mismatches} --no-coverage-search --segment-length {params.segment_length} --min-coverage-intron {params.min_coverage_intron} --max-coverage-intron {params.max_coverage_intron} --min-segment-intron {params.min_segment_intron} --max-segment-intron {params.max_segment_intron} --b2-sensitive {params.idx_path}/{basename} {input.fastq}) &> {log}
        """

rule select_unique_mappers:
    input:
        "alignment/{sample}/accepted_hits.bam"
    output:
        f"alignment/{{sample}}_{ASSAY}-uniquemappers.bam"
    threads: config["threads"]
    log: "logs/select_unique_mappers/select_unique_mappers-{sample}.log"
    shell: """
        (samtools view -buh -q 50 -@ {threads} {input} | samtools sort -T .{wildcards.sample} -@ {threads} - > {output}) &> {log}
        """

rule remove_PCR_duplicates:
    input:
        f"alignment/{{sample}}_{ASSAY}-uniquemappers.bam"
    output:
        f"alignment/{{sample}}_{ASSAY}-noPCRduplicates.bam"
    log: "logs/remove_PCR_duplicates/remove_PCR_duplicates-{sample}.log"
    shell: """
        (python scripts/removePCRdupsFromBAM.py {input} {output}) &> {log}
        (rm {input}) &>> {log}
        """

#indexing is required for separating species by samtools view
rule index_bam:
    input:
        f"alignment/{{sample}}_{ASSAY}-noPCRduplicates.bam" if config["random-hexamer"] else f"alignment/{{sample}}_{ASSAY}-uniquemappers.bam"
    output:
        f"alignment/{{sample}}_{ASSAY}-noPCRduplicates.bam.bai" if config["random-hexamer"] else f"alignment/{{sample}}_{ASSAY}-uniquemappers.bam.bai"
    log : "logs/index_bam/index_bam-{sample}.log"
    shell: """
        (samtools index {input}) &> {log}
        """

rule bam_separate_species:
    input:
        bam = f"alignment/{{sample}}_{ASSAY}-noPCRduplicates.bam" if config["random-hexamer"] else f"alignment/{{sample}}_{ASSAY}-uniquemappers.bam",
        bai = f"alignment/{{sample}}_{ASSAY}-noPCRduplicates.bam.bai" if config["random-hexamer"] else f"alignment/{{sample}}_{ASSAY}-uniquemappers.bam.bai",
        fasta = "{directory}/{bn}.fa".format(directory = os.path.split(config["genome"]["fasta"])[0], bn=basename) if SISAMPLES else [],
    output:
        f"alignment/{{sample}}_{ASSAY}-noPCRduplicates-{{species}}.bam" if config["random-hexamer"] else f"alignment/{{sample}}_{ASSAY}-uniquemappers-{{species}}.bam"
    params:
        filterprefix = lambda wc: config["spike_in"]["name"] if wc.species=="experimental" else config["genome"]["name"],
        prefix = lambda wc: config["genome"]["name"] if wc.species=="experimental" else config["spike_in"]["name"]
    threads: config["threads"]
    log: "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    shell: """
        (samtools view -h {input.bam} $(faidx {input.fasta} -i chromsizes | grep {params.prefix}_ | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') | grep -v -e 'SN:{params.filterprefix}_' | sed 's/{params.prefix}_//g' | samtools view -bh -@ {threads} -o {output} -) &> {log}
        """

