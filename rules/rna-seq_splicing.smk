#!/usr/bin/env python

localrules:
    filter_intron_alignments,
    get_intron_counts,
    cat_intron_counts

rule filter_intron_alignments:
    input:
        intron_annotation = config["splicing"]["intron_annotation"],
        blacklist_annotation = config["splicing"]["blacklist_annotation"] if config["splicing"]["blacklist_annotation"] else [],
        alignment = lambda wc:  {True:
                                    {   True : f"alignment/{wc.sample}_{ASSAY}-noPCRduplicates-experimental.bam",
                                        False: f"alignment/{wc.sample}_{ASSAY}-noPCRduplicates.bam"
                                    },
                                 False:
                                    {   True: f"alignment/{wc.sample}_{ASSAY}-uniquemappers-experimental.bam",
                                        False: f"alignment/{wc.sample}_{ASSAY}-uniquemappers.bam"
                                    }
                                }.get(config["molecular-barcode"]).get(len(SISAMPLES)>0)
    output:
        temp("splicing/{sample}_filtered_intron_alignments.bam")
    params:
        strand_flag = "-s" if config["sequence-from-5prime"] else "-S"
    log:
        "logs/filter_intron_alignments/filter_intron_alignments_{sample}.log"
    run:
        if config["splicing"]["blacklist_annotation"]:
            shell("(bedtools intersect {params.strand_flag} -a {input.alignment} -b {input.intron_annotation} | \
                    bedtools intersect {params.strand_flag} -v -a stdin -b {input.blacklist_annotation} > \
                    {output}) &> {log}")
        else:
            shell("(bedtools intersect {params.strand_flag} -a {input.alignment} -b {input.intron_annotation} > \
                    {output}) &> {log}")

rule get_intron_counts:
    input:
        intron_annotation = config["splicing"]["intron_annotation"],
        alignment = "splicing/{sample}_filtered_intron_alignments.bam"
    output:
        "splicing/{sample}_intron_counts.tsv"
    params:
        swap_alignment_strand = "" if config["sequence-from-5prime"] else "--swap_alignment_strand",
        closed_intervals = "--open_intervals" if config["splicing"]["open_intervals"] else "",
        min_overhang = config["splicing"]["min_overhang"],
        group_id = lambda wc: SAMPLES[wc.sample]["group"]
    log:
        "logs/get_intron_counts/get_intron_counts_{sample}.log"
    conda:
        "../envs/htseq.yaml"
    shell: """
        (python scripts/get_intron_counts.py \
                --introns {input.intron_annotation} \
                --alignment {input.alignment} \
                {params.swap_alignment_strand} \
                {params.closed_intervals} \
                --sample_id {wildcards.sample} \
                --group_id {params.group_id} \
                --min_overhang {params.min_overhang} \
                --output {output}) &> {log}
        """

rule cat_intron_counts:
    input:
        expand("splicing/{sample}_intron_counts.tsv", sample=PASSING)
    output:
        "splicing/allsamples_intron_counts.tsv.gz"
    log:
        "logs/cat_intron_counts.log"
    shell: """
        (awk 'NR==1 || FNR !=1' {input} | \
            pigz --stdout > {output}) &> {log}
        """

rule compare_intron_retention:
    input:
        data = "splicing/allsamples_intron_counts.tsv.gz"
    output:
        credible_intervals = "splicing/{condition}-v-{control}/{condition}-v-{control}_intron_retention_credible_intervals.tsv",
        qcplots = "splicing/{condition}-v-{control}/{condition}-v-{control}_intron_retention_qcplots.svg",
        results = "splicing/{condition}-v-{control}/{condition}-v-{control}_intron_retention_results.tsv",
        credible_intervals_plot = "splicing/{condition}-v-{control}/{condition}-v-{control}_intron_retention_credible_intervals.svg",
        volcano_plot = "splicing/{condition}-v-{control}/{condition}-v-{control}_intron_retention_volcano.svg"
    params:
        n_trials = config["splicing"]["n_trials"],
        credible_interval_level = config["splicing"]["credible_interval_level"],
        fdr_cutoff = config["splicing"]["fdr_cutoff"],
    conda:
        "../envs/intron_retention.yaml"
    script:
        "../scripts/intron_retention_empirical_bayes.R"

