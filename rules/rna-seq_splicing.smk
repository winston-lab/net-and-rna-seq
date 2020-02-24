#!/usr/bin/env python

localrules:

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
