#!/usr/bin/env python

rule map_to_windows:
    input:
        bg = "coverage/{norm}/{sample}_{assay}-5end-{norm}-SENSE.bedgraph",
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        temp("qual_ctrl/scatter_plots/{assay}_{sample}-{norm}-window-{windowsize}.bedgraph")
    log: "logs/map_to_windows/map_to_windows-{norm}_{sample}_{windowsize}-{assay}.log"
    shell: """
        (bedtools makewindows -g <(faidx {input.fasta} -i chromsizes | awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2; print $1"-minus", $2}}' | LC_COLLATE=C sort -k1,1) -w {wildcards.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output}) &> {log}
        """

rule join_window_counts:
    input:
        lambda wc: expand(f"qual_ctrl/scatter_plots/{ASSAY}_{{sample}}-{wc.norm}-window-{wc.windowsize}.bedgraph", sample=(SAMPLES if wc.norm=="libsizenorm" else SISAMPLES))
    output:
        "qual_ctrl/scatter_plots/{assay}_union-bedgraph-{norm}-window-{windowsize}-allsamples.tsv.gz"
    params:
        names = lambda wc: list(SAMPLES.keys()) if wc.norm=="libsizenorm" else list(SISAMPLES.keys())
    log: "logs/join_window_counts/join_window_counts-{norm}-{windowsize}-{assay}.log"
    shell: """
        (bedtools unionbedg -i {input} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz -f > {output}) &> {log}
        """

rule plot_scatter_plots:
    input:
        "qual_ctrl/scatter_plots/{assay}_union-bedgraph-{norm}-window-{windowsize}-allsamples.tsv.gz"
    output:
        "qual_ctrl/scatter_plots/{condition}-v-{control}/{status}/{condition}-v-{control}_{assay}-{norm}-scatterplots-{status}-window-{windowsize}.svg"
    params:
        pcount = lambda wc: 0.01*int(wc.windowsize),
        samplelist = lambda wc: get_samples(wc.status, wc.norm, [wc.condition, wc.control]),
        assay = {"rnaseq": "RNA-seq",
                 "netseq": "NET-seq"}.get(ASSAY)
    conda: "../envs/tidyverse.yaml"
    script:
        "../scripts/plot_scatter_plots.R"

# rule pca_and_cluster:
#     input:
#         "coverage/{norm}/union-bedgraph-window-{windowsize}-{norm}.tsv.gz",
#     output:
#         scree = "qual_ctrl/{status}/{condition}-v-{control}/{condition}-v-{control}-netseq-{status}-window-{windowsize}-{norm}-pca_scree.svg",
#         pca = "qual_ctrl/{status}/{condition}-v-{control}/{condition}-v-{control}-netseq-{status}-window-{windowsize}-{norm}-pca.svg",
#         dist = "qual_ctrl/{status}/{condition}-v-{control}/{condition}-v-{control}-netseq-{status}-window-{windowsize}-{norm}-euclidean_distances.svg",
#     params:
#         samplelist = plotcorrsamples,
#         grouplist = lambda wc: [SAMPLES[sample]["group"] for sample in plotcorrsamples(wc)]
#     script:
#         "scripts/pca_and_dist.R"


