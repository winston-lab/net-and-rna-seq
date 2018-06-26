#!/usr/bin/env python

rule map_to_windows:
    input:
        bg = "coverage/{norm}/{sample}-netseq-{norm}-5end-SENSE.bedgraph",
        chrsizes = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
    output:
        exp = temp("coverage/{norm}/{sample}-netseq-window-{windowsize}-coverage-{norm}.bedgraph"),
    shell: """
        bedtools makewindows -g {input.chrsizes} -w {wildcards.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output.exp}
        """

rule join_window_counts:
    input:
        exp = expand("coverage/{{norm}}/{sample}-netseq-window-{{windowsize}}-coverage-{{norm}}.bedgraph", sample=SAMPLES),
    output:
        exp = "coverage/{norm}/union-bedgraph-window-{windowsize}-{norm}.tsv.gz",
    params:
        names = list(SAMPLES.keys())
    log: "logs/join_window_counts/join_window_counts-{norm}.log"
    shell: """
        (bedtools unionbedg -i {input.exp} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz -f > {output.exp}) &> {log}
        """

rule plotcorrelations:
    input:
        "coverage/{norm}/union-bedgraph-window-{windowsize}-{norm}.tsv.gz",
    output:
        "qual_ctrl/{status}/{condition}-v-{control}/{condition}-v-{control}-netseq-{status}-window-{windowsize}-{norm}-correlations.svg"
    params:
        pcount = lambda wc: 0.01*int(wc.windowsize),
        samplelist = plotcorrsamples
    script:
        "scripts/plotcorr.R"

rule pca_and_cluster:
    input:
        "coverage/{norm}/union-bedgraph-window-{windowsize}-{norm}.tsv.gz",
    output:
        scree = "qual_ctrl/{status}/{condition}-v-{control}/{condition}-v-{control}-netseq-{status}-window-{windowsize}-{norm}-pca_scree.svg",
        pca = "qual_ctrl/{status}/{condition}-v-{control}/{condition}-v-{control}-netseq-{status}-window-{windowsize}-{norm}-pca.svg",
        dist = "qual_ctrl/{status}/{condition}-v-{control}/{condition}-v-{control}-netseq-{status}-window-{windowsize}-{norm}-euclidean_distances.svg",
    params:
        samplelist = plotcorrsamples,
        grouplist = lambda wc: [SAMPLES[sample]["group"] for sample in plotcorrsamples(wc)]
    script:
        "scripts/pca_and_dist.R"


