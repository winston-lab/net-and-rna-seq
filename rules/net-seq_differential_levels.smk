#!/usr/bin/env python

localrules:
    map_counts_to_transcripts,
    combine_transcript_counts

rule map_counts_to_transcripts:
    input:
        bed = lambda wc: config["genome"]["transcripts"] if wc.species=="experimental" else config["genome"]["spikein-transcripts"],
        bg = lambda wc: f"coverage/counts/{wc.sample}_netseq-5end-counts-SENSE.bedgraph" if wc.species=="experimental" else f"coverage/sicounts/{wc.sample}_netseq-5end-sicounts-SENSE.bedgraph"
    output:
        temp("diff_exp/{condition}-v-{control}/{sample}_{species}-transcript-counts.tsv")
    log: "logs/map_counts_to_transcripts/map_counts_to_transcripts_{condition}-v-{control}_{sample}-{species}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input.bed} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $4"~"$1"~"$2"~"$3, $7}}' &> {output}) &> {log}
        """

rule combine_transcript_counts:
    input:
        lambda wc : ["diff_exp/{condition}-v-{control}/".format(**wc) + x + f"_{wc.species}-transcript-counts.tsv" for x in get_samples("passing", "libsizenorm", [wc.control, wc.condition])]
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{species}-allsample-transcript-counts.tsv"
    params:
        n = lambda wc: 2*len(get_samples("passing", "libsizenorm", [wc.control, wc.condition])),
        names = lambda wc: "\t".join(get_samples("passing", "libsizenorm", [wc.control, wc.condition]))
    log: "logs/combine_transcript_counts/combine_transcript_counts-{condition}-v-{control}-{species}.log"
    shell: """
        (paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) | cat <(echo -e "name\t" "{params.names}" ) - > {output}) &> {log}
        """

rule call_diffexp_transcripts:
    input:
        expcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-experimental-allsample-transcript-counts.tsv",
        sicounts = lambda wc: [] if wc.norm=="libsizenorm" else "diff_exp/{condition}-v-{control}/{condition}-v-{control}-spikein-allsample-transcript-counts.tsv".format(**wc)
    output:
        results_all = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-results-{norm}-all.tsv",
        results_up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-results-{norm}-up.tsv",
        results_down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-results-{norm}-down.tsv",
        results_unch = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-results-{norm}-unch.tsv",
        bed_all = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-results-{norm}-all.bed",
        bed_up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-results-{norm}-up.bed",
        bed_down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-results-{norm}-down.bed",
        bed_unch = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-results-{norm}-unch.bed",
        normcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-counts-sfnorm-{norm}.tsv",
        rldcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-counts-rlog-{norm}.tsv",
        qcplots = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-netseq-qcplots-{norm}.svg"
    params:
        samples = lambda wc : get_samples("passing", wc.norm, [wc.control, wc.condition]),
        groups = lambda wc : [PASSING[x]["group"] for x in get_samples("passing", wc.norm, [wc.control, wc.condition])],
        alpha = config["deseq"]["fdr"],
        lfc = log2(config["deseq"]["fold-change-threshold"])
    conda: "../envs/diff_exp.yaml"
    script:
        "../scripts/call_diffexp_transcripts.R"

