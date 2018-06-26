#!/usr/bin/env python

rule map_counts_to_transcripts:
    input:
        bed = lambda wc: config["genome"]["transcripts"] if wc.type=="exp" else config["genome"]["spikein-transcripts"],
        bg = lambda wc: "coverage/counts/" + wc.sample + "-netseq-counts-5end-SENSE.bedgraph" if wc.type=="exp" else "coverage/sicounts/" + wc.sample + "-netseq-sicounts-5end-SENSE.bedgraph"
    output:
        temp("diff_exp/{condition}-v-{control}/{sample}-{type}-alltranscriptcounts.tsv")
    log: "logs/map_counts_to_transcripts/map_counts_to_transcripts-{condition}-v-{control}-{sample}-{type}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input.bed} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $4"~"$1"~"$2"~"$3, $7}}' &> {output}) &> {log}
        """

rule get_transcript_counts:
    input:
        lambda wc : ["diff_exp/" + wc.condition + "-v-" + wc.control + "/" + x + "-" + wc.type + "-alltranscriptcounts.tsv" for x in getsamples(wc.control, wc.condition)]
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{type}-transcript-counts.tsv"
    params:
        n = lambda wc: 2*len(getsamples(wc.control, wc.condition)),
        names = lambda wc: "\t".join(getsamples(wc.control, wc.condition))
    log: "logs/get_transcript_counts/get_transcript_counts-{condition}-v-{control}-{type}.log"
    shell: """
        (paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) | cat <(echo -e "name\t" "{params.names}" ) - > {output}) &> {log}
        """

rule call_de_transcripts:
    input:
        expcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-exp-transcript-counts.tsv",
        sicounts = lambda wc: "diff_exp/" + wc.condition + "-v-" + wc.control + "/" + wc.condition + "-v-" + wc.control + "-si-transcript-counts.tsv" if wc.norm=="spikenorm" else "diff_exp/" + wc.condition + "-v-" + wc.control + "/" + wc.condition + "-v-" + wc.control + "-exp-transcript-counts.tsv"
    params:
        samples = lambda wc : getsamples(wc.control, wc.condition),
        groups = lambda wc : [PASSING[x]["group"] for x in getsamples(wc.control, wc.condition)],
        alpha = config["deseq"]["fdr"],
        lfc = log2(config["deseq"]["fold-change-threshold"])
    output:
        results_all = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv",
        results_up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-up.tsv",
        results_down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-down.tsv",
        results_unch = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-unch.tsv",
        bed_all = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.bed",
        bed_up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-up.bed",
        bed_down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-down.bed",
        bed_unch = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-unch.bed",
        normcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-counts-sfnorm-{norm}.tsv",
        rldcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-counts-rlog-{norm}.tsv",
        qcplots = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-qcplots-{norm}.svg"
    script:
        "scripts/call_de_transcripts.R"



