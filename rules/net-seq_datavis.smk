#!/usr/bin/env python

# rule make_stranded_annotations:
#     input:
#         lambda wc : FIGURES[wc.figure]["annotations"][wc.annotation]["path"]
#     output:
#         "{annopath}/stranded/{figure}_{annotation}-STRANDED.{ext}"
#     log : "logs/make_stranded_annotations/make_stranded_annotations-{figure}_{annotation}_{ext}.log"
#     shell: """
#         (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
#         """

rule compute_matrix:
    input:
        # annotation = lambda wc: os.path.dirname(FIGURES[wc.figure]["annotations"][wc.annotation]["path"]) + "/stranded/" + wc.figure + "_" + wc.annotation + "-STRANDED" + os.path.splitext(FIGURES[wc.figure]["annotations"][wc.annotation]["path"])[1],
        bw = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.bw"
    output:
        dtfile = temp("datavis/{figure}/{norm}/{annotation}_{sample}_{norm}-{readtype}-{strand}.mat.gz"),
        matrix = temp("datavis/{figure}/{norm}/{annotation}_{sample}_{norm}-{readtype}-{strand}.tsv"),
        melted = temp("datavis/{figure}/{norm}/{annotation}_{sample}_{norm}-{readtype}-{strand}-melted.tsv.gz")
    params:
        group = lambda wc : SAMPLES[wc.sample]["group"],
        refpoint = lambda wc: "TSS" if FIGURES[wc.figure]["parameters"]["type"]=="scaled" else FIGURES[wc.figure]["parameters"]["refpoint"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"] + FIGURES[wc.figure]["parameters"]["binsize"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["dnstream"] + FIGURES[wc.figure]["parameters"]["binsize"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        binsize = lambda wc: FIGURES[wc.figure]["parameters"]["binsize"],
        binstat = lambda wc: FIGURES[wc.figure]["parameters"]["binstat"],
        nan_afterend = lambda wc: [] if FIGURES[wc.figure]["parameters"]["type"]=="scaled" or not FIGURES[wc.figure]["parameters"]["nan_afterend"] else "--nanAfterEnd",
        anno_label = lambda wc: FIGURES[wc.figure]["annotations"][wc.annotation]["label"]
    threads : config["threads"]
    log: "logs/compute_matrix/compute_matrix-{figure}_{annotation}_{sample}_{norm}-{strand}.log"
    run:
        if FIGURES[wildcards.figure]["parameters"]["type"]=="absolute":
            shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} {params.nan_afterend} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {params.scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {params.refpoint} -g {params.group} -s {wildcards.sample} -a {params.anno_label} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        lambda wc: expand("datavis/{figure}/{norm}/{annotation}_{sample}_{norm}-{readtype}-{strand}-melted.tsv.gz", annotation=[k for k,v in FIGURES[wc.figure]["annotations"].items()], sample=SAMPLES, figure=wc.figure, norm=wc.norm, strand=wc.strand, readtype=wc.readtype)
    output:
        "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-{norm}-{readtype}-{strand}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{figure}_{norm}-{readtype}-{strand}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_figures:
    input:
        matrices = expand("datavis/{{figure}}/{{norm}}/{{figure}}-allsamples-allannotations-{{norm}}-{{readtype}}-{strand}.tsv.gz", strand=["SENSE", "ANTISENSE"]),
        annotations = lambda wc: [v["path"] for k,v in FIGURES[wc.figure]["annotations"].items()]
    output:
        heatmap_sample_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-heatmap-bysample-bothstrands.svg",
        heatmap_sample_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-heatmap-bysample-sense.svg",
        heatmap_sample_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-heatmap-bysample-antisense.svg",
        heatmap_group_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-heatmap-bygroup-bothstrands.svg",
        heatmap_group_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-heatmap-bygroup-sense.svg",
        heatmap_group_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-heatmap-bygroup-antisense.svg",
        metagene_sample_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bysample-bothstrands.svg",
        metagene_sample_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bysample-sense.svg",
        metagene_sample_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bysample-antisense.svg",
        metagene_sample_overlay_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-sampleoverlay-bothstrands.svg",
        metagene_sample_overlay_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-sampleoverlay-sense.svg",
        metagene_sample_overlay_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-sampleoverlay-antisense.svg",
        metagene_group_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bygroup-bothstrands.svg",
        metagene_group_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bygroup-sense.svg",
        metagene_group_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-bygroup-antisense.svg",
        metagene_sampleanno_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byanno-sample-bothstrands.svg",
        metagene_sampleanno_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byanno-sample-sense.svg",
        metagene_sampleanno_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byanno-sample-antisense.svg",
        metagene_groupanno_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byanno-group-bothstrands.svg",
        metagene_groupanno_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byanno-group-sense.svg",
        metagene_groupanno_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byanno-group-antisense.svg",
        metagene_sampleclust_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byclustersample-bothstrands.svg",
        metagene_sampleclust_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byclustersample-sense.svg",
        metagene_sampleclust_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byclustersample-antisense.svg",
        metagene_groupclust_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byclustergroup-bothstrands.svg",
        metagene_groupclust_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byclustergroup-sense.svg",
        metagene_groupclust_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{readtype}/netseq-{figure}-{norm}-{status}_{condition}-v-{control}_{readtype}-metagene-byclustergroup-antisense.svg",
    params:
        # abusing snakemake a bit here...using params as output paths in order to use lambda functions
        annotations_out = lambda wc: ["datavis/" + wc.figure + "/" + wc.norm + "/" + wc.condition + "-v-" + wc.control + "/" + wc.status + "/" + wc.readtype + "/" + annotation + "_cluster-" + str(cluster) + ".bed" for annotation in FIGURES[wc.figure]["annotations"] for cluster in range(1, FIGURES[wc.figure]["annotations"][annotation]["n_clusters"]+1)],
        clusters_out = lambda wc: ["datavis/" + wc.figure + "/" + wc.norm + "/" + wc.condition + "-v-" + wc.control + "/" + wc.status + "/" + wc.readtype + "/" + annotation + ".pdf" for annotation in FIGURES[wc.figure]["annotations"]],
        samplelist = plotcorrsamples,
        plottype = lambda wc: FIGURES[wc.figure]["parameters"]["type"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["dnstream"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        pct_cutoff = lambda wc: FIGURES[wc.figure]["parameters"]["pct_cutoff"],
        log_transform = lambda wc: str(FIGURES[wc.figure]["parameters"]["log_transform"]).upper(),
        pcount = lambda wc: 0 if not FIGURES[wc.figure]["parameters"]["log_transform"] else FIGURES[wc.figure]["parameters"]["pseudocount"],
        meta_spread = lambda wc: FIGURES[wc.figure]["parameters"]["meta_spread"],
        trim_pct = lambda wc: FIGURES[wc.figure]["parameters"]["trim_pct"],
        refpointlabel = lambda wc: FIGURES[wc.figure]["parameters"]["refpointlabel"],
        endlabel = lambda wc:  "HAIL SATAN" if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["three_prime_label"],
        cmap = lambda wc: FIGURES[wc.figure]["parameters"]["heatmap_colormap"],
        sortmethod = lambda wc: FIGURES[wc.figure]["parameters"]["arrange"],
        cluster_scale = lambda wc: "FALSE" if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else str(FIGURES[wc.figure]["parameters"]["cluster_scale"]).upper(),
        cluster_samples = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else cluster_samples(wc.status, wc.norm, FIGURES[wc.figure]["parameters"]["cluster_conditions"], FIGURES[wc.figure]["parameters"]["cluster_strands"]),
        cluster_five = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_five"],
        cluster_three = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_three"],
        k = lambda wc: [v["n_clusters"] for k,v in FIGURES[wc.figure]["annotations"].items()],
    script:
        "scripts/plot_netseq_figures.R"

