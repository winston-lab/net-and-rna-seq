#!/usr/bin/env python

import os
from math import log2
import itertools

configfile: "config.yaml"

subworkflow build_annotations:
    workdir: config["genome"]["annotation_workflow"]

configfile: build_annotations("config.yaml")

ASSAY = config["assay"]

SAMPLES = config["samples"]
SISAMPLES = {k:v for k,v in SAMPLES.items() if v["spikein"]}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]}
SIPASSING = {k:v for k,v in PASSING.items() if v["spikein"]}

controlgroups = list(itertools.chain(*[d.values() for d in config["comparisons"]["libsizenorm"]]))
conditiongroups = list(itertools.chain(*[d.keys() for d in config["comparisons"]["libsizenorm"]]))

comparisons_si = config["comparisons"]["spikenorm"]

if comparisons_si:
    controlgroups_si = list(itertools.chain(*[d.values() for d in config["comparisons"]["spikenorm"]]))
    conditiongroups_si = list(itertools.chain(*[d.keys() for d in config["comparisons"]["spikenorm"]]))

COUNTTYPES = ["counts", "sicounts"] if SISAMPLES else ["counts"]
NORMS = ["libsizenorm", "spikenorm"] if SISAMPLES else ["libsizenorm"]

CATEGORIES = ["genic", "antisense", "convergent", "divergent", "intergenic"]

FIGURES = config["figures"]

wildcard_constraints:
    sample = "|".join(re.escape(x) for x in list(SAMPLES.keys())),
    group = "|".join(set(re.escape(v["group"]) for k,v in SAMPLES.items())),
    control = "|".join(set(re.escape(x) for x in controlgroups + (controlgroups_si if comparisons_si else []) + ["all"])),
    condition = "|".join(set(re.escape(x) for x in conditiongroups + (conditiongroups_si if comparisons_si else []) + ["all"])),
    species = "experimental|spikein",
    read_status = "raw|cleaned|aligned|unaligned",
    figure = "|".join(re.escape(x) for x in list(FIGURES.keys())),
    annotation = "|".join(re.escape(x) for x in set(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES]))),
    status = "all|passing",
    counttype= "counts|sicounts",
    norm = "counts|sicounts|libsizenorm|spikenorm",
    readtype = "5end|wholeread",
    windowsize = "\d+",
    direction = "all|up|unchanged|down",
    assay = ASSAY

status_norm_sample_dict = {
    "all":
        {   "libsizenorm" : SAMPLES,
            "spikenorm" : SISAMPLES
        },
    "passing":
        {   "libsizenorm" : PASSING,
            "spikenorm" : SIPASSING
        }
    }

def get_samples(status, norm, groups):
    if "all" in groups:
        return(list(status_norm_sample_dict[status][norm].keys()))
    else:
        return([k for k,v in status_norm_sample_dict[status][norm].items() if v["group"] in groups])

def cluster_samples(status, norm, cluster_groups, cluster_strands):
    ll = []
    for group, strand in zip(cluster_groups, cluster_strands):
        sublist = [k for k,v in status_norm_sample_dict[status][norm].items() if v["group"] in cluster_groups]
        if strand in ["sense", "both"]:
            ll.append([f"{sample}-sense" for sample in sublist])
        if strand in ["antisense", "both"]:
            ll.append([f"{sample}-antisense" for sample in sublist])
    return(list(itertools.chain(*ll)))

include: "rules/net-seq_clean_reads.smk"
include: "rules/net-seq_alignment.smk"
include: "rules/net-seq_genome_coverage.smk"
include: "rules/net-seq_fastqc.smk"
include: "rules/net-seq_library_processing_summary.smk"
include: "rules/net-seq_sample_similarity.smk"
include: "rules/net-seq_datavis.smk"
include: "rules/net-seq_differential_levels.smk"
include: "rules/net-seq_transcript_annotation.smk"
include: "rules/net-seq_transcript_classification.smk"

localrules:
    all

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

rule all:
    input:
        #require config file so that it gets archived
        "config.yaml",
        #FastQC
        f'qual_ctrl/fastqc/{ASSAY}-per_base_sequence_content.svg',
        #alignment
        expand(f"alignment/{{sample}}_{ASSAY}-noPCRduplicates.bam", sample=SAMPLES) if config["random-hexamer"] else expand(f"alignment/{{sample}}_{ASSAY}-uniquemappers.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/{sample}_{assay}-{readtype}-{norm}-{strand}.bw", norm=["counts","libsizenorm"], sample=SAMPLES, readtype=["5end", "wholeread"], strand=["SENSE", "ANTISENSE", "plus", "minus"], assay=ASSAY),
        expand("coverage/{norm}/{sample}_{assay}-{readtype}-{norm}-{strand}.bw", norm=["sicounts","spikenorm"], sample=SISAMPLES, readtype=["5end", "wholeread"], strand=["SENSE", "ANTISENSE", "plus", "minus"], assay=ASSAY),
        #quality control
        f"qual_ctrl/read_processing/{ASSAY}_read-processing-loss.svg",
        expand(f"qual_ctrl/spikein/{ASSAY}_spikein-plots-{{status}}.svg", status=["all", "passing"]) if SISAMPLES else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{assay}}-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all", "passing"], windowsize=config["scatterplot_binsizes"], assay=ASSAY),
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{assay}}-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all", "passing"], windowsize=config["scatterplot_binsizes"], assay=ASSAY) if SISAMPLES and comparisons_si else [],
        #datavis
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{readtype}}/{{assay}}-{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bygroup-sense.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), figure=FIGURES, status=["all","passing"], readtype=["5end", "wholeread"], assay=ASSAY) if config["plot_figures"] and SISAMPLES and comparisons_si else [],
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/{{assay}}-{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bygroup-sense.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, status=["all","passing"], readtype=["5end", "wholeread"], assay=ASSAY) if config["plot_figures"] else [],
        expand(f"diff_exp/{{condition}}-v-{{control}}/libsizenorm/{{condition}}-v-{{control}}_{ASSAY}-libsizenorm-diffexp-results-all.tsv", zip, condition=conditiongroups, control=controlgroups),
        expand(f"diff_exp/{{condition}}-v-{{control}}/spikenorm/{{condition}}-v-{{control}}_{ASSAY}-spikenorm-diffexp-results-all.tsv", zip, condition=conditiongroups_si, control=controlgroups_si) if SISAMPLES and comparisons_si else [],
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_{{assay}}-libsizenorm-diffexp-results-{{category}}-{{direction}}.tsv", zip, condition=conditiongroups, control=controlgroups), category=CATEGORIES, assay=ASSAY, direction=["all", "up", "down", "unchanged"]),
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_{{assay}}-spikenorm-diffexp-results-{{category}}-{{direction}}.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), category=CATEGORIES, assay=ASSAY, direction=["all", "up", "down", "unchanged"]) if SISAMPLES and comparisons_si else []

# rule make_ratio_annotation:
#     input:
#         lambda wc: config["ratios"][wc.ratio]["path"]
#     params:
#         totalsize = lambda wc: config["ratios"][wc.ratio]["numerator"]["upstream"] + config["ratios"][wc.ratio]["numerator"]["dnstream"] + config["ratios"][wc.ratio]["denominator"]["upstream"] + config["ratios"][wc.ratio]["denominator"]["dnstream"],
#     output:
#         "ratios/{ratio}/{ratio}.bed"
#     log: "logs/make_ratio_annotation/make_ratio_annotation-{ratio}.log"
#     shell:  """
#         (bash scripts/makeStrandedBed.sh {input} | awk 'BEGIN{{FS=OFS="\t"}} ($3-$2)>={params.totalsize}' > {output}) &> {log}
#         """

# rule ratio_counts:
#     input:
#         annotation = "ratios/{ratio}/{ratio}.bed",
#         bw = "coverage/libsizenorm/{sample}-netseq-libsizenorm-5end-SENSE.bw"
#     output:
#         dtfile = temp("ratios/{ratio}/{ratio}_{fractype}_{sample}.mat.gz"),
#         matrix = temp("ratios/{ratio}/{ratio}_{fractype}_{sample}.tsv"),
#         melted = temp("ratios/{ratio}/{ratio}_{fractype}_{sample}-melted.tsv.gz"),
#     params:
#         group = lambda wc : SAMPLES[wc.sample]["group"],
#         upstream = lambda wc: config["ratios"][wc.ratio][wc.fractype]["upstream"],
#         dnstream = lambda wc: config["ratios"][wc.ratio][wc.fractype]["dnstream"],
#         refpoint = lambda wc: config["ratios"][wc.ratio][wc.fractype]["refpoint"]
#     threads: config["threads"]
#     log: "logs/ratio_counts/ratio_counts-{ratio}-{fractype}-{sample}.log"
#     shell: """
#         (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize $(echo {params.upstream} + {params.dnstream} | bc) --averageTypeBins sum -p {threads}) &> {log}
#         (Rscript scripts/melt_matrix.R -i {output.matrix} -r TSS --group {params.group} -s {wildcards.sample} -a none -b $(echo {params.upstream} + {params.dnstream} | bc) -u {params.upstream} -o {output.melted}) &>> {log}
#         """

# rule cat_ratio_counts:
#     input:
#         expand("ratios/{{ratio}}/{{ratio}}_{{fractype}}_{sample}-melted.tsv.gz", sample=SAMPLES)
#     output:
#         "ratios/{ratio}/allsamples_{ratio}_{fractype}.tsv.gz"
#     log: "logs/cat_ratio_counts/cat_ratio_counts-{ratio}-{fractype}.log"
#     shell: """
#         (cat {input} > {output}) &> {log}
#         """

# def ratiosamples(wc):
#     dd = SAMPLES if wc.status=="all" else PASSING
#     if wc.condition=="all":
#         return list(dd.keys())
#     else:
#         return [k for k,v in dd.items() if v["group"] in [wc.control, wc.condition]]

# rule plot_ratios:
#     input:
#         numerator = "ratios/{ratio}/allsamples_{ratio}_numerator.tsv.gz",
#         denominator = "ratios/{ratio}/allsamples_{ratio}_denominator.tsv.gz",
#     output:
#         violin = "ratios/{ratio}/{condition}-v-{control}/netseq-{ratio}_{status}_{condition}-v-{control}_violin.svg",
#         ecdf = "ratios/{ratio}/{condition}-v-{control}/netseq-{ratio}_{status}_{condition}-v-{control}_ecdf.svg"
#     params:
#         num_size = lambda wc: config["ratios"][wc.ratio]["numerator"]["upstream"] + config["ratios"][wc.ratio]["numerator"]["dnstream"],
#         den_size = lambda wc: config["ratios"][wc.ratio]["denominator"]["upstream"] + config["ratios"][wc.ratio]["denominator"]["dnstream"],
#         pcount = 1e-3,
#         samplelist = ratiosamples,
#         ratio_label = lambda wc: config["ratios"][wc.ratio]["ratio_name"],
#         num_label = lambda wc: config["ratios"][wc.ratio]["numerator"]["region_label"],
#         den_label = lambda wc: config["ratios"][wc.ratio]["denominator"]["region_label"],
#         annotation_label = lambda wc: config["ratios"][wc.ratio]["label"]
#     script:
#         "scripts/ratio.R"

# rule direction_counts:
#     input:
#         annotation = lambda wc: os.path.dirname(config["directionality"][wc.annotation]["path"]) + "/stranded/" + wc.annotation + "-STRANDED" + os.path.splitext(config["directionality"][wc.annotation]["path"])[1],
#         bw = "coverage/libsizenorm/{sample}-netseq-libsizenorm-5end-{strand}.bw"
#     output:
#         dtfile = temp("directionality/{annotation}/{annotation}_{sample}-{strand}.mat.gz"),
#         matrix = temp("directionality/{annotation}/{annotation}_{sample}-{strand}.tsv"),
#         melted = temp("directionality/{annotation}/{annotation}_{sample}-{strand}-melted.tsv.gz"),
#     params:
#         group = lambda wc : SAMPLES[wc.sample]["group"],
#         upstream = lambda wc: config["directionality"][wc.annotation]["distance"] + 1 if wc.strand=="ANTISENSE" else 0,
#         dnstream = lambda wc: config["directionality"][wc.annotation]["distance"] + 1 if wc.strand=="SENSE" else 0,
#         refpoint = lambda wc: config["directionality"][wc.annotation]["refpoint"]
#     threads: config["threads"]
#     log: "logs/direction_counts/direction_counts-{annotation}-{sample}-{strand}.log"
#     shell: """
#         (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize 1 --averageTypeBins mean -p {threads}) &> {log}
#         (Rscript scripts/melt_matrix.R -i {output.matrix} -r {params.refpoint} --group {params.group} -s {wildcards.sample} -a {wildcards.strand} -b 1 -u $(echo {params.upstream}-1 | bc) -o {output.melted}) &>> {log}
#         """

# rule cat_direction_counts:
#     input:
#         expand("directionality/{{annotation}}/{{annotation}}_{sample}-{strand}-melted.tsv.gz", sample=SAMPLES, strand=["SENSE", "ANTISENSE"])
#     output:
#         "directionality/{annotation}/allsamples_{annotation}.tsv.gz"
#     log: "logs/cat_direction_counts/cat_direction_counts-{annotation}.log"
#     shell: """
#         (cat {input} > {output}) &> {log}
#         """

