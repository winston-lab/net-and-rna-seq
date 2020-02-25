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

comparisons = config["comparisons"]["libsizenorm"]
if comparisons:
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
    control = "|".join(set(re.escape(x) for x in (controlgroups if comparisons else []) + (controlgroups_si if comparisons_si else []) + ["all"])),
    condition = "|".join(set(re.escape(x) for x in (conditiongroups if comparisons else []) + (conditiongroups_si if comparisons_si else []) + ["all"])),
    species = "experimental|spikein",
    read_status = "raw|cleaned|aligned|unaligned",
    figure = "|".join(re.escape(x) for x in list(FIGURES.keys())),
    annotation = "|".join(re.escape(x) for x in set(list(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES])) + list(config["differential_expression"]["annotations"].keys() if config["differential_expression"]["annotations"] else []) + ["transcripts"])),
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
include: "rules/rna-seq_splicing.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules: all

def statuscheck(dict1, dict2):
    return(["passing"] if dict1 == dict2 else ["all", "passing"])

def conditioncheck(conditionlist):
    return(conditionlist if len(conditionlist)==1 else conditionlist + ["all"])

rule all:
    input:
        #require config file so that it gets archived
        "config.yaml",
        #FastQC
        f'qual_ctrl/fastqc/{ASSAY}-per_base_sequence_content.svg',
        #alignment
        expand(f"alignment/{{sample}}_{ASSAY}-noPCRduplicates.bam", sample=SAMPLES) if config["molecular-barcode"] else expand(f"alignment/{{sample}}_{ASSAY}-uniquemappers.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/{sample}_{assay}-{readtype}-{norm}-{strand}.bw", norm=["counts","libsizenorm"], sample=SAMPLES, readtype=["5end", "wholeread"], strand=["SENSE", "ANTISENSE", "plus", "minus"], assay=ASSAY),
        expand("coverage/{norm}/{sample}_{assay}-{readtype}-{norm}-{strand}.bw", norm=["sicounts","spikenorm"], sample=SISAMPLES, readtype=["5end", "wholeread"], strand=["SENSE", "ANTISENSE", "plus", "minus"], assay=ASSAY),
        #quality control
        f"qual_ctrl/read_processing/{ASSAY}_read-processing-loss.svg",
        expand(f"qual_ctrl/spikein/{ASSAY}_spikein-plots-{{status}}.svg", status=statuscheck(SISAMPLES, SIPASSING)) if SISAMPLES else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{assay}}-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), status=statuscheck(SAMPLES, PASSING), windowsize=config["scatterplot_binsizes"], assay=ASSAY) if comparisons else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{assay}}-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), status=statuscheck(SISAMPLES, SIPASSING), windowsize=config["scatterplot_binsizes"], assay=ASSAY) if SISAMPLES and comparisons_si else [],
        #datavis
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{readtype}}/{{assay}}-{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bygroup-sense.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), figure=FIGURES, status=statuscheck(SISAMPLES, SIPASSING), readtype=["5end", "wholeread"], assay=ASSAY) if config["plot_figures"] and SISAMPLES and comparisons_si else [],
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/{{assay}}-{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bygroup-sense.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), figure=FIGURES, status=statuscheck(SAMPLES, PASSING), readtype=["5end", "wholeread"], assay=ASSAY) if config["plot_figures"] and comparisons else [],
        #differential expression
        expand(expand("diff_exp/transcripts/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_{{assay}}-libsizenorm-transcripts-diffexp-results-{{category}}-{{direction}}.tsv", zip, condition=conditiongroups, control=controlgroups), category=CATEGORIES, assay=ASSAY, direction=["all", "up", "down", "unchanged"]) if comparisons else [],
        expand(expand("diff_exp/transcripts/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_{{assay}}-spikenorm-transcripts-diffexp-results-{{category}}-{{direction}}.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), category=CATEGORIES, assay=ASSAY, direction=["all", "up", "down", "unchanged"]) if SISAMPLES and comparisons_si else [],
        expand(expand("diff_exp/{{annotation}}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_{{assay}}-libsizenorm-{{annotation}}-diffexp-results-{{direction}}.tsv", zip, condition=conditiongroups, control=controlgroups), assay=ASSAY, direction=["all", "up", "down", "unchanged"], annotation=list(config["differential_expression"]["annotations"].keys()) if config["differential_expression"]["annotations"] else []) if comparisons else [],
        expand(expand("diff_exp/{{annotation}}/{condition}-v-{control}/spikenorm/{condition}-v-{control}_{{assay}}-spikenorm-{{annotation}}-diffexp-results-{{direction}}.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), assay=ASSAY, direction=["all", "up", "down", "unchanged"], annotation=list(config["differential_expression"]["annotations"].keys()) if config["differential_expression"]["annotations"] else []) if SISAMPLES and comparisons_si else [],
        #differential expression summary
        expand(expand("diff_exp/transcripts/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_{{assay}}-libsizenorm-diffexp-{{plot}}.svg", zip, condition=conditiongroups, control=controlgroups), plot = ["mosaic", "maplot", "volcano"], assay=ASSAY) if comparisons else [],
        expand(expand("diff_exp/transcripts/{condition}-v-{control}/spikenorm/{condition}-v-{control}_{{assay}}-spikenorm-diffexp-{{plot}}.svg", zip, condition=conditiongroups_si, control=controlgroups_si), plot = ["mosaic", "maplot", "volcano"], assay=ASSAY) if SISAMPLES and comparisons_si else [],
        #splicing
        expand("splicing/{condition}-v-{control}/{condition}-v-{control}_intron_retention_results.tsv", zip, condition=conditiongroups, control=controlgroups) if config["assay"]=="rnaseq" and config["analyze_splicing"] and comparisons else [],
        expand("splicing/{condition}-v-{control}/{condition}-v-{control}_intron_retention_results.tsv", zip, condition=conditiongroups_si, control=controlgroups_si) if config["assay"]=="rnaseq" and config["analyze_splicing"] and comparisons_si else [],



