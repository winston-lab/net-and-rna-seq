#!/usr/bin/env python

localrules:
    gene_ontology

rule gene_ontology:
    input:
        universe = config["genome"]["transcript_annotation"],
        diffexp_path = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-{category}-{direction}.tsv",
        go_anno_path = config["gene_ontology_mapping_file"]
    output:
        results = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-{category}-{direction}-gene-ontology-results.tsv",
        enriched_combined = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-{category}-{direction}-gene-ontology-enriched-all.svg",
        depleted_combined = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-{category}-{direction}-gene-ontology-depleted-all.svg",
        enriched_facet = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-{category}-{direction}-gene-ontology-enriched-facetted.svg",
        depleted_facet = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-{category}-{direction}-gene-ontology-depleted-facetted.svg",
    params:
        direction = lambda wc: {"up": "upregulated",
                                "down": "downregulated"}.get(wc.direction, wc.direction)
    conda:
        "../envs/gene_ontology.yaml"
    script:
        "../scripts/gene_ontology.R"

