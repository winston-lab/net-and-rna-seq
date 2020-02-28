#!/usr/bin/env python

localrules:
    gene_ontology_called_transcripts,
    gene_ontology_custom_annotation,
    gsea

rule gene_ontology_called_transcripts:
    input:
        universe = config["genome"]["transcript_annotation"],
        diffexp_path = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-{category}-{direction}.tsv",
        go_anno_path = config["gene_ontology_mapping_file"]
    output:
        results = "gene_ontology/transcripts/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-transcripts-{category}-{direction}-gene-ontology-results.tsv",
        enriched_combined = "gene_ontology/transcripts/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-transcripts-{category}-{direction}-gene-ontology-enriched-all.svg",
        depleted_combined = "gene_ontology/transcripts/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-transcripts-{category}-{direction}-gene-ontology-depleted-all.svg",
        enriched_facet = "gene_ontology/transcripts/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-transcripts-{category}-{direction}-gene-ontology-enriched-facetted.svg",
        depleted_facet = "gene_ontology/transcripts/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-transcripts-{category}-{direction}-gene-ontology-depleted-facetted.svg",
    params:
        direction = lambda wc: {"up": "upregulated",
                                "down": "downregulated"}.get(wc.direction, wc.direction),
        category = lambda wc: wc.category,
    conda:
        "../envs/gene_ontology.yaml"
    script:
        "../scripts/gene_ontology.R"

rule gene_ontology_custom_annotation:
    input:
        universe = lambda wc: config["differential_expression"]["annotations"][wc.annotation]["annotation"],
        diffexp_path = "diff_exp/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{annotation}-diffexp-results-{direction}.tsv",
        go_anno_path = lambda wc: config["differential_expression"]["annotations"][wc.annotation]["gene_ontology_mapping_file"],
    output:
        results = "gene_ontology/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{annotation}-{direction}-gene-ontology-results.tsv",
        enriched_combined = "gene_ontology/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{annotation}-{direction}-gene-ontology-enriched-all.svg",
        depleted_combined = "gene_ontology/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{annotation}-{direction}-gene-ontology-depleted-all.svg",
        enriched_facet = "gene_ontology/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{annotation}-{direction}-gene-ontology-enriched-facetted.svg",
        depleted_facet = "gene_ontology/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{annotation}-{direction}-gene-ontology-depleted-facetted.svg",
    params:
        direction = lambda wc: {"up": "upregulated",
                                "down": "downregulated"}.get(wc.direction, wc.direction),
        category = "genic",
    conda:
        "../envs/gene_ontology.yaml"
    script:
        "../scripts/gene_ontology.R"

rule gsea:
    input:
        diffexp_path = lambda wc: "diff_exp/transcripts/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-genic-all.tsv" if wc.annotation=="transcripts" else "diff_exp/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{annotation}-diffexp-results-all.tsv",
        go_mapping_path = lambda wc: config["gene_ontology_mapping_file"] if wc.annotation=="transcripts" else config["differential_expression"]["annotations"][wc.annotation]["gene_ontology_mapping_file"],
    output:
        results = "gene_ontology/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{annotation}-gsea-results.tsv",
        volcano = "gene_ontology/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{annotation}-gsea-volcano.svg",
        dotplot = "gene_ontology/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{annotation}-gsea-dotplot.svg",
    params:
        min_go_group_size = config["gsea"]["min_group_size"],
        max_go_group_size = config["gsea"]["max_group_size"],
        n_permutations = config["gsea"]["n_permutations"],
        fdr_cutoff = config["gsea"]["fdr_cutoff"],
    conda:
        "../envs/fgsea.yaml"
    script:
        "../scripts/fgsea.R"

