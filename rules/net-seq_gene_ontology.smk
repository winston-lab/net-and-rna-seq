#!/usr/bin/env python

localrules:
    gene_ontology_called_transcripts

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

# rule gsea:
#     input:
#         diffexp_path = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-{category}-all.tsv",
#         go_mapping_path = config["gene_ontology"]["gene_ontology_mapping_file"],
#     output:
#         results = "gene_ontology/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-{category}-{direction}-gene-ontology-results.tsv",
#     params:
#         min_go_group_size = config["gene_ontology"]["gsea_min_group_size"],
#         max_go_group_size = config["gene_ontology"]["gsea_max_group_size"],
#         n_permutations = config["gene_ontology"]["gsea_n_permutations"],
