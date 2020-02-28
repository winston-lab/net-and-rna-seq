library(tidyverse)
library(magrittr)
library(fgsea)
library(GO.db)

main = function(diffexp_results_path = "Spn1-IAA-v-Spn1-DMSO_rnaseq-spikenorm-verified-coding-genes-diffexp-results-all.tsv",
                go_mapping_path = "test_combined.tsv",
                min_go_group_size = 15,
                max_go_group_size = Inf,
                n_permutations = 1e3,
                results_out="results.tsv",
                volcano_out="volcano.png",
                dotplot_out="dotplot.png",
                fdr_cutoff=0.1){
    goterm_lookup = GOTERM %>%
        as.data.frame() %>%
        magrittr::set_colnames(c("go_id",
                                 "junk",
                                 "term",
                                 "ontology",
                                 "definition",
                                 "synonym",
                                 "secondary")) %>%
        as_tibble() %>%
        distinct(go_id, term, ontology)

    diffexp_results = read_tsv(diffexp_results_path)

    feature_stats = diffexp_results %>%
        pull(log2_foldchange) %>%
        set_names(diffexp_results %>%
                      pull(name))

    go_mapping = read_tsv(go_mapping_path)

    go_mapping_list = split(go_mapping[["feature_name"]],
                            go_mapping[["go_category"]])

    fgsea_results = fgsea(pathways=go_mapping_list,
                          stats=feature_stats,
                          minSize=min_go_group_size,
                          maxSize=max_go_group_size,
                          nperm=n_permutations)

    df = fgsea_results %>%
        as_tibble() %>%
        left_join(goterm_lookup,
                  by=c("pathway"="go_id")) %>%
        transmute(go_id=pathway,
                  enrichment_score=ES,
                  normalized_enrichment_score=NES,
                  log_pval=-log10(pval),
                  log_padj=-log10(padj),
                  n_random_more_extreme=nMoreExtreme,
                  feature_set_size=size,
                  ontology,
                  go_description=term,
                  leading_edge_features=(leadingEdge)) %>%
        unnest(leading_edge_features) %>%
        group_by(go_id, enrichment_score, normalized_enrichment_score,
                 log_pval, log_padj, n_random_more_extreme, feature_set_size,
                 ontology, go_description) %>%
        summarize(leading_edge_features=paste(leading_edge_features, collapse=",")) %>%
        ungroup() %>%
        arrange(desc(log_padj),
                normalized_enrichment_score)

    df %>%
        mutate_if(is.numeric,
                  ~signif(., digits=4)) %>%
        write_tsv(results_out)

    df %<>%
        mutate(label=if_else(is.na(go_description),
                             go_id,
                             go_description),
               label=fct_inorder(label))

    volcano = ggplot(data=df %>%
               mutate(ontology=factor(ontology,
                                      levels=c("BP",
                                               "MF",
                                               "CC"),
                                      labels=c("biological process",
                                               "molecular function",
                                               "cellular compartment"))),
           aes(x=normalized_enrichment_score,
               y=log_padj,
               color=ontology)) +
        geom_vline(xintercept=0,
                   color="gray70",
                   size=0.2) +
        geom_point(alpha=0.8,
                   size=1,
                   shape=16) +
        scale_x_continuous(name="normalized enrichment score") +
        scale_y_continuous(name=expression(-"log"[10] ~ "p"["adj"])) +
        scale_color_viridis_d(end=0.83,
                              na.value="gray50") +
        ggtitle("gene set enrichment") +
        theme_light() +
        theme(panel.grid=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              text=element_text(color="black"),
              axis.text=element_text(color="black"))

    ggsave(volcano_out,
           plot=volcano,
           width=16,
           height=9,
           units="cm")

    dotplot = ggplot(data=df %>%
               filter(log_padj > -log10(fdr_cutoff)),
           aes(y=label,
               x=log_padj,
               size=abs(normalized_enrichment_score),
               color=normalized_enrichment_score)) +
        geom_point() +
        scale_color_distiller(palette="PRGn",
                              limits=c(-max(abs(df[["normalized_enrichment_score"]]), na.rm=TRUE),
                                       max(abs(df[["normalized_enrichment_score"]]), na.rm=TRUE)),

                              breaks=scales::pretty_breaks(4),
                              name="normalized\nenrichment score") +
        scale_size_continuous(name=expression(bgroup("|",
                                                     atop("normalized", "enrichment score"),
                                                     "|"))) +
        scale_x_continuous(name=expression(-"log"[10] ~ "p"["adj"])) +
        theme_light() +
        theme(axis.title.y=element_blank(),
              text=element_text(color="black"),
              axis.text=element_text(color="black"))

    ggsave(dotplot_out,
           plot=dotplot,
           width=16*1.5,
           height=9*1.5,
           units="cm")
}

main(diffexp_results_path = snakemake@input[["diffexp_path"]],
     go_mapping_path = snakemake@input[["go_mapping_path"]],
     min_go_group_size = as.numeric(snakemake@params[["min_go_group_size"]]),
     max_go_group_size = as.numeric(snakemake@params[["max_go_group_size"]]),
     n_permutations = as.numeric(snakemake@params[["n_permutations"]]),
     results_out= snakemake@output[["results"]],
     volcano_out= snakemake@output[["volcano"]],
     dotplot_out= snakemake@output[["dotplot"]],
     fdr_cutoff=snakemake@params[["fdr_cutoff"]])

