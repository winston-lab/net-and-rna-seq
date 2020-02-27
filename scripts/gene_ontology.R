library(tidyverse)
library(magrittr)
library(goseq)
library(ggthemes)
library(cowplot)

main = function(universe_path="Scer_polIItranscripts-adjustedTSS.bed",
                diffexp_path="Spn1-IAA-v-Spn1-DMSO_rnaseq-spikenorm-transcripts-diffexp-results-down.tsv",
                go_anno_path="sgd_go_slim_mapping_2018_2_7.tsv",
                ttype="genic",
                diffexp_direction="downregulated",
                condition, control, results_out, e_combined_path, d_combined_path, e_facet_path, d_facet_path){
    universe = read_tsv(universe_path,
                    col_names=c('chrom', 'start', 'end', 'name', 'score', 'strand'))

    diffexp_results = read_tsv(diffexp_path)

    if (nrow(diffexp_results) > 0){
        if (ttype=="genic"){
            id_column=quo(name)
        } else {
            id_column=quo(transcript_name)
        }

        diffexp_results %<>%
            rename(feature_name=!!id_column) %>%
            group_by(feature_name) %>%
            arrange(desc(log10_padj),
                    .by_group=TRUE) %>%
            dplyr::slice(1) %>%
            ungroup()

        all_genes = universe %>%
            left_join(diffexp_results,
                      by=c('name'='feature_name')) %>%
            dplyr::select(name, log10_padj) %>%
            mutate(log10_padj = if_else(is.na(log10_padj), 0, 1))

        all_genes_vector = all_genes[['log10_padj']]
        names(all_genes_vector) = all_genes[['name']]

        lengths_vector = universe %>%
            transmute(length=end-start) %>%
            pull(length)

        go_anno = read_tsv(go_anno_path,
                           col_names = c('sys_name',
                                         'common_name',
                                         'go_id')) %>%
            filter(!(is.na(go_id))) %>%
            mutate(id = if_else(is.na(common_name),
                                sys_name,
                                common_name)) %>%
            dplyr::select(id, go_id) %>%
            filter(id %in% universe[['name']]) %>%
            as.data.frame()

        pwf = nullp(DEgenes = all_genes_vector,
                    bias.data=lengths_vector,
                    plot.fit=FALSE)

        results = goseq(pwf,
                        gene2cat=go_anno) %>%
            as_tibble() %>%
            mutate(over_represented_pvalue = p.adjust(over_represented_pvalue, method="BH"),
                   under_represented_pvalue = p.adjust(under_represented_pvalue, method="BH")) %>%
            write_tsv(results_out) %>%
            mutate(term = if_else(is.na(term), category, term),
                   ontology = case_when(ontology=="BP" ~ "biological process",
                                        ontology=="MF" ~ "molecular function",
                                        ontology=="CC" ~ "cellular compartment",
                                        is.na(ontology) ~ "other"))

        results_enriched = results %>%
            rename(pval = over_represented_pvalue) %>%
            filter(pval < 0.2) %>%
            arrange(pval) %>%
            mutate(term = fct_rev(fct_inorder(term, ordered=TRUE)))

        results_depleted = results %>%
            rename(pval = under_represented_pvalue) %>%
            filter(pval < 0.2) %>%
            arrange(pval) %>%
            mutate(term = fct_rev(fct_inorder(term, ordered=TRUE)))

        combined_plot = function(df,
                                 enrichment_direction){
            ggplot(data = df, aes(y=-log10(pval), x=term, fill=ontology)) +
                geom_hline(yintercept = 1, linetype="dashed") +
                geom_col(color="black", size=0.3) +
                coord_flip() +
                scale_fill_manual(values = ptol_pal()(4)) +
                ylab(expression(-"log"[10] ~ "p"["adj"])) +
                ggtitle(paste(enrichment_direction, "gene ontology terms"),
                        subtitle = paste(ttype, "transcripts",
                                         diffexp_direction,
                                         "in", condition, "vs.", control)) +
                theme_light() +
                theme(text = element_text(size=12, color="black"),
                      axis.title.y = element_blank(),
                      axis.text.y = element_text(size=10, color="black"),
                      axis.text.x = element_text(size=12, color="black"),
                      legend.text = element_text(size=12))
        }

        e_combined_out = combined_plot(results_enriched, "enriched")
        d_combined_out = combined_plot(results_depleted, "depleted")

        ggplot2::ggsave(e_combined_path, plot=e_combined_out, width=28,
               height=9.5+.24*nrow(results_enriched), units="cm", limitsize=FALSE)
        ggplot2::ggsave(d_combined_path, plot=d_combined_out, width=28,
               height=9.5+.24*nrow(results_depleted), units="cm", limitsize=FALSE)

        filtered_plot = function(df,
                                 ontology_type,
                                 enrichment_direction){
            ggplot(data = df %>%
                       filter(ontology==ontology_type),
               aes(y=-log10(pval),
                   x=term,
                   fill=ontology)) +
            geom_hline(yintercept = 1,
                       linetype="dashed",
                       size=0.5) +
            geom_col(fill="#114477") +
            coord_flip() +
            ylab(expression(-"log"[10] ~ "p"["adj"])) +
            ggtitle(paste0(enrichment_direction, " ", ontology_type,
                          if(ontology_type=="other") {" misc. categories"}
                          else if(substr(ontology_type, nchar(ontology_type), nchar(ontology_type))=="s"){"es"}
                          else {"s"}),
                    subtitle = paste(ttype, "transcripts",
                                     diffexp_direction,
                                     "in", condition, "vs.", control)) +
            theme_light() +
            theme(text = element_text(size=12, color="black"),
                  axis.title.y = element_blank(),
                  axis.text.y = element_text(size=10, color="black"),
                  axis.text.x = element_text(size=12, color="black"),
                  legend.text = element_text(size=12),
                  plot.title = element_text(size=12),
                  plot.subtitle = element_text(size=10))
        }

        e_plotlist = list(filtered_plot(results_enriched, "biological process", "enriched"),
                          filtered_plot(results_enriched, "cellular compartment", "enriched"),
                          filtered_plot(results_enriched, "molecular function", "enriched"),
                          filtered_plot(results_enriched, "other", "enriched"))

        d_plotlist = list(filtered_plot(results_depleted, "biological process", "depleted"),
                          filtered_plot(results_depleted, "cellular compartment", "depleted"),
                          filtered_plot(results_depleted, "molecular function", "depleted"),
                          filtered_plot(results_depleted, "other", "depleted"))

        get_heights = function(df){
            df %>%
                count(ontology) %>%
                complete(ontology = c("biological process", "cellular compartment",
                                      "molecular function", "other"),
                         fill=list(n=0)) %>%
                pull(n) %>%
                return()
        }

        e_heights = get_heights(results_enriched)
        d_heights = get_heights(results_depleted)

        e_facet_out = plot_grid(plotlist = e_plotlist,
                                align="v",
                                axis="l",
                                ncol=1,
                                rel_heights = e_heights/sum(1e-3+e_heights)+0.45)

        d_facet_out = plot_grid(plotlist = d_plotlist,
                                align="v",
                                axis="l",
                                ncol=1,
                                rel_heights = d_heights/sum(1e-3+d_heights)+0.45)

        ggplot2::ggsave(e_facet_path, plot=e_facet_out, width=20,
                        height=8+0.45*sum(e_heights), units="cm", limitsize=FALSE)
        ggplot2::ggsave(d_facet_path, plot=d_facet_out, width=20,
                        height=8+0.45*sum(d_heights), units="cm", limitsize=FALSE)
    } else {
        file.create(results_out,
                    e_combined_path,
                    d_combined_path,
                    e_facet_path,
                    d_facet_path)
    }
}

main(universe_path = snakemake@input[["universe"]],
     diffexp_path = snakemake@input[["diffexp_path"]],
     go_anno_path = snakemake@input[["go_anno_path"]],
     ttype = snakemake@wildcards[["category"]],
     diffexp_direction= snakemake@params[["direction"]],
     condition= snakemake@wildcards[["condition"]],
     control= snakemake@wildcards[["control"]],
     results_out = snakemake@output[["results"]],
     e_combined_path = snakemake@output[["enriched_combined"]],
     d_combined_path = snakemake@output[["depleted_combined"]],
     e_facet_path = snakemake@output[["enriched_facet"]],
     d_facet_path = snakemake@output[["depleted_facet"]])

