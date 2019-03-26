library(tidyverse)
library(ggrepel)
library(ggpmisc)

import = function(df, path, alpha, label_col_id, category){
    df = read_tsv(path) %>%
        distinct(chrom, start, end, transcript_id, score, strand, .keep_all = TRUE) %>%
        select(chrom, start, end, label=label_col_id, score, strand, log2_foldchange, lfc_SE, stat,
               log10_pval, log10_padj, mean_expr, condition_expr, control_expr) %>%
        mutate(sig = if_else(log10_padj > -log10(alpha), TRUE, FALSE),
               category = category) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(in_all, in_genic, in_anti,
                in_conv, in_div, in_inter,
                assay, condition, control, lfc, alpha,
                out_ma, out_volcano, out_volcano_free, out_mosaic, out_summary_table){
    df = tibble() %>%
        import(in_all, alpha=alpha, label_col_id="name", category="all") %>%
        import(in_genic, alpha=alpha, label_col_id="name", category="genic") %>%
        import(in_anti, alpha=alpha, label_col_id="transcript_name", category="antisense") %>%
        import(in_conv, alpha=alpha, label_col_id="transcript_name", category="convergent") %>%
        import(in_div, alpha=alpha, label_col_id="transcript_name", category="divergent") %>%
        import(in_inter, alpha=alpha, label_col_id="name", category="intergenic") %>%
        mutate(category = fct_inorder(category, ordered=TRUE))

    min_x = quantile(df[["mean_expr"]], .2)

    maplot = ggplot(data = df, aes(x=mean_expr, y=log2_foldchange)) +
                geom_hline(yintercept = 0, linetype="dashed") +
                geom_point(aes(color=sig), shape=1, alpha=0.4, size=.4) +
                scale_color_manual(values = c("grey40", "#440154FF"), guide=FALSE) +
                stat_dens2d_filter(data = df %>% filter(sig & mean_expr>min_x & log2_foldchange > 0),
                                   geom = "text_repel",
                                   aes(label = label),
                                   keep.number = 5,
                                   point.padding = unit(0.1, "lines"),
                                   box.padding = unit(0.1, "lines"),
                                   nudge_y = .3,
                                   size=1.5) +
                stat_dens2d_filter(data = df %>% filter(sig & mean_expr>min_x & log2_foldchange < 0),
                                   geom = "text_repel",
                                   aes(label = label),
                                   keep.number = 5,
                                   point.padding = unit(0.1, "lines"),
                                   box.padding = unit(0.1, "lines"),
                                   nudge_y = -.3,
                                   size=1.5) +
                scale_x_log10(name="mean expression level", limits=c(1, NA),
                              expand=c(0,0)) +
                scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
                ylab(bquote(bold(log[2]~frac(.(condition), .(control))))) +
                facet_wrap(~category) +
                ggtitle(paste(assay, "MA plots:", condition, "vs.", control),
                        subtitle = bquote("DESeq2: |" ~ log[2] ~ "fold-change" ~ "| >" ~ log[2](.(lfc)) ~ "@ FDR" ~ .(alpha)))  +
                theme_bw() +
                theme(text = element_text(size=12, face="bold"),
                      axis.text = element_text(size=10, color="black"),
                      axis.title.y = element_text(angle=0, vjust=0.5),
                      strip.background = element_blank(),
                      strip.text = element_text(size=12, face="bold", color="black"))

    ggsave(out_ma, maplot, height=16, width=22, units="cm")

    volcano = ggplot(data = df, aes(x=log2_foldchange, y=log10_padj))+
                geom_point(aes(color=sig), shape=1, alpha=0.4, size=.4) +
                scale_color_manual(values = c("grey40", "#440154FF"), guide=FALSE) +
                stat_dens2d_filter(data = df %>% filter(sig & log2_foldchange < 0),
                                   geom = "text_repel",
                                   aes(label = label),
                                   keep.number = 5,
                                   point.padding = unit(0.1, "lines"),
                                   box.padding = unit(0.1, "lines"),
                                   nudge_x = -0.3,
                                   size=1.5) +
                stat_dens2d_filter(data = df %>% filter(sig & log2_foldchange > 0),
                                   geom = "text_repel",
                                   aes(label = label),
                                   keep.number = 5,
                                   point.padding = unit(0.1, "lines"),
                                   box.padding = unit(0.1, "lines"),
                                   nudge_x = 0.3,
                                   size=1.5) +
                ylab(expression(bold(-log[10] ~ p[adj]))) +
                xlab(bquote(bold(log[2]~frac(condition,control)))) +
                ggtitle(paste(assay, "volcano plots:", condition, "vs.", control),
                        subtitle = bquote("DESeq2: |" ~ log[2] ~ "fold-change" ~ "| >" ~ log[2](.(lfc)) ~ "@ FDR" ~ .(alpha)))  +
                theme_bw() +
                theme(text = element_text(size=12, face="bold"),
                      axis.text = element_text(size=10, color="black"),
                      axis.title.y = element_text(angle=0, vjust=0.5),
                      strip.background = element_blank(),
                      strip.text = element_text(size=12, face="bold", color="black"))

    ggsave(out_volcano, volcano + facet_wrap(~category) , height=16, width=20, units="cm")
    ggsave(out_volcano_free, volcano + facet_wrap(~category, scales="free_y") , height=16, width=20, units="cm")

    count_df = df %>%
        filter(category != "all") %>%
        mutate(change = if_else(sig, if_else(log2_foldchange > 0, "up", "down"), "unchanged")) %>%
        count(category, change) %>%
        write_tsv(out_summary_table) %>%
        filter(! is.na(change)) %>%
        mutate(xmax = cumsum(n), xmin=cumsum(n)-n) %>%
        group_by(category) %>%
        mutate(ymax=cumsum(n), ymin=cumsum(n)-n,
               xmax=max(xmax), xmin=min(xmin)) %>%
        mutate_at(vars(ymin, ymax), ~(./max(ymax))) %>%
        ungroup() %>%
        mutate_at(vars(xmin, xmax), ~(./max(xmax))) %>%
        mutate(x=(xmin+xmax)/2,
               y=(ymin+ymax)/2)

    mosaic = ggplot() +
        geom_rect(data = count_df,
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=change),
                  color = "white", size=1.5) +
        geom_text(data = count_df,
                  aes(x=x, y=y, label=n),
                  size=12/75*25.4, color="black", fontface="bold") +
        geom_text(data = count_df %>% group_by(category) %>% summarise(x=first(x)),
                  aes(x=x, label=category), y=1.04, angle=30, hjust=.1,
                  size=12/75*25.4, fontface="bold") +
        scale_fill_brewer(palette = "Set1",
                          guide=guide_legend(reverse=TRUE)) +
        scale_x_continuous(limits = c(NA, 1.08), expand=c(0,0)) +
        scale_y_continuous(limits = c(NA, 1.15)) +
        theme_void() +
        ggtitle(paste(assay, "differential expression:", condition, "vs.", control),
                subtitle = bquote("DESeq2: |" ~ log[2] ~ "fold-change" ~ "| >" ~ log[2](.(lfc)) ~ "@ FDR" ~ .(alpha)))  +
        theme(legend.text = element_text(size=12, face="bold", color="black"),
              legend.title = element_blank(),
              legend.position = c(.96, .5),
              legend.justification = "left",
              legend.key.size = unit(1, "cm"),
              plot.margin = unit(c(.5, 3.25, 0, 0.25), "cm"),
              plot.title = element_text(size=12, face="bold", color="black"),
              plot.subtitle = element_text(size=10))

    ggsave(out_mosaic, mosaic, height=10, width=20, units="cm")
}

main(in_all = snakemake@input[["total"]],
     in_genic = snakemake@input[["genic"]],
     in_anti = snakemake@input[["antisense"]],
     in_conv = snakemake@input[["convergent"]],
     in_div = snakemake@input[["divergent"]],
     in_inter = snakemake@input[["intergenic"]],
     assay = ifelse(snakemake@wildcards[["assay"]]=="rnaseq",
                    "RNA-seq", "NET-seq"),
     condition = snakemake@wildcards[["condition"]],
     control = snakemake@wildcards[["control"]],
     lfc = snakemake@params[["lfc"]],
     alpha = snakemake@params[["alpha"]],
     out_ma = snakemake@output[["maplot"]],
     out_volcano = snakemake@output[["volcano"]],
     out_volcano_free = snakemake@output[["volcano_free"]],
     out_mosaic = snakemake@output[["mosaic"]],
     out_summary_table = snakemake@output[["summary_table"]])

