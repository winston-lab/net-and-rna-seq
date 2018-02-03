library(tidyverse)
library(forcats)
library(ggthemes)

main = function(in_numerator, in_denominator, samplelist, num_size, den_size, pcount=1e-3,
                rlabel, nlabel, dlabel, alabel, violin_out, ecdf_out){
    df = read_tsv(in_numerator,
                     col_names=c('group', 'sample', 'index', 'position', 'signal')) %>%
        filter(sample %in% samplelist) %>% 
        select(-position) %>% 
        left_join(read_tsv(in_denominator,
                           col_names=c('group', 'sample', 'index', 'position', 'signal')) %>%
                      select(-position) %>% filter(sample %in% samplelist),
                  by=c("group", "sample", "index"), suffix=c("_num", "_den")) %>% 
        mutate(group = fct_inorder(group, ordered=TRUE),
               sample = fct_inorder(sample, ordered=TRUE),
               signal_num = signal_num/num_size,
               signal_den = signal_den/den_size) %>% 
        mutate(ratio = log2((signal_num+pcount)/(signal_den+pcount)))
    
    n_anno = n_distinct(df[["index"]])
    n_samples = length(samplelist)
    
    theme_default = theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              plot.title = element_text(size=12),
              plot.subtitle = element_text(size=10, face="plain"),
              legend.text = element_text(size=12))
    
    violin_plot = ggplot(data = df, aes(x=sample, y=ratio, fill=group)) +
        geom_hline(yintercept = 0, color="grey65", size=1) +
        geom_violin(draw_quantiles = c(0.5), color="black") +
        scale_fill_ptol(guide=guide_legend(label.position="top"), name=NULL) +
        scale_color_ptol() +
        ylab(bquote(bold(log[2] ~ frac(.(nlabel),.(dlabel))))) +
        ggtitle(paste0("NET-seq \"", rlabel, "\""),
                subtitle = paste0(n_anno," ",alabel, " >", (num_size+den_size)/1000, "kb")) +
        theme_default +
        theme(legend.position="top",
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle=30, hjust=0.9),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))
    
    ecdf_plot = ggplot(data = df, aes(x=ratio, group=sample, color=group)) +
        geom_vline(xintercept = 0, color="grey65", size=0.4) +
        stat_ecdf(geom="step", size=1) +
        scale_color_ptol(name=NULL) +
        xlab(bquote(bold(log[2] ~ frac(.(nlabel),.(dlabel))))) +
        ylab("cumulative probability") +
        ggtitle(paste0("eCDF of NET-seq \"", rlabel, "\""),
                subtitle = paste0(n_anno," ",alabel, " >", (num_size+den_size)/1000, "kb")) +
        theme_default +
        theme(legend.position = c(0.025,0.95),
              legend.justification = c(0,1),
              legend.background = element_blank(),
              axis.title.y = element_text(size=10, face="plain"),
              axis.text.y = element_text(size=10, face="plain"),
              legend.key = element_blank())
    
    ggsave(violin_out, plot=violin_plot, width=6+2.5*n_samples, height=12, units="cm")
    ggsave(ecdf_out, plot=ecdf_plot, width=14, height=10, units="cm")
    
    # group_df = df %>% group_by(group, index) %>% 
    #     summarise(n = mean(signal_num),
    #               d = mean(signal_den)) %>% 
    #     mutate(r = log2((n+pcount)/(d+pcount))) %>% 
    #     select(-c(n,d))
    # 
    # ggplot(data = group_df, aes(x=group, y=r, group=index)) +
    #     geom_line(size=0.1, alpha=0.2)
}

main(in_numerator = snakemake@input[["numerator"]],
    in_denominator = snakemake@input[["denominator"]],
    num_size = snakemake@params[["num_size"]],
    den_size = snakemake@params[["den_size"]],
    pcount = snakemake@params[["pcount"]],
    samplelist = snakemake@params[["samplelist"]],
    rlabel = snakemake@params[["ratio_label"]],
    nlabel = snakemake@params[["num_label"]],
    dlabel = snakemake@params[["den_label"]],
    alabel = snakemake@params[["annotation_label"]],
    violin_out = snakemake@output[["violin"]],
    ecdf_out = snakemake@output[["ecdf"]])
