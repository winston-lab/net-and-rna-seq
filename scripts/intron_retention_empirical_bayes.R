library(tidyverse)
library(magrittr)
library(VGAM)

theme_default = theme_light() +
    theme(text = element_text(size=10, color="black"),
          axis.text = element_text(size=10, color="black"),
          axis.title.y = element_text(angle=0, vjust=0.5),
          panel.grid.major.x = element_line(color="grey95"),
          strip.background = element_blank(),
          strip.text = element_text(size=10, color="black"))

main = function(input_path="input.tsv",
                n_trials=1e4,
                cred_int_level=0.95,
                fdr_cutoff=0.05,
                control_group_id = "Spn1-DMSO",
                query_group_id = "Spn1-IAA",
                credible_intervals_out = "credible_intervals.tsv",
                qcplots_out = "qcplots.png",
                results_out = "results.tsv",
                credible_intervals_plot_out = "credible_intervals.png",
                volcano_plot_out = "volcano.png"){

    df = read_tsv(input_path) %>%
        filter(group_id %in% c(control_group_id, query_group_id)) %>%
        group_by(chrom, start, end, name, strand, group_id) %>%
        summarize_at(vars(spliced, junction_5, junction_3, intronic, ambiguous),
                     ~sum(.)) %>%
        mutate(unspliced = junction_5 + junction_3 + intronic,
               total = unspliced + spliced) %>%
        group_by(chrom, start, end, name, strand) %>%
        mutate(never_spliced = sum(spliced) == 0)

    dispersion_plot = ggplot(data=df,
           aes(x=spliced + unspliced,
               y=unspliced / (spliced + unspliced),
               color=never_spliced)) +
        geom_point(alpha=0.85) +
        scale_y_continuous(name=expression(textstyle(frac("unspliced",
                                                          "spliced + unspliced")))) +
        scale_color_brewer(palette="Set1",
                           name="never spliced",
                           direction=-1) +
        theme_default

     intronic_v_junction_plot = ggplot(data=df %>%
                                           filter(!never_spliced),
                                       aes(x=junction_5 + junction_3,
                                           y=intronic)) +
         geom_abline(slope=1,
                     intercept=0,
                     color="gray70",
                     size=0.2) +
         geom_point(alpha=0.8) +
         scale_x_log10(name="junction alignments") +
         scale_y_log10(name="intronic\nalignments") +
         theme_default

    nll = function(alpha, beta){
        -sum(dbetabinom.ab(df %>%
                               filter(!never_spliced) %>%
                               pull(unspliced),
                           df %>%
                               filter(!never_spliced) %>%
                               pull(total),
                           alpha, beta,
                           log=TRUE))
    }

    m = mle(nll,
            start=list(alpha = 1,
                       beta = 1),
            lower=list(alpha=0, beta=0),
            method="L-BFGS-B")
    alpha_0 = coef(m)[["alpha"]]
    beta_0 = coef(m)[["beta"]]

    beta_binom_fit = ggplot(data = df,
           aes(x=unspliced/total)) +
        geom_histogram(binwidth=0.01,
                       aes(y=..density..)) +
        # geom_density(bw=0.01) +
        stat_function(fun = function(x) dbeta(x, shape1=alpha_0, shape2=beta_0),
                      color="red", size=1) +
        scale_y_continuous(expand=c(0,0),
                           breaks = scales::pretty_breaks(3)) +
        scale_x_continuous(name = expression(textstyle(frac("unspliced", "unspliced + spliced"))),
                           expand=c(0,0)) +
        theme_default

    df %<>%
        mutate(alpha_1 = alpha_0 + unspliced,
               beta_1 = beta_0 + spliced,
               ir_est = alpha_1/(alpha_1 + beta_1),
               ir_cred_int_low = qbeta((1-cred_int_level)/2, alpha_1, beta_1),
               ir_cred_int_high = qbeta((1+cred_int_level)/2, alpha_1, beta_1))

    df %>%
        mutate_if(is.numeric,
                  ~signif(., digits=4)) %>%
        write_tsv(credible_intervals_out)


    shrinkage_plot = ggplot(data = df,
           aes(x = unspliced/total,
               y = ir_est,
               color = total)) +
        geom_hline(yintercept = alpha_0/(alpha_0 + beta_0),
                   linetype = "dashed") +
        geom_abline(slope=1,
                    intercept=0) +
        geom_jitter(alpha = 0.6,
                   size = 0.5,
                   width = 0.01,
                   height = 0.01) +
        scale_color_viridis_c(limits = c(NA, quantile(df[["total"]], 0.9)),
                              oob = scales::squish,
                              name="unspliced + spliced",
                              guide=guide_colorbar(title.position="top",
                                                   title.hjust=0.5,
                                                   barwidth=8,
                                                   barheight=0.5)) +
        xlab(expression(bgroup("(", textstyle(frac("unspliced", "unspliced + spliced")), ")")["raw"])) +
        ylab(expression(bgroup("(", textstyle(frac("unspliced", "unspliced + spliced")), ")")["shrunken"])) +
        theme_default +
        theme(legend.position = c(0.05,0.95),
              legend.justification = c(0,1),
              legend.direction="horizontal")

    qc_plots = gridExtra::arrangeGrob(intronic_v_junction_plot,
                                      dispersion_plot,
                                      beta_binom_fit,
                                      shrinkage_plot,
                                      ncol=2,
                                      widths=c(0.8, 1))
    ggsave(qcplots_out,
           plot=qc_plots,
           width=16*2,
           height=9*2,
           units="cm")

    df_proptest = df %>%
        ungroup() %>%
        group_by_all() %>%
        mutate(rbeta = list(rbeta(n_trials, alpha_1, beta_1))) %>%
        group_by(chrom, start, end, name, strand) %>%
        do(tibble(rbeta_condition = .$rbeta[.$group_id == query_group_id],
                  rbeta_controls = list(.$rbeta[.$group_id != query_group_id]))) %>%
        mutate(ir_est_condition = rbeta_condition %>%
                   unlist() %>%
                   median(),
               ir_low_condition = rbeta_condition %>%
                   unlist() %>%
                   quantile(1-cred_int_level),
               ir_high_condition = rbeta_condition %>%
                   unlist() %>%
                   quantile(cred_int_level),
               ir_est_controls = rbeta_controls %>%
                   unlist() %>%
                   median(),
               ir_low_controls = rbeta_controls %>%
                   unlist() %>%
                   quantile(1-cred_int_level),
               ir_high_controls = rbeta_controls %>%
                   unlist() %>%
                   quantile(cred_int_level),
               rbeta_control_max = rbeta_controls %>%
                   unlist() %>%
                   matrix(ncol=n_trials, byrow=TRUE) %>%
                   apply(2, max) %>%
                   list(),
               rbeta_control_min = rbeta_controls %>%
                   unlist() %>%
                   matrix(ncol=n_trials, byrow=TRUE) %>%
                   apply(2, min) %>%
                   list()) %>%
        group_by(chrom, start, end, name, strand,
                 ir_est_condition, ir_low_condition, ir_high_condition,
                 ir_est_controls, ir_low_controls, ir_high_controls) %>%
        do(tibble(prob_increase = mean(.$rbeta_condition[[1]] > .$rbeta_control_max[[1]]),
                  prob_decrease = mean(.$rbeta_condition[[1]] < .$rbeta_control_min[[1]]))) %>%
        ungroup() %>%
        arrange(desc(prob_increase)) %>%
        mutate(qvalue_increase = cummean(1-prob_increase)) %>%
        arrange(desc(prob_decrease)) %>%
        mutate(qvalue_decrease = cummean(1-prob_decrease),
               category = case_when(qvalue_decrease < fdr_cutoff ~ "decreased",
                                    qvalue_increase < fdr_cutoff ~ "increased",
                                    TRUE ~ "not significant"),
               category = ordered(category, levels=c("decreased", "not significant", "increased"))) %>%
        arrange(pmin(qvalue_increase, qvalue_decrease))

    df_proptest %>%
        mutate_if(is.numeric,
                  ~signif(., digits = 4)) %>%
        write_tsv(results_out)

    intron_order = c(df_proptest %>%
                         filter(category=="decreased") %>%
                         arrange(qvalue_decrease,
                                 ir_est_controls - ir_est_condition) %>%
                         pull(name),
                     df_proptest %>%
                         filter(category=="not significant" &
                                   ir_est_condition < ir_est_controls) %>%
                         arrange(qvalue_decrease,
                                 ir_est_controls - ir_est_condition) %>%
                         pull(name),
                     df_proptest %>%
                         filter(category=="not significant" &
                                   ir_est_condition >= ir_est_controls) %>%
                         arrange(desc(qvalue_increase),
                                 ir_est_condition - ir_est_controls) %>%
                         pull(name),
                     df_proptest %>%
                         filter(category=="increased") %>%
                         arrange(desc(qvalue_increase),
                                 ir_est_condition - ir_est_controls) %>%
                         pull(name)
    )

    df_proptest %<>%
        mutate(name = ordered(name, levels=intron_order))

    credible_intervals_plot = ggplot(data = df_proptest) +
        geom_segment(aes(y=ir_low_controls,
                         yend=ir_high_controls,
                         x=name, xend=name),
                     color="black",
                     size=0.8,
                     # position=position_nudge(x=0.05),
                     alpha=0.8) +
        geom_point(aes(y=ir_est_controls,
                       x=name),
                   color="black",
                   size=0.4) +
        geom_segment(aes(y=ir_low_condition,
                         yend=ir_high_condition,
                         x=name,
                         xend=name),
                     color="red",
                     size=0.3,
                     # position=position_nudge(x=-0.05),
                     alpha=0.8) +
        geom_point(aes(y=ir_est_condition,
                       x=name),
                   color="red",
                   size=0.4) +
        scale_y_continuous(name = expression(textstyle(frac("unspliced","unspliced + spliced"))),
                           expand = c(0,0.01)) +
        xlab("introns") +
        facet_grid(.~category, scales="free_x", space="free_x") +
        theme_default +
        theme(axis.text.x = element_text(size=2, color="black", angle=90, hjust=1, vjust=0.5))

    ggsave(credible_intervals_plot_out,
           plot=credible_intervals_plot,
           width=16*1.5,
           height=9*1.5, units="cm")

    volcano = ggplot(data=df_proptest,
           aes(x=ir_est_condition - ir_est_controls,
               y=ifelse(ir_est_condition >= ir_est_controls,
                        -log10(qvalue_increase),
                        -log10(qvalue_decrease)),
               color=(category=="not significant")))  +
        geom_vline(xintercept=0,
                   color="grey70",
                   size=0.2) +
        geom_point() +
        xlab(expression("IR"["condition"] - "IR"["controls"])) +
        scale_y_continuous(name = expression("-log"[10]("q-value"))) +
        theme_default +
        theme(legend.position="none",
              axis.title=element_text(size=16),
              axis.title.y=element_text(angle=90, vjust=1),
              axis.text=element_text(size=16),
              plot.margin=margin(r=8, unit="pt"),
              panel.grid=element_blank())

    ggsave(volcano_plot_out,
           plot=volcano,
           width=16,
           height=9,
           units="cm")
}

main(input_path=snakemake@input[["data"]],
     n_trials=as.numeric(snakemake@params[["n_trials"]]),
     cred_int_level=snakemake@params[["credible_interval_level"]],
     fdr_cutoff=snakemake@params[["fdr_cutoff"]],
     control_group_id = snakemake@wildcards[["control"]],
     query_group_id = snakemake@wildcards[["condition"]],
     credible_intervals_out = snakemake@output[["credible_intervals"]],
     qcplots_out = snakemake@output[["qcplots"]],
     results_out = snakemake@output[["results"]],
     credible_intervals_plot_out = snakemake@output[["credible_intervals_plot"]],
     volcano_plot_out = snakemake@output[["volcano_plot"]])

