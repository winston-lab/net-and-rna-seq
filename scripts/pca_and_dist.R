library(tidyverse)
library(forcats)
library(stringr)
# library(broom)
# library(ggfortify)
library(seriation)
library(viridis)
library(ggdendro)
library(ggrepel)
library(ggthemes)
library(gridExtra)
library(gmodels)

main = function(intable, samplelist, grouplist, binsize,
                scree_out, pca_out, cluster_out){
    df = read_tsv(intable) %>% 
        gather(key=sample, value=signal, -name) %>% 
        filter(sample %in% samplelist) %>%
        group_by(name) %>% 
        filter(any(signal>0)) %>%
        mutate(sample = fct_inorder(sample, ordered=TRUE)) %>% 
        spread(key=name, value=signal)
    
    n_samples = length(samplelist)
    
    pcobject = df %>% select(-sample) %>% fast.prcomp(center=FALSE, scale=FALSE)
    
    data = bind_cols(df["sample"], pcobject[["x"]] %>% as_tibble(), tibble(group=grouplist)) %>% 
        mutate_if(is.numeric, funs(.8*./max(abs(.)))) %>% 
        mutate(group = fct_inorder(group, ordered=TRUE))
    
    rotations = pcobject[["rotation"]] %>% as_tibble(rownames="location") %>% 
        select(1:3) %>% 
        mutate(norm = sqrt(PC1^2+PC2^2)) %>% 
        arrange(desc(norm))
    
    pct_var = pcobject[["sdev"]] %>% as_tibble() %>% rowid_to_column(var="pc") %>% 
        mutate(variance = value^2,
               pct_variance = variance/sum(variance),
               cum_variance = cumsum(pct_variance))
    
    pct_var_plot = ggplot(data = pct_var, aes(x=pc, y=pct_variance)) +
        geom_line() +
        geom_point() +
        xlab("principal component") +
        ylab("variance explained") +
        theme_light() +
        theme(text = element_text(size=12),
              axis.text = element_text(size=12, color="black"))
    
    cum_var_plot = ggplot(data = pct_var, aes(x=pc, y=cum_variance)) +
        geom_hline(yintercept = 1, linetype="dashed", color="grey75") +
        geom_line() +
        geom_point() +
        xlab("principal component") +
        ylab("cumulative variance explained") +
        theme_light() +
        theme(text = element_text(size=12),
              axis.text = element_text(size=12, color="black"))
    
    ggsave(scree_out, plot=arrangeGrob(pct_var_plot, cum_var_plot, ncol=1),
           width=16, height=16, units="cm")
    
    pca_plot = ggplot() +
        # geom_segment(data = rotations %>% filter(norm>.15*max(norm)),
        #              aes(xend=0, yend=0, x=PC1, y=PC2)) +
        # geom_text(data = rotations %>% filter(norm>.15*max(norm)),
        #           aes(x=PC1, y=PC2, label=location, angle=atan(PC2/PC1)*360/(2*pi)),
        #           size=1) +
        geom_text_repel(data = data, aes(x=PC1, y=PC2, label=sample), 
                        size=10/72*25.4,
                        fontface="bold",
                        box.padding = unit(0.4, "lines"),
                        point.padding = unit(1e-6, "lines")) +
        geom_point(data = data, aes(x=PC1, y=PC2, color=group)) +
        xlab(paste0("PC1: ", pct_var %>% slice(1) %>% pull(pct_variance) %>% "*"(100) %>% round(1), "% variance explained")) +
        ylab(paste0("PC2: ", pct_var %>% slice(2) %>% pull(pct_variance) %>% "*"(100) %>% round(1), "% variance explained")) +
        scale_color_ptol(guide=FALSE) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"))
    
    pca_loadings_plot = ggplot() +
        geom_segment(data = rotations %>% filter(norm>.5*max(norm)),
                     aes(xend=0, yend=0, x=PC1, y=PC2)) +
        geom_text_repel(data = data, aes(x=PC1, y=PC2, label=sample),
                        size=8/72*25.4) +
        geom_point(data = data, aes(x=PC1, y=PC2, color=group)) +
        geom_label_repel(data = rotations %>% filter(norm>.5*max(norm)) %>% 
                      mutate(location = str_replace_all(location, c("-minus-"="-:", "-plus-"="+:"))),
                  aes(x=PC1, y=PC2, label=location), 
                  point.padding=unit(0, "lines"), label.padding=unit(0.15, "lines"), box.padding=unit(0.05, "lines"),
                  size=8/72*25.4,
                  alpha=0.9) +
        xlab(paste0("PC1: ", pct_var %>% slice(1) %>% pull(pct_variance) %>% "*"(100) %>% round(1), "% variance explained")) +
        ylab(paste0("PC2: ", pct_var %>% slice(2) %>% pull(pct_variance) %>% "*"(100) %>% round(1), "% variance explained")) +
        scale_color_ptol(guide=FALSE) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"))
    
    ggsave(pca_out, plot=arrangeGrob(pca_plot, pca_loadings_plot, ncol = 1), width=16, height=16, units="cm")
    #####
    ### hierarchical clustering on Euclidean distances
    distances = df %>% remove_rownames() %>% 
        column_to_rownames(var="sample") %>% dist()
    
    dend = distances %>% reorder(hclust(distances), .) %>% dendro_data()
        
    dist_df = distances %>% 
        as.matrix() %>% as_tibble(rownames="sample") %>% 
        left_join(dend[["labels"]], by=c("sample"="label")) %>% 
        arrange(x) %>% 
        mutate(sample = fct_inorder(sample, ordered=TRUE)) %>% 
        select(-c(x,y)) %>%  
        gather(key=sample2, value=distance, -sample) %>% 
        mutate(sample2 = ordered(sample2, levels=levels(sample)))
    
    dist_plot = ggplot() +
        geom_raster(data = dist_df, aes(x=sample, y=sample2, fill=distance)) +
        geom_segment(data = dend[["segments"]] %>%
                         rowid_to_column(var="index") %>% 
                         gather(end, value, c(y,yend)) %>% 
                         mutate(value = scales::rescale(value, to=c(0,n_samples/2))) %>% 
                         spread(end, value),
                     aes(yend=x, xend=-(2+y), y=xend, x=-(2+yend)),
                     lineend="square") +
        geom_label(data = dend[["labels"]],
                  aes(x=y, y=x, label=label),
                  size=12/72*25.4, fontface="bold", hjust=1,
                  label.size=0) +
        scale_fill_viridis(option="inferno",
                           name=paste0("Euclidean distance,\nNET-seq signal in ", binsize, "nt bins"),
                           breaks = scales::pretty_breaks(n=2),
                           guide=guide_colorbar(barwidth=4+.5*n_samples, barheight=1, title.vjust=0.8)) +
        theme_minimal() +
        theme(panel.grid.major.y = element_blank(),
              text = element_text(size=12, color="black", face="bold"),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=12, color="black", angle=30, hjust=0.9),
              axis.title = element_blank(),
              legend.position="top")
    
    ggsave(cluster_out, plot=dist_plot, width=2.5*n_samples, height=2+1.5*n_samples, units="cm")
}

main(intable = snakemake@input[[1]],
     samplelist = snakemake@params[["samplelist"]],
     grouplist = snakemake@params[["grouplist"]],
     binsize= snakemake@wildcards[["windowsize"]],
     scree_out = snakemake@output[["scree"]],
     pca_out = snakemake@output[["pca"]],
     cluster_out = snakemake@output[["dist"]])

