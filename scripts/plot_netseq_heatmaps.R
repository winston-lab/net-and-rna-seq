#!/usr/bin/env Rscript
library(argparse)
library(tidyverse)
library(forcats)
library(viridis)
library(dendsort)

parser = ArgumentParser()
parser$add_argument('-i', dest='input', type='character')
parser$add_argument('-s', dest='samplelist', type='character', nargs='+')
parser$add_argument('-t', dest='type', type='character')
parser$add_argument('-u', dest='upstream', type='integer')
parser$add_argument('-d', dest='downstream', type='integer')
parser$add_argument('-c', dest='pct_cutoff', type='double')
parser$add_argument('-z', dest='cluster', type='character')
parser$add_argument('-k', dest='k', type='integer')
parser$add_argument('-r', dest='refptlabel', type='character', nargs='+')
parser$add_argument('-f', dest='strand', type='character')
parser$add_argument('-l', dest='scaled_length', type='integer')
parser$add_argument('-e', dest='endlabel', type='character', nargs='+')
parser$add_argument('-y', dest='ylabel', type='character', nargs='+')
parser$add_argument('-m', dest='cmap', type='character')
parser$add_argument('-o', dest='samples_out', type='character')
parser$add_argument('-p', dest='group_out', type='character')

args = parser$parse_args()

format_xaxis = function(refptlabel, upstream, dnstream){
    function(x){
        if (first(upstream)>500 | first(dnstream)>500){
            return(if_else(x==0, refptlabel, as.character(x)))    
        }    
        else {
            return(if_else(x==0, refptlabel, as.character(x*1000)))
        }
    }
}

main = function(intable, samplelist, type, upstream, dnstream, pct_cutoff,
                cluster, k, refptlabel, strand, scaled_length, endlabel, ylabel, cmap, samples_out, group_out){
    hmap = function(df, refptlabel="refpt", scaled_length=0, endlabel="endpt"){
        #pseudocount for log-transform
        pcount = .01

        heatmap_base = ggplot(data = df) +
            geom_raster(aes(x=position, y=index, fill=log2(cpm+pcount))) +
            scale_y_reverse(name=paste(nindices, ylabel), expand=c(0.02, 0)) +
            scale_fill_viridis(option = cmap,
                               na.value="FFFFFF00",
                               name=bquote(bold(log[2] ~ .(strand) ~ "NET-seq signal")),
                               guide=guide_colorbar(title.position="top",
                                                    barwidth=15, barheight=1, title.hjust=0.5)) +
            theme_minimal() +
            theme(text = element_text(size=12, face="bold", color="black"),
                  legend.position = "top",
                  legend.title = element_text(size=12, face="bold", color="black"),
                  legend.text = element_text(size=10, face="plain"),
                  legend.margin = margin(0,0,0,0),
                  legend.box.margin = margin(0,0,0,0),
                  strip.text = element_text(size=12, face="bold", color="black"),
                  axis.text.y = element_blank(),
                  axis.text.x = element_text(size=12, face="bold", color="black", margin = unit(c(0,0,0,0),"cm")),
                  panel.grid.major.x = element_line(color="black"),
                  panel.grid.minor.x = element_line(color="black"),
                  panel.grid.major.y = element_line(color="black"),
                  panel.grid.minor.y = element_blank(),
                  panel.spacing.x = unit(.5, "cm"))
        if(type=="absolute"){
            heatmap_base = heatmap_base +
                scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                                   labels=format_xaxis(refptlabel=refptlabel,
                                                       upstream=upstream,
                                                       dnstream=dnstream),
                                   name=paste("distance from", refptlabel,
                                              if_else(upstream>500 | dnstream>500, "(kb)", "(nt)")),
                                   expand=c(0.05,0))
        }
        else {
            heatmap_base = heatmap_base +
                scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                                   labels=c(refptlabel, "", endlabel),
                                   name="scaled distance",
                                   expand=c(0.05,0))
            
        }
        return(heatmap_base)
    }
    
    raw = read_tsv(intable, col_names=c("group", "sample", "index", "position","cpm")) %>%
        filter(sample %in% samplelist & !is.na(cpm)) %>% 
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE)))
        
    nindices = max(raw$index, na.rm=TRUE)
    nsamples = length(samplelist)
    ngroups = length(fct_unique(raw$group))

    #clustering
    if (cluster=="yes"){
        #first k-means clustering on NOTE: unscaled but log transformed data
        pcount = 0.1
        rr = raw %>% mutate(cpm = log2(cpm+pcount)) %>% select(-group) %>% unite(cid, c(sample, position), sep="~") %>%
                spread(cid, cpm, fill=0) %>% select(-index)
        clust = kmeans(rr, centers=k)
       
        #then hierarchical clustering on the k-means centers and sorting
        #of the resulting dendrogram
        centerclust = clust$centers %>% dist() %>% hclust() %>% dendsort(isReverse=TRUE)
        
        reorder = clust$cluster %>% as_tibble() %>% rename(cluster=value) %>%
            mutate(cluster = factor(cluster, levels=centerclust$order, ordered=TRUE),
                   og_index=row_number()) %>%
            arrange(cluster,og_index) %>%
            mutate(new_index=row_number())
        raw = raw %>% left_join(reorder, by=c("index"="og_index")) %>%
            select(-index) %>% rename(index=new_index)
    } 
    
    #get replicate info for sample facetting
    repl_df = raw %>% select(group, sample) %>% distinct() %>% group_by(group) %>%
        mutate(replicate=row_number()) %>% ungroup() %>% select(-group)
    
    #plot heatmap facetted by sample and group
    df_sample = raw %>% left_join(repl_df, by="sample")
    sample_cutoff = quantile(df_sample$cpm, probs=pct_cutoff, na.rm=TRUE)
    df_sample = df_sample %>%
        mutate(cpm = pmin(sample_cutoff, cpm))
    
    hmap_sample = hmap(df=df_sample, refptlabel=refptlabel, scaled_length=scaled_length,
                       endlabel=endlabel) +
        facet_grid(replicate~group) +
        theme(strip.text.y=element_text(angle=0))
 
    hmap.width = max(12, (.0008*(upstream+scaled_length+dnstream)+3.4)*ngroups) 
    ggsave(samples_out, plot=hmap_sample, height=(.0005*nindices+7.5)*max(repl_df$replicate),
           width=hmap.width, units="cm", limitsize=FALSE)
    rm(hmap_sample, df_sample, repl_df)
    
    df_group = raw %>% group_by(group, index, position) %>% summarise(cpm = mean(cpm))
    rm(raw)
    group_cutoff = quantile(df_group$cpm, probs=pct_cutoff, na.rm=TRUE)
    df_group = df_group %>% 
        mutate(cpm = pmin(group_cutoff, cpm))

    hmap_group = hmap(df=df_group, refptlabel=refptlabel, scaled_length=scaled_length,
                      endlabel=endlabel) +
        facet_wrap(~group, ncol=ngroups)
    ggsave(group_out, plot = hmap_group, height= .0009*nindices+11.5,
           width=hmap.width, units="cm")
}

main(intable= args$input,
     samplelist = args$samplelist,
     type= args$type,
     upstream = args$upstream,
     dnstream= args$downstream,
     pct_cutoff = args$pct_cutoff,
     cluster = args$cluster,
     k = args$k,
     refptlabel = paste(args$refptlabel, collapse=" "),
     strand = tolower(args$strand),
     scaled_length = args$scaled_length,
     endlabel = paste(args$endlabel, collapse=" "),
     ylabel = paste(args$ylabel, collapse=" "),
     cmap = args$cmap,
     samples_out = args$samples_out,
     group_out = args$group_out)
