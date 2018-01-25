#!/usr/bin/env Rscript
library(argparse)
library(psych)
library(tidyverse)
library(forcats)
library(viridis)
library(ggthemes)

parser = ArgumentParser()
parser$add_argument('-i', dest='input', type='character')
parser$add_argument('-s', dest='samplelist', type='character', nargs='+')
parser$add_argument('-t', dest='type', type='character')
parser$add_argument('-f', dest='strand', type='character')
parser$add_argument('-u', dest='upstream', type='integer')
parser$add_argument('-d', dest='downstream', type='integer')
parser$add_argument('-c', dest='trim_pct', type='double')
parser$add_argument('-r', dest='refptlabel', type='character', nargs='+')
parser$add_argument('-l', dest='scaled_length', type='integer')
parser$add_argument('-e', dest='endlabel', type='character', nargs='+')
parser$add_argument('-y', dest='ylabel', type='character', nargs='+')
parser$add_argument('--out1', dest='meta_sample', type='character')
parser$add_argument('--out2', dest='meta_sample_overlay', type='character')
parser$add_argument('--out3', dest='meta_heatmap_sample', type='character')
parser$add_argument('--out4', dest='meta_group', type='character')
parser$add_argument('--out5', dest='meta_group_overlay', type='character')
parser$add_argument('--out6', dest='meta_heatmap_group', type='character')

args = parser$parse_args()

##stat_stepribbon by user 'felasa':
##https://groups.google.com/forum/?fromgroups=#!topic/ggplot2/9cFWHaH1CPs
#stairstepn <- function( data, direction="hv", yvars="y" ) {
#    direction <- match.arg( direction, c( "hv", "vh" ) )
#    data <- as.data.frame( data )[ order( data$x ), ]
#    n <- nrow( data )
#    
#    if ( direction == "vh" ) {
#        xs <- rep( 1:n, each = 2 )[ -2 * n ]
#        ys <- c( 1, rep( 2:n, each = 2 ) )
#    } else {
#        ys <- rep( 1:n, each = 2 )[ -2 * n ]
#        xs <- c( 1, rep( 2:n, each = 2))
#    }
#    
#    data.frame(
#        x = data$x[ xs ]
#        , data[ ys, yvars, drop=FALSE ]
#        , data[ xs, setdiff( names( data ), c( "x", yvars ) ), drop=FALSE ]
#    ) 
#}
#
#stat_stepribbon <- 
#    function(mapping = NULL, data = NULL, geom = "ribbon", position = "identity", inherit.aes = TRUE) {
#        ggplot2::layer(
#            stat = Stepribbon, mapping = mapping, data = data, geom = geom, 
#            position = position, inherit.aes = inherit.aes
#        )
#    }
#
#Stepribbon <- 
#    ggproto("stepribbon", Stat,
#            compute_group = function(., data, scales, direction = "hv", yvars = c( "ymin", "ymax" ), ...) {
#                stairstepn( data = data, direction = direction, yvars = yvars )
#            },                        
#            required_aes = c( "x", "ymin", "ymax" )
#    )
#
#format_xaxis = function(refptlabel, upstream, dnstream){
#    function(x){
#        if (first(upstream)>500 | first(dnstream)>500){
#            return(if_else(x==0, refptlabel, as.character(x)))    
#        }    
#        else {
#            return(if_else(x==0, refptlabel, as.character(x*1000)))
#        }
#    }
#}

theme_default = theme_light() +
    theme(text = element_text(size=12, color="black", face="bold"),
          strip.background = element_blank(),
          strip.text = element_text(size=12, color="black", face="bold"),
          strip.text.y = element_text(angle=-180),
          axis.text.x = element_text(size=12, color="black", face="bold"),
          axis.text.y = element_text(size=10, color="black", face="plain"),
          axis.title = element_text(size=10, face="plain"),
          plot.title = element_text(size=12),
          plot.subtitle = element_text(size=12, face="plain"))
    
main = function(intable, samplelist, type, strand, upstream, dnstream, trim_pct,
                refptlabel, scaled_length, endlabel, ylabel, meta_sample_out, 
                meta_sample_overlay_out, meta_heatmap_sample_out, meta_group_out,
                meta_group_overlay_out, meta_heatmap_group_out){
    x_label = function(ggp){
        if(type=="absolute"){
            ggp = ggp +
                scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                                   labels=format_xaxis(refptlabel=refptlabel,
                                                       upstream=upstream,
                                                       dnstream=dnstream),
                                   name=paste("distance from", refptlabel,
                                              if_else(upstream>500 | dnstream>500, "(kb)", "(nt)")),
                                   limits = c(-upstream/1000, (scaled_length+dnstream)/1000),
                                   expand=c(0,0))
        }
        else {
            ggp = ggp +
                scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                                   labels=c(refptlabel, "", endlabel),
                                   name="scaled distance",
                                   limits = c(-upstream/1000, (scaled_length+dnstream)/1000),
                                   expand=c(0,0))
            
        }
        return(ggp)
    }

    meta = function(df){
        metagene_base = ggplot(data = df, aes(x=position, y=mean, ymin = mean-1.96*sem,
                                    ymax=mean+1.96*sem)) +
        geom_vline(xintercept = c(0, scaled_length/1000), size=1, color="grey65")
        if(type=="scaled"){
            metagene_base = metagene_base +
                geom_vline(xintercept=scaled_length/1000, size=1, color="grey65")
        }
        metagene_base = metagene_base +
            geom_ribbon(fill="#114477", alpha=0.4, size=0) +
            geom_line(color="#114477", alpha=0.9) +
            scale_y_continuous(limits=c(0, NA), "normalized counts") +
            ggtitle(paste("mean", strand, "NET-seq signal"),
                    subtitle = paste(nindices, ylabel)) +
            theme_default
        return(x_label(metagene_base))
    }
    meta_overlay = function(df){
        metagene_base = ggplot(data = df,
                      aes(x=position, y=mean, ymin = mean-1.96*sem,
                          ymax=mean+1.96*sem, fill=group, color=group)) +
        geom_vline(xintercept = c(0, scaled_length/1000), size=1, color="grey65")
        if(type=="scaled"){
            metagene_base = metagene_base +
                geom_vline(xintercept=scaled_length/1000, size=1, color="grey65")
        }
        metagene_base = metagene_base +
            scale_fill_ptol(name=NULL) +
            scale_color_ptol(name=NULL) +
            scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                               labels=c(refptlabel, "", endlabel),
                               name="scaled distance",
                               limits = c(-upstream/1000, (dnstream+scaled_length)/1000),
                               expand=c(0,0)) +
            scale_y_continuous(limits=c(0, NA), name="normalized counts") +
            ggtitle(paste("mean", strand, "NET-seq signal"),
                    subtitle = paste(nindices, ylabel)) +
            theme_default +
            theme(legend.text = element_text(size=12))
        return(x_label(metagene_base))
    }
    
    meta_heatmap = function(df){
        metagene_base = ggplot(data = df, aes(x=position, y=0, fill=log2(mean+pcount))) +
            geom_raster() +
            scale_fill_viridis(option="inferno", name=expression(bold(log[2] ~ signal)),
                               guide=guide_colorbar(barheight=8, barwidth=1)) +
            scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                               labels=c(refptlabel, "", endlabel),
                               name="scaled distance",
                               limits = c(-upstream/1000, (dnstream+scaled_length)/1000),
                               expand=c(0,0)) +
            scale_y_continuous(breaks=0, expand=c(0,0), name=NULL) +
            ggtitle(paste("mean", strand, "NET-seq signal"),
                    subtitle = paste(nindices, ylabel)) +
            theme_default +
            theme(axis.text.y = element_blank(),
                  strip.text.y = element_text(hjust=1),
                  legend.title = element_text(size=10, face="plain"),
                  axis.ticks.x = element_line(color="black", size=1),
                  axis.ticks.y = element_blank(),
                  panel.border = element_blank())
        return(x_label(metagene_base))
    }
    
    raw = read_tsv(intable, col_names=c("group", "sample", "index", "position","cpm")) %>%
        filter(sample %in% samplelist & !is.na(cpm)) %>% 
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE)))
        
    nindices = max(raw$index, na.rm=TRUE)
    nsamples = length(samplelist)
    ngroups = length(fct_unique(raw$group))

    pcount=0.1
    #get replicate info for sample facetting
    repl_df = raw %>% select(group, sample) %>% distinct() %>% group_by(group) %>%
        mutate(replicate=row_number()) %>% ungroup() %>% select(-group)
    
    #plot heatmap facetted by sample and group
    df_sample = raw %>% group_by(group, sample, position) %>% 
        summarize(mean = winsor.mean(cpm, trim = trim_pct),
                  sem = winsor.sd(cpm, trim = trim_pct)/sqrt(n())) %>% 
        left_join(repl_df, by="sample")
    
    meta_sample = meta(df_sample) +
        facet_grid(replicate~group) +
        theme(strip.text.y = element_text(angle=0))
    ggsave(meta_sample_out, plot=meta_sample, height=2+4.5*max(repl_df$replicate),
           width=7*ngroups, units="cm", limitsize=FALSE)
    rm(meta_sample)
    
    meta_sample_overlay = meta_overlay(df_sample) +
        geom_ribbon(aes(group=sample), alpha=0.3, size=0) +
        geom_line(aes(group=sample), alpha=0.9)
    ggsave(meta_sample_overlay_out, plot=meta_sample_overlay,
           height=8, width=14, units="cm")
    rm(meta_sample_overlay)
    
    meta_heatmap_sample = meta_heatmap(df_sample) +
        facet_grid(sample~., switch="y")
    ggsave(meta_heatmap_sample_out, plot=meta_heatmap_sample,
           height=2.5+1.25*nsamples, width=14, units="cm")
    rm(meta_heatmap_sample, df_sample, repl_df)
    
    df_group = raw %>% group_by(group, position) %>%
        summarise(mean = winsor.mean(cpm, trim=trim_pct),
                  sem = winsor.sd(cpm, trim=trim_pct)/sqrt(n()))
    rm(raw)
    
    meta_group = meta(df_group) +
        facet_grid(.~group)
    ggsave(meta_group_out, plot=meta_group, height=8,
           width=7*ngroups, units="cm", limitsize=FALSE)
    rm(meta_group_out)
    
    meta_group_overlay = meta_overlay(df_group) +
        geom_ribbon(aes(group=group), alpha=0.3, size=0) +
        geom_line(aes(group=group), alpha=0.9)
    ggsave(meta_group_overlay_out, plot=meta_group_overlay,
           height=8, width=14, units="cm")
    rm(meta_group_overlay)
    
    meta_heatmap_group = meta_heatmap(df_group) +
        facet_grid(group~., switch="y")
    ggsave(meta_heatmap_group_out, meta_heatmap_group,
           height=2.5+1.5*ngroups, width=14, units="cm")
}

main(intable= args$input,
     samplelist = args$samplelist,
     type= args$type,
     strand = tolower(args$strand),
     upstream = args$upstream,
     dnstream= args$downstream,
     trim_pct = args$trim_pct,
     refptlabel = paste(args$refptlabel, collapse=" "),
     scaled_length = args$scaled_length,
     endlabel = paste(args$endlabel, collapse=" "),
     ylabel = paste(args$ylabel, collapse=" "),
     meta_sample_out = args$meta_sample,
     meta_sample_overlay_out = args$meta_sample_overlay,
     meta_heatmap_sample_out = args$meta_heatmap_sample,
     meta_group_out = args$meta_group,
     meta_group_overlay_out = args$meta_group_overlay,
     meta_heatmap_group_out = args$meta_heatmap_group)
