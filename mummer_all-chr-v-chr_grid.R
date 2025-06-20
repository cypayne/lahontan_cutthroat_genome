## Matrix of pairwise MUMmer alignments visualized as dot plots between assemblies
## Comparing Lahontan Cutthroat Trout genome assembly with other salmonid genome assemblies
## Code creates an all chromosomes versus all chromosomes matrix of dotplots between 
## genome assemblies. These figures appear in Payne & Escalona et al 2025. JHeredity.
## The following code is adapted from code written by @jmonlong 
## https://jmonlong.github.io/Hippocamplus/2017/09/19/mummerplots-with-ggplot2/

library(patchwork)
library(dplyr)
library(magrittr)
library(knitr)
library(ggplot2)
library(tidyr)

###  functions and global variables
# read in filtered delta file
readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1] 
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

# filter out poor alignments
# i.e. keep contigs with at least one aligned segment larger than a threshold size
# smaller alignments in those contigs are kept if they're in the same range as the 
# large aligned segments
filterMum <- function(df, minl=1000, flanks=1e4){
  coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
    summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
    ungroup %>% arrange(desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
  merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
    mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

# set number to divide the position by. 
scalar <- 1000000 # i.e. units will be in Mb


### plots

## Choose the comparison to make by setting the path to the alignment delta files
deltas_path <- "LCT_v_WCT_mummer/"
#deltas_path <- "LCT_v_OmyA_mummer/"
#deltas_path <- "LCT_v_OmyS_mummer/"
#deltas_path <- "LCT_v_SsaSally_mummer/"
#deltas_path <- "OmyA_v_SsaSally_mummer/"

# read in the order you want the chromosomes on each axis
x_chr_order <- read.table(paste0(deltas_path,"fOncCla1.0.p_chr-name2num.txt"),header=T)
#x_chr_order <- read.table(paste0(deltas_path,"OmyArlee_chr-name2num.txt"),header=T)[1:32,]
#y_chr_order <- read.table("LCT_v_WCT_mummer/WCT-UVic2024_chr-name2num.txt",header=T)
#y_chr_order <- read.table(paste0(deltas_path,"OmyArlee_chr-name2num.txt"),header=T)[1:32,]
y_chr_order <- read.table(paste0(deltas_path,"OmySwanson_chr-name2num.txt"),header=T)[1:29,]
#y_chr_order <- read.table(paste0(deltas_path,"Ssalar.Sally_chr-name2num.txt"),header=T)[1:29,]

# for each chrom combo, grab the delta.m file
plot_list <- list()
chr_x_length_vector <- list()
chr_y_length_vector <- list()
counter <- 1
chr_x_list <- sort(x_chr_order$chr_num)
chr_y_list <- sort(y_chr_order$chr_num)
for(chr_y in chr_y_list) {
  for(chr_x in chr_x_list) {
    # get both chromosome names and lengths
    chr_x_name <- x_chr_order[x_chr_order$chr_num==chr_x,]$chr_name 
    chr_y_name <- y_chr_order[y_chr_order$chr_num==chr_y,]$chr_name
    chr_x_length <- x_chr_order[x_chr_order$chr_num==chr_x,]$chr_length / scalar
    chr_y_length <- y_chr_order[y_chr_order$chr_num==chr_y,]$chr_length / scalar
    
    delta_str <- paste0(deltas_path,chr_x_name,"_v_",chr_y_name,"_nucmer.delta.m")
    mumgp = readDelta(delta_str) # read the delta.m output
    #mumgp.filt = filterMum(mumgp, minl=18000) # filter the delta.m output
    mumgp.filt = filterMum(mumgp, minl=10000) # filter the delta.m output
    
    
    if(chr_x == chr_x_list[1] & chr_y == chr_y_list[length(chr_y_list)]) {
      # if the first chromosome in assembly X AND the last chromosome in 
      # assembly Y, add labels to the X AND Y axes
      dotplot <-
        ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
        geom_point(alpha=0.1) + scale_y_reverse() +
        theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                           strip.text.y=element_text(angle=180, size=5),
                           legend.position="none",
                           #legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                           panel.grid.minor = element_blank(),
                           strip.background=element_blank(),
                           axis.text = element_blank(),
                           axis.title = element_text(size = 60),
                           #axis.text.x = element_text(size = 14),
                           axis.title.x = element_text(size = 60)) +
        xlim(0,chr_x_length) + ylim(chr_y_length,0) +
        scale_colour_brewer(palette='Set1') +
        xlab(chr_x) + ylab(chr_y)
    } else if(chr_x == chr_x_list[1]) {
      # if the first chromosome in assembly X, add labels to the Y axis
      dotplot <-
        ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
        geom_point(alpha=0.1) + scale_y_reverse() +
        theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                         strip.text.y=element_text(angle=180, size=5),
                         legend.position="none",
                         #legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                         panel.grid.minor = element_blank(),
                         strip.background=element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.text = element_blank(),
                         axis.title = element_blank(),
                         #axis.text.y = element_text(size = 14),
                         axis.title.y = element_text(size = 60)) +
        xlim(0,chr_x_length) + ylim(chr_y_length,0) +
        scale_colour_brewer(palette='Set1') +
        ylab(chr_y)
    } else if(chr_y == chr_y_list[length(chr_y_list)]) {
        # if the last chromosome in assembly Y, add labels to the X axis
        dotplot <-
          ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
          geom_point(alpha=0.1) + scale_y_reverse() +
          theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                           strip.text.y=element_text(angle=180, size=5),
                           legend.position="none",
                           #legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                           panel.grid.minor = element_blank(),
                           strip.background=element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.text = element_blank(),
                           axis.title = element_blank(),
                           #axis.text.x = element_text(size = 14),
                           axis.title.x = element_text(size = 60)) +
          xlim(0,chr_x_length) + ylim(chr_y_length,0) +
          scale_colour_brewer(palette='Set1') +
          xlab(chr_x) 
    } else {
        # otherwise, add no labels or ticks to either axis
        dotplot <-
          ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
          geom_point(alpha=0.1) + scale_y_reverse() +
          theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                         strip.text.y=element_text(angle=180, size=5),
                         legend.position="none",
                         #legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                         panel.grid.minor = element_blank(),
                         strip.background=element_blank(),
                         axis.ticks = element_blank(),
                         axis.text = element_blank(),
                         axis.title = element_blank()) +
          xlim(0,chr_x_length) + ylim(chr_y_length,0) +
          scale_colour_brewer(palette='Set1') 
          #xlab(chr_x) + ylab(chr_y)
    }
    dotplot
    plot_list[[counter]] <- dotplot
    counter <- counter + 1
  }
}

# create vectors with the sizes of each chromosome
x_chr_order <- x_chr_order %>% arrange(chr_num) 
chr_x_length_vector <- x_chr_order$chr_length
y_chr_order <- y_chr_order %>% arrange(chr_num) 
chr_y_length_vector <- y_chr_order$chr_length

# arrange the plots so they're stacked on top of each other
#big_grid_plot <- wrap_plots(plot_list, nrow = 32, ncol = 32, widths=chr_x_length_vector,heights=chr_y_length_vector)
big_grid_plot <- wrap_plots(plot_list, nrow = 29, ncol = 32, widths=chr_x_length_vector,heights=chr_y_length_vector)
ggsave(paste0(deltas_path,"LCT_v_OmyS-10000_nucmer-grid-plot.png"), big_grid_plot, width = 49, height = 49, units = "in", dpi = 150)

# test
#big_grid_plot <- wrap_plots(plot_list, nrow = 10, ncol = 10, widths=chr_x_length_vector,heights=chr_y_length_vector)
#big_grid_plot
#ggsave(paste0(deltas_path,"LCT_v_WCT_nucmer-grid-plot_10x10.png"), big_grid_plot, width = 48, height = 48, units = "in", dpi = 150)
