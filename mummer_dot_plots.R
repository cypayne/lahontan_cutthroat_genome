## MUMmer alignments visualized as dot plots
## Comparing Lahontan Cutthroat Trout genome assembly with other salmonid genome assemblies
## Code creates dotplot figures that appear in Payne & Escalona et al 2025. JHeredity.
## The following code is adapted from code written by @jmonlong 
## https://jmonlong.github.io/Hippocamplus/2017/09/19/mummerplots-with-ggplot2/

library(dplyr)
library(magrittr)
library(knitr)
library(ggplot2)
library(tidyr)

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



#### Lahontan vs Westslope Cutthroat ####

## Figure 6 C-i
## OclaPP08 (JBEFNN010000002.1) vs OclaCL08 + OclaCL28
## fusion
## OclaCL28 0-49.7 Mb = OclaPP08 0-52.4 Mb
## OclaCL08 0-42.0 Mb = OclaPP08 52.9-95.8 Mb
deltas_path <- "LCT_v_WCT_mummer/fission_fusion_alignments/"
mumgp1 = readDelta(paste0(deltas_path,"JBEFNN010000002.1_v_CM099244.1_nucmer.delta.m"))
mumgp.filt1 = filterMum(mumgp1, minl=18000)
mumgp2 = readDelta(paste0(deltas_path,"JBEFNN010000002.1_v_CM099264.1_nucmer.delta.m"))
mumgp.filt2 = filterMum(mumgp2, minl=18000)

mumgp.filt <- mumgp.filt1 %>% bind_rows(mumgp.filt2) %>% mutate(qid = recode(qid, CM099244.1 = "chr08",CM099264.1 = "chr28"),rid = recode(rid, JBEFNN010000002.1 = "chr08"))
mumgp.filt$qid <- relevel(mumgp.filt$qid,ref="chr08")

dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(qid~., scales='free', space='free', switch='y') + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_blank(),
                     strip.text.y=element_text(size=14),
                     strip.placement="outside",
                     legend.position=c(.1,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16),
                     plot.margin = margin(t=0.2,r=0.4,b=0.2,l=0.2, "cm")) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10),limits=c(0,100),expand=c(0.02,0)) + 
  scale_colour_brewer(palette='Set1') +
  ylab('OclaCL (Mb)') + xlab('OclaPP chr08 (Mb)') 
dotplot
ggsave(dotplot, filename='dotplot_nucmer_OclaPP08_OclaCL08_OclaCL28.pdf',  width=8.5, height=8, bg="transparent")


## Figure 6 C-ii
## OclaPP14 (JBEFNN010000006.1) vs OclaCL14 + OclaCL32
## fusion
# OclaCL14 0-42.7 Mb = OclaPP14 0-45.4 Mb
# OclaCL32 0-40.6 Mb = OclaPP14 45.9-88.9 Mb
deltas_path <- "LCT_v_WCT_mummer/fission_fusion_alignments/"
mumgp1 = readDelta(paste0(deltas_path,"JBEFNN010000006.1_v_CM099250.1_nucmer.delta.m"))
mumgp.filt1 = filterMum(mumgp1, minl=18000)
mumgp2 = readDelta(paste0(deltas_path,"JBEFNN010000006.1_v_CM099268.1_nucmer.delta.m"))
mumgp.filt2 = filterMum(mumgp2, minl=18000)

mumgp.filt <- mumgp.filt1 %>% bind_rows(mumgp.filt2) %>% mutate(qid = recode(qid, CM099250.1 = "chr14", CM099268.1 = "chr32"),rid = recode(rid, JBEFNN010000006.1 = "chr14"))

dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(qid~., scales='free', space='free', switch='y') + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_blank(),
                     strip.text.y=element_text(size=14),
                     strip.placement="outside",
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16),
                     plot.margin = margin(t=0.2,r=0.4,b=0.2,l=0.2, "cm")) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10),limits=c(0,90),expand=c(0.02,0)) + 
  scale_colour_brewer(palette='Set1') +
  ylab('OclaCL (Mb)') + xlab('OclaPP chr14 (Mb)') 
dotplot
ggsave(dotplot, filename='dotplot_nucmer_OclaPP14_OclaCL14_OclaCL32.pdf',  width=8.5, height=8, bg="transparent")


## Figure 6 C-iii
## OclaPP05 (JBEFNN010000021.1) + OclaPP28 (JBEFNN010000023.1) vs OclaCL05
## fission
## OcaPP05 0-51.9 Mb =	OclaCL05 0-49.2 Mb
## OcaPP28 0-50.4 Mb =	OclaCL05 49.6-92.6 Mb
deltas_path <- "LCT_v_WCT_mummer/fission_fusion_alignments/"
mumgp1 = readDelta(paste0(deltas_path,"JBEFNN010000021.1_v_CM099241.1_nucmer.delta.m"))
mumgp.filt1 = filterMum(mumgp1, minl=18000)
mumgp2 = readDelta(paste0(deltas_path,"JBEFNN010000023.1_v_CM099241.1_nucmer.delta.m"))
mumgp.filt2 = filterMum(mumgp2, minl=18000)

mumgp.filt <- mumgp.filt1 %>% bind_rows(mumgp.filt2) %>% mutate(qid = recode(qid, CM099241.1 = "chr05"),rid = recode(rid, JBEFNN010000021.1 = "chr05", JBEFNN010000023.1 = "chr28"))
#mumgp.filt$rid <- relevel(mumgp.filt$rid,ref="chr05")

dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(.~rid, scales='free', space='free', switch='x') + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(size=14),
                     strip.text.y=element_blank(),
                     strip.placement="outside",
                     legend.position=c(.99,.87), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_y_reverse(breaks=scales::pretty_breaks(n=10),limits=c(94,0),expand=c(0.02,0)) + 
  scale_colour_brewer(palette='Set1') +
  ylab('OclaCL chr05 (Mb)') + xlab('OclaPP (Mb)') 
dotplot
ggsave(dotplot, filename='dotplot_nucmer_OclaPP05_OclaCL28_OclaCL05.pdf',  width=8.5, height=8, bg="transparent")


## Figure 4E
## OclaPP22 (JBEFNN010000017.1) vs OclaCL22
## translocated inversion: OclaPP22 47.9-50.8 Mb = OclaCL22 14.5-17.3 Mb
## the segment it translocated around is OclaPP22 41-47.9Mb = OcaCL22 7.6-14.3 ~7mb
deltas_path <- "LCT_v_WCT_mummer/"
mumgp = readDelta(paste0(deltas_path,"JBEFNN010000017.1_v_CM099258.1_nucmer.delta.m"))
mumgp.filt = filterMum(mumgp, minl=18000)

mumgp.filt <- mumgp.filt %>% mutate(qid = recode(qid, CM099258.1 = "chr22"),rid = recode(rid, JBEFNN010000017.1 = "chr22"))

dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(size=14),
                     strip.text.y=element_blank(),
                     strip.placement="outside",
                     legend.position=c(.1,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=6)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=6)) + 
  scale_colour_brewer(palette='Set1') + 
  ylab('OclaCL chr22 (Mb)') + xlab('OclaPP chr22 (Mb)') 
dotplot
ggsave(dotplot, filename='dotplot_nucmer_OclaPP22_OclaCL22.pdf',  width=8.5, height=8, bg="transparent")

# zoom in
dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.09,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') + 
  #xlim(40,49)+ ylim(18,5) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10), limits=c(40,52)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=10), limits=c(18,6)) + 
  ylab('OclaCL chr22 (Mb)') + xlab('OclaPP chr22 (Mb)') 
dotplot
ggsave(dotplot, filename='dotplot_nucmer_OclaCL22_OclaPP22_zoomin.pdf', width=8, height=8, bg="transparent")


## OclaPP13 (JBEFNN010000027.1) vs OclaCL13
# inversion? may not be real: OclaPP13 41.8-48.0 Mb = OclaCL13 5.5-15.4 Mb
mumgp = readDelta("OclaCL13_OclaPP27_nucmer.delta.m") 
mumgp.filt = filterMum(mumgp, minl=1e4)

dotplot <-
  ggplot(mumgp.filt, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.13,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') +
  ylab('OclaCL chr13 (Mb)') + xlab('OclaPP chr13 (Mb)') 

dotplot
ggsave(dotplot, filename='dotplot_nucmer_OclaCL13_OclaPP13.pdf', width=8, height=8, bg="transparent")

# zoom in
dotplot <-
  ggplot(mumgp.filt, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.08,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') + 
  #xlim(40,49)+ ylim(18,5) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10), limits=c(40,49)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=10), limits=c(17,5)) + 
  ylab('OclaCL chr13 (Mb)') + xlab('OclaPP chr13 (Mb)') 
dotplot
ggsave(dotplot, filename='dotplot_nucmer_OclaCL13_OclaPP13_zoomin.pdf', width=8, height=8, bg="transparent")



#### Lahontan Cutthroat vs Rainbow Arlee ####

## Closer inspection of alignments of fusions and fissions 
## OclaPP04 0-49.7Mb = OmyA04 45.8-0Mb (-)
## OclaPP04 56.2-89.4Mb = OmyA30 43.3-11.8Mb (-)
mumgp = readDelta("LCT_v_OmyA_mummer/JBEFNN010000001.1_v_NC_048568.1_nucmer.delta.m") # OmyA04
mumgp = readDelta("LCT_v_OmyA_mummer/JBEFNN010000001.1_v_NC_050570.1_nucmer.delta.m") # OmyA30
## OclaPP14 0-45.0Mb = OmyA14 42.8-0Mb (-)
## OclaPP14 46.3-89.0Mb = OmyA32 41.5-0Mb (-)
mumgp = readDelta("LCT_v_OmyA_mummer/JBEFNN010000006.1_v_NC_048578.1_nucmer.delta.m") # OmyA14
mumgp = readDelta("LCT_v_OmyA_mummer/JBEFNN010000006.1_v_NC_050572.1_nucmer.delta.m") # OmyA32
## OclaPP11 0-39.8Mb = OmyA11 38.6-0Mb (-)
## OclaPP30 0-38.0Mb = OmyA11 39.8-74.7Mb (+)
mumgp = readDelta("LCT_v_OmyA_mummer/JBEFNN010000030.1_v_NC_048575.1_nucmer.delta.m") # OclaPP11
mumgp = readDelta("LCT_v_OmyA_mummer/JBEFNN010000022.1_v_NC_048575.1_nucmer.delta.m") # OclaPP30
## OclaPP15 0-28.1Mb = OmyA15 38.2-14.4Mb (-)
## OclaPP32 0-25.0 = OmyA15 39.4-62.8Mb (+)
mumgp = readDelta("LCT_v_OmyA_mummer/JBEFNN010000029.1_v_NC_048579.1_nucmer.delta.m") # OclaPP15
mumgp = readDelta("LCT_v_OmyA_mummer/JBEFNN010000032.1_v_NC_048579.1_nucmer.delta.m") # OclaPP32

## OclaPP05 (JBEFNN010000021.1) + OclaPP28 (JBEFNN010000023.1) vs OmyA05
## fission
## OclaPP05 0-51.7Mb = OmyA05 49.6-0Mb (-)
## OclaPP28 0-50.4Mb = OmyA05 100.8-50.1Mb (-)
mumgp = readDelta("Omy05_Ocla21_Ocla23_nucmer.delta.m")
mumgp.filt = filterMum(mumgp, minl=1e4)

mumgp.filt <- mumgp.filt %>% mutate(qid = recode(qid, JBEFNN010000023.1 = "chr28", JBEFNN010000021.1 = "chr05" ))
mumgp.filt$qid <- relevel(mumgp.filt$qid,ref="chr05")

dotplot <-
  ggplot(mumgp.filt, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(.~qid, scales='free', space='free', switch='x') + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(size=14),
                     strip.text.y=element_blank(),
                     strip.placement="outside",
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') +
  ylab('OmyA chr05 (Mb)') + xlab('OclaPP (Mb)') 
dotplot
ggsave(dotplot, filename='dotplot_nucmer_OmyA05_OclaPP05_Ocla28.pdf',  width=8.5, height=8, bg="transparent")


## Figure 4B
## OclaPP20 (JBEFNN010000004.1) vs OmyA20 + OmyA28 
## fusion
## OclaPP20 0-47.5 Mb = OmyA20 46.6-0 Mb (-)
## OclaPP20 48.3-93.4 Mb = OmyA28 0-43.7 Mb (+)
## inversion: OclaPP20 26.5-37.4Mb = OmyA20 6.5-17.3Mb
## translocated inversion: OclaPP20 42.7-44.9Mb = OmyA20 4.7-6.5Mb
## translocation: OclaPP20 44.9-47.1Mb = OmyA20 19.0-21.2Mb
mumgp = readDelta("Omy20_Omy28_Ocla4_nucmer.delta.m")
mumgp.filt = filterMum(mumgp, minl=18000)

mumgp.filt <- mumgp.filt %>% mutate(rid = recode(rid, NC_048584.1 = "chr20", NC_048592.1 = "chr28" ))

# pull out structural rearrangement coordinates
# doesn't help haha
subset(mumgp.filt,strand=="+" & qs>25000000 & rs>6000000 & rs<40000000 & strand=="+")
max(subset(mumgp.filt,rid=="chr20")$re)

dotplot <-
  ggplot(mumgp.filt, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(rid~., scales='free', space='free', switch='y') + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_blank(),
                     strip.text.y=element_text(size=14),
                     strip.placement="outside",
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16),
                     panel.grid.minor=element_blank()) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=8)) + 
  scale_colour_brewer(palette='Set1') +
  ylab('OmyA (Mb)') + xlab('OclaPP chr20 (Mb)')  
dotplot
ggsave(dotplot, filename='dotplot_nucmer_OmyA20_OmyA28_OclaPP20.pdf',  width=8.5, height=8, bg="transparent")

# zoom in
ocla20_v_omy20 <- subset(mumgp.filt,rid=="chr20")
dotplot <-
  ggplot(ocla20_v_omy20, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_blank(),
                     strip.text.y=element_text(size=14),
                     strip.placement="outside",
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 13),
                     axis.title = element_text(size = 16)) +
#  xlim(25,48) + ylim(20,0) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=20), limits=c(25.5,47.0)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=20), limits=c(22,0.5)) + 
  scale_colour_brewer(palette='Set1') +
  ylab('OmyA chr20 (Mb)') + xlab('OclaPP chr20 (Mb)')  
dotplot
ggsave(dotplot, filename='dotplot_nucmer_OmyA20_OclaPP20_zoomin.pdf',  width=8, height=8, bg="transparent")


## Figure 4D
## OclaPP29 (JBEFNN010000020.1) vs OmyA29
## inversion: OmyA29 6.6-14.2Mb = OclaPP29 6.3-14.0Mb
## translocated inversion: OmyA29 3.0-5.1Mb = OclaPP29 0.8-2.9Mb
## translocated around a ~3Mb segment
mumgp = readDelta("Omy29_Ocla20_nucmer.delta.m") 
mumgp.filt = filterMum(mumgp, minl=18000)

dotplot <-
  ggplot(mumgp.filt, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') +
  ylab('OmyA chr29 (Mb)') + xlab('OclaPP chr29 (Mb)') 

dotplot
ggsave(dotplot, filename='dotplot_nucmer_Omy29_Ocla29.pdf', width=8, height=8, bg="transparent")

# zoom in
dotplot <-
  ggplot(mumgp.filt, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=20), limits=c(0,16)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=20), limits=c(16,0)) + 
  ylab('OmyA chr29 (Mb)') + xlab('OclaPP chr29 (Mb)') 

dotplot
ggsave(dotplot, filename='dotplot_nucmer_Omy29_Ocla29_zoomin.pdf', width=8, height=8, bg="transparent")


## Figure 4C
## OclaPP22 (JBEFNN010000017.1) vs OmyA22 (NC_048586.1)
## translocated inversion: OmyA22 7.3-12.6Mb = OclaPP22 34.8-40.5Mb  
## the segment it translocated around is OclaPP22 41-50.8 = OmyA 12.6-22.1 ~10mb
mumgp = readDelta("LCT_v_OmyA_mummer/JBEFNN010000017.1_v_NC_048586.1_nucmer.delta.m") 
mumgp.filt = filterMum(mumgp, minl=18000)

mumgp.filt <- mumgp.filt %>% mutate(qid = recode(qid, NC_048586.1 = "chr22"),rid = recode(rid, JBEFNN010000017.1 = "chr22"))

dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.1,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=6)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=6),limits=c(52,0)) + 
  scale_colour_brewer(palette='Set1') +
  ylab('OmyA chr22 (Mb)') + xlab('OclaPP chr22 (Mb)') 

dotplot
ggsave(dotplot, filename='dotplot_nucmer_OmyA22_OclaPP22.pdf', width=8, height=8, bg="transparent")

# zoom in
dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.09,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=20), limits=c(33.5,51.5)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=20), limits=c(24.5,6.5)) + 
  ylab('OmyA chr22 (Mb)') + xlab('OclaPP chr22 (Mb)') 

dotplot
ggsave(dotplot, filename='dotplot_nucmer_Omy22_Ocla22_zoomin.pdf', width=8, height=8, bg="transparent")


## Figure 4A
## OclaPP17 (JBEFNN010000013.1) vs OmyA17
## inversion: OclaPP17 19.0-26.0Mb = OmyA17 36.3-43.1Mb
mumgp = readDelta("LCT_v_OmyA_mummer/JBEFNN010000013.1_v_NC_048581.1_nucmer.delta.m") 
mumgp.filt = filterMum(mumgp, minl=18000)

dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16),
                     panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=8)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=10), limits=c(95,0)) + 
  scale_colour_brewer(palette='Set1') +
  ylab('OmyA chr17 (Mb)') + xlab('OclaPP chr17 (Mb)') 

dotplot
ggsave(dotplot, filename='dotplot_nucmer_OmyA17_OclaPP17.pdf', width=8, height=8, bg="transparent")



#### Lahontan Cutthroat vs Atlantic Salmon ####

## Figure S7A
## OclaPP20 (JBEFNN010000004.1) vs Ssa03 + Ssa28
mumgp = readDelta("Ssa03_Ssa28_Ocla4_nucmer.delta.m")
mumgp.filt = filterMum(mumgp, minl=1e4)

mumgp.filt <- mumgp.filt %>% mutate(rid = recode(rid, NC_027302.1 = "chr03", NC_027327.1 = "chr28" ))
mumgp.filt <- subset(mumgp.filt,rs<45000000)
dotplot <-
  ggplot(mumgp.filt, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(rid~., scales='free', space='free', switch='y') + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_blank(),
                     strip.text.y=element_text(size=14),
                     strip.placement="outside",
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') +
  ylab('SsaSally (Mb)') + xlab('OclaPP chr20 (Mb)') +
  expand_limits(y = 0)
dotplot
ggsave(dotplot, filename='dotplot_nucmer_Ssa03_Ssa28_OclaPP20.pdf',  width=8.5, height=8, bg="transparent")


## OclaPP22 (JBEFNN010000017.1) vs Ssa21 (NC_027320.1)
## inversion: Ssa21 9.1-19Mb = OclaPP22 41-50.9Mb  
mumgp = readDelta("LCT_v_SsaSally_mummer/JBEFNN010000017.1_v_NC_027320.1_nucmer.delta.m") 
mumgp.filt = filterMum(mumgp, minl=1000)

mumgp.filt <- mumgp.filt %>% mutate(qid = recode(qid, NC_027320.1 = "chr21"),rid = recode(rid, JBEFNN010000017.1 = "chr22"))

dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.1,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=6)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=6)) + 
  scale_colour_brewer(palette='Set1') +
  ylab('SsaSally chr21 (Mb)') + xlab('OclaPP chr22 (Mb)') 

dotplot
ggsave(dotplot, filename='dotplot_nucmer_SsaS21_OclaPP22.pdf', width=8, height=8, bg="transparent")

# zoom in
dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.09,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=20), limits=c(40,55)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=20), limits=c(20,5)) + 
  ylab('SsaSally chr21 (Mb)') + xlab('OclaPP chr22 (Mb)') 

dotplot
ggsave(dotplot, filename='dotplot_nucmer_Ssa21_OclaPP22_zoomin.pdf', width=8, height=8, bg="transparent")


## OclaPP29 (JBEFNN010000020.1) vs Ssa11 (NC_027310.1)
## inversion: OmyA29 6.6-14.2Mb = OclaPP29 6.3-14.0Mb
## translocated inversion: OmyA29 3.0-5.1Mb = OclaPP29 0.8-2.9Mb
## translocated around a ~3Mb segment
mumgp = readDelta("LCT_v_SsaSally_mummer/JBEFNN010000020.1_v_NC_027310.1_nucmer.delta.m") 
mumgp.filt = filterMum(mumgp, minl=1000)

dotplot <-
  ggplot(mumgp.filt, aes(x=rs/scalar, xend=re/scalar, y=qs/scalar, yend=qe/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') +
  ylab('OmyA chr29 (Mb)') + xlab('OclaPP chr29 (Mb)') 

dotplot
ggsave(dotplot, filename='dotplot_nucmer_Ssa11_Ocla29.pdf', width=8, height=8, bg="transparent")

# zoom in
dotplot <-
  ggplot(mumgp.filt, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_text(angle=180, size=5),
                     strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=20), limits=c(0,16)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=20), limits=c(16,0)) + 
  ylab('OmyA chr29 (Mb)') + xlab('OclaPP chr29 (Mb)') 

dotplot
ggsave(dotplot, filename='dotplot_nucmer_Omy29_Ocla29_zoomin.pdf', width=8, height=8, bg="transparent")


#### Rainbow Trout Arlee vs Atlantic Salmon ####

## Figure S7B
## Ssa03+Ssa28 vs Omy20+Omy28
mumgp = readDelta("Ssa03_Ssa28_Omy20_Omy28_nucmer.delta.m")
mumgp.filt = filterMum(mumgp, minl=1e4)

mumgp.filt <- mumgp.filt %>% mutate(rid = recode(rid, NC_027302.1 = "chr03", NC_027327.1 = "chr28" ))
mumgp.filt <- mumgp.filt %>% mutate(qid = recode(qid, NC_048584.1 = "chr20", NC_048592.1 = "chr28" ))
mumgp.filt$qid <- relevel(mumgp.filt$qid,ref="chr20")
mumgp.filt <- subset(mumgp.filt,rs<45000000)

dotplot <-
  ggplot(mumgp.filt, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(rid~qid, scales='free', space='free', switch='both') + scale_y_reverse() +
  theme_bw() + theme(strip.text=element_text(size=14),
                     strip.placement="outside",
                     legend.position=c(.1,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 16)) +
  scale_colour_brewer(palette='Set1') +
  ylab('SsaSally (Mb)') + xlab('OmyA (Mb)') +
  expand_limits(y = 0)
dotplot 
ggsave(dotplot, filename='~/Desktop/trout_projects/cutthroat_assembly/dotplot_nucmer_Ssa03+Ssa28_Omy20+Omy28.pdf',  width=8.5, height=8.5, bg="transparent")

submum <- subset(mumgp.filt,rid=="chr28")
dotplot <-
  ggplot(submum, aes(x=qs/scalar, xend=qe/scalar, y=rs/scalar, yend=re/scalar, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + scale_y_reverse() +
  theme_bw() + theme(strip.text.x=element_blank(),
                     strip.text.y=element_text(size=14),
                     strip.placement="outside",
                     legend.position=c(.99,.88), legend.justification=c(1,0), legend.text=element_text(size = 14),
                     strip.background=element_blank(),
                     axis.text = element_text(size = 13),
                     axis.title = element_text(size = 16)) +
  #  xlim(25,48) + ylim(20,0) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=20), limits=c(0,30)) + 
  scale_y_reverse(breaks=scales::pretty_breaks(n=20), limits=c(25,0)) + 
  scale_colour_brewer(palette='Set1') +
  ylab('OmyA chr20 (Mb)') + xlab('OclaPP chr20 (Mb)')  
dotplot



