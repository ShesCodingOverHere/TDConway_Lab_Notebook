library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)
library(ggplot2)
library(tidyr)


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

mumgp = readDelta("C:/Users/tdcon/Desktop/Pseudo_Male.delta")

mumgp %>% head %>% kable



filterMum <- function(df, minl=1000, flanks=1e4){
  coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
    summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
    ungroup %>% arrange(desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
  merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
    mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

mumgp$rid[which(PAF$RefID=="NC_046603.1")]="mtDNA"
mumgp$rid[which(PAF$RefID=="NC_046679.1")]="Chr2_MullerE"
mumgp$rid[which(PAF$RefID=="NC_046680.1")]="Chr3_MullerC"
mumgp$rid[which(PAF$RefID=="NC_046681.1")]="Chr4_MullerB"
mumgp$rid[which(PAF$RefID=="NC_046682.1")]="Chr5_MullerF"
mumgp$rid[which(PAF$RefID=="NC_046683.1")]="ChrX_MullerAD"

mumgp.filt = filterMum(mumgp, minl=1e4)
mumgp.filt %>% head %>% kable

table(mumgp$rid)


## Whole mummerplot
ggplot(mumgp.filt, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(qid~., scales='free', space='free', switch='y') +
  theme_bw() + theme(strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.01), legend.justification=c(1,0),
                     strip.background=element_blank(),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1') +
  ggtitle("Pseudoobscura-Male Affinis Main Matches")


## Diagonalize
diagMum <- function(df){
  ## Find best qid order
  rid.o = df %>% group_by(qid, rid) %>% summarize(base=sum(abs(qe-qs)),
                                                  rs=weighted.mean(rs, abs(qe-qs))) %>%
    ungroup %>% arrange(desc(base)) %>% group_by(qid) %>% do(head(., 1)) %>%
    ungroup %>% arrange(desc(rid), desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid)))
  ## Find best qid strand
  major.strand = df %>% group_by(qid) %>%
    summarize(major.strand=ifelse(sum(sign(qe-qs)*abs(qe-qs))>0, '+', '-'),
              maxQ=max(c(qe, qs)))
  merge(df, major.strand) %>% mutate(qs=ifelse(major.strand=='-', maxQ-qs, qs),
                                     qe=ifelse(major.strand=='-', maxQ-qe, qe),
                                     qid=factor(qid, levels=levels(rid.o$qid)))
}

mumgp.filt.diag = diagMum(mumgp.filt)

ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
  geom_segment() + geom_point(alpha=.4) + theme_bw() + 
  facet_grid(qid~., scales='free', space='free', switch='y') +
  theme(strip.text.y=element_text(angle=180, size=6), strip.background=element_blank(),
        legend.position=c(.99,.01), legend.justification=c(1,0),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1') +
  ggtitle("Pseudoobscura-Male Affinis Diagonalized Matches")

## multiple reference regions
ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
  geom_segment() + geom_point(alpha=.4) + theme_bw() + 
  facet_grid(qid~rid, scales='free', space='free', switch='y') +
  theme(strip.text.y=element_text(angle=180, size=6), strip.background=element_blank(),
        legend.position=c(.99,.01), legend.justification=c(1,0),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1') +
  ggtitle("Pseudoobscura Aligned with Male Affinis")
