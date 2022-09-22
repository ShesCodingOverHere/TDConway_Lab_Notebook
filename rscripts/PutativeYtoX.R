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

mumgp = readDelta("C:/Users/tdcon/Desktop/New_Female_to_Male.delta")

mumgp %>% head %>% kable

table(mumgp$qid)


filterMum <- function(df, minl=750, flanks=1e4){
  coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
    summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
    ungroup %>% arrange(desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
  merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
    mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

mumgp.filt = filterMum(mumgp, minl=1e4)
mumgp.filt %>% head %>% kable

# PGA_scaffold11__133_contigs__length_7237977 
# PGA_scaffold33__196_contigs__length_9964445
# PGA_scaffold31__9_contigs__length_634783
# Consensus_Consensus_Consensus_disjointig_37_pilon_pilon_pilon_pilon_pilon__unscaffolded
# Consensus_Consensus_Consensus_disjointig_566_pilon_pilon_pilon_pilon_pilon__unscaffolded

# potential y 1
potY1=mumgp.filt[which(mumgp.filt$qid=="PGA_scaffold11__133_contigs__length_7237977"),]
potY1.filt=potY1[which(potY1$rid %in% names(tail(sort(table(potY1$rid)),5))),]

sort(table(potY1$rid))


ggplot(potY1.filt, aes(x=qs, xend=qe, y=rs, yend=re, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(rid~., scales='free', space='free', switch='y') +
  theme_bw() + theme(strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.01), legend.justification=c(1,0),
                     strip.background=element_blank(),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('Assembly') + ylab('reference sequence') + scale_colour_brewer(palette='Set1') +
  ggtitle("Female Affinis to Male Affinis: PGA_scaffold11__133")

head(potY1.filt)

# gives the number of base matches in potY1.filt
tapply(abs(potY1.filt$qs-potY1.filt$qe),potY1.filt$rid,sum)

# gives the mtDNA for potY1.filt
potY1.filt[which(potY1.filt$rid=="mtDNA"),]

# We made a new variable to filter everything on scaffold 11 < 1500000
potY1.notmt=potY1.filt[which(potY1.filt$qs<1500000),]

# then plotted it
ggplot(potY1.notmt, aes(x=qs, xend=qe, y=rs, yend=re, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(rid~., scales='free', space='free', switch='y') +
  theme_bw() + theme(strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.01), legend.justification=c(1,0),
                     strip.background=element_blank(),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('Assembly') + ylab('reference sequence') + scale_colour_brewer(palette='Set1') +
  ggtitle("Female Affinis to Male Affinis: PGA_scaffold11__133")

# made an empty plot
plot(NA,xlim=c(0,7237977),ylim=c(0,2))

# plotted to that empty plot to see the mtDNA
segments(potY1.filt[which(potY1.filt$rid=="mtDNA"),]$qs,1,potY1.filt[which(potY1.filt$rid=="mtDNA"),]$qe,1)



# potential y 2
potY2=mumgp.filt[which(mumgp.filt$qid=="PGA_scaffold33__196_contigs__length_9964445"),]
potY2.filt=potY2[which(potY2$rid %in% names(tail(sort(table(potY2$rid)),5))),]

sort(table(potY2$rid))


ggplot(potY2.filt, aes(x=qs, xend=qe, y=rs, yend=re, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(rid~., scales='free', space='free', switch='y') +
  theme_bw() + theme(strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.01), legend.justification=c(1,0),
                     strip.background=element_blank(),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1') +
  ggtitle("Female Affinis to Male Affinis: PGA_scaffold33__196")




# Potential Y 3
potY3=mumgp.filt[which(mumgp.filt$qid=="PGA_scaffold31__9_contigs__length_634783"),]
potY3.filt=potY3[which(potY3$qid %in% names(tail(sort(table(potY3$qid)),5))),]

sort(table(potY3$rid))


ggplot(potY3.filt, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(rid~., scales='free', space='free', switch='y') +
  theme_bw() + theme(strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.01), legend.justification=c(1,0),
                     strip.background=element_blank(),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1') +
  ggtitle("Female Affinis to Male Affinis: PGA_scaffold31__9")



# Potential Y 4
potY4=mumgp.filt[which(mumgp.filt$qid=="Consensus_Consensus_Consensus_disjointig_37_pilon_pilon_pilon_pilon_pilon__unscaffolded"),]
potY4.filt=potY4[which(potY4$qid %in% names(tail(sort(table(potY4$qid)),5))),]

sort(table(potY4$rid))


ggplot(potY4.filt, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(rid~., scales='free', space='free', switch='y') +
  theme_bw() + theme(strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.01), legend.justification=c(1,0),
                     strip.background=element_blank(),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1') +
  ggtitle("Female Affinis to Male Affinis: disjointig_37")


# potential y 5

potY5=mumgp.filt[which(mumgp.filt$qid=="Consensus_Consensus_Consensus_disjointig_566_pilon_pilon_pilon_pilon_pilon__unscaffolded"),]
potY5.filt=potY5[which(potY5$qid %in% names(tail(sort(table(potY5$qid)),5))),]

sort(table(potY5$rid))


ggplot(potY5.filt, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(rid~., scales='free', space='free', switch='y') +
  theme_bw() + theme(strip.text.y=element_text(angle=180, size=5),
                     legend.position=c(.99,.01), legend.justification=c(1,0),
                     strip.background=element_blank(),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1') +
  ggtitle("Female Affinis to Male Affinis: disjointig_566")

