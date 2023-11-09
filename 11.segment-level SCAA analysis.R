# PREREQUISITIS: load section 4 data 

# "#440154FF", "#39568CFF", "#55C667FF", "#FDE725FF"

# DATA =====
wd <- do.call(rbind, segments)
wd <- wd[c(1,7)]
wd <- wd[!duplicated(wd),]
mean(table(wd$patient))
min(table(wd$patient))
max(table(wd$patient))

wd <- do.call(rbind, segments)
wd <- wd[which(wd$clonalityCNA>0.5), colnames(wd) %!in% c("sample","cn","countSame")]
wd <- wd[which(wd$n.peaks>0),]
wd <- wd[!duplicated(wd),]
wd <- wd[which(wd$cna!="LOH"),]
wd <- na.omit(wd)

wd$coverage <- wd$n.peaks/wd$size

t.test(wd$coverage[which(wd$cna=="diploid")], wd$coverage[which(wd$cna=="gain")])
t.test(wd$coverage[which(wd$cna=="diploid")], wd$coverage[which(wd$cna=="loss")])
t.test(wd$coverage[which(wd$cna=="gain")], wd$coverage[which(wd$cna=="loss")])


# log windows so floor to nearest million
wd$sizeGroup <- "100-150"
wd$sizeGroup[wd$size<=100000000] <- "75-100"
wd$sizeGroup[wd$size<=75000000] <- "50-75"
wd$sizeGroup[wd$size<=50000000] <- "25-50"
wd$sizeGroup[wd$size<=25000000] <- "1-25"
wd$sizeGroup[wd$size<=1000000] <- "0-1"
wd$sizeGroup <- factor(wd$sizeGroup, levels = c("0-1","1-25","25-50","50-75","75-100","100-150"))

# centromere data
library(data.table)
hg38.centromere <- read.csv("~/Documents/SCAA/Data/hg38centromeres.csv")
hg38centromere <- data.frame(aggregate(hg38.centromere$chromStart, list(hg38.centromere$chrom), FUN=min))
colnames(hg38centromere) <- c("chr","centroStart")
hg38centromere$centroStop <- unlist(aggregate(hg38.centromere$chromEnd, list(hg38.centromere$chrom), FUN=max)[2])

# telomere data
library("DescTools")
library(rtracklayer)
hg38.length <- SeqinfoForUCSCGenome("hg38")
hg38.length <- data.frame("chr"=substr(hg38.length@seqnames,4,nchar(hg38.length@seqnames)), "length"=hg38.length@seqlengths)
hg38.length <- hg38.length[c(1:22),]

# assign to each peak
dist.centro <- list()
dist.telo <- list()
for (i in 1:nrow(wd)) {
  print(paste(i,"/",nrow(wd)))
  chr <- paste("chr",wd$chr[i], sep = "")
  start <- as.numeric(wd$cnaStart[i])
  stop <- as.numeric(wd$cnaStop[i])

  # centromere per peak
  cen <- hg38centromere[hg38centromere$chr == chr,]
  dist.centro[[i]] <- min(abs(stop - cen$centroStop), abs(start - cen$centroStart))
  
  # telomere per peak
  tel <- hg38.length[hg38.length$chr == substr(chr, 4, nchar(chr)),]
  dist.telo[[i]] <- min(abs(start - 1), abs(stop - tel$length))
  
}
wd$dist.centro <- unlist(dist.centro)
wd$dist.telo <- unlist(dist.telo)

wd$pga <- scaaPP[match(wd$patient, scaaPP$patient), 8]
wd$coverage <- wd$n.peaks/wd$size
wd$MSI <- scaaPP[match(wd$patient, scaaPP$patient), 7]
wd$sizenorm.prop.SCAA <- wd$prop.SCAA/wd$size



# SCAA BURDEN ON DIPLOID/GAIN/LOSS (VIOLIN) =====

library(ggpubr)

my_comparisons <- list( c("diploid", "gain"), c("diploid", "loss"),c("gain", "loss"))

# exc 3, min 2 peak
#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & wd$n.peaks>1),], aes(y=log(as.numeric(prop.SCAA)), x=(cna), fill=cna)) +
  ggtitle("Exc. C516/C542/C543") +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitterdodge(jitter.width = 0.2,
                                                                          jitter.height = 0.1,
                                                                          dodge.width = 0.6,)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,2.4),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(comparisons = my_comparisons, label="p.format", size=8) +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), axis.title.x = element_blank())
dev.off()


mean(wd$prop.SCAA[which((wd$patient %!in% c("C516","C542","C543")) & wd$n.peaks>1 & wd$cna=="loss")])
mean(wd$prop.SCAA[which((wd$patient %!in% c("C516","C542","C543")) & wd$n.peaks>1 & wd$cna=="gain")])
mean(wd$prop.SCAA[which((wd$patient %!in% c("C516","C542","C543")) & wd$n.peaks>1 & wd$cna=="diploid")])

# all patients, & min 2 peak
#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot2::ggplot(data=wd[which(wd$n.peaks>1),], aes(y=log(as.numeric(prop.SCAA)), x=(cna), fill=cna)) +
  ggtitle("All patients") +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitterdodge(jitter.width = 0.2,
                                                                          jitter.height = 0.1,
                                                                          dodge.width = 0.6,)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,2.4),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(comparisons = my_comparisons, label="p.format", size=8) +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), axis.title.x = element_blank())
dev.off()

# split by MSI
#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & wd$n.peaks>1 &wd$MSI=="MSI"),], aes(y=log(as.numeric(prop.SCAA)), x=(cna), fill=cna)) +
  ggtitle("MSI patients") +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitterdodge(jitter.width = 0.2,
                                                                          jitter.height = 0.1,
                                                                          dodge.width = 0.6,)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,2.4),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(comparisons = my_comparisons, label="p.format", size=8) +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), axis.title.x = element_blank())
dev.off()

#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & wd$n.peaks>1 &wd$MSI=="MSS"),], aes(y=log(as.numeric(prop.SCAA)), x=(cna), fill=cna)) +
  ggtitle("MSS patients") +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitterdodge(jitter.width = 0.2,
                                                                          jitter.height = 0.1,
                                                                          dodge.width = 0.6,)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,2.4),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(comparisons = my_comparisons, label="p.format", size=8) +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), axis.title.x = element_blank())
dev.off()


# SCAA BURDEN ON AMPLIFIED REGIONS =====
#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & (wd$cna=="gain")),], aes(x=as.numeric(aveCNperCNA), y=log(as.numeric(prop.SCAA)))) +
  geom_point(shape=21, size=3, fill="#D95F02") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,2.4),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  xlab("Average amplification of gains") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank())
dev.off()


ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & (wd$cna=="gain") & (wd$n.peaks>1)),], aes(x=as.numeric(aveCNperCNA), y=log(as.numeric(prop.SCAA)))) +
  geom_point(shape=21, size=3, fill="#D95F02") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,2.4),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  xlab("Average amplification of gains") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# normalised for size as amps on small, more than 1 peak
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & (wd$cna=="gain") & wd$n.peaks>1),], aes(x=as.numeric(aveCNperCNA), y=log(as.numeric(sizenorm.prop.SCAA)))) +
  geom_point(size=3, alpha = .5, shape=21, fill="#D95F02") +
  scale_y_continuous(name="Normalised SCAA burden \n(SCAA burden / segment size)", 
                     breaks = c(-25,-20,-15,-10,-5),
                     labels = c(round(exp(-25),11),round(exp(-20),9),round(exp(-15),7),round(exp(-10),5),round(exp(-5),4))) +
  xlab("Average amplification of gains") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
dev.off()



# if SCAAs are controlling exp of amps, what genes are in amps
z <- wd[which((wd$patient %!in% c("C516","C542","C543")) & (wd$cna=="gain") & (wd$prop.SCAA>0) & (wd$n.peaks>1) & wd$aveCNperCNA>=5 & wd$sizeGroup=="0-1"),]
output <- list()
l <- 1
for (i in 1:nrow(z)) {
  patient <- z$patient[i]
  chr <- z$chr[i]
  cnaStart <- z$cnaStart[i]
  cnaStop <- z$cnaStop[i]
  aveCNperCNA <- z$aveCNperCNA[i]
  n.peaks <- z$n.peaks[i]
  n.SCAA <- z$n.SCAA[i]
  prop.SCAA <- z$prop.SCAA[i]
  sizenorm.prop.SCAA <- z$sizenorm.prop.SCAA[i]
  symbols <- unique(scaa.long$symbol[which(scaa.long$chr==z$chr[i] & scaa.long$cnaStart==z$cnaStart[i] & scaa.long$cnaStop==z$cnaStop[i])])
  for (j in 1:length(symbols)) {
    output[[l]] <- setNames(c(patient, chr, cnaStart, cnaStop, aveCNperCNA, n.peaks,n.SCAA,prop.SCAA, sizenorm.prop.SCAA, symbols[j]),
                            c("patient","chr","cnaStart","cnaStop","aveCNperCNA","n.peaks","n.SCAA","prop.SCAA","sizenorm.prop.SCAA", "symbol"))
    l <- l + 1
  }
}
y <- data.frame(do.call(rbind,output))
y$n.peaks <- as.numeric(y$n.peaks)
y$n.SCAA <- as.numeric(y$n.SCAA)
y$prop.SCAA <- as.numeric(y$prop.SCAA)
y$sizenorm.prop.SCAA <- as.numeric(y$sizenorm.prop.SCAA)
y <- y[which(y$prop.SCAA>0.2),]


x <- z
x$bin <- 1:nrow(z)
z <- x[,c(23,3,4,5)]
colnames(z) <- c("bin","chr","start","stop")
GOlist.hg38 <- read.delim("~/Documents/SCAA/Data/hg38.1_all_gene_GO_annotations.txt")
GOlist.hg38 <- GOlist.hg38[c(8,11,5,3,4,1)]
GOlist.hg38 <- GOlist.hg38[!duplicated(GOlist.hg38),]
geneAnno_amp <- gene_anno(z, GOlist.hg38)
geneAnno_amp$start <- z[match(geneAnno_amp$bin, z$bin), 3]
geneAnno_amp$stop <- z[match(geneAnno_amp$bin, z$bin), 4]
geneAnno_amp$symbol <- mapIds(org.Hs.eg.db, keys = as.character(geneAnno_amp$entrezgene),
                                  column = c('SYMBOL'), keytype = 'ENTREZID')

output <- list()
for (i in 1:nrow(x)) {
  y <- geneAnno_amp[which(geneAnno_amp$start==x$cnaStart[i]),]
  if(nrow(y)==0) {
    output[[i]] <- NA
    next
  }
  output[[i]] <- paste0(c(y$symbol), collapse=", ")
  #output[[i]] <- y$symbol
}
x$geneMatch <- unlist(output)
x <- x[!is.na(x$geneMatch),]
x <- x[which(x$prop.SCAA>0),]

table(unlist(output))

z <- scaa.long[which(scaa.long$cnaStart==127191793),]

# PER PATIENT, SCAA PROPORTION LOOKS DIFFERENT ON SMALL SEGMENTS IN DIPLOID/GAIN/LOSS  =====
# scatter (per patient): Y=proportion of SCAAs, X=size of diploid segment 

# diploid
library(ggpubr)
dip.figs <- list()
gain.figs <- list()
loss.figs <- list()
for (p in 1:length(patients)) {
  dip.figs[[p]] <- ggplot(data=wd[which(wd$patient==patients[p] & wd$cna=="diploid" & wd$n.peaks>1),], aes(x=as.numeric(size)/1000000, y=as.numeric(prop.SCAA))) +
    geom_point(color="#336600", size=3) +
    ggtitle(patients[p]) +
    stat_smooth(method = "glm", color='black') +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
             label.x.npc = "left", size=4, color="black")+
    theme_custom() +
    theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size=18)) 
  
  gain.figs[[p]] <-  ggplot(data=wd[which(wd$patient==patients[p] & wd$cna=="gain" & wd$n.peaks>1),], aes(x=(as.numeric(size))/1000000, y=as.numeric(prop.SCAA))) +
    geom_point(color="#CC0033", size=3) +
    ggtitle(patients[p]) +
    stat_smooth(method = "glm", color='black') +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
             label.x.npc = "left", size=4, color="black")+
    theme_custom() +
    theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size=18)) 
  
  loss.figs[[p]] <- ggplot(data=wd[which(wd$patient==patients[p] & wd$cna=="loss" & wd$n.peaks>1),], aes(x=as.numeric(size)/1000000, y=as.numeric(prop.SCAA))) +
    geom_point(color="#39568CFF", size=3) +
    ggtitle(patients[p]) +
    stat_smooth(method = "glm", color='black') +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
             label.x.npc = "left", size=4, color="black")+
    theme_custom() +
    theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size=18)) 
}

dip.figs <- dip.figs[lengths(dip.figs) != 0]
gain.figs <- gain.figs[lengths(gain.figs) != 0]
loss.figs <- loss.figs[lengths(loss.figs) != 0]


jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*8.76))
annotate_figure(cowplot::plot_grid(plotlist = dip.figs, ncol=6), 
                left = text_grob("Proportion of SCAAs in diploid segments", rot = 90, size = 28, vjust = 0.5),
                bottom = text_grob("Mb of diploid segment", size = 28, vjust = 0))
dev.off()

jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*8.76))
annotate_figure(cowplot::plot_grid(plotlist = gain.figs, ncol=6), 
                left = text_grob("Proportion of SCAAs in gain segments", rot = 90, size = 28, vjust = 0.5),
                bottom = text_grob("Mb of gain segment", size = 28, vjust = 0))
dev.off()

jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*8.76))
annotate_figure(cowplot::plot_grid(plotlist = loss.figs, ncol=6), 
                left = text_grob("Proportion of SCAAs in loss segments", rot = 90, size = 28, vjust = 0.5),
                bottom = text_grob("Mb of loss segment", size = 28, vjust = 0))
dev.off()

# ACROSS PATIENT, SCAA PROPORTION LOOKS DIFFERENT ON SMALL SEGMENTS IN DIPLOID/GAIN =====
# scatter for all, and with the 3 removed. X=size, Y=prop SCAA

# diploid
x <- wd[which(wd$cna=="diploid" & wd$n.peaks>1),]
jpeg('tempfig.jpeg', width = (50), height = (80), units = "cm", res = 300)
ggplot(data=x, aes(x=(as.numeric(size)), y=(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#336600",0.4), size=3, shape=21) + 
  facet_wrap(~patient, ncol = 4, scales="free_y") +
  scale_y_continuous(name="SCAA burden") +
  scale_x_continuous(name="Size of diploid") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(),
        strip.text = element_text(size=24)) 
dev.off()

d1 <- ggplot(data=x, aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#336600",0.4), size=3, shape=21) + 
  ggtitle("All patients") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-8,-6,-4,-2,0),
                     limits = c(-8,0),
                     labels = c(round(exp(-8),5),round(exp(-6),4),round(exp(-4),2),round(exp(-2),2),round(exp(0),2))) +
  scale_x_continuous(name="Size of diploid", 
                     limits = c(7.5,20),
                     breaks = c(10,12.5,15,17.5),
                     labels = c(signif(round(exp(10),0),1),signif(round(exp(12.5),0),1), signif(round(exp(15),0),1), signif(round(exp(17.5),0),1))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(),
        strip.text = element_text(size=24)) 

d2 <- ggplot(data=x[-which(x$patient %in% c("C543","C542","C516")),], aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#336600",0.4), size=3, shape=21) +
  ggtitle("exc C516/542/543") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-8,-6,-4,-2,0),
                     limits = c(-8,0),
                     labels = c(round(exp(-8),5),round(exp(-6),4),round(exp(-4),2),round(exp(-2),2),round(exp(0),2))) +
  scale_x_continuous(name="Size of diploid", 
                     limits = c(6,19),
                     breaks = c(6,10,14,18),
                     labels = c(signif(round(exp(6),0),1), signif(round(exp(10),0),1), signif(round(exp(12),0),1), signif(round(exp(18),0),1))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

d3 <- ggplot(data=x[which(x$patient %in% c("C543","C542","C516")),], aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#336600",0.4), size=3, shape=21) +
  ggtitle("C516/542/543") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-3,-2,-1,0),
                     limits = c(-3,0),
                     labels = c(round(exp(-3),4),round(exp(-2),2),round(exp(-1),2),round(exp(0),2))) +
  scale_x_continuous(name="Size of diploid", 
                     limits = c(10,19),
                     breaks = c(10,14,18),
                     labels = c(signif(round(exp(10),0),1), signif(round(exp(12),0),1), signif(round(exp(18),0),1))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

#jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (50), height = (15), units = "cm", res = 300)
annotate_figure(cowplot::plot_grid(d1, d2, d3, ncol = 3, labels = c("A","B","C"), label_size = 34),left = text_grob('SCAA burden', size = 28, rot = 90))
dev.off()

# gain 
x <- wd[which(wd$cna=="gain" & wd$n.peaks>1),]

g1 <- ggplot(data=x, aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#CC0033",0.4), size=3, shape=21) + 
  ggtitle("All patients") +
  facet_wrap(~patient, ncol = 4) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-8,-6,-4,-2,0),
                     limits = c(-8,0),
                     labels = c(round(exp(-8),5),round(exp(-6),4),round(exp(-4),2),round(exp(-2),2),round(exp(0),2))) +
  scale_x_continuous(name="Size of gain", 
                     limits = c(5.5,19),
                     breaks = c(6,10,14,18),
                     labels = c(signif(round(exp(6),0),1), signif(round(exp(10),0),1), signif(round(exp(12),0),1), signif(round(exp(18),0),1))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

g2 <- ggplot(data=x[-which(x$patient %in% c("C543","C542","C516")),], aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#CC0033",0.4), size=3, shape=21) +
  ggtitle("exc C516/542/543") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-8,-6,-4,-2,0),
                     limits = c(-8,0),
                     labels = c(round(exp(-8),5),round(exp(-6),4),round(exp(-4),2),round(exp(-2),2),round(exp(0),2))) +
  scale_x_continuous(name="Size of gain", 
                     limits = c(5.5,19),
                     breaks = c(6,10,14,18),
                     labels = c(signif(round(exp(6),0),1), signif(round(exp(10),0),1), signif(round(exp(12),0),1), signif(round(exp(18),0),1))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

g3 <- ggplot(data=x[which(x$patient %in% c("C543","C542","C516")),], aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#CC0033",0.4), size=3, shape=21) +
  ggtitle("C516/542/543") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-3,-2,-1,0),
                     limits = c(-3,0),
                     labels = c(round(exp(-3),4),round(exp(-2),2),round(exp(-1),2),round(exp(0),2))) +
  scale_x_continuous(name="Size of gain", 
                     limits = c(10,19),
                     breaks = c(10,14,18),
                     labels = c(signif(round(exp(10),0),1), signif(round(exp(12),0),1), signif(round(exp(18),0),1))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

#jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (50), height = (15), units = "cm", res = 300)
annotate_figure(cowplot::plot_grid(g1, g2, g3, ncol = 3, labels = c("D","E","F"), label_size = 34),left = text_grob('SCAA burden', size = 28, rot = 90))
dev.off()

# loss 39568CFF
x <- wd[which(wd$cna=="loss" & wd$n.peaks>1),]

l1 <- ggplot(data=x, aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#39568CFF",0.4), size=3, shape=21) + 
  ggtitle("All patients") +
  facet_wrap(~patient, ncol = 4) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-8,-6,-4,-2,0),
                     limits = c(-8,0),
                     labels = c(round(exp(-8),5),round(exp(-6),4),round(exp(-4),2),round(exp(-2),2),round(exp(0),2))) +
  scale_x_continuous(name="Size of loss", 
                     limits = c(10,19),
                     breaks = c(10,14,18),
                     labels = c(signif(round(exp(10),0),1), signif(round(exp(12),0),1), signif(round(exp(18),0),1))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

l2 <- ggplot(data=x[-which(x$patient %in% c("C543","C542","C516")),], aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#39568CFF",0.4), size=3, shape=21) +
  ggtitle("exc C516/542/543") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-8,-6,-4,-2,0),
                     limits = c(-8,0),
                     labels = c(round(exp(-8),5),round(exp(-6),4),round(exp(-4),2),round(exp(-2),2),round(exp(0),2))) +
  scale_x_continuous(name="Size of loss", 
                     limits = c(10,19),
                     breaks = c(10,14,18),
                     labels = c(signif(round(exp(10),0),1), signif(round(exp(12),0),1), signif(round(exp(18),0),1))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

l3 <- ggplot(data=x[which(x$patient %in% c("C543","C542","C516")),], aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#39568CFF",0.4), size=3, shape=21) +
  ggtitle("C516/542/543") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-1.5,-1.2,-0.9),
                     limits = c(-1.8,-0.5),
                     labels = c(round(exp(-1.5),2),round(exp(-1.2),2),round(exp(-0.9),2))) +
  scale_x_continuous(name="Size of loss", 
                     limits = c(10,18),
                     breaks = c(12,14,16),
                     labels = c(signif(round(exp(12),0),1), signif(round(exp(14),0),1), signif(round(exp(16),0),1))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

#jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (50), height = (15), units = "cm", res = 300)
annotate_figure(cowplot::plot_grid(l1, l2, l3, ncol = 3, labels = c("G","H","I"), label_size = 34),left = text_grob('SCAA burden', size = 28, rot = 90))
dev.off()

jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*4.19))
annotate_figure(cowplot::plot_grid(d2, g2, l2, ncol = 3),left = text_grob('SCAA burden', size = 28, rot = 90))
dev.off()



# ACROSS PATIENT, SCAA PROPORTION LOOKS DIFFERENT ON SMALL SEGMENTS IN DIPLOID/GAIN (VIOLIN) =====


probSCAA.dip$size <- as.numeric(probSCAA.dip$size)
max(probSCAA.dip$size)
probSCAA.dip$group <- "100-150"
probSCAA.dip$group[probSCAA.dip$size<=100000000] <- "75-100"
probSCAA.dip$group[probSCAA.dip$size<=75000000] <- "50-75"
probSCAA.dip$group[probSCAA.dip$size<=50000000] <- "25-50"
probSCAA.dip$group[probSCAA.dip$size<=25000000] <- "1-25"
probSCAA.dip$group[probSCAA.dip$size<=1000000] <- "0-1"
probSCAA.dip$group <- factor(probSCAA.dip$group, levels = c("0-1","1-25","25-50","50-75","75-100","100-150"))

probSCAA.gain$size <- as.numeric(probSCAA.gain$size)
probSCAA.gain$group <- "100-150"
probSCAA.gain$group[probSCAA.gain$size<=100000000] <- "75-100"
probSCAA.gain$group[probSCAA.gain$size<=75000000] <- "50-75"
probSCAA.gain$group[probSCAA.gain$size<=50000000] <- "25-50"
probSCAA.gain$group[probSCAA.gain$size<=25000000] <- "1-25"
probSCAA.gain$group[probSCAA.gain$size<=1000000] <- "0-1"
probSCAA.gain$group <- factor(probSCAA.gain$group, levels = c("0-1","1-25","25-50","50-75","75-100","100-150"))

probSCAA.loss$size <- as.numeric(probSCAA.loss$size)
max(probSCAA.loss$size)
probSCAA.loss$group <- "100-150"
probSCAA.loss$group[probSCAA.loss$size<=100000000] <- "75-100"
probSCAA.loss$group[probSCAA.loss$size<=75000000] <- "50-75"
probSCAA.loss$group[probSCAA.loss$size<=50000000] <- "25-50"
probSCAA.loss$group[probSCAA.loss$size<=25000000] <- "1-25"
probSCAA.loss$group[probSCAA.loss$size<=1000000] <- "0-1"
probSCAA.loss$group <- factor(probSCAA.loss$group, levels = c("0-1","1-25","25-50","50-75","75-100","100-150"))

my_comparisons <- list( c("0-1", "1-25"), c("0-1", "25-50"), c("0-1", "50-75"), c("0-1", "75-100"), c("0-1", "100-150"),
                        c("1-25", "25-50"), c("1-25", "50-75"), c("1-25", "75-100"), c("1-25", "100-150"),
                        c("25-50", "50-75"), c("25-50", "75-100"), c("25-50", "100-150"),
                        c("50-75", "75-100"), c("50-75", "100-150"),
                        c("75-100", "100-150"))

fig1 <- ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & wd$cna=="diploid"),], aes(y=log(as.numeric(prop.SCAA)), x=(sizeGroup), fill=as.factor(sizeGroup))) +
  geom_violin(aes(fill=as.factor(sizeGroup)), trim = F, alpha = .4, position=position_dodge(0.2)) +
  geom_jitter(size=2, shape=21, alpha=.4, position=position_jitterdodge(jitter.width = 0.7, jitter.height = 0, dodge.width = 0.2)) +
  geom_boxplot(width = 0.3, alpha = .8,  show.legend = FALSE,position=position_dodge(0.2)) +
  ggtitle("Diploid") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  xlab("Segment size (Mb)") +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size=6) +
  theme(legend.position = "none",axis.text.x = element_blank(), axis.title = element_blank())
  
fig2 <- ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & wd$cna=="gain"),], aes(y=log(as.numeric(prop.SCAA)), x=(sizeGroup), fill=as.factor(sizeGroup))) +
  geom_violin(aes(fill=as.factor(sizeGroup)), trim = F, alpha = .4, position=position_dodge(0.2)) +
  geom_jitter(size=2, shape=21, alpha=.4, position=position_jitterdodge(jitter.width = 0.7, jitter.height = 0, dodge.width = 0.2)) +
  geom_boxplot(width = 0.3, alpha = .8,  show.legend = FALSE,position=position_dodge(0.2)) +
  ggtitle("Gain") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  xlab("Segment size (Mb)") +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size=6) +
  theme(legend.position = "none",axis.text.x = element_blank(), axis.title = element_blank())

fig3 <- ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & wd$cna=="loss"),], aes(y=log(as.numeric(prop.SCAA)), x=(sizeGroup), fill=as.factor(sizeGroup))) +
  geom_violin(aes(fill=as.factor(sizeGroup)), trim = F, alpha = .4, position=position_dodge(0.2)) +
  geom_jitter(size=2, shape=21, alpha=.4, position=position_jitterdodge(jitter.width = 0.7, jitter.height = 0, dodge.width = 0.2)) +
  geom_boxplot(width = 0.3, alpha = .8,  show.legend = FALSE,position=position_dodge(0.2)) +
  ggtitle("Loss") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  xlab("Segment size (Mb)") +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size=6) +
  theme(legend.position = "none", axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), axis.title = element_blank())


jpeg('tempfig.jpeg', width = (3*37.795*10), height = (3*37.795*11))
annotate_figure(cowplot::plot_grid(fig1, fig2, fig3, ncol=1, rel_heights = c(1,1,1.3), align = "v"), 
                left = text_grob("Proportion of SCAAs (log)", rot = 90, size = 28, vjust = 0.5),
                bottom = text_grob("Segment size (Mb)", size = 28, vjust = 0))
dev.off()



# SCAA ENRICHMENT ON GAIN V DIPLOID GROUPED BY SIZE (PAIRED VIOLIN) =====
wd$sizeGroup <- factor(wd$sizeGroup, levels = c("100-150","75-100","50-75","25-50","1-25","0-1"))

# diploid v gain
library(ggpubr)
#jpeg('tempfig.jpeg', height = (3*37.795*7.5), width = (3*37.795*6))
jpeg('tempfig.jpeg', height = (15), width = (40), units = "cm", res = 300)
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & (wd$cna=="diploid" | wd$cna=="gain")),], aes(y=log(as.numeric(prop.SCAA)), x=(sizeGroup), fill=cna)) +
  #ggtitle("Diploid V Gain") +
  geom_split_violin(width=0.7, alpha = .4, trim = FALSE, scale = "width") +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitterdodge(jitter.width = 0.3,
                                                                          jitter.height = 0,
                                                                          dodge.width = 0.6,)) +
  geom_boxplot(width = .3, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_y_continuous(name="SCAA burden", 
                     limits = c(-10,2.3),
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  xlab("Segment size (Mb)") +
  scale_fill_manual(values=c("#1B9E77","#D95F02")) +
  theme_custom() +
  #coord_flip() +
  stat_compare_means(aes(group = cna), label="p.format", label.y = 2,  size=9, method = "wilcox.test") +
  theme(legend.position = "right", legend.text = element_text(size=28), legend.title = element_blank())
dev.off()

# diploid v loss
library(ggpubr)
#jpeg('tempfig.jpeg', height = (3*37.795*7.5), width = (3*37.795*6))
jpeg('tempfig.jpeg', height = (15), width = (40), units = "cm", res = 300)
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & (wd$cna=="diploid" | wd$cna=="loss")),], aes(y=log(as.numeric(prop.SCAA)), x=(sizeGroup), fill=cna)) +
  #ggtitle("Diploid V Loss") +
  geom_split_violin(width=0.7, alpha = .4, trim = FALSE, scale = "width") +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitterdodge(jitter.width = 0.3,
                                                                          jitter.height = 0,
                                                                          dodge.width = 0.6,)) +
  geom_boxplot(width = .3, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_y_continuous(name="SCAA burden", 
                     limits = c(-10,2.3),
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  xlab("Segment size (Mb)") +
  scale_fill_manual(values=c("#1B9E77","#7570B3")) +
  theme_custom() +
  #coord_flip() +
  stat_compare_means(aes(group = cna), label="p.format", label.y = 2,  size=9, method = "wilcox.test") +
  theme(legend.position = "right", legend.text = element_text(size=28), legend.title = element_blank())
dev.off()

# SCAA POPRTION BY CLONALITY =====

jpeg('tempfig.jpeg', width = (3*37.795*7), height = (3*37.795*5))
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) ),], aes(y=log(as.numeric(prop.SCAA)), x=log(clonalityCNA))) +
  geom_point() +
  scale_fill_brewer(palette = "Dark2") +
  stat_smooth(color='black') +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x.npc = "left", size=4, color="black")+
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
dev.off()

# SCAA POPRTION BY SEGMENT TELOMERE/CENTROMERE DISTRANCE =====


jpeg('tempfig.jpeg', width = (3*37.795*7), height = (3*37.795*5))
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) ),], aes(y=log(as.numeric(prop.SCAA)), x=log(dist.telo))) +
  geom_point() +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme_custom() +
  stat_smooth(color='black') +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x.npc = "left", size=4, color="black")+
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
dev.off()

jpeg('tempfig.jpeg', width = (3*37.795*7), height = (3*37.795*5))
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) ),], aes(y=log(as.numeric(prop.SCAA)), x=log(dist.centro))) +
  geom_point() +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme_custom() +
  stat_smooth(color='black') +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x.npc = "left", size=4, color="black")+
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
dev.off()

# MODELLING NUMBER OF SCAAS ====


model4 <- lm(prop.SCAA ~ log(size) + cna, data = wd[wd$patient %!in% c("C516","C542","C543"),])
summary(model4)

model5 <- lm(prop.SCAA ~ log(size) + cna + log(size)*cna, data = wd[wd$patient %!in% c("C516","C542","C543"),])
summary(model5)

model6 <- lm(prop.SCAA ~ log(size) + cna + clonalityCNA, data = wd[wd$patient %!in% c("C516","C542","C543"),])
summary(model6)

model7 <- lm(prop.SCAA ~ log(size) + cna + clonalityCNA + dist.centro, data = wd[wd$patient %!in% c("C516","C542","C543"),])
summary(model7)

model8 <- lm(prop.SCAA ~ log(size) + cna + clonalityCNA, data = wd[wd$patient %!in% c("C516","C542","C543"),])
summary(model8)

library(caret)
set.seed(123)
wd2 <- wd[which(wd$patient %!in% c("C516","C542","C543")),]
wd2 <- wd2[which(wd2$prop.SCAA!=0),]
cut <- createDataPartition(y = wd2$prop.SCAA, p = .75, list = FALSE)
training <- wd2[ cut,]
testing  <- wd2[-cut,]

model <- lm(log(prop.SCAA) ~ log(size) + cna + clonalityCNA, data = training)
summary(model)
predictions <- data.frame(prediction=predict(model, newdata = testing), actual=log(testing$prop.SCAA))
ggplot2::ggplot(data=predictions, aes(y=prediction, x=actual)) +
  geom_point() +
  scale_fill_brewer(palette = "Dark2") +
  stat_smooth(color='black') +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
           label.x.npc = "left", size=4, color="black")+
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

wd2 <- wd[which(wd$patient %!in% c("C516","C542","C543")),]
wd2 <- wd2[which(wd2$prop.SCAA!=0),]
set.seed(123)
train.control <- trainControl(method = "LOOCV")
model <- train(log(prop.SCAA) ~ log(size) + cna + clonalityCNA +clonalityCNA*log(size), data = wd2, method = "lm", trControl = train.control)
summary(model)  


# ARE SCAAS AT END OF LARGE SEGMENTS =====

peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")]

scaa.bins <- data.frame(chr = substr(sub("\\:.*", "", rownames(scaa)),4,nchar(sub("\\:.*", "", rownames(scaa)))),
                        peakStart = as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", rownames(scaa))),
                        peakStop = as.numeric(sub(".*-", "", rownames(scaa))) )

z <- wd[which(wd$sizeGroup=="25-50"),]
prop.SCAA <- list()
l <- 1
for (i in 1:nrow(z)) {
  patient <- z$patient[i]
  cna <- z$cna[i]
  segID <- z$segID[i]
  chr <- z$chr[i]
  cnaStart <- z$cnaStart[i]
  cnaStop <- z$cnaStop[i]
  
  # create buffer zeon on each side
  zoneStart <- cnaStart-500000
  if (zoneStart<1) {
    zoneStart <- 1
  }
  zoneStop <- cnaStop+500000
  if (zoneStop>hg38.length$length[hg38.length$chr==chr]) {
    zoneStop <- hg38.length$length[hg38.length$chr==chr]
  }
  
  # split into roughly windows
  zoneSize <- zoneStop-zoneStart
  n.windows <- ceiling(zoneSize / 1000000)
  windowSize <- floor(zoneSize / n.windows)
  windows <- seq(zoneStart,zoneStop,by=windowSize)
  
  # for each window...
  windowMidpoint <- list()
  for (j in 2:length(windows)) {
    print(paste(i,"/",nrow(z),"-",j,"/",length(windows)))
    windowStart <- windows[j-1]
    windowStop <- windows[j]
    windowMidpoint <- j/length(windows)
    
    # ...pull scaa proportion
    x <- scaa[which(scaa.bins$chr==chr & (scaa.bins$peakStart>=(windowStart-250)) & (scaa.bins$peakStop<=(windowStop+250)) ), substr(colnames(scaa),1,4) == patient]
    prop.SCAA[[l]] <- setNames(c(patient, cna, segID, windowMidpoint, (sum(x==1) / length(x))), 
                               c("patient","cna", "segID", "location", "prop.SCAA"))
    l <- l + 1
  }
  
  
}

x <- data.frame(do.call(rbind, prop.SCAA))
x$location <- as.numeric(x$location)
x$prop.SCAA <- as.numeric(x$prop.SCAA)
ggplot(data=x[which(x$patient==unique(x$patient)[1]),], aes(x=location, y=(prop.SCAA), group=segID, colour=as.factor(segID))) + 
  geom_line() + 
  geom_point()

# CHECK THAT GAINS ARENT ASS WITH OPEN SCAAS (IE COVERAGE ISSUE, CAN ONLY CHECK FOR SCAAS IN PATIENTS WITH RNA DATA) ====
openclose <- readxl::read_excel("~/Documents/SCAA/Data/41586_2022_5202_MOESM9_ESM.xlsx")

# cant tell which patients had open closed, so only consider where 100%
openclose <- openclose[which(openclose$peak.gain==0 | openclose$peak.loss==0),]
openclose <- openclose[which(openclose$peak.recurrence>10),]

x <- scaa.long[which(scaa.long$peak %in% openclose$peak.peak),]
x$peak_event.type <- openclose[match(x$peak, openclose$peak.peak), 7]
x <- x[which(x$value==1),]
x <- x[which(x$mainClonality>0.5),]

y <- data.frame(diploid=c(nrow(x[which(x$cna=="diploid" & x$peak_event.type$peak.event_type=="gain"),]), nrow(x[which(x$cna=="diploid" & x$peak_event.type$peak.event_type=="loss"),])),
                gain=c(nrow(x[which(x$cna=="gain" & x$peak_event.type$peak.event_type=="gain"),]), nrow(x[which(x$cna=="gain" & x$peak_event.type$peak.event_type=="loss"),])),
                loss=c(nrow(x[which(x$cna=="loss" & x$peak_event.type$peak.event_type=="gain"),]), nrow(x[which(x$cna=="loss" & x$peak_event.type$peak.event_type=="loss"),])))

rownames(y) <- c("open","close")

chisq <- chisq.test(y)
enrich <- chisq$observed/chisq$expected


# pull cn norm ATAC data for each SCAA. per patient, does increase in openness correlate with clonality of gain?
# use extent of gain to give clairty on corr

vst <- assay(deseqVst.dWGS.cohort)
#x <- scaa.long[which(scaa.long$value==1),]
x <- scaa.long
#x <- x[which(x$peak %in% openclose$peak.peak),]
outut <- list()
for (i in 1:nrow(x)) {
  print(paste(i,"/",nrow(x)))
  patient <- x$patient[i]
  y <- vst[which(rownames(vst)==x$peak[i]), grep(patient, colnames(vst))]
  if(length(y)==0) {
    output[[i]] <- setNames(c(x$peak[i],x$value[i],patient, NA, NA),
                            c("peak","scaa","patient","p","estimate"))
    next
  }
  z <- dwgs.perpatient[[patient]]
  z <- z[which(z$chromosome==x$chr[i] & z$start==x$cnaStart[i] & z$stop==x$cnaStop[i]),-c(1:3)]
  z <- z[names(z) %in% names(y)]
  if(sum(!is.na(z))<=1) {
    output[[i]] <- setNames(c(x$peak[i],x$value[i],patient,NA, NA),
                            c("peak","scaa","patient","p","estimate"))
    next
  }
  correlation <- cor.test(as.numeric(z), as.numeric(y))
  output[[i]] <- setNames(c(x$peak[i],x$value[i],patient,correlation$p.value, correlation$estimate),
                          c("peak","scaa","patient","p","estimate"))
}

a <- data.frame(do.call(rbind, output))
a$p <- as.numeric(a$p)
a$p.adj <- p.adjust(a$p, method = "fdr")
a$sig <- ifelse(a$p.adj<0.01, TRUE, FALSE)
a$sigSameDir <- ifelse(a$sig==TRUE & a$estimate>0, TRUE, FALSE)
sum(a$sigSameDir==TRUE, na.rm = T)/nrow(a)

# per patient/peak approach
vst <- assay(deseqVst.dWGS.cohort)
x <- scaa.long
output <- list()
for (i in 1:nrow(x)) {
  print(paste(i,"/",nrow(x)))
  patient <- x$patient[i]
  y <- vst[which(rownames(vst)==x$peak[i]), grep(patient, colnames(vst))]
  if(length(y)==0) {
    next
  }
  z <- dwgs.perpatient[[patient]]
  z <- z[which(z$chromosome==x$chr[i] & z$start==x$cnaStart[i] & z$stop==x$cnaStop[i]),-c(1:3)]
  z <- z[names(z) %in% names(y)]
  if(sum(!is.na(z))<=1) {
    next
  }
  wd <- data.frame(t(rbind(y,z)))
  colnames(wd) <- c("peakSignal","cn")
  wd$sample <- rownames(wd)
  wd$patient <- patient
  wd$peak <- x$peak[i]
  output[[i]] <- wd
}
b <- data.frame(do.call(rbind, output))
b$patient <- substr(b$sample,1,4)
cor.test(as.numeric(b$peakSignal), as.numeric(b$cn))
mod <- lm(peakSignal ~ cn, data=b[which(b$patient=="C555"),])
mod <- lm(peakSignal ~ cn, data=b[which(b$peak=="chr1:904506-905006"),])
mod <- lm(peakSignal ~ cn, data=b)
summary(mod)

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data=b[which(b$patient=="C524"),], aes(x=as.numeric(peakSignal), y=as.numeric(cn))) +
  geom_point(size=3, shape=21, fill=alpha("#333399",0.85)) +
  ggtitle(patients[p]) +
  stat_smooth(method = "glm", color='black') +
  ylab("Copy number") +
  xlab("Peak signal") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 
dev.off()

