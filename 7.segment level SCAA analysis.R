# PREREQUISITIS: load section 4 data 

# "#440154FF", "#39568CFF", "#55C667FF", "#FDE725FF"

# DATA =====

wd <- do.call(rbind, segments)
wd <- wd[which(wd$clonalityCNA>0.5), colnames(wd) %!in% c("sample","cn","countSame")]
wd <- wd[which(wd$n.peaks>0),]
wd <- wd[!duplicated(wd),]
wd <- wd[which(wd$cna!="LOH"),]
wd <- na.omit(wd)

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

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543"))),], aes(y=log(as.numeric(prop.SCAA)), x=(cna), fill=cna)) +
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

ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & wd$n.peaks>1),], aes(y=log(as.numeric(prop.SCAA)), x=(cna), fill=cna)) +
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

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot2::ggplot(data=wd, aes(y=log(as.numeric(prop.SCAA)), x=(cna), fill=cna)) +
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
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
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
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
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

# normalised for size as amps on small
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & (wd$cna=="gain")),], aes(x=as.numeric(aveCNperCNA), y=log(as.numeric(sizenorm.prop.SCAA)))) +
  geom_point(shape=21, size=3, fill="#D95F02") +
  scale_y_continuous(name="Normalised SCAA burden \n(proportion of SCAAs / segment size)", 
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

# investigate
z <- wd[which((wd$patient %!in% c("C516","C542","C543")) & (wd$cna=="gain") & (wd$n.peaks>1) & wd$prop.SCAA>=0.5),]

# PER PATIENT, SCAA PROPORTION LOOKS DIFFERENT ON SMALL SEGMENTS IN DIPLOID/GAIN/LOSS  =====
# scatter (per patient): Y=proportion of SCAAs, X=size of diploid segment 

# diploid
library(ggpubr)
dip.figs <- list()
gain.figs <- list()
loss.figs <- list()
for (p in 1:length(patients)) {
  dip.figs[[p]] <- ggplot(data=wd[which(wd$patient==patients[p] & wd$cna=="diploid"),], aes(x=as.numeric(size)/1000000, y=as.numeric(prop.SCAA))) +
    geom_point(color="#336600", size=3) +
    ggtitle(patients[p]) +
    stat_smooth(method = "glm", color='black') +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
             label.x.npc = "left", size=4, color="black")+
    theme_custom() +
    theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size=18)) 
  
  gain.figs[[p]] <-  ggplot(data=wd[which(wd$patient==patients[p] & wd$cna=="gain"),], aes(x=(as.numeric(size))/1000000, y=as.numeric(prop.SCAA))) +
    geom_point(color="#CC0033", size=3) +
    ggtitle(patients[p]) +
    stat_smooth(method = "glm", color='black') +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
             label.x.npc = "left", size=4, color="black")+
    theme_custom() +
    theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size=18)) 
  
  loss.figs[[p]] <- ggplot(data=wd[which(wd$patient==patients[p] & wd$cna=="loss"),], aes(x=as.numeric(size)/1000000, y=as.numeric(prop.SCAA))) +
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
x <- wd[which(wd$cna=="diploid"),]

d1 <- ggplot(data=x, aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#336600",0.4), size=3, shape=21) + 
  ggtitle("All patients") +
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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*4.19))
annotate_figure(cowplot::plot_grid(d1, d2, d3, ncol = 3),left = text_grob('SCAA burden', size = 28, rot = 90))
dev.off()

# gain 
x <- wd[which(wd$cna=="gain"),]

g1 <- ggplot(data=x, aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#CC0033",0.4), size=3, shape=21) + 
  ggtitle("All patients") +
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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*4.19))
annotate_figure(cowplot::plot_grid(g1, g2, g3, ncol = 3),left = text_grob('SCAA burden', size = 28, rot = 90))
dev.off()

# loss 39568CFF
x <- wd[which(wd$cna=="loss"),]

l1 <- ggplot(data=x, aes(x=log(as.numeric(size)), y=log(as.numeric(prop.SCAA)))) +
  geom_point(fill=alpha("#39568CFF",0.4), size=3, shape=21) + 
  ggtitle("All patients") +
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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank(), axis.title.y = element_blank()) 

jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*4.19))
annotate_figure(cowplot::plot_grid(l1, l2, l3, ncol = 3),left = text_grob('SCAA burden', size = 28, rot = 90))
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
jpeg('tempfig.jpeg', height = (3*37.795*7.5), width = (3*37.795*6))
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & (wd$cna=="diploid" | wd$cna=="gain")),], aes(y=log(as.numeric(prop.SCAA)), x=(sizeGroup), fill=cna)) +
  ggtitle("Diploid V Gain") +
  geom_split_violin(width=0.7, alpha = .4, trim = FALSE, scale = "width") +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitterdodge(jitter.width = 0.3,
                                                                          jitter.height = 0,
                                                                          dodge.width = 0.6,)) +
  geom_boxplot(width = .3, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_y_continuous(name="SCAA burden", 
                     limits = c(-10,4),
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  xlab("Segment size (Mb)") +
  scale_fill_manual(values=c("#1B9E77","#D95F02")) +
  theme_custom() +
  coord_flip() +
  stat_compare_means(aes(group = cna), label="p.format", label.y = 3,  size=8, method = "wilcox.test") +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank())
dev.off()

# diploid v loss
library(ggpubr)
jpeg('tempfig.jpeg', height = (3*37.795*7.5), width = (3*37.795*6))
ggplot2::ggplot(data=wd[which((wd$patient %!in% c("C516","C542","C543")) & (wd$cna=="diploid" | wd$cna=="loss")),], aes(y=log(as.numeric(prop.SCAA)), x=(sizeGroup), fill=cna)) +
  ggtitle("Diploid V Loss") +
  geom_split_violin(width=0.7, alpha = .4, trim = FALSE, scale = "width") +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitterdodge(jitter.width = 0.3,
                                                                          jitter.height = 0,
                                                                          dodge.width = 0.6,)) +
  geom_boxplot(width = .3, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_y_continuous(name="SCAA burden", 
                     limits = c(-10,4),
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  xlab("Segment size (Mb)") +
  scale_fill_manual(values=c("#1B9E77","#7570B3")) +
  theme_custom() +
  coord_flip() +
  stat_compare_means(aes(group = cna), label="p.format", label.y = 3,  size=8, method = "wilcox.test") +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank())
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

# SCAA POPRTION BY SEGMENTO TELOMERE/CENTROMERE DISTRANCE =====


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

set.seed(123)
train.control <- trainControl(method = "LOOCV")
model <- train(log(prop.SCAA) ~ log(size) + cna + clonalityCNA, data = wd2, method = "lm", trControl = train.control)
summary(model)  


# =-=-=-=-=-====
# MONTE CARLO APPROACH --(WIP)-- =====
# scatter: Y=proportion SCAA, X=size of segment the 1Mb loci was sampled from, per patient

# pull the gains with appropriate clonality
gains <- peakperCNA[peakperCNA$CNA.clonality > 0.5 & peakperCNA$cna == "gain",]

# remove segments which have 20% overlap with the regions missing peaks due to 'other reasons' (acro/centro)
peakGapOverlap <- list()
l <- 1
for (i in 1:nrow(gains)) {
  print(i)
  for (j in 1:nrow(peakGaps)) {
    if ( ((Overlap(c(gains$cnaStart[i], gains$cnaStop[i]), c(peakGaps$gapStart[j], peakGaps$gapStop[j])) / gains$size[i]) > 0.2) & gains$chr[i]==peakGaps$chr[j] ) {
      peakGapOverlap[[i]] <- TRUE
      break
    } else { 
      peakGapOverlap[[i]] <- FALSE
    }
  }
}
gains$peakGapOverlap <- unlist(peakGapOverlap)
gains <- gains[gains$peakGapOverlap==FALSE,]

# really small gains cant be subsampled because tehe MC region size would have to become very small, which would miss peaks
gains <- gains[gains$size>25000000,]

# larger gains should have more random MC loci. Each segment should be sampled at least once
gains$weight <- gains$size/sum(gains$size)
gains$n <- ceiling(gains$weight*100000)

# for each gain >1Mb, generate n 1Mb subsampled regions (randomly choose a start site in segment)
MC.regions <- list()
for (i in 1:nrow(gains)) {
  MC.regions[[i]] <- data.frame(patient=gains$patient[i], chr=gains$chr[i], segSize=gains$size[i], segPeakCount=gains$peakCount[i], sampleStart=floor(runif(gains$n[i], min=gains$cnaStart[i], max=gains$cnaStop[i]-25000000)))
}
MC.regions <- data.frame(do.call(rbind, MC.regions))
MC.regions$sampleStop <- MC.regions$sampleStart + 25000000

# pull proportion SCAAs
peaks <- scaaCNA
peaks$chr <- as.numeric(substr(sub("\\:.*", "", peaks$peak),4,nchar(sub("\\:.*", "", peaks$peak))))
peaks$peakStart <- as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", peaks$peak))
peaks$peakStop <- as.numeric(sub(".*-", "", peaks$peak)) 

scaaProp <- list()
for (i in 1:nrow(MC.regions)) {
  print(paste(i,"/",nrow(MC.regions)))
  patient <- MC.regions$patient[i]
  chr <- MC.regions$chr[i]
  start <- MC.regions$sampleStart[i]
  stop <- MC.regions$sampleStop[i]
  
  # pull peaks/SCAA
  peak.wd <- peaks[peaks$patient==patient & peaks$chr==chr & peaks$peakStart>=start & peaks$peakStop<=stop,]
  
  # calculate SCAA prop
  scaaProp[[i]] <- sum(peak.wd$scaaTRUE==1)/nrow(peak.wd)
}

MC.regions$propSCAA <- unlist(scaaProp)
MC.regions <- na.omit(MC.regions)
x <- table(MC.regions$patient)

saveRDS(MC.regions, "~/Documents/SCAA/Data/MC.regions.25Mb.100000x.rds") #25kb x100,000


# plot
jpeg('tempfig.jpeg', width = (3*37.795*7), height = (3*37.795*5))
ggplot(MC.regions, aes(x=(segSize), y=(propSCAA))) + 
  geom_point(size=2, alpha=0.5) +
  ylab("SCAA proportion") +
  xlab("Size of segment sampled") +
  #facet_grid(vars(patient)) +
  #scale_color_manual(values=c("#FF6600", "#55C667FF","#440154FF")) +
  stat_smooth(method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 0, label.y = 'top', size=8) +
  theme_custom() +
  theme(legend.text = element_text(size=20), legend.title = element_blank())
dev.off()


# AT EACH SIZE GROUP, DO GAINS HAVE HIGHER SCAA PROPORTION THAN DIPLOIDS?

