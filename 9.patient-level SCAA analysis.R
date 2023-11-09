
# PCA TO DETECT OUTLIERS =====
pca_input <- t(scaa)
pca <- prcomp(pca_input)
percentVar <- round(100 * summary(pca)$importance[2,])
df <- cbind(data.frame(patient=colnames(scaa)), pca$x)

#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(df, aes(x=PC1, y=PC2, fill=patient)) + 
  geom_point(size=7, shape=21, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_text_repel(aes(label=substr(patient,1,4)), size=8,box.padding = 1) +
  coord_fixed() +
  theme_custom() +
  theme(legend.position = "none", plot.margin = margin(0.5, 0,0,0, "cm"))
dev.off()

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot2::ggplot(scaaPP, aes(x=reorder(patient, pga_recentre), group=1, fill=patient)) +
  #geom_bar(aes(y=(propSCAA)), fill=alpha("#FF9933",0.6), color="#FF9933", stat="identity") +
  geom_bar(aes(y=(propSCAA)), stat="identity", alpha=0.85, color="black") +
  scale_y_continuous(name="SCAA burden") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), axis.title.x = element_blank(), axis.title.y = element_text(vjust = 0),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "none", 
        strip.text = element_text(size=28), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 3, l = 0))
dev.off()


# SCAA PROPORTION AND WGD/MSI =====
# scatter: Y=proportion of SCAAs, X=ploidy mean per patient
# finding: WGD patients have similar numbers of SCAAS per patient, background SCAA prop is highly patient specific
pvWGD <- round(wilcox.test(as.numeric(scaaPP$propSCAA[which(scaaPP$WGD==">2.5")]), as.numeric(scaaPP$propSCAA[which(scaaPP$WGD=="<2.5")]), alternative = "two.sided")$p.value,3)
pvMSI <- round(wilcox.test(as.numeric(scaaPP$propSCAA[which(scaaPP$type=="MSI")]), as.numeric(scaaPP$propSCAA[which(scaaPP$type=="MSS")]), alternative = "two.sided")$p.value,3)

library(scales)
p1 <- ggplot2::ggplot(data=scaaPP, aes(y=log(propSCAA), x=as.factor(WGD), fill=as.factor(WGD))) +
  geom_violin(trim = F, alpha=0.4) +
  geom_boxplot(width=0.1, outlier.size = 0) +
  geom_jitter(size=5, shape=21, position=position_jitter(0.2, seed=1), alpha=0.8) +
  scale_fill_manual(values = c("#E6AB02","#660033")) +
  xlab("Average ploidy per patient") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  geom_text_repel(aes(label=patient), position=position_jitter(0.2, seed=1),min.segment.length = 0, size=8, box.padding = 1) +
  #stat_compare_means(aes(group = WGD), label="p.format", label.y = 3,  size=8, method = "wilcox.test") +
  annotate('text', y=c(1.1), x=c(1.2), 
           label=c(paste("p", "==", pvWGD, sep = ""), paste("p", "==", pvWGD, sep = "")), size=9, parse=TRUE, hjust = 0) +
  annotate('segment', x = c(1.1), xend = c(1.8), y = c(0.6), yend = c(0.6)) +
  theme_custom() +
  theme(legend.position = "none")

p2 <- ggplot2::ggplot(data=scaaPP, aes(y=log(propSCAA), x=as.factor(type), fill=as.factor(type))) +
  geom_violin(trim = F, alpha=0.4) +
  geom_boxplot(width=0.1, outlier.size = 0) +
  geom_jitter(size=5, shape=21, position=position_jitter(0.2, seed=1), alpha=0.8) +
  scale_fill_manual(values = c("#A6CEE3","#666633")) +
  xlab("MSI status") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  geom_text_repel(aes(label=patient), position=position_jitter(0.2, seed=1),min.segment.length = 0, size=8, box.padding = 1) +
  annotate('text', y=c(1.1), x=c(1.2), 
           label=c(paste("p", "==", pvMSI, sep = ""), paste("p", "==", pvMSI, sep = "")), size=9, parse=TRUE, hjust = 0) +
  annotate('segment', x = c(1.1), xend = c(1.8), y = c(0.6), yend = c(0.6)) +
  theme_custom() +
  theme(legend.position = "none")

jpeg('tempfig.jpeg', height = (3*37.795*4.30),  width = (3*37.795*8.38))
cowplot::plot_grid(p1, p2, ncol = 2)
dev.off()

t.test(log(scaaPP$propSCAA[which(scaaPP$WGD==">2.5" & 
                                        scaaPP$patient %!in% c("C516","C542","C543"))]), 
            log(scaaPP$propSCAA[which(scaaPP$WGD=="<2.5" &
                                        scaaPP$patient %!in% c("C516","C542","C543"))]))

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
p1
dev.off()

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
p2
dev.off()

# PGA VS N.SEGMENTS =====
library(ggrepel)
#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (28), height = (20), units = "cm", res = 300)
ggplot(data=scaaPP, aes(x=pga, y=n.segments, fill=type, color=type)) +
  #geom_point(size=3, shape=21, color="black") +
  geom_point(size=5, color="black", aes(shape=WGD)) +
  stat_smooth(method = "glm") +
  scale_fill_manual(values = c("#1F78B4", "#33A02C")) +
  scale_color_manual(values = c("#1F78B4", "#33A02C")) +
  scale_shape_manual(values = c(21,24)) +
  ylab("Number of segments") +
  xlab("PGA") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  geom_text_repel(aes(label=patient), color="black", size=7, box.padding = 0.7) +
  theme_custom() +
  theme(legend.text = element_text(size=28), legend.title = element_text(size=28), legend.key.size =unit(1, "cm") )
dev.off()

# SCAA PROPORTION AND PGA/SEGMENTS CORRELATE =====

# all patients
p1 <- ggplot(data=scaaPP, aes(x=log(pga), y=log(as.numeric(propSCAA)))) +
  geom_point(size=5, shape=21, aes(fill=as.factor(remove))) +
  geom_text_repel(aes(label=patient), size=7, box.padding = 0.7) +
  scale_fill_manual(values = c(alpha("#333399",0.85), alpha("#CC9966",0.8))) +
  ggtitle("All patients") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_x_continuous(name="PGA", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "none", legend.text = element_text(size=24), legend.title = element_blank()) 

p2 <- ggplot(data=scaaPP, aes(x=(n.segments), y=log(as.numeric(propSCAA)))) +
  geom_point(size=5, shape=21, aes(fill=as.factor(remove))) +
  geom_text_repel(aes(label=patient), size=7, box.padding = 0.7) +
  scale_fill_manual(values = c(alpha("#333399",0.85), alpha("#CC9966",0.8))) +
  #ggtitle("All patients") +
  ggtitle("") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_x_continuous(name="Number of segments") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "none", legend.text = element_text(size=24), legend.title = element_blank()) 

#jpeg('tempfig.jpeg', height = (3*37.795*4.30),  width = (3*37.795*8.38))
jpeg('tempfig.jpeg', width = (40), height = (20), units = "cm", res = 300)
cowplot::plot_grid(p1, p2, ncol = 2)
dev.off()


# exc 3
p1 <- ggplot(data=scaaPP[which(scaaPP$patient %!in% c("C543","C542","C516")),], aes(x=log(pga), y=log(as.numeric(propSCAA)))) +
  geom_point(size=5, shape=21, fill=alpha("#333399",0.85)) +
  ggtitle("Exc. C516/542/543") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_x_continuous(name="PGA", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 

p2 <- ggplot(data=scaaPP[which((scaaPP$patient %!in% c("C543","C542","C516"))),], aes(x=(n.segments), y=log(as.numeric(propSCAA)))) +
  geom_point(size=5, shape=21, fill=alpha("#333399",0.85)) +
  #ggtitle("Exc. C516/542/543") +
  ggtitle("") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_x_continuous(name="Number of segments") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 

#jpeg('tempfig.jpeg', height = (3*37.795*4.30),  width = (3*37.795*8.38))
jpeg('tempfig.jpeg', width = (40), height = (20), units = "cm", res = 300)
cowplot::plot_grid(p1, p2, ncol = 2)
dev.off()

# exc 3 + sep by MSI
p1 <- ggplot(data=scaaPP[which(scaaPP$patient %!in% c("C543","C542","C516")),], aes(x=log(pga), y=log(as.numeric(propSCAA)), fill=type, color=type)) +
  geom_point(size=5, shape=21, color="black") +
  ggtitle("Exc. C516/542/543") +
  scale_color_manual(values=c("#6699CC","#666633")) +
  scale_fill_manual(values=c("#6699CC","#666633")) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_x_continuous(name="PGA", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  stat_smooth(method = "glm") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 

p2 <- ggplot(data=scaaPP[which(scaaPP$patient %!in% c("C543","C542","C516")),], aes(x=n.segments, y=log(as.numeric(propSCAA)), fill=type, color=type)) +
  geom_point(size=5, shape=21, color="black") +
  #ggtitle("Exc. C516/542/543") +
  ggtitle("") +
  scale_color_manual(values=c("#6699CC","#666633")) +
  scale_fill_manual(values=c("#6699CC","#666633")) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_x_continuous(name="Number of segments") +
  stat_smooth(method = "glm") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 

#jpeg('tempfig.jpeg', height = (3*37.795*4.30),  width = (3*37.795*8.38))
jpeg('tempfig.jpeg', width = (40), height = (20), units = "cm", res = 300)
cowplot::plot_grid(p1, p2, ncol = 2)
dev.off()

# NUMBER OF SNVS WITH SCAA BURDEN =====

# c=clonal, sc=subclonal, it=indel/trunctating, snv=snv
cm_snv <- readxl::read_excel("~/Documents/SCAA/Data/CM_SNVs.xlsx")
x <- data.frame(frac=(colSums(!is.na(cm_snv[-1])))/117)
x$pga <- scaaPP[match(rownames(x), scaaPP$patient), 8]
x$propSCAA <- scaaPP[match(rownames(x), scaaPP$patient), 2]
x$patient <- rownames(x)
x$ploidy <- scaaPP[match(x$patient, scaaPP$patient), 6]
x$MSI <- scaaPP[match(x$patient, scaaPP$patient), 7]
totalSNVcount <- read.delim("~/Documents/SCAA/Data/total_number_of_mutations.txt")
x$totalSNVcount <- totalSNVcount[match(rownames(x), totalSNVcount$Patient), 2]
x <- x[which(!is.na(x$pga)),]
x <- x[which(x$patient %in% scaaPP$patient),]

jpeg('tempfig.jpeg', width = (3*37.795*4.19), height = (3*37.795*4.35))
ggplot2::ggplot(x, aes(y=log(as.numeric(propSCAA)), x=as.numeric(frac), fill=MSI, color=MSI)) +
  geom_point(size=3, shape=21, color="black") +
  xlab("Fraction of chromatin \nmodifiers with SNV/indel") +
  scale_fill_manual(values=c("#6699CC","#666633")) +
  scale_color_manual(values=c("#6699CC","#666633")) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  stat_smooth(method = "glm") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 
dev.off()

# the relationship just with SN.indel burden
#jpeg('tempfig.jpeg', width = (3*37.795*4.19), height = (3*37.795*4.35))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot2::ggplot(x[which(x$patient %!in% c("C516","C542","C543")),], aes(y=log(as.numeric(propSCAA)), x=as.numeric(totalSNVcount), fill=MSI, color=MSI)) +
  geom_point(size=5, shape=21, color="black") +
  xlab("Total SNV/indel cout") +
  ggtitle("Exc. C516/542/543") +
  #ggtitle("All patients") +
  scale_fill_manual(values=c("#6699CC","#666633")) +
  scale_color_manual(values=c("#6699CC","#666633")) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  stat_smooth(method = "glm") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 
dev.off()

# CNAS IN CHROMATIN MODIFIERS (WITH PLOIDY RECENTRE) =====

pgaplot <- ggplot2::ggplot(scaaPP, aes(x=reorder(patient, pga_recentre), group=1)) +
  geom_area(aes(y=pga_recentre), fill="#69b3a2", alpha=0.4) +
  geom_line(aes(y=pga_recentre), color="#69b3a2", size=2) +
  geom_point(aes(y=pga_recentre), size=3, color="#69b3a2") +
  ylab("PGA \n(recentered)") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  theme_custom() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(vjust = -10),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "top", 
        strip.text = element_text(size=28), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 3, l = 0))

scaaplot <- ggplot2::ggplot(scaaPP, aes(x=reorder(patient, pga_recentre), group=1)) +
  geom_bar(aes(y=(propSCAA)), fill=alpha("#FF9933",0.6), color="#FF9933", stat="identity") +
  scale_y_continuous(name="SCAA \nburden") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  theme_custom() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(vjust = -10),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "top", 
        strip.text = element_blank(), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 3, l = 0))

# create chromatin gene list
chromatinGenes <- read.table('~/Documents/CNA/Github/Data/GeneLists/chromatinGenes', header = T, sep = "\t")
DNAmodGenes <- read.table('~/Documents/CNA/Github/Data/GeneLists/DNAmodGenes', header = T, sep = "\t")
histoneGenes <- read.table('~/Documents/CNA/Github/Data/GeneLists/histoneGenes', header = T, sep = "\t")

library(tidyverse)
hallmarkGenes <- readRDS('~/Documents/CNA/Github/Data/GeneLists/hallmarkGenes.rds')
cosmic <- read_csv('~/Documents/CNA/Github/Data/GeneLists/COSMIC 11_12_03 2022.csv')
#cm_snv <- readxl::read_excel("~/Documents/SCAA/Data/CM_SNVs.xlsx")
GOlist <- read.delim("~/Documents/CNA/Data/hg38.1_all_gene_GO_annotations.txt")
GOlist <- GOlist[GOlist$Chromosome.scaffold.name %in% c(seq(1,22,1),"X","Y"),]

chromatin.ass.genes <- unique(c(chromatinGenes$Symbol, DNAmodGenes$Symbol, histoneGenes$Symbol))
chromatin.ass.genes <- GOlist[GOlist$HGNC.symbol %in% chromatin.ass.genes,]
chromatin.ass.genes <- chromatin.ass.genes[,c(1,3,4,5,8)]
chromatin.ass.genes <- chromatin.ass.genes[!duplicated(chromatin.ass.genes),]

# per patient see if have clonal/subclonal gain/loss in each chromatin gene
output <- list()
l <- 1
for (i in 1:length(segments_ploidyrecenterd)) { # for this patient...
  wd <- segments_ploidyrecenterd[[i]]
  
  for (j in 1:nrow(chromatin.ass.genes)) { #...do they have this gene?
    print(paste(i,":", j,"/",nrow(chromatin.ass.genes)))
    gene <- chromatin.ass.genes$HGNC.symbol[j]
    geneChr <- chromatin.ass.genes$Chromosome.scaffold.name[j]
    geneStart <- chromatin.ass.genes$Gene.start..bp.[j]
    geneStop <- chromatin.ass.genes$Gene.end..bp.[j]
    wd2 <- wd[wd$chr==geneChr,]
    
    if (nrow(wd2)==0) { #dont have the sex chr
      output[[l]] <- "noChr"
      l <- l + 1
      next
    }
    
    overlapIndex <- apply(wd2, 1, function(x) DescTools::Overlap(as.numeric(c(x[5], x[6])), as.numeric(c(geneStart,geneStop)))) > 0.9*(geneStop-geneStart)
    wd2 <- wd2[overlapIndex,]
    wd2 <- wd2[!is.na(wd2$cna),]
    
    if ( sum(wd2$cna=="diploid") == nrow(wd2) ) {
      output[[l]] <- "diploid"
      l <- l + 1
    } else if ( sum(wd2$cna=="gain") == nrow(wd2) ) {
      output[[l]] <- "clonal.gain"
      l <- l + 1
    } else if ( sum(wd2$cna=="loss") == nrow(wd2) ) {
      output[[l]] <- "clonal.loss"
      l <- l + 1
    } else if ( (sum(wd2$cna=="loss") > 1) & (sum(wd2$cna=="gain") > 1) ) {
      output[[l]] <- "subclonal.gain.loss"
      l <- l + 1
    } else if ( (sum(wd2$cna=="loss") > 1)) {
      output[[l]] <- "subclonal.loss"
      l <- l + 1
    } else if ( (sum(wd2$cna=="gain") > 1)) {
      output[[l]] <- "subclonal.gain" 
      l <- l + 1
    } else {
      output[[l]] <- NA 
      l <- l + 1
    }
  }
}

chromatin.cnas.pp <- data.frame(matrix(unlist(output), nrow(chromatin.ass.genes), byrow = F)) # fill col first
rownames(chromatin.cnas.pp) <- chromatin.ass.genes$HGNC.symbol
colnames(chromatin.cnas.pp) <- unique(substr(names(dwgs.perpatient), 1, 4))
chromatin.cnas.pp <- chromatin.cnas.pp[,which(colnames(chromatin.cnas.pp) %in% scaaPP$patient)]
library(ggplot2)

# plot gains
chromatin.cnas.pp_gain <- chromatin.cnas.pp
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="clonal.loss"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="noChr"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="diploid"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="noChr"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="subclonal.loss"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="subclonal.gain.loss"] <- "subclonal.gain"
x <- chromatin.cnas.pp_gain[,colSums(!is.na(chromatin.cnas.pp_gain))>0]

keep <- chromatin.cnas.pp_gain[rowSums(!is.na(chromatin.cnas.pp_gain), na.rm = T) > 12,]
keep$geneID <- rownames(keep)

x <- tidyr::gather(keep, key="key", value="value", -geneID)
x$type <- ifelse(x$key %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
x$pga <- scaaPP[match(x$key, scaaPP$patient),8]
x$pga_recentre <- scaaPP[match(x$key, scaaPP$patient),12]
x$propSCAA <- scaaPP[match(x$key, scaaPP$patient),2]

gainplot <- ggplot(x, aes(y=geneID, x=reorder(key, pga_recentre), fill= value)) + 
  geom_tile(color="black") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  scale_fill_manual(values = c("#CC0033",alpha("#CC0033",0.5)), na.value = "#FFFFFF") +
  theme_custom() +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size=24),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "none", 
        strip.text = element_blank(), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 3, l = 0))

# plot losses
chromatin.cnas.pp_loss <- chromatin.cnas.pp
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="clonal.gain"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="noChr"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="diploid"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="noChr"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="subclonal.gain"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="subclonal.gain.loss"] <- "subclonal.loss"
x <- chromatin.cnas.pp_loss[,colSums(!is.na(chromatin.cnas.pp_loss))>0]

keep <- chromatin.cnas.pp_loss[rowSums(!is.na(chromatin.cnas.pp_loss), na.rm = T) > 12,]
keep$geneID <- rownames(keep)

x <- tidyr::gather(keep, key="key", value="value", -geneID)
x$type <- ifelse(x$key %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
x$pga <- scaaPP[match(x$key, scaaPP$patient),8]
x$pga_recentre <- scaaPP[match(x$key, scaaPP$patient),12]
x$propSCAA <- scaaPP[match(x$key, scaaPP$patient),2]

lossplot <- ggplot(x, aes(y=geneID, x=reorder(key, pga_recentre), fill= value)) + 
  geom_tile(color="black") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  scale_fill_manual(values = c("#3300FF",alpha("#3300FF",0.5)), na.value = "#FFFFFF") +
  theme_custom() +
  theme(axis.text.x = element_text(size=24, angle=90, vjust = 0.5), axis.title = element_blank(), axis.text.y = element_text(size=24),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "none",
        strip.text = element_blank(), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 3, l = 0))



#jpeg('tempfig.jpeg', width = (3*37.795*7), height = (3*37.795*16))
jpeg('tempfig.jpeg', width = (40), height = (60), units = "cm", res = 300)
cowplot::plot_grid(pgaplot, scaaplot, gainplot, lossplot, ncol = 1, align = 'v', rel_heights = c(0.12, 0.1, 0.7, 0.15))
dev.off()







x <- chromatin.cnas.pp
x[x=="noChr"] <- "diploid"
x[is.na(x)] <- "diploid"
x <- x[rowSums(x!="diploid")>0,]
x <- x[,colnames(x) %in% scaaPP$patient[which(scaaPP$type=="MSS")]]

sum(colSums(x=="clonal.gain" | x=="subclona.gain" |x=="subclonal.gain.loss")>0)
sum(colSums(x=="clonal.gain")>0)
colSums(x=="clonal.gain")
mean(colSums(x=="clonal.gain"))

sum(colSums(x=="clonal.loss" | x=="subclona.loss" |x=="subclonal.gain.loss")>0)
sum(colSums(x=="clonal.loss")>0)
colSums(x=="clonal.loss")
mean(colSums(x=="clonal.loss"))

sum(colSums(x!="diploid")>0)

x <- t(x)
x <- x[,colSums(x!="diploid")>18]

y <- data.frame(rowSums(x=="clonal.gain"))
y <- data.frame(rowSums(x=="clonal.loss"))

# NUMBER OF CNA IN CM WITH PGA AND SCAA BURDEN =====

# just SWI/SNF
swisnf <-  c("ARID1A","ARID1B","ARID2","PBRM1","SMARCA4","SMARCB1", "SMARCA2", "SMARCC1", "SMARCA5", "SMARCC2", "SMARCD1", "SMARCD2", "SMARCD3")
x <- chromatin.cnas.pp[which(rownames(chromatin.cnas.pp) %in% swisnf),]
x <- chromatin.cnas.pp[grep("KDM", rownames(chromatin.cnas.pp)),]
x <- data.frame(clonal=colSums(x=="clonal.loss" | x=="clonal.gain", na.rm = T ),
                subclonal=colSums(x=="subclonal.loss" | x=="subclonal.gain" | x=="subclonal.gain.loss", na.rm = T ))
x$any <- x$clonal+x$subclonal
x$pga <- scaaPP[match(rownames(x), scaaPP$patient), 8]
x$pga_recentre <- scaaPP[match(rownames(x), scaaPP$patient), 12]
x$propSCAA <- scaaPP[match(rownames(x), scaaPP$patient), 2]
x$frac.clonal <- x$clonal/259
x$frac.subclonal <- x$subclonal/259
x$frac.any <- x$any/259
x$WGD <- scaaPP[match(rownames(x), scaaPP$patient),6]
x$type <- scaaPP[match(rownames(x), scaaPP$patient),7]
x$remove <- scaaPP[match(rownames(x), scaaPP$patient),11]
x$patient <- rownames(x)
x <- x[which(x$patient %in% scaaPP$patient),]

# all CM
x <- data.frame(clonal=colSums(chromatin.cnas.pp=="clonal.loss" | chromatin.cnas.pp=="clonal.gain", na.rm = T ),
                subclonal=colSums(chromatin.cnas.pp=="subclonal.loss" | chromatin.cnas.pp=="subclonal.gain" | chromatin.cnas.pp=="subclonal.gain.loss", na.rm = T ))
x$any <- x$clonal+x$subclonal
x$pga <- scaaPP[match(rownames(x), scaaPP$patient), 8]
x$pga_recentre <- scaaPP[match(rownames(x), scaaPP$patient), 12]
x$propSCAA <- scaaPP[match(rownames(x), scaaPP$patient), 2]
x$frac.clonal <- x$clonal/259
x$frac.subclonal <- x$subclonal/259
x$frac.any <- x$any/259
x$WGD <- scaaPP[match(rownames(x), scaaPP$patient),6]
x$type <- scaaPP[match(rownames(x), scaaPP$patient),7]
x$remove <- scaaPP[match(rownames(x), scaaPP$patient),11]
x$patient <- rownames(x)
x <- x[which(x$patient %in% scaaPP$patient),]



#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot2::ggplot(x[which(x$patient %!in% c("C516","C542","C543")),], aes(x=as.numeric(pga), y=as.numeric(frac.any))) +
  geom_point(size=5, shape=21, aes(fill=as.factor(remove))) +
  geom_text_repel(aes(label=patient), size=8, box.padding = 1) +
  scale_fill_manual(values = c(alpha("#333399",0.85), alpha("#CC9966",0.8))) +
  scale_y_continuous(name="Fraction of chromatin \nmodifiers with CNA", breaks = seq(0,1,0.25)) +
  scale_x_continuous(name="PGA", breaks = seq(0,1,0.25)) +
  ggtitle("") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "none",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
dev.off()

#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot2::ggplot(x[which(x$patient %!in% c("C516","C542","C543")),], aes(x=log(as.numeric(propSCAA)), y=as.numeric(frac.any))) +
  geom_point(size=5, shape=21, aes(fill=as.factor(remove))) +
  geom_text_repel(aes(label=patient), size=8, box.padding = 1) +
  scale_fill_manual(values = c(alpha("#333399",0.85), alpha("#CC9966",0.8))) +
  scale_y_continuous(name="Fraction of chromatin \nmodifiers with CNA", breaks = seq(0,1,0.25)) +
  scale_x_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  ggtitle("Exc. C516/542/543") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  theme_custom() +
  theme(legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "none",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
dev.off()


# ANY SCAA INDIVIDUALLY CORRELATES WITH PGA =====

# pull all peaks associated with promo/enh which are SCAAs across patient exc 3
peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")]
scaa <- scaa[,substr(colnames(scaa),1,4) %!in% c("C516","C542","C543")]

scaa.peaks <- rownames(scaa[rowSums(scaa)>2,])
scaa.peaks <- scaa.peaks[scaa.peaks %in% c(promoter,enhancer)]
scaa.long$pga <- scaaPP[match(scaa.long$patient, scaaPP$patient), 8]
scaa.long$pga_recentre <- scaaPP[match(scaa.long$patient, scaaPP$patient), 12]

# split 
set.seed(1234)
training.samples <- sample(substr(colnames(scaa),1,4), 21)
testing.samples  <- substr(colnames(scaa),1,4)[substr(colnames(scaa),1,4) %!in% training.samples]
#training.samples <- substr(colnames(scaa),1,4)

# for each scaa.peak, pull pag of patients with and without scaa, t.test to compare
library(lsr)
p <- list()
direction <- list()
d <- list()
nSCAA <- list()
peak <- list()
for (i in 1:length(scaa.peaks)) {
  print(i)
  wd <- scaa.long[which(scaa.long$peak==scaa.peaks[i] ),]
  wd$pga <- scaaPP[match(wd$patient, scaaPP$patient), 8]
  #wd$pga <- scaaPP[match(wd$patient, scaaPP$patient), 9]
  if(sum(wd$value)<=1) {
    next
  }
  peak[[i]] <- scaa.peaks[i]
  nSCAA[[i]] <- sum(wd$value)
  p[[i]] <- t.test((wd$pga[which(wd$value==1)]), (wd$pga[which(wd$value==0)]))$p.value
  if( t.test((wd$pga[which(wd$value==1)]), (wd$pga[which(wd$value==0)]))$estimate[1] >  t.test((wd$pga[which(wd$value==1)]), (wd$pga[which(wd$value==0)]))$estimate[2] ) {
    direction[[i]] <- "positive"
  } else {
    direction[[i]] <- "negative"
  }
  d[[i]] <- cohensD(pga ~ value, data = wd)
}

# for each peak, add FDR and gene associated
x <- data.frame(peak=unlist(peak), nSCAA=unlist(nSCAA), direction=unlist(direction), d=unlist(d), pvalue=unlist(p))
x$FDR <- p.adjust(x$pvalue, "fdr")
x$symbol <- peakAnn_all[match(x$peak, peakAnn_all$peak), 17]
x$entrez <- peakAnn_all[match(x$peak, peakAnn_all$peak), 14]

x <- x[which(x$FDR<0.05),] 
openClose <- readxl::read_excel("~/Documents/SCAA/Data/41586_2022_5202_MOESM9_ESM.xlsx")
x$SCAA <- openClose[match(x$peak, openClose$peak.peak), 7]
x$cosmic <- ifelse(x$symbol %in% cosmic$`Gene Symbol`, TRUE, FALSE)

#overlap with rep genes 
betaregRepGenes <- readRDS("~/Documents/CNA/Github/Data/betaregRepGenes.rds")
x$betaregRepGenes <- ifelse(x$symbol %in% betaregRepGenes$symbol, TRUE, FALSE)

# boxplot for a gene
y <- scaa.long[which(scaa.long$symbol=="FOXK2"),]
ggplot2::ggplot(data=y[which((y$patient %!in% c("C516","C542","C543"))),], aes(y=(pga), x=as.factor(value), fill=as.factor(value))) +
  geom_violin(trim = T, alpha=0.4) +
  geom_jitter(size=3, shape=21, position=position_jitter(0.2, 0.01, seed=1), alpha=0.8) +
  geom_boxplot(width=0.1, outlier.size = 0) +
  #scale_fill_manual(values = c("#E6AB02","#660033")) +
  scale_x_discrete(name=NULL, labels = c("0" = "Not SCAA",
                                         "1" = "SCAA")) +
  stat_compare_means(method = "t.test", label.y = 1, label.x.npc = 0.2, size=8) +
  ggtitle(unique(y$symbol)) +
  scale_y_continuous(name="PGA") +
  geom_text_repel(aes(label=patient), position=position_jitter(0.2, 0.01, seed=1),min.segment.length = 0, size=5) +
  theme_custom() +
  theme(legend.position = "none")

# GO on the FDR<0.05 peaks shows nothing
GOx <- limma::goana(list(positive=x$entrez[which(x$direction=="positive")], negative=x$entrez[which(x$direction=="negative")]))
GOx$p.adj.pos = p.adjust(GOx$P.positive, method = "fdr")
GOx$p.adj.neg = p.adjust(GOx$P.negative, method = "fdr")

data <- GOx[which(GOx$p.adj.neg < 0.05),]
data$ID <- rownames(data)
simMatrix <- rrvgo::calculateSimMatrix(data$ID, orgdb="org.Hs.eg.db", ont = c("BP"), method="Rel")
scores <- setNames(-log10(data$p.adj.neg), data$ID)
reducedTerms <- rrvgo::reduceSimMatrix(simMatrix, scores = scores, threshold=0.9, orgdb="org.Hs.eg.db")
parentGO <- data.frame(value=tapply(reducedTerms$score, reducedTerms$parentTerm, FUN=sum))
parentGO$term <- rownames(parentGO)
parentGO$fraction <- ( (parentGO$value/sum(parentGO$value))*100 )


# above doesnt work with recnetered data, maybe just detecting scaa association with WGD
# do chi for scaa and wgd
library(lsr)
p <- list()
direction <- list()
SCAA.WGD <- list()
SCAA.notWGD <- list()
nSCAA <- list()
peak <- list()
for (i in 1:length(scaa.peaks)) {
  print(i)
  wd <- scaa.long[which(scaa.long$peak==scaa.peaks[i] ),]
  wd$pga_recentre <- scaaPP[match(wd$patient, scaaPP$patient), 12]
  wd$WGD <- scaaPP[match(wd$patient, scaaPP$patient), 6]
  if(sum(wd$value)<=1) {
    next
  }
  peak[[i]] <- scaa.peaks[i]
  nSCAA[[i]] <- sum(wd$value)
  x <- data.frame(NotSCAA=c(nrow(wd[which(wd$value==0 & wd$WGD=="<2.5"),]), nrow(wd[which(wd$value==0 & wd$WGD==">2.5"),])), 
                  SCAA=c(nrow(wd[which(wd$value==1 & wd$WGD=="<2.5"),]), nrow(wd[which(wd$value==1 & wd$WGD==">2.5"),])))
  rownames(x) <- c("<2.5",">2.5")
  chisq <- chisq.test(x)
  p[[i]] <- chisq$p.value
  SCAA.WGD[[i]] <- (chisq$observed/chisq$expected)[2,2]
  SCAA.notWGD[[i]] <- (chisq$observed/chisq$expected)[1,2]
}

x <- data.frame(peak=unlist(peak), nSCAA=unlist(nSCAA), enrich.SCAA.WGD=unlist(SCAA.WGD), enrich.SCAA.notWGD=unlist(SCAA.notWGD), pvalue=unlist(p))
x$FDR <- p.adjust(x$pvalue, "fdr")
x$symbol <- peakAnn_all[match(x$peak, peakAnn_all$peak), 17]
x$entrez <- peakAnn_all[match(x$peak, peakAnn_all$peak), 14]
x <- x[which(x$FDR<0.05),] 
openClose <- readxl::read_excel("~/Documents/SCAA/Data/41586_2022_5202_MOESM9_ESM.xlsx")
x$SCAA <- openClose[match(x$peak, openClose$peak.peak), 7]















# model
#data <- scaa.long[which(scaa.long$peak %in% x$peak),c(10,11,22,28)]
data <- scaa.long[which(scaa.long$peak %in% x$peak),c(10,11,22,28)]
data$value <- as.factor(data$value) 
data_wide <- spread(data, peak, value)

training.data <- data_wide[which(data_wide$patient %in% training.samples),]
testing.data <- data_wide[which(data_wide$patient %in% testing.samples),]

model <- lm(pga ~., data = data_wide[-1])
stepMod <- stepAIC(model, direction="both")
summary(stepMod)
predictions <- stepMod %>% predict(testing.data)
data.frame( R2 = R2(predictions, testing.data$pga),
            RMSE = RMSE(predictions, testing.data$pga),
            MAE = MAE(predictions, testing.data$pga))

train.control <- trainControl(method = "LOOCV")
model <- train(pga ~., data = data_wide[-1], method = "lmStepAIC",
               trControl = train.control)
model$results
model$finalModel
summary(model$finalModel)

# plot predicted 
testing.data$predicted <- model %>% predict(testing.data)
data_wide$predicted <- model %>% predict(data_wide)
ggplot(testing.data, aes(x=pga, y=predicted)) +
  geom_point(shape=21) +
  geom_text_repel(aes(label=patient)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, )


par(mfrow = c(2, 2))
plot(model)






# ANY CM SNV CORRELATE WITH SCAA BURDEN =====

# an cm snv corr with scaa burden
cm_snv <- readxl::read_excel("~/Documents/SCAA/Data/CM_SNVs.xlsx")
cm_snv[is.na(cm_snv)] <- "none"
cm_snv[cm_snv=="c.it"] <- "none"
cm_snv[cm_snv=="sc.it"] <- "none"
gene <- list()
nMut <- list()
p <- list()
direction <- list()
for (i in 1:nrow(cm_snv)) {
  wd <- data.frame(group=t(cm_snv[i,-1]))
  wd$propSCAA <- scaaPP[match(rownames(wd), scaaPP$patient), 2]
  wd <- na.omit(wd)
  if (sum(wd$group!="none")<2) {
    next
  }
  gene[[i]] <- cm_snv$GeneID[i]
  nMut[[i]] <- sum(wd$group!="none")
  p[[i]] <- t.test((wd$propSCAA[which(wd$group!="none")]), (wd$propSCAA[which(wd$group=="none")]))$p.value
  if( t.test((wd$propSCAA[which(wd$group!="none")]), (wd$propSCAA[which(wd$group=="none")]))$estimate[1] >  t.test((wd$propSCAA[which(wd$group!="none")]), (wd$propSCAA[which(wd$group=="none")]))$estimate[2] ) {
    direction[[i]] <- "positive"
  } else {
    direction[[i]] <- "negative"
  }
}

x <- data.frame(peak=unlist(gene), nSCAA=unlist(nMut), direction=unlist(direction), pvalue=unlist(p))
x$FDR <- p.adjust(x$pvalue, "fdr")

# specific to clonal trunc
cm_snv <- readxl::read_excel("~/Documents/SCAA/Data/CM_SNVs.xlsx")
cm_snv[is.na(cm_snv)] <- "none"
#swisnf <-  c("ARID1A","ARID1B","ARID2","PBRM1","SMARCA4","SMARCB1", "SMARCA2", "SMARCC1", "SMARCA5", "SMARCC2", "SMARCD1", "SMARCD2", "SMARCD3")
#swisnf <- GOlist[GOlist$HGNC.symbol %in% swisnf,]
#swisnf <- swisnf[,c(1,3,4,5,8)]
#swisnf <- swisnf[!duplicated(swisnf),]
#cm_snv <- cm_snv[which(cm_snv$GeneID %in% swisnf$HGNC.symbol),]

x<- data.frame(any=colSums(cm_snv[-1]!="none"), 
               it=colSums(cm_snv[-1]=="c.it" | cm_snv[-1]=="sc.it"),
               clonal.it=colSums(cm_snv[-1]=="c.it"),
               clonal.it.swisnf=colSums(cm_snv[which(cm_snv$GeneID %in% swisnf$HGNC.symbol),-1]=="c.it"),
               snv=colSums(cm_snv[-1]=="c.snv" | cm_snv[-1]=="sc.snv"),
               binary.any=(colSums(cm_snv[-1]!="none")>0),
               binary.it=(colSums(cm_snv[-1]=="c.it" | cm_snv[-1]=="sc.it"))>0,
               binary.snv=(colSums(cm_snv[-1]=="c.snv" | cm_snv[-1]=="sc.snv"))>0,
               binary.clonal.it=(colSums(cm_snv[-1]=="c.it"))>0,
               binary.clonal.it.swinf=(colSums(cm_snv[which(cm_snv$GeneID %in% swisnf$HGNC.symbol),-1]=="c.it"))>0)

x$propSCAA <- scaaPP[match(rownames(x), scaaPP$patient), 2]
x$pga <- scaaPP[match(rownames(x), scaaPP$patient), 8]
x$n.segments <- scaaPP[match(rownames(x), scaaPP$patient), 9]
x$MSI <- scaaPP[match(rownames(x), scaaPP$patient), 7]
x <- na.omit(x)
t.test((x$propSCAA[which(x$binary.any==TRUE)]), (x$propSCAA[which(x$binary.any==FALSE)]))
t.test((x$propSCAA[which(x$binary.it==TRUE)]), (x$propSCAA[which(x$binary.it==FALSE)]))
t.test((x$propSCAA[which(x$binary.snv==TRUE)]), (x$propSCAA[which(x$binary.snv==FALSE)]))
t.test((x$propSCAA[which(x$binary.clonal.it==TRUE)]), (x$propSCAA[which(x$binary.clonal.it==FALSE)]))
t.test((x$propSCAA[which(x$binary.clonal.it.swinf==TRUE)]), (x$propSCAA[which(x$binary.clonal.it.swinf==FALSE)]))
summary(lm(propSCAA ~ clonal.it, data=x))
summary(lm(propSCAA ~ clonal.it.swisnf, data=x))

t.test((x$pga[which(x$binary.any==TRUE)]), (x$pga[which(x$binary.any==FALSE)]))
t.test((x$pga[which(x$binary.it==TRUE)]), (x$pga[which(x$binary.it==FALSE)])) #
t.test((x$pga[which(x$binary.snv==TRUE)]), (x$pga[which(x$binary.snv==FALSE)]))
t.test((x$pga[which(x$binary.clonal.it==TRUE)]), (x$pga[which(x$binary.clonal.it==FALSE)])) #
t.test((x$pga[which(x$binary.clonal.it.swinf==TRUE)]), (x$pga[which(x$binary.clonal.it.swinf==FALSE)])) ##
t.test((x$n.segments[which(x$binary.clonal.it==TRUE)]), (x$n.segments[which(x$binary.clonal.it==FALSE)])) #
t.test((x$n.segments[which(x$binary.clonal.it.swinf==TRUE)]), (x$n.segments[which(x$binary.clonal.it.swinf==FALSE)])) #

t.test((x$propSCAA[which(x$MSI=="MSS" & x$binary.clonal.it==TRUE)]), (x$propSCAA[which(x$MSI=="MSS" & x$binary.clonal.it==FALSE)]))
t.test((x$pga[which(x$MSI=="MSS" & x$binary.clonal.it==TRUE)]), (x$pga[which(x$MSI=="MSS" & x$binary.clonal.it==FALSE)]))
t.test((x$propSCAA[which(x$MSI=="MSS" & x$binary.clonal.it.swinf==TRUE)]), (x$propSCAA[which(x$MSI=="MSS" & x$binary.clonal.it.swinf==FALSE)]))
t.test((x$pga[which(x$MSI=="MSS" & x$binary.clonal.it.swinf==TRUE)]), (x$pga[which(x$MSI=="MSS" & x$binary.clonal.it.swinf==FALSE)]))
t.test((x$n.segments[which(x$MSI=="MSS" & x$binary.clonal.it==TRUE)]), (x$n.segments[which(x$MSI=="MSS" & x$binary.clonal.it==FALSE)]))


ggplot2::ggplot(data=x, aes(y=(pga), x=as.factor(binary.clonal.it), fill=as.factor(binary.clonal.it))) +
  geom_violin(trim = T, alpha=0.4) +
  geom_jitter(size=3, shape=21, position=position_jitter(0.2, 0.01, seed=1), alpha=0.8) +
  geom_boxplot(width=0.1, outlier.size = 0) +
  scale_fill_manual(values = c("#E6AB02","#660033")) +
  scale_x_discrete(name=NULL, labels = c("0" = "Not SCAA",
                                         "1" = "SCAA")) +
  stat_compare_means(method = "t.test", label.y = 1, label.x.npc = 0.2, size=8) +
  scale_y_continuous(name="PGA") +
  geom_text_repel(aes(label=rownames(x)), position=position_jitter(0.2, 0.01, seed=1),min.segment.length = 0, size=5) +
  theme_custom() +
  theme(legend.position = "none")

# ======
# ANY SPECIFIC CNA IN CM WITH SCAA BURDEN =====

# pull all peaks associated with promo/enh which are SCAAs across patient exc 3
peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")]
scaa <- scaa[,substr(colnames(scaa),1,4) %!in% c("C516","C542","C543")]

scaa.peaks <- rownames(scaa[rowSums(scaa)>2,])
scaa.peaks <- scaa.peaks[scaa.peaks %in% c(promoter,enhancer)]
scaa.long$propSCAA <- scaaPP[match(scaa.long$patient, scaaPP$patient), 2]

# split 
set.seed(1234)
training.samples <- sample(substr(colnames(scaa),1,4), 21)
testing.samples  <- substr(colnames(scaa),1,4)[substr(colnames(scaa),1,4) %!in% training.samples]
#training.samples <- substr(colnames(scaa),1,4)

# for each scaa.peak, pull pag of patients with and without scaa, t.test to compare
library(lsr)
p <- list()
direction <- list()
d <- list()
nSCAA <- list()
peak <- list()
for (i in 1:length(scaa.peaks)) {
  print(i)
  wd <- scaa.long[which(scaa.long$peak==scaa.peaks[i] ),]
  wd$propSCAA <- scaaPP[match(wd$patient, scaaPP$patient), 2]
  if(sum(wd$value)<=1) {
    next
  }
  peak[[i]] <- scaa.peaks[i]
  nSCAA[[i]] <- sum(wd$value)
  p[[i]] <- t.test((wd$propSCAA[which(wd$value==1)]), (wd$propSCAA[which(wd$value==0)]))$p.value
  if( t.test((wd$propSCAA[which(wd$value==1)]), (wd$propSCAA[which(wd$value==0)]))$estimate[1] >  t.test((wd$propSCAA[which(wd$value==1)]), (wd$propSCAA[which(wd$value==0)]))$estimate[2] ) {
    direction[[i]] <- "positive"
  } else {
    direction[[i]] <- "negative"
  }
  d[[i]] <- cohensD(propSCAA ~ value, data = wd)
}

# for each peak, add FDR and gene associated
x <- data.frame(peak=unlist(peak), nSCAA=unlist(nSCAA), direction=unlist(direction), d=unlist(d), pvalue=unlist(p))
x$FDR <- p.adjust(x$pvalue, "fdr")
x$symbol <- peakAnn_all[match(x$peak, peakAnn_all$peak), 17]
x$entrez <- peakAnn_all[match(x$peak, peakAnn_all$peak), 14]

x <- x[which(x$FDR<0.05),] 
openClose <- readxl::read_excel("~/Documents/SCAA/Data/41586_2022_5202_MOESM9_ESM.xlsx")
x$SCAA <- openClose[match(x$peak, openClose$peak.peak), 7]

# boxplot for a gene
y <- scaa.long[which(scaa.long$symbol=="FOXK2"),]
ggplot2::ggplot(data=y[which((y$patient %!in% c("C516","C542","C543"))),], aes(y=(pga), x=as.factor(value), fill=as.factor(value))) +
  geom_violin(trim = T, alpha=0.4) +
  geom_jitter(size=3, shape=21, position=position_jitter(0.2, 0.01, seed=1), alpha=0.8) +
  geom_boxplot(width=0.1, outlier.size = 0) +
  #scale_fill_manual(values = c("#E6AB02","#660033")) +
  scale_x_discrete(name=NULL, labels = c("0" = "Not SCAA",
                                         "1" = "SCAA")) +
  stat_compare_means(method = "t.test", label.y = 1, label.x.npc = 0.2, size=8) +
  ggtitle(unique(y$symbol)) +
  scale_y_continuous(name="PGA") +
  geom_text_repel(aes(label=patient), position=position_jitter(0.2, 0.01, seed=1),min.segment.length = 0, size=5) +
  theme_custom() +
  theme(legend.position = "none")

# GO on the FDR<0.05 peaks shows nothing
GOx <- limma::goana(list(positive=x$entrez[which(x$direction=="positive")], negative=x$entrez[which(x$direction=="negative")]))
GOx$p.adj.pos = p.adjust(GOx$P.positive, method = "fdr")
GOx$p.adj.neg = p.adjust(GOx$P.negative, method = "fdr")




# DO HIGH SCAA PATIENTS HAVE CNA IN CHROMATIN PREDICTORS? =====

betaregRepGenes <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/betaregRepGenes.rds")
chromatin.cnas.pp_repgenes <- chromatin.cnas.pp[rownames(chromatin.cnas.pp) %in% betaregRepGenes$symbol,]

# plot gains
chromatin.cnas.pp_gain <- chromatin.cnas.pp_repgenes
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="clonal.loss"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="noChr"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="diploid"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="noChr"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="subclonal.loss"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="subclonal.gain.loss"] <- "subclonal.gain"
hist(rowSums(!is.na(chromatin.cnas.pp_gain), na.rm = T))

keep <- chromatin.cnas.pp_gain[rowSums(!is.na(chromatin.cnas.pp_gain), na.rm = T) > 0,]
keep$geneID <- rownames(keep)

x <- tidyr::gather(keep, key="key", value="value", -geneID)
x$type <- ifelse(x$key %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
x$pga <- scaaPP[match(x$key, scaaPP$patient),6]

gainplot <- ggplot(x, aes(y=geneID, x=reorder(key, pga), fill= value)) + 
  geom_tile(color="black") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  scale_fill_manual(values = c("#CC0033",alpha("#CC0033",0.5)), na.value = "#FFFFFF") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5), axis.title = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "bottom", 
        strip.text = element_blank(), panel.spacing = unit(2, "lines"), panel.border = element_rect(size=1, fill=NA))

jpeg('tempfig.jpeg', width = (3*37.795*10), height = (3*37.795*15))
cowplot::plot_grid(pgaplot, scaaplot, gainplot, ncol = 1, align = 'v', rel_heights = c(0.12, 0.1, 1))
dev.off()


# plot losses
chromatin.cnas.pp_loss <- chromatin.cnas.pp_repgenes
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="clonal.gain"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="noChr"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="diploid"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="noChr"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="subclonal.gain"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="subclonal.gain.loss"] <- "subclonal.loss"
hist(rowSums(!is.na(chromatin.cnas.pp_loss), na.rm = T))

keep <- chromatin.cnas.pp_loss[rowSums(!is.na(chromatin.cnas.pp_loss), na.rm = T) > 4,]
keep$geneID <- rownames(keep)

x <- tidyr::gather(keep, key="key", value="value", -geneID)
x$type <- ifelse(x$key %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
x$pga <- scaaPP[match(x$key, scaaPP$patient),6]

lossplot <- ggplot(x, aes(y=geneID, x=reorder(key, pga), fill= value)) + 
  geom_tile(color="black") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  scale_fill_manual(values = c("#3300FF",alpha("#3300FF",0.5)), na.value = "#FFFFFF") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5), axis.title = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "bottom",
        strip.text = element_blank(), panel.spacing = unit(2, "lines"), panel.border = element_rect(size=1, fill=NA))

jpeg('tempfig.jpeg', width = (3*37.795*10), height = (3*37.795*15))
cowplot::plot_grid(pgaplot, scaaplot, lossplot, ncol = 1, align = 'v', rel_heights = c(0.12, 0.1, 1))
dev.off()



# DO THE 3 HAVE HIGH PEAK OR GENETIC ITH =====
stdev.per.region <- readRDS("~/Documents/SCAA/Data/stdev.per.region_cn.norm.peaks.rds")
patientwise.pITH <- readRDS("~/Documents/SCAA/Data/patientwise.pITH_cn.norm.peaks.rds")
regionwise.pITH <- readRDS("~/Documents/SCAA/Data/regionwise.pITH_cn.norm.peaks.rds")
cohortwise.pIPH <- readRDS("~/Documents/SCAA/Data/cohortwise.pIPH_cn.norm.peaks.rds")

hist(as.numeric(patientwise.pITH$pITH))

patientwise.pITH$propSCAA <- scaaPP[match(rownames(patientwise.pITH), scaaPP$patient), 2]
patientwise.pITH$remove <- scaaPP[match(rownames(patientwise.pITH), scaaPP$patient), 11]
patientwise.pITH$pga <- scaaPP[match(rownames(patientwise.pITH), scaaPP$patient), 8]
patientwise.pITH$patient <- rownames(patientwise.pITH)

ggplot(data=patientwise.pITH[which(patientwise.pITH$remove=="FALSE" & patientwise.pITH$patient!="C519"),], aes(x=(as.numeric(pITH)), y=log(as.numeric(propSCAA)))) +
  geom_point(size=3, shape=21, aes(fill=as.factor(remove)))+
  geom_text_repel(aes(label=patient), size=5) +
  scale_fill_manual(values = c(alpha("#333399",0.85), alpha("#CC9966",0.8))) +
  ggtitle("All patients") +
  #scale_y_continuous(name="SCAA burden", 
  #                   breaks = c(-10,-7.5,-5,-2.5,0),
  #                   limits = c(-10,0),
  #                   labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  #scale_x_continuous(name="PGA", 
  #                   breaks = c(-10,-7.5,-5,-2.5,0),
  #                   labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none", legend.text = element_text(size=24), legend.title = element_blank()) 

# DE ON HIGHER V LOW PGA =====

dds <- readRDS('~/Documents/SCAA/Data/EPICC/EPICC_expression/allgenes.dds.ensembl.rds')
res <- results(dds)
vst <- readRDS('~/Documents/SCAA/Data/EPICC/EPICC_expression/allgenes.vsd.ensembl.rds')

plotPCA(vst, intgroup=c("Patient")) +
  geom_label(aes(label = Patient))

res <- results(dds, contrast=c("Patient","C539","C555"))

par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=NULL, xlim=c(-3,3), cex=1, cex.axis=2, cex.lab=2, cex.main=2))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


res.sig <- as.data.frame(res)
res.sig <- res.sig[which(res.sig$padj<0.05),]

library("org.Hs.eg.db")
res.sig$ensemble <- gsub("\\..*","", rownames(res.sig))
res.sig$gene_id <- mapIds(org.Hs.eg.db, keys = res.sig$ensemble,
                          column = c('SYMBOL'), keytype = 'ENSEMBL')
res.sig$entrez <- mapIds(org.Hs.eg.db, keys = res.sig$ensemble,
                         column = c('ENTREZID'), keytype = 'ENSEMBL')

hallmarkGenes <- readRDS('~/Documents/CNA/Github/Data/GeneLists/hallmarkGenes.rds')
cosmic <- read_csv('~/Documents/CNA/Github/Data/GeneLists/COSMIC 11_12_03 2022.csv')
chromatinGenes <- read.table('~/Documents/CNA/Github/Data/GeneLists/chromatinGenes', header = T, sep = "\t")
DNAmodGenes <- read.table('~/Documents/CNA/Github/Data/GeneLists/DNAmodGenes', header = T, sep = "\t")
histoneGenes <- read.table('~/Documents/CNA/Github/Data/GeneLists/histoneGenes', header = T, sep = "\t")

res.sig$hallmark <- ifelse(res.sig$gene_id %in% hallmarkGenes$gene_symbol, TRUE, FALSE)
res.sig$cosmic <- ifelse(res.sig$gene_id %in% cosmic$`Gene Symbol`, TRUE, FALSE)
res.sig$chromatin <- ifelse(res.sig$gene_id %in% chromatinGenes$Symbol, TRUE, FALSE)
res.sig$DNAmodGenes <- ifelse(res.sig$gene_id %in% DNAmodGenes$Symbol, TRUE, FALSE)
res.sig$histoneGenes <- ifelse(res.sig$gene_id %in% histoneGenes$Symbol, TRUE, FALSE)


# for each chromatin gene, find those who exp corr with SCAA
cm_snv <- readxl::read_excel("~/Documents/SCAA/Data/CM_SNVs.xlsx")
chromatinGenes <- read.table('~/Documents/CNA/Github/Data/GeneLists/chromatinGenes', header = T, sep = "\t")
DNAmodGenes <- read.table('~/Documents/CNA/Github/Data/GeneLists/DNAmodGenes', header = T, sep = "\t")
histoneGenes <- read.table('~/Documents/CNA/Github/Data/GeneLists/histoneGenes', header = T, sep = "\t")

chromatinGenesAll <- data.frame(symbol=unique(c(cm_snv$GeneID, chromatinGenes$Symbol, DNAmodGenes$Symbol, histoneGenes$Symbol)))
chromatinGenesAll <- data.frame(symbol=c("ARID1A","ARID1B","ARID2","PBRM1","SMARCA4","SMARCB1", "SMARCA2", "SMARCC1", "SMARCA5", "SMARCC2", "SMARCD1", "SMARCD2", "SMARCD3"))
chromatinGenesAll$ENSEMBL <- mapIds(org.Hs.eg.db, keys = chromatinGenesAll$symbol,
                         column = c('ENSEMBL'), keytype = 'SYMBOL')

geneexp <- as.data.frame(assay(vst))

output <- list()
for (i in 1:nrow(chromatinGenesAll)) {
  ensembl <- chromatinGenesAll$ENSEMBL[i]
  if (sum(!is.na(ensembl))==0) {
    next
  }
  wd <- data.frame(t(geneexp[which(rownames(geneexp)==ensembl),]))
  wd$scaaPP <- scaaPP[match(substr(rownames(wd), 1, 4), scaaPP$patient), 2]
  colnames(wd) <- c("exp","scaaPP")
  mod <- lm(scaaPP ~ exp, data=wd)
  summary(mod)[13]
  output[[i]] <- setNames(c(chromatinGenesAll$symbol[i], summary(mod)[9], lmp(mod)),
                          c("symbol","adjR2","pvalue"))
}
x <- data.frame(do.call(rbind, output))
x$adjR2 <- as.numeric(x$adjR2)
x$pvalue <- as.numeric(x$pvalue)
x$p.adjust <- p.adjust(x$pvalue, method="fdr")

expressed_genes_by_patient <- readRDS('~/Documents/SCAA/Data/expressed_genes_by_patient.rds')

# CNAS IN SWI/SNF =====

GOlist <- read.delim("~/Documents/CNA/Data/hg38.1_all_gene_GO_annotations.txt")
GOlist <- GOlist[GOlist$Chromosome.scaffold.name %in% c(seq(1,22,1),"X","Y"),]

swisnf <-  c("ARID1A","ARID1B","ARID2","PBRM1","SMARCA4","SMARCB1", "SMARCA2", "SMARCC1", "SMARCA5", "SMARCC2", "SMARCD1", "SMARCD2", "SMARCD3")
swisnf <- GOlist[GOlist$HGNC.symbol %in% swisnf,]
swisnf <- swisnf[,c(1,3,4,5,8)]
swisnf <- swisnf[!duplicated(swisnf),]

output <- list()
l <- 1
for (i in 1:length(segments)) { # for this patient...
  wd <- segments[[i]]
  
  for (j in 1:nrow(swisnf)) { #...do they have this gene?
    print(paste(i,":", j,"/",nrow(swisnf)))
    gene <- swisnf$HGNC.symbol[j]
    geneChr <- swisnf$Chromosome.scaffold.name[j]
    geneStart <- swisnf$Gene.start..bp.[j]
    geneStop <- swisnf$Gene.end..bp.[j]
    wd2 <- wd[wd$chr==geneChr,]
    
    if (nrow(wd2)==0) { #dont have the sex chr
      output[[l]] <- "noChr"
      l <- l + 1
      next
    }
    
    overlapIndex <- apply(wd2, 1, function(x) DescTools::Overlap(as.numeric(c(x[5], x[6])), as.numeric(c(geneStart,geneStop)))) > 0.9*(geneStop-geneStart)
    wd2 <- wd2[overlapIndex,]
    wd2 <- wd2[!is.na(wd2$cna),]
    
    if ( sum(wd2$cna=="diploid") == nrow(wd2) ) {
      output[[l]] <- "diploid"
      l <- l + 1
    } else if ( sum(wd2$cna=="gain") == nrow(wd2) ) {
      output[[l]] <- "clonal.gain"
      l <- l + 1
    } else if ( sum(wd2$cna=="loss") == nrow(wd2) ) {
      output[[l]] <- "clonal.loss"
      l <- l + 1
    } else if ( (sum(wd2$cna=="loss") > 1) & (sum(wd2$cna=="gain") > 1) ) {
      output[[l]] <- "subclonal.gain.loss"
      l <- l + 1
    } else if ( (sum(wd2$cna=="loss") > 1)) {
      output[[l]] <- "subclonal.loss"
      l <- l + 1
    } else if ( (sum(wd2$cna=="gain") > 1)) {
      output[[l]] <- "subclonal.gain" 
      l <- l + 1
    } else {
      output[[l]] <- NA 
      l <- l + 1
    }
  }
}

swisnf.cnas.pp <- data.frame(matrix(unlist(output), nrow(swisnf), byrow = F)) # fill col first
rownames(swisnf.cnas.pp) <- swisnf$HGNC.symbol
colnames(swisnf.cnas.pp) <- unique(substr(names(dwgs.perpatient), 1, 4))


# plot gains
swisnf.cnas.pp_gain <- swisnf.cnas.pp
swisnf.cnas.pp_gain[swisnf.cnas.pp_gain=="clonal.loss"] <- NA
swisnf.cnas.pp_gain[swisnf.cnas.pp_gain=="noChr"] <- NA
swisnf.cnas.pp_gain[swisnf.cnas.pp_gain=="diploid"] <- NA
swisnf.cnas.pp_gain[swisnf.cnas.pp_gain=="noChr"] <- NA
swisnf.cnas.pp_gain[swisnf.cnas.pp_gain=="subclonal.loss"] <- NA
swisnf.cnas.pp_gain[swisnf.cnas.pp_gain=="subclonal.gain.loss"] <- "subclonal.gain"

keep <- swisnf.cnas.pp_gain[rowSums(!is.na(swisnf.cnas.pp_gain), na.rm = T) > 0,]
keep$geneID <- rownames(keep)

x <- tidyr::gather(keep, key="key", value="value", -geneID)
x$type <- ifelse(x$key %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
x$pga <- scaaPP[match(x$key, scaaPP$patient),8]
x$propSCAA <- scaaPP[match(x$key, scaaPP$patient),2]

ggplot(x, aes(y=geneID, x=reorder(key, pga), fill= value)) + 
  geom_tile(color="black") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  scale_fill_manual(values = c("#CC0033",alpha("#CC0033",0.5)), na.value = "#FFFFFF") +
  theme_custom() +
  theme(axis.text.x = element_text(size=24, angle=90, vjust = 0.5), axis.title = element_blank(), axis.text.y = element_text(size=24),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "bottom", 
        strip.text = element_blank(), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA))

# plot losses
swisnf.cnas.pp_loss <- swisnf.cnas.pp
swisnf.cnas.pp_loss[swisnf.cnas.pp_loss=="clonal.gain"] <- NA
swisnf.cnas.pp_loss[swisnf.cnas.pp_loss=="noChr"] <- NA
swisnf.cnas.pp_loss[swisnf.cnas.pp_loss=="diploid"] <- NA
swisnf.cnas.pp_loss[swisnf.cnas.pp_loss=="noChr"] <- NA
swisnf.cnas.pp_loss[swisnf.cnas.pp_loss=="subclonal.gain"] <- NA
swisnf.cnas.pp_loss[swisnf.cnas.pp_loss=="subclonal.gain.loss"] <- "subclonal.loss"

keep <- swisnf.cnas.pp_loss[rowSums(!is.na(swisnf.cnas.pp_loss), na.rm = T) > 0,]
keep$geneID <- rownames(keep)

x <- tidyr::gather(keep, key="key", value="value", -geneID)
x$type <- ifelse(x$key %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
x$pga <- scaaPP[match(x$key, scaaPP$patient),8]
x$propSCAA <- scaaPP[match(x$key, scaaPP$patient),2]

ggplot(x, aes(y=geneID, x=reorder(key, pga), fill= value)) + 
  geom_tile(color="black") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  scale_fill_manual(values = c("#3300FF",alpha("#3300FF",0.5)), na.value = "#FFFFFF") +
  theme_custom() +
  theme(axis.text.x = element_text(size=24, angle=90, vjust = 0.5), axis.title = element_blank(), axis.text.y = element_text(size=24),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "bottom",
        strip.text = element_blank(), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA))

# number of swisnf cna and scaa burden
x <- data.frame(clonal.gain=apply(swisnf.cnas.pp_gain, 2, function(x) sum(x=="clonal.gain", na.rm = T)),
                subclonal.gain=apply(swisnf.cnas.pp_gain, 2, function(x) sum(x=="subclonal.gain", na.rm = T)),
                any.gain=apply(swisnf.cnas.pp_gain, 2, function(x) sum(!is.na(x), na.rm = T)))
x$pga <- scaaPP[match(rownames(x), scaaPP$patient), 8]
x$propSCAA <- scaaPP[match(rownames(x), scaaPP$patient), 2]
x$frac <- x$any.gain/259
x$WGD <- scaaPP[match(rownames(x), scaaPP$patient),6]
x$type <- scaaPP[match(rownames(x), scaaPP$patient),7]
x$patient <- rownames(x)
x <- x[which(x$patient %!in% c("C516","C542","C543")),]

x <- data.frame(clonal.gain=apply(swisnf.cnas.pp_loss, 2, function(x) sum(x=="clonal.loss", na.rm = T)),
                subclonal.gain=apply(swisnf.cnas.pp_loss, 2, function(x) sum(x=="subclonal.loss", na.rm = T)),
                any.gain=apply(swisnf.cnas.pp_loss, 2, function(x) sum(!is.na(x), na.rm = T)))
x$pga <- scaaPP[match(rownames(x), scaaPP$patient), 8]
x$propSCAA <- scaaPP[match(rownames(x), scaaPP$patient), 2]
x$frac <- x$any.gain/259
x$WGD <- scaaPP[match(rownames(x), scaaPP$patient),6]
x$type <- scaaPP[match(rownames(x), scaaPP$patient),7]
x$patient <- rownames(x)
x <- x[which(x$patient %!in% c("C516","C542","C543")),]

ggplot2::ggplot(x, aes(x=as.numeric(propSCAA), y=as.numeric(frac))) +
  geom_point(size=3, shape=21, fill="#CC0033") +
  ggtitle("Exc. C516/542/543") +
  ylab("Fraction of chromatin \nmodifiers gained") +
  xlab("PGA") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "none",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))



# NUMBER OF MUTS IN SWISNF AND SCAA BURDEN =====
swisnf <-  c("ARID1A","ARID1B","ARID2","PBRM1","SMARCA4","SMARCB1", "SMARCA2", "SMARCC1", "SMARCA5", "SMARCC2", "SMARCD1", "SMARCD2", "SMARCD3")
cm_snv <- readxl::read_excel("~/Documents/SCAA/Data/CM_SNVs.xlsx")
swisnf_snv <- cm_snv[which(cm_snv$GeneID %in% swisnf),]
swisnf_snv <- swisnf_snv[,colSums(!is.na(swisnf_snv))>0]