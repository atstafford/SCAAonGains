
# DATA =====
peakAnn_all <- readRDS("~/Documents/SCAA/Data/peakAnn_all.rds")
promoter <- peakAnn_all$peak[which(peakAnn_all$Promoter==TRUE)]
enhancer <- peakAnn_all$peak[which(!is.na(peakAnn_all$elementType))]

patients <- names(dwgs.perpatient)
scaaPP <- list()
scaaPP_promo <- list()
scaaPP_notPromo <- list()
for (i in 1:length(patients)) {
  wd <- scaa.long[scaa.long$patient==patients[i],]
  scaaPP[[i]] <- sum(wd$value==1, na.rm = T) / nrow(wd)
  wd2 <- wd[which(wd$peak %in% promoter),]
  scaaPP_promo[[i]] <- sum(wd2$value==1, na.rm = T) / nrow(wd2)
  wd3 <- wd[which(wd$peak %!in% promoter),]
  scaaPP_notPromo[[i]] <- sum(wd3$value==1, na.rm = T) / nrow(wd3)
}
scaaPP <- data.frame(patient=patients, propSCAA=unlist(scaaPP), propSCAA.promo=unlist(scaaPP_promo), propSCAA.notPromo=unlist(scaaPP_notPromo))
scaaPP$meanPloidy <- unlist(ploidy.perpatient[match(scaaPP$patient, ploidy.perpatient$patient), 2])
scaaPP$WGD <- ifelse(scaaPP$meanPloidy>2.5, ">2.5", "<2.5")
scaaPP$type <- ifelse(scaaPP$patient %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")

# calculate PGA per sample and average for each patient
output <- list()
for (i in 1:length(dwgs.perpatient) ) {
  wd <- dwgs.perpatient[[i]]
  wd[,-c(1:3)] <- round(wd[,-c(1:3)])
  weights <- wd$stop-wd$start
  weights <- weights / sum(weights)
  
  y <- apply(wd[,-c(1:3)], 2, function(x) {
    weights * ifelse(x!=2, 1, 0)
  })
  
  output[[i]] <- mean(colSums(y, na.rm = T))
}

pga <- data.frame(patient=names(dwgs.perpatient), pga=unlist(output))

# calculate number of segmewnts pp
output <- list()
for (i in 1:length(dwgs.perpatient) ) {
  wd <- dwgs.perpatient[[i]]
  wd$x <- rowMeans(wd[,-c(1:3)])
  output[[i]] <- length((rle(wd$x))[[1]])
}

pga$n.segments <- unlist(output)

scaaPP$pga <- pga[match(scaaPP$patient, pga$patient),2]
scaaPP$propSCAA[is.na(scaaPP$propSCAA)] <- 0
scaaPP$n.segments <- pga[match(scaaPP$patient, pga$patient),3]
scaaPP$pga.ploidyRC <- pga[match(scaaPP$patient, pga.ploidyRC$patient),2]

# SCAA PROPORTION AND WGD/MSI =====
# scatter: Y=proportion of SCAAs, X=ploidy mean per patient
# finding: WGD patients have similar numbers of SCAAS per patient, background SCAA prop is highly patient specific
library(scales)
p1 <- ggplot2::ggplot(data=scaaPP, aes(y=log(propSCAA), x=as.factor(WGD), fill=as.factor(WGD))) +
  geom_violin(trim = F, alpha=0.4) +
  geom_jitter(size=3, shape=21, position=position_jitter(0.2, seed=1), alpha=0.8) +
  geom_boxplot(width=0.1, outlier.size = 0) +
  scale_fill_manual(values = c("#E6AB02","#660033")) +
  xlab("Average ploidy per patient") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  geom_text_repel(aes(label=patient), position=position_jitter(0.2, seed=1),min.segment.length = 0, size=5) +
  theme_custom() +
  theme(legend.position = "none")

p2 <- ggplot2::ggplot(data=scaaPP, aes(y=log(propSCAA), x=as.factor(type), fill=as.factor(type))) +
  geom_violin(trim = F, alpha=0.4) +
  geom_jitter(size=3, shape=21, position=position_jitter(0.2, seed=1), alpha=0.8) +
  geom_boxplot(width=0.1, outlier.size = 0) +
  scale_fill_manual(values = c("#A6CEE3","#666633")) +
  xlab("MSI status") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  geom_text_repel(aes(label=patient), position=position_jitter(0.2, seed=1),min.segment.length = 0, size=5) +
  theme_custom() +
  theme(legend.position = "none")

jpeg('tempfig.jpeg', width = (3*37.795*4.30),  height = (3*37.795*8.38))
cowplot::plot_grid(p1, p2, ncol = 1)
dev.off()

wilcox.test((z$log[which(z$WGD=="MSI")]), 
            (z$log[which(z$WGD=="MSS")]) )

wilcox.test((z$log[which(z$WGD=="<2.5")]), 
            (z$log[which(z$WGD==">2.5")]) )

# PGA VS N.SEGMENTS =====
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data=scaaPP, aes(x=pga, y=n.segments, fill=type, color=type)) +
  geom_point(size=3, shape=21, color="black") +
  geom_text_repel(aes(label=patient), color="black") +
  stat_smooth(method = "glm", color='black') +
  scale_fill_manual(values = c("#1F78B4", "#33A02C")) +
  scale_color_manual(values = c("#1F78B4", "#33A02C")) +
  ylab("Number of segments") +
  xlab("PGA") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom()
dev.off()



# SCAA PROPORTION AND PGA CORRELATE =====

# all patients
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data=scaaPP, aes(x=log(pga), y=log(as.numeric(propSCAA)))) +
  geom_point(size=3, color="#39568CFF") +
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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 
dev.off()

# exc 3
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data=scaaPP[which(scaaPP$patient %!in% c("C543","C542","C516")),], aes(x=log(pga), y=log(as.numeric(propSCAA)))) +
  geom_point(size=3, color="#39568CFF") +
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
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 
dev.off()

# exc 3 + sep by MSI
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data=scaaPP[which(scaaPP$patient %!in% c("C543","C542","C516")),], aes(x=(pga.ploidyRC), y=log(as.numeric(propSCAA)), fill=type, color=type)) +
  geom_point(size=3, aes(shape=WGD), color="black") +
  ggtitle("Exc. C516/542/543") +
  scale_color_manual(values=c("#6699CC","#666633")) +
  scale_fill_manual(values=c("#6699CC","#666633")) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  #scale_x_continuous(name="PGA", 
  #                   breaks = c(-10,-7.5,-5,-2.5,0),
  #                   labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  stat_smooth(method = "glm") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 
dev.off()


# SCAA PROPORTION AND NUMBER OF SEGMENTS CORRELATE =====

# all patients
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data=scaaPP, aes(x=(n.segments), y=log(as.numeric(propSCAA)))) +
  geom_point(size=3, color="#39568CFF") +
  ggtitle("All patients") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_x_continuous(name="Number of segments") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 
dev.off()

# exc 3
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data=scaaPP[which((scaaPP$patient %!in% c("C543","C542","C516"))),], aes(x=(n.segments), y=log(as.numeric(propSCAA)))) +
  geom_point(size=3, color="#39568CFF") +
  ggtitle("Exc. C516/542/543") +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_x_continuous(name="Number of segments") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 
dev.off()

# exc 3 + sep by MSI
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data=scaaPP[which(scaaPP$patient %!in% c("C543","C542","C516")),], aes(x=n.segments, y=log(as.numeric(propSCAA)), fill=type, color=type)) +
  geom_point(size=3, shape=21, color="black") +
  ggtitle("Exc. C516/542/543") +
  scale_color_manual(values=c("#6699CC","#666633")) +
  scale_fill_manual(values=c("#6699CC","#666633")) +
  scale_y_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_x_continuous(name="Number of segments", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  stat_smooth(method = "glm") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 
dev.off()

# NO SCAA INDIVIDUALLY CORRELATES WITH PGA =====

# pull all peaks which are SCAAs
scaa.peaks <- rownames(scaa[rowSums(scaa)>1,])
scaa.long$pga <- scaaPP[match(scaa.long$patient, scaaPP$patient), 8]

# for each scaa.peak, pull pag of patients with and without scaa
test <- list()
for (i in 1:length(scaa.peaks)) {
  print(i)
  wd <- scaa.long[which(scaa.long$peak==scaa.peaks[i]),]
  wd$pga <- scaaPP[match(wd$patient, scaaPP$patient), 8]
  test[[i]] <- wilcox.test(log(wd$pga[which(wd$value==1)]), log(wd$pga[which(wd$value==0)]))$p.value
}

x <- data.frame(peak=scaa.peaks, pvalue=unlist(test))
x$FDR <- p.adjust(x$pvalue, "fdr")


# CNAS IN CHROMATIN MODIFIERS (WITHOUT PLOIDY RECENTRE) =====

pgaplot <- ggplot2::ggplot(scaaPP, aes(x=reorder(patient, pga), group=1)) +
  geom_area(aes(y=pga), fill="#69b3a2", alpha=0.4) +
  geom_line(aes(y=pga), color="#69b3a2", size=2) +
  geom_point(aes(y=pga), size=3, color="#69b3a2") +
  ylab("PGA") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  theme_custom() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(vjust = -10),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "top", 
        strip.text = element_text(size=28), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

scaaplot <- ggplot2::ggplot(scaaPP, aes(x=reorder(patient, pga), group=1)) +
  geom_bar(aes(y=(propSCAA)), fill=alpha("#FF9933",0.6), color="#FF9933", stat="identity") +
  scale_y_continuous(name="SCAA \nburden") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  theme_custom() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(vjust = -10),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "top", 
        strip.text = element_blank(), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

# create chromatin gene list
chromatinGenes <- read.table('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/chromatinGenes', header = T, sep = "\t")
DNAmodGenes <- read.table('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/DNAmodGenes', header = T, sep = "\t")
histoneGenes <- read.table('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/histoneGenes', header = T, sep = "\t")

library(tidyverse)
hallmarkGenes <- readRDS('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/hallmarkGenes.rds')
cosmic <- read_csv('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/COSMIC 11_12_03 2022.csv')

GOlist <- read.delim("~/Documents/CNA/Data/hg38.1_all_gene_GO_annotations.txt")
GOlist <- GOlist[GOlist$Chromosome.scaffold.name %in% c(seq(1,22,1),"X","Y"),]

chromatin.ass.genes <- unique(c(chromatinGenes$Symbol, DNAmodGenes$Symbol, histoneGenes$Symbol))
#chromatin.ass.genes <- unique(c(hallmarkGenes$gene_symbol, cosmic$`Gene Symbol`))
chromatin.ass.genes <- GOlist[GOlist$HGNC.symbol %in% chromatin.ass.genes,]
chromatin.ass.genes <- chromatin.ass.genes[,c(1,3,4,5,8)]
chromatin.ass.genes <- chromatin.ass.genes[!duplicated(chromatin.ass.genes),]

# per patient see if have clonal/subclonal gain/loss in each chromatin gene
output <- list()
l <- 1
for (i in 1:length(segments)) { # for this patient...
  wd <- segments[[i]]
  
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
    
    overlapIndex <- apply(wd2, 1, function(x) DescTools::Overlap(as.numeric(c(x[6], x[7])), as.numeric(c(geneStart,geneStop)))) > 0.9*(geneStop-geneStart)
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
library(ggplot2)

# plot gains
chromatin.cnas.pp_gain <- chromatin.cnas.pp
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="clonal.loss"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="noChr"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="diploid"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="noChr"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="subclonal.loss"] <- NA
chromatin.cnas.pp_gain[chromatin.cnas.pp_gain=="subclonal.gain.loss"] <- "subclonal.gain"

keep <- chromatin.cnas.pp_gain[rowSums(!is.na(chromatin.cnas.pp_gain), na.rm = T) > 16,]
keep$geneID <- rownames(keep)

x <- tidyr::gather(keep, key="key", value="value", -geneID)
x$type <- ifelse(x$key %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
x$pga <- scaaPP[match(x$key, scaaPP$patient),8]
x$propSCAA <- scaaPP[match(x$key, scaaPP$patient),2]

gainplot <- ggplot(x, aes(y=geneID, x=reorder(key, pga), fill= value)) + 
  geom_tile(color="black") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  scale_fill_manual(values = c("#CC0033",alpha("#CC0033",0.5)), na.value = "#FFFFFF") +
  theme_custom() +
  theme(axis.text.x = element_text(size=24, angle=90, vjust = 0.5), axis.title = element_blank(), axis.text.y = element_text(size=24),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "bottom", 
        strip.text = element_blank(), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA))

jpeg('tempfig.jpeg', width = (3*37.795*7), height = (3*37.795*9.07))
cowplot::plot_grid(pgaplot, scaaplot, gainplot, ncol = 1, align = 'v', rel_heights = c(0.12, 0.1, 0.6))
dev.off()

# plot losses
chromatin.cnas.pp_loss <- chromatin.cnas.pp
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="clonal.gain"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="noChr"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="diploid"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="noChr"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="subclonal.gain"] <- NA
chromatin.cnas.pp_loss[chromatin.cnas.pp_loss=="subclonal.gain.loss"] <- "subclonal.loss"

keep <- chromatin.cnas.pp_loss[rowSums(!is.na(chromatin.cnas.pp_loss), na.rm = T) > 4,]
keep$geneID <- rownames(keep)

x <- tidyr::gather(keep, key="key", value="value", -geneID)
x$type <- ifelse(x$key %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")
x$pga <- scaaPP[match(x$key, scaaPP$patient),8]
x$propSCAA <- scaaPP[match(x$key, scaaPP$patient),2]

lossplot <- ggplot(x, aes(y=geneID, x=reorder(key, pga), fill= value)) + 
  geom_tile(color="black") +
  facet_grid(cols = vars(type), scales = "free", space = "free") +
  scale_fill_manual(values = c("#3300FF",alpha("#3300FF",0.5)), na.value = "#FFFFFF") +
  theme_custom() +
  theme(axis.text.x = element_text(size=24, angle=90, vjust = 0.5), axis.title = element_blank(), axis.text.y = element_text(size=24),
        legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "bottom",
        strip.text = element_blank(), panel.spacing = unit(5, "pt"), panel.border = element_rect(size=1, fill=NA))

jpeg('tempfig.jpeg', width = (3*37.795*7), height = (3*37.795*9.07))
cowplot::plot_grid(pgaplot, scaaplot, lossplot, ncol = 1, align = 'v', rel_heights = c(0.12, 0.1, 1))
dev.off()

# NUMBER OF GAINS/LOSSES IN CM WITH PGA AND SCAA BURDEN =====
x <- data.frame(clonal.gain=apply(chromatin.cnas.pp_gain, 2, function(x) sum(x=="clonal.gain", na.rm = T)),
                subclonal.gain=apply(chromatin.cnas.pp_gain, 2, function(x) sum(x=="subclonal.gain", na.rm = T)),
                any.gain=apply(chromatin.cnas.pp_gain, 2, function(x) sum(!is.na(x), na.rm = T)))
x$pga <- scaaPP[match(rownames(x), scaaPP$patient), 8]
x$propSCAA <- scaaPP[match(rownames(x), scaaPP$patient), 2]
x$frac <- x$any.gain/259
x$WGD <- scaaPP[match(rownames(x), scaaPP$patient),6]
x$type <- scaaPP[match(rownames(x), scaaPP$patient),7]

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot2::ggplot(x, aes(x=as.numeric(pga), y=as.numeric(frac))) +
  geom_point(size=3, shape=21, fill="#CC0033") +
  ylab("Fraction of chromatin \nmodifiers gained") +
  xlab("PGA") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "none", 
        strip.text = element_text(size=28), panel.spacing = unit(2, "lines"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
dev.off()

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot2::ggplot(x, aes(x=log(as.numeric(propSCAA)), y=as.numeric(frac))) +
  geom_point(size=3, shape=21, fill="#CC0033") +
  ylab("Fraction of chromatin \nmodifiers gained") +
  scale_x_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "none", 
        strip.text = element_text(size=28), panel.spacing = unit(2, "lines"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
dev.off()


x <- data.frame(clonal.loss=apply(chromatin.cnas.pp_loss, 2, function(x) sum(x=="clonal.loss", na.rm = T)),
                subclonal.loss=apply(chromatin.cnas.pp_loss, 2, function(x) sum(x=="subclonal.loss", na.rm = T)),
                any.loss=apply(chromatin.cnas.pp_loss, 2, function(x) sum(!is.na(x), na.rm = T)))
x$pga <- scaaPP[match(rownames(x), scaaPP$patient), 8]
x$propSCAA <- scaaPP[match(rownames(x), scaaPP$patient), 2]
x$frac <- x$any.loss/259
x$WGD <- scaaPP[match(rownames(x), scaaPP$patient),6]
x$type <- scaaPP[match(rownames(x), scaaPP$patient),7]

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot2::ggplot(x, aes(x=as.numeric(pga), y=as.numeric(any.loss)/nrow(x))) +
  geom_point(size=3, shape=21, fill="#3300FF") +
  ylab("Fraction of chromatin \nmodifiers lost") +
  xlab("PGA") +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "none", 
        strip.text = element_text(size=28), panel.spacing = unit(2, "lines"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
dev.off()

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot2::ggplot(x, aes(x=log(as.numeric(propSCAA)), y=as.numeric(frac))) +
  geom_point(size=3, shape=21, fill="#3300FF") +
  ylab("Fraction of chromatin \nmodifiers lost") +
  scale_x_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     limits = c(-10,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.title = element_blank(), legend.text = element_text(size=28), legend.position = "none", 
        strip.text = element_text(size=28), panel.spacing = unit(2, "lines"), panel.border = element_rect(size=1, fill=NA),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
dev.off()

# NUMBER OF SNVS CM WITH PGA AND SCAA BURDEN =====

# c=clonal, sc=subclonal, it=indel/trunctating, snv=snv
cm_snv <- readxl::read_excel("~/Documents/SCAA/Data/CM_SNVs.xlsx")
x <- data.frame(frac=(colSums(!is.na(cm_snv[-1])))/117)
x$pga <- scaaPP[match(rownames(x), scaaPP$patient), 8]
x$propSCAA <- scaaPP[match(rownames(x), scaaPP$patient), 2]
x$patient <- rownames(x)
x$ploidy <- scaaPP[match(x$patient, scaaPP$patient), 6]
x$MSI <- scaaPP[match(x$patient, scaaPP$patient), 7]

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.35))
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

# ======
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


