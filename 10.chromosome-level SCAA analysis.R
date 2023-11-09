# DATA =====
patients <- unique(scaa.long$patient)
propSCAA_chr <- list()
l <- 1
for (i in 1:length(patients)) {
  patient <- patients[i]
  for (chr in 1:22) {
    wd <- scaa.long[which(scaa.long$patient==patient & scaa.long$chr==chr),]
    propSCAA_chr[[l]] <- setNames(c(patient, chr, sum(wd$value==1)/nrow(wd)),
                                  c("patient", "chr", "propSCAA"))
    l <- l + 1
  }
}
propSCAA_chr <- data.frame(do.call(rbind, propSCAA_chr))
propSCAA_chr$propSCAAPP <- scaaPP[match(propSCAA_chr$patient, scaaPP$patient), 2]
propSCAA_chr$propSCAA <- as.numeric(propSCAA_chr$propSCAA)
propSCAA_chr$chrSCAAratio <- propSCAA_chr$propSCAA /propSCAA_chr$propSCAAPP
propSCAA_chr$MSI <- scaaPP[match(propSCAA_chr$patient, scaaPP$patient), 7]

# calculate PGA per chr for each patient
output <- list()
l <- 1
for (i in 1:length(dwgs.perpatient) ) {
  for (j in 1:22) {
    wd <- dwgs.perpatient[[i]]
    wd <- wd[which(wd$chromosome==j),]
    wd[,-c(1:3)] <- round(wd[,-c(1:3)])
    weights <- wd$stop-wd$start
    weights <- weights / sum(weights)
    
    if (ncol(wd)==4) {
      y <- weights * ifelse(wd$C527.B1_G7!=2, 1, 0)
      y <- mean(y)
    } else {
      y <- apply(wd[,-c(1:3)], 2, function(x) {
        weights * ifelse(x!=2, 1, 0)
      })
      y <- mean(colSums(y, na.rm = T))
    }
   
    if (nrow(wd)>0) {
      output[[l]] <- setNames(c(names(dwgs.perpatient)[i], j, y),
                              c("patient","chr","chr_pga")) 
      l <- l + 1
    } else {
      output[[l]] <- setNames(c(names(dwgs.perpatient)[i], j, NA),
                              c("patient","chr","chr_pga")) 
      l <- l + 1
    }
  }
}

chr_pga <- data.frame(do.call(rbind, output))

# calculate number of segmewnts per chr
output <- list()
l <- 1
for (i in 1:length(dwgs.perpatient) ) {
  for (j in 1:22) {
    wd <- dwgs.perpatient[[i]]
    wd <- wd[which(wd$chromosome==j),]
    x <- apply(wd[-c(1:3)],2, function(x) {length((rle(x))[[1]])})
    output[[l]] <- mean(x)
    l <- l + 1
  }
}

chr_pga$chr_n.segments <- unlist(output)

# calculate average copy number per chr
output <- list()
l <- 1
aveChrCN <- list()
for (i in 1:length(dwgs.perpatient) ) {
  for (j in 1:22) {
    wd <- dwgs.perpatient[[i]]
    wd <- wd[which(wd$chromosome==j),]
    weights <- wd$stop-wd$start
    sumweight <- sum(weights)
    weightCN <- apply(wd[-c(1:3)], 2, function(x) {x*weights})
    sumweightCN <- colSums(weightCN)
    chr.weightedCN <- sumweightCN/sumweight
    ave.chr.weightedCN <- mean(chr.weightedCN)
    aveChrCN[[l]] <- ave.chr.weightedCN
    l <- l + 1
  }
}

chr_pga$aveChrCN <- unlist(aveChrCN)

# GC content per chromosome (https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-019-4137-z)
GCperChr <- data.frame(chr=seq(1,22,1), 
                       GC=c(41.72,40.23,39.67,38.24,39.51,39.61,40.70,
                            40.16,41.28,41.54,41.54,40.77,38.55,40.83,
                            42.03,44.58,45.32,39.78,47.94,43.80,40.94,47))

# protein coding genes per chromosome (https://en.wikipedia.org/wiki/Human_genome)
genesperChr <- data.frame(chr=seq(1,22,1), 
                          genes=c(2058,1309,1078,752,876,1048,989,677,786,733,1298,1034,327,830,613,873,1197,270,1472,544,234,488))


# add
propSCAA_chr <- merge(propSCAA_chr, chr_pga)
propSCAA_chr$pga <- scaaPP[match(propSCAA_chr$patient, scaaPP$patient), 8]
propSCAA_chr$n.segments <- scaaPP[match(propSCAA_chr$patient, scaaPP$patient), 9]
propSCAA_chr$n.chr <- scaaPP[match(propSCAA_chr$patient, scaaPP$patient), 10]
propSCAA_chr$chrPGAratio <- as.numeric(propSCAA_chr$chr_pga) /propSCAA_chr$pga
propSCAA_chr$chrSEGratio <- (as.numeric(propSCAA_chr$chr_n.segments)-1) / (propSCAA_chr$n.segments-1)

propSCAA_chr <- propSCAA_chr[c(1,2,6,3,4,5,9,7,10,13,12,8,11,14)]
library("DescTools")
library(rtracklayer)
hg38.length <- SeqinfoForUCSCGenome("hg38")
hg38.length <- data.frame("chr"=substr(hg38.length@seqnames,4,nchar(hg38.length@seqnames)), "length"=hg38.length@seqlengths)
hg38.length <- hg38.length[c(1:22),]
propSCAA_chr$chrLength <- hg38.length[match(propSCAA_chr$chr, hg38.length$chr), 2]
propSCAA_chr$GCperChr <- GCperChr[match(propSCAA_chr$chr, GCperChr$chr), 2]
propSCAA_chr$genesperChr <- genesperChr[match(propSCAA_chr$chr, genesperChr$chr), 2]
propSCAA_chr$chr_pga <- as.numeric(propSCAA_chr$chr_pga)

#add number of TSG/onco per chr
cosmic <- read_csv('~/Documents/CNA/Github/Data/GeneLists/COSMIC 11_12_03 2022.csv')
#cosmic <- cosmic[which(cosmic$Tier=="1"),]
onco <- data.frame(table(sub("\\:.*", "", cosmic$`Genome Location`[grep("oncogene", cosmic$`Role in Cancer`)])))
names(onco) <- c("chr","onco")
tsg <- data.frame(table(sub("\\:.*", "", cosmic$`Genome Location`[grep("TSG", cosmic$`Role in Cancer`)])))
names(tsg) <- c("chr","tsg")
ntsgonco <- merge(onco,tsg)

propSCAA_chr$n.tsg <- ntsgonco[match(propSCAA_chr$chr, ntsgonco$chr), 3]
propSCAA_chr$n.onco <- ntsgonco[match(propSCAA_chr$chr, ntsgonco$chr), 2]
propSCAA_chr$n.tsg[is.na(propSCAA_chr$n.tsg)] <- 0
propSCAA_chr$n.onco[is.na(propSCAA_chr$n.onco)] <- 0


# SCAA BURDEN PER CHR =====

#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*10.35))
jpeg('tempfig.jpeg', width = (40), height = (60), units = "cm", res = 300)
#ggplot(propSCAA_chr[which(propSCAA_chr$patient %!in% c("C516","C542","C543")),], aes(x=as.numeric(chr), y=log(chrSCAAratio, base = 2))) +
ggplot(propSCAA_chr, aes(x=as.numeric(chr), y=log(chrSCAAratio, base = 2))) +
  #geom_jitter(shape=21, position=position_jitter(0), fill="#333399") +
  geom_line(color="#333399") +
  geom_jitter(shape=21, size=5, position=position_jitter(0), aes(fill=log(chrSCAAratio, base = 2))) +
  scale_fill_gradient2(low = "blue", mid="yellow", high = "red") +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_x_continuous(breaks=seq(1,22,2), name = "Chromosome") +
  scale_y_continuous(name = "SCAA burden per chromosome, relative to per patient SCAA burden (logFold2)") +
  #facet_grid(rows = vars(patient)) +
  facet_wrap(~patient, ncol = 3) +
  theme_custom() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(fill=NA), strip.text = element_text(size=22))
dev.off()

library(stats)
library(tidyverse)
library(ggpubr)
library(rstatix)
pwc <- propSCAA_chr %>%
  pairwise_t_test(chrSCAAratio ~ chr, p.adjust.method = "bonferroni")
pwc <- pwc[which(pwc$p.adj.signif <0.05),]

# relationship between rel SCAA ratio and SCAA burden
propSCAA_chr$dist <- abs(1-propSCAA_chr$chrSCAAratio)
x <- aggregate(propSCAA_chr$dist, list(propSCAA_chr$patient), FUN=mean)
x$scaaPP <- propSCAA_chr[match(x$Group.1, propSCAA_chr$patient), 5]
x$note <- ifelse(x$Group.1 %in% c("C516","C542","C543"), "TRUE","FALSE")

jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x[which(x$note==FALSE),], aes(x=log(scaaPP), y=(x))) +
  #geom_jitter(shape=21, position=position_jitter(0), fill="#333399") +
  #scale_x_continuous(name = "SCAA burden") +
  scale_x_continuous(name="SCAA burden", 
                     breaks = c(-10,-7.5,-5,-2.5,0),
                     labels = c(round(exp(-10),5),round(exp(-7.5),4),round(exp(-5),2),round(exp(-2.5),2),round(exp(0),2))) +
  scale_y_continuous(name = "Mean fold change in \nSCAA burden per chromosome") +
  scale_fill_manual(values = c(alpha("#333399",0.85), alpha("#CC9966",0.8))) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=9, ) +
  
  geom_jitter(shape=21, position=position_jitter(0), aes(fill=note), size=5) +
  theme_custom() +
  theme(legend.position = "none", plot.margin = unit(c(0.5,0.5,0,0), "cm"))
dev.off()

# MODEL CHR SCAA BURDEN =====
# do those chr with abnoramlly high SCAA for the patient, have abnormally high PGA for the patient?
ggplot(propSCAA_chr[which(propSCAA_chr$patient %!in% c("C516","C542","C543")),], aes(x=log(chrSCAAratio), y=log(chrPGAratio))) +
  geom_point(shape=21, aes(fill=chr)) +
  geom_text_repel(aes(label=chr)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, )


library(MASS)
library(caret)

cor.test(propSCAA_chr$chr_pga, propSCAA_chr$chr_n.segments)
cor.test(propSCAA_chr$chr_pga, propSCAA_chr$aveChrCN)
cor.test(propSCAA_chr$chr_n.segments, propSCAA_chr$aveChrCN)
cor.test(propSCAA_chr$genesperChr, propSCAA_chr$GCperChr)
cor.test(scaaPP$pga, scaaPP$n.segments)

# variance in scaa and pga per chrom
data <- propSCAA_chr[which(propSCAA_chr$patient %!in% c("C516","C542","C543") & propSCAA_chr$MSI=="MSS"),]
scaa.chr_var <- data.frame(aggregate(data$chrSCAAratio, list(data$chr), FUN=sd))
colnames(scaa.chr_var) <- c("chr","scaa.var")
pga.chr_var <- data.frame(aggregate(data$chrPGAratio, list(data$chr), FUN=sd))
colnames(pga.chr_var) <- c("chr","pga.var")
x <- merge(scaa.chr_var,pga.chr_var)
cor.test(x$scaa.var, x$pga.var)

# all
data <- propSCAA_chr[which(propSCAA_chr$patient %!in% c("C516","C542","C543")) , colnames(propSCAA_chr) %in% c("chrSCAAratio","chr_n.segments","chr_pga", "chrLength", "GCperChr","genesperChr","chr","n.tsg","n.onco")]
data$chr <- as.factor(data$chr)
data$chr <- factor(data$chr, levels = c(2, 1, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22))
fullMod <- lm(chrSCAAratio ~0+., data=data)
summary(fullMod)
stepMod1 <- stepAIC(fullMod, direction="both")
summary(stepMod1)
vif <- data.frame(car::vif(stepMod1))


data <- propSCAA_chr[which(propSCAA_chr$patient %!in% c("C516","C542","C543")) , colnames(propSCAA_chr) %in% c("chrSCAAratio","chr_pga","chr_n.segments","aveChrCN", "chrLength", "GCperChr","genesperChr","chr","chrPGAratio")]
data$predicted <- predict(stepMod1, newdata = data)
data <- data[which(data$chr!=20),]
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.8))
ggplot(data, aes(x=chrSCAAratio, y=predicted)) +
  geom_point(shape=21, aes(fill=chr)) +
  geom_text_repel(aes(label=chr)) +
  ylab("Predictied chrSCAA ratio") +
  xlab("Actual chrSCAA ratio") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none")
dev.off()

par(mfrow = c(2, 2))
plot(stepMod1)

# remove chr
data <- propSCAA_chr[which(propSCAA_chr$patient %!in% c("C516","C542","C543")) , colnames(propSCAA_chr) %in% c("chrSCAAratio","chr_n.segments","chr_pga", "chrLength", "GCperChr","genesperChr","n.tsg","n.onco")]
fullMod <- lm(chrSCAAratio ~0+., data=data)
summary(fullMod)
stepMod2 <- stepAIC(fullMod, direction="both")
summary(stepMod2)

data <- propSCAA_chr[which(propSCAA_chr$patient %!in% c("C516","C542","C543")) , colnames(propSCAA_chr) %in% c("chr", "chrSCAAratio","chr_n.segments","chr_pga", "chrLength", "GCperChr","genesperChr","n.tsg","n.onco")]
data$predicted <- predict(stepMod2, newdata = data)
#data <- data[which(data$chr!=20),]
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.8))
ggplot(data, aes(x=chrSCAAratio, y=predicted)) +
  geom_point(shape=21, aes(fill=chr)) +
  geom_text_repel(aes(label=chr)) +
  ylab("Predictied chrSCAA ratio") +
  xlab("Actual chrSCAA ratio") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  stat_smooth(method = "glm", color='black') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=8, ) +
  theme_custom() +
  theme(legend.position = "none")
dev.off()


par(mfrow = c(2, 2))
plot(stepMod2)

# just using CIN
data <- propSCAA_chr[which(propSCAA_chr$patient %!in% c("C516","C542","C543")) , colnames(propSCAA_chr) %in% c("chrSCAAratio","chr_n.segments","chr_pga")]
fullMod <- lm(chrSCAAratio ~0+., data=data)
summary(fullMod)

par(mfrow = c(2, 2))
plot(fullMod)





