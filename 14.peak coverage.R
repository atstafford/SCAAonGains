# PREREQUISITIS: load section 4 data 

# DATA =====

peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")]

peaks <- data.frame(peak=unique(rownames(scaa)))
peaks$chr <- sub("\\:.*", "", peaks$peak)

library(data.table)
centromere <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", 
                    col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
centromere <- centromere[ , .(length = sum(chromEnd - chromStart)), 
                          by = .(chrom, arm = substring(name, 1, 1)) ]
centromere <- centromere[centromere$arm=="p"]

library("DescTools")
library(rtracklayer)
hg38.length <- SeqinfoForUCSCGenome("hg38")
hg38.length <- data.frame("chr"=substr(hg38.length@seqnames,4,nchar(hg38.length@seqnames)), "length"=hg38.length@seqlengths)
hg38.length <- hg38.length[c(1:22),]

peaks$centromere <- unlist(centromere[match(peaks$chr, centromere$chrom),3])
peaks$chr <- as.numeric(substr(peaks$chr,4,nchar(peaks$chr)))
peaks$mid <- (as.numeric(sub(".*-", "", peaks$peak)) - as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", peaks$peak)))/2 + as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", peaks$peak))
peaks$chr.start <- 1
peaks$chr.end <- hg38.length[match(peaks$chr, hg38.length$chr),2]

# NUMBER OF SEGMENTS PP =====

# pull unique segments per patient
x <- do.call(rbind, segments)
x <- x[ ,colnames(x) %in% c("patient","segID")]
x <- x[!duplicated(x),]
length(unique(x$patient))

# number of segments per person
segPP <- data.frame(table(x$patient))
mean(segPP$Freq)
min(segPP$Freq)
max(segPP$Freq)

# N CHANGES WITH CLONALITY CUTOFF OF SEGMENTS ====

# bar chart showing how n chages with clonality cutoff
x <- do.call(rbind, segments)
x <- x[ ,colnames(x) %in% c("patient","segID","cna","clonalityCNA")]
x <- x[!duplicated(x),]

dip0 <- nrow(x[which(x$cna=="diploid" & x$clonality>0),])
dip0.5 <- nrow(x[which(x$cna=="diploid" & x$clonality>0.5),])
dip0.7 <- nrow(x[which(x$cna=="diploid" & x$clonality>0.7),])

gain0 <- nrow(x[which(x$cna=="gain" & x$clonality>0),])
gain0.5 <- nrow(x[which(x$cna=="gain" & x$clonality>0.5),])
gain0.7 <- nrow(x[which(x$cna=="gain" & x$clonality>0.7),])

loss0 <- nrow(x[which(x$cna=="loss" & x$clonality>0),])
loss0.5 <- nrow(x[which(x$cna=="loss" & x$clonality>0.5),])
loss0.7 <- nrow(x[which(x$cna=="loss" & x$clonality>0.7),])

data <- data.frame(category=c('dip','dip','dip','gain','gain','gain','loss','loss','loss'), 
                   clonality=c('>0','>0.5','>0.7','>0','>0.5','>0.7','>0','>0.5','>0.7'),
                   count=c(dip0, dip0.5, dip0.7, gain0, gain0.5, gain0.7, loss0, loss0.5, loss0.7))

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data, aes(x=(category), y=count, fill=clonality)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c(alpha("#39568CFF",0.5),  alpha("#39568CFF",0.7), alpha("#39568CFF",1))) +
  ylab("Count") +
  xlab("Segment category") +
  theme_custom() +
  theme(legend.text = element_text(size=24), legend.title = element_text(size=24))
dev.off()



# DISTRIBUTION OF PEAKS ACROSS GENOME AND PEAK GAP REGIONS =====

# line: Y=binary display of peak, X=genome region
# finding: roughly even, missing peaks on some small arms of acrocentric. 

library("DescTools")
library(rtracklayer)
hg38.length <- SeqinfoForUCSCGenome("hg38")
hg38.length <- data.frame("chr"=substr(hg38.length@seqnames,4,nchar(hg38.length@seqnames)), "length"=hg38.length@seqlengths)
hg38.length <- hg38.length[c(1:22),]

wd2 <- tidyr::gather(peaks[-1], key="key", value="value", -chr)
wd2$key[wd2$key=="chr.start"] <- "chr.end"

#jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*7.87))
jpeg('tempfig.jpeg', width = (40), height = (40), units="cm", res=300)
ggplot() +
  geom_hline(yintercept = c(seq(1,22,1)), color="#666666") +
  geom_point(data=wd2[wd2$key!="mid",], aes(y=chr, x=value, shape=key), color="#CC3300") +
  scale_shape_manual(values=c(16,15)) +
  geom_point(data=wd2[wd2$key=="mid",], aes(y=chr, x=value), color=alpha("#660099",0.5), shape="|") +
  xlab("Genomic loci") +
  scale_y_continuous(name="Chromosome", breaks = seq(1,22,1), labels = seq(1,22,1)) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=24), legend.title = element_blank()) 
dev.off()

# FIND LARGE PEAK GAPS =====

output <- list()
l <- 1
for (i in 1:22) {
  chr <- i
  x <- peaks[peaks$chr==chr,]
  x <- x[order(x$mid,decreasing=F),]
  y <- c(1,unlist(x$mid),x$chr.end[1])
  gapIndex <- which(diff(y)>8000000)
  if (length(gapIndex)==0) {
    next
  }
  for (j in 1:length(gapIndex)) {
    output[[l]] <- setNames(c(x$chr[gapIndex[j]], y[gapIndex[j]], y[gapIndex[j]+1]),
                          c("chr","gapStart","gapStop"))
    l <- l + 1
  }
}
peakGaps <- data.frame(do.call(rbind,output))
peakGaps$gap <- peakGaps$gapStop - peakGaps$gapStart

# overlap of peakGaps with hg38gaps
gap <- read.csv("~/Documents/SCAA/Data/UCSC_search/gap.csv")

gapOverlap <- list()
gapType <- list()
for (i in 1:nrow(peakGaps)) {
  print(i)
  wd <- gap[which(gap$chrom == paste("chr",peakGaps$chr[i],sep="")),]
  overlapIndex <- apply(wd, 1, function(x) DescTools::Overlap(as.numeric(c(x[3], x[4])), as.numeric(c(peakGaps$gapStart[i],peakGaps$gapStop[i])))) 
  if (length(overlapIndex)>0) {
    gapOverlap[[i]] <- sum(overlapIndex, na.rm = T) / peakGaps$gap[i]
    gapType[[i]] <- wd$type[which(overlapIndex == max(overlapIndex))]
  } else { # no overlaps
    gapOverlap[[i]] <- FALSE
    gapType[[i]] <- FALSE
  }
}
  
peakGaps$gapOverlap <- unlist(gapOverlap)
peakGaps$gapType <- unlist(gapType)


# STILL SEGMENTS MISSING PEAKS WITH LARGE GAP REMOVAL =====

# so as to not double count because some segments can be a mix of g/l/d per patient, only consider segments with >0.5 clonality
x <- do.call(rbind, segments)
x <- x[, colnames(x) %in% c("patient","chr","cnaStart","cnaStop","size","segID","cna","clonalityCNA","n.peaks","pcArm")]
x <- x[x$clonalityCNA > 0.5,]
x <- x[!duplicated(x),]

# ignore LOH
x <- x[which(x$cna!="LOH"),]
x$group <- ifelse(x$n.peaks>0, "withPeaks", "withoutPeaks")
x$coverage <- x$n.peaks/x$size

# remove segments which overlap with the regions missing peaks due to 'other reasons' (acro/centro)
peakGapOverlap <- list()
l <- 1
for (i in 1:nrow(x)) {
  print(i)
  for (j in 1:nrow(peakGaps)) {
    if ( ((Overlap(c(x$cnaStart[i], x$cnaStop[i]), c(peakGaps$gapStart[j], peakGaps$gapStop[j])) / x$size[i]) > 0) & x$chr[i]==peakGaps$chr[j] ) {
      peakGapOverlap[[i]] <- TRUE
      break
    } else { 
      peakGapOverlap[[i]] <- FALSE
    }
  }
}
x$peakGapOverlap <- unlist(peakGapOverlap)
#x <- x[x$peakGapOverlap==FALSE,]

# number of segments with no peaks
nrow(x[which(x$cna=="diploid" & x$n.peaks==0),])/nrow(x[which(x$cna=="diploid"),])
nrow(x[which(x$cna=="loss" & x$n.peaks==0),])/nrow(x[which(x$cna=="loss"),])
nrow(x[which(x$cna=="gain" & x$n.peaks==0),])/nrow(x[which(x$cna=="gain"),])
nrow(x[which(x$n.peaks==0),])/nrow(x)

# find more overlaps with gaps
gapOverlap <- list()
gapType <- list()
for (i in 1:nrow(x)) {
  print(i)
  wd <- gap[which(gap$chrom == paste("chr",x$chr[i],sep="")),]
  overlapIndex <- apply(wd, 1, function(y) DescTools::Overlap(as.numeric(c(y[3], y[4])), as.numeric(c(x$cnaStart[i],x$cnaStop[i])))) 
  if (sum(overlapIndex)>0) {
    gapOverlap[[i]] <- sum(overlapIndex, na.rm = T) / x$size[i]
    gapType[[i]] <- wd$type[which(overlapIndex == max(overlapIndex))][1]
  } else { # no overlaps
    gapOverlap[[i]] <- FALSE
    gapType[[i]] <- FALSE
  }
}

x$gapOverlap <- unlist(gapOverlap)
x$gapType <- unlist(gapType)

ggplot(data=x, aes(y=log(gapOverlap), x=(group))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.2)) +
  geom_jitter(size=2, shape=21, alpha=.4) +
  geom_boxplot(width = 0.3, alpha = .8,  show.legend = FALSE,position=position_dodge(0.2))


y <- x[which(x$group=="withoutPeaks" & x$peakGapOverlap==FALSE),]
sum(y$gapOverlap==0)/nrow(y)

# PROPORTION OF SEGMENTS WITH AT LEAST 1 PEAK VARIES BY SIZE =====

# histogram showing number of peaks per segment
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(x[which(x$gapType==FALSE),], aes(x=(n.peaks))) + 
  geom_histogram(alpha=0.8, position="identity", color="black", fill="#39568CFF", bins = 50) +
  #scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) +
  ylab("Count") +
  xlab("Number of peaks per segment") +
  #facet_grid(rows=vars(cna), scales = "free") +
  theme_custom() +
  theme(legend.text = element_text(size=24),
        legend.title = element_blank(),strip.text = element_text(size=24),
        panel.spacing = unit(1, "cm"))
dev.off()

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(x, aes(x=log(n.peaks), fill=cna)) + 
  geom_histogram(alpha=0.5, position="identity", color="black", 
                 bins = 50) +
  scale_fill_manual(values=c("#440154FF", "#FDE725FF", "#39568CFF")) +
  ylab("Fraction") +
  xlab("Number of peaks per segment (log)") +
  facet_grid(rows=vars(cna)) +
  theme_custom() +
  theme(legend.text = element_text(size=24),
        legend.title = element_blank())
dev.off()


# PEAKS MISSING BY SIZE =====

# histogram showing number of peaks per segment with/without peaks by size
z <- x[which(x$gapType==FALSE),]

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot() + 
  geom_histogram(data=z, aes(x=(size), fill="allSegments"), alpha=0.5, position="identity", color="black", bins = 50) +
  geom_histogram(data=z[z$group=="withPeaks",], aes(x=(size), fill="withPeaks"), alpha=0.5, position="identity", color="black", bins = 50) +
  scale_fill_manual(values=c("#440154FF","#FDE725FF"))+
  ylab("Count") +
  xlab("Size of segment") +
  theme_custom() +
  theme(legend.text = element_text(size=24),
        legend.title = element_blank(),
        legend.position = "none") +
  guides(fill=guide_legend())
dev.off()

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot() + 
  geom_histogram(data=x, aes(x=log(size), fill="allSegments"), alpha=0.5, position="identity", color="black", bins = 50) +
  geom_histogram(data=x[x$group=="withPeaks",], aes(x=log(size), fill="withPeaks"), alpha=0.5, position="identity", color="black", bins = 50) +
  scale_fill_manual(values=c("#440154FF","#FDE725FF"))+
  ylab("Count") +
  xlab("Size of segment (log)") +
  theme_custom() +
  theme(legend.text = element_text(size=24),
        legend.title = element_blank(),
        legend.position = "none") +
  guides(fill=guide_legend())
dev.off()

# does size predict proportion of peaks?
# need to categorise size to output a proportion, but n will be different for each
size <- c(seq(10000,100000,10000),seq(200000,1000000,100000),seq(2000000,150000000,1000000))
gain <- list()
loss <- list()
diploid <- list()
for (i in 2:length(size)) {
  wd <- x[x$size>size[i-1] & x$size<=size[i] & x$cna=="gain", ]
  gain[[i]] <- sum(wd$n.peaks>0, na.rm = T)/nrow(wd)
  wd <- x[x$size>size[i-1] & x$size<=size[i] & x$cna=="loss", ]
  loss[[i]] <- sum(wd$n.peaks>0, na.rm = T)/nrow(wd)
  wd <- x[x$size>size[i-1] & x$size<=size[i] & x$cna=="diploid", ]
  diploid[[i]] <- sum(wd$n.peaks>0, na.rm = T)/nrow(wd)
}

gain <- data.frame(size=as.numeric(size[-1]), proportion=as.numeric(unlist(gain)), cna="gain")
diploid <- data.frame(size=as.numeric(size[-1]), proportion=as.numeric(unlist(diploid)), cna="diploid")
loss <- data.frame(size=as.numeric(size[-1]), proportion=as.numeric(unlist(loss)), cna="loss")
output <- rbind(gain, diploid, loss)
output <- output[!is.nan(output$proportion),] 
output$cna <- as.factor(output$cna)

jpeg('tempfig.jpeg', width = (3*37.795*7), height = (3*37.795*5))
ggplot(output, aes(x=log(size), y=as.numeric(proportion), color=as.factor(cna))) + 
  geom_point(size=2, alpha=0.5) +
  ylab("Proportion of segments with >0 peaks") +
  xlab("Size of segment (log)") +
  scale_color_manual(values=c("#FF6600", "#55C667FF","#440154FF")) +
  stat_smooth(method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 0, label.y = 'top', size=8) +
  theme_custom() +
  theme(legend.text = element_text(size=20), legend.title = element_blank())
dev.off()

# relationship between size and proportion of segments with a peak is not due to clonality
clonality <- c(seq(0,1,0.1))
gain <- list()
loss <- list()
diploid <- list()
for (i in 2:length(clonality)) {
  wd <- x[x$clonalityCNA>clonality[i-1] & x$clonalityCNA<=clonality[i] & x$cna=="gain", ]
  gain[[i]] <- sum(wd$n.peaks>0, na.rm = T)/nrow(wd)
  wd <- x[x$clonalityCNA>clonality[i-1] & x$clonalityCNA<=clonality[i] & x$cna=="loss", ]
  loss[[i]] <- sum(wd$n.peaks>0, na.rm = T)/nrow(wd)
  wd <- x[x$clonalityCNA>clonality[i-1] & x$clonalityCNA<=clonality[i] & x$cna=="diploid", ]
  diploid[[i]] <- sum(wd$n.peaks>0, na.rm = T)/nrow(wd)
}

gain <- data.frame(clonality=as.numeric(clonality[-1]), proportion=as.numeric(unlist(gain)), cna="gain")
diploid <- data.frame(clonality=as.numeric(clonality[-1]), proportion=as.numeric(unlist(diploid)), cna="diploid")
loss <- data.frame(clonality=as.numeric(clonality[-1]), proportion=as.numeric(unlist(loss)), cna="loss")
output <- rbind(gain, diploid, loss)
output <- output[!is.nan(output$proportion),] 
output$cna <- as.factor(output$cna)

jpeg('tempfig.jpeg', width = (3*37.795*7), height = (3*37.795*5))
ggplot(output, aes(x=(clonality), y=as.numeric(proportion), color=as.factor(cna))) + 
  geom_point(size=2, alpha=0.5) +
  ylab("Proportion of segments with >0 peaks") +
  xlab("Clonality") +
  scale_color_manual(values=c("#FF6600", "#55C667FF","#440154FF")) +
  stat_smooth(method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 0, label.y = 'top', size=8) +
  theme_custom() +
  theme(legend.text = element_text(size=20), legend.title = element_blank())
dev.off()



# PEAK COVERAGE IS HIGHER ON SMALLER SEGMENTS (STILL UNRELATED TO CLOANLITY) =====

# higher peak coverage on smaller segments
jpeg('tempfig.jpeg', width = (3*37.795*8.38), height = (3*37.795*7.87))
ggplot(x, aes(x=log(size), y=log(coverage), fill=as.factor(cna))) + 
  geom_point(size=3, alpha=0.5, shape=21) +
  scale_y_continuous(name="Coverage (n.peak/segment size)", 
                     breaks = c(-15,-12.5,-10,-7.5),
                     labels = c((round(exp(-15),7)), (round(exp(-12.5),6)), (round(exp(-10),5)),(round(exp(-7.5),4)))) +
  scale_x_continuous(name="Size of segment", 
                     breaks = c(5,10,15),
                     labels = c(signif(round(exp(5),0),1), signif(round(exp(10),0),1), signif(round(exp(15),0),1))) +
  geom_hline(yintercept = log(0.0005), linetype="dashed") +
  scale_fill_manual(values=c("#336600", "#CC0033","#39568CFF")) +
  stat_smooth(method = "glm", color="black") +
  facet_grid(rows=vars(cna)) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = "left", label.y = 1, size=8) +
  theme_custom() +
  theme(legend.text = element_text(size=20), legend.title = element_blank(),panel.spacing = unit(2, "lines"), strip.text = element_text(size=28), legend.position = "none" )
dev.off()

mod <- lm((coverage) ~ log(size), data=x)
summary(mod)

# pull the segments with abnormally high pea coverage
highCoverageSegs <- x[which(x$coverage>0.00055),]

jpeg('tempfig.jpeg', width = (3*37.795*8.38), height = (3*37.795*7.87))
ggplot(x[which(x$n.peaks>1),], aes(x=log(size), y=log(coverage), fill=as.factor(cna))) + 
  geom_point(size=3, alpha=0.5, shape=21) +
  scale_y_continuous(name="Coverage (n.peak/segment size)", 
                     breaks = c(-15,-12.5,-10,-7.5),
                     labels = c((round(exp(-15),7)), (round(exp(-12.5),6)), (round(exp(-10),5)),(round(exp(-7.5),4)))) +
  scale_x_continuous(name="Size of segment", 
                     breaks = c(5,10,15),
                     labels = c(signif(round(exp(5),0),1), signif(round(exp(10),0),1), signif(round(exp(15),0),1))) +
  geom_hline(yintercept = log(0.0005), linetype="dashed") +
  scale_fill_manual(values=c("#336600", "#CC0033","#39568CFF")) +
  stat_smooth(method = "glm", color="black") +
  facet_grid(rows=vars(cna)) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = "left", label.y = 1, size=8) +
  theme_custom() +
  theme(legend.text = element_text(size=20), legend.title = element_blank(),panel.spacing = unit(2, "lines"), strip.text = element_text(size=28), legend.position = "none" )
dev.off()



# peak coverage on gains is explained a bit by clonality
jpeg('tempfig.jpeg', width = (3*37.795*6), height = (3*37.795*5))
ggplot(x, aes(x=(clonalityCNA), y=log(coverage), color=as.factor(cna))) + 
  geom_point(size=2, alpha=0.5) +
  ylab("Coverage (n.peak/segment size, log)") +
  xlab("Clonality of segment") +
  scale_color_manual(values=c("#FF6600", "#55C667FF","#440154FF")) +
  stat_smooth(method = "lm") +
  facet_grid(rows=vars(cna)) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 0, label.y = 'top', size=7) +
  theme_custom() +
  theme(legend.text = element_text(size=20), legend.title = element_blank())
dev.off()

# the relationship between size and clonality is positive but weak
jpeg('tempfig.jpeg', width = (3*37.795*6), height = (3*37.795*5))
ggplot(x, aes(x=(clonalityCNA), y=(size), color=as.factor(cna))) + 
  geom_point(size=2, alpha=0.5) +
  ylab("Coverage (n.peak/segment size, log)") +
  xlab("Clonality of segment (log") +
  scale_color_manual(values=c("#FF6600", "#55C667FF","#440154FF")) +
  stat_smooth(method = "lm") +
  facet_grid(rows=vars(cna)) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 0, label.y = 'top', size=7) +
  theme_custom() +
  theme(legend.text = element_text(size=20), legend.title = element_blank())
dev.off()

x$cna <- factor(x$cna, levels = c("diploid","gain","loss"))
model <- lm((coverage) ~ cna + log(size) + clonalityCNA, data = x)
summary(model)

# PEAK COVERAGE ANEUPLOIDY =====

my_comparisons <- list( c("diploid", "gain"), c("diploid", "loss"),c("gain", "loss"))

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot2::ggplot(data=x, aes(y=log(as.numeric(coverage)), x=(cna), fill=cna)) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitterdodge(jitter.width = 0.2,
                                                                          jitter.height = 0,
                                                                          dodge.width = 0.6,)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_y_continuous(name="Peak coverage \n(number of peaks/segment size)", 
                     breaks = c(-15,-12.5,-10,-7.5,-5),
                     labels = c(round(exp(-15),7),round(exp(-12.5),6),round(exp(-10),5),round(exp(-7.5),4),round(exp(5),0))) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(comparisons = my_comparisons, label="p.format", size=8) +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), axis.title.x = element_blank())
dev.off()

ggplot(data=x, aes(x=log(as.numeric(coverage)), fill=cna)) + 
  geom_histogram(alpha=0.4, position="identity", bins = 100) +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_continuous(name="Peak coverage \n(number of peaks/segment size)", 
                     breaks = c(-15,-12.5,-10,-7.5,-5),
                     limits = c(-15,-5),
                     labels = c(round(exp(-15),7),round(exp(-12.5),6),round(exp(-10),5),round(exp(-7.5),4),round(exp(5),0))) +
  theme_custom()

x <- x[which(x$coverage!=0),]
1000000*median(x$coverage)
t.test(x$coverage[which(x$cna=="diploid")], x$coverage[which(x$cna=="gain")])
t.test(x$coverage[which(x$cna=="diploid")], x$coverage[which(x$cna=="loss")])
t.test(x$coverage[which(x$cna=="gain")], x$coverage[which(x$cna=="loss")])
