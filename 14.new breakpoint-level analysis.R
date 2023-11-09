# GET BREAKPOINT DATA =====
gap <- read.csv("~/Documents/SCAA/Data/UCSC_search/gap.csv")

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

peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")] # keep only carcinomas ('pure')
peak.location <- data.frame(chr=sub("\\:.*", "", rownames(scaa)), start=sub(".*[:]([^.]+)[-].*", "\\1", rownames(scaa)), stop=sub(".*-", "", rownames(scaa)))

# unique breakpoints
breakpoint.data <- breakpoints_dWGS[,-which(colnames(breakpoints_dWGS) %in% c("sample"))]
breakpoint.data <- breakpoint.data[,-which(colnames(breakpoint.data) %in% c("border"))]
breakpoint.data <- breakpoint.data[!duplicated(breakpoint.data),]
breakpoint.data <- breakpoint.data[which(breakpoint.data$chr!="chrX" & breakpoint.data$chr!="chrY"),]

# only ones on a dip/loss border 
#breakpoint.data <- breakpoint.data[which(breakpoint.data$border=="dip_loss"),]
#breakpoint.data <- breakpoint.data[,-which(colnames(breakpoint.data) %in% c("border"))]

# only breakpoints with high enough clonality (>0.5), and small enough break region
breakpoint.data <- breakpoint.data[which(breakpoint.data$BPclonality>0.9),]
breakpoint.data <- breakpoint.data[which(breakpoint.data$breakpointRegion_size< 20000),]
hist(breakpoint.data$breakpointRegion_size)

# add arm and distance to centromere/telomere
breakpoint.data$centro_start <- hg38.centromere[match(breakpoint.data$chr, hg38.centromere$chrom), 2]
breakpoint.data$centro_stop <- hg38.centromere[match(breakpoint.data$chr, hg38.centromere$chrom), 3]
arm <-list()
centro_dist <- list()
for (i in 1:nrow(breakpoint.data)) {
  if (breakpoint.data$breakpoint[i] < breakpoint.data$centro_start[i]) {
    arm[[i]] <- "p"
    centro_dist[[i]] <- breakpoint.data$centro_start[i] - breakpoint.data$breakpoint[i]
  } else if (breakpoint.data$breakpoint[i] > breakpoint.data$centro_stop[i]) {
    arm[[i]] <- "q"
    centro_dist[[i]] <- breakpoint.data$breakpoint[i] - breakpoint.data$centro_stop[i]
  } else {
    arm[[i]] <- "centro"
    centro_dist[[i]] <- min(abs(breakpoint.data$centro_start[i] - breakpoint.data$breakpoint[i]),
                            abs(breakpoint.data$breakpoint[i] - breakpoint.data$centro_stop))
  }
}
breakpoint.data$arm <- unlist(arm)
breakpoint.data$centro_dist <- unlist(centro_dist)

# calculate distance to nearest BP within same patient
nearestBP <- list()
for (i in 1:nrow(breakpoint.data)) {
  print(paste(i,"/",nrow(breakpoint.data)))
  wd <- breakpoint.data[which(breakpoint.data$patient==breakpoint.data$patient[i] & breakpoint.data$chr==breakpoint.data$chr[i] & breakpoint.data$breakpoint!=breakpoint.data$breakpoint[i]),]
  nearestBP[[i]] <- min(abs(wd$breakpoint - breakpoint.data$breakpoint[i]))
}
breakpoint.data$nearestBP <- unlist(nearestBP)
breakpoint.data$nearestBP[!is.finite(breakpoint.data$nearestBP)] <- NA
hist(breakpoint.data$nearestBP, breaks = 100)

# add per patient scaa prop
breakpoint.data$propSCAApp <- scaaPP[match(breakpoint.data$patient, scaaPP$patient), 2]

# CLUSTER NEARBY BREAKPOINTS =====

# cluster per patient
clusteredBreakpoints <- list()
l <- 1
for (i in 1:length(unique(breakpoint.data$patient))) {
  print(paste(i, '/', length(unique(breakpoint.data$patient))))
  patient <- unique(breakpoint.data$patient)[i]
  wd <- breakpoint.data[which(breakpoint.data$patient==patient),]
  
  for (k in 1:length(unique(wd$chr))) {
    chr <- unique(wd$chr)[k]
    wd2 <- wd[which(wd$chr==chr), ]
    dat <- data.frame(breakpoint=unique(wd2$breakpoint))
    
    if (nrow(dat) == 1) { # only one break on chromosome
      dat$cluster <- c(1)
      dat$accept <- TRUE
      dat$patient <- patient
      dat$chr <- chr
      clusteredBreakpoints[[l]] <- dat
      l <- l + 1
      next
    } else if (nrow(dat) == 2 & abs(dat$breakpoint[2] - dat$breakpoint[1])>1000000) { # only two BP and the base distance i over 1mb
      dat$cluster <- c(1:2)
      dat$accept <- TRUE
      dat$patient <- patient
      dat$chr <- chr
      clusteredBreakpoints[[l]] <- dat
      l <- l + 1
      next
    }
    
    # cluster bases
    dist <- dist(dat$breakpoint, method = "euclidean")
    hc <- hclust(dist, method = 'complete')
    x <- cophenetic(hc)
    cor(dist, x)
    
    # choose height that keeps within cluster <1000000
    for (j in length(hc$height):1) {
      h <- hc$height[j]
      cut <- cutree(hc, h=h)
      dat$cluster <- unlist(cut)
      
      dat <- dat %>%
        group_by(cluster) %>%
        mutate(
          accept = ifelse( (max(breakpoint, na.rm = T) -  min(breakpoint, na.rm = T) )<=1000000, TRUE, FALSE )
        ) 
      
      if (sum(dat$accept) == nrow(dat)) { # Are they all TRUE? then h is correct. if not, 1 or more of the clusters is too large
        dat$patient <- patient
        dat$chr <- chr
        clusteredBreakpoints[[l]] <- dat
        l <- l + 1
        break
      } else if ( (sum(dat$accept) != nrow(dat)) & j==1) { # even the lowest height is yielding too big a cluster
        problemCluster_length <- length(dat$cluster[dat$accept==FALSE])
        if (problemCluster_length > 2) {
          print("fuck")
        }
        problemCluster <- unique(dat$cluster[dat$accept==FALSE])
        dat$cluster[dat$cluster==problemCluster][2] <- max(dat$cluster)+1
        dat$patient <- patient
        dat$chr <- chr
        clusteredBreakpoints[[l]] <- dat
        l <- l + 1
      } else { # try a lower height
        next
      }
    }
  }
}

# make clusters unique across chromosomes
for (i in 1:length(clusteredBreakpoints)) {
  if (i==1) {
    max <- max(clusteredBreakpoints[[i]]$cluster)
    next
  } else {
    clusteredBreakpoints[[i]]$cluster <- clusteredBreakpoints[[i]]$cluster + max
    max <- max(clusteredBreakpoints[[i]]$cluster)
  }
}

x <- do.call(rbind, clusteredBreakpoints)
breakpoint.data <- merge(x=breakpoint.data, y=x[-3], by=c("patient","breakpoint","chr"), all.x=TRUE)
sum(is.na(breakpoint.data$cluster))



for (i in 1:length(unique((breakpoint.data$cluster)))) {
  wd <- breakpoint.data[which(breakpoint.data$cluster==unique(breakpoint.data$cluster)[i]),]
  if (length(unique(wd$patient))>1) {
    print("patient")
  }
  if (length(unique(wd$arm))>1) {
    print("arm")
  }
  if (length(unique(wd$chr))>1) {
    print("chr")
  }
}





# GET PEAKS AND SCAAS AROUND BREAKPOINTS AND RANDOM POINTS =====

windows <- c(1000, 5000, 10000, 50000, 100000, 500000, 1000000)

# pull peaks and scaa per breakpoint and random point
distFromRandtoBP <- 2000000
l <- 1
breakpoint_peak.sum <- list()
breakpoint_scaa.sum <- list()
randompoint_peak.sum <- list()
randompoint_scaa.sum <- list()
randompoint <- list()
for (i in 1:nrow(breakpoint.data)) {
  patient <- breakpoint.data$patient[i]
  chr <- breakpoint.data$chr[i]
  arm <-  breakpoint.data$arm[i]
  breakpoint <- breakpoint.data$breakpoint[i]
  
  # randomly generate a locus on same chromosome arm, <2Mb away, not in 10000 telomeres (allow 1mb into centromere)
  centro_start <- hg38.centromere$start[hg38.centromere$chrom==chr]
  centro_stop <- hg38.centromere$stop[hg38.centromere$chrom==chr]
  chr_length <- hg38.length$length[hg38.length==substr(chr,4,nchar(chr))]

  drop <- TRUE
  while(drop==TRUE) {
    if(arm=="p") {
      rando <- sample(10000:centro_start, 1)
    }
    else if(arm=="q") {
      rando <- sample(centro_stop:chr_length-10000, 1)
    }
    
    wd <- breakpoint.data[which(breakpoint.data$patient==patient & breakpoint.data$chr==chr),]
    drop <- ifelse(min(abs(wd$breakpoint - rando)) <= distFromRandtoBP, TRUE, FALSE)
  }
  
  # get scaa and peak locations for chromosome/patient
  peak.wd <- peak.location[which(peak.location$chr==chr),]
  scaa.wd <- peak.location[scaa[,substr(colnames(scaa), 1, 4) == patient]==1,]
  scaa.wd <- scaa.wd[which(scaa.wd$chr==chr),]
  
  # pull number of peaks and scaas around breakpoint and random point
  for (j in 1:length(windows)) {
    print(paste(i,"/",nrow(breakpoint.data),":" ,j,"/",length(windows)))
    breakpoint_peak.sum[[l]] <- sum(apply(peak.wd, 1, function(x) DescTools::Overlap(as.numeric(c(x[2], x[3])), as.numeric(c(breakpoint-windows[j],breakpoint+windows[j]))))>1)
    breakpoint_scaa.sum[[l]] <- sum(apply(scaa.wd, 1, function(x) DescTools::Overlap(as.numeric(c(x[2], x[3])), as.numeric(c(breakpoint-windows[j],breakpoint+windows[j]))))>1)
    randompoint_peak.sum[[l]] <- sum(apply(peak.wd, 1, function(x) DescTools::Overlap(as.numeric(c(x[2], x[3])), as.numeric(c(rando-windows[j],rando+windows[j]))))>1)
    randompoint_scaa.sum[[l]] <- sum(apply(scaa.wd, 1, function(x) DescTools::Overlap(as.numeric(c(x[2], x[3])), as.numeric(c(rando-windows[j],rando+windows[j]))))>1)
    randompoint[[l]] <- rando
    l <- l + 1
  }
}

# SAVE AND PUT IN MATRIX =====
saveRDS(breakpoint_peak.sum, "~/Documents/SCAA/Data/breakpoint_peak.sum_diploss.rds")
saveRDS(breakpoint_scaa.sum, "~/Documents/SCAA/Data/breakpoint_scaa.sum_diploss.rds")
saveRDS(randompoint_peak.sum, "~/Documents/SCAA/Data/randompoint_peak.sum_diploss.rds")
saveRDS(randompoint_scaa.sum, "~/Documents/SCAA/Data/randompoint_scaa.sum_diploss.rds")

# put in matrices
breakpoint_peak.df <- data.frame(matrix(unlist(breakpoint_peak.sum), nrow = nrow(breakpoint.data), byrow = TRUE))
breakpoint_scaa.df <- data.frame(matrix(unlist(breakpoint_scaa.sum), nrow = nrow(breakpoint.data), byrow = TRUE))
randompoint_peak.df <- data.frame(matrix(unlist(randompoint_peak.sum), nrow = nrow(breakpoint.data), byrow = TRUE))
randompoint_scaa.df <- data.frame(matrix(unlist(randompoint_scaa.sum), nrow = nrow(breakpoint.data), byrow = TRUE))
colnames(breakpoint_peak.df) <- colnames(breakpoint_scaa.df) <- colnames(randompoint_peak.df) <- colnames(randompoint_scaa.df) <- windows

breakpoint_propSCAA.df <- breakpoint_scaa.df/breakpoint_peak.df
randompoint_propSCAA.df <- randompoint_scaa.df/randompoint_peak.df

breakpoint_propSCAA.df$location <- breakpoint.data$location
randompoint_propSCAA.df$location <- breakpoint.data$location
breakpoint_scaa.df$location <- breakpoint.data$location
randompoint_scaa.df$location <- breakpoint.data$location

saveRDS(breakpoint_peak.df, "~/Documents/SCAA/Data/breakpoint_peak.df.rds")
saveRDS(breakpoint_scaa.df, "~/Documents/SCAA/Data/breakpoint_scaa.df.rds")
saveRDS(randompoint_peak.df, "~/Documents/SCAA/Data/randompoint_peak.df.rds")
saveRDS(randompoint_scaa.df, "~/Documents/SCAA/Data/randompoint_scaa.df.rds")

breakpoint_peak.df <- readRDS( "~/Documents/SCAA/Data/breakpoint_peak.df.rds")
breakpoint_scaa.df <- readRDS( "~/Documents/SCAA/Data/breakpoint_scaa.df.rds")
randompoint_peak.df <- readRDS( "~/Documents/SCAA/Data/randompoint_peak.df.rds")
randompoint_scaa.df <- readRDS( "~/Documents/SCAA/Data/randompoint_scaa.df.rds")

# VIOLIN OF CUMULATIVE SCAA =====

# onyl consider BP that are >2Mb away from other BP (only consider up to 1Mb distance cumulative)
# remove BP (and their random pairs) that dont have any scaas within 1Mb
lonelyBP <- breakpoint_scaa.df[which(breakpoint.data$nearestBP>2000000 & breakpoint_scaa.df$`1e+06`!=0),c(1:7)]
lonelyBP <- data.frame(t(apply(lonelyBP, 1, function(c) ecdf(c)(c))))
colnames(lonelyBP) <- windows[c(1:7)]
lonelyBP$location <- breakpoint.data$location[which(breakpoint.data$nearestBP>2000000 & breakpoint_scaa.df$`1e+06`!=0)]
lonelyBP <- melt(lonelyBP, na.rm = TRUE, value.name = 'value', id="location")
colnames(lonelyBP) <- c("location","window","value")
lonelyBP$window <- as.numeric(as.character(lonelyBP$window))
lonelyBP$group <- "break"

lonelyRP <- randompoint_scaa.df[which(breakpoint.data$nearestBP>2000000 & breakpoint_scaa.df$`1e+06`!=0),c(1:7)]
lonelyRP <- data.frame(t(apply(lonelyRP, 1, function(c) ecdf(c)(c))))
colnames(lonelyRP) <- windows[c(1:7)]
lonelyRP$location <- breakpoint.data$location[which(breakpoint.data$nearestBP>2000000 & breakpoint_scaa.df$`1e+06`!=0)]
lonelyRP <- melt(lonelyRP, na.rm = TRUE, value.name = 'value', id="location")
colnames(lonelyRP) <- c("location","window","value")
lonelyRP$window <- as.numeric(as.character(lonelyRP$window))
lonelyRP$group <- "random"

x <- rbind(lonelyRP, lonelyBP)

jpeg('tempfig.jpeg', width = (3*37.795*12), height = (3*37.795*5))
ggplot(x, aes(x=factor(window, levels = windows), y=(value), fill=as.factor(group))) + 
  #geom_violin(trim = T, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=4, alpha = .4, shape=21, position=position_jitterdodge(0.2,0.03)) +
  geom_boxplot(width = .5, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_fill_manual(values = c("#2c728eff","#20a486ff")) +
  ylab("Cumulative frequency") +
  theme_custom() +
  stat_compare_means(aes(group = group), label="p.format", label.y = 1.1,  size=7, method = "wilcox.test") +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
dev.off()

ggplot(x[which(substr(sub("\\:.*", "", x$location),4,nchar(sub("\\:.*", "", z$location)))==1),], aes(x=as.factor(window), y=(value), color=group)) + 
  geom_line(aes(group=interaction(group, location)), alpha=1) +
  geom_hline(yintercept = 0) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=5))

# DO CLUSTERS OF BREAKPOINTS HAVE A HIGHER SCAA =====
clusters <- data.frame(table(breakpoint.data$cluster))
clusters <- clusters[which(clusters$Freq>1),]


friendlyBP <- breakpoint.data[which(breakpoint.data$cluster %in% clusters$Var1),]
friendlyBP$clusterFreq <- clusters[match(friendlyBP$cluster, clusters$Var1), 2]
friendlyBP <- friendlyBP %>%
  group_by(cluster) %>%
  mutate(
    clusterSize = max(breakpoint, na.rm = T) -  min(breakpoint, na.rm = T) 
  ) 

hist(clusters$Freq) #most in pairs, some v crowded
hist(unique(friendlyBP$clusterSize[which(friendlyBP$clusterSize<200000)])) # most are small patches

# grab the busy patches (may be chromothripsis)
superFriendlyBP <- friendlyBP[which(friendlyBP$clusterFreq>3),]
x <- data.frame(table(superFriendlyBP$patient, superFriendlyBP$chr, superFriendlyBP$cluster))

patient <- list()
chr <- list()
arm <- list()
clusterBreak_start <- list()
clusterBreak_stop <- list()
clusterSize <- list()
clusterFreq <- list()
clusterBPdensity <- list()
for (i in 1:length(unique(superFriendlyBP$cluster))) { # for each super-friendly patch
  wd <- superFriendlyBP[which(superFriendlyBP$cluster == unique(superFriendlyBP$cluster)[i]),]
  patient[[i]] <- unique(wd$patient)
  chr[[i]] <- unique(wd$chr)
  arm[[i]] <-  unique(wd$arm)
  clusterBreak_start[[i]] <- min(wd$breakpoint)
  clusterBreak_stop[[i]] <- max(wd$breakpoint)
  clusterSize[[i]] <- unique(wd$clusterSize)
  clusterFreq[[i]] <- unique(wd$clusterFreq)
  clusterBPdensity[[i]] <- unique(wd$clusterFreq) / unique(wd$clusterSize)
}

superFriendlyClusters <- data.frame(patient=unlist(patient), chr=unlist(chr), arm=unlist(arm),
                                    clusterBreak_start=unlist(clusterBreak_start), clusterBreak_stop=unlist(clusterBreak_stop),
                                    clusterSize=unlist(clusterSize), clusterFreq=unlist(clusterFreq), clusterBPdensity=unlist(clusterBPdensity))
  

clusterBreak_peak.sum <- list()
clusterBreak_scaa.sum <- list()
clusterRandom_peak.sum <- list()
clusterRandom_scaa.sum <- list()
randomCluster_start <- list()
randomCluster_stop <- list()
l <- 1
for (i in 1:nrow(superFriendlyClusters)) {
  print(paste(i,"/",nrow(superFriendlyClusters)))
  patient <- superFriendlyClusters$patient[i]
  chr <- superFriendlyClusters$chr[i]
  arm <-  superFriendlyClusters$arm[i]
  clusterBreak_start <- superFriendlyClusters$clusterBreak_start[i]
  clusterBreak_stop <- superFriendlyClusters$clusterBreak_stop[i]
  #size <- superFriendlyClusters$clusterSize[i]
  size <- 1000000
  
  # generate a random patch of same size with no breakpoints on same arm/patient (no overlap with centro/telo)
  centro_start <- hg38.centromere$start[hg38.centromere$chrom==chr]
  centro_stop <- hg38.centromere$stop[hg38.centromere$chrom==chr]
  chr_length <- hg38.length$length[hg38.length==substr(chr,4,nchar(chr))]
  
  drop <- TRUE
  while(drop==TRUE) {
    if(arm=="p") {
      rando_start <- sample(10000:centro_start-size, 1)
    }
    else if(arm=="q") {
      rando_start <- sample(centro_stop:chr_length-(10000+size), 1)
    }
    
    rando_stop <- rando_start+size
    
    wd <- superFriendlyClusters[which(superFriendlyClusters$patient==patient & superFriendlyClusters$chr==chr),]
    drop <- ifelse(sum(apply(wd, 1, function(x) DescTools::Overlap(as.numeric(c(x[4], x[5])), as.numeric(c(rando_start, rando_stop)))))>1, TRUE, FALSE)
  }
  
  # get scaa and peak locations for chromosome/patient
  peak.wd <- peak.location[which(peak.location$chr==chr),]
  scaa.wd <- peak.location[scaa[,substr(colnames(scaa), 1, 4) == patient]==1,]
  scaa.wd <- scaa.wd[which(scaa.wd$chr==chr),]
  
  # pull number of peaks and scaas in the breakpoint and random patch
  clusterBreak_peak.sum[[i]] <- sum(apply(peak.wd, 1, function(x) DescTools::Overlap(as.numeric(c(x[2], x[3])), as.numeric(c(clusterBreak_start-500000, clusterBreak_stop+500000))))>1)
  clusterBreak_scaa.sum[[i]] <- sum(apply(scaa.wd, 1, function(x) DescTools::Overlap(as.numeric(c(x[2], x[3])), as.numeric(c(clusterBreak_start-500000, clusterBreak_stop+500000))))>1)
  clusterRandom_peak.sum[[i]] <- sum(apply(peak.wd, 1, function(x) DescTools::Overlap(as.numeric(c(x[2], x[3])), as.numeric(c(rando_start-500000, rando_stop+500000))))>1)
  clusterRandom_scaa.sum[[i]] <- sum(apply(scaa.wd, 1, function(x) DescTools::Overlap(as.numeric(c(x[2], x[3])), as.numeric(c(rando_start-500000, rando_stop+500000))))>1)
  randomCluster_start[[i]] <- rando_start
  randomCluster_stop[[i]] <- rando_stop
}
  
superFriendlyClusters$n.peak <- unlist(clusterBreak_peak.sum)
superFriendlyClusters$n.scaa <- unlist(clusterBreak_scaa.sum)
superFriendlyClusters$random_start <- unlist(randomCluster_start)
superFriendlyClusters$random_stop <- unlist(randomCluster_stop)
superFriendlyClusters$random_n.peak <- unlist(clusterRandom_peak.sum)
superFriendlyClusters$random_n.scaa <- unlist(clusterRandom_scaa.sum)
  
x <- data.frame(prop=superFriendlyClusters$n.scaa/superFriendlyClusters$n.peak, group="BP")
y <- data.frame(prop=superFriendlyClusters$random_n.scaa/superFriendlyClusters$random_n.peak, group="RP")
x <- rbind(x,y)
ggplot(data=x, aes(x=group, y=(prop))) +
  geom_violin() +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0) +
  stat_compare_means(aes(group = group), label="p.format",  size=7, method = "wilcox.test") +
  theme_custom()

  # get SCAA prop within patch
  # get SCAA prop within patch plus 50kb
  # get SCAA prop within patch plus 1Mb
  # from patch centre, pull SCAA number at windows


  





# look at location on loney- near whole arm charges vs focal. near centro? near telo?
# plot dist betwene clonal enough BP, grab v close BP, near centro/telo?

# 
# cluster brealpoints (max size 1mb)
# record size of cluster and number of scaas in breakpoint cluster
# find patch of same size with no breakpoint and pull scaa number
# also only consider BP patch with no BP for same window either side and compar scaa against agjacent patch
# cum freq moving away from centre of patch
# near=scaa induced BP, farther=cna induce BP
# ecdf on all BP window 26bp










# ======
# VIOLIN OF PROP SCAA =====
bp <- breakpoint_propSCAA.df[which(breakpoint.data$nearestBP>2000000 & breakpoint_scaa.df$`1e+06`!=0),]
bp <- melt(bp, na.rm = FALSE, value.name = 'value', id="location")
colnames(bp) <- c("location","window","SCAA_burden")
bp$group <- "break"

rp <- randompoint_propSCAA.df[which(breakpoint.data$nearestBP>2000000 & breakpoint_scaa.df$`1e+06`!=0),]
rp <- melt(rp, na.rm = FALSE, value.name = 'value', id="location")
colnames(rp) <- c("location","window","SCAA_burden")
rp$group <- "random"

x <- rbind(bp,rp)

ggplot(x, aes(x=as.factor(window), y=log(SCAA_burden), fill=as.factor(group))) + 
  #geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=2, alpha = .4, shape=21, position=position_jitterdodge()) +
  geom_boxplot(width = .5, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("#2c728eff","#20a486ff")) +
  stat_compare_means(aes(group = group), label="p.format", label.y = 3,  size=8, method = "wilcox.test") +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=5))

# of a single location
ggplot(x[which(x$location=="chr5:85137418"),], aes(x=as.numeric(window), y=(value), fill=as.factor(group))) + 
  geom_point(shape=21) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=5))
