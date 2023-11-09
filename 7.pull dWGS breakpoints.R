# PREREQUISITIS: load section 1 data 

breakpoints_dWGS <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/breakpoints_dWGS.rds")
breakpoints_dWGS_highres10000 <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/breakpoints_dWGS_highres10000.rds")
randomPoints <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/randomPoints.rds")
library(GWASTools)

# PULL BREAKPOINTS ===== 

# find centromeres
hg38.centromere <- read.csv("~/Documents/SCAA/Data/UCSC_search/hg38.centromere")
library(dplyr)
start <- hg38.centromere %>% 
  group_by(chrom) %>% 
  summarise(start = min(chromStart))
end <- hg38.centromere %>% 
  group_by(chrom) %>% 
  summarise(stop = max(chromEnd))
hg38.centromere <- cbind(start, end[2])

# the breakpoint region may be large which reduces ability to know where 'real' breakpoint is
#BPregionMax <- 10000

# only consider patients with SCAA data
peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig 
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")] # keep only carcinomas ('pure')
dat <- dwgs[which(substr(names(dwgs),1,4) %in% substr(colnames(scaa),1,4))] # use raw as distinguishes diff levels of amp

sample <- list()
patient <- list()
location <- list()
chr <- list()
breakpoint <- list()
breakpointRegion_size <- list()
border <- list()
l <- 1
for (i in 1:length(dat)) {
  print(paste(i, "/", length(dat)))
  wd <- dat[[i]]
  wd <- wd[!is.na(wd$CNt),]
  
  for (j in 2:nrow(wd)) {
    if( (wd$CNt[j] != wd$CNt[j-1]) & #a breakpoint
        (wd$chromosome[j] == wd$chromosome[j-1]) & #on the same chromosome
        (DescTools::Overlap(c(wd$end.pos[j-1], wd$start.pos[j]), c(hg38.centromere$start[hg38.centromere$chrom==wd$chromosome[j]], hg38.centromere$stop[hg38.centromere$chrom==wd$chromosome[j]])) == 0)  #on same arm (aka dont overlap centromere)
         ) { 
      
      sample[[l]] <- names(dat)[i]
      patient[[l]] <- substr(names(dat)[i],1,4)
      chr[[l]] <- wd$chromosome[j]
      breakpoint[[l]] <- floor(( as.numeric(wd$start.pos[j]) - as.numeric(wd$end.pos[j-1]) ) / 2 ) + as.numeric(wd$start.pos[j])
      breakpointRegion_size[[l]] <- (as.numeric(wd$start.pos[j]) - as.numeric(wd$end.pos[j-1]))
      location[[l]] <- paste(wd$chromosome[j-1], breakpoint[[l]], sep = ":")
      
      if ( (wd$CNt[j]<2 | wd$CNt[j-1]<2) & (wd$CNt[j]==2 | wd$CNt[j-1]==2) ) {
        border[[l]] <- "dip_loss"
      } else if ( (wd$CNt[j]>2 | wd$CNt[j-1]>2) & (wd$CNt[j]==2 | wd$CNt[j-1]==2) ) {
        border[[l]] <- "dip_gain"
      } else if ( (wd$CNt[j]>2 | wd$CNt[j-1]>2) & (wd$CNt[j]<2 | wd$CNt[j-1]<2) ) {
        border[[l]] <- "gain_loss"
      } else if ( (wd$CNt[j]>2 | wd$CNt[j-1]>2) & (wd$CNt[j]>2 | wd$CNt[j-1]>2) ) {
        border[[l]] <- "gain_gain"
      } else if ( (wd$CNt[j]<2 | wd$CNt[j-1]<2) & (wd$CNt[j]<2 | wd$CNt[j-1]<2) ) {
        border[[l]] <- "loss_loss"
      }
      
      
      l <- l + 1
    } else {
      next
    }
  }
}

breakpoints_dWGS <- data.frame(sample = unlist(sample), patient = unlist(patient), location=unlist(location), 
                               chr=unlist(chr), breakpoint=unlist(breakpoint), breakpointRegion_size=unlist(breakpointRegion_size),
                               border=unlist(border))

# count number of samples per patient to work out clonality of breakpoint
breakpoints_dWGS <- breakpoints_dWGS %>%
  group_by(patient) %>%
  mutate(
    nSamples = length(unique(sample)) 
  ) 
breakpoints_dWGS <- breakpoints_dWGS %>%
  group_by(patient, location) %>%
  mutate(
    nBP = length(unique(sample)) 
  ) 
breakpoints_dWGS$BPclonality <- breakpoints_dWGS$nBP/breakpoints_dWGS$nSamples

# STOP HERE FOR NEW ANALYSIS =====
# CLUSTER NEARBY BREAKPOINTS =====

# only consider patients with SCAA data
# cluster size chosen based on what produces <2000 clusters to analyse
chrs <- unique(breakpoints_dWGS$chr)

clusteredBreakpoints <- list()
for (i in 1:length(chrs)) {
  print(paste(i, '/', length(chrs)))
  wd <- breakpoints_dWGS[which(breakpoints_dWGS$chr==chrs[i]), ]
  dat <- data.frame(base=unique(wd$break_base))
  
  # cluster bases
  dist <- dist(dat$base, method = "euclidean")
  hc <- hclust(dist, method = 'complete')
  x <- cophenetic(hc)
  cor(dist, x)
  
  # choose height that keeps within cluster <250000
  for (j in 1:length(hc$height)) {
    h <- hc$height[j]
    cut <- cutree(hc, h=h)
    dat$cluster <- unlist(cut)
    
    dat <- dat %>%
      group_by(cluster) %>%
      mutate(
        accept = ifelse( (max(base, na.rm = T) -  min(base, na.rm = T) )<=250000, TRUE, FALSE )
      ) 
    
    if (sum(dat$accept) != nrow(dat)) { # there is a FALSE, so revert
      h <- hc$height[j-1]
      cut <- cutree(hc, h=h)
      dat$cluster <- unlist(cut)
      
      dat <- dat %>%
        group_by(cluster) %>%
        mutate(
          accept = ifelse( (max(base, na.rm = T) -  min(base, na.rm = T) )<=250000, TRUE, FALSE )
        ) 
      
      clusteredBreakpoints[[i]] <- dat
      break
    }
  }
}

names(clusteredBreakpoints) <- chrs

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

# ASSIGN BREAKPOINTS TO A CLUSTER =====
# assign
output <- list()
for (i in 1:nrow(breakpoints_dWGS)) {
  chr <- breakpoints_dWGS$break_chr[i]
  wd <- clusteredBreakpoints[[which(names(clusteredBreakpoints)==chr)]]
  
  output[[i]] <- wd$cluster[which(wd$base ==  breakpoints_dWGS$break_base[i])]
}

breakpoints_dWGS$cluster <- unlist(output)

#find mid point per cluster
breakpoints_dWGS <- breakpoints_dWGS %>%
  group_by(cluster) %>%
  mutate(
    centrePoint = round(mean(break_base, na.rm = T))
  ) 


# CLASSIFY RECURRING BREAKPOINTS ------------------------------------------------------------------------------------------------------------------------------
breakpoints_dWGS <- breakpoints_dWGS %>%
  group_by(cluster) %>%
  mutate(
    clusterFrequency = length(unique(patient)) / length(unique(breakpoints_dWGS$patient))
  ) 

breakpoints_dWGS <- breakpoints_dWGS %>%
  group_by(cluster) %>%
  mutate(
    clusterFrequency2 = length(unique(patient)) 
  ) 

# CLASSIFY CLONAL BREAKPOINTS =====
# count number of cluster/patient
breakpoints_dWGS <- breakpoints_dWGS %>%
  group_by(patient, cluster) %>%
  mutate(
    clusterClonality = length(unique(sample))
  ) 

# divide by number of samples per patient
breakpoints_dWGS$clusterClonality <- breakpoints_dWGS$clusterClonality / breakpoints_dWGS$nSamples

# GENERATE RANDOM POINTS THAT ARNT IN HG38 GAPS =====

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

# sample 10,000x per chromosome
output <- list()
for (i in 1:nrow(hg38.length)) {
  output[[i]] <- paste(i, sample(1:hg38.length$length[i], 100000), sep=":")
}
randomPoints <- unlist(output)

# remove randompoints in gap regions
library(dplyr)
gap <- read.csv("~/Documents/SCAA/Data/UCSC_search/gap.csv")
drop <- list()
for (i in 1:length(randomPoints)) {
  print(paste(i,"/",length(randomPoints)))
  chr <- paste("chr", sub("\\:.*", "", randomPoints[i]), sep="")
  point <- as.numeric(sub(".*:", "", randomPoints[i])) 
  wd <- gap[which(gap$chrom==chr),]
  for (j in 1:nrow(wd)) {
    if( between(point, wd$chromStart[j], wd$chromEnd[j])  ) {
      drop[[i]] <- TRUE
      break
    } else {
      drop[[i]] <- FALSE
    }
  }
}
sum(drop==TRUE)
sum(drop==FALSE)
randomPoints <- randomPoints[drop==FALSE]

# ONLY KEEP BREAKPPOINTS THAT ARNT IN HG38 GAPS =====
drop <- list()
for (i in 1:nrow(breakpoints_dWGS)) {
  print(paste(i,"/",nrow(breakpoints_dWGS)))
  chr <- breakpoints_dWGS$break_chr[i]
  point <- breakpoints_dWGS$break_base[i]
  wd <- gap[which(gap$chrom==chr),]
  for (j in 1:nrow(wd)) {
    if( between(point, wd$chromStart[j], wd$chromEnd[j])  ) {
      drop[[i]] <- TRUE
      break
    } else {
      drop[[i]] <- FALSE
    }
  }
}
sum(drop==TRUE)
sum(drop==FALSE)
breakpoints_dWGS <- breakpoints_dWGS[drop==FALSE,]

# ALSO SELECT NON BREAKPOINTS =====
x <- sample(randomPoints,1000000)
drop <- list()
for (i in 1:length(x)) {
  print(paste(i,"/",length(x)))
  chr <- paste("chr", sub("\\:.*", "", x[i]), sep="")
  point <- as.numeric(sub(".*:", "", x[i])) 
  wd <- breakpoints_dWGS[which(breakpoints_dWGS$break_chr==chr),]
  for (j in 1:nrow(wd)) {
    if( abs(point-wd$break_base[j])<2000000  ) {
      drop[[i]] <- TRUE
      break
    } else {
      drop[[i]] <- FALSE
    }
  }
}
sum(drop==TRUE)
sum(drop==FALSE)
randomPoints_2Mbaway <- randomPoints[drop==FALSE]

# SAVE =====
saveRDS(breakpoints_dWGS, "~/Documents/SCAA/Data/EPICC/SCAA/breakpoints_dWGS.rds")
saveRDS(randomPoints, "~/Documents/SCAA/Data/EPICC/SCAA/randomPoints.rds")





