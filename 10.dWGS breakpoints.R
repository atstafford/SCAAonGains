# PREREQUISITIS: load section 1 data 

breakpoints_dWGS <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/breakpoints_dWGS.rds")
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
BPregionMax <- 10000

# only consider patients with SCAA data
peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig 
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")] # keep only carcinomas ('pure')
dat <- dwgs[which(substr(names(dwgs),1,4) %in% substr(colnames(scaa),1,4))] # use raw as distinguishes diff levels of amp

sample <- list()
patient <- list()
location <- list()
l <- 1
for (i in 1:length(dat)) {
  print(i)
  wd <- dat[[i]]
  wd <- wd[!is.na(wd$CNt),]
  
  for (j in 2:nrow(wd)) {
    if( (wd$CNt[j] != wd$CNt[j-1]) & #a breakpoint
        (wd$chromosome[j] == wd$chromosome[j-1]) & #on the same chromosome
        (DescTools::Overlap(c(wd$end.pos[j-1], wd$start.pos[j]), c(hg38.centromere$start[hg38.centromere$chrom==wd$chromosome[j]], hg38.centromere$stop[hg38.centromere$chrom==wd$chromosome[j]])) == 0) & #on same arm (aka dont overlap centromere)
        (as.numeric(wd$start.pos[j]) - as.numeric(wd$end.pos[j-1]) < BPregionMax) ) { #near enough to call the breakpoint with some certainty
      
      breakpoint <- (( as.numeric(wd$start.pos[j]) - as.numeric(wd$end.pos[j-1]) ) / 2 ) + as.numeric(wd$start.pos[j])
      
      sample[[l]] <- names(dat)[i]
      patient[[l]] <- substr(names(dat)[i],1,4)
      location[[l]] <- paste(wd$chromosome[j-1], breakpoint, sep = ":")
      
      l <- l + 1
    } else {
      next
    }
  }
}

breakpoints_dWGS <- data.frame(sample = unlist(sample), patient = unlist(patient))

# count number of samples per patient
breakpoints_dWGS <- breakpoints_dWGS %>%
  group_by(patient) %>%
  mutate(
    nSamples = length(unique(sample)) 
  ) 

breakpoints_dWGS$location = unlist(location)
breakpoints_dWGS$break_chr <- sub("\\:.*", "",breakpoints_dWGS$location)
breakpoints_dWGS$break_base <- round(as.numeric(sub(".*:", "", breakpoints_dWGS$location)))



# DO NUMBER OF BREAKPOINTS CORRELATE WITH PGA? --(WIP)-- =====

BP.persample <- breakpoints_dWGS %>% group_by(sample) %>% 
  summarise(count=length(break_base))
BP.persample$patient <- substr(BP.persample$sample,1,4)

BP.perpatient <- BP.persample %>% group_by(patient) %>% 
  summarise(mean=mean(count, na.rm = T))

scaaPP$meanBP <- BP.perpatient[match(scaaPP$patient, BP.perpatient$patient), 2]

# CLUSTER NEARBY BREAKPOINTS =====

# only consider patients with SCAA data
# cluster size chosen based on what produces <2000 clusters to analyse
chrs <- unique(breakpoints_dWGS$break_chr)

clusteredBreakpoints <- list()
for (i in 1:length(chrs)) {
  print(paste(i, '/', length(chrs)))
  wd <- breakpoints_dWGS[which(breakpoints_dWGS$break_chr==chrs[i]), ]
  dat <- data.frame(base=unique(wd$break_base))
  
  # cluster bases
  dist <- dist(dat$base, method = "euclidean")
  hc <- hclust(dist, method = 'complete')
  x <- cophenetic(hc)
  cor(dist, x)
  
  # choose height that keeps within cluster <300000
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

# CLASSIFY CLONAL BREAKPOINTS =====
# count number of cluster/patient
breakpoints_dWGS <- breakpoints_dWGS %>%
  group_by(patient, cluster) %>%
  mutate(
    clusterClonality = length(unique(sample))
  ) 

# divide by number of samples per patient
breakpoints_dWGS$clusterClonality <- breakpoints_dWGS$clusterClonality / breakpoints_dWGS$nSamples

# SAVE =====
saveRDS(breakpoints_dWGS, "~/Documents/SCAA/Data/EPICC/SCAA/breakpoints_dWGS_highres10000.rds")
saveRDS(dat, "~/Documents/SCAA/Data/EPICC/SCAA/dat.rds")



