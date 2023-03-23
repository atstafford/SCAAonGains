# all analysis on carcinoma glands only

# GAP FILL =====
# read in files
myFiles <- list.files("~/Documents/SCAA/Data/EPICC/deep/Segments", full.names=T)
dwgs <- lapply(myFiles, read.delim)

# extract sample names and correct format
sampleNames <- list.files("~/Documents/SCAA/Data/EPICC/deep/Segments", full.names=F)
sampleNames <- substr(sampleNames, 7, 16)
sampleNames <- sub("_", ".", sampleNames)
names(dwgs) <- sampleNames

# keep only gland data
type <- substr(sampleNames,9,9)
dwgs <- dwgs[type == "G"]

# subset based on adenoma/carcinoma: adenoma=F-G (and C/D in C516), normal=E, carcinoma=A-D
normal <- grep("\\.E", names(dwgs))
adenomas <- c( grep("C516.C", names(dwgs)), grep("C516.D", names(dwgs)), grep("\\.F", names(dwgs)), grep("\\.G", names(dwgs)), grep("\\.H", names(dwgs)) )
carcinomas <- c( grep("\\.A", names(dwgs)), grep("\\.B", names(dwgs)), grep("\\.C", names(dwgs)), grep("\\.D", names(dwgs)) ) 
carcinomas <- carcinomas[-which(carcinomas %in% adenomas)]
dwgs <- dwgs[carcinomas]

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

# fill in gaps
for (j in 1:length(dwgs)) {
  start <- list()
  stop <- list()
  wd <- dwgs[[j]]
  
  for (i in 1:nrow(wd)) {
    chr <- wd$chromosome[i]
    if (i == 1) { # top row
      start[[i]] <- wd$start.pos[i]
    } else if ( wd$chromosome[i] != wd$chromosome[i-1]  ) {  # first row of chr
      start[[i]] <- wd$start.pos[i]
    } else {
      diff <- ceiling( (wd$start.pos[i] - wd$end.pos[i-1]) / 2 ) - 1
      start[[i]] <- wd$start.pos[i] - diff
    }
    
    
    if (i == nrow(wd)) { # last row pf chr
      stop[[i]] <- wd$end.pos[i]
    } else if (wd$chromosome[i] != wd$chromosome[i+1] ) { 
      stop[[i]] <- wd$end.pos[i]
    } else {
      diff <- floor( (wd$start.pos[i+1] - wd$end.pos[i]) / 2 )
      stop[[i]] <- wd$end.pos[i] + diff
    }
  }
  
  dwgs[[j]]$start <- unlist(start)
  dwgs[[j]]$stop <- unlist(stop)
}

# remove duplicates
index <- !duplicated(names(dwgs))
dwgs <- dwgs[index]

# PLOIDY RECENTRE =====

# assess ploidy per sample and recentre 
sample <- names(dwgs)
dwgs.ploidyRecentre <- list()
ploidy.persample <- list()
for (i in 1:length(dwgs)) {
  weights <- dwgs[[i]]$end.pos - dwgs[[i]]$start.pos
  sumweight <- sum(weights)
  ploidy <- round(sum(dwgs[[i]]$CNt*weights, na.rm = T)/sumweight, 1)
  ploidy.persample[[i]] <- ploidy
  
  if ( ploidy > 2.5 ) {
    distance <- abs(2.5 - ploidy)
    print(paste(i, ploidy))
    dwgs.ploidyRecentre[[i]] <- dwgs[[i]]
    dwgs.ploidyRecentre[[i]]$CNt <- dwgs.ploidyRecentre[[i]]$CNt - distance
  }
  else {
    dwgs.ploidyRecentre[[i]] <- dwgs[[i]]
  }
  dwgs.ploidyRecentre[[i]] <- dwgs.ploidyRecentre[[i]][which(dwgs.ploidyRecentre[[i]]$chr!="X" | dwgs.ploidyRecentre[[i]]$chr!="Y"),]
}
names(dwgs.ploidyRecentre) <- names(dwgs)

# only carcinomas are ploidy
ploidy.persample <- data.frame(patient=substr(names(dwgs),1,4), sample=names(dwgs), ploidy=unlist(ploidy.persample))

# PLOIDY PER PATIENT (CARCINOMA ONLY) ----------------------------------------------------------------
library(dplyr)

ploidy.perpatient <- ploidy.persample %>% group_by(patient) %>% 
  summarise(mean_ploidy=mean(ploidy),
            SD_ploidy=sd(ploidy),
            max_ploidy=max(ploidy),
            min_ploidy=min(ploidy))

# PATIENT PER DATAFRAME (CARCINOMA ONLY)----------------------------------------------------------------

# non recentered
sample <- names(dwgs)
patient <- sub("\\..*", "", sample)

dwgs.perpatient <- list()

for (i in 1:length(unique(patient))) {
  index <- grep(unique(patient)[i], sample)
  list <- list()
  
  for (j in 1:length(index)) {
    list[[j]] <- data.frame(chromosome = substr(dwgs[[index[j]]]$chromosome, 4, nchar(dwgs[[index[j]]]$chromosome)),
                            start = dwgs[[index[j]]]$start,
                            stop = dwgs[[index[j]]]$stop,
                            cna = dwgs[[index[j]]]$CNt )
    colnames(list[[j]])[4] <- sample[index[j]]
  }
  
  dwgs.perpatient[[i]] <- powerjoin::power_full_join(list, by = c("chromosome","start","stop"))
}
names(dwgs.perpatient) <- unique(patient)

# recentered
sample <- names(dwgs.ploidyRecentre)
patient <- sub("\\..*", "", sample)

dwgs.ploidyRecentre.perpatient <- list()

for (i in 1:length(unique(patient))) {
  index <- grep(unique(patient)[i], sample)
  list <- list()
  
  for (j in 1:length(index)) {
    list[[j]] <- data.frame(chromosome = substr(dwgs.ploidyRecentre[[index[j]]]$chromosome, 4, nchar(dwgs.ploidyRecentre[[index[j]]]$chromosome)),
                            start = dwgs.ploidyRecentre[[index[j]]]$start,
                            stop = dwgs.ploidyRecentre[[index[j]]]$stop,
                            cna = dwgs.ploidyRecentre[[index[j]]]$CNt )
    colnames(list[[j]])[4] <- sample[index[j]]
  }
  
  dwgs.ploidyRecentre.perpatient[[i]] <- powerjoin::power_full_join(list, by = c("chromosome","start","stop"))
}
names(dwgs.ploidyRecentre.perpatient) <- unique(patient)


# SAVE -----------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(dwgs, "~/Documents/SCAA/Data/dwgs.rds")
saveRDS(dwgs.ploidyRecentre, "~/Documents/SCAA/Data/dwgs.ploidyRecentre.rds")
saveRDS(ploidy.persample, "~/Documents/SCAA/Data/ploidy.persample.rds")
saveRDS(ploidy.perpatient, "~/Documents/SCAA/Data/ploidy.perpatient.rds")
saveRDS(dwgs.perpatient, "~/Documents/SCAA/Data/dwgs.perpatient.rds")
saveRDS(dwgs.ploidyRecentre.perpatient, "~/Documents/SCAA/Data/dwgs.ploidyRecentre.perpatient.rds")
