# PREREQUISITIS: load section 1 data 


# PULL CNA DATA PER PEAK (diploid may be LOH) -------------------------------------------------------

# pull regions with a scaa
peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")]

# for non recentered
# list with element per scaa detailing any cnas
cnaPerScaa <- list
for (i in 1:nrow(scaa)) {
  print(paste(i, "/",nrow(scaa)))
  peak <- rownames(scaa)[i]
  
  # pull cna data for the patient
  list <- list
  for (j in 1:ncol(scaa)) {
    patient <- substr(colnames(scaa)[j],1,4)
    scaaTRUE <- scaa[i,j]
    dwgs.wd <- dwgs.perpatient[[patient]]
    
    # pull the cnas overlapping the peak
    dwgs.wd <- dwgs.wd[dwgs.wd$chr==substr(sub("\\:.*", "", peak),4,nchar(sub("\\:.*", "", peak))), ]
    for (k in 1:nrow(dwgs.wd)) {
      start <- as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", peak))
      stop <- as.numeric(sub(".*-", "", peak))
      if ( DescTools::Overlap(c(start,stop), c(dwgs.wd$start[k], dwgs.wd$stop[k])) >0 ) {
        cnaStart <- dwgs.wd$start[k]
        cnaStop <- dwgs.wd$stop[k]
        segSize <- dwgs.wd$stop[k]- dwgs.wd$start[k]
        n <- ncol(dwgs.wd)-3
        nDip <- sum(dwgs.wd[k, -c(1:3)]==2)
        nLoss <- sum(dwgs.wd[k, -c(1:3)]<2)
        nGain <- sum(dwgs.wd[k, -c(1:3)]>2)
        fracDip <- round(nDip/n, 2)
        fracLoss <- round(nLoss/n, 2)
        fracGain <- round(nGain/n, 2)
        distancetoBP <- min((start-cnaStart), (cnaStop-stop))
        list[[j]] <- setNames(c(patient, peak, scaaTRUE,cnaStart, cnaStop, segSize, nDip, nLoss, nGain, fracDip, fracLoss, fracGain, distancetoBP), 
                              c("patient", "peak", "scaaTRUE","cnaStart", "cnaStop", "segSize", "nDip", "nLoss", "nGain", "fracDip", "fracLoss", "fracGain", "distancetoBP"))
        break
      } else {
        next
      }
    }
  }
  cnaPerScaa[[i]] <- do.call(rbind, list)
}
scaaCNA <- do.call(rbind, cnaPerScaa)
scaaCNA <- data.frame(scaaCNA)
scaaCNA <- na.omit(scaaCNA)
scaaCNA[ , c(4:13)] <- apply(scaaCNA[ , c(4:13)], 2, function(x) as.numeric(as.character(x)))


# for recentered
# list with element per scaa detailing any cnas
cnaPerScaa.ploidyRecenter <- list
for (i in 1:nrow(scaa)) {
  print(i)
  peak <- rownames(scaa)[i]
  
  # pull cna data for the patient
  list <- list
  for (j in 1:ncol(scaa)) {
    patient <- substr(colnames(scaa)[j],1,4)
    scaaTRUE <- scaa[i,j]
    dwgs.wd <- dwgs.ploidyRecentre.perpatient[[patient]]
    
    # pull the cnas overlapping the peak
    dwgs.wd <- dwgs.wd[dwgs.wd$chr==substr(sub("\\:.*", "", peak),4,nchar(sub("\\:.*", "", peak))), ]
    for (k in 1:nrow(dwgs.wd)) {
      start <- as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", peak))
      stop <- as.numeric(sub(".*-", "", peak))
      if ( DescTools::Overlap(c(start,stop), c(dwgs.wd$start[k], dwgs.wd$stop[k])) >0 ) {
        cnaStart <- dwgs.wd$start[k]
        cnaStop <- dwgs.wd$stop[k]
        segSize <- dwgs.wd$stop[k]- dwgs.wd$start[k]
        n <- ncol(dwgs.wd)-3
        nDip <- sum(dwgs.wd[k, -c(1:3)]==2)
        nLoss <- sum(dwgs.wd[k, -c(1:3)]<2)
        nGain <- sum(dwgs.wd[k, -c(1:3)]>2)
        fracDip <- round(nDip/n, 2)
        fracLoss <- round(nLoss/n, 2)
        fracGain <- round(nGain/n, 2)
        distancetoBP <- min((start-cnaStart), (cnaStop-stop))
        list[[j]] <- setNames(c(patient, peak, scaaTRUE,cnaStart, cnaStop, segSize, nDip, nLoss, nGain, fracDip, fracLoss, fracGain, distancetoBP), 
                              c("patient", "peak", "scaaTRUE","cnaStart", "cnaStop", "segSize", "nDip", "nLoss", "nGain", "fracDip", "fracLoss", "fracGain", "distancetoBP"))
        break
      } else {
        next
      }
    }
  }
  cnaPerScaa.ploidyRecenter[[i]] <- do.call(rbind, list)
}

scaaCNA.ploidyRecenter <- do.call(rbind, cnaPerScaa.ploidyRecenter)
scaaCNA.ploidyRecenter <- data.frame(scaaCNA.ploidyRecenter)
scaaCNA.ploidyRecenter <- na.omit(scaaCNA.ploidyRecenter)
scaaCNA.ploidyRecenter[ , c(4:13)] <- apply(scaaCNA.ploidyRecenter[ , c(4:13)], 2, function(x) as.numeric(as.character(x)))



#cnaPerScaa <- readRDS("~/Documents/SCAA/Data/cnaPerScaa.rds")
cnaPerScaa.ploidyRecenter <- cnaPerScaa 
scaaCNA.ploidyRecenter <- scaaCNA 

# DEFINE FOCAL CNAS WITH PEAKS -------------------------------------------------------

# long segments will appear to be enriched because of chance
# you can compare the relative enrichment for SCAAs (as compared to diploid) between long vs short gains
# need to decide what level of clonality to consider a region a gain or diploid

# label CNAs if focal or arm length
library(data.table)
x <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", 
           col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
x <- x[ , .(length = sum(chromEnd - chromStart)), 
        by = .(chrom, arm = substring(name, 1, 1)) ]

# for non recentered
pcArm <- list()
arm <- list()
for (i in 1:nrow(scaaCNA)) {
  print(i)
  chr <- sub("\\:.*", "", scaaCNA$peak[i])
  start <- as.numeric(scaaCNA$cnaStart[i])
  stop <-  as.numeric(scaaCNA$cnaStop[i])
  length <- stop-start
  wd <- x[which(x$chrom==chr)]
  
  if (start < wd$length[which(wd$arm == "p")] & stop < wd$length[which(wd$arm == "p")]) { # if CNA is in p
    pcArm[[i]] <- length / wd$length[wd$arm=="p"]
    arm[[i]] <- "p"
  } else if (start > wd$length[which(wd$arm == "p")] & stop > wd$length[which(wd$arm == "p")]) { # if CNA is in q
    pcArm[[i]] <- length / wd$length[wd$arm=="q"]
    arm[[i]] <- "q"
  } else if (start < wd$length[which(wd$arm == "p")] & stop > wd$length[which(wd$arm == "p")]) {
    pcArm[[i]] <- "centrospan"
    arm[[i]] <- "chr"
  }
}

scaaCNA$pcArm <- unlist(pcArm)
scaaCNA$arm <- unlist(arm)


# for recentered
pcArm.ploidyRecenter <- list()
arm.ploidyRecenter <- list()
for (i in 1:nrow(scaaCNA.ploidyRecenter)) {
  print(i)
  chr <- sub("\\:.*", "", scaaCNA.ploidyRecenter$peak[i])
  start <- as.numeric(scaaCNA.ploidyRecenter$cnaStart[i])
  stop <-  as.numeric(scaaCNA.ploidyRecenter$cnaStop[i])
  length <- stop-start
  wd <- x[which(x$chrom==chr)]
  
  if (start < wd$length[which(wd$arm == "p")] & stop < wd$length[which(wd$arm == "p")]) { # if CNA is in p
    pcArm.ploidyRecenter[[i]] <- length / wd$length[wd$arm=="p"]
    arm.ploidyRecenter[[i]] <- "p"
  } else if (start > wd$length[which(wd$arm == "p")] & stop > wd$length[which(wd$arm == "p")]) { # if CNA is in q
    pcArm.ploidyRecenter[[i]] <- length / wd$length[wd$arm=="q"]
    arm.ploidyRecenter[[i]] <- "q"
  } else if (start < wd$length[which(wd$arm == "p")] & stop > wd$length[which(wd$arm == "p")]) {
    pcArm.ploidyRecenter[[i]] <- "centrospan"
    arm.ploidyRecenter[[i]] <- "chr"
  }
}

scaaCNA.ploidyRecenter$pcArm <- unlist(pcArm.ploidyRecenter)
scaaCNA.ploidyRecenter$arm <- unlist(arm.ploidyRecenter)



# FIND TRUE DIPLOIDS (ALLELIC CHECK) ----------------------------------------------------

# concious of LOH and allelic imbalance
sample <- names(dwgs)
patient <- sub("\\..*", "", sample)

trueDiploid.perpatient <- list()

for (i in 1:length(unique(patient))) {
  index <- grep(unique(patient)[i], sample)
  list <- list()
  
  for (j in 1:length(index)) {
    list[[j]] <- data.frame(chromosome = substr(dwgs[[index[j]]]$chromosome, 4, nchar(dwgs[[index[j]]]$chromosome)),
                            start = dwgs[[index[j]]]$start,
                            stop = dwgs[[index[j]]]$stop,
                            cna = rowSums(dwgs[[index[j]]][ ,c(11:12)]==1) )
    colnames(list[[j]])[4] <- sample[index[j]]
  }
  
  trueDiploid.perpatient[[i]] <- powerjoin::power_full_join(list, by = c("chromosome","start","stop"))
}
names(trueDiploid.perpatient) <- unique(patient)

# only keep the true diploid
trueDiploid.perpatient <- lapply(trueDiploid.perpatient, function(x) {
  x[,-c(1:3)] <- ifelse(x[,-c(1:3)]==2,2,0)
  colnames(x)[1] <- "chr"
  x
})



# PULL TRUE DIPLOID PER PEAK ----------------------------------------------------------------

# pull regions with a scaa
peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")]

# list with element per scaa detailing any cnas
diploidPerScaa <- list
for (i in 1:nrow(scaa)) {
  print(i)
  peak <- rownames(scaa)[i]
  
  # pull cna data for the patient
  list <- list
  for (j in 1:ncol(scaa)) {
    patient <- substr(colnames(scaa)[j],1,4)
    scaaTRUE <- scaa[i,j]
    dwgs.wd <- trueDiploid.perpatient[[patient]]
    
    # pull the cnas overlapping the peak
    dwgs.wd <- dwgs.wd[dwgs.wd$chr==substr(sub("\\:.*", "", peak),4,nchar(sub("\\:.*", "", peak))), ]
    for (k in 1:nrow(dwgs.wd)) {
      start <- as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", peak))
      stop <- as.numeric(sub(".*-", "", peak))
      if ( DescTools::Overlap(c(start,stop), c(dwgs.wd$start[k], dwgs.wd$stop[k])) >0 ) {
        cnaStart <- dwgs.wd$start[k]
        cnaStop <- dwgs.wd$stop[k]
        segSize <- dwgs.wd$stop[k]- dwgs.wd$start[k]
        n <- ncol(dwgs.wd)-3
        nDip <- sum(dwgs.wd[k, -c(1:3)]==2)
        fracDip <- round(nDip/n, 2)
        list[[j]] <- setNames(c(patient, peak, scaaTRUE,cnaStart, cnaStop, segSize, nDip, fracDip), 
                              c("patient", "peak", "scaaTRUE","cnaStart", "cnaStop", "segSize", "nDip", "fracDip"))
        break
      } else {
        next
      }
    }
  }
  diploidPerScaa[[i]] <- do.call(rbind, list)
}
scaaCNA_dip <- do.call(rbind, diploidPerScaa)
scaaCNA_dip <- data.frame(scaaCNA_dip)
scaaCNA_dip <- na.omit(scaaCNA_dip)
scaaCNA_dip <- scaaCNA_dip[scaaCNA_dip$fracDip>0,]
scaaCNA_dip[ , -c(1:3)] <- apply(scaaCNA_dip[ , -c(1:3)], 2, function(x) as.numeric(as.character(x)))

# CHARACTERISE TRUE DIPLOIDS ----------------------------------------------------------------

# find length of diploid
diploidInfo<- list()
for (i in 1:length(trueDiploid.perpatient)) {
  wd <- trueDiploid.perpatient[[i]]
  
  # list subclonal diploid per patient
  diploid <- tidyr::gather(wd, key = "sample", value = "cn", 4:ncol(wd))
  diploid <- diploid[diploid$cn == 2, ]
  
  # find length
  diploid$size <- diploid$stop - diploid$start
  
  # some diploid may occur in >1 sample (althouhg still subclonal)
  freq <- setDT(diploid[,c(1:3,5)])[,list(countSame=.N),names(diploid[,c(1:3,5)])]
  
  # for a given subclonal gain, whats the average copy number and SD
  aveCN <- aggregate(diploid$cn, by=list(chr=diploid$chr, start=diploid$start, stop=diploid$stop), FUN=mean)
  
  diploid <- merge(diploid,freq)
  diploid <- merge(diploid,aveCN)
  colnames(diploid)[ncol(diploid)] <- "aveCN"
  colnames(diploid)[2] <- "cnaStart"
  colnames(diploid)[3] <- "cnaStop"
  diploid$nsample <- ncol(trueDiploid.perpatient[[i]])-3
  diploid$patient <- names(trueDiploid.perpatient)[i]
  
  diploid$cna <- "trueDiploid"
  freq <- setDT(diploid[,c(1:3,11)])[,list(countSameCNA=.N),names(diploid[,c(1:3,11)])]
  diploid <- merge(diploid,freq)
  diploid <- diploid[,c(11,6,10,1,2,3,7,4,12,5,8,9)]
  
  diploidInfo[[i]] <- diploid
}
names(diploidInfo) <- unique(patient) 
x <- do.call(rbind, diploidInfo)

# DEFINE FOCAL ON TRUE DIPLOIDS -------------------------------------------------------

# long segments will appear to be enriched because of chance
# you can compare the relative enrichment for SCAAs (as compared to diploid) between long vs short gains
# need to decide what level of clonality to consider a region a gain or diploid

# label CNAs if focal or arm length
library(data.table)
x <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", 
           col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
x <- x[ , .(length = sum(chromEnd - chromStart)), 
        by = .(chrom, arm = substring(name, 1, 1)) ]

# for non recentered
for (p in 1:length(diploidInfo)) {
  print(p)
  pcArm <- list()
  arm <- list()
  for (i in 1:nrow(diploidInfo[[p]])) {
    y <- diploidInfo[[p]]
    chr <- y$chr[i]
    start <- as.numeric(y$cnaStart[i])
    stop <-  as.numeric(y$cnaStop[i])
    length <- stop-start
    wd <- x[which(substr(x$chrom,4,nchar(x$chrom))==chr)]
    
    if (start < wd$length[which(wd$arm == "p")] & stop < wd$length[which(wd$arm == "p")]) { # if CNA is in p
      pcArm[[i]] <- length / wd$length[wd$arm=="p"]
      arm[[i]] <- "p"
    } else if (start > wd$length[which(wd$arm == "p")] & stop > wd$length[which(wd$arm == "p")]) { # if CNA is in q
      pcArm[[i]] <- length / wd$length[wd$arm=="q"]
      arm[[i]] <- "q"
    } else if (start < wd$length[which(wd$arm == "p")] & stop > wd$length[which(wd$arm == "p")]) {
      pcArm[[i]] <- "centrospan"
      arm[[i]] <- "chr"
    }
  }
  diploidInfo[[p]]$pcArm <- unlist(pcArm)
  diploidInfo[[p]]$arm <- unlist(arm)
}


# CHARACTERISE ALL SEGS (although diploid may be LOH)----------------------------------------------------------------

library(data.table)
# find length of X
segInfo<- list()
for (i in 1:length(dwgs.perpatient)) {
  wd <- dwgs.perpatient[[i]]
  
  # list subclonal X per patient
  x <- tidyr::gather(wd, key = "sample", value = "cn", 4:ncol(wd))
  
  # find length
  x$size <- x$stop - x$start
  
  # some X may occur in >1 sample (althouhg still subclonal)
  freq <- setDT(x[,c(1:3,5)])[,list(countSame=.N),names(x[,c(1:3,5)])]
  
  # for a given subclonal gain, whats the average copy number and SD
  aveCN <- aggregate(x$cn, by=list(chr=x$chr, start=x$start, stop=x$stop), FUN=mean)
  
  x <- merge(x,freq)
  x <- merge(x,aveCN)
  colnames(x)[ncol(x)] <- "aveCN"
  colnames(x)[2] <- "cnaStart"
  colnames(x)[3] <- "cnaStop"
  x$nsample <- ncol(dwgs.perpatient[[i]])-3
  x$patient <- names(dwgs.perpatient)[i]
  
  x$cna <- ifelse(x$cn>2, "gain", ifelse(x$cn<2, "loss", "diploid"))
  freq <- setDT(x[,c(1:3,11)])[,list(countSameCNA=.N),names(x[,c(1:3,11)])]
  x <- merge(x,freq)
  x <- x[,c(11,6,10,1,2,3,7,4,12,5,8,9)]
    
  segInfo[[i]] <- x
}
names(segInfo) <- unique(patient) 
x <- do.call(rbind,segInfo)

# and with recentre
# find length of X
segInfo.ploidyRecentre<- list()
for (i in 1:length(dwgs.ploidyRecentre.perpatient)) {
  wd <- dwgs.ploidyRecentre.perpatient[[i]]
  
  # list subclonal X per patient
  x <- tidyr::gather(wd, key = "sample", value = "cn", 4:ncol(wd))
  
  # find length
  x$size <- x$stop - x$start
  
  # some X may occur in >1 sample (althouhg still subclonal)
  freq <- setDT(x[,c(1:3,5)])[,list(countSame=.N),names(x[,c(1:3,5)])]
  
  # for a given subclonal gain, whats the average copy number and SD
  aveCN <- aggregate(x$cn, by=list(chr=x$chr, start=x$start, stop=x$stop), FUN=mean)
  
  x <- merge(x,freq)
  x <- merge(x,aveCN)
  colnames(x)[ncol(x)] <- "aveCN"
  x <- x[c(8,1,2,3,4,5,6,7,9)]
  colnames(x)[2] <- "cnaStart"
  colnames(x)[3] <- "cnaStop"
  x$nsample <- ncol(dwgs.ploidyRecentre.perpatient[[i]])-3
  x$patient <- names(dwgs.ploidyRecentre.perpatient)[i]
  
  x$cna <- ifelse(x$cn>2, "gain", ifelse(x$cn<2, "loss", "diploid"))
  freq <- setDT(x[,c(1:3,11)])[,list(countSameCNA=.N),names(x[,c(1:3,11)])]
  x <- merge(x,freq)
  x <- x[,c(11,6,10,1,2,3,7,4,12,5,8,9)]
  
  segInfo.ploidyRecentre[[i]] <- x
}
names(segInfo.ploidyRecentre) <- unique(patient) 
x <- do.call(rbind,segInfo.ploidyRecentre)

# DEFINE FOCAL ON ALL SEGMENTS -------------------------------------------------------

# long segments will appear to be enriched because of chance
# you can compare the relative enrichment for SCAAs (as compared to diploid) between long vs short gains
# need to decide what level of clonality to consider a region a gain or diploid

# label CNAs if focal or arm length
library(data.table)
x <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", 
           col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
x <- x[ , .(length = sum(chromEnd - chromStart)), 
        by = .(chrom, arm = substring(name, 1, 1)) ]

# for non recentered
for (p in 1:length(segInfo)) {
  print(p)
  pcArm <- list()
  arm <- list()
for (i in 1:nrow(segInfo[[p]])) {
  y <- segInfo[[p]]
  chr <- y$chr[i]
  start <- as.numeric(y$cnaStart[i])
  stop <-  as.numeric(y$cnaStop[i])
  length <- stop-start
  wd <- x[which(substr(x$chrom,4,nchar(x$chrom))==chr)]
  
  if (start < wd$length[which(wd$arm == "p")] & stop < wd$length[which(wd$arm == "p")]) { # if CNA is in p
    pcArm[[i]] <- length / wd$length[wd$arm=="p"]
    arm[[i]] <- "p"
  } else if (start > wd$length[which(wd$arm == "p")] & stop > wd$length[which(wd$arm == "p")]) { # if CNA is in q
    pcArm[[i]] <- length / wd$length[wd$arm=="q"]
    arm[[i]] <- "q"
  } else if (start < wd$length[which(wd$arm == "p")] & stop > wd$length[which(wd$arm == "p")]) {
    pcArm[[i]] <- "centrospan"
    arm[[i]] <- "chr"
  }
}
  segInfo[[p]]$pcArm <- unlist(pcArm)
  segInfo[[p]]$arm <- unlist(arm)
}

# for  recentered
for (p in 1:length(segInfo.ploidyRecentre)) {
  print(p)
  pcArm <- list()
  arm <- list()
  for (i in 1:nrow(segInfo.ploidyRecentre[[p]])) {
    y <- segInfo.ploidyRecentre[[p]]
    chr <- y$chr[i]
    start <- as.numeric(y$cnaStart[i])
    stop <-  as.numeric(y$cnaStop[i])
    length <- stop-start
    wd <- x[which(substr(x$chrom,4,nchar(x$chrom))==chr)]
    
    if (start < wd$length[which(wd$arm == "p")] & stop < wd$length[which(wd$arm == "p")]) { # if CNA is in p
      pcArm[[i]] <- length / wd$length[wd$arm=="p"]
      arm[[i]] <- "p"
    } else if (start > wd$length[which(wd$arm == "p")] & stop > wd$length[which(wd$arm == "p")]) { # if CNA is in q
      pcArm[[i]] <- length / wd$length[wd$arm=="q"]
      arm[[i]] <- "q"
    } else if (start < wd$length[which(wd$arm == "p")] & stop > wd$length[which(wd$arm == "p")]) {
      pcArm[[i]] <- "centrospan"
      arm[[i]] <- "chr"
    }
  }
  segInfo.ploidyRecentre[[p]]$pcArm <- unlist(pcArm)
  segInfo.ploidyRecentre[[p]]$arm <- unlist(arm)
}


# SAVE -----------
saveRDS(dwgs.perpatient, "~/Documents/SCAA/Data/dwgs.perpatient.rds")
saveRDS(dwgs.ploidyRecentre.perpatient, "~/Documents/SCAA/Data/dwgs.ploidyRecentre.perpatient.rds")

saveRDS(ploidy.persample, "~/Documents/SCAA/Data/ploidy.persample.rds")
saveRDS(ploidy.perpatient, "~/Documents/SCAA/Data/ploidy.perpatient.rds")

saveRDS(diploidPerScaa,"~/Documents/SCAA/Data/diploidPerScaa.rds")
saveRDS(cnaPerScaa,"~/Documents/SCAA/Data/cnaPerScaa.rds")
saveRDS(cnaPerScaa.ploidyRecenter,"~/Documents/SCAA/Data/cnaPerScaa.ploidyRecenter.rds")

saveRDS(focal, "~/Documents/SCAA/Data/focal.rds")
saveRDS(arm, "~/Documents/SCAA/Data/arm.rds")
saveRDS(pcArm.ploidyRecenter, "~/Documents/SCAA/Data/pcArm.ploidyRecenter.rds")
saveRDS(arm.ploidyRecenter, "~/Documents/SCAA/Data/arm.ploidyRecenter.rds")

saveRDS(scaaCNA_dip, "~/Documents/SCAA/Data/scaaCNA_dip.rds")
saveRDS(scaaCNA, "~/Documents/SCAA/Data/scaaCNA.rds")
saveRDS(scaaCNA.ploidyRecenter, "~/Documents/SCAA/Data/scaaCNA.ploidyRecenter.rds")

saveRDS(gainInfo, "~/Documents/SCAA/Data/gainInfo.rds")
saveRDS(diploidInfo, "~/Documents/SCAA/Data/diploidInfo.rds")
saveRDS(lossInfo, "~/Documents/SCAA/Data/lossInfo.rds")
saveRDS(segInfo, "~/Documents/SCAA/Data/segInfo.rds")
saveRDS(segInfo.ploidyRecentre, "~/Documents/SCAA/Data/segInfo.ploidyRecentre.rds")

saveRDS(trueDiploid.perpatient, "~/Documents/SCAA/Data/trueDiploid.perpatient.rds")
saveRDS(trueDiploids, "~/Documents/SCAA/Data/trueDiploids.rds")



