

# FIND TRUE DIPLOIDS (ALLELIC CHECK) =====

# conscious of LOH and allelic imbalance
# per sample create list of the segments which are true diploid and not LOH (appear diploid)

trueDiploid.persample <- list()
for (j in 1:length(dwgs)) {
  x <- data.frame(patient = substr(names(dwgs)[j], 1, 4),
                  sample = names(dwgs)[j],
                  chromosome = substr(dwgs[[j]]$chromosome, 4, nchar(dwgs[[j]]$chromosome)),
                  start = dwgs[[j]]$start,
                  stop = dwgs[[j]]$stop,
                  cna = rowSums(dwgs[[j]][ ,c(11:12)]==1) )
  trueDiploid.persample[[j]] <- na.omit(x[x$cna==2,])
}
names(trueDiploid.persample) <- names(dwgs)
trueDiploids <- do.call(rbind, trueDiploid.persample)
trueDiploids$segID <- paste(trueDiploids$chromosome, ":", trueDiploids$start, "-", trueDiploids$stop, sep = "")
trueDiploids <- trueDiploids[c(7,1:6)]

# PULL PEAK/SCAA LOCI =====

# for carcinomas only
peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")]

scaa.bins <- data.frame(chr = substr(sub("\\:.*", "", rownames(scaa)),4,nchar(sub("\\:.*", "", rownames(scaa)))),
         peakStart = as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", rownames(scaa))),
         peakStop = as.numeric(sub(".*-", "", rownames(scaa))) )

# CHARACTERISE ALL SEGS =====

library(data.table)

segments <- list()
for (i in 1:length(dwgs.perpatient)) {
  wd <- dwgs.perpatient[[i]]
  
  # list subclonal X per patient
  x <- tidyr::gather(wd, key = "sample", value = "cn", 4:ncol(wd))
  
  # find length
  x$size <- x$stop - x$start

  # some X may occur in >1 sample (althouhg still subclonal)
  freq <- setDT(x[,c(1:3,5)])[,list(countSame=.N),names(x[,c(1:3,5)])]
  
  # for a given subclonal gain, whats the average copy number and SD
  aveCN <- aggregate(x$cn, by=list(chromosome=x$chromosome, start=x$start, stop=x$stop), FUN=mean, na.rm = TRUE)
  
  x <- merge(x,freq)
  x <- merge(x,aveCN)
  colnames(x)[ncol(x)] <- "aveCN"
  
  colnames(x)[1] <- "chr"
  colnames(x)[2] <- "cnaStart"
  colnames(x)[3] <- "cnaStop"
  
  x$nsample <- ncol(dwgs.perpatient[[i]])-3
  x$patient <- names(dwgs.perpatient)[i]
  x$segID <- paste(x$chr, ":", x$cnaStart, "-", x$cnaStop, sep = "")
  
  x$cna <- ifelse(x$cn>2, "gain", ifelse(x$cn<2, "loss", "diploid"))
  
  aveCNperCNA <- aggregate(x$cn, by=list(chr=x$chr, cnaStart=x$cnaStart, cnaStop=x$cnaStop, cna=x$cna), FUN=mean, na.rm = TRUE)
  x <- merge(x,aveCNperCNA)
  colnames(x)[ncol(x)] <- "aveCNperCNA"
  
  # tag the true diploids
  list <- list()
  for (j in 1:nrow(x)) {
    sample <- x$sample[j]
    wd2 <- trueDiploids[trueDiploids$sample==sample,]
    list[[j]] <- ifelse(x$cna[j] == "diploid" & x$segID[j] %in% wd2$segID, "diploid", 
                    ifelse(x$cna[j] == "diploid" & x$segID[j] %!in% wd2$segID, "LOH", 
                         ifelse(x$cna[j] == "gain", "gain", ifelse(x$cna[j] == "loss","loss", NA))))
  }
  
  x$cna <- unlist(list)
  
  freq <- setDT(x[,c(1:4)])[,list(countSameCNA=.N),names(x[,c(1:4)])]
  x <- merge(x,freq)
  x <- x[,c(11,6,10,1,2,3,12,7,4,5,8,9,14,13)]
  x$clonalityCNA <- x$countSameCNA / x$nsample 
  
  # find number of peaks/scaa per segment (as long as over half peak in segment)
  n.peaks <- list()
  n.SCAA <- list()
  
  for (j in 1:nrow(x)) {
    print(paste(i,"/",length(dwgs.perpatient),"-",j,"/",nrow(x)))
    chr <- x$chr[j]
    start <- x$cnaStart[j]
    stop <- x$cnaStop[j]
    
    wd <- scaa[which(scaa.bins$chr==chr & (scaa.bins$peakStart>=(start-250)) & (scaa.bins$peakStop<=(stop+250)) ), substr(colnames(scaa),1,4) == unique(x$patient)]
    n.peaks[[j]] <- length(wd)
    
    # calculate scaa per segment
    n.SCAA[[j]] <- sum(wd==1)
    
    #index.overlap <- list()
    #for (k in 1:nrow(scaa.bins)) {
    #  if ( scaa.bins$chr[k] == chr & DescTools::Overlap(c(start,stop), c(scaa.bins$peakStart[k], scaa.bins$peakStop[k])) >= 250 ) {
    #    index.overlap[[k]] <- TRUE
    #  } else {
    #    index.overlap[[k]] <- FALSE
    #  }
    #}
    
    #n.peaks[[j]] <- sum(index.overlap==TRUE)
    
    # calculate scaa per segment
    #z <- scaa[index.overlap==TRUE, substr(colnames(scaa),1,4) == unique(x$patient)]
    #n.SCAA[[j]] <- sum(z==1)
  }
  
  x$n.peaks <- unlist(n.peaks)
  x$n.SCAA <- unlist(n.SCAA)
  
  x$prop.SCAA <- x$n.SCAA / x$n.peaks
  
  segments[[i]] <- x
}

names(segments) <- names(dwgs.perpatient)
x <- do.call(rbind, segments)
sum(x$cna=="diploid", na.rm = T)


# DEFINE FOCAL CNAS =====

# long segments will appear to be enriched because of chance
# you can compare the relative enrichment for SCAAs (as compared to diploid) between long vs short gains
# need to decide what level of clonality to consider a region a gain or diploid

# label CNAs if focal or arm length
library(data.table)
armLength <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", 
           col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
armLength <- armLength[ , .(length = sum(chromEnd - chromStart)), 
        by = .(chrom, arm = substring(name, 1, 1)) ]

# for non recentered
for (j in 1:length(segments)) {
  print(paste(j,"/",length(segments)))
  wd1 <- segments[[j]]
  pcArm <- list()
  arm <- list()
  
  for (i in 1:nrow(wd1)) {
    chr <- paste("chr", wd1$chr[i], sep = "")
    start <- as.numeric(wd1$cnaStart[i])
    stop <-  as.numeric(wd1$cnaStop[i])
    length <- as.numeric(wd1$size[i])
    wd2 <- armLength[which(armLength$chrom==chr)]
    
    if (start < wd2$length[which(wd2$arm == "p")] & stop < wd2$length[which(wd2$arm == "p")]) { # if CNA is in p
      pcArm[[i]] <- length / wd2$length[wd2$arm=="p"]
      arm[[i]] <- "p"
    } else if (start > wd2$length[which(wd2$arm == "p")] & stop > wd2$length[which(wd2$arm == "p")]) { # if CNA is in q
      pcArm[[i]] <- length / wd2$length[wd2$arm=="q"]
      arm[[i]] <- "q"
    } else if (start < wd2$length[which(wd2$arm == "p")] & stop > wd2$length[which(wd2$arm == "p")]) {
      pcArm[[i]] <- "centrospan"
      arm[[i]] <- "chr"
    }
  }
  segments[[j]]$pcArm <- unlist(pcArm)
  segments[[j]]$arm <- unlist(arm)
}
names(segments) <- names(dwgs.perpatient)

# SEGMENTS ON PLOIDY RECENTERED =====


library(data.table)

segments_ploidyrecenterd <- list()
for (i in 1:length(dwgs.ploidyRecentre.perpatient)) {
  wd <- dwgs.ploidyRecentre.perpatient[[i]]
  
  # list subclonal X per patient
  x <- tidyr::gather(wd, key = "sample", value = "cn", 4:ncol(wd))
  
  # find length
  x$size <- x$stop - x$start
  
  # some X may occur in >1 sample (althouhg still subclonal)
  freq <- setDT(x[,c(1:3,5)])[,list(countSame=.N),names(x[,c(1:3,5)])]
  
  # for a given subclonal gain, whats the average copy number and SD
  aveCN <- aggregate(x$cn, by=list(chromosome=x$chromosome, start=x$start, stop=x$stop), FUN=mean, na.rm = TRUE)
  
  x <- merge(x,freq)
  x <- merge(x,aveCN)
  colnames(x)[ncol(x)] <- "aveCN"
  
  colnames(x)[1] <- "chr"
  colnames(x)[2] <- "cnaStart"
  colnames(x)[3] <- "cnaStop"
  
  x$nsample <- ncol(dwgs.ploidyRecentre.perpatient[[i]])-3
  x$patient <- names(dwgs.ploidyRecentre.perpatient)[i]
  x$segID <- paste(x$chr, ":", x$cnaStart, "-", x$cnaStop, sep = "")
  
  x$cna <- ifelse(x$cn>2, "gain", ifelse(x$cn<2, "loss", "diploid"))
  
  aveCNperCNA <- aggregate(x$cn, by=list(chr=x$chr, cnaStart=x$cnaStart, cnaStop=x$cnaStop, cna=x$cna), FUN=mean, na.rm = TRUE)
  x <- merge(x,aveCNperCNA)
  colnames(x)[ncol(x)] <- "aveCNperCNA"
  
  # tag the true diploids
  list <- list()
  for (j in 1:nrow(x)) {
    sample <- x$sample[j]
    wd2 <- trueDiploids[trueDiploids$sample==sample,]
    list[[j]] <- ifelse(x$cna[j] == "diploid" & x$segID[j] %in% wd2$segID, "diploid", 
                        ifelse(x$cna[j] == "diploid" & x$segID[j] %!in% wd2$segID, "LOH", 
                               ifelse(x$cna[j] == "gain", "gain", ifelse(x$cna[j] == "loss","loss", NA))))
  }
  
  x$cna <- unlist(list)
  
  freq <- setDT(x[,c(1:4)])[,list(countSameCNA=.N),names(x[,c(1:4)])]
  x <- merge(x,freq)
  x <- x[,c(11,6,10,1,2,3,12,7,4,5,8,9,14,13)]
  x$clonalityCNA <- x$countSameCNA / x$nsample 
  
  # find number of peaks/scaa per segment (as long as over half peak in segment)
  n.peaks <- list()
  n.SCAA <- list()
  
  for (j in 1:nrow(x)) {
    print(paste(i,"/",length(dwgs.ploidyRecentre.perpatient),"-",j,"/",nrow(x)))
    chr <- x$chr[j]
    start <- x$cnaStart[j]
    stop <- x$cnaStop[j]
    
    wd <- scaa[which(scaa.bins$chr==chr & (scaa.bins$peakStart>=(start-250)) & (scaa.bins$peakStop<=(stop+250)) ), substr(colnames(scaa),1,4) == unique(x$patient)]
    n.peaks[[j]] <- length(wd)
    
    # calculate scaa per segment
    n.SCAA[[j]] <- sum(wd==1)
  }
  
  x$n.peaks <- unlist(n.peaks)
  x$n.SCAA <- unlist(n.SCAA)
  
  x$prop.SCAA <- x$n.SCAA / x$n.peaks
  
  segments_ploidyrecenterd[[i]] <- x
}

names(segments_ploidyrecenterd) <- names(dwgs.ploidyRecentre.perpatient)
x <- do.call(rbind, segments_ploidyrecenterd)
sum(x$cna=="diploid", na.rm = T)


# SAVE -----------
saveRDS(trueDiploids, "~/Documents/SCAA/Data/trueDiploids.rds")
saveRDS(segments, "~/Documents/SCAA/Data/segments.rds")
saveRDS(segments_ploidyrecenterd, "~/Documents/SCAA/Data/segments_ploidyrecenterd.rds")
saveRDS(scaa, "~/Documents/SCAA/Data/scaa.rds")


