# ADD IN CYTOBAND, TSS ENRICHMENT, AND TELOMERE/CENTROMERE DISTANCE, and SEGMENT PER SCAA =====

# cytoband data
hg38.chromosomeBand <- read.csv("~/Documents/SCAA/Data/UCSC_search/hg38.chromosomeBands.csv")

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

# collect peak locations
peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")]

scaa.bins <- data.frame(chr = substr(sub("\\:.*", "", rownames(scaa)),4,nchar(sub("\\:.*", "", rownames(scaa)))),
                        peakStart = as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", rownames(scaa))),
                        peakStop = as.numeric(sub(".*-", "", rownames(scaa))) )

# TSS data
mySession = browserSession("UCSC")
genome(mySession) <- "hg38"
range <- GRanges(paste("chr", scaa.bins$chr[c(1:100)], sep = ""), IRanges(scaa.bins$peakStart[c(1:100)], scaa.bins$peakStop[c(1:100)]))
track.names <- trackNames(ucscTableQuery(mySession))
tableNames(ucscTableQuery(mySession, track="geneHancer"))
TSS <- getTable(ucscTableQuery(mySession, track="geneHancer",
                               range=range, table="geneHancerGenesDoubleElite"))

# assign to each peak
stain <- list()
bandName <- list()
dist.centro <- list()
dist.telo <- list()

for (i in 1:nrow(scaa.bins)) {
  print(paste(i,"/",nrow(scaa.bins)))
  chr <- paste("chr",scaa.bins$chr[i], sep = "")
  start <- as.numeric(scaa.bins$peakStart[i])
  stop <- as.numeric(scaa.bins$peakStop[i])
  
  # cytoband per peak
  stain[i] <- "no data"
  bandName[i] <- "no data"
  wd <- hg38.chromosomeBand[hg38.chromosomeBand$X.chrom==chr,]
  for (j in 1:nrow(wd)) {
    if( nrow(wd)==0 ) {
      stain[i] <- "no data"
      bandName[i] <- "no data"
      break
    }
    if( start>=wd$chromStart[j] & stop<=wd$chromEnd[j] ) {
      stain[i] <- wd$gieStain[j]
      bandName[i] <- wd$name[j]
      break
    } else {
      next
    }
  }
  
  # centromere per peak
  wd <- hg38centromere[hg38centromere$chr == chr,]
  dist.centro[[i]] <- min(abs(stop - wd$centroStop), abs(start - wd$centroStart))
  
  # telomere per peak
  wd <- hg38.length[hg38.length$chr == substr(chr, 4, nchar(chr)),]
  dist.telo[[i]] <- min(abs(start - 1), abs(stop - wd$length))
  
}
scaa.bins$geiStain <- unlist(stain)
scaa.bins$bandName <- unlist(bandName)
scaa.bins$dist.centro <- unlist(dist.centro)
scaa.bins$dist.telo <- unlist(dist.telo)

# add TSS per peak
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9450098/ says promoter is 3kb +/- from TSS but didnt make differe ce
TSS.n <- list()
for (i in 1:nrow(scaa.bins)) {
  print(paste(i,"/",nrow(scaa.bins)))
  chr <- paste("chr", scaa.bins$chr[i], sep = "")
  start <- scaa.bins$peakStart[i]-3000
  stop <- scaa.bins$peakStop[i]+3000
  
  wd <- TSS[TSS$chrom==chr,]
  for (j in 1:nrow(wd)) {
    count <- 0
    if (nrow(wd)==0) {
      break
    }
    if( (wd$chromStart[j])>=(start) & (wd$chromEnd[j])<=(stop)) {
      count <- count + 1
    } else {
      next
    }
  }
  TSS.n[[i]] <- count
}
scaa.bins$TSS.n <- unlist(TSS.n)

# CpG data
library(rtracklayer)
mySession = browserSession("UCSC")
genome(mySession) <- "hg38"
range <- GRanges(paste("chr", scaa.bins$chr[c(1:10)], sep = ""), IRanges(scaa.bins$peakStart[c(1:10)], scaa.bins$peakStop[c(1:10)]))
CpG <- getTable(ucscTableQuery(mySession, track="cpgIslandExt",
                               range=range, table="cpgIslandExt"))

# add CpG island coverage per peak
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9450098/ says promoter is 3kb +/- from TSS but didnt make differe ce
CpG.n <- list()
for (i in 1:nrow(scaa.bins)) {
  print(paste(i,"/",nrow(scaa.bins)))
  chr <- paste("chr", scaa.bins$chr[i], sep = "")
  start <- scaa.bins$peakStart[i]-3000
  stop <- scaa.bins$peakStop[i]+3000
  
  wd <- CpG[CpG$chrom==chr,]
  for (j in 1:nrow(wd)) {
    count <- 0
    if (nrow(wd)==0) {
      break
    }
    if( (wd$chromStart[j])>=(start) & (wd$chromEnd[j])<=(stop)) {
      count <- count + 1
    } else {
      next
    }
  }
  CpG.n[[i]] <- count
}
scaa.bins$CpG.n <- unlist(CpG.n)

scaa.long <- cbind(scaa.bins, scaa)
scaa.long <- gather(scaa.long, key="patient", value = "value", -chr, -peakStart, -peakStop, -geiStain,-bandName, -dist.centro, -dist.telo, -TSS.n, -CpG.n)
scaa.long$patient <- substr(scaa.long$patient, 1, 4)

# add segment per peak (needs to be patient specific)
cnaStart <- list()
cnaStop <- list()
segSize <- list()
n.samples <- list()
fracDip <- list()
fracLoss <- list()
fracGain <- list()
distancetoBP <- list()
for (i in 1:nrow(scaa.long)) {
  print(paste(i,"/",nrow(scaa.long)))
  patient <- scaa.long$patient[i]
  
  chr <- scaa.long$chr[i]
  peakStart <- as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", scaa.long$peakStart[i]))
  peakStop <- as.numeric(sub(".*-", "", scaa.long$peakStop[i]))
  
  dwgs.wd <- dwgs.perpatient[[patient]]
  dwgs.wd <- dwgs.wd[dwgs.wd$chr == chr, ]
  
  cnaStart[[i]] <- "no data"
  cnaStop[[i]] <- "no data"
  segSize[[i]] <- "no data"
  n.samples[[i]] <- "no data"
  fracDip[[i]] <- "no data"
  fracLoss[[i]] <- "no data"
  fracGain[[i]] <- "no data"
  distancetoBP[[i]] <- "no data"
  
  if(nrow(dwgs.wd)==0) { # no chromosome data
    next
  }
  
  for (k in 1:nrow(dwgs.wd)) {
    if ( DescTools::Overlap(c(peakStart, peakStop), c(dwgs.wd$start[k], dwgs.wd$stop[k])) >= 250 ) { # choose segment which majority of peak sits on
      cnaStart[[i]] <- dwgs.wd$start[k]
      cnaStop[[i]] <- dwgs.wd$stop[k]
      segSize[[i]] <- dwgs.wd$stop[k] - dwgs.wd$start[k]
      n.samples[[i]] <- ncol(dwgs.wd)-3
      fracDip[[i]] <- round( (sum(dwgs.wd[k, -c(1:3)]==2, na.rm = T)) / (ncol(dwgs.wd)-3), 2 )
      fracLoss[[i]] <- round( (sum(dwgs.wd[k, -c(1:3)]<2, na.rm = T)) / (ncol(dwgs.wd)-3), 2 )
      fracGain[[i]] <- round( (sum(dwgs.wd[k, -c(1:3)]>2, na.rm = T)) / (ncol(dwgs.wd)-3), 2 )
      peakMid <- peakStart + 250
      distancetoBP[[i]] <- min((peakMid-dwgs.wd$start[k]), (dwgs.wd$stop[k]-peakMid))
      break
    } else {
      next
    }
  }
}

scaa.long$cnaStart <- as.numeric(unlist(cnaStart))
scaa.long$cnaStop <- as.numeric(unlist(cnaStop))
scaa.long$segSize <- as.numeric(unlist(segSize))
scaa.long$n.samples <- as.numeric(unlist(n.samples))
scaa.long$fracDip <- as.numeric(unlist(fracDip))
scaa.long$fracLoss <- as.numeric(unlist(fracLoss))
scaa.long$fracGain <- as.numeric(unlist(fracGain))
scaa.long$distancetoBP <- as.numeric(unlist(distancetoBP))

x <- scaa.long[which(scaa.long$fracGain == 0.5 & scaa.long$fracLoss ==0.5),] # no peaks are half gain/half loss
scaa.long$cna <- ifelse(scaa.long$fracGain>=0.5, "gain", ifelse(scaa.long$fracLoss>=0.5, "loss", "diploid"))
scaa.long$mainClonality <- pmax(scaa.long$fracDip, scaa.long$fracLoss, scaa.long$fracGain)

peakAnn_all <- readRDS("~/Documents/SCAA/Data/peakAnn_all.rds")
promoter <- peakAnn_all$peak[which(peakAnn_all$Promoter==TRUE)]
enhancer <- peakAnn_all$peak[which(!is.na(peakAnn_all$elementType))]
scaa.long$peak <- paste("chr",scaa.long$chr,":",scaa.long$peakStart,"-",scaa.long$peakStop, sep = "")
scaa.long$promoter <- ifelse(scaa.long$peak %in% promoter, "TRUE","FALSE")
scaa.long$enhancer <- ifelse(scaa.long$peak %in% enhancer, "TRUE","FALSE")

scaa.long$sizeGroup <- "100-150"
scaa.long$sizeGroup[scaa.long$segSize<=100000000] <- "75-100"
scaa.long$sizeGroup[scaa.long$segSize<=75000000] <- "50-75"
scaa.long$sizeGroup[scaa.long$segSize<=50000000] <- "25-50"
scaa.long$sizeGroup[scaa.long$segSize<=25000000] <- "1-25"
scaa.long$sizeGroup[scaa.long$segSize<=1000000] <- "0-1"
scaa.long$sizeGroup <- factor(scaa.long$sizeGroup, levels = c("0-1","1-25","25-50","50-75","75-100","100-150"))

scaa.long$geneID <- peakAnn_all[match(scaa.long$peak, peakAnn_all$peak), 14]
scaa.long$symbol <- peakAnn_all[match(scaa.long$peak, peakAnn_all$peak), 17]
library("org.Hs.eg.db")
scaa.long$entrez <- mapIds(org.Hs.eg.db, keys = scaa.long$symbol,
                           column = c('ENTREZID'), keytype = 'SYMBOL')

# length
library(rtracklayer)
hg38.length <- SeqinfoForUCSCGenome("hg38")
hg38.length <- data.frame("chr"=substr(hg38.length@seqnames,4,nchar(hg38.length@seqnames)), "length"=hg38.length@seqlengths)
hg38.length <- hg38.length[c(1:22),]

dist.centroPC <- list()
dist.teloPC <- list()
for (i in 1:nrow(scaa.long)) {
  print(i)
  chr1 <- paste("chr",scaa.long$chr[i], sep="")
  chr <- scaa.long$chr[i]
  centroStart <- min(hg38.centromere$chromStart[which(hg38.centromere$chrom==chr1)])
  centroStop <- min(hg38.centromere$chromEnd[which(hg38.centromere$chrom==chr1)])
  length <- hg38.length$length[which(hg38.length$chr==chr)]
  if (scaa.long$peakStart[i]<centroStart) {
    dist.centroPC[[i]] <- scaa.long$dist.centro[i]/centroStart
    dist.teloPC[[i]] <- scaa.long$dist.telo[i]/length
  } else {
    dist.centroPC[[i]] <- scaa.long$dist.centro[i]/(length-centroStop)
    dist.teloPC[[i]] <- scaa.long$dist.telo[i]/length
  }
}

scaa.long$dist.centroPC <- unlist(dist.centroPC)
scaa.long$dist.teloPC <- unlist(dist.teloPC)
# save
saveRDS(scaa.long, "~/Documents/SCAA/Data/scaa.long.rds")
