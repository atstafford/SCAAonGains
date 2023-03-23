# PREREQUISITIS: load section 4 data 

# epicc showed that dirving genetic alterations were sort of rare, but loads of SCAAs in cancer genes were recurrent and recurrent SCAAs in other gene not ass with cancer (yet...)
# most scaas are clonal

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

# TSS data
mySession = browserSession("UCSC")
genome(mySession) <- "hg38"
range <- GRanges(paste("chr", scaa.bins$chr[c(1:100)], sep = ""), IRanges(scaa.bins$peakStart[c(1:100)], scaa.bins$peakStop[c(1:100)]))
track.names <- trackNames(ucscTableQuery(mySession))
tableNames(ucscTableQuery(mySession, track="geneHancer"))
TSS <- getTable(ucscTableQuery(mySession, track="geneHancer",
                               range=range, table="geneHancerGenesDoubleElite"))

# collect peak locations
peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")]

scaa.bins <- data.frame(chr = substr(sub("\\:.*", "", rownames(scaa)),4,nchar(sub("\\:.*", "", rownames(scaa)))),
                        peakStart = as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", rownames(scaa))),
                        peakStop = as.numeric(sub(".*-", "", rownames(scaa))) )

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

# SCAA ENRICHMENT BY SEGMENT SIZE =====

unique(scaa.long$sizeGroup)
x <- data.frame(SCAA= c(nrow(scaa.long[scaa.long$sizeGroup=="0-1" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$sizeGroup=="1-25" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$sizeGroup=="25-50" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$sizeGroup=="50-75" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$sizeGroup=="75-100" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$sizeGroup=="100-150" & scaa.long$value==1,])),
                notSCAA= c(nrow(scaa.long[scaa.long$sizeGroup=="0-1" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$sizeGroup=="1-25" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$sizeGroup=="25-50" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$sizeGroup=="50-75" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$sizeGroup=="75-100" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$sizeGroup=="100-150" & scaa.long$value==0,])))
rownames(x) <- c("0-1","1-25","25-50","50-75","75-100","100-150")

chisq <- chisq.test(x)
enrich <- chisq$observed/chisq$expected

jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
corrplot::corrplot(chisq$residuals, is.cor = FALSE)
dev.off()

round(rowSums(x)/sum(rowSums(x)),3)

model <- glm(value ~ sizeGroup, data=scaa.long, family = "binomial")
summary(model)

# SCAA ENRICHMENT BY CNA =====
x <- data.frame(SCAA=c( nrow(scaa.long[scaa.long$value==1 & scaa.long$cna=="gain",]), nrow(scaa.long[scaa.long$value==1 & scaa.long$cna=="loss",]), nrow(scaa.long[scaa.long$value==1 & scaa.long$cna=="diploid",])  ),
                noSCAA=c( nrow(scaa.long[scaa.long$value==0 & scaa.long$cna=="gain",]), nrow(scaa.long[scaa.long$value==0 & scaa.long$cna=="loss",]), nrow(scaa.long[scaa.long$value==0 & scaa.long$cna=="diploid",]) ))
rownames(x) <- c("gain","loss","diploid")

chisq <- chisq.test(x)
enrich <- chisq$observed/chisq$expected

jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
corrplot::corrplot(chisq$residuals, is.cor = FALSE)
dev.off()

round(rowSums(x)/sum(rowSums(x)),2)

model <- glm(value ~ cna, data=scaa.long, family = "binomial")
summary(model)


model <- glm(value ~ geiStain + dist.telo + dist.centro + cna, data=scaa.long, family = "binomial")
summary(model)


# SCAA ENRICHMENT BY G-BAND =====

# promoters / non promoters peaks across G-bands

x <- data.frame("Promoter/enhancer peak"= c(nrow(scaa.long[scaa.long$geiStain=="gneg" & (scaa.long$promoter=="TRUE" | scaa.long$enhancer=="TRUE"),]),
                        nrow(scaa.long[scaa.long$geiStain=="gvar" &  (scaa.long$promoter=="TRUE" | scaa.long$enhancer=="TRUE"),]),
                        nrow(scaa.long[scaa.long$geiStain=="acen" &  (scaa.long$promoter=="TRUE" | scaa.long$enhancer=="TRUE"),]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos25" &  (scaa.long$promoter=="TRUE" | scaa.long$enhancer=="TRUE"),]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos50" &  (scaa.long$promoter=="TRUE" | scaa.long$enhancer=="TRUE"),]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos75" &  (scaa.long$promoter=="TRUE" | scaa.long$enhancer=="TRUE"),]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos100" &  (scaa.long$promoter=="TRUE" | scaa.long$enhancer=="TRUE"),])),
                "Non-promoter/enhancer peak"= c(nrow(scaa.long[scaa.long$geiStain=="gneg" &  (scaa.long$promoter=="FALSE" | scaa.long$enhancer=="FALSE"),]),
                           nrow(scaa.long[scaa.long$geiStain=="gvar" & (scaa.long$promoter=="FALSE" | scaa.long$enhancer=="FALSE"),]),
                           nrow(scaa.long[scaa.long$geiStain=="acen" & (scaa.long$promoter=="FALSE" | scaa.long$enhancer=="FALSE"),]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos25" & (scaa.long$promoter=="FALSE" | scaa.long$enhancer=="FALSE"),]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos50" & (scaa.long$promoter=="FALSE" | scaa.long$enhancer=="FALSE"),]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos75" & (scaa.long$promoter=="FALSE" | scaa.long$enhancer=="FALSE"),]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos100" & (scaa.long$promoter=="FALSE" | scaa.long$enhancer=="FALSE"),])))
rownames(x) <- c("gneg","gvar","acen","gpos25","gpos50","gpos75","gpos100")
x$Stain <- rownames(x)
x <- gather(x, key = "key", value="Count", -Stain)
x$fraction <- paste(round(100*(x$Count/sum(x$Count)),1), "%", sep="")
jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
ggplot(x, aes(y=Count, x=Stain, fill= key)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label=fraction), hjust = 0.5, position = position_dodge(width = 1), size=5) +
  coord_flip() +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  theme(
        legend.text = element_text(size=20), legend.title = element_blank(), legend.position = "top")
dev.off()



# SCAA / non SCAA peaks across G-bands
library(ggrepel)
round(rowSums(x)/sum(rowSums(x)),3)
round(colSums(x)/sum(colSums(x)),3)

x <- data.frame(SCAA_peak= c(nrow(scaa.long[scaa.long$geiStain=="gneg" & (scaa.long$value==1),]),
                                            nrow(scaa.long[scaa.long$geiStain=="gvar" &  (scaa.long$value==1),]),
                                            nrow(scaa.long[scaa.long$geiStain=="acen" &  (scaa.long$value==1),]),
                                            nrow(scaa.long[scaa.long$geiStain=="gpos25" &   (scaa.long$value==1),]),
                                            nrow(scaa.long[scaa.long$geiStain=="gpos50" &   (scaa.long$value==1),]),
                                            nrow(scaa.long[scaa.long$geiStain=="gpos75" &   (scaa.long$value==1),]),
                                            nrow(scaa.long[scaa.long$geiStain=="gpos100" &   (scaa.long$value==1),])),
                NonSCAA_peak= c(nrow(scaa.long[scaa.long$geiStain=="gneg" &   (scaa.long$value==0),]),
                                                nrow(scaa.long[scaa.long$geiStain=="gvar" & (scaa.long$value==0),]),
                                                nrow(scaa.long[scaa.long$geiStain=="acen" & (scaa.long$value==0),]),
                                                nrow(scaa.long[scaa.long$geiStain=="gpos25" & (scaa.long$value==0),]),
                                                nrow(scaa.long[scaa.long$geiStain=="gpos50" & (scaa.long$value==0),]),
                                                nrow(scaa.long[scaa.long$geiStain=="gpos75" & (scaa.long$value==0),]),
                                                nrow(scaa.long[scaa.long$geiStain=="gpos100" & (scaa.long$value==0),])))
rownames(x) <- c("gneg","gvar","acen","gpos25","gpos50","gpos75","gpos100")
x$Stain <- rownames(x)
x <- gather(x, key = "key", value="Count", -Stain)

jpeg('tempfig.jpeg', width = (3*37.795*6), height = (3*37.795*4.19))
ggplot(x, aes(y=Count, x=Stain, fill= key)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label=Count), hjust = 0.5, position = position_dodge(width = 1), size=5) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 900000), breaks = c(0,400000,800000)) +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  theme(
    legend.text = element_text(size=20), legend.title = element_blank(), legend.position = "top",
    axis.title.y = element_blank())
dev.off()

# across patients cytoband
unique(scaa.long$geiStain)
x <- data.frame(SCAA= c(nrow(scaa.long[scaa.long$geiStain=="gneg" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gvar" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="acen" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos25" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos50" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos75" & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos100" & scaa.long$value==1,])),
                notSCAA= c(nrow(scaa.long[scaa.long$geiStain=="gneg" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="gvar" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="acen" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos25" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos50" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos75" & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos100" & scaa.long$value==0,])))
rownames(x) <- c("gneg","gvar","acen","gpos25","gpos50","gpos75","gpos100")

chisq <- chisq.test(x)
enrich <- chisq$observed/chisq$expected

jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
corrplot::corrplot(chisq$residuals, is.cor = FALSE)
dev.off()


# but per patient cytoband analysis shows high variability in data
output <- list()
corr <- list()
for (i in 1:length(unique(scaa.long$patient))) {
  patient <- unique(scaa.long$patient)[i]
  x <- scaa.long[which(scaa.long$patient == patient),]
  x <- data.frame(SCAA= c(nrow(x[x$geiStain=="gneg" & x$value==1,]),
                          nrow(x[x$geiStain=="gvar" & x$value==1,]),
                          nrow(x[x$geiStain=="acen" & x$value==1,]),
                          nrow(x[x$geiStain=="gpos25" & x$value==1,]),
                          nrow(x[x$geiStain=="gpos50" & x$value==1,]),
                          nrow(x[x$geiStain=="gpos75" & x$value==1,]),
                          nrow(x[x$geiStain=="gpos100" & x$value==1,])),
                  notSCAA= c(nrow(x[x$geiStain=="gneg" & x$value==0,]),
                             nrow(x[x$geiStain=="gvar" & x$value==0,]),
                             nrow(x[x$geiStain=="acen" & x$value==0,]),
                             nrow(x[x$geiStain=="gpos25" & x$value==0,]),
                             nrow(x[x$geiStain=="gpos50" & x$value==0,]),
                             nrow(x[x$geiStain=="gpos75" & x$value==0,]),
                             nrow(x[x$geiStain=="gpos100" & x$value==0,])))
  rownames(x) <- c("gneg","gvar","acen","gpos25","gpos50","gpos75","gpos100")
  
  chisq <- chisq.test(x)
  enrich <- chisq$observed/chisq$expected
  
  output[[i]] <- setNames(c(patient, chisq$p.value, enrich[,1]),
                          c("patient", "chiPvalue", names(enrich[,1])))
}

gstainPP <- data.frame(do.call(rbind, output))
gstainPP[-1] <- data.frame(apply(gstainPP[-1], 2, function(x) as.numeric(as.character(x))))
gstainPP$FDR <- p.adjust(gstainPP$chiPvalue, "fdr")
gstainPP <- gstainPP[gstainPP$FDR<=0.05,]
gstainPP <- gstainPP[,-c(2,10)]
gstainPP <- gather(gstainPP, key="stain", value="value", -patient)
gstainPP$value <- round(gstainPP$value, 2)
gstainPP$group <- ifelse(gstainPP$value==0,"0",ifelse(gstainPP$value==1,"1",ifelse(gstainPP$value<1,"less",ifelse(gstainPP$value>3,"more3",ifelse(gstainPP$value>2,"more2",ifelse(gstainPP$value>1,"more1",NA))))))

jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*7))
ggplot(gstainPP, aes(y=patient, x=stain, fill= group)) + 
  geom_tile(color="black") +
  scale_fill_manual(values = c("#FFFFFF","#CCCCCC", alpha("#3333FF",0.6),alpha("#CC0066",0.5),alpha("#CC0066",0.7),alpha("#CC0066",0.9))) +
  geom_text(aes(label=round(value,2))) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
dev.off()

scaa.long$geiStain <- factor(scaa.long$geiStain, levels = c("gneg","gpos100","acen","gpos25","gpos50","gpos75","no data","gvar"))
model <- glm(value ~ geiStain, data=scaa.long[which(scaa.long$geiStain!="no data"),], family = "binomial")
summary(model)

# where are the gvar peaks?
x <- scaa.long[which(scaa.long$geiStain=="gvar"),]
x <- na.omit(x)
x <- x[,c(1:7,23,25)]
x <- x[!duplicated(x),]
x$name <- paste(x$chr, x$bandName,sep="")

jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
ggplot(x, aes(x=(dist.telo)/1000000, fill=name)) + 
  geom_histogram(alpha=0.5, position="identity", color="black", 
                 bins = 50) +
  ylab("Count") +
  xlab("Distance to telomere (Mb)") +
  facet_grid(rows = vars(name)) +
  scale_x_continuous(breaks = c(10,25,50,75,100)) +
  scale_y_continuous(breaks = c(0,20,40)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  theme(legend.text = element_text(size=24),
        legend.title = element_blank(), legend.position = "none")
dev.off()

jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
ggplot(x, aes(x=(dist.centro)/1000000, fill=name)) + 
  geom_histogram(alpha=0.5, position="identity", color="black", 
                 bins = 50) +
  ylab("Count") +
  xlab("Distance to centromere (Mb)") +
  facet_grid(rows = vars(name)) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  scale_y_continuous(breaks = c(0,4,8,12)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  theme(legend.text = element_text(size=24),
        legend.title = element_blank(), legend.position = "none")
dev.off()

# what specific bands hold the peaks
unique(paste(x$chr,x$bandName,sep="."))


# SCAA ENRICHMENT BY TELOMERE DISTANCE =====

my_comparisons <- list( c("0-1", "1-25"), c("0-1", "25-50"), c("0-1", "50-75"), c("0-1", "75-100"), c("0-1", "100-150"),
                        c("1-25", "25-50"), c("1-25", "50-75"), c("1-25", "75-100"), c("1-25", "100-150"),
                        c("25-50", "50-75"), c("25-50", "75-100"), c("25-50", "100-150"),
                        c("50-75", "75-100"), c("50-75", "100-150"),
                        c("75-100", "100-150"))

jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
ggplot2::ggplot(data=scaa.long, aes(x=as.factor(value), y=(dist.telo), fill=as.factor(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to telomere (bp)") +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 125000000, label.x=1.4, size=7, method = "wilcox.test") +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
         axis.title.x = element_blank())
dev.off()

# by cna
jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
ggplot2::ggplot(data=scaa.long[!is.na(scaa.long$cna),], aes(x=as.factor(cna), y=log(dist.telo), fill=as.factor(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to telomere (bp)") +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 125000000, label.x=1.4, size=7, method = "wilcox.test") +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank())
dev.off()

wilcox.test(scaa.long$dist.telo[which(scaa.long$value==0 & scaa.long$cna=="diploid")], scaa.long$dist.telo[(scaa.long$value==1 & scaa.long$cna=="diploid")], 
            alternative = "greater")

model <- glm(value ~ dist.telo + cna + segSize, data=scaa.long, family = "binomial")
summary(model)

# by size
my_comparisons <- list( c("0-1", "1-25"), c("0-1", "25-50"), c("0-1", "50-75"), c("0-1", "75-100"), c("0-1", "100-150"),
                        c("1-25", "25-50"), c("1-25", "50-75"), c("1-25", "75-100"), c("1-25", "100-150"),
                        c("25-50", "50-75"), c("25-50", "75-100"), c("25-50", "100-150"),
                        c("50-75", "75-100"), c("50-75", "100-150"),
                        c("75-100", "100-150"))

jpeg('tempfig.jpeg', width = (3*37.795*10), height = (3*37.795*5))
ggplot2::ggplot(data=scaa.long, aes(x=as.factor(sizeGroup), y=(dist.telo), fill=as.factor(value))) +
  geom_split_violin(width=0.7, alpha = .4, trim = FALSE, scale = "width") +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to telomere (bp)") +
  xlab("Segment size (Mb)") +
  scale_fill_brewer(palette = "Dark2", labels=c("not SCAA","SCAA")) +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 150000000, label.x=1.4, size=7, method = "wilcox.test") +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
dev.off()



# SCAA ENRICHMENT BY CENTROMERE DISTANCE =====

jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
ggplot2::ggplot(data=scaa.long, aes(x=as.factor(value), y=(dist.centro), fill=as.factor(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to centromere (bp)") +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 125000000, label.x=1.4, size=7, method = "wilcox.test") +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank())
dev.off()

t.test(scaa.long$dist.centro[scaa.long$value==0], scaa.long$dist.centro[scaa.long$value==1])

model <- glm(value ~ dist.centro, data=scaa.long, family = "binomial")
summary(model)

# by cna
jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
ggplot2::ggplot(data=scaa.long[!is.na(scaa.long$cna),], aes(x=as.factor(cna), y=(dist.centro), fill=as.factor(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to centromere (bp)") +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 150000000, label.x=1.4, size=7, method = "wilcox.test") +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank())
dev.off()

wilcox.test(scaa.long$dist.centro[which(scaa.long$value==0 & scaa.long$cna=="diploid")], scaa.long$dist.centro[(scaa.long$value==1 & scaa.long$cna=="diploid")], 
            alternative = "greater")

# by size
jpeg('tempfig.jpeg', width = (3*37.795*10), height = (3*37.795*5))
ggplot2::ggplot(data=scaa.long, aes(x=as.factor(sizeGroup), y=(dist.centro), fill=as.factor(value))) +
  geom_split_violin(width=0.7, alpha = .4, trim = FALSE, scale = "width") +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to centromere (bp)") +
  xlab("Segment size (Mb)") +
  scale_fill_brewer(palette = "Dark2", labels=c("not SCAA","SCAA")) +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 150000000, label.x=1.4, size=7, method = "wilcox.test") +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
dev.off()

# SCAA ENRICHMENT BY BP DISTANCE =====

jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
ggplot2::ggplot(data=scaa.long, aes(x=as.factor(value), y=(distancetoBP), fill=as.factor(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to BP (bp)") +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 125000000, label.x=1.4, size=7, method = "wilcox.test") +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank())
dev.off()

t.test(scaa.long$dist.centro[scaa.long$value==0], scaa.long$dist.centro[scaa.long$value==1])

model <- glm(value ~ dist.centro, data=scaa.long, family = "binomial")
summary(model)

# by cna
jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
ggplot2::ggplot(data=scaa.long[!is.na(scaa.long$cna),], aes(x=as.factor(cna), y=(distancetoBP), fill=as.factor(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to centromere (bp)") +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 150000000, label.x=1.4, size=7, method = "wilcox.test") +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank())
dev.off()

wilcox.test(scaa.long$distancetoBP[which(scaa.long$value==0 & scaa.long$cna=="diploid")], scaa.long$distancetoBP[(scaa.long$value==1 & scaa.long$cna=="diploid")], 
            alternative = "greater")

# SAVE =====
saveRDS(scaaPP, "~/Documents/SCAA/Data/scaaPP.rds")
saveRDS(scaa.long, "~/Documents/SCAA/Data/scaa.long.rds")
saveRDS(scaa.bins, "~/Documents/SCAA/Data/scaa.bins.rds")

scaaPP <- readRDS("~/Documents/SCAA/Data/scaaPP.rds")
scaa.long <- readRDS("~/Documents/SCAA/Data/scaa.long.rds")
scaa.bins <- readRDS("~/Documents/SCAA/Data/scaa.bins.rds")
