# PREREQUISITIS: load section 1

# LOAD ATAC AND PURITY/PLOIDY DATA --------------------------------------------------------------
atac_raw <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/atacseq_counts_per_peak.rds")
colnames(atac_raw) <- sub("_", ".", substr(colnames(atac_raw),7,nchar(colnames(atac_raw))-3) ) # correct format of IDs (patient.region_gland)

atac_pur <- readxl::read_excel("~/Documents/SCAA/Data/EPICC/Table_S2.xlsx") # s2 atac seq. app purity is est from atac data. not great corr with est from muts
atac_pur$Sample <- sub("_", ".", substr(atac_pur$sample_barcode,7,nchar(atac_pur$sample_barcode)-3) )

wgs_pur <- readxl::read_excel("~/Documents/SCAA/Data/EPICC/Table_S3.xlsx") # s3 wgs purity from seq (ignore). true est from muts is best for wgs samples
wgs_pur <- wgs_pur[which(wgs_pur$tissue_type=="cancer"), ]
wgs_pur$Sample <- sub("_", ".", substr(wgs_pur$sample_barcode,7,nchar(wgs_pur$sample_barcode)-3) )

# KEEP SAMPLES THAT HAVE BOTH ATAC AND DWGS DATA ----------------------------------------------------------------------------------------------------------------------------------

keep <- intersect(names(dwgs), colnames(atac_raw))
keep <- unique(keep)

dwgs_sub <- dwgs[keep]
atac_sub <- atac_raw[ ,keep] # same order

# ALIGN BINS BETWEEN ATAC AND WGS DATA -----------------------------------------------------------------------------------------------------------------------------------------------
atacbins <- rownames(atac_sub)
atacbins <- data.frame(bin=as.numeric(c(1:length(atacbins))),
                       chr=(substr(sub("\\:.*", "", atacbins),4,nchar(sub("\\:.*", "", atacbins)))),
                       start=as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", atacbins)),
                       stop=as.numeric(sub(".*-", "", atacbins)) 
)
atacbins <- na.omit(atacbins)

binmatch_dwgs<- list()
for (i in 1:length(dwgs_sub)) {
  print(i)
  wd <- dwgs_sub[[i]]
  wgsbins <- data.frame(bin=as.numeric(c(1:nrow(wd))),
                        chr=(substr(wd$chromosome, 4, nchar(wd$chromosome))),
                        start=as.numeric(wd$start),
                        stop=as.numeric(wd$stop) 
  )
  wgsbins <- na.omit(wgsbins)
  
  binmatch <- SimpleAlignBins(targetBins = wgsbins, incomingBins = atacbins)
  binmatch <- cbind(atacbins, data.frame(matched_wgs_bin = unlist(binmatch)))
  binmatch$matched_wgs_chr <- wgsbins[match(binmatch$matched_wgs_bin, wgsbins$bin), 2]
  binmatch$matched_wgs_start <- wgsbins[match(binmatch$matched_wgs_bin, wgsbins$bin), 3]
  binmatch$matched_wgs_stop <- wgsbins[match(binmatch$matched_wgs_bin, wgsbins$bin), 4]
  binmatch_dwgs[[i]] <- binmatch
}

# put atac data into dataframe per sample 
atac_sub.list <- list()
for (i in 1:ncol(atac_sub)) {
  atac_sub.list[[i]] <- atac_sub[,i]
  atac_sub.list[[i]] <- data.frame(atac_sub.list[[i]])
  colnames(atac_sub.list[[i]]) <- colnames(atac_sub)[i]
}

# ASSEMBLY MATCHING PURITY PLOIDY LIST FOR EACH ATAC SAMPLE -----------------------------------------------------------------------------------------------------------------------------------------------

# assemble purity and ploidy lists
purity <- list()
ploidy <- list()

for ( i in 1:length(atac_sub.list) ) {
  print(i)
  wd <- atac_sub.list[[i]]
  sample <- colnames(wd)
  
  # take apparent_purity from atac data. Not as good as from muts in wgs, but best practice 
  if ( sample %in% atac_pur$Sample) {
    purity[[i]] <- as.numeric(atac_pur[match(sample, atac_pur$Sample), 12])
  }
  else {
    purity[[i]] <- NA
  }
  
  # take average of ploidies associated with region from wgs
  if ( substr(sample,1,nchar(sample)-3) %in% substr(wgs_pur$Sample,1,nchar(wgs_pur$Sample)-3) ) {
    x <- wgs_pur$ploidy[which(substr(wgs_pur$Sample,1,nchar(wgs_pur$Sample)-3) == substr(sample,1,nchar(sample)-3))]
    ploidy[[i]] <- as.numeric(mean(x))
  }
  else {
    ploidy[[i]] <- NA
  }
}

# COPY NUMBER CORRECTION -----------------------------------------------------------------------------------------------------------------------------------------------
# correct for CN: S = Sn * [2(1 - purity) + CNA*purity]/[2(1 - purity) + ploidy*purity]
cn.correction <- list()
for(s in 1:length(atac_sub.list)) { # for a given sample
  corrected <- list()
  print(s)
  for (b in 1:nrow(atac_sub.list[[s]])) { # for a given bin
    
    # check we have purity, ploidy, and CNA data
    if ( is.na(purity[[s]]) | is.na(ploidy[[s]]) | is.na(binmatch_dwgs[[s]]$matched_wgs_chr[b]) ) {
      corrected[[b]] <- NA
      next
    }
    
    # pull cn for that sample/bin
    wgsbin <- as.numeric(binmatch_dwgs[[s]]$matched_wgs_bin[b])
    atacsample <- substr(colnames(atac_sub.list[[s]]), 1, nchar(colnames(atac_sub.list[[s]]))-4)
    wgscol <- which(atacsample == substr(names(dwgs_sub), 1, 6))
    
    if (length(wgscol)==1) {
      CNA <- dwgs_sub[[wgscol]]$CNt[wgsbin]
    } else {
      CNA <- list()
      for(i in 1:length(wgscol)) {
        CNA[[i]] <- dwgs_sub[[wgscol[i]]]$CNt[wgsbin]
      }
      CNA <- mean(unlist(CNA))
    }
    
    A <- atac_sub.list[[s]][b,]
    corrected[[b]] <- A * (  ( (2*(1-purity[[s]])) + (CNA*purity[[s]]) ) / ( (2*(1-purity[[s]])) + (ploidy[[s]]*purity[[s]]) )  )
  }
  cn.correction[[s]] <- matrix(unlist(corrected))
}

names(cn.correction) <- keep

# add back in peak names and remove missing bins
for (i in 1:length(cn.correction)) {
  rownames(cn.correction[[i]]) <- rownames(atac_sub.list[[i]])
  #cn.correction[[i]] <- na.omit(cn.correction[[i]])
}
cn.correction <- cn.correction[lapply(cn.correction,length)>0]

atac_cn.deepWGSnorm <- do.call(cbind, cn.correction)
colnames(atac_cn.deepWGSnorm) <- keep

# remove samples with no data
atac_cn.deepWGSnorm <- atac_cn.deepWGSnorm[,colSums(is.na(atac_cn.deepWGSnorm)) < nrow(atac_cn.deepWGSnorm)]

# some samples were missing dWGS data so peaks couldnt be normalised. DESeq2 cant handle missing values so these peaks need to be remove across all samples
atac_cn.deepWGSnorm <- atac_cn.deepWGSnorm[rowSums(is.na(atac_cn.deepWGSnorm)) == 0,]

# SAVE -----------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(atac_raw, "~/Documents/SCAA/Data/atac_raw.rds")
saveRDS(atac_pur, "~/Documents/SCAA/Data/atac_pur.rds")
saveRDS(wgs_pur, "~/Documents/SCAA/Data/wgs_pur.rds")
saveRDS(binmatch_dwgs, "~/Documents/SCAA/Data/binmatch_dwgs.rds")
saveRDS(cn.correction, "~/Documents/SCAA/Data/cn.correction.rds")
saveRDS(atac_cn.deepWGSnorm, "~/Documents/SCAA/Data/atac_cn.deepWGSnorm.rds")



