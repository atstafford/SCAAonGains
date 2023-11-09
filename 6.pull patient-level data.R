
# DATA =====
peakAnn_all <- readRDS("~/Documents/SCAA/Data/peakAnn_all.rds")
promoter <- peakAnn_all$peak[which(peakAnn_all$Promoter==TRUE)]
enhancer <- peakAnn_all$peak[which(!is.na(peakAnn_all$elementType))]

patients <- names(dwgs.perpatient)
scaaPP <- list()
scaaPP_promo <- list()
scaaPP_notPromo <- list()
for (i in 1:length(patients)) {
  wd <- scaa.long[scaa.long$patient==patients[i],]
  scaaPP[[i]] <- sum(wd$value==1, na.rm = T) / nrow(wd)
  wd2 <- wd[which(wd$peak %in% enhancer),]
  scaaPP_promo[[i]] <- sum(wd2$value==1, na.rm = T) / nrow(wd2)
  wd3 <- wd[which(wd$peak %!in% enhancer),]
  scaaPP_notPromo[[i]] <- sum(wd3$value==1, na.rm = T) / nrow(wd3)
}
scaaPP <- data.frame(patient=patients, propSCAA=unlist(scaaPP), propSCAA.promo=unlist(scaaPP_promo), propSCAA.notPromo=unlist(scaaPP_notPromo))
scaaPP$meanPloidy <- unlist(ploidy.perpatient[match(scaaPP$patient, ploidy.perpatient$patient), 2])
scaaPP$WGD <- ifelse(scaaPP$meanPloidy>2.5, ">2.5", "<2.5")
scaaPP$type <- ifelse(scaaPP$patient %in% c("C552","C548","C518","C562","C516","C536"), "MSI", "MSS")

scaaPP <- scaaPP[which(scaaPP$patient %in% substr(colnames(scaa), 1,4)),]

# calculate PGA per sample and average for each patient
output <- list()
for (i in 1:length(dwgs.perpatient) ) {
  wd <- dwgs.perpatient[[i]]
  wd[,-c(1:3)] <- round(wd[,-c(1:3)])
  weights <- wd$stop-wd$start
  weights <- weights / sum(weights)
  
  if (ncol(wd)==4) {
    y <- weights * ifelse(wd$C527.B1_G7!=2, 1, 0)
    output[[i]] <- mean(y)
  } else {
    y <- apply(wd[,-c(1:3)], 2, function(x) {
      weights * ifelse(x!=2, 1, 0)
    })
    output[[i]] <- mean(colSums(y, na.rm = T))
  }
  
}

pga <- data.frame(patient=names(dwgs.perpatient), pga=unlist(output))

# calculate PGA per sample and average for each patient (with recenter)
output <- list()
for (i in 1:length(dwgs.ploidyRecentre.perpatient) ) {
  wd <- dwgs.ploidyRecentre.perpatient[[i]]
  wd[,-c(1:3)] <- round(wd[,-c(1:3)])
  weights <- wd$stop-wd$start
  weights <- weights / sum(weights)
  
  if (ncol(wd)==4) {
    y <- weights * ifelse(wd$C527.B1_G7!=2, 1, 0)
    output[[i]] <- mean(y)
  } else {
    y <- apply(wd[,-c(1:3)], 2, function(x) {
      weights * ifelse(x!=2, 1, 0)
    })
    output[[i]] <- mean(colSums(y, na.rm = T))
  }
}

pga_recentre <- data.frame(patient=names(dwgs.ploidyRecentre.perpatient), pga_recentre=unlist(output))

# calculate number of segmewnts pp
output <- list()
outputc <- list()
for (i in 1:length(dwgs.perpatient) ) {
  wd <- dwgs.perpatient[[i]]
  x <- apply(wd[-c(1:3)],2, function(x) {length((rle(x))[[1]])})
  output[[i]] <- mean(x)
  outputc[[i]] <- length(unique(wd$chromosome))
}

pga$n.segments <- unlist(output)
pga$n.chr <- unlist(outputc)

scaaPP$pga <- pga[match(scaaPP$patient, pga$patient),2]
scaaPP$propSCAA[is.na(scaaPP$propSCAA)] <- 0
scaaPP$n.segments <- pga[match(scaaPP$patient, pga$patient),3]
scaaPP$n.chr <- pga[match(scaaPP$patient, pga$patient),4]
scaaPP$remove <- ifelse(scaaPP$patient %in% c("C516","C542","C543"), TRUE, FALSE)
scaaPP$pga_recentre <- pga_recentre[match(scaaPP$patient, pga_recentre$patient),2]

# save
saveRDS(scaaPP, "~/Documents/SCAA/Data/scaaPP.rds")
