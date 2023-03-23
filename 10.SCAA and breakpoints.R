# PREREQUISITIS: load section 7 data 

# FREQ/DISTANCE OF SCAAS AROUND BREAKPOINTS  =====

peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")] # keep only carcinomas ('pure')
#scaa <- scaa[which(rownames(scaa) %in% peakAnn$peak[which(peakAnn$Promoter == "TRUE" | !is.na(peakAnn$elementType))]),]
#scaa <- scaa[which(rowSums(scaa)>=(ncol(scaa)*0.2)),]
scaa.location <- data.frame(chr=sub("\\:.*", "", rownames(scaa)), start=sub(".*[:]([^.]+)[-].*", "\\1", rownames(scaa)), stop=sub(".*-", "", rownames(scaa)))

nearest_scaa <- list()
distance <- list()
n.scaa_within1Mb <- list()
for (i in 1:nrow(breakpoints_dWGS)) {
  print(paste(i,"/",nrow(breakpoints_dWGS)))
  patient <- breakpoints_dWGS$patient[i]
  chr <- sub("\\:.*", "", breakpoints_dWGS$location[i])
  
  scaa.wd <- cbind(scaa.location, scaa[which(substr(colnames(scaa),1,4)==patient)]) # pull SCAAs for the patient
  scaa.wd <- scaa.wd[which(sub("\\:.*", "",rownames(scaa))==chr), ] # pull SCAAs for the chr
  scaa.wd <- scaa.wd[which(scaa.wd[4]==1),]
  scaa.wd$mid.peak <- (as.numeric(scaa.wd$stop)-as.numeric(scaa.wd$start)) + as.numeric(scaa.wd$start)
  
  nearest_scaa[[i]] <- scaa.wd[which.min(abs(breakpoints_dWGS$break_base[i]-scaa.wd$mid.peak)),-4]
  
  if(nrow(nearest_scaa[[i]])<1) { # no SCAA on chromosome
    nearest_scaa[[i]] <- data.frame(chr=NA,start=NA,stop=NA,mid.peak=NA)
    distance[[i]] <- NA
    n.scaa_within1Mb[[i]] <- NA
  }
  distance[[i]] <- abs(breakpoints_dWGS$break_base[i]-nearest_scaa[[i]]$mid.peak)
  n.scaa_within1Mb[[i]] <- nrow(scaa.wd[which(abs(breakpoints_dWGS$break_base[i]-scaa.wd$mid.peak)<=1000000),-4])
}

scaaPerBreakpoint <- breakpoints_dWGS 
scaaPerBreakpoint$nearest_scaa <- do.call("rbind",nearest_scaa) 
scaaPerBreakpoint$distance <- unlist(distance)
scaaPerBreakpoint$n.scaa_within1Mb <- unlist(n.scaa_within1Mb)

hist(scaaPerBreakpoint$n.scaa_within1Mb, breaks = 200)
hist(scaaPerBreakpoint$distance, breaks = 2000, xlim = c(0,100000000))

# GENERATE RANDOM POINTS =====

# pull min/max start/stop per chromosome
chrs <- unique(dat[[1]]$chromosome)

output1 <- list()
output2 <- list()
for (i in 1:length(dat)) {
  wd1 <- dat[[i]]
  wd1 <- wd1[!is.na(wd1$CNt),]
  
  start <- list()
  stop <- list()
  for (j in 1:length(chrs)) {
    chr <- chrs[j]
    wd2 <- wd1[which(wd1$chromosome == chr),]
    
    if (nrow(wd2)==0) {
      start[[j]] <- NA
      stop[[j]] <- NA
    } else {
      start[[j]] <- min(wd2$start.pos)
      stop[[j]] <- max(wd2$end.pos)
    }
  }
  output1[[i]] <- data.frame("minStart"=unlist(start))
  output2[[i]] <- data.frame("maxStop"=unlist(stop))
}

Chrstarts <- as.matrix(do.call("cbind", output1))
Chrstops <- as.matrix(do.call("cbind", output2))

ChrEnds <- data.frame("start"=rowMaxs(Chrstarts, na.rm = T), "stop"=rowMins(Chrstops, na.rm = T)) # choose latest start and earliest stop so random points always generated within appropriate range

# sample 10,000x per chromosome
output <- list()
for (i in 1:nrow(ChrEnds)) {
  output[[i]] <- paste(i, sample(ChrEnds$start[i]:ChrEnds$stop[i], 10000), sep=":")
}
randomPoints <- unlist(output)

# FREQ/DISTANCE OF SCAAS AROUND RANDOM POINTS =====
randomPoints_sample <- data.frame("location"=sample(randomPoints, size=nrow(breakpoints_dWGS), replace = FALSE))
randomPoints_sample$break_base <- sub(".*:", "", randomPoints_sample$location)
randomPoints_sample$patient_sub <- print(sample(unique(substr(names(dat),1,4)), nrow(randomPoints_sample),replace=TRUE))

nearest_scaa <- list()
distance <- list()
n.scaa_within1Mb <- list()
for (i in 1:nrow(randomPoints_sample)) {
  print(paste(i,"/",nrow(randomPoints_sample)))
  patient <- randomPoints_sample$patient_sub[i]
  chr <- paste("chr",sub("\\:.*", "", randomPoints_sample$location[i]), sep = "")
  
  scaa.wd <- cbind(scaa.location, scaa[which(substr(colnames(scaa),1,4)==patient)])
  scaa.wd <- scaa.wd[which(sub("\\:.*", "",rownames(scaa))==chr), ]
  scaa.wd <- scaa.wd[which(scaa.wd[4]==1),]
  scaa.wd$mid.peak <- (as.numeric(scaa.wd$stop)-as.numeric(scaa.wd$start)) + as.numeric(scaa.wd$start)
  
  nearest_scaa[[i]] <- scaa.wd[which.min(abs(as.numeric(randomPoints_sample$break_base[i])-scaa.wd$mid.peak)),-4]
  if(nrow(nearest_scaa[[i]])<1) {
    nearest_scaa[[i]] <- data.frame(chr=NA,start=NA,stop=NA,mid.peak=NA)
  }
  distance[[i]] <- abs(as.numeric(randomPoints_sample$break_base[i])-nearest_scaa[[i]]$mid.peak)
  n.scaa_within1Mb[[i]] <- nrow(scaa.wd[which(abs(as.numeric(randomPoints_sample$break_base[i])-scaa.wd$mid.peak)<=1000000),-4])
}

scaaPerRandomPoints <- randomPoints_sample
scaaPerRandomPoints$nearest_scaa <- do.call("rbind",nearest_scaa)
scaaPerRandomPoints$distance <- unlist(distance)
scaaPerRandomPoints$n.scaa_within1Mb <- unlist(n.scaa_within1Mb)

hist(scaaPerRandomPoints$n.scaa_within1Mb, breaks = 200)
hist(scaaPerRandomPoints$distance, breaks = 2000, xlim = c(0,100000000))

# COMPARISON BASED ON SCAA WITHIN 1MB =====
# comparing with random
ks.test(scaaPerBreakpoint$n.scaa_within1Mb,scaaPerRandomPoints$n.scaa_within1Mb)
plot(ecdf(x=scaaPerBreakpoint$n.scaa_within1Mb))
lines(ecdf(x=scaaPerRandomPoints$n.scaa_within1Mb), col=2)

data <- rbind(data.frame(Log_value=log(scaaPerBreakpoint$n.scaa_within1Mb),group="Breakpoints"),
              data.frame(Log_value=log(scaaPerRandomPoints$n.scaa_within1Mb),group="Not Breakpoints"))

ggplot(data, aes(x=group, y=Log_value)) +
  geom_boxplot(width=0.4) +
  geom_violin(aes(fill=group)) +
  geom_jitter(shape=16, position=position_jitter(0.4)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(fill=NA))

# investigate bps in super scaa rich areas
t.test(scaaPerBreakpoint$n.scaa_within1Mb, scaaPerRandomPoints$n.scaa_within1Mb)

# COMPARISON BASED ON DISTANCE TO NEAREST SCAA =====
# comparing with random
ks.test(scaaPerBreakpoint$distance,scaaPerRandomPoints$distance)

plot(ecdf(x=scaaPerBreakpoint$distance))
lines(ecdf(x=scaaPerRandomPoints$distance), col=2)

data <- rbind(data.frame(Log_value=log(scaaPerBreakpoint$distance),group="Breakpoints"),
              data.frame(Log_value=log(scaaPerRandomPoints$distance),group="Not Breakpoints"))

ggplot(data, aes(x=group, y=Log_value)) +
  geom_violin(aes(fill=group)) +
  geom_jitter(shape=16, position=position_jitter(0.2), fill=alpha("black",0.1), col=alpha("black",0.1), size=1) +
  geom_boxplot(width=0.1,aes(fill=group)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(fill=NA))

# is there a min distance?
t.test(scaaPerBreakpoint$distance, scaaPerRandomPoints$distance)

# SAVE =====
saveRDS(scaaPerBreakpoint, "~/Documents/SCAA/Data/EPICC/SCAA/scaaPerBreakpoint.rds")
saveRDS(scaaPerRandomPoints, "~/Documents/SCAA/Data/EPICC/SCAA/scaaPerRandomPoints.rds")
saveRDS(scaaPerNonBreakpoint, "~/Documents/SCAA/Data/EPICC/SCAA/scaaPerNonBreakpoint.rds")
saveRDS(randomPoints, "~/Documents/SCAA/Data/EPICC/SCAA/randomPoints.rds")


scaaPerBreakpoint <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/scaaPerBreakpoint.rds")
scaaPerRandomPoints <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/scaaPerRandomPoints.rds")

