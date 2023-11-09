# PREREQUISITIS: load section 7 data 

# DATA =====
breakpoints_dWGS <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/breakpoints_dWGS.rds")
breakpoints_dWGS_highres10000 <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/breakpoints_dWGS_highres10000.rds")
randomPoints <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/randomPoints.rds")

# FREQ/DISTANCE OF SCAAS AROUND BREAKPOINTS  =====

peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")] # keep only carcinomas ('pure')
scaa.location <- data.frame(chr=sub("\\:.*", "", rownames(scaa)), start=sub(".*[:]([^.]+)[-].*", "\\1", rownames(scaa)), stop=sub(".*-", "", rownames(scaa)))

nearest_scaa <- list()
distance_scaa <- list()
n.scaa_within100kb <- list()
n.scaa_within250kb <- list()
n.scaa_within500kb <- list()
n.scaa_within750kb <- list()
n.scaa_within1Mb <- list()
n.scaa_within1.5Mb <- list()
n.scaa_within5Mb <- list()
n.scaa_within20Mb <- list()

nearest_peak <- list()
distance_peak <- list()
n.peaks_within100kb <- list()
n.peaks_within250kb <- list()
n.peaks_within500kb <- list()
n.peaks_within750kb <- list()
n.peaks_within1Mb <- list()
n.peaks_within1.5Mb <- list()
n.peaks_within5Mb <- list()
n.peaks_within20Mb <- list()

for (i in 1:nrow(breakpoints_dWGS_sub)) {
  print(paste(i,"/",nrow(breakpoints_dWGS_sub)))
  patient <- breakpoints_dWGS_sub$patient[i]
  chr <- sub("\\:.*", "", breakpoints_dWGS_sub$location[i])
  
  # find SCAAs
  scaa.wd <- cbind(scaa.location, scaa[which(substr(colnames(scaa),1,4)==patient)]) # pull SCAAs for the patient
  scaa.wd <- scaa.wd[which(sub("\\:.*", "",rownames(scaa))==chr), ] # pull SCAAs for the chr
  scaa.wd <- scaa.wd[which(scaa.wd[4]==1),]
  scaa.wd$mid.peak <- (as.numeric(scaa.wd$stop)-as.numeric(scaa.wd$start)) + as.numeric(scaa.wd$start)
  
  nearest_scaa[[i]] <- rownames(scaa.wd[which.min(abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak)),-4])
  
  if(length(nearest_scaa[[i]])<1) { # no SCAA on chromosome
    nearest_scaa[[i]] <- NA
    distance_scaa[[i]] <- NA
    n.scaa_within100kb[[i]] <- NA
    n.scaa_within250kb[[i]] <- NA
    n.scaa_within500kb[[i]] <- NA
    n.scaa_within750kb[[i]] <- NA
    n.scaa_within1Mb[[i]] <- NA
    n.scaa_within1.5Mb[[i]] <- NA
    n.scaa_within5Mb[[i]] <- NA
    n.scaa_within20Mb[[i]] <- NA
  }
  distance_scaa[[i]] <- abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak[which.min(abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak))])
  n.scaa_within100kb[[i]] <- nrow(scaa.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak)<=100000),-4])
  n.scaa_within250kb[[i]] <- nrow(scaa.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak)<=250000),-4])
  n.scaa_within500kb[[i]] <- nrow(scaa.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak)<=500000),-4])
  n.scaa_within750kb[[i]] <- nrow(scaa.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak)<=750000),-4])
  n.scaa_within1Mb[[i]] <- nrow(scaa.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak)<=1000000),-4])
  n.scaa_within1.5Mb[[i]] <- nrow(scaa.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak)<=1500000),-4])
  n.scaa_within5Mb[[i]] <- nrow(scaa.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak)<=5000000),-4])
  n.scaa_within20Mb[[i]] <- nrow(scaa.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-scaa.wd$mid.peak)<=20000000),-4])
  
  # find peaks
  peak.wd <- cbind(scaa.location, scaa[which(substr(colnames(scaa),1,4)==patient)])
  peak.wd <- peak.wd[which(sub("\\:.*", "",rownames(scaa))==chr), ] # pull peaks for the chr
  peak.wd$mid.peak <- (as.numeric(peak.wd$stop)-as.numeric(peak.wd$start)) + as.numeric(peak.wd$start)
  
  nearest_peak[[i]] <- rownames(peak.wd[which.min(abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak)),-4])
  
  if(length(nearest_peak[[i]])<1) { # no SCAA on chromosome
    nearest_peak[[i]] <- NA
    distance_peak[[i]] <- NA
    n.peaks_within1Mb[[i]] <- NA
  }
  distance_peak[[i]] <- abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak[which.min(abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak))])
  n.peaks_within100kb[[i]] <- nrow(peak.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak)<=100000),-4])
  n.peaks_within250kb[[i]] <- nrow(peak.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak)<=250000),-4])
  n.peaks_within500kb[[i]] <- nrow(peak.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak)<=500000),-4])
  n.peaks_within750kb[[i]] <- nrow(peak.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak)<=750000),-4])
  n.peaks_within1Mb[[i]] <- nrow(peak.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak)<=1000000),-4])
  n.peaks_within1.5Mb[[i]] <- nrow(peak.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak)<=1500000),-4])
  n.peaks_within5Mb[[i]] <- nrow(peak.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak)<=5000000),-4])
  n.peaks_within20Mb[[i]] <- nrow(peak.wd[which(abs(breakpoints_dWGS_sub$breakpoint[i]-peak.wd$mid.peak)<=20000000),-4])
}

# load into dataframe----
scaaPerBreakpoint <- breakpoints_dWGS_sub

scaaPerBreakpoint$nearest_scaa <- do.call("rbind",nearest_scaa) 
distance_scaa[lengths(distance_scaa) == 0] <- NA
scaaPerBreakpoint$distance_scaa <- unlist(distance_scaa)
scaaPerBreakpoint$n.scaa_within100kb <- unlist(n.scaa_within100kb)
scaaPerBreakpoint$n.scaa_within250kb <- unlist(n.scaa_within250kb)
scaaPerBreakpoint$n.scaa_within500kb <- unlist(n.scaa_within500kb)
scaaPerBreakpoint$n.scaa_within750kb <- unlist(n.scaa_within750kb)
scaaPerBreakpoint$n.scaa_within1Mb <- unlist(n.scaa_within1Mb)
scaaPerBreakpoint$n.scaa_within1.5Mb <- unlist(n.scaa_within1.5Mb)
scaaPerBreakpoint$n.scaa_within5Mb <- unlist(n.scaa_within5Mb)
scaaPerBreakpoint$n.scaa_within20Mb <- unlist(n.scaa_within20Mb)

scaaPerBreakpoint$nearest_peak <- do.call("rbind",nearest_peak) 
distance_peak[lengths(distance_peak) == 0] <- NA
scaaPerBreakpoint$distance_peak <- unlist(distance_peak)
scaaPerBreakpoint$n.peaks_within100kb <- unlist(n.peaks_within100kb)
scaaPerBreakpoint$n.peaks_within250kb <- unlist(n.peaks_within250kb)
scaaPerBreakpoint$n.peaks_within500kb <- unlist(n.peaks_within500kb)
scaaPerBreakpoint$n.peaks_within750kb <- unlist(n.peaks_within750kb)
scaaPerBreakpoint$n.peaks_within1Mb <- unlist(n.peaks_within1Mb)
scaaPerBreakpoint$n.peaks_within1.5Mb <- unlist(n.peaks_within1.5Mb)
scaaPerBreakpoint$n.peaks_within5Mb <- unlist(n.peaks_within5Mb)
scaaPerBreakpoint$n.peaks_within20Mb <- unlist(n.peaks_within20Mb)


scaaPerBreakpoint$propSCAA_100kb <- scaaPerBreakpoint$n.scaa_within100kb / scaaPerBreakpoint$n.peaks_within100kb
scaaPerBreakpoint$propSCAA_250kb <- scaaPerBreakpoint$n.scaa_within250kb / scaaPerBreakpoint$n.peaks_within250kb
scaaPerBreakpoint$propSCAA_500kb <- scaaPerBreakpoint$n.scaa_within500kb / scaaPerBreakpoint$n.peaks_within500kb
scaaPerBreakpoint$propSCAA_750kb <- scaaPerBreakpoint$n.scaa_within750kb / scaaPerBreakpoint$n.peaks_within750kb
scaaPerBreakpoint$propSCAA_1Mb <- scaaPerBreakpoint$n.scaa_within1Mb / scaaPerBreakpoint$n.peaks_within1Mb
scaaPerBreakpoint$propSCAA_1.5Mb <- scaaPerBreakpoint$n.scaa_within1.5Mb / scaaPerBreakpoint$n.peaks_within1.5Mb
scaaPerBreakpoint$propSCAA_5Mb <- scaaPerBreakpoint$n.scaa_within5Mb / scaaPerBreakpoint$n.peaks_within5Mb
scaaPerBreakpoint$propSCAA_20Mb <- scaaPerBreakpoint$n.scaa_within20Mb / scaaPerBreakpoint$n.peaks_within20Mb

scaaPerBreakpoint$propSCAA_pp <- scaaPP[match(scaaPerBreakpoint$patient, scaaPP$patient),2]

# FREQ/DISTANCE OF SCAAS AROUND RANDOM POINTS =====
set.seed(12345)
randomPoints_sample <- data.frame("location"=sample(randomPoints, size=nrow(scaaPerBreakpoint), replace = FALSE))
randomPoints_sample$break_base <- as.numeric(sub(".*:", "", randomPoints_sample$location))
set.seed(12345)
#randomPoints_sample$patient <- print(sample(unique(scaaPerBreakpoint$patient)[unique(scaaPerBreakpoint$patient) %!in% c("C516","C542","C543")], nrow(randomPoints_sample),replace=TRUE))
randomPoints_sample$patient <- print(sample(scaaPerBreakpoint$patient, nrow(randomPoints_sample),replace=FALSE))

nearest_scaa <- list()
distance_scaa <- list()
n.scaa_within100kb <- list()
n.scaa_within250kb <- list()
n.scaa_within500kb <- list()
n.scaa_within750kb <- list()
n.scaa_within1Mb <- list()
n.scaa_within1.5Mb <- list()
n.scaa_within5Mb <- list()
n.scaa_within20Mb <- list()

nearest_peak <- list()
distance_peak <- list()
n.peaks_within100kb <- list()
n.peaks_within250kb <- list()
n.peaks_within500kb <- list()
n.peaks_within750kb <- list()
n.peaks_within1Mb <- list()
n.peaks_within1.5Mb <- list()
n.peaks_within5Mb <- list()
n.peaks_within20Mb <- list()

for (i in 1:nrow(randomPoints_sample)) {
  print(paste(i,"/",nrow(randomPoints_sample)))
  patient <- randomPoints_sample$patient[i]
  chr <- paste("chr",sub("\\:.*", "", randomPoints_sample$location[i]), sep = "")
  
  # find SCAAs
  scaa.wd <- cbind(scaa.location, scaa[which(substr(colnames(scaa),1,4)==patient)]) # pull SCAAs for the patient
  scaa.wd <- scaa.wd[which(sub("\\:.*", "",rownames(scaa))==chr), ] # pull SCAAs for the chr
  scaa.wd <- scaa.wd[which(scaa.wd[4]==1),]
  scaa.wd$mid.peak <- (as.numeric(scaa.wd$stop)-as.numeric(scaa.wd$start)) + as.numeric(scaa.wd$start)
  
  nearest_scaa[[i]] <- rownames(scaa.wd[which.min(abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak)),-4])
  
  if(length(nearest_scaa[[i]])<1) { # no SCAA on chromosome
    nearest_scaa[[i]] <- NA
    distance_scaa[[i]] <- NA
    n.scaa_within100kb[[i]] <- NA
    n.scaa_within250kb[[i]] <- NA
    n.scaa_within500kb[[i]] <- NA
    n.scaa_within750kb[[i]] <- NA
    n.scaa_within1Mb[[i]] <- NA
    n.scaa_within1.5Mb[[i]] <- NA
    n.scaa_within5Mb[[i]] <- NA
    n.scaa_within20Mb[[i]] <- NA
  }
  distance_scaa[[i]] <- abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak[which.min(abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak))])
  n.scaa_within100kb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak)<=100000),-4])
  n.scaa_within250kb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak)<=250000),-4])
  n.scaa_within500kb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak)<=500000),-4])
  n.scaa_within750kb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak)<=750000),-4])
  n.scaa_within1Mb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak)<=1000000),-4])
  n.scaa_within1.5Mb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak)<=1500000),-4])
  n.scaa_within5Mb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak)<=5000000),-4])
  n.scaa_within20Mb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_sample$break_base[i]-scaa.wd$mid.peak)<=20000000),-4])
  
  # find peaks
  peak.wd <- cbind(scaa.location, scaa[which(substr(colnames(scaa),1,4)==patient)])
  peak.wd <- peak.wd[which(sub("\\:.*", "",rownames(scaa))==chr), ] # pull peaks for the chr
  peak.wd$mid.peak <- (as.numeric(peak.wd$stop)-as.numeric(peak.wd$start)) + as.numeric(peak.wd$start)
  
  nearest_peak[[i]] <- rownames(peak.wd[which.min(abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak)),-4])
  
  if(length(nearest_peak[[i]])<1) { # no SCAA on chromosome
    nearest_peak[[i]] <- NA
    distance_peak[[i]] <- NA
    n.peaks_within1Mb[[i]] <- NA
  }
  distance_peak[[i]] <- abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak[which.min(abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak))])
  n.peaks_within100kb[[i]] <- nrow(peak.wd[which(abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak)<=100000),-4])
  n.peaks_within250kb[[i]] <- nrow(peak.wd[which(abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak)<=250000),-4])
  n.peaks_within500kb[[i]] <- nrow(peak.wd[which(abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak)<=500000),-4])
  n.peaks_within750kb[[i]] <- nrow(peak.wd[which(abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak)<=750000),-4])
  n.peaks_within1Mb[[i]] <- nrow(peak.wd[which(abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak)<=1000000),-4])
  n.peaks_within1.5Mb[[i]] <- nrow(peak.wd[which(abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak)<=1500000),-4])
  n.peaks_within5Mb[[i]] <- nrow(peak.wd[which(abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak)<=5000000),-4])
  n.peaks_within20Mb[[i]] <- nrow(peak.wd[which(abs(randomPoints_sample$break_base[i]-peak.wd$mid.peak)<=20000000),-4])
}

# load into dataframe----
scaaPerRandomPoints <- randomPoints_sample

scaaPerRandomPoints$nearest_scaa <- do.call("rbind",nearest_scaa) 
distance_scaa[lengths(distance_scaa) == 0] <- NA
scaaPerRandomPoints$distance_scaa <- unlist(distance_scaa)
scaaPerRandomPoints$n.scaa_within100kb <- unlist(n.scaa_within100kb)
scaaPerRandomPoints$n.scaa_within250kb <- unlist(n.scaa_within250kb)
scaaPerRandomPoints$n.scaa_within500kb <- unlist(n.scaa_within500kb)
scaaPerRandomPoints$n.scaa_within750kb <- unlist(n.scaa_within750kb)
scaaPerRandomPoints$n.scaa_within1Mb <- unlist(n.scaa_within1Mb)
scaaPerRandomPoints$n.scaa_within1.5Mb <- unlist(n.scaa_within1.5Mb)
scaaPerRandomPoints$n.scaa_within5Mb <- unlist(n.scaa_within5Mb)
scaaPerRandomPoints$n.scaa_within20Mb <- unlist(n.scaa_within20Mb)

scaaPerRandomPoints$nearest_peak <- do.call("rbind",nearest_peak) 
distance_peak[lengths(distance_peak) == 0] <- NA
scaaPerRandomPoints$distance_peak <- unlist(distance_peak)
scaaPerRandomPoints$n.peaks_within100kb <- unlist(n.peaks_within100kb)
scaaPerRandomPoints$n.peaks_within250kb <- unlist(n.peaks_within250kb)
scaaPerRandomPoints$n.peaks_within500kb <- unlist(n.peaks_within500kb)
scaaPerRandomPoints$n.peaks_within750kb <- unlist(n.peaks_within750kb)
scaaPerRandomPoints$n.peaks_within1Mb <- unlist(n.peaks_within1Mb)
scaaPerRandomPoints$n.peaks_within1.5Mb <- unlist(n.peaks_within1.5Mb)
scaaPerRandomPoints$n.peaks_within5Mb <- unlist(n.peaks_within5Mb)
scaaPerRandomPoints$n.peaks_within20Mb <- unlist(n.peaks_within20Mb)

scaaPerRandomPoints$propSCAA_100kb <- scaaPerRandomPoints$n.scaa_within100kb / scaaPerRandomPoints$n.peaks_within100kb
scaaPerRandomPoints$propSCAA_250kb <- scaaPerRandomPoints$n.scaa_within250kb / scaaPerRandomPoints$n.peaks_within250kb
scaaPerRandomPoints$propSCAA_500kb <- scaaPerRandomPoints$n.scaa_within500kb / scaaPerRandomPoints$n.peaks_within500kb
scaaPerRandomPoints$propSCAA_750kb <- scaaPerRandomPoints$n.scaa_within750kb / scaaPerRandomPoints$n.peaks_within750kb
scaaPerRandomPoints$propSCAA_1Mb <- scaaPerRandomPoints$n.scaa_within1Mb / scaaPerRandomPoints$n.peaks_within1Mb
scaaPerRandomPoints$propSCAA_1.5Mb <- scaaPerRandomPoints$n.scaa_within1.5Mb / scaaPerRandomPoints$n.peaks_within1.5Mb
scaaPerRandomPoints$propSCAA_5Mb <- scaaPerRandomPoints$n.scaa_within5Mb / scaaPerRandomPoints$n.peaks_within5Mb
scaaPerRandomPoints$propSCAA_20Mb <- scaaPerRandomPoints$n.scaa_within20Mb / scaaPerRandomPoints$n.peaks_within20Mb

scaaPerRandomPoints$propSCAA_pp <- scaaPP[match(scaaPerRandomPoints$patient, scaaPP$patient),2]

# stop here-----
# FREQ/DISTANCE OF SCAAS AROUND RANDOM POINTS 2Mb AWAY FROM BP =====
set.seed(1234)
randomPoints_2Mbaway_sample <- data.frame("location"=sample(randomPoints_2Mbaway, size=nrow(breakpoints_dWGS), replace = FALSE))
randomPoints_2Mbaway_sample$break_base <- as.numeric(sub(".*:", "", randomPoints_2Mbaway_sample$location))
set.seed(1234)
randomPoints_2Mbaway_sample$patient_sub <- print(sample(unique(scaaPerBreakpoint$patient)[unique(scaaPerBreakpoint$patient) %!in% c("C516","C542","C543")], nrow(randomPoints_sample),replace=TRUE))

nearest_scaa <- list()
distance_scaa <- list()
n.scaa_within1Mb <- list()
n.scaa_within25Mb <- list()
n.scaa_within500kb <- list()
nearest_peak <- list()
distance_peak <- list()
n.peaks_within1Mb <- list()
n.peaks_within25Mb <- list()
n.peaks_within500kb <- list()
for (i in 1:nrow(randomPoints_2Mbaway_sample)) {
  print(paste(i,"/",nrow(randomPoints_2Mbaway_sample)))
  patient <- randomPoints_2Mbaway_sample$patient_sub[i]
  chr <- paste("chr",sub("\\:.*", "", randomPoints_2Mbaway_sample$location[i]), sep = "")
  
  # find SCAAs
  scaa.wd <- cbind(scaa.location, scaa[which(substr(colnames(scaa),1,4)==patient)]) # pull SCAAs for the patient
  scaa.wd <- scaa.wd[which(sub("\\:.*", "",rownames(scaa))==chr), ] # pull SCAAs for the chr
  scaa.wd <- scaa.wd[which(scaa.wd[4]==1),]
  scaa.wd$mid.peak <- (as.numeric(scaa.wd$stop)-as.numeric(scaa.wd$start)) + as.numeric(scaa.wd$start)
  
  nearest_scaa[[i]] <- scaa.wd[which.min(abs(randomPoints_2Mbaway_sample$break_base[i]-scaa.wd$mid.peak)),-4]
  
  if(nrow(nearest_scaa[[i]])<1) { # no SCAA on chromosome
    nearest_scaa[[i]] <- data.frame(chr=NA,start=NA,stop=NA,mid.peak=NA)
    distance_scaa[[i]] <- NA
    n.scaa_within1Mb[[i]] <- NA
  }
  distance_scaa[[i]] <- abs(randomPoints_2Mbaway_sample$break_base[i]-nearest_scaa[[i]]$mid.peak)
  n.scaa_within1Mb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_2Mbaway_sample$break_base[i]-scaa.wd$mid.peak)<=1000000),-4])
  n.scaa_within25Mb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_2Mbaway_sample$break_base[i]-scaa.wd$mid.peak)<=25000000),-4])
  n.scaa_within500kb[[i]] <- nrow(scaa.wd[which(abs(randomPoints_2Mbaway_sample$break_base[i]-scaa.wd$mid.peak)<=500000),-4])
  
  # find peaks
  peak.wd <- cbind(scaa.location, scaa[which(substr(colnames(scaa),1,4)==patient)])
  peak.wd <- peak.wd[which(sub("\\:.*", "",rownames(scaa))==chr), ] # pull peaks for the chr
  peak.wd$mid.peak <- (as.numeric(peak.wd$stop)-as.numeric(peak.wd$start)) + as.numeric(peak.wd$start)
  
  nearest_peak[[i]] <- peak.wd[which.min(abs(randomPoints_2Mbaway_sample$break_base[i]-peak.wd$mid.peak)),-4]
  
  if(nrow(nearest_peak[[i]])<1) { # no SCAA on chromosome
    nearest_peak[[i]] <- data.frame(chr=NA,start=NA,stop=NA,mid.peak=NA)
    distance_peak[[i]] <- NA
    n.peaks_within1Mb[[i]] <- NA
  }
  distance_peak[[i]] <- abs(randomPoints_2Mbaway_sample$break_base[i]-nearest_peak[[i]]$mid.peak)
  n.peaks_within1Mb[[i]] <- nrow(peak.wd[which(abs(randomPoints_2Mbaway_sample$break_base[i]-peak.wd$mid.peak)<=1000000),-4])
  n.peaks_within25Mb[[i]] <- nrow(peak.wd[which(abs(randomPoints_2Mbaway_sample$break_base[i]-peak.wd$mid.peak)<=25000000),-4])
  n.peaks_within500kb[[i]] <- nrow(peak.wd[which(abs(randomPoints_2Mbaway_sample$break_base[i]-peak.wd$mid.peak)<=500000),-4])
  
}

scaaPerRandomPoints2Mbaway <- randomPoints_2Mbaway_sample

scaaPerRandomPoints2Mbaway$nearest_scaa <- do.call("rbind",nearest_scaa) 
scaaPerRandomPoints2Mbaway$distance_scaa <- unlist(distance_scaa)
scaaPerRandomPoints2Mbaway$n.scaa_within1Mb <- unlist(n.scaa_within1Mb)
scaaPerRandomPoints2Mbaway$n.scaa_within25Mb <- unlist(n.scaa_within25Mb)
scaaPerRandomPoints2Mbaway$n.scaa_within500kb <- unlist(n.scaa_within500kb)

scaaPerRandomPoints2Mbaway$nearest_peak <- do.call("rbind",nearest_peak) 
scaaPerRandomPoints2Mbaway$distance_peak <- unlist(distance_peak)
scaaPerRandomPoints2Mbaway$n.peaks_within1Mb <- unlist(n.peaks_within1Mb)
scaaPerRandomPoints2Mbaway$n.peaks_within25Mb <- unlist(n.peaks_within25Mb)
scaaPerRandomPoints2Mbaway$n.peaks_within500kb <- unlist(n.peaks_within500kb)

scaaPerRandomPoints$propSCAA_1Mb <- scaaPerRandomPoints$n.scaa_within1Mb / scaaPerRandomPoints$n.peaks_within1Mb
scaaPerRandomPoints$propSCAA_25Mb <- scaaPerRandomPoints$n.scaa_within25Mb / scaaPerRandomPoints$n.peaks_within25Mb
scaaPerRandomPoints$propSCAA_500kb <- scaaPerRandomPoints$n.scaa_within500kb / scaaPerRandomPoints$n.peaks_within500kb



# DISTANCE TO NEAREST BP =====

peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")] # keep only carcinomas ('pure')
scaa.location <- data.frame(chr=sub("\\:.*", "", rownames(scaa)), start=sub(".*[:]([^.]+)[-].*", "\\1", rownames(scaa)), stop=sub(".*-", "", rownames(scaa)))

output <- list()
l <- 1
for (i in 1:ncol(scaa)) {
  print(paste(i,"/",ncol(scaa)))
  patient <- substr(colnames(scaa),1,4)[i]
  
  for (j in 1:nrow(scaa)) {
    chr <- sub("\\:.*", "", rownames(scaa)[j])
    wd <- breakpoints_dWGS[which(breakpoints_dWGS$patient==patient & breakpoints_dWGS$break_chr==chr),]
    peak_base <- as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", rownames(scaa)[j])) +250
    index <- which.min(apply(wd[6], 2, function(x) {abs(x-scaa_base)}))
    clonality <- wd[index,10]
    nearest_BP <- wd[index,6]
    distance_BP <- abs(nearest_BP-peak_base)
    output[[l]] <- setNames(c(patient, rownames(scaa)[j], chr, peak_base, scaa[j,i], nearest_BP, distance_BP, clonality),
                            c("patient", "peak", "chr", "peak_base", "scaa","nearest_BP", "distance_BP", "BP_clonality"))
    l <- l + 1
  }
}
  
breakpointPerScaa <- data.frame(do.call(rbind, output))
breakpointPerScaa$distance_BP <- as.numeric(breakpointPerScaa$distance_BP)
breakpointPerScaa$BP_clonality <- as.numeric(breakpointPerScaa$BP_clonality)
breakpointPerScaa$promoter <- peakAnn_all[match(breakpointPerScaa$peak, peakAnn_all$peak), 4]


# SAVE =====
saveRDS(scaaPerBreakpoint, "~/Documents/SCAA/Data/EPICC/SCAA/scaaPerBreakpoint.rds")
saveRDS(scaaPerRandomPoints, "~/Documents/SCAA/Data/EPICC/SCAA/scaaPerRandomPoints.rds")
saveRDS(breakpointPerScaa, "~/Documents/SCAA/Data/EPICC/SCAA/breakpointPerScaa.rds")






