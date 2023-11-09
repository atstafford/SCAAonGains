# Output from section 1 (dWGS data setup) =============================================
dwgs <- readRDS("~/Documents/SCAA/Data/dwgs.rds")
dwgs.ploidyRecentre <- readRDS("~/Documents/SCAA/Data/dwgs.ploidyRecentre.rds")
ploidy.persample <- readRDS("~/Documents/SCAA/Data/ploidy.persample.rds")
ploidy.perpatient <- readRDS("~/Documents/SCAA/Data/ploidy.perpatient.rds")
dwgs.perpatient <- readRDS("~/Documents/SCAA/Data/dwgs.perpatient.rds")
dwgs.ploidyRecentre.perpatient <- readRDS("~/Documents/SCAA/Data/dwgs.ploidyRecentre.perpatient.rds")

# Output from section 2 (CN normalisation) =============================================
atac_raw <- readRDS("~/Documents/SCAA/Data/atac_raw.rds")
atac_pur <- readRDS("~/Documents/SCAA/Data/atac_pur.rds")
wgs_pur <- readRDS("~/Documents/SCAA/Data/wgs_pur.rds")
atac_cn.deepWGSnorm <- readRDS("~/Documents/SCAA/Data/atac_cn.deepWGSnorm.rds")

# Output from section 3 (deseq2 of dWGS normalised across cohort data, ignore that it says per patient) =========
metaData.dWGS.cohort <- readRDS("~/Documents/SCAA/Data/metaData.dWGS.cohort.rds")
dds.dWGS.cohort <- readRDS("~/Documents/SCAA/Data/dds.dWGS.cohort.rds")
deseq.dWGS.cohort <- readRDS("~/Documents/SCAA/Data/deseq.dWGS.cohort.rds")
deseqVst.dWGS.cohort <- readRDS("~/Documents/SCAA/Data/deseqVst.dWGS.cohort.rds")

metaData.perpatient <- readRDS("~/Documents/SCAA/Data/metaData.perpatient.rds") #old
dds_rlog.perpatient <- readRDS("~/Documents/SCAA/Data/dds_rlog.perpatient.rds") #old

# Output from section 4 (segment data) ==============
trueDiploids <- readRDS("~/Documents/SCAA/Data/trueDiploids.rds")
segments <- readRDS("~/Documents/SCAA/Data/segments.rds")
segments_ploidyrecenterd <- readRDS("~/Documents/SCAA/Data/segments_ploidyrecenterd.rds")
scaa <- readRDS("~/Documents/SCAA/Data/scaa.rds")

# Output from section 5 (peak-level data) ==============
scaa.long <- readRDS("~/Documents/SCAA/Data/scaa.long.rds")

# Output from section 6 (patient-level data) ==============
scaaPP <- readRDS("~/Documents/SCAA/Data/scaaPP.rds")

# Output from section 7 (breakpoint data) ==============
breakpoints_dWGS <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/breakpoints_dWGS.rds")
randomPoints <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/randomPoints.rds")

# Output from section 8 (SCAAs at breakpoint data) ==============
scaaPerBreakpoint <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/scaaPerBreakpoint.rds")
scaaPerRandomPoints <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/scaaPerRandomPoints.rds")
breakpointPerScaa <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/breakpointPerScaa.rds")









# Output from section 4 (CNA data per peak) ==============
ploidy.persample <- readRDS("~/Documents/SCAA/Data/ploidy.persample.rds")
ploidy.perpatient <- readRDS("~/Documents/SCAA/Data/ploidy.perpatient.rds")

dwgs.perpatient <- readRDS("~/Documents/SCAA/Data/dwgs.perpatient.rds")
dwgs.ploidyRecentre.perpatient <- readRDS("~/Documents/SCAA/Data/dwgs.ploidyRecentre.perpatient.rds")

#cnaPerScaa <- readRDS("~/Documents/SCAA/Data/cnaPerScaa.rds")
#focal <- readRDS("~/Documents/SCAA/Data/focal.rds")
#arm <- readRDS("~/Documents/SCAA/Data/arm.rds")
scaaCNA <- readRDS("~/Documents/SCAA/Data/scaaCNA.rds")
gainInfo <- readRDS("~/Documents/SCAA/Data/gainInfo.rds")

diploidPerScaa <- readRDS("~/Documents/SCAA/Data/diploidPerScaa.rds")
cnaPerScaa <- readRDS("~/Documents/SCAA/Data/cnaPerScaa.rds")
cnaPerScaa.ploidyRecenter <- readRDS("~/Documents/SCAA/Data/cnaPerScaa.ploidyRecenter.rds")

focal <- readRDS("~/Documents/SCAA/Data/focal.rds")
arm <- readRDS("~/Documents/SCAA/Data/arm.rds")
pcArm.ploidyRecenter <-readRDS("~/Documents/SCAA/Data/pcArm.ploidyRecenter.rds")
arm.ploidyRecenter <- readRDS("~/Documents/SCAA/Data/arm.ploidyRecenter.rds")

scaaCNA_dip <- readRDS("~/Documents/SCAA/Data/scaaCNA_dip.rds")
scaaCNA <- readRDS("~/Documents/SCAA/Data/scaaCNA.rds")
scaaCNA.ploidyRecenter <- readRDS("~/Documents/SCAA/Data/scaaCNA.ploidyRecenter.rds")

diploidInfo <- readRDS("~/Documents/SCAA/Data/diploidInfo.rds")
segInfo <- readRDS("~/Documents/SCAA/Data/segInfo.rds")
segInfo.ploidyRecentre <- readRDS("~/Documents/SCAA/Data/segInfo.ploidyRecentre.rds")

trueDiploid.perpatient <- readRDS("~/Documents/SCAA/Data/trueDiploid.perpatient.rds")
trueDiploids <- readRDS("~/Documents/SCAA/Data/trueDiploids.rds")

dwgs.ploidyRecentre.perpatient <- readRDS("~/Documents/SCAA/Data/dwgs.ploidyRecentre.perpatient.rds")
#cnaPerScaa.ploidyRecenter <- readRDS("~/Documents/SCAA/Data/cnaPerScaa.ploidyRecenter.rds")
#pcArm.ploidyRecenter <- readRDS("~/Documents/SCAA/Data/pcArm.ploidyRecenter.rds")
#arm.ploidyRecenter <- readRDS("~/Documents/SCAA/Data/arm.ploidyRecenter.rds")
scaaCNA.ploidyRecenter <- readRDS("~/Documents/SCAA/Data/scaaCNA.ploidyRecenter.rds")

# Output from section 5 (peak coverage) ==============
peakperCNA <- readRDS("~/Documents/SCAA/Data/peakperCNA.rds")

# Output from section 6 (SCAA cause) ==============
scaaPP <- readRDS("~/Documents/SCAA/Data/scaaPP.rds")

# Output from section 7 (breakpoints) ==============
saveRDS(breakpoints_dWGS, "~/Documents/SCAA/Data/breakpoints_dWGS.rds")
saveRDS(dat, "~/Documents/SCAA/Data/dat.rds")