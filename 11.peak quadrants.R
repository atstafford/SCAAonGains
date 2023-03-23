# PREREQUISITIS: load section 1, 7 data 

# Idea from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6984959/

# SET QUADRANTS ACORSS ALL PEAKS =====

# Splitting intra-tumour RNA heterogeneity and inter-tumour RNA heterogeneity by their respective average (mean)
# value generates RNA heterogeneity quadrants
mean.pITH <- mean(regionwise.pITH$pITH)
mean.pIPH <- mean(cohortwise.pIPH$pIPH)

yhist.pITH <- ggplot(regionwise.pITH, aes(x=pITH)) + 
  geom_histogram(aes(y=..density..), color="black", fill="white", bins = 100) + 
  geom_vline(aes(xintercept=mean(pITH)), color="black", linetype="dashed", size=1) +
  geom_density(alpha=.2, fill="#FF6666") +
  xlim(0+0.25,(2*mean.pITH)-0.25) +
  coord_flip() +
  scale_y_reverse() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(fill=NA),
        plot.margin=margin(t=0,r=0,b=0,l=0,"cm"))

xhist.pIPH <- ggplot(cohortwise.pIPH, aes(x=pIPH)) + 
  geom_histogram(aes(y=..density..), color="black", fill="white", bins = 100) + 
  geom_vline(aes(xintercept=mean(pIPH)), color="black", linetype="dashed", size=1) +
  geom_density(alpha=.2, fill="#FF6666") +
  xlim(0,(2*mean.pIPH)) +
  scale_y_reverse() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(fill=NA),
        plot.margin=margin(t=0,r=0,b=0,l=0,"cm"))

low.pITH <- rownames(regionwise.pITH)[which(regionwise.pITH$pITH < mean.pITH)]
high.pITH <- rownames(regionwise.pITH)[which(regionwise.pITH$pITH > mean.pITH)]
low.pIPH <- rownames(cohortwise.pIPH)[which(cohortwise.pIPH$pIPH < mean.pIPH)]
high.pIPH <- rownames(cohortwise.pIPH)[which(cohortwise.pIPH$pIPH > mean.pIPH)]
x <- rownames(regionwise.pITH)[which(regionwise.pITH$pITH == mean.pITH)] # none equal mean
x <- rownames(cohortwise.pIPH)[which(cohortwise.pIPH$pIPH == mean.pIPH)]

q2.peaks <- intersect(low.pIPH, low.pITH)
q1.peaks <- intersect(low.pIPH, high.pITH) 
q4.peaks <- intersect(high.pIPH, low.pITH) 
q3.peaks <- intersect(high.pIPH, high.pITH) 

nrow(regionwise.pITH)-length(q1.peaks)-length(q2.peaks)-length(q3.peaks)-length(q4.peaks)

# SET QUADRANTS ACORSS PROMOTERS ONLY =====

PE.peaks <- peakAnn_all$peak[which(peakAnn_all$Promoter == "TRUE" | !is.na(peakAnn_all$elementType))]

regionwise.pITH_PE <- data.frame(pITH=regionwise.pITH[which(rownames(regionwise.pITH) %in% PE.peaks),], row.names = rownames(regionwise.pITH)[which(rownames(regionwise.pITH) %in% PE.peaks)])
cohortwise.pIPH_PE <- data.frame(pIPH=cohortwise.pIPH[which(rownames(cohortwise.pIPH) %in% PE.peaks),], row.names = rownames(cohortwise.pIPH)[which(rownames(cohortwise.pIPH) %in% PE.peaks)])

mean.pITH_PE <- mean(regionwise.pITH_PE$pITH)
mean.pIPH_PE <- mean(cohortwise.pIPH_PE$pIPH)

yhist.pITH_PE <- ggplot(regionwise.pITH_PE, aes(x=pITH)) + 
  geom_histogram(aes(y=..density..), color="black", fill="white", bins = 100) + 
  geom_vline(aes(xintercept=mean(pITH)), color="black", linetype="dashed", size=1) +
  geom_density(alpha=.2, fill="#FF6666") +
  xlim(0+0.25,(2*mean.pITH_PE)-0.25) +
  coord_flip() +
  scale_y_reverse() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(fill=NA),
        plot.margin=margin(t=0,r=0,b=0,l=0,"cm"))

xhist.pIPH_PE <- ggplot(cohortwise.pIPH_PE, aes(x=pIPH)) + 
  geom_histogram(aes(y=..density..), color="black", fill="white", bins = 100) + 
  geom_vline(aes(xintercept=mean(pIPH)), color="black", linetype="dashed", size=1) +
  geom_density(alpha=.2, fill="#FF6666") +
  xlim(0,(2*mean.pIPH_PE)) +
  scale_y_reverse() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(fill=NA),
        plot.margin=margin(t=0,r=0,b=0,l=0,"cm"))

low.pITH <- rownames(regionwise.pITH_PE)[which(regionwise.pITH_PE$pITH < mean.pITH_PE)]
high.pITH <- rownames(regionwise.pITH_PE)[which(regionwise.pITH_PE$pITH > mean.pITH_PE)]
low.pIPH <- rownames(cohortwise.pIPH_PE)[which(cohortwise.pIPH_PE$pIPH < mean.pIPH_PE)]
high.pIPH <- rownames(cohortwise.pIPH_PE)[which(cohortwise.pIPH_PE$pIPH > mean.pIPH_PE)]
x <- rownames(regionwise.pITH_PE)[which(regionwise.pITH_PE$pITH == mean.pITH_PE)] # none equal mean
x <- rownames(cohortwise.pIPH_PE)[which(cohortwise.pIPH_PE$pIPH == mean.pIPH_PE)]

q2.peaks_PE <- intersect(low.pIPH, low.pITH)
q1.peaks_PE <- intersect(low.pIPH, high.pITH) 
q4.peaks_PE <- intersect(high.pIPH, low.pITH) 
q3.peaks_PE <- intersect(high.pIPH, high.pITH) 

nrow(regionwise.pITH_PE)-length(q1.peaks_PE)-length(q2.peaks_PE)-length(q3.peaks_PE)-length(q4.peaks_PE)

# PEAK ANNOTATION OF QUADRANTS --------------------------------------------------------------------------------------------------------------------

# annotate peaks per quarter
q1.peakAnn <- peakAnn_all[which(peakAnn_all$peak %in% q1.peaks),]
q2.peakAnn <- peakAnn_all[which(peakAnn_all$peak %in% q2.peaks),]
q3.peakAnn <- peakAnn_all[which(peakAnn_all$peak %in% q3.peaks),]
q4.peakAnn <- peakAnn_all[which(peakAnn_all$peak %in% q4.peaks),]

# annotate peaks per quarter PE only
q1.peakAnn_PE <- peakAnn_all[which(peakAnn_all$peak %in% q1.peaks_PE),]
q2.peakAnn_PE <- peakAnn_all[which(peakAnn_all$peak %in% q2.peaks_PE),]
q3.peakAnn_PE <- peakAnn_all[which(peakAnn_all$peak %in% q3.peaks_PE),]
q4.peakAnn_PE <- peakAnn_all[which(peakAnn_all$peak %in% q4.peaks_PE),]

# PLOTS WITH QUADRANTS BASED ON ALL PEAKS =====

# plot all peaks
outsidebox_size <- c(length(q3.peaks), length(q1.peaks), length(q2.peaks), length(q4.peaks))
insidebox_size <- NULL
outsidebox_label <- c(length(q3.peaks), length(q1.peaks), length(q2.peaks), length(q4.peaks))
insidebox_label <- NULL
quad.plot_allpeaks <- four_quadrant(outsidebox_size = outsidebox_size, insidebox_size = insidebox_size, 
                                    outsidebox_label = outsidebox_label, insidebox_label = insidebox_label, p.value = NA)

p1 <- plot_grid(yhist.pITH, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(quad.plot_allpeaks, xhist.pIPH, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# occurrence of promo/enh in each quarter: using both promoter and GH column
x <- quadEnrichmentChi(enrichedPeakList = PE.peaks, q1.peaks = q1.peaks, q2.peaks = q2.peaks, q3.peaks =q3.peaks, q4.peaks = q4.peaks)
p1 <- plot_grid(yhist.pITH, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x$quad.plot, xhist.pIPH, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# all SCAAs: enriched in Q1/4
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")] #keep 'pures'
scaa <- scaa[which(rownames(scaa) %in% peakAnn_all$peak),] #remove scaas missing cna data
scaa <- scaa[which(rowSums(scaa)>=1),] #keep peaks that are a SCAA in at least one patient
scaaPeaks <- rownames(scaa)

x <- quadEnrichmentChi(enrichedPeakList = scaaPeaks, q1.peaks = q1.peaks, q2.peaks = q2.peaks, q3.peaks =q3.peaks, q4.peaks = q4.peaks)
p1 <- plot_grid(yhist.pITH, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x$quad.plot, xhist.pIPH, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# recurring (>20%) SCAAs: enriched in Q3/4
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")] #keep 'pures'
scaa <- scaa[which(rownames(scaa) %in% peakAnn_all$peak),] #remove scaas missing cna data
scaa <- scaa[which(rowSums(scaa)>=(ncol(scaa)*0.2)),] #keep peaks that are a SCAA in at least 20% patient
scaaPeaks <- rownames(scaa)

x <- quadEnrichmentChi(enrichedPeakList = scaaPeaks, q1.peaks = q1.peaks, q2.peaks = q2.peaks, q3.peaks =q3.peaks, q4.peaks = q4.peaks)
p1 <- plot_grid(yhist.pITH, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x$quad.plot, xhist.pIPH, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# PLOTS WITH QUADRANTS BASED ON PE PEAKS =====

# plot all PE peaks
outsidebox_size <- c(length(q3.peaks_PE), length(q1.peaks_PE), length(q2.peaks_PE), length(q4.peaks_PE))
insidebox_size <- NULL
outsidebox_label <- c(length(q3.peaks_PE), length(q1.peaks_PE), length(q2.peaks_PE), length(q4.peaks_PE))
insidebox_label <- NULL
quad.plot <- four_quadrant(outsidebox_size = outsidebox_size, insidebox_size = insidebox_size, 
                           outsidebox_label = outsidebox_label, insidebox_label = insidebox_label, p.value = NA)

p1 <- plot_grid(yhist.pITH_PE, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(quad.plot, xhist.pIPH_PE, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# Hallmark evenly dist
HM.PEpeaks <- peakAnn$peak[which(peakAnn$SYMBOL %in% hallmarkGenes$gene_symbol)]
x <- quadEnrichmentChi(enrichedPeakList = HM.PEpeaks, q1.peaks = q1.peaks_PE, q2.peaks = q2.peaks_PE, q3.peaks =q3.peaks_PE, q4.peaks = q4.peaks_PE)
p1 <- plot_grid(yhist.pITH_PE, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x, xhist.pIPH_PE, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# cosmic cancer genes are overrep in Q2
CS.PEpeaks <- peakAnn$peak[which(peakAnn$SYMBOL %in% cosmic$`Gene Symbol`)]
x <- quadEnrichmentChi(enrichedPeakList = CS.PEpeaks, q1.peaks = q1.peaks_PE, q2.peaks = q2.peaks_PE, q3.peaks =q3.peaks_PE, q4.peaks = q4.peaks_PE)
p1 <- plot_grid(yhist.pITH_PE, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x, xhist.pIPH_PE, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# chromatin cancer genes are overrep in Q3
CG.PEpeaks <- peakAnn$peak[which(peakAnn$SYMBOL %in% chromatinGenes$Symbol)]
x <- quadEnrichmentChi(enrichedPeakList = CG.PEpeaks, q1.peaks = q1.peaks_PE, q2.peaks = q2.peaks_PE, q3.peaks =q3.peaks_PE, q4.peaks = q4.peaks_PE)
p1 <- plot_grid(yhist.pITH_PE, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x, xhist.pIPH_PE, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# DNAmodGenes genes evenly dist
DG.PEpeaks <- peakAnn$peak[which(peakAnn$SYMBOL %in% DNAmodGenes$Symbol)]
x <- quadEnrichmentChi(enrichedPeakList = DG.PEpeaks, q1.peaks = q1.peaks_PE, q2.peaks = q2.peaks_PE, q3.peaks =q3.peaks_PE, q4.peaks = q4.peaks_PE)
p1 <- plot_grid(yhist.pITH_PE, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x, xhist.pIPH_PE, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# histone genes are overrep in Q3
HG.PEpeaks <- peakAnn$peak[which(peakAnn$SYMBOL %in% histoneGenes$Symbol)]
x <- quadEnrichmentChi(enrichedPeakList = HG.PEpeaks, q1.peaks = q1.peaks_PE, q2.peaks = q2.peaks_PE, q3.peaks =q3.peaks_PE, q4.peaks = q4.peaks_PE)
p1 <- plot_grid(yhist.pITH_PE, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x, xhist.pIPH_PE, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# stable genes are evenly dist
stable.genes <- data.frame(gene=c("G6PD","ACTB","GAPDH","HPRT1","B2M","GUSB","IPO8","UBC","PPIA","TBP","HMBS","YWHAZ","PGK1","TFRC"),
                           ID=c("2539","60","2597","3251","567","2990","10526","7316","5478","6908","3145","7534","5230","7037"))
SG.PEpeaks <- peakAnn$peak[which(peakAnn$geneID %in% stable.genes$ID)]
x <- quadEnrichmentChi(enrichedPeakList = SG.PEpeaks, q1.peaks = q1.peaks_PE, q2.peaks = q2.peaks_PE, q3.peaks =q3.peaks_PE, q4.peaks = q4.peaks_PE)
p1 <- plot_grid(yhist.pITH_PE, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x, xhist.pIPH_PE, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# all SCAAs: enriched in Q1/2/3
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")] #keep 'pures'
scaa <- scaa[which(rownames(scaa) %in% peakAnn$peak),] #remove scaas missing cna data
scaa <- scaa[which(rowSums(scaa)>=1),] #keep peaks that are a SCAA in at least one patient
scaaPeaks <- rownames(scaa)

x <- quadEnrichmentChi(enrichedPeakList = scaaPeaks, q1.peaks = q1.peaks_PE, q2.peaks = q2.peaks_PE, q3.peaks =q3.peaks_PE, q4.peaks = q4.peaks_PE)
p1 <- plot_grid(yhist.pITH_PE, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x$quad.plot, xhist.pIPH_PE, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3

# recurring (>20%) SCAAs: enriched in Q1/2/4
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")] #keep 'pures'
scaa <- scaa[which(rownames(scaa) %in% peakAnn$peak),] #remove scaas missing cna data
scaa <- scaa[which(rowSums(scaa)>=(ncol(scaa)*0.2)),] #keep peaks that are a SCAA in at least 20% patient
scaaPeaks <- rownames(scaa)

x <- quadEnrichmentChi(enrichedPeakList = scaaPeaks, q1.peaks = q1.peaks_PE, q2.peaks = q2.peaks_PE, q3.peaks =q3.peaks_PE, q4.peaks = q4.peaks_PE)
p1 <- plot_grid(yhist.pITH_PE, NULL, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p2 <- plot_grid(x$quad.plot, xhist.pIPH_PE, ncol=1, align = "v", rel_heights = c(1,0.5), rel_widths = c(1))
p3 <- plot_grid(p1, p2, ncol=2, align = "v", axis = "lr",  rel_heights = c(1), rel_widths = c(0.5,1))
p3


# SAVE -----------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(q1.peaks, "~/Documents/SCAA/Data/q1.peaks.rds")
saveRDS(q2.peaks, "~/Documents/SCAA/Data/q2.peaks.rds")
saveRDS(q3.peaks, "~/Documents/SCAA/Data/q3.peaks.rds")
saveRDS(q4.peaks, "~/Documents/SCAA/Data/q4.peaks.rds")

saveRDS(q1.peaks_PE, "~/Documents/SCAA/Data/q1.peaks_PE.rds")
saveRDS(q2.peaks_PE, "~/Documents/SCAA/Data/q2.peaks_PE.rds")
saveRDS(q3.peaks_PE, "~/Documents/SCAA/Data/q3.peaks_PE.rds")
saveRDS(q4.peaks_PE, "~/Documents/SCAA/Data/q4.peaks_PE.rds")

saveRDS(q1.peakAnn, "~/Documents/SCAA/Data/q1.peakAnn.rds")
saveRDS(q2.peakAnn, "~/Documents/SCAA/Data/q2.peakAnn.rds")
saveRDS(q3.peakAnn, "~/Documents/SCAA/Data/q3.peakAnn.rds")
saveRDS(q4.peakAnn, "~/Documents/SCAA/Data/q4.peakAnn.rds")

saveRDS(q1.peakAnn_PE, "~/Documents/SCAA/Data/q1.peakAnn_PE.rds")
saveRDS(q2.peakAnn_PE, "~/Documents/SCAA/Data/q2.peakAnn_PE.rds")
saveRDS(q3.peakAnn_PE, "~/Documents/SCAA/Data/q3.peakAnn_PE.rds")
saveRDS(q4.peakAnn_PE, "~/Documents/SCAA/Data/q4.peakAnn_PE.rds")

saveRDS(PE.peaks, "~/Documents/SCAA/Data/PE.peaks.rds")
saveRDS(regionwise.pITH_PE, "~/Documents/SCAA/Data/regionwise.pITH_PE.rds")

saveRDS(xhist.pIPH, "~/Documents/SCAA/Data/xhist.pIPH.rds")
saveRDS(yhist.pITH, "~/Documents/SCAA/Data/yhist.pITH.rds")

