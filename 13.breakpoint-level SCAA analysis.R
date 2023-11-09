# DATA =====
scaaPerBreakpoint <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/scaaPerBreakpoint.rds")
scaaPerRandomPoints <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/scaaPerRandomPoints.rds")

# with over 45 peaks
x <- scaaPerBreakpoint[which(scaaPerBreakpoint$gap==FALSE & scaaPerBreakpoint$n.peaks_within1Mb>0 ),c(1,11:32)]
y <- scaaPerRandomPoints[scaaPerRandomPoints$n.peaks_within1Mb>0, c(3,5:26)]
#y <- y[c(2,1,3,4,5,6)]
#colnames(y)[1] <- c("patient")
x$group <- "Breakpoints"
y$group <- "Random loci"
#x <- x[!duplicated(x),]
#set.seed(1234)
#y <- y[sample(nrow(y), nrow(x)),]

data <- rbind(x,y)
data$pga <- scaaPP[match(data$patient, scaaPP$patient),8]

# COMPARISON BASED ON SCAA WITHIN 1MB ----------------------------------------------------------------------------------------------

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data, aes(x=group, y=log(propSCAA_20Mb), fill=group)) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitter(0.15)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  scale_y_continuous(name="SCAA burden within 1Mb of point", 
                     breaks = c(-4,-2,0),
                     labels = c(round(exp(-4),2),round(exp(-2),2),round(exp(-0),2))) +
  scale_x_discrete(labels = c('Breakpoints','Random loci')) +
  scale_fill_manual(values = c("#2c728eff","#20a486ff")) +
  theme_custom() +
  stat_compare_means(method = "wilcox.test", label.x.npc = 0.3,label.y = 1,  size=8) +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank())
dev.off()



# paried violin
#x <- scaaPerBreakpoint[which(scaaPerBreakpoint$gap==FALSE & scaaPerBreakpoint$n.peaks_within100kb>5 ),c(30:37)]
x <- scaaPerBreakpoint[which(scaaPerBreakpoint$gap==FALSE & scaaPerBreakpoint$n.peaks_within1Mb>49 ),c(30:37)]
y <- scaaPerRandomPoints[scaaPerRandomPoints$n.peaks_within1Mb>49, c(24:31)]
#y <- scaaPerRandomPoints[scaaPerRandomPoints$n.peaks_within100kb>5, c(24:31)]

x$group <- "Breakpoint"
y$group <- "Random"

data <- rbind(x,y)

wilcox.test(log(data$propSCAA_100kb[which(data$group=="Breakpoint")]), log(data$propSCAA_100kb[which(data$group=="Random")])) 
mean(data$propSCAA_1.5Mb[which(data$group=="Breakpoint")])
mean(data$propSCAA_1.5Mb[which(data$group=="Random")])

colnames(data) <- c("100kb","250kb","500kb","750kb","1Mb","1.5Mb","5Mb","20Mb","group")
data <- gather(data, key="key", value="value", -group)
data$key <- factor(data$key, levels = c("100kb","250kb","500kb","750kb","1Mb","1.5Mb","5Mb","20Mb"))

library(ggpubr)

jpeg('tempfig.jpeg', width = (40), height = (20), units = "cm", res = 300)
ggplot2::ggplot(data[-which(data$key %in% c("5Mb","20Mb")),], aes(y=log(value), x=(key), fill=group)) +
  geom_split_violin(width=0.95, alpha = .4, trim = FALSE, scale = "width") +
  geom_boxplot(width = .4, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  geom_jitter(size=2, alpha = .4, shape=21, position=position_jitterdodge(jitter.width = 0.1,
                                                                          jitter.height = 0,
                                                                          dodge.width = 0.6,)) +
  scale_y_continuous(name="SCAA burden", 
                     limits = c(-9,2),
                     breaks = c(-6,-3,0),
                     labels = c(round(exp(-6),5),round(exp(-3),4),round(exp(0),2))) +
  xlab("Distance from locus (Mb)") +
  scale_fill_manual(values=c("#2c728eff","#20a486ff")) +
  theme_custom() +
  #coord_flip() +
  stat_compare_means(aes(group = group), label="p.format", label.y = 2,  size=9, method = "wilcox.test") +
  theme(legend.position = "right", legend.text = element_text(size=28), legend.title = element_blank())
dev.off()


t.test(log(data$propSCAA_100kb[which(data$group=="Breakpoints" & data$propSCAA_100kb>0)]), log(data$propSCAA_100kb[which(data$group=="Random loci" & data$propSCAA_100kb>0)])) 
t.test(log(data$propSCAA_500kb[which(data$group=="Breakpoints" & data$propSCAA_500kb>0)]), log(data$propSCAA_500kb[which(data$group=="Random loci" & data$propSCAA_500kb>0)])) 
t.test(log(data$propSCAA_750kb[which(data$group=="Breakpoints" & data$propSCAA_750kb>0)]), log(data$propSCAA_750kb[which(data$group=="Random loci" & data$propSCAA_750kb>0)])) 
t.test(log(data$propSCAA_1Mb[which(data$group=="Breakpoints" & data$propSCAA_1Mb>0)]), log(data$propSCAA_1Mb[which(data$group=="Random loci" & data$propSCAA_1Mb>0)])) 
t.test(log(data$propSCAA_5Mb[which(data$group=="Breakpoints" & data$propSCAA_5Mb>0)]), log(data$propSCAA_5Mb[which(data$group=="Random loci" & data$propSCAA_5Mb>0)])) 


# per patient
data <- rbind(data.frame(value=(scaaPerBreakpoint$propSCAA_20Mb),
                         group="Breakpoints",
                         patient=scaaPerBreakpoint$patient),
              data.frame(value=(scaaPerRandomPoints$propSCAA_20Mb),
                         group="Not Breakpoints",
                         patient=scaaPerRandomPoints$patient))

#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*25))
jpeg('tempfig2.jpeg', width = (10), height = (60), units = "cm", res = 300)
ggplot(data, aes(x=group, y=log(value))) +
  geom_boxplot(width=1) +
  #geom_jitter(shape=16, position=position_jitter(0.4)) +
  #geom_violin(aes(fill=group)) +
  facet_grid(rows=vars(patient)) +
  stat_compare_means(method = "wilcox.test", label.x.npc = 0.3,label.y = -2,  size=4) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(fill=NA))
dev.off()

# what proportion of SCAA burden per patient is near a BP?

# COMPARISON BASED ON DISTANCE TO NEAREST SCAA ----------------------------------------------------------------------------------------------

# remove 3 patients
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(data, aes(x=group, y=log(distance_scaa), fill=group)) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitter(0.15)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  #geom_hline(yintercept = log(1000000)) +
  scale_y_continuous(name="Distance to nearest SCAA", 
                     breaks = c(5,10,15,20),
                     #labels = c(signif(exp(5),2),signif(exp(10),2),signif(exp(15),2),signif(exp(20),6))) +
                     labels = c("150bp","22kb","3.3Mb","485Mb")) +
  scale_x_discrete(labels = c('Breakpoints','Random loci')) +
  scale_fill_manual(values = c("#2c728eff","#20a486ff")) +
  theme_custom() +
  stat_compare_means(method = "t.test", label.x.npc = 0.3,label.y = log(500000000),  size=8) +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank())
dev.off()

ks.test(data$value[which(data$group=="Breakpoints")], data$value[which(data$group=="Not Breakpoints")]) 

# per patient
data <- rbind(data.frame(value=(scaaPerBreakpoint$distance_scaa[which(scaaPerBreakpoint$patient=="C531")]),group="Breakpoints"),
              data.frame(value=(scaaPerRandomPoints$distance_scaa[which(scaaPerBreakpoint$patient=="C531")]),group="Not Breakpoints"))

ggplot(data, aes(x=group, y=(value))) +
  geom_boxplot(width=0.4) +
  geom_jitter(shape=16, position=position_jitter(0.4)) +
  geom_violin(aes(fill=group)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(fill=NA))

ks.test(data$value[which(data$group=="Breakpoints")], data$value[which(data$group=="Not Breakpoints")]) 

# COMPARISON BASED ON DISTANCE TO NEAREST BP ----------------------------------------------------------------------------------------------
x <- breakpointPerScaa[!is.na(breakpointPerScaa$distance_BP),]

# remove 3 patients
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
ggplot(x[which(x$BP_clonality>0.5 & x$promoter==FALSE),], aes(x=as.character(scaa), y=log(distance_BP), fill=scaa)) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitter(0.15)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  #geom_hline(yintercept = log(1000000)) +
  #scale_y_continuous(name="Distance to nearest SCAA", 
  #                   breaks = c(5,10,15,20),
  #                   #labels = c(signif(exp(5),2),signif(exp(10),2),signif(exp(15),2),signif(exp(20),6))) +
  #                   labels = c("150bp","22kb","3.3Mb","485Mb")) +
  #scale_x_discrete(labels = c('Breakpoints','Random loci')) +
  #scale_fill_manual(values = c("#2c728eff","#20a486ff")) +
  theme_custom() +
  #stat_compare_means(method = "t.test", label.x.npc = 0.3,label.y = log(500000000),  size=8) +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank())
dev.off()

t.test(log(x$distance_BP[which(x$scaa==1)]), 
        log(x$distance_BP[which(x$scaa==0)]) )

x <- breakpointPerScaa[which(breakpointPerScaa$scaa==0),]

# per patient
data <- rbind(data.frame(value=(scaaPerBreakpoint$distance_scaa[which(scaaPerBreakpoint$patient=="C531")]),group="Breakpoints"),
              data.frame(value=(scaaPerRandomPoints$distance_scaa[which(scaaPerBreakpoint$patient=="C531")]),group="Not Breakpoints"))

ggplot(data, aes(x=group, y=(value))) +
  geom_boxplot(width=0.4) +
  geom_jitter(shape=16, position=position_jitter(0.4)) +
  geom_violin(aes(fill=group)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(fill=NA))

ks.test(data$value[which(data$group=="Breakpoints")], data$value[which(data$group=="Not Breakpoints")]) 

# SCAA PROP AS DISTANCE FROM BP INCREASES =====

# on a per patient basis
x <- scaaPerBreakpoint[,colnames(scaaPerBreakpoint) %in% c("patient","group","propSCAA_100kb","propSCAA_500kb","propSCAA_750kb","propSCAA_1Mb","propSCAA_20Mb","propSCAA_pp")]
x$group <- "breakpoint"
y <- scaaPerRandomPoints[,colnames(scaaPerRandomPoints) %in% c("patient","group","propSCAA_100kb","propSCAA_500kb","propSCAA_750kb","propSCAA_1Mb","propSCAA_20Mb","propSCAA_pp")]
y$group <- "random"

x <- rbind(x,y)
#x <- x[which(x$patient=="C524"),]
x$propSCAA_100kb <- x$propSCAA_100kb/x$propSCAA_pp
x$propSCAA_500kb <- x$propSCAA_500kb/x$propSCAA_pp
x$propSCAA_750kb <- x$propSCAA_750kb/x$propSCAA_pp
x$propSCAA_1Mb <- x$propSCAA_1Mb/x$propSCAA_pp
x$propSCAA_20Mb <- x$propSCAA_20Mb/x$propSCAA_pp
#x$propSCAA_25Mb <- x$propSCAA_25Mb/x$propSCAA_pp
#x$propSCAA_50Mb <- x$propSCAA_50Mb/x$propSCAA_pp
x <- x[-c(1)]

x <- gather(x, key="key", value="value", -group)
x$key <- factor(x$key, levels = c("propSCAA_100kb","propSCAA_500kb","propSCAA_750kb","propSCAA_1Mb","propSCAA_20Mb"))
ggplot(x, aes(x=key, y=log(value,base=2), fill=as.factor(group))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
 # geom_jitter(size=3, alpha = .4, shape=21, position=position_jitter(0.15)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  geom_hline(yintercept = 0) +
  #scale_y_continuous(name="SCAA burden within 1Mb of point", 
  #                   breaks = c(-4,-2,0),
  #                   labels = c(round(exp(-4),2),round(exp(-2),2),round(exp(-0),2))) +
  #scale_x_discrete(labels = c('Breakpoints','Random loci')) +
  #scale_fill_manual(values = c("#2c728eff","#20a486ff")) +
  theme_custom() +
  #stat_compare_means(method = "t.test", label.x.npc = 0.3,label.y = 1,  size=8) +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=5))

x <- scaaPerRandomPoints[which(scaaPerRandomPoints$patient_sub=="C552"),6:14]
x <- gather(x, key="key", value="value")
x$key <- factor(x$key, levels = c("n.scaa_within100kb","n.scaa_within500kb","n.scaa_within750kb","n.scaa_within1Mb","n.scaa_within5Mb","n.scaa_within10Mb","n.scaa_within20Mb","n.scaa_within25Mb","n.scaa_within50Mb"))
ggplot(x, aes(x=key, y=(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitter(0.15)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  geom_hline(yintercept = 0) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=5))

x <- scaaPerRandomPoints[which(scaaPerRandomPoints$patient_sub=="C552"),17:25]
x <- gather(x, key="key", value="value")
x$key <- factor(x$key, levels = c("n.peaks_within100kb","n.peaks_within500kb","n.peaks_within750kb","n.peaks_within1Mb","n.peaks_within5Mb","n.peaks_within10Mb","n.peaks_within20Mb","n.peaks_within25Mb","n.peaks_within50Mb"))
ggplot(x, aes(x=key, y=(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitter(0.15)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  geom_hline(yintercept = 0) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=5))

x <- scaaPerRandomPoints[which(scaaPerRandomPoints$patient_sub=="C552"),26:34]
x <- gather(x, key="key", value="value")
x$key <- factor(x$key, levels = c("propSCAA_100kb","propSCAA_500kb","propSCAA_750kb","propSCAA_1Mb","propSCAA_5Mb","propSCAA_10Mb","propSCAA_20Mb","propSCAA_25Mb","propSCAA_50Mb"))
ggplot(x, aes(x=key, y=(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_jitter(size=3, alpha = .4, shape=21, position=position_jitter(0.15)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  geom_hline(yintercept = 0) +
  theme_custom() +
  theme(legend.position = "top", legend.text = element_text(size=20), legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=5))

# =====
# RECURRING BREAKPOINTS =====
cor.test(breakpoints_dWGS$clusterFrequency, breakpoints_dWGS$clusterClonality)
x <- breakpoints_dWGS[which(breakpoints_dWGS$clusterFrequency>0.5 & breakpoints_dWGS$clusterClonality>0.5),]
length(unique(x$break_base))

for (i in 1:length(unique(x$cluster))) {
  cluster <- unique(x$cluster)[i]
  wd <- x[which(x$cluster==cluster),]
  chr <- substr(unique(wd$break_chr),4,nchar(unique(wd$break_chr)))
  clusterBreak <- unique(wd$centrePoint)
  patients <- unique(wd$patient)
  
  wd2 <- scaa.long[which(scaa.long$chr==chr & 
                            scaa.long$patient %in% patients &
                            (scaa.long$peakStart+250)>=(clusterBreak-1000000) & 
                            (scaa.long$peakStart+250)<=(clusterBreak+1000000)),]
  
  table(wd2$patient, wd2$value)
    
} 