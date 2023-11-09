# PREREQUISITIS: load section 4 data 

# epicc showed that dirving genetic alterations were sort of rare, but loads of SCAAs in cancer genes were recurrent and recurrent SCAAs in other gene not ass with cancer (yet...)
# most scaas are clonal


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

# promoters / non promoters peaks across G-bands ----
x <- scaa.long[,c(1:4,23,24)]
x <- x[!duplicated(x),]

x <- data.frame("Promoter/enhancer peak"= c(nrow(x[x$geiStain=="gneg" & (x$promoter=="TRUE" | x$enhancer=="TRUE"),]),
                        nrow(x[x$geiStain=="gvar" &  (x$promoter=="TRUE" | x$enhancer=="TRUE"),]),
                        nrow(x[x$geiStain=="acen" &  (x$promoter=="TRUE" | x$enhancer=="TRUE"),]),
                        nrow(x[x$geiStain=="gpos25" &  (x$promoter=="TRUE" | x$enhancer=="TRUE"),]),
                        nrow(x[x$geiStain=="gpos50" &  (x$promoter=="TRUE" | x$enhancer=="TRUE"),]),
                        nrow(x[x$geiStain=="gpos75" &  (x$promoter=="TRUE" | x$enhancer=="TRUE"),]),
                        nrow(x[x$geiStain=="gpos100" &  (x$promoter=="TRUE" | x$enhancer=="TRUE"),]),
                        nrow(x[x$geiStain=="no data" &  (x$promoter=="TRUE" | x$enhancer=="TRUE"),])),
                "Non-promoter/enhancer peak"= c(nrow(x[x$geiStain=="gneg" &  (x$promoter=="FALSE" & x$enhancer=="FALSE"),]),
                           nrow(x[x$geiStain=="gvar" & (x$promoter=="FALSE" & x$enhancer=="FALSE"),]),
                           nrow(x[x$geiStain=="acen" & (x$promoter=="FALSE" & x$enhancer=="FALSE"),]),
                           nrow(x[x$geiStain=="gpos25" & (x$promoter=="FALSE" & x$enhancer=="FALSE"),]),
                           nrow(x[x$geiStain=="gpos50" & (x$promoter=="FALSE" & x$enhancer=="FALSE"),]),
                           nrow(x[x$geiStain=="gpos75" & (x$promoter=="FALSE" & x$enhancer=="FALSE"),]),
                           nrow(x[x$geiStain=="gpos100" & (x$promoter=="FALSE" & x$enhancer=="FALSE"),]),
                           nrow(x[x$geiStain=="no data" & (x$promoter=="FALSE" & x$enhancer=="FALSE"),])))
sum(x)
rownames(x) <- c("gneg","gvar","acen","gpos25","gpos50","gpos75","gpos100", "no data")
x <- x[which(rownames(x)!="no data"),]

chisq <- chisq.test(x)
enrich <- chisq$observed/chisq$expected

x$Stain <- rownames(x)
x <- gather(x, key = "key", value="Count", -Stain)
x$fraction <- ifelse(x$key=="Promoter.enhancer.peak", 
                     round(100*(x$Count/sum(x$Count[which(x$key=="Promoter.enhancer.peak")])),1),
                     round(100*(x$Count/sum(x$Count[which(x$key=="Non.promoter.enhancer.peak")])),1) )
sum(x$fraction)

library(ggrepel)

#jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=key, y=fraction, fill=Stain)) + 
  geom_bar(position="stack", stat="identity", alpha=0.5, width = 0.4) +
  geom_text_repel(aes(label = paste(Stain, paste(fraction,"%", sep = ""), sep = ": ")), position = position_stack(vjust = 0.5),
                  vjust=1, size = 7, box.padding = 0.5) +
  scale_x_discrete(name=NULL, labels = c("Not promoter/\nenhancer peak","Promoter/\nenhancer peak")) +
  scale_fill_brewer(palette = "Dark2") +
  ylab("Percent") +
  theme_custom() +
  theme(
        legend.text = element_text(size=28), legend.title = element_blank(), legend.position = "top")
dev.off()



# SCAA / non SCAA peaks across G-bands ----
library(ggrepel)
library(ggbreak) 
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

chisq <- chisq.test(x)
enrich <- chisq$observed/chisq$expected

x$Stain <- rownames(x)
x <- gather(x, key = "key", value="Count", -Stain)
x$fraction <- ifelse(x$key=="SCAA_peak", 
                     round(100*(x$Count/sum(x$Count[which(x$key=="SCAA_peak")])),1),
                     round(100*(x$Count/sum(x$Count[which(x$key=="NonSCAA_peak")])),1) )
sum(x$fraction)

#jpeg('tempfig.jpeg', height = (3*37.795*6), width = (3*37.795*6))
jpeg('tempfig.jpeg', width = (20), height = (20), units = "cm", res = 300)
ggplot(x, aes(x=key, y=fraction, fill=Stain)) + 
  geom_bar(position="stack", stat="identity", alpha=0.5, width = 0.4) +
  geom_text_repel(aes(label = paste(Stain, paste(fraction,"%", sep = ""), sep = ": ")), position = position_stack(vjust = 0.5),
                  vjust=1, size = 7, box.padding = 0.5) +
  scale_x_discrete(name=NULL, labels = c("Not SCAA peak","SCAA peak")) +
  scale_fill_brewer(palette = "Dark2") +
  ylab("Percent") +
  theme_custom() +
  theme(
    legend.text = element_text(size=28), legend.title = element_blank(), legend.position = "top")
dev.off()



ggplot(x, aes(y=Count, x=Stain, fill= key)) + 
  geom_bar(position="dodge", stat="identity",color="black") +
  geom_text(aes(label=Count), hjust = 0.5, vjust=1, position = position_dodge(width = 1), size=5) +
  scale_y_continuous(limits = c(0, 880000), breaks = c(0,2000,4000,6000,8000,50000,125000,200000,860000,880000)) +
  scale_y_break(c(8000, 20000),  scales = 0.5) + 
  scale_y_break(c(200000, 850000),  scales = 0.3) + 
  scale_fill_manual(labels = c( "peak","SCAA"), values=c("#d8b365","#5ab4ac")) +
  theme_custom() +
  theme(
    legend.text = element_text(size=20), legend.title = element_blank(), legend.position = "top",
    axis.title.y = element_blank(), axis.text.y.right = element_blank(), 
    axis.title.x = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1))
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


itrack <- IdeogramTrack(genome = "hg38", chromosome = "chr3")
plotTracks(itrack, from=94004246, to=98593978)

itrack <- IdeogramTrack(genome = "hg38", chromosome = "chr16")
plotTracks(itrack, from=46554286, to=46973978)

itrack <- IdeogramTrack(genome = "hg38", chromosome = "chr19")
plotTracks(itrack, from=19900625, to=31479397)

itrack <- IdeogramTrack(genome = "hg38", chromosome = "chr13")
plotTracks(itrack, from=16010054, to=16201673)

itrack <- IdeogramTrack(genome = "hg38", chromosome = "chr22")
plotTracks(itrack, from=11973990, to=11974490)

itrack <- IdeogramTrack(genome = "hg38", chromosome = "chr21")
plotTracks(itrack, from=8765041, to=10616677)


z <- scaa.long[which(scaa.long$geiStain=="gvar"),]

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

# do patient with the gvar peaks have WGD?
scaa.long$ploidy <- ploidy.perpatient[match(scaa.long$patient, ploidy.perpatient$patient),2]
scaa.long$WGD <- ifelse(scaa.long$ploidy>2.5, TRUE, FALSE)

x <- data.frame(WGD= c(nrow(scaa.long[scaa.long$geiStain=="gneg" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gvar" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="acen" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos25" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos50" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos75" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gpos100" & scaa.long$WGD==TRUE & scaa.long$value==1,])),
                diploid= c(nrow(scaa.long[scaa.long$geiStain=="gneg" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                           nrow(scaa.long[scaa.long$geiStain=="gvar" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                           nrow(scaa.long[scaa.long$geiStain=="acen" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos25" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos50" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos75" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos100" & scaa.long$WGD==FALSE & scaa.long$value==1,])))

x <- data.frame(SCAA= c(nrow(scaa.long[scaa.long$geiStain=="gneg" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                       nrow(scaa.long[scaa.long$geiStain=="gvar" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                       nrow(scaa.long[scaa.long$geiStain=="acen" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                       nrow(scaa.long[scaa.long$geiStain=="gpos25" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                       nrow(scaa.long[scaa.long$geiStain=="gpos50" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                       nrow(scaa.long[scaa.long$geiStain=="gpos75" & scaa.long$WGD==FALSE & scaa.long$value==1,]),
                       nrow(scaa.long[scaa.long$geiStain=="gpos100" & scaa.long$WGD==FALSE & scaa.long$value==1,])),
                notSCAA= c(nrow(scaa.long[scaa.long$geiStain=="gneg" & scaa.long$WGD==FALSE & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="gvar" & scaa.long$WGD==FALSE & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="acen" & scaa.long$WGD==FALSE & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos25" & scaa.long$WGD==FALSE & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos50" & scaa.long$WGD==FALSE & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos75" & scaa.long$WGD==FALSE & scaa.long$value==0,]),
                           nrow(scaa.long[scaa.long$geiStain=="gpos100" & scaa.long$WGD==FALSE & scaa.long$value==0,])))

x <- data.frame(SCAA= c(nrow(scaa.long[scaa.long$geiStain=="gneg" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                      nrow(scaa.long[scaa.long$geiStain=="gvar" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                      nrow(scaa.long[scaa.long$geiStain=="acen" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                      nrow(scaa.long[scaa.long$geiStain=="gpos25" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                      nrow(scaa.long[scaa.long$geiStain=="gpos50" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                      nrow(scaa.long[scaa.long$geiStain=="gpos75" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                      nrow(scaa.long[scaa.long$geiStain=="gpos100" & scaa.long$WGD==TRUE & scaa.long$value==1,])),
              notSCAA= c(nrow(scaa.long[scaa.long$geiStain=="gneg" & scaa.long$WGD==TRUE & scaa.long$value==0,]),
                         nrow(scaa.long[scaa.long$geiStain=="gvar" & scaa.long$WGD==TRUE & scaa.long$value==0,]),
                         nrow(scaa.long[scaa.long$geiStain=="acen" & scaa.long$WGD==TRUE & scaa.long$value==0,]),
                         nrow(scaa.long[scaa.long$geiStain=="gpos25" & scaa.long$WGD==TRUE & scaa.long$value==0,]),
                         nrow(scaa.long[scaa.long$geiStain=="gpos50" & scaa.long$WGD==TRUE & scaa.long$value==0,]),
                         nrow(scaa.long[scaa.long$geiStain=="gpos75" & scaa.long$WGD==TRUE & scaa.long$value==0,]),
                         nrow(scaa.long[scaa.long$geiStain=="gpos100" & scaa.long$WGD==TRUE & scaa.long$value==0,])))
rownames(x) <- c("gneg","gvar","acen","gpos25","gpos50","gpos75","gpos100")

chisq <- chisq.test(x)
enrich <- chisq$observed/chisq$expected


x <- data.frame(Gvar_T= c(nrow(scaa.long[scaa.long$geiStain=="gvar" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                        nrow(scaa.long[scaa.long$geiStain=="gvar" & scaa.long$WGD==FALSE & scaa.long$value==1,])),
                Gvar_F= c(nrow(scaa.long[scaa.long$geiStain!="gvar" & scaa.long$WGD==TRUE & scaa.long$value==1,]),
                           nrow(scaa.long[scaa.long$geiStain!="gvar" & scaa.long$WGD==FALSE & scaa.long$value==1,])))
rownames(x) <- c("WGD","diploid")
chisq <- chisq.test(x)
enrich <- chisq$observed/chisq$expected

# SCAA ENRICHMENT BY TELOMERE DISTANCE =====
library(ggpubr)
my_comparisons <- list( c("0-1", "1-25"), c("0-1", "25-50"), c("0-1", "50-75"), c("0-1", "75-100"), c("0-1", "100-150"),
                        c("1-25", "25-50"), c("1-25", "50-75"), c("1-25", "75-100"), c("1-25", "100-150"),
                        c("25-50", "50-75"), c("25-50", "75-100"), c("25-50", "100-150"),
                        c("50-75", "75-100"), c("50-75", "100-150"),
                        c("75-100", "100-150"))

#jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
jpeg('tempfig.jpeg', width = (60), height = (40), units="cm", res=300)
ggplot2::ggplot(data=scaa.long, aes(x=as.factor(value), y=(dist.telo), fill=as.factor(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to telomere (bp)") +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  #facet_wrap(~chr) +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 150000000, label.x=1.4, size=9, method = "wilcox.test") +
  theme(legend.position = "none", legend.text = element_text(size=20), legend.title = element_blank(),
         axis.title.x = element_blank())
dev.off()


t.test(scaa.long$dist.telo[which(scaa.long$value==0)], scaa.long$dist.telo[(scaa.long$value==1)])

# by cna
jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
ggplot2::ggplot(data=scaa.long[!is.na(scaa.long$cna),], aes(x=as.factor(cna), y=(dist.telo), fill=as.factor(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to telomere (bp)") +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 125000000, label.x=1.4, size=9, method = "wilcox.test") +
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

#jpeg('tempfig.jpeg', width = (3*37.795*5), height = (3*37.795*5))
jpeg('tempfig.jpeg', width = (20), height = (20), units="cm", res=300)
ggplot2::ggplot(data=scaa.long, aes(x=as.factor(value), y=(dist.centroPC), fill=as.factor(value))) +
  geom_violin(trim = F, alpha = .4, position=position_dodge(0.6)) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE,position=position_dodge(0.6)) +
  ylab("Distance to centromere (bp)") +
  scale_x_discrete(labels = c("1" = "SCAA", "0" = "not SCAA")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_custom() +
  stat_compare_means(aes(group = as.factor(value)), label="p.format", label.y = 1, label.x=1.4, size=9, method = "wilcox.test") +
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

