# PREREQUISITIS: load section 1, 4 data 

# TO DO
# hypo: patients with most separation have altered immune/metabolic expression compared to normal
# more outside forces driving change in ATAC

# BOTTOM LINE
# Regions can be separated based on their ATAC data (PC1 & PC2), however loading in PC1/2 are unique to patient
# Regional sep may be based on purity

library("RColorBrewer")
library("pheatmap")
library("nnet")
library("caret")

# DATA QUALITY ASSESSMENT (optional) =====

# Data quality assessment and quality control (i.e. the removal of insufficiently good data) are essential steps of any data analysis. 
# These steps should typically be performed very early in the analysis of a new data set, preceding or in parallel to the differential 
# expression testing.

# We define the term quality as fitness for purpose. Our purpose is the detection of differentially expressed genes, and we are looking in 
# particular for samples whose experimental treatment suffered from an anormality that renders the data points obtained from these particular 
# samples detrimental to our purpose.

# heatmap of the transformed count matrix
select <- order(rowMeans(counts(dds_deseq.perpatient[[i]],normalized=TRUE)),
                decreasing=TRUE)[1:200]
df <- as.data.frame(colData(dds_deseq.perpatient[[i]])[,c("subtissue","purity")])
pheatmap(assay(dds_rlog.perpatient[[i]])[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# HEATMAPS OF SAMPLE-TO-SAMPLE RELATIONSHIPS =====
# heatmap of the sample-to-sample euclidean distances
sample2sample_eucl_hm <- list()
for(i in 1:length(dds_rlog.perpatient)) {
  sampleDists <- dist(t(assay(dds_rlog.perpatient[[i]])))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(dds_rlog.perpatient[[i]]$subtissue, dds_rlog.perpatient[[i]]$purity, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  sample2sample_eucl_hm[[i]] <- pheatmap(sampleDistMatrix,
                                         clustering_distance_rows=sampleDists,
                                         clustering_distance_cols=sampleDists,
                                         col=colors)
}

# heatmap of the sample-to-sample pearson correlation 
sample2sample_pearson_hm <- list()
for(i in 1:length(dds_rlog.perpatient)) {
  mat <- assay(dds_rlog.perpatient[[i]])
  mat_cor <- cor(mat)
  sample2sample_pearson_hm[[i]] <- pheatmap(mat_cor)
}

# Pearson heatmap and assocxiated HC of i=1 show regional separation
# In these comparisons, we used ATAC-seq data without any filtering (except for basic filtering to exclude unreliable signals) 
# to directly compare the nature of the data themselves.
# try: First, the hierarchical clustering method was made with the R heatmap.2 function using ward.D2 and the Pearson correlation 
# as a measure of distance.

# PCA AND PREDICT SUBTISSUE USING TRAINING/TESTING SPLIT (maybe with region removal with low samples) =====

# new pca
pca <- list()
maxCentroidDist <- list()
pca_plot.perpatient <- list()
patient <- list()
mnlm <- list()
mnlm_confu <- list()
mnlm_accur <- list()
mnlm_kappa <- list()
mnlm.cv <- list()
mnlm.cv_confu <- list()
mnlm.cv_accur <- list()

set.seed(123)
for(i in 1:length(dds_rlog.perpatient)) {
  # remove regions with <2 samples
  wd.metaData <- metaData.perpatient[[i]][!metaData.perpatient[[i]]$subtissue %in% names(which(table(metaData.perpatient[[i]]$subtissue) <2)), ]
  wd.rlog <- assay(dds_rlog.perpatient[[i]])[ ,which(colnames(assay(dds_rlog.perpatient[[i]])) %in% rownames(wd.metaData))]
  
  rv <- rowVars(wd.rlog)
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  pca_input <- t(wd.rlog[select,])
  
  # names
  patient[[i]] <- unique(substr(rownames(wd.metaData),1,4))
  
  # save pca on all samples
  pca[[i]] <- prcomp(pca_input)
  percentVar <- round(100 * summary(pca[[i]])$importance[2,])
  df <- cbind(wd.metaData, pca[[i]]$x)
  pca_plot.perpatient[[i]] <- ggplot(df, aes(x=PC1, y=PC2, color = subtissue)) + 
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    stat_ellipse()
  
  # find centroid for PC1/2
  centroid <- as.data.frame(pca[[i]]$x)
  centroid$groups <- substr(rownames(centroid), 6,6)
  centroid <- aggregate(centroid[,1:ncol(centroid)-1], list(Type = centroid$groups), mean)[,c(1:3)]
  
  # find distance between all centroids
  combinations <- expand.grid(first = unique(centroid$Type), second = unique(centroid$Type))
  dist <- list()
  for (j in 1:nrow(combinations)) {
    dist[[j]] <- setNames( c(dist(rbind(centroid[centroid$Type == combinations$first[j],2:3],centroid[centroid$Type == combinations$second[j],2:3]), method = "euclidean")),
                           c(paste(combinations$first[j], combinations$second[j])) )
  }
  
  # keep maximum distance (if there are two with same, keep one)
  dist <- unlist(dist)
  maxCentroidDist[[i]] <- dist[dist == max(dist)][1]
  
  # save the model and associated performance metrics
  df$subtissue <- relevel(as.factor(df$subtissue), ref = unique(df$subtissue)[1])
  mnlm[[i]] <- multinom(subtissue ~ PC1+PC2, data = df, maxit=1000) # build model
  p <- predict(mnlm[[i]], df)
  mnlm_confu[[i]] <- confusionMatrix(p, factor(df$subtissue, levels = levels(p)))
  mnlm_accur[[i]] <- mnlm_confu[[i]]$overall[1]
  mnlm_kappa[[i]] <- mnlm_confu[[i]]$overall[2]
  
  # repeated  cross val
  fit.control <- trainControl(method = "repeatedcv", number = 5, repeats = 10)
  set.seed(123)  
  mnlm.cv[[i]] <- train(subtissue ~ PC1+PC2, data = df, method = "multinom", trControl = fit.control, trace = TRUE)
  mnlm.cv_confu[[i]] <- confusionMatrix(mnlm.cv[[i]])
  mnlm.cv_accur[[i]] <- sum(diag(mnlm.cv_confu[[i]]$table))/100
}

mnlm.accuracies <- data.frame(patient = unlist(patient), accuracy = unlist(mnlm_accur)  , accuracy.cv = unlist(mnlm.cv_accur))
hist(mnlm.accuracies$accuracy)
hist(mnlm.accuracies$accuracy.cv)
mean(mnlm.accuracies$accuracy)
mean(mnlm.accuracies$accuracy.cv)

# drop any patient who couldnt be predicted well (cv<0.7)
names(maxCentroidDist) <- unlist(patient)
names(pca_plot.perpatient) <- unlist(patient)
index <- which(mnlm.accuracies$accuracy.cv <0.7)
maxCentroidDist_sub <- maxCentroidDist[-index]
pca_plot.perpatient_sub <- pca_plot.perpatient[-index]

i=10
pca_plot.perpatient[[i]]
maxCentroidDist_sub[[i]]

# DE per patient
dds <- readRDS('~/Documents/SCAA/Data/EPICC/EPICC_expression/allgenes.dds.ensembl.rds')
res <- results(dds)
vst <- readRDS('~/Documents/SCAA/Data/EPICC/EPICC_expression/allgenes.vsd.ensembl.rds')
geneexp <- as.data.frame(assay(vst))

# CONSIDER PCA LOADINGS =====
top_peaks <- list()
top_loadings <- list()
loadings_plot <- list()

for (i in 1:length(pca)) {
  loadings <- as.data.frame(pca[[i]]$rotation)
  loadings$peak <- rownames(loadings)
  
  top_peaks[[i]] <- loadings %>% 
    select(peak, PC1, PC2) %>% # select only the PCs we are interested in
    pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% # convert to a "long" format
    group_by(PC) %>% # for each PC
    arrange(desc(abs(loading))) %>% # arrange by descending order of loading
    dplyr::slice(1:10) %>% # take the 10 top rows
    pull(peak) %>%# pull the gene column as a vector
    unique() # ensure only unique genes are retained
  
  top_loadings[[i]] <- loadings %>% 
    filter(peak %in% top_peaks[[i]])
  
  loadings_plot[[i]] <- ggplot(data = top_loadings[[i]]) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                 arrow = arrow(length = unit(0.1, "in")),
                 colour = "brown") +
    geom_text(aes(x = PC1, y = PC2, label = peak),
              nudge_y = 0.005, size = 3) +
    scale_x_continuous(expand = c(0.02, 0.02))
}

# the peaks defining a region arent always high in promoters
#list <- list()
#for (i in 1:length(top_peaks)) {
#  list[[i]] <- length(intersect(top_peaks[[i]], PE.peaks)) / length(top_peaks[[i]])
#}
#hist(unlist(list))
#mean(unlist(list))

# Are some differential peaks shared across cohort, with grouping nearby peaks
diff_peaks <- data.frame(location=unlist(top_peaks))
diff_peaks$chr <- sub("\\:.*", "", diff_peaks$location)
diff_peaks$start <- sub(".*[:]([^.]+)[-].*", "\\1", diff_peaks$location)
diff_peaks$stop <- sub(".*-", "", diff_peaks$location)

peak.position <- data.frame(location=rownames(assay(dds_rlog.perpatient[[1]])), bin=1:nrow(assay(dds_rlog.perpatient[[1]])))
diff_peaks$position <- peak.position[match(diff_peaks$location, peak.position$location),2]
diff_peaks$rank <- rank(diff_peaks$position, ties.method = "random")
diff_peaks <- diff_peaks[order(diff_peaks$rank),]

group <- list()
for(i in 1:nrow(diff_peaks)) {
  if (i==1) {
    group[[i]] <- 1
  } else {
    if(diff_peaks$chr[i]==diff_peaks$chr[i-1] & as.numeric(diff_peaks$start[i])-as.numeric(diff_peaks$stop[i-1])<10000) {
      group[[i]] <- group[[i-1]]
    } else {
      group[[i]] <- group[[i-1]]+1
    }
  }
}

diff_peaks$group <- unlist(group)

# any shared pathways?
diff_peaks$geneID <- peakAnn[match(diff_peaks$location, peakAnn_all$peak), 14]
GOanno <- goana(list(diff_peaks=diff_peaks$geneID))
GOanno$P.diff_peaks.adj = p.adjust(GOanno$P.diff_peaks, method = "fdr")
data <- GOanno[which(GOanno$P.diff_peaks.adj <= 0.05),]
data$ID <- rownames(data)
simMatrix_q1 <- calculateSimMatrix(data$ID, orgdb="org.Hs.eg.db", ont = c("BP"), method="Rel")
scores <- setNames(-log10(data$P.diff_peaks.adj), data$ID)
reducedTerms_q1 <- reduceSimMatrix(simMatrix_q1, scores = scores, threshold=0.6, orgdb="org.Hs.eg.db")
whichparent <- data.frame(tapply(reducedTerms_q1$score, reducedTerms_q1$parentTerm, FUN=sum))
include <- rownames(whichparent)[which(whichparent$tapply.reducedTerms_q1.score..reducedTerms_q1.parentTerm..FUN...sum.>0)]
adjusted <- reducedTerms_q1[which(reducedTerms_q1$parentTerm %in% include),]
adj_reducedTerms_q1 <- data.frame(group=adjusted$parentTerm, subgroup=adjusted$term, value=adjusted$score)

jpeg('tempfig1.jpeg', width = 1500, height = 1500)
treemap::treemap(adj_reducedTerms_q1,
                 index=c("group","subgroup"),vSize="value",type="index",
                 fontsize.labels=c(10,5), fontcolor.labels=c("black","white"), fontface.labels=c(2,1),  
                 bg.labels=c("transparent"), align.labels=list(c("center", "center"), c("center", "center")),
                 overlap.labels=1, inflate.labels=F, border.col=c("black","white"), border.lwds=c(3,1),
                 palette = "Set3")
dev.off()

# plot all PCAs
jpeg('tempfig1.jpeg', width = 1500, height = 1500)
plot_grid(pca_plot.perpatient[[1]], pca_plot.perpatient[[2]], pca_plot.perpatient[[3]], pca_plot.perpatient[[4]], pca_plot.perpatient[[5]], 
          pca_plot.perpatient[[6]], pca_plot.perpatient[[7]], pca_plot.perpatient[[8]], pca_plot.perpatient[[9]], pca_plot.perpatient[[10]], 
          pca_plot.perpatient[[11]], pca_plot.perpatient[[12]], pca_plot.perpatient[[13]], pca_plot.perpatient[[14]], pca_plot.perpatient[[15]], 
          pca_plot.perpatient[[16]], pca_plot.perpatient[[17]], pca_plot.perpatient[[18]], pca_plot.perpatient[[19]], pca_plot.perpatient[[20]], 
          pca_plot.perpatient[[21]], pca_plot.perpatient[[22]], pca_plot.perpatient[[23]], pca_plot.perpatient[[24]], 
          nrow = 3, align = "hv")
dev.off()

# HIERARCHICAL CLUSTERING TO COMPARE TO GENETIC PHYLOGENIES =====

hc <- list()
for(i in 1:length(dds_rlog.perpatient)) {
  mat <- assay(dds_rlog.perpatient[[i]])
  dist <- dist(t(mat), method = 'euclidean')
  clusters <- hclust(dist, method = "complete")
  plot(clusters)
}


# SAVE =====
saveRDS(pca_plot.perpatient, "~/Documents/SCAA/Data/pca_plot.perpatient.rds")
saveRDS(top_peaks, "~/Documents/SCAA/Data/top_peaks_cn.norm.peaks.rds")
saveRDS(top_loadings, "~/Documents/SCAA/Data/top_loadings_cn.norm.peaks.rds")
saveRDS(loadings_plot, "~/Documents/SCAA/Data/loadings_plot_cn.norm.peaks.rds")
saveRDS(diff_peaks, "~/Documents/SCAA/Data/diff_peaks.rds")

saveRDS(mnlm, "~/Documents/SCAA/Data/mnlm_cn.norm.peaks.rds")
saveRDS(mnlm_confu, "~/Documents/SCAA/Data/mnlm_confu_cn.norm.peaks.rds")
saveRDS(mnlm_accur, "~/Documents/SCAA/Data/mnlm_accur_cn.norm.peaks.rds")
saveRDS(mnlm_kappa, "~/Documents/SCAA/Data/mnlm_kappa_cn.norm.peaks.rds")
saveRDS(mnlm.cv, "~/Documents/SCAA/Data/mnlm.cv_cn.norm.peaks.rds")
saveRDS(mnlm.cv_confu, "~/Documents/SCAA/Data/mnlm.cv_confu_cn.norm.peaks.rds")
saveRDS(mnlm.cv_accur, "~/Documents/SCAA/Data/mnlm.cv_accur_cn.norm.peaks.rds")
saveRDS(mnlm.accuracies, "~/Documents/SCAA/Data/mnlm.accuracies_cn.norm.peaks.rds")









