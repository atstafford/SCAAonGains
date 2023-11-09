# PREREQUISITIS: load section 1, 2 data 

# PULL GLANDS INTO PER PATIENT LIST -----------------------------------------------------------------------------------------------------------------------------------------------

# as theyve already been CN normalised using gland dWGS, these are all glands but may contain some adenomas/normal

# Extract string before ':' in rowname as chr
chr <- sub("\\:.*", "", rownames(atac_cn.deepWGSnorm))
chr <- substr(chr, 4, nchar(chr))

# Extract string inbetween ([^.]+)  ':' and '-' in rownames as start
start <- sub(".*[:]([^.]+)[-].*", "\\1", rownames(atac_cn.deepWGSnorm))

# Extract string after '-'
stop <- sub(".*-", "", rownames(atac_cn.deepWGSnorm))

# Extract identifiers
sample <- unique(colnames(atac_cn.deepWGSnorm))
patient <- sub("\\..*", "", sample)
type <- as.factor(sub("^([[:alpha:]]*).*", "\\1", sub(".*_", "", sample)))

# Check only regions A-D (E is distant normal, F-H is adenoma) c516.c/d were adenoma too
unique(substr(colnames(atac_cn.deepWGSnorm),6,6))
grep("C516.C", colnames(atac_cn.deepWGSnorm))
grep("C516.D", colnames(atac_cn.deepWGSnorm))

# Keep only samples which have purity data (should be all)
purity <- atac_pur$Sample[!is.na(atac_pur$apparent_purity)]
atac_cn.deepWGSnorm <- atac_cn.deepWGSnorm[ ,which(colnames(atac_cn.deepWGSnorm) %in% purity)]

# NORMALIZE FOR LIBRARY SIZE ACROSS ALL PATIENTS USING CN NORMALIZED DATA -------------------------------------------------------------------------------

# input for DESEQ2
counts <- round(atac_cn.deepWGSnorm)
counts <- counts+1 #remove zeros to allow for geomatric calcs
metaData.dWGS.cohort <- data.frame(subtissue = substr(colnames(counts),6,6), 
                                      row.names = colnames(counts))
metaData.dWGS.cohort$purity <- unlist(atac_pur[match(rownames(metaData.dWGS.cohort), atac_pur$Sample),12])
metaData.dWGS.cohort$patient <- substr(rownames(metaData.dWGS.cohort),1,4)
metaData.dWGS.cohort$patient.subtissue <- substr(rownames(metaData.dWGS.cohort),1,6)
dds.dWGS.cohort <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                             colData = metaData.dWGS.cohort, 
                                             design = ~ purity + patient.subtissue)

# run
library(DESeq2)
deseq.dWGS.cohort <- DESeq(dds.dWGS.cohort)
deseqVst.dWGS.cohort <- vst(deseq.dWGS.cohort, blind = F)

# SAVE -----------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(metaData.dWGS.cohort, "~/Documents/SCAA/Data/metaData.dWGS.cohort.rds")
saveRDS(dds.dWGS.cohort, "~/Documents/SCAA/Data/dds.dWGS.cohort.rds")
saveRDS(deseq.dWGS.cohort, "~/Documents/SCAA/Data/deseq.dWGS.cohort.rds")
saveRDS(deseqVst.dWGS.cohort, "~/Documents/SCAA/Data/deseqVst.dWGS.cohort.rds")






