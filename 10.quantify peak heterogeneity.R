# PREREQUISITIS: load section 3 data 

# Idea from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6984959/

# QUANTIFY GENE-SPECIFIC PATIENT-SPECFIC PEAK-ITH =====

# For a given tumour, the standard deviation of expression values for a particular gene across tumour regions was 
# calculated yielding a gene-specific, patient-specific measure of peak-ITH (σg,p). 
# This was repeated for all genes, then all tumours, generating a matrix of σg,p values
# The peaks with high ITH per patient may not be in PC1/2 as pca isnt guided by groups 

# dont use mean per region because hetero within region. Instead consider SD across all glands and samples
patients <- unique(substr(colnames(assay(dds_vst.dWGSperpatient)),1,4))

stdev.per.region <- list()
names <- list()
for (i in 1:length(patients)) {
  wd <- assay(dds_vst.dWGSperpatient)[ ,which(substr(colnames(assay(dds_vst.dWGSperpatient)),1,4) == patients[i])]
  stdev.per.region[[i]] <- data.frame(rowSds(wd))
  names(stdev.per.region[[i]]) <- patients[[i]]
  rownames(stdev.per.region[[i]]) <- rownames(assay(dds_vst.dWGSperpatient))
}
stdev.per.region <- do.call("cbind",stdev.per.region)

hist(stdev.per.region$C516)

# CALCULATE PEAK-ITH PER PATIENT =====
# Patient-wise peak-ITH values are summarised as the average (median) value per tumour across all expressed genes (σp)
patientwise.pITH <- data.frame(pITH = colMedians(as.matrix(stdev.per.region)))
rownames(patientwise.pITH) <- patients

# CALCULATE PEAK-ITH PER REGION =====
# Region-wise peak-ITH values are summarised as the average (median) value per peak across all tumours in the cohort (σg)
regionwise.pITH <- data.frame(pITH = rowMedians(as.matrix(stdev.per.region)))
rownames(regionwise.pITH) <- rownames(assay(dds_vst.dWGSperpatient))
hist(regionwise.pITH$pITH)

# CALCULATE PEAK-IPH PER COHORT =====
# inter-tumour heterogeneity is derived for each peak by randomly sampling one region per patient and taking the st dev across the 
# resulting single-biopsy cohort, then repeating this process 10 times to take the average score across iteration
cohortwise.pIPH <- list()
for (j in 1:100) {
  print(j)
  sample <- list()
  for (i in 1:length(patients)) {
    set.seed(i*j)
    wd <- assay(dds_vst.dWGSperpatient)[ ,which(substr(colnames(assay(dds_vst.dWGSperpatient)),1,4) == patients[i])] 
    sample[[i]] <- wd[ ,sample(1:ncol(wd), 1)]  
  }
  sample <- do.call("cbind", sample)
  cohortwise.pIPH[[j]] <- data.frame(rowSds(sample))
}

cohortwise.pIPH <- do.call("cbind", cohortwise.pIPH)
cohortwise.pIPH <- data.frame(pIPH = rowMedians(as.matrix(cohortwise.pIPH)))
rownames(cohortwise.pIPH) <- rownames(assay(dds_vst.dWGSperpatient))
hist(cohortwise.pIPH$pIPH)

# SAVE =====
saveRDS(stdev.per.region, "~/Documents/SCAA/Data/stdev.per.region_cn.norm.peaks.rds")
saveRDS(patientwise.pITH, "~/Documents/SCAA/Data/patientwise.pITH_cn.norm.peaks.rds")
saveRDS(regionwise.pITH, "~/Documents/SCAA/Data/regionwise.pITH_cn.norm.peaks.rds")
saveRDS(cohortwise.pIPH, "~/Documents/SCAA/Data/cohortwise.pIPH_cn.norm.peaks.rds")
