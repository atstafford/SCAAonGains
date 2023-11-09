# PREREQUISITIS: load section  data 

# =====
# HYPO: LOSS HAVE PEAK OPEN< AND GAIN HAVE PEAK CLOSED
x <- readxl::read_excel("~/Documents/SCAA/Data/41586_2022_5202_MOESM9_ESM.xlsx")
open <- x[which(x$peak.event_type=="gain"),]
closed <- x[which(x$peak.event_type=="loss"),]

y <- scaa.long[which(scaa.long$peak %in% open$peak.peak & (scaa.long$value==1)),]
z <- scaa.long[which(scaa.long$peak %in% closed$peak.peak & (scaa.long$value==1)),]

table(y$cna)/nrow(y)
table(z$cna)/nrow(z)

# GO OF SCAAS IN HIGH PGA BUT NOT IN LOW PGA =====

# find scaa specific for high GO


# ARE SOME SCAA ASSOCIATED WITH DIPLOID?GAIN?LOSS SEGMENTS =====

# GO all SCAA has some metabolic parent terms
x <- unique((scaa.long$peak[which(scaa.long$value==1)]))
x <- unique(unlist(scaa.long$entrez[which(scaa.long$value==1)]))

# no sig GO terms for SCAAs with sig resulting RNA seq (even if consider peak gain/loss only)
x <- readxl::read_excel("~/Documents/SCAA/Data/41586_2022_5202_MOESM9_ESM.xlsx")
library("org.Hs.eg.db")
x$entrez <- mapIds(org.Hs.eg.db, keys = x$symbol,
                           column = c('ENTREZID'), keytype = 'SYMBOL')
x <- unique(unlist(x$entrez[which(x$p_deseq<=0.05)]))

# no sig GO terms for SCAAs with sig resulting RNA seq, only considering peaks on aneu segments
x <- readxl::read_excel("~/Documents/SCAA/Data/41586_2022_5202_MOESM9_ESM.xlsx")
library("org.Hs.eg.db")
x$entrez <- mapIds(org.Hs.eg.db, keys = x$symbol,
                   column = c('ENTREZID'), keytype = 'SYMBOL')
x <- unique(unlist(x$entrez[which(x$p_deseq<=0.05 & (x$peak.peak %in% scaa.long$peak[which(scaa.long$cna=="gain" | scaa.long$cna=="loss")]))]))
x <- unique(unlist(x$entrez[which(x$p_deseq<=0.05 & (x$peak.peak %in% scaa.long$peak[which(scaa.long$cna=="diploid")]))]))

# manual look at SCAAs on aneuplod segments that influence RNA
x <- x[which(x$p_deseq<=0.05 & (x$peak.peak %in% scaa.long$peak[which(scaa.long$cna=="gain" | scaa.long$cna=="loss")])),]

# GRAMD4 helps induce apop in absense of p53. 12 show loss of enhancer, 10 of which have less GRAMD4.  
# EPB41L3 tumour suppressor. 

GOx <- limma::goana(x)
GOx$p.adj = p.adjust(GOx$P.DE, method = "fdr")
data <- GOx[which(GOx$p.adj <= 0.05),]
data$ID <- rownames(data)
simMatrix <- rrvgo::calculateSimMatrix(data$ID, orgdb="org.Hs.eg.db", ont = c("BP"), method="Rel")
scores <- setNames(-log10(data$p.adj), data$ID)
reducedTerms <- rrvgo::reduceSimMatrix(simMatrix, scores = scores, threshold=0.9, orgdb="org.Hs.eg.db")
parentGO <- data.frame(value=tapply(reducedTerms$score, reducedTerms$parentTerm, FUN=sum))
parentGO$term <- rownames(parentGO)
parentGO$fraction <- ( (parentGO$value/sum(parentGO$value))*100 )
parentGO <- parentGO %>% add_row(term = "other", value=NA, fraction = sum(parentGO$fraction[which(parentGO$fraction<=4)]))
parentGO <- parentGO[which(parentGO$fraction>4),]


# RNA =====
vst <- readRDS("~/Documents/SCAA/Data/EPICC/EPICC_expression/allgenes.vsd.ensembl.rds")
geneexp <- as.data.frame(assay(vst))

z <- readRDS("~/Documents/SCAA/Data/EPICC/EPICC_expression/allgenes.dds.ensembl.rds")
