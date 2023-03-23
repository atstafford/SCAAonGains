# GLM MODEL SCAAS =====

# model g-band (univariate)
model <- glm(as.factor(value) ~ geiStain, data = scaa.long[scaa.long$geiStain!="no data",], family="binomial")
summary(model)

# model BP distance (univariate)
hist(scaa.long$distancetoBP)
hist(scaa.long$distancetoBP[scaa.long$value==1])
hist(scaa.long$distancetoBP[scaa.long$value==0])
ks.test(scaa.long$distancetoBP[scaa.long$value==1], scaa.long$distancetoBP[scaa.long$value==0])
plot(ecdf(x = scaa.long$distancetoBP[scaa.long$value==0]), main = "ECDF of x and y")
lines(ecdf(x = scaa.long$distancetoBP[scaa.long$value==1]), col = 2)

ggplot(scaa.long, aes(x=as.factor(value), y=(distancetoBP)))+
  geom_boxplot()

model <- glm(as.factor(value) ~ distancetoBP, data = scaa.long, family="binomial")
summary(model)

# model telomere/centromere (univariate)
hist(scaa.long$dist.centro[scaa.long$value==1])
hist(scaa.long$dist.centro[scaa.long$value==0])
ks.test(scaa.long$dist.centro[scaa.long$value==1], scaa.long$dist.centro[scaa.long$value==0])
plot(ecdf(x = scaa.long$dist.centro[scaa.long$value==0]), main = "ECDF of x and y")
lines(ecdf(x = scaa.long$dist.centro[scaa.long$value==1]), col = 2)

ggplot(scaa.long, aes(x=as.factor(value), y=(dist.centro)))+
  geom_boxplot()

model <- glm(as.factor(value) ~ dist.centro, data = scaa.long, family="binomial")
summary(model)

hist(scaa.long$dist.telo[scaa.long$value==1])
hist(scaa.long$dist.telo[scaa.long$value==0])
ks.test(scaa.long$dist.telo[scaa.long$value==1], scaa.long$dist.telo[scaa.long$value==0])
plot(ecdf(x = scaa.long$dist.telo[scaa.long$value==0]), main = "ECDF of x and y")
lines(ecdf(x = scaa.long$dist.telo[scaa.long$value==1]), col = 2)

ggplot(scaa.long, aes(x=as.factor(value), y=(dist.telo)))+
  geom_boxplot()

model <- glm(as.factor(value) ~ dist.telo, data = scaa.long, family="binomial")
summary(model)

# model cna/clonality (univariate)
model <- glm(as.factor(value) ~ fracGain, data = scaa.long, family="binomial")
summary(model)
model <- glm(as.factor(value) ~ fracLoss, data = scaa.long, family="binomial")
summary(model)
model <- glm(as.factor(value) ~ fracDip, data = scaa.long, family="binomial")
summary(model)
model <- glm(as.factor(value) ~ fracGain + fracLoss + fracDip, data = scaa.long, family="binomial")
summary(model)


model <- glm(as.factor(value) ~ cna, data = scaa.long, family="binomial")
summary(model)
model <- glm(as.factor(value) ~ mainClonality, data = scaa.long, family="binomial")
summary(model)
model <- glm(as.factor(value) ~ cna + mainClonality, data = scaa.long, family="binomial")
summary(model)

# model multivariate g-band + location
model <- glm(as.factor(value) ~ geiStain + dist.telo + dist.centro, data = scaa.long[scaa.long$geiStain!="no data",], family="binomial")
summary(model)

model <- glm(as.factor(value) ~ geiStain + dist.telo + dist.centro + distancetoBP, data = scaa.long[scaa.long$geiStain!="no data",], family="binomial")
summary(model)

model <- glm(as.factor(value) ~ geiStain + dist.telo + dist.centro + distancetoBP + cna + mainClonality, data = scaa.long[scaa.long$geiStain!="no data",], family="binomial")
summary(model)

model <- glm(as.factor(value) ~ dist.telo + dist.centro + distancetoBP + cna + mainClonality, data = scaa.long, family="binomial")
summary(model)

# multivariate based model diagnostics
library(tidyr)
library(dplyr)
library(ggplot2)
data <- scaa.long[, which(colnames(scaa.long) %in% c("value","distancetoBP","geiStain","dist.telo","dist.centro","cna","mainClonality"))]
data <- na.omit(data)
data <- data[data$geiStain!="no data",]
data$dist.telo2 <- data$dist.telo^2
model <- glm(as.factor(value) ~ ., data = data, family="binomial")
summary(model)
probabilities <- predict(model, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "SCAA", "noSCAA")

mydata <- data[colnames(data) %in% c("distancetoBP","dist.telo","dist.centro","mainClonality","dist.telo2")]
mydata <- data[colnames(data) %in% c("distancetoBP","dist.telo","dist.centro","mainClonality")]
mydata <- mydata %>%
  mutate(logit = log(probabilities/(1-probabilities))) %>%
  gather(key = "predictors", value = "predictor.value", -logit)

x <- sample_n(mydata, 1000)
jpeg('tempfig.jpeg', width = (3*37.795*7), height = (3*37.795*7))
ggplot(x, aes(logit, predictor.value)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw() +  
  facet_wrap(~predictors, scales = "free_y")
dev.off()

# CARET
data <- scaa.long[, which(colnames(scaa.long) %in% c("value","distancetoBP","geiStain","dist.telo","dist.centro","cna","mainClonality"))]
data <- na.omit(data)
data <- data[data$geiStain!="no data",]
data$value <- as.factor(data$value)

library(caret)
set.seed(430)
part <- createDataPartition(data$value, p = 0.75, list = FALSE)
train <- data[part, ]
test <- data[-part, ]

glm_mod = train(
  form = value ~ .,
  data = train,
  trControl = trainControl(method = "cv", number = 5),
  method = "glm",
  family = "binomial"
)

glm_mod
glm_mod$results
glm_mod$finalModel
summary(glm_mod)

calc_acc = function(actual, predicted) {
  mean(actual == predicted)
}

head(predict(glm_mod, newdata = test))

calc_acc(actual = test$value,
         predicted = predict(glm_mod, newdata = test))

x <- (predict(glm_mod, newdata = train, type = "prob"))


# DATA =====
library(caret)
library(ROSE)
library(rpart)
library(naivebayes)

#Will need to use balancing
data <- scaa.long[, which(colnames(scaa.long) %in% c("value","distancetoBP","geiStain","dist.telo","dist.centro","cna","mainClonality"))]
data <- na.omit(data)
data <- data[data$geiStain!="no data",]
data$value <- as.factor(data$value)
data$cna <- as.factor(data$cna)
data$geiStain <- as.factor(data$geiStain)

# split
set.seed(430)
part <- createDataPartition(data$value, p = 0.75, list = FALSE)
train <- data[part, ]
test <- data[-part, ]

# balance
data_balanced_over <- ovun.sample(value ~ ., data = train, method = "over",N = 2100000)$data #over
table(data_balanced_over$value)
prop.table(table(data_balanced_over$value))

data_balanced_under <- ovun.sample(value ~ ., data = train, method = "under",N = 44796*2)$data #under
table(data_balanced_under$value)
prop.table(table(data_balanced_under$value))

data_balanced_both <- ovun.sample(value ~ ., data = train, method = "both",N = 1000000, seed=1)$data #both
table(data_balanced_both$value)
prop.table(table(data_balanced_both$value))

data.rose <- ROSE(value ~ ., data = train, seed = 1)$data #synthetic
table(data.rose$value)
prop.table(table(data.rose$value))

# DECISION TREE =====
# As the rules are learned sequentially, from trunk to leaf, a decision tree requires high quality, 
# clean data from the outset of training, or the branches may become over-fitted or skewed.
# assumption of independence among predictors.

# imbalanced
table(train$value)
prop.table(table(train$value))

treeimb <- rpart(value ~ ., data = train)
pred.treeimb <- predict(treeimb, newdata = test)

accuracy.meas(test$value, pred.treeimb[,2]) # preicion 1=low false pos, recall 1=low false negs

roc.curve(test$value, pred.treeimb[,2], plotit = F) # AUC

#build decision tree models
tree.rose <- rpart(value ~ ., data = data.rose) 
tree.over <- rpart(value ~ ., data = data_balanced_over)
tree.under <- rpart(value ~ ., data = data_balanced_under)
tree.both <- rpart(value ~ ., data = data_balanced_both)

#make predictions on unseen data
pred.tree.rose <- predict(tree.rose, newdata = test)
pred.tree.over <- predict(tree.over, newdata = test)
pred.tree.under <- predict(tree.under, newdata = test)
pred.tree.both <- predict(tree.both, newdata = test)

#AUC
roc.curve(test$value, pred.tree.rose[,2])
roc.curve(test$value, pred.tree.over[,2])
roc.curve(test$value, pred.tree.under[,2])
roc.curve(test$value, pred.tree.both[,2])

# check the model accuracy using holdout and bagging method. This helps us to ensure that our resultant predictions doesn’t suffer from high variance.
ROSE.holdout <- ROSE.eval(cls ~ ., data = train, learner = rpart, method.assess = "holdout", extr.pred = function(obj)obj[,2], seed = 1)

# NAIVE BAYES =====

# assumption of independence among predictors.
# useful for very large data sets

model <- naive_bayes(value ~ ., data = data_balanced_over, usekernel = T) 
p <- predict(model, test)
tab <- table(p, test$value)
1 - sum(diag(tab)) / sum(tab) #misclassification rate is high (0.5)

model <- naive_bayes(value ~ ., data = data_balanced_under, usekernel = T) 
p <- predict(model, test)
tab <- table(p, test$value)
1 - sum(diag(tab)) / sum(tab) #misclassification rate is high (0.8)

model <- naive_bayes(value ~ ., data = data_balanced_both, usekernel = T) 
p <- predict(model, test)
tab <- table(p, test$value)
1 - sum(diag(tab)) / sum(tab) #misclassification rate is high (0.7)

model <- naive_bayes(value ~ ., data = data.rose, usekernel = T) 
p <- predict(model, test)
tab <- table(p, test$value)
1 - sum(diag(tab)) / sum(tab) #misclassification rate is low (0.12)

# =====
# cytoband data
hg38.chromosomeBand <- read.csv("~/Documents/SCAA/Data/UCSC_search/hg38.chromosomeBands.csv")

# centromere data
library(data.table)
centromere <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", 
                    col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
centromere <- centromere[ , .(length = sum(chromEnd - chromStart)), 
                          by = .(chrom, arm = substring(name, 1, 1)) ]
centromere <- centromere[centromere$arm=="p"]

# telomere data
library("DescTools")
library(rtracklayer)
hg38.length <- SeqinfoForUCSCGenome("hg38")
hg38.length <- data.frame("chr"=substr(hg38.length@seqnames,4,nchar(hg38.length@seqnames)), "length"=hg38.length@seqlengths)
hg38.length <- hg38.length[c(1:22),]

# TSS data
mySession = browserSession("UCSC")
genome(mySession) <- "hg38"
range <- GRanges(paste("chr", scaa.bins$chr[c(1:100)], sep = ""), IRanges(scaa.bins$peakStart[c(1:100)], scaa.bins$peakStop[c(1:100)]))
track.names <- trackNames(ucscTableQuery(mySession))
tableNames(ucscTableQuery(mySession, track="geneHancer"))
TSS <- getTable(ucscTableQuery(mySession, track="geneHancer",
                               range=range, table="geneHancerGenesDoubleElite"))

# collect peak locations
peak_analysis <- readRDS("~/Documents/SCAA/Data/EPICC/SCAA/analysis_results.rds")
scaa <- peak_analysis$sig
scaa <- scaa[ ,which(substr(colnames(scaa),6,6) == "p")]

scaa.bins <- data.frame(chr = substr(sub("\\:.*", "", rownames(scaa)),4,nchar(sub("\\:.*", "", rownames(scaa)))),
                        peakStart = as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", rownames(scaa))),
                        peakStop = as.numeric(sub(".*-", "", rownames(scaa))) )

# assign to each peak
stain <- list()
dist.centro <- list()
dist.telo <- list()

for (i in 1:nrow(scaa.bins)) {
  print(paste(i,"/",nrow(scaa.bins)))
  chr <- paste("chr",scaa.bins$chr[i], sep = "")
  start <- as.numeric(scaa.bins$peakStart[i])
  stop <- as.numeric(scaa.bins$peakStop[i])
  
  # cytoband per peak
  stain[i] <- "no data"
  wd <- hg38.chromosomeBand[hg38.chromosomeBand$X.chrom==chr,]
  for (j in 1:nrow(wd)) {
    if( nrow(wd)==0 ) {
      stain[i] <- "no data"
      break
    }
    if( start>=wd$chromStart[j] & stop<=wd$chromEnd[j] ) {
      stain[i] <- wd$gieStain[j]
      break
    } else {
      next
    }
  }
  
  # centromere per peak
  wd <- centromere[centromere$chrom == chr,]
  dist.centro[[i]] <- min(abs(stop - wd$length), abs(start - wd$length))
  
  # telomere per peak
  wd <- hg38.length[hg38.length$chr == substr(chr, 4, nchar(chr)),]
  dist.telo[[i]] <- min(abs(start - 1), abs(stop - wd$length))
  
}
scaa.bins$geiStain <- unlist(stain)
scaa.bins$dist.centro <- unlist(dist.centro)
scaa.bins$dist.telo <- unlist(dist.telo)

# add TSS per peak
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9450098/ says promoter is 3kb +/- from TSS but didnt make differe ce
TSS.n <- list()
for (i in 1:nrow(scaa.bins)) {
  print(paste(i,"/",nrow(scaa.bins)))
  chr <- paste("chr", scaa.bins$chr[i], sep = "")
  start <- scaa.bins$peakStart[i]-3000
  stop <- scaa.bins$peakStop[i]+3000
  
  wd <- TSS[TSS$chrom==chr,]
  for (j in 1:nrow(wd)) {
    count <- 0
    if (nrow(wd)==0) {
      break
    }
    if( (wd$chromStart[j])>=(start) & (wd$chromEnd[j])<=(stop)) {
      count <- count + 1
    } else {
      next
    }
  }
  TSS.n[[i]] <- count
}
scaa.bins$TSS.n <- unlist(TSS.n)

# CpG data
library(rtracklayer)
mySession = browserSession("UCSC")
genome(mySession) <- "hg38"
range <- GRanges(paste("chr", scaa.bins$chr[c(1:10)], sep = ""), IRanges(scaa.bins$peakStart[c(1:10)], scaa.bins$peakStop[c(1:10)]))
CpG <- getTable(ucscTableQuery(mySession, track="cpgIslandExt",
                               range=range, table="cpgIslandExt"))

# add CpG island coverage per peak
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9450098/ says promoter is 3kb +/- from TSS but didnt make differe ce
CpG.n <- list()
for (i in 1:nrow(scaa.bins)) {
  print(paste(i,"/",nrow(scaa.bins)))
  chr <- paste("chr", scaa.bins$chr[i], sep = "")
  start <- scaa.bins$peakStart[i]-3000
  stop <- scaa.bins$peakStop[i]+3000
  
  wd <- CpG[CpG$chrom==chr,]
  for (j in 1:nrow(wd)) {
    count <- 0
    if (nrow(wd)==0) {
      break
    }
    if( (wd$chromStart[j])>=(start) & (wd$chromEnd[j])<=(stop)) {
      count <- count + 1
    } else {
      next
    }
  }
  CpG.n[[i]] <- count
}
scaa.bins$CpG.n <- unlist(CpG.n)



scaa.long <- cbind(scaa.bins, scaa)
scaa.long <- gather(scaa.long, key="patient", value = "value", -chr, -peakStart, -peakStop, -geiStain, -dist.centro, -dist.telo, -TSS.n, -CpG.n)
scaa.long$patient <- substr(scaa.long$patient, 1, 4)



# add segment per peak (needs to be patient specific)
cnaStart <- list()
cnaStop <- list()
segSize <- list()
n.samples <- list()
fracDip <- list()
fracLoss <- list()
fracGain <- list()
distancetoBP <- list()
for (i in 1:nrow(scaa.long)) {
  print(paste(i,"/",nrow(scaa.long)))
  patient <- scaa.long$patient[i]
  
  chr <- scaa.long$chr[i]
  peakStart <- as.numeric(sub(".*[:]([^.]+)[-].*", "\\1", scaa.long$peakStart[i]))
  peakStop <- as.numeric(sub(".*-", "", scaa.long$peakStop[i]))
  
  dwgs.wd <- dwgs.perpatient[[patient]]
  dwgs.wd <- dwgs.wd[dwgs.wd$chr == chr, ]
  
  cnaStart[[i]] <- "no data"
  cnaStop[[i]] <- "no data"
  segSize[[i]] <- "no data"
  n.samples[[i]] <- "no data"
  fracDip[[i]] <- "no data"
  fracLoss[[i]] <- "no data"
  fracGain[[i]] <- "no data"
  distancetoBP[[i]] <- "no data"
  
  if(nrow(dwgs.wd)==0) { # no chromosome data
    next
  }
  
  for (k in 1:nrow(dwgs.wd)) {
    if ( DescTools::Overlap(c(peakStart, peakStop), c(dwgs.wd$start[k], dwgs.wd$stop[k])) >= 250 ) { # choose segment which majority of peak sits on
      cnaStart[[i]] <- dwgs.wd$start[k]
      cnaStop[[i]] <- dwgs.wd$stop[k]
      segSize[[i]] <- dwgs.wd$stop[k] - dwgs.wd$start[k]
      n.samples[[i]] <- ncol(dwgs.wd)-3
      fracDip[[i]] <- round( (sum(dwgs.wd[k, -c(1:3)]==2, na.rm = T)) / (ncol(dwgs.wd)-3), 2 )
      fracLoss[[i]] <- round( (sum(dwgs.wd[k, -c(1:3)]<2, na.rm = T)) / (ncol(dwgs.wd)-3), 2 )
      fracGain[[i]] <- round( (sum(dwgs.wd[k, -c(1:3)]>2, na.rm = T)) / (ncol(dwgs.wd)-3), 2 )
      peakMid <- peakStart + 250
      distancetoBP[[i]] <- min((peakMid-dwgs.wd$start[k]), (dwgs.wd$stop[k]-peakMid))
      break
    } else {
      next
    }
  }
}

scaa.long$cnaStart <- as.numeric(unlist(cnaStart))
scaa.long$cnaStop <- as.numeric(unlist(cnaStop))
scaa.long$segSize <- as.numeric(unlist(segSize))
scaa.long$n.samples <- as.numeric(unlist(n.samples))
scaa.long$fracDip <- as.numeric(unlist(fracDip))
scaa.long$fracLoss <- as.numeric(unlist(fracLoss))
scaa.long$fracGain <- as.numeric(unlist(fracGain))
scaa.long$distancetoBP <- as.numeric(unlist(distancetoBP))

x <- scaa.long[which(scaa.long$fracGain == 0.5 & scaa.long$fracLoss ==0.5),] # no peaks are half gain/half loss
scaa.long$cna <- ifelse(scaa.long$fracGain>=0.5, "gain", ifelse(scaa.long$fracLoss>=0.5, "loss", "diploid"))
scaa.long$mainClonality <- pmax(scaa.long$fracDip, scaa.long$fracLoss, scaa.long$fracGain)