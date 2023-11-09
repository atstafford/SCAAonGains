# FUNCTIONS -------------------------------------------------------------------------------------------------

SimpleAlignBins <- function(targetBins, incomingBins) {
  
  output <- list()
  # Bin incoming dataset 
  # per bin
  for ( b in 1:nrow(incomingBins) ) {
    chr <- incomingBins$chr[b]
    start <- incomingBins$start[b] 
    stop <- incomingBins$stop[b]
    bin <- incomingBins$bin[b]
    
    wd <- targetBins[targetBins$chr == chr,]
    wd <- wd[order(wd$start),]
    
    # no chromosome match as its sex
    if(nrow(wd)==0) {
      output[[b]] <- "autosome"
      next
    }
    
    for ( r in 1:nrow(wd) ) {
      
      ## IF START IS BEFORE ROW R 
      if ( start < wd$start[r] ) {
        
        #1.3              # stop is also before start of first row of wd
        if ( stop < wd$start[r] ) {
          output[[b]] <- "too early"
          break
        }
        #1.1              # if stop is in row r of wd, or predictor doesnt reach next bin, or this is the last row of wd...
        else if ( dplyr::between(stop, wd$start[r], wd$stop[r]) | stop < wd$start[r+1] | r == nrow(wd) ) { 
          output[[b]] <- wd$bin[r]
          break
        }
        
        #1.2              # ...or else stop is beyond row r of wd
        else if ( stop > wd$stop[r] ) {
          
          fraction <- list()
          f <- 1 
          
          # find number of extra rows required
          for ( e in 0:(nrow(wd)-r) ) {
            
            # if stop is before row r+e of wd, or stop doesnt reach next bin, or this is the last row of wd...
            if ( stop <= wd$stop[r+e] | stop < wd$start[r+e+1] | r+e == nrow(wd) ) {
              fraction[[f]] <- ( min(wd$stop[r+e], stop)-wd$start[r+e])/(stop-start) # we know it always covers wd$start
              break}
            
            else if (stop > wd$stop[r+e]) {
              fraction[[f]] <- ( min(wd$stop[r+e], stop)-wd$start[r+e])/(stop-start) # take fraction and move onto next
              f <- f + 1
              next}
          }
          
          # take max fraction
          output[[b]] <- wd$bin[r] + ((which(unlist(fraction) == max(unlist(fraction))))-1)
          break
          
        }
        
      }
      
      ## IF START IS IN ROW R
      else if ( dplyr::between(start, wd$start[r], wd$stop[r]) ) { 
        
        #2.1              # if stop is in row r of wd, or predictor doesnt reach next bin, or this is the last row of wd...
        if ( dplyr::between(stop, wd$start[r], wd$stop[r]) | stop < wd$start[r+1] | r == nrow(wd) ) { 
          output[[b]] <- wd$bin[r]
          break
        }
        
        #2.2              # ...or else stop is beyond row r of wd
        else if ( stop > wd$stop[r] ) {
          
          fraction <- list()
          f <- 1 
          
          # find number of extra rows required
          for ( e in 0:(nrow(wd)-r) ) {
            
            # if stop is before row r+e of wd, or stop doesnt reach next bin, or this is the last row of wd...
            if ( stop <= wd$stop[r+e] | stop < wd$start[r+e+1] | r+e == nrow(wd) ) {
              fraction[[f]] <- ( min(wd$stop[r+e], stop) - max(wd$start[r+e], start))/(stop-start) # we dont know it always covers wd$start
              break}
            
            else if (stop > wd$stop[r+e]) {
              fraction[[f]] <- ( min(wd$stop[r+e], stop) - max(wd$start[r+e], start))/(stop-start) # take fraction and move onto next
              f <- f + 1
              next}
          }
          
          # take max fraction
          output[[b]] <- wd$bin[r] + ((which(unlist(fraction) == max(unlist(fraction))))-1)
          break
          
        }
      }
      
      ## IF START IS IN A LATER ROW AND THERE ARE LATER ROWS
      else if ( start > wd$stop[r] & nrow(wd)>r ){
        next
      }
      
      ## IF START IS IN A LATER ROW AND THERE ARE NO LATER ROWS
      else if ( start > wd$stop[r] & nrow(wd)==r ){
        output[[b]] <- "too late"
        break
      }
    }
  }
  
  return(output)
}

PullDataInfo <- function(rawdata) {
  
  # Dataframe identifying start and stop codon, and chromosome for each bin
  start.stop <- rawdata[,c(1:3)]
  start.stop$bin <- 1:nrow(start.stop)
  start.stop <- start.stop[, c(4, 1, 2, 3)]
  
  # Vector holding sample ID
  sampleIDs <- colnames(rawdata)[-c(1:3)] 
  
  # Vector holding patient identifiers, ie the string before the '.' in sample IDs
  patientIDs <- unique(sub("\\..*", "", colnames(rawdata)))[-c(1:3)] 
  
  # Number of samples
  noSamples <- length(sampleIDs)
  
  # Number of patients
  noPatients <- length(patientIDs)
  
  # list of number of samples per patient
  sampPerPatient <- list()
  for ( i in 1:noPatients ) {
    wd <- rawdata[, which(sub("\\..*", "", colnames(rawdata))==patientIDs[i])]
    sampPerPatient[[i]] <- ncol(wd)
  }
  
  # Number of bins
  noBins <- length(start.stop$bin)
  
  # Visualisation may require a vector to identify of chromosome ends and chromosome midpoints 
  chr.ends <- cumsum(table(start.stop$chr))
  list <- list()
  l <- 1
  for ( i in 1:length(chr.ends) ) {  
    if ( i == 1 ) { 
      list[[l]] <- chr.ends[i]/2 
      l <- l+1
    }
    else { 
      list[[l]] <- chr.ends[i-1] + ((chr.ends[i]-chr.ends[i-1])/2)
      l <- l+1
    }
  }
  chr.mid <- unlist(list)
  chr.ends <- data.frame(start=c(0,chr.ends[-22]), end=(chr.ends), col=c('G','W'))
  
  #average number of bins per patient
  binsPerPatient <- list()
  for (i in 1:length(patientIDs)) {
    wd <- rawdata[,which(sub("\\..*", "", colnames(rawdata))==patientIDs[i])]
    binsPerPatient[[i]] <- round(mean(colSums(!is.na(wd))))
  }
  
  
  # Return data
  newData <- list('start.stop'=start.stop, 'sampleIDs'=sampleIDs, 'patientIDs'=patientIDs, 'noSamples'=noSamples, 'noPatients'=noPatients,
                  'sampPerPatient'=sampPerPatient, 'noBins'=noBins,'chr.mid'=chr.mid,  'chr.end'=chr.ends, 'binsPerPatient'=binsPerPatient)
  return(newData)
}
PullDataClonality <- function(rawdata, dataInfo) {
  
  # A dataframe with a bin per row and a patient per column, with values indicating clonality. 
  # 0=notCNA, 1=subclonalCNA, 2=clonalCNA
  clonal.data <- rawdata[,-c(1:3)]
  l <- 1
  clonal <- list()
  for ( k in 1:nrow(clonal.data) ) {
    for ( i in 1:length(dataInfo$patientIDs) ) {
      wd <- clonal.data[k, which(sub("\\..*", "", colnames(clonal.data))==dataInfo$patientIDs[i])]
      wd <- wd[, is.na(wd)!=TRUE]
      
      if ( 1 %in% wd | 3 %in% wd ) { #if one of the samples has a mutation, proceed
        if ( length(unique(t(wd)))==1 ) { #if all the same then clonal
          clonal[[l]] <- 2
          l <- l + 1
        }
        else { #different = subclonal
          clonal[[l]] <- 1
          l <- l + 1
        }
      }
      else if ( length(wd)==0 ) { #all MR for that bit are NA
        clonal[[l]] <- NA
        l <- l + 1
      }
      else { #neither sample has a mutation
        clonal[[l]] <- 0 
        l <- l + 1
      }
    }
  }
  
  clonal.data <- data.frame(t(matrix(unlist(clonal), ncol=dataInfo$noBins)))
  colnames(clonal.data) <- dataInfo$patientIDs
  clonal.data[] <- lapply(clonal.data, factor, levels=unique(unlist(clonal.data)))
  
  
  # Dataframe detailing the counts of gains/losses and whether they are subclonal or clonal
  CNA.clo.counts <- data.frame(bin = dataInfo$start.stop$bin, chr1 = NA, chr2 = NA,
                               clonal.aneu = NA, subclonal.aneu = NA, gain = NA, loss = NA,
                               clonal.gain = NA, clonal.loss = NA, clonal.noCNA = NA, 
                               subclonal.gain = NA, subclonal.loss = NA)
  
  CNA.clo.counts$chr1 <- as.numeric(dataInfo$start.stop[match(CNA.clo.counts$bin, dataInfo$start.stop$bin), 2]) 
  CNA.clo.counts$chr2 <- as.numeric(dataInfo$start.stop[match(CNA.clo.counts$bin, dataInfo$start.stop$bin), 2]) 
  CNA.clo.counts$chr2 <- factor(CNA.clo.counts$chr2,levels = rev(seq(1:22))) 
  
  data <- rawdata[,-c(1:3)]
  for ( k in 1:nrow(data) ) { # for a bin
    clonal.all <- clonal.gain <- clonal.loss <- clonal.noCNA <- subclonal.all <- subclonal.gain <- subclonal.loss <- subclonal.noCNA <- 0
    
    for ( i in 1:dataInfo$noPatients ) { # for a patient
      
      # If its NA
      if ( is.na(clonal.data[k,i]) == TRUE ) {
        next
      }
      
      # If its clonal
      if ( clonal.data[k,i] == 2 ) {
        wd <- data[k, which(sub("\\..*", "", colnames(data))==dataInfo$patientIDs[i])]
        wd <- wd[, is.na(wd)!=TRUE]
        
        if ( 3 %in% wd ) { # if both are gains
          clonal.gain <- clonal.gain + 1
        }
        else if ( 1 %in% wd ) { # if both are losses
          clonal.loss <- clonal.loss + 1
        }
      }
      
      # if its subclonal
      else if ( clonal.data[k,i] == 1 ) {
        wd <- data[k, which(sub("\\..*", "", colnames(data))==dataInfo$patientIDs[i])]
        wd <- wd[, is.na(wd)!=TRUE]
        
        if ( 3 %in% wd ) { # if one is a gain
          subclonal.gain <- subclonal.gain + 1
        }
        if ( 1 %in% wd ) { # if one is a loss
          subclonal.loss <- subclonal.loss + 1
        }
      }
      
      # if its no CNA
      else if ( clonal.data[k,i] == 0 ) {
        clonal.noCNA <- clonal.noCNA + 1
      }
    }
    
    CNA.clo.counts$clonal.gain[k] <- clonal.gain
    CNA.clo.counts$clonal.loss[k] <- clonal.loss
    CNA.clo.counts$clonal.noCNA[k] <- clonal.noCNA
    CNA.clo.counts$subclonal.gain[k] <- subclonal.gain
    CNA.clo.counts$subclonal.loss[k] <- subclonal.loss
  }
  
  CNA.clo.counts$clonal.aneu <- CNA.clo.counts$clonal.gain + CNA.clo.counts$clonal.loss
  CNA.clo.counts$subclonal.aneu <- CNA.clo.counts$subclonal.gain + CNA.clo.counts$subclonal.loss
  CNA.clo.counts$gain <- CNA.clo.counts$clonal.gain + CNA.clo.counts$subclonal.gain
  CNA.clo.counts$loss <- CNA.clo.counts$clonal.loss + CNA.clo.counts$subclonal.loss
  
  # dataframe showing: bin | countGain/NoPatient | countLoss/NoPatient | 
  CloFreq <- cbind(bin=CNA.clo.counts$bin, gain=CNA.clo.counts$gain/dataInfo$noPatients, loss=CNA.clo.counts$loss/dataInfo$noPatients)
  
  # dataframe showing what percent of gain and loss are subclonal. On a patient basis: bin | gain | loss
  pcSubclonal <- data.frame(bin=1:dataInfo$noBins, gain=CNA.clo.counts$subclonal.gain / CNA.clo.counts$gain, loss=CNA.clo.counts$subclonal.loss / CNA.clo.counts$loss)
  
  # Count of noCNA, subclonalCNA, and clonalCNA by patient
  patientClo <- as.data.frame(t(sapply(clonal.data, table))) 
  if( nrow(patientClo)==3) {
    patientClo <- patientClo[, c('0','1','2')]
    colnames(patientClo) <- c('noCNA','subclonal','clonal')
    patientClo$CNA <- patientClo$subclonal + patientClo$clonal
  } else if ( "2" %!in% colnames(patientClo) ) { # no clonal
    patientClo <- patientClo[, c('0','1')]
    colnames(patientClo) <- c('noCNA','subclonal')
    patientClo$CNA <- patientClo$subclonal
  }
  
  
  patientClo$patient <- rownames(patientClo)
  
  # Return data
  newData <- list('clonal.data'=clonal.data, 'CNA.clo.counts'=CNA.clo.counts, 
                  'CloFreq'=CloFreq, 'pcSubclonal'=pcSubclonal, 'patientClo'=patientClo)
  return(newData)
}
PullDataDiversity <- function(rawdata, dataInfo) {
  
  # Cannot use averaged raw CN across the MR samples in case of variable number of samples per patient.
  # Will therefore can use a PIC.frac for each patient. 
  # This is a continuous, not binary, measure of clonality
  # Define the maximum pic score depending on the number of samples (up to 13)
  max.pics <- list()
  for ( d in 1:50 ) {
    if ( d %% 3 == 0 ) {
      n1 <- n2 <- n3 <- d/3
    }
    else if ( d %% 3 == 1 ) {
      n1 <- ((d-1)/3) + 1
      n2 <- ((d-1)/3)
      n3 <- ((d-1)/3)
    }
    else if ( d %% 3 == 2 ) {
      n1 <- ((d-2)/3) + 1
      n2 <- ((d-2)/3) + 1
      n3 <- ((d-2)/3)
    }
    max.pics[[d]] <- pic.score <- 1 - ( (n1/d)^2 + (n2/d)^2 + (n3/d)^2 )
  }
  
  # A dataframe with a bin per row and a patient per column, with values indicating pic score
  # And a dataframe showing pic.frac per patient
  pic <- list()
  pic.frac <- list()
  for ( i in 1:length(dataInfo$patientIDs)) {
    # Set working data as the cols holding the samples
    wd <- rawdata[, which(sub("\\..*", "", colnames(rawdata))==dataInfo$patientIDs[i])]
    
    # Record the number of samples
    upto <- ncol(wd)
    
    # Use PIC function on wd
    pic[[i]] <- PIC(wd, upto, c(1:upto))
    pic[[i]] <- na.omit(pic[[i]])
    
    # Define the max possible diversity given the number of sample, as maxPIC*number of bins
    # max.ITH <- max.pics[[upto]] * length(pic[[i]])
    max.ITH <- max.pics[[upto]] * dataInfo$binsPerPatient[[i]]
    
    # Store as a dataframe
    pic.frac[[i]] <- data.frame(pic.frac = sum(pic[[i]], na.rm = TRUE)/max.ITH)
  }
  
  pic.data <- as.data.frame(do.call('cbind',pic))
  colnames(pic.data) <- dataInfo$patientIDs
  pic.frac <- do.call('rbind',pic.frac)
  
  # Calulate average PIC per bin as a measure of average bin hetero across patients 
  ave.pic <- dataInfo$start.stop[,c(1:2)]
  ave.pic$avePic <- rowSums(pic.data)/(ncol(pic.data))
  
  # A dataframe of the propotion of genome gained/lost per sample, alongside ith
  pga <- data.frame(t(rawdata[,-c(1:3)]), check.names = FALSE)
  pga <- data.frame(prop.gain=apply(pga,1,function(x) sum(x == 3, na.rm = TRUE)/ncol(pga)),
                    prop.loss=apply(pga,1,function(x) sum(x == 1, na.rm = TRUE)/ncol(pga)))
  pga$prop.aneu <- pga$prop.gain + pga$prop.loss
  pga <- cbind(as.data.frame(lapply(pic.frac, rep, dataInfo$sampPerPatient)),
               pga)
  
  # Return data
  newData <- list('pic.data'=pic.data, 'pic.frac'=pic.frac, 'ave.pic'=ave.pic, 'pga'=pga)
  return(newData)
}
PIC <- function(data, sample.no, index) { 
  # based on table of unique patient nos
  # returns PIC per row (bin), calc across cols as given by index
  # PIC formula: 1- ((CN1/n)^2 + (CN2/n)^2 + (CN3/n)^2), where CN1 is no. of counts of copy number1
  PIC <- 1 - ((rowSums(data[,index]==1, na.rm = TRUE)/sample.no)^2 + 
                (rowSums(data[,index]==2, na.rm = TRUE)/sample.no)^2 + 
                (rowSums(data[,index]==3, na.rm = TRUE)/sample.no)^2)
  return(PIC)
}
'%!in%' <- function(x,y)!('%in%'(x,y))

grouper <- function(df, n) {
  
  # create a random number for each row
  random <- sample(1:nrow(df), replace = FALSE, nrow(df))
  
  # divide the random number by the group size
  df$group_number <- ceiling(random / (nrow(df) / n))
  
  return(df)  
}

getChrLength <- function(genome = "BSgenome.Hsapiens.UCSC.hg38"){
  g <- getBSgenome(genome, masked=FALSE)
  data.frame(chrom=1:24, length=seqlengths(g)[1:24])
}
.chrAsNum <- function(tbl){
  tbl$chrom <- gsub("chr", "", tbl$chrom)
  tbl$chrom[tbl$chrom=="X"] <- 23
  tbl$chrom[tbl$chrom=="Y"] <- 24
  tbl$chrom <- as.numeric(tbl$chrom)
  tbl[order(tbl$chrom),]
}
getCentromeres <- function( genome="hg38" ){
  mySession <- try(browserSession("UCSC"), silent=TRUE)
  # In case of failure, try another mirror
  if(inherits(mySession, "try-error"))
    mySession <- browserSession("UCSC",
                                url="http://genome-euro.ucsc.edu/cgi-bin/")
  genome(mySession) <- genome
  obj <- ucscTableQuery(mySession, table="gap")
  tbl <- getTable(obj)
  tbl <- tbl[tbl$type=="centromere", c("chrom", "chromStart", "chromEnd")]
  colnames(tbl)[2:3] <- c("centromerStart", "centromerEnd")
  .chrAsNum(tbl)
}
makeHg38 <- function(){
  tbl <- merge(getChrLength(), getCentromeres(), by="chrom")
  cumlen <- c(0, cumsum(as.numeric(tbl$length))[-nrow(tbl)])
  cbind.data.frame(tbl, cumlen=cumlen)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

four_quadrant <- function(outsidebox_size, insidebox_size, 
                          outsidebox_label, insidebox_label, 
                          col_quad1=c(alpha("#FF3366",0.4),alpha("#FF9900",0.6),alpha("#9900CC",0.4),alpha("#66CCCC",0.7)), 
                          col_quad2=c(alpha("#FF3366",0.7),alpha("#FF9900",0.9),alpha("#9900CC",0.7),alpha("#66CCCC",1)), 
                          col_text1="black", col_text2="white",
                          p.value) {    
  nx <- length(outsidebox_size)
  sqx <- sqrt(outsidebox_size) 
  outside_df <- data.frame(x=c(sqx[1],-sqx[2],-sqx[3],sqx[4])/2, 
                           y=c(sqx[1],sqx[2],-sqx[3],-sqx[4])/2, 
                           size=sqx, label=paste("n=",outsidebox_label, sep = ""))
  outside_df$label.yloc <- outside_df$y+(outside_df$size/2.3)
  outside_df$label.yloc[which(outside_df$y<0)] <- outside_df$y[which(outside_df$y<0)]-(outside_df$size[which(outside_df$y<0)]/2.3)
  
  if( !is.null(insidebox_size) ) {
    ny <- length(insidebox_size)
    sqy <- sqrt(insidebox_size) 
    inside_df <- data.frame(x=c(sqy[1],-sqy[2],-sqy[3],sqy[4])/2, 
                            y=c(sqy[1],sqy[2],-sqy[3],-sqy[4])/2, 
                            size=sqy, label=paste(round(insidebox_label, digits = 2),"X",sep = ""))
  }
  
  mm <- max(outside_df$size)*1.1
  
  p <- ggplot() +
    geom_tile(data=outside_df, aes(x=x, y=y, width=size, height=size, 
                                   group=factor(size)), fill=col_quad1) +
    geom_text(data=outside_df, aes(x=x, y=label.yloc, label=label), col=col_text1, size=5, fontface="italic")
  
  if( !is.null(insidebox_size) ) {
    p <- p +
      geom_tile(data=inside_df, aes(x=x, y=y, width=size, height=size, 
                                    group=factor(size)), fill=col_quad2) +
      geom_text(data=inside_df, aes(x=x, y=y, label=label), col=col_text2, size=5) 
  }
  
  p <- p + geom_hline(aes(yintercept=0), size=0.8, linetype="dashed") +
    geom_vline(aes(xintercept=0), size=0.8, linetype="dashed") +
    xlim(c(-mm,mm)) + ylim(c(-mm,mm)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.border = element_rect(fill=NA),
          plot.margin=margin(t=0,r=0,b=0,l=0,"cm")) 
  p <- p + annotate("text", x=c(-mm,-mm,mm,mm), y=c(mm,-mm,mm,-mm), label=c("Q1","Q2","Q3","Q4"), fontface = "bold")
  
  if( !is.na(p.value) ) {
    p <- p + annotate("text", x=c(-mm+mm*1.1), y=c(-mm), label=paste("p=",round(p.value, digits = 3), sep = ""), fontface = "italic")
  }
  
  p <- p
}

quadEnrichmentChi <- function(enrichedPeakList, q1.peaks, q2.peaks, q3.peaks, q4.peaks) {
  
  totalTarget <- length(enrichedPeakList)
  
  q1yes <- length(intersect(q1.peaks, enrichedPeakList))
  q1no <- length(q1.peaks) - q1yes
  q2yes <- length(intersect(q2.peaks, enrichedPeakList))
  q2no <- length(q2.peaks) - q2yes
  q3yes <- length(intersect(q3.peaks, enrichedPeakList))
  q3no <- length(q3.peaks) - q3yes
  q4yes <- length(intersect(q4.peaks, enrichedPeakList))
  q4no <- length(q4.peaks) - q4yes
  
  df <- data.frame(yes = c(q1yes,q2yes,q3yes,q4yes), no = c(q1no,q2no,q3no,q4no))
  chisq <- chisq.test(df)
  
  corplot <- corrplot(chisq$residuals, is.corr = F, method="number")
  corrpv <- chisq$p.value
  
  obs <- data.frame(chisq$observed)
  exp <- data.frame(chisq$expected)
  
  enrich.q1 <- (obs$yes[1]/totalTarget) / (exp$yes[1]/totalTarget)
  enrich.q2 <- (obs$yes[2]/totalTarget) / (exp$yes[2]/totalTarget)
  enrich.q3 <- (obs$yes[3]/totalTarget) / (exp$yes[3]/totalTarget)
  enrich.q4 <- (obs$yes[4]/totalTarget) / (exp$yes[4]/totalTarget)
  
  outsidebox_size <- c(length(q3.peaks), length(q1.peaks), length(q2.peaks), length(q4.peaks))
  insidebox_size <- c(q3yes, q1yes, q2yes, q4yes)
  outsidebox_label <- c(length(q3.peaks), length(q1.peaks), length(q2.peaks), length(q4.peaks))
  insidebox_label <- c(enrich.q3, enrich.q1, enrich.q2, enrich.q4)
  quad.plot <- four_quadrant(outsidebox_size = outsidebox_size, insidebox_size = insidebox_size, 
                             outsidebox_label = outsidebox_label, insidebox_label = insidebox_label, p.value = corrpv)
  
  list <- list("quad.plot"=quad.plot,"enrich.q1"=enrich.q1, "enrich.q2"=enrich.q2, "enrich.q3"=enrich.q3, "enrich.q4"=enrich.q4)
  return(list)
}

# Function to align bin with new bins and use weight average copy number
newAlignBins <- function(bins, cn.list) {
  # bins needs to be a dataframe holding: bin | chr | start | stop, for the bins to align to
  # cn.list is the output from ploidyRecentre with skipcol=3
  
  # Create dataframe for each patient and put into list
  cnBinned.list <- rep(list(bins), length(cn.list))
  
  # Convert to numeric
  cnBinned.list <- lapply(cnBinned.list, function (x) {
    x[] <- apply(x,2,as.numeric)
    x
  })
  
  # Add empty columns to hold output
  for ( i in 1:length(cnBinned.list) ) {
    sampleNo <- ncol(cn.list[[i]])-3 #
    cnBinned.list[[i]] <- do.call("cbind", list(cnBinned.list[[i]], rep(list(NA), sampleNo)))
    colnames(cnBinned.list[[i]])[-c(1:4)] <- colnames(cn.list[[i]])[-c(1:3)]
  }
  
  # Bin incoming dataset (cn.list)
  for ( p in 1:length(cnBinned.list) ) {
    print(paste(p,"/",length(cnBinned.list)))
    # per bin
    for ( b in 1:nrow(cnBinned.list[[p]]) ) {
      chr <- cnBinned.list[[p]]$chr[b]
      start <- cnBinned.list[[p]]$start[b] 
      stop <- cnBinned.list[[p]]$stop[b]
      bin <- cnBinned.list[[p]]$bin[b]
      
      wd <- cn.list[[p]][cn.list[[p]]$chr == chr,]
      wd <- wd[order(wd$start),]
      
      for ( r in 1:nrow(wd) ) {
        
        ## IF START IS BEFORE ROW R 
        if ( start < wd$start[r] ) {
          
          #1.3              # stop is also before start of first row of wd
          if ( stop < wd$start[r] ) {
            cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- "too early"
            break
          }
          #1.1              # if stop is in row r of wd, or predictor doesnt reach next bin, or this is the last row of wd...
          else if ( dplyr::between(stop, wd$start[r], wd$stop[r]) | stop < wd$start[r+1] | r == nrow(wd) ) { 
            cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- wd[r,-c(1:3)]
            break
          }
          
          #1.2              # ...or else stop is beyond row r of wd
          else if ( stop > wd$stop[r] ) {
            
            fraction <- list()
            f <- 1 
            
            # find number of extra rows required
            for ( e in 0:(nrow(wd)-r) ) {
              
              # if stop is before row r+e of wd, or stop doesnt reach next bin, or this is the last row of wd...
              if ( stop <= wd$stop[r+e] | stop < wd$start[r+e+1] | r+e == nrow(wd) ) {
                fraction[[f]] <- ( min(wd$stop[r+e], stop)-wd$start[r+e])/(stop-start) # we know it always covers wd$start
                break}
              
              else if (stop > wd$stop[r+e]) {
                fraction[[f]] <- ( min(wd$stop[r+e], stop)-wd$start[r+e])/(stop-start) # take fraction and move onto next
                f <- f + 1
                next}
            }
            
            # take weighted mean
            cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- colSums( (data.frame(wd[c(r:(r+e)),-c(1:3)]) *unlist(fraction)), na.rm=T) / sum(unlist(fraction))
            break
            
          }
          
        }
        
        ## IF START IS IN ROW R
        else if ( dplyr::between(start, wd$start[r], wd$stop[r]) ) { 
          
          #2.1              # if stop is in row r of wd, or predictor doesnt reach next bin, or this is the last row of wd...
          if ( dplyr::between(stop, wd$start[r], wd$stop[r]) | stop < wd$start[r+1] | r == nrow(wd) ) { 
            cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- wd[r,-c(1:3)]
            break
          }
          
          #2.2              # ...or else stop is beyond row r of wd
          else if ( stop > wd$stop[r] ) {
            
            fraction <- list()
            f <- 1 
            
            # find number of extra rows required
            for ( e in 0:(nrow(wd)-r) ) {
              
              # if stop is before row r+e of wd, or stop doesnt reach next bin, or this is the last row of wd...
              if ( stop <= wd$stop[r+e] | stop < wd$start[r+e+1] | r+e == nrow(wd) ) {
                fraction[[f]] <- ( min(wd$stop[r+e], stop) - max(wd$start[r+e], start))/(stop-start) # we dont know it always covers wd$start
                break}
              
              else if (stop > wd$stop[r+e]) {
                fraction[[f]] <- ( min(wd$stop[r+e], stop) - max(wd$start[r+e], start))/(stop-start) # take fraction and move onto next
                f <- f + 1
                next}
            }
            
            # take weighted mean
            cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- colSums( (data.frame(wd[c(r:(r+e)),-c(1:3)]) *unlist(fraction)), na.rm=T) / sum(unlist(fraction))
            break
            
          }
        }
        
        ## IF START IS IN A LATER ROW
        else if ( start > wd$stop[r] ){
          next
        }
        
      }
    }
  }
  return(cnBinned.list)
}

theme_custom <- function(){
  theme_classic() %+replace%  
    theme(
      plot.margin=margin(t=0.2,r=0.5,b=0.5,l=0.5,"cm"),
      panel.background = element_blank(),
      plot.title = element_text(size=28, colour='black',face='bold', hjust=0),
      axis.line = element_line(size = 0.5, colour = "black"),
      axis.title = element_text(size=28, colour='black'),
      axis.text = element_text(size=28, colour='black'),
      axis.ticks.length=unit(0.2, "cm")
    )
}
theme_legend <- function(){
  theme_classic() %+replace%  
    theme(
      plot.margin = unit(c(t=-100,r=0,b=-100,l=0), "cm"),
      legend.position = "top",legend.direction="horizontal",
      legend.title = element_text(size=28, colour='black',face='bold'),
      legend.margin = margin(grid::unit(c(t=-100,r=0,b=-100,l=0),"cm")),
      legend.text = element_text(size=28, colour='black'),
      legend.key.height = grid::unit(0.8,"cm"),
      legend.key.width = grid::unit(1.4,"cm")
    )
}
# Function to plot frequency and clonality of gains and losses
cloFreqPlot <- function(clonalityData, dataInfo, annotation, title=NULL, ylab="Fraction of patients with aberration", xlab="Chromosome", colourChoice=c("#1B7837","#5AAE61","#A6DBA0","#D9F0D3","#E7D4E8","#C2A5CF","#9970AB","#762A83")) {
  
  # pull CloFreq
  CloFreq <- data.frame(clonalityData$CloFreq)
  
  # Make losses negative to create a mirror plot
  CloFreq$loss <- CloFreq$loss*-1 
  
  # Gather into long form
  CloFreq <- gather(CloFreq,CNA,freq,-bin) 
  
  # Look up the percent of the gain/loss that is subclonal
  CloFreq$pcSubclonal <- NA
  
  for ( i in 1: nrow(CloFreq) ) {
    bin <- CloFreq$bin[i]
    CNA <- CloFreq$CNA[i]
    
    if ( CNA == 'gain' ) {
      CloFreq$pcSubclonal[i] <- clonalityData$pcSubclonal$gain[bin] *100
    }
    if ( CNA == 'loss' ) {
      CloFreq$pcSubclonal[i] <- clonalityData$pcSubclonal$loss[bin] *100
    }
  }
  
  # Position label along y axis
  for ( i in 1:nrow(annotation) ) {
    bin <- annotation$x[i]
    data <- CloFreq[which(CloFreq$bin==bin),]
    annotation$y[i] <- data$freq[which.max(abs(data$freq))]
  }
  
  # Create plots
  fig <- ggplot() + 
    geom_rect(data = dataInfo$chr.end[which(dataInfo$chr.end$col=='W'),], 
              aes(NULL,NULL,xmin=start, xmax=end),
              fill = alpha("#CCCCCC", 0.1),
              ymin = -1,
              ymax = 1) +
    
    geom_bar(data=CloFreq, aes(fill=pcSubclonal, y=freq, x=bin),
             position="stack", stat="identity", width = 1) +
    scale_fill_gradientn(colours = colourChoice) +
    
    ggtitle(title) +
    scale_x_continuous(expand = c(0,0), name=xlab, breaks=dataInfo$chr.mid, labels = c(1:18,'\n19','20','\n21','22')) +
    scale_y_continuous(expand = c(0,0), name=ylab, limits = c(-1,1), breaks = c(seq(-1,1,0.2)), 
                       labels = c(1,0.8,0.6,0.4,0.2,0,0.2,0.4,0.6,0.8,1)) +
    
    geom_hline(yintercept = 0, size=0.3, color="black") +
    
    geom_point(data = annotation, aes(x = x, y = y),
               shape = 18, size = 0) +
    geom_text_repel(data = subset(annotation, y >= 0),
                    aes(x = x, y = y, label = label),
                    size = 7,
                    nudge_y = 0.2,direction = 'x',
                    arrow = arrow(length = unit(0.015, "npc"))) +
    geom_text_repel(data = subset(annotation, y < 0),
                    aes(x = x, y = y, label = label),
                    size = 7,
                    nudge_y = -0.2,direction = 'x',
                    arrow = arrow(length = unit(0.015, "npc"))) +
    annotate('text', y=c(0.9,-0.9), x=c(50,50), 
             label=c('italic(Gains)','italic(Losses)'), size=10, parse=TRUE, hjust = 0) +
    theme_custom() +
    theme(legend.position = "none", plot.title = element_text(hjust=0))
  
  # Create legend
  galo.legend <- as_ggplot(cowplot::get_legend(fig + 
                                                 guides(fill = guide_colorbar(title="% subclonal", label.position = "bottom",
                                                                              title.position = "left", title.vjust = 0.9)) +
                                                 theme_legend() ))
  
  # remove negative for loss to allow for calculation of correlation
  CloFreq$freq[which(CloFreq$freq<0)] <- CloFreq$freq[which(CloFreq$freq<0)]*-1
  return(list(CloFreqDF = CloFreq, plot = fig, legend=galo.legend))
}

# Function to plot a split violin plot
library(tidyverse)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# annotate genes
gene_anno <- function (top.bins, results) {
  #top.bins data must have bin=col1, chr=col2, start=col3, stop=col4
  #results (gene annotation) must have hgnc_symbol=col1, entrez=col2, chr=col3, start=col4, stop=col5, geneID=col6
  k <- 1
  i <- 1
  l <- 1
  g <- c <- sta <- sto <- id <- b <- NULL
  g.l <- c.l <- sta.l <- sto.l <- id.l <- b.l <- list() #to match genes to bins
  
  for ( k in 1:nrow(top.bins) ) {
    bin <- top.bins[k,1]
    chr <- top.bins[k,2]
    start <- top.bins[k,3]
    stop <- top.bins[k,4]
    
    wd <- results[which(results[,3] == chr),]
    
    for ( i in 1:nrow(wd) ) {
      if ( (chr == wd[i,3]) && ( between(wd[i,4], start, stop) || between(wd[i,5], start, stop)) ) {
        g <- append(g,wd[i,2])
        c <- append(c,wd[i,3])
        sta <- append(sta,wd[i,4])
        sto <- append(sto,wd[i,5])
        id <- append(id,wd[i,6])
        b <- append(b,bin)
      }
    }
    g.l[[l]] <- g
    c.l[[l]] <- c
    sta.l[[l]] <- sta
    sto.l[[l]] <- sto
    id.l[[l]] <- id
    b.l[[l]] <- b
    l <- l+1
    g <- c <- sta <- sto <- id <- b <- NULL
  }
  top.genes <- data.frame(bin=unlist(b.l), entrezgene=unlist(g.l), chr=unlist(c.l), gene.start=unlist(sta.l), gene.stop=unlist(sto.l), geneID=unlist(id.l))
  #top.genes[top.genes==""] <- NA
  top.genes <- na.omit(top.genes) #remove rows with blanks
  #colnames(top.genes) <- c("bin",colnames(wd)[2])
  #top.genes <- merge(top.genes, results, by = colnames(wd)[2]) #add in data from gene annotation data
  #colnames(top.genes)[4] <- "chr"
  return(top.genes)
}

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}