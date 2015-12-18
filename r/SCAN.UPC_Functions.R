require(SCAN.UPC)
require(rtracklayer)
require(caTools)

#Function to read in and create .bw files
geo2bw<-function(GSM, baseDir){
  #Create temporary directory
  tempDir<-file.path(baseDir, "GEOtemp")
  dir.create(tempDir)
  
  #Download all data files from GEO for this experiment
  getGEOSuppFiles(GSM, makeDirectory=FALSE, baseDir=tempDir)
  
  #Get list of wig files and zips
  wigFiles<-list.files(tempDir, pattern = ".wig")
  
  #Create new directory to put GSM files in
  newDir<-file.path(baseDir, GSM)
  dir.create(newDir)
  
  #Copy the .wig and zips to new directory
  for(i in wigFiles){
    file.copy(file.path(tempDir, i), newDir)}
  
  #Delete the temporary directory and files
  unlink(tempDir, recursive=TRUE, force=TRUE)
  
  #Create bigWig files from wig files
  wigFiles<-list.files(newDir, pattern=".wig")
  
  #Get .soft file on this data
  getGEOfile(GSM, file.path(newDir))
  
  for(i in wigFiles){
    
    #Get Assembly information from .soft file
    build<-"empty"
    soft<-readLines(file.path(newDir, paste(GSM, ".soft", sep="")))
    iSub<-sub(".gz", "", i)
    end<-FALSE
    for(j in 1:length(soft)){
      line<-soft[j]
      if(grepl(iSub, line) & end==FALSE){
        flag<-TRUE
        q<-j+1
        while(flag==TRUE & q<length(soft)){
          q<-q+1
          lineSub<-soft[q]
          if(grepl("GENOME_ASSEMBLY", lineSub)){
            build<-lineSub
            flag<-FALSE
            end<-TRUE}
          if(grepl("[*******]", lineSub)){
            flag<-FALSE}}}}
    
    #Set Seq to the correct  build.  Otherwise you can't write the bigWig file. 
    Seq<-NULL
    if(grepl("hg19", build)){
      build<-"hg19"
      Seq<-Seqinfo(genome=build)}
    if(grepl("hg18", build)){
      build<-"hg18"
      Seq<-Seqinfo(genome=build)}
    if(grepl("hg38", build)){
      build<-"hg38"
      Seq<-Seqinfo(genome=build)}
    
    if(!is.null(Seq)){
      wigToBigWig(file.path(newDir, i), seqinfo = Seq, dest = paste(file_path_sans_ext(file.path(newDir, i), TRUE),"_", build, ".bw", sep=""))} else {
        print(paste("Could not create bigWig file from, ", i," because genome build is unknown."))}
  }
}


#Function to run UPC on BigWig files
UPC_bigWig<-function(BWfile, binWidth=150, overlap = 25){
  
  bw<-import.bw(BWfile)
  
  
  #Get coordinates
  startCoord <- bw@ranges@start
  
  #Get consensus bin sizes
  bin<-startCoord[2:length(startCoord)]-startCoord[1:(length(startCoord)-1)]
  bin<-table(sample(bin, 1000))
  bin<-as.numeric(names(which.max(bin)))
  print(paste("bigWig file binned to ", bin, " bp.", sep=""))
  
  #Find number of adjacent bins to combine
  binCount<-ceiling(binWidth/bin)
  ActualBinWidth<-binCount*bin
  print(paste("Creating signal object. BinWidth=", ActualBinWidth, " bp.", sep=""))
  
  #Find bin overlap
  overCount<-ceiling(overlap/bin)
  ActualOverlap<-overCount*bin
  
  if(overCount>=binCount) stop(paste("Overlap size too large. For BinWidth = ", binWidth, " bp the maximum overlap is ", binWidth-bin, " bp.", sep=""))
  print(paste("Adjacent regions will overlap by ", ActualOverlap, "bp on either side.", sep=""))
  
  #Get chromosome list
  Chrs<-as.character(bw@seqnames@values)
  Length<-bw@seqnames@lengths
  
  #Generate vector of chromosome assignments
  chrVec<-NULL
  for(i in 1:length(Chrs)){
    chrVec<-c(chrVec, rep(Chrs[i], Length[i]))}
  
  endCoord <- startCoord[(binCount):length(startCoord)]+(bin-1)
  startCoord<-startCoord[1:(length(startCoord)-(binCount-1))]
  binMeans<-runmean(bw@elementMetadata@listData$score, k=binCount, alg="fast", endrule = "trim")
  
  #Remove GIANT bigWig object
  rm(bw)
  
  #Create data frame of coordinates and bin means
  binStats<-data.frame(chr = chrVec[1:length(startCoord)], start = startCoord, end = endCoord, binMeans = binMeans)
  
  #Reove GIANT vectors
  rm(startCoord, endCoord, binMeans, chrVec)
  
  #Flag rows that don't have width = ActualBinWidth
  BinCheck <- binStats$end - binStats$start + 1 == ActualBinWidth
  
  #Create vector of overlap keepers
  overKeep <-c(TRUE, rep(FALSE, binCount - overCount - 1))
  overKeep <-rep_len(overKeep, length.out = nrow(binStats))
  
  
  #Subset for overlap keepers and correct bin width
  binStats<-binStats[BinCheck & overKeep,]
  
  UPCout <- UPC_Generic(binStats$binMeans)
  binStats$UPCout<-UPCout
  
  
  return(binStats)
}

#Alternate function to run UPC on BigWig files.  This will bin all data file from the same genome build the same.
UPC_bigWig_NO<-function(BWfile, 
                        binWidth = 150, 
                        build = "hg19", 
                        chrs = NULL){
  
  #Tile the genome for specified width
  if(is.character(build)) Seq = Seqinfo(genome=build) else Seq = build
  
  bins <- tileGenome(Seq, tilewidth=binWidth,
                     cut.last.tile.in.chrom=TRUE)
  bins <- keepSeqlevels(bins, chrs)
  
  #Read in BigWig file as RleList
  bw <- import.bw(BWfile, as = "RleList")
  
  #Order the RleList to match the tiled genome (bins)
  bw <- bw[names(bw) %in% seqlevels(bins)]
  bw <- bw[order(match(names(bw), seqlevels(bins)))]
  
  #Find average signal in each bin
  bwBinned <- binnedAverage(bins, bw, "signal")
  if(!is.null(chrs)) bwBinned <- bwBinned[bwBinned@seqnames%in%chrs,]
  
  #Create data from from binned averages
  bwFram <- as.data.frame(bwBinned)
  bwFram$UPCout <- 0
  
  bwFram$UPCout[bwFram$signal!=0] <- UPC_Generic(bwFram$signal[bwFram$signal!=0])
  
  return(bwFram)
}

#Alternate function to run UPC on BigWig files.  This requires an annotation file in GRanges format.
UPC_bigWig_annot<-function(BWfile,
                           annot, 
                           chrs = NULL,
                           max = FALSE){
  
  #Read in BigWig file as RleList
  bw <- import.bw(BWfile, as = "RleList")
  
  #Order the RleList to match the tiled genome (bins)
  bw <- bw[names(bw) %in% seqlevels(annot)]
  bw <- bw[order(match(names(bw), seqlevels(annot)))]
  
  #function (bins, numvar, varname) 
  binnedMax<-function(bins, numvar, varname) {
    if (!is(bins, "GRanges")) 
      stop("'x' must be a GRanges object")
    if (!is(numvar, "RleList")) 
      stop("'numvar' must be an RleList object")
    if (!identical(seqlevels(bins), names(numvar))) 
      stop("'seqlevels(bin)' and 'names(numvar)' must be identical")
    viewMeans2 <- function(v) {
      means <- viewMaxs(v)
      w0 <- width(v)
      w1 <- width(trim(v))
      means <- means * w1/w0
      means[w0 != 0L & w1 == 0L] <- 0
      means
    }
    bins_per_chrom <- split(ranges(bins), seqnames(bins))
    means_list <- lapply(names(numvar), function(seqname) {
      v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
      viewMeans2(v)
    })
    new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
    mcols(bins)[[varname]] <- new_mcol
    bins
    }
  
  #Find average signal in each bin
  if(max==FALSE){
    bwBinned <- binnedAverage(annot, bw, "signal")
  } else {bwBinned <- binnedMax(annot, bw, "signal")}
  if(!is.null(chrs)) bwBinned <- bwBinned[bwBinned@seqnames%in%chrs,]
  
  #Create data from from binned averages
  bwFram <- as.data.frame(bwBinned)
  
  bwFram$UPCout<-UPC_Generic(bwFram$signal, lengths = bwFram$width)
  colnames(bwFram)[1]<-"chr"
  
  return(bwFram)
}