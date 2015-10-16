#Function to read in and create .bw files
require(SCAN.UPC)
require(rtracklayer)
require(caTools)

#Set base directly to read and write files
baseDir<-getwd()
GSM<-"GSM409307"
BWdir<-file.path(baseDir, GSM)
BWfile<-list.files(BWdir, pattern=".bw", full.names = TRUE)[1]


UPC_bigWig<-function(BWfile, binWidth=160, overlap = 40){
  
  bw<-import.bw(BWfile)
  
  #Check that the genome is binned the same throughout
  bin<-unique(bw@ranges@width)
  if(length(bin)>1) stop(paste("Different bin widths.", sep=""))
  print(paste("bigWig file is binned to ", bin, " bp.", sep=""))
  
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
  
  #Get coordinates and bin means
  startCoord = bw@ranges@start[1:(length(bw@ranges@start)-(binCount-1))]
  endCoord = bw@ranges@start[(binCount):length(bw@ranges@start)]+(bin-1)
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

start<-Sys.time()
UPCout<-UPC_bigWig(BWfile)
end<-Sys.time()
end-start
