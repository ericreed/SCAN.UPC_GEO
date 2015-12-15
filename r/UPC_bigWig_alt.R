#Function to read in and create .bw files
require(SCAN.UPC)
require(rtracklayer)
require(caTools)

#Set base directly to read and write files
baseDir<-getwd()
GSM<-"GSM409307"
BWdir<-file.path(baseDir, GSM)
BWfile<-list.files(BWdir, pattern=".bw", full.names = TRUE)[1]



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

start<-Sys.time()
UPCout<-UPC_bigWig_NO(BWfile, chrs = paste("chr", c(1:22, "X", "Y"), sep=""))
end<-Sys.time()
end-start
