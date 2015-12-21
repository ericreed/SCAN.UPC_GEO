#Function to read in and create .bw files
require(SCAN.UPC)
require(rtracklayer)
require(caTools)

#Set base directly to read and write files
baseDir<-getwd()
GSM<-"GSM772735"
BWdir<-file.path(baseDir, GSM)
BWfile<-list.files(BWdir, pattern=".bw", full.names = TRUE)[1]

#Read in promoter annotation
regDir<-"/Users/ericreed/Desktop/JohnsonLab/SCAN.UPC_GEO/EnsReg"
prom<-read.table(file.path(regDir, "Promoters_EnsReg_hg19.txt"), header=T, stringsAsFactors = FALSE)

#Create GenomicRanges object
chrU<-unique(prom$chr)

chrTab<-as.data.frame(table(prom$chr))
colnames(chrTab)<-c("chrs", "lengths")
chrTab<-chrTab[order(match(chrTab$chrs, chrU)),]

annot <- GRanges(seqnames =
                   Rle(chrTab$chrs, chrTab$lengths),
                 ranges =
                   IRanges(prom$start, prom$end))
values(annot)$id <- prom$name

#Alternate function to run UPC on BigWig files.  This will bin all data file from the same genome build the same.
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

start<-Sys.time()
UPCout<-UPC_bigWig_annot(BWfile, annot, max=FALSE)
end<-Sys.time()
end-start
