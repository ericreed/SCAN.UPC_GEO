require(R.utils, quietly = TRUE)

args=(commandArgs(TRUE))

#qsub -P combat -pe single_node 8 -v GSM="GSM772735" -v annot="Promoters_EnsReg_hg19.txt" -v baseDir="/restricted/projectnb/combat/UPC_EGRM/data/UPCtest" -v Seq="GenomBuilds.RData" -v build="hg19" RunUPC.sh

GSM = args[[1]]
annotFile  = args[[2]]
baseDir = args[[3]]
seqData = args[[4]]
build = args[[5]]

print(GSM)

####Run UPC on Input
suppressMessages(source("SCAN.UPC_Functions.R"))


test<-try(seqData, silent = TRUE)
if(class(test)!="try-error"){
  load(seqData)
  SeqList<-list(Seq19, Seq18, Seq38)
  names(SeqList)<-c("Seq19", "Seq18", "Seq38")}

#Read in annotations of EndReg file
prom<-read.table(annotFile, header=T, stringsAsFactors = FALSE)
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


wrapFunc<-function(GSM, baseDir, test){
  if(class(test)=="try-error"){
    geo2bw(GSM, baseDir)} else {
      geo2bw(GSM, baseDir, SeqList = SeqList)}
  BWdir<-file.path(baseDir, GSM)
  BWfile<-list.files(BWdir, pattern=".bw", full.names = TRUE)[1]
  BWfile<-BWfile[grepl(build, BWfile)]
  if(!is.na(BWfile)){
    UPCout<-UPC_bigWig_annot(BWfile, annot, max=FALSE)
    write.table(UPCout, paste(file_path_sans_ext(BWfile), "_UPCout.txt", sep=""), row.names=FALSE, col.names = TRUE, quote=FALSE)
  }}

#######NEED TO ADD BINNED AVERAGE FUNCTION
binnedAverage <- function (bins, numvar, varname) 
{
  if (!is(bins, "GRanges")) 
    stop("'x' must be a GRanges object")
  if (!is(numvar, "RleList")) 
    stop("'numvar' must be an RleList object")
  if (!identical(seqlevels(bins), names(numvar))) 
    stop("'seqlevels(bin)' and 'names(numvar)' must be identical")
  viewMeans2 <- function(v) {
    means <- viewMeans(v)
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
#######


UPCout<-wrapFunc(GSM, baseDir, test)


