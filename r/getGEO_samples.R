#Function to read in and create .bw files
require(SCAN.UPC)
require(rtracklayer)
require(tools)

#Set base directly to read and write files
baseDir<-getwd()

GSM<-"GSM772735"

geo2bw<-function(GSM, baseDir, Seq = NULL){
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
    if(is.null(Seq)){
      if(grepl("hg19", build)){
        build<-"hg19"
        Seq<-Seqinfo(genome=build)}
      if(grepl("hg18", build)){
        build<-"hg18"
        Seq<-Seqinfo(genome=build)}
      if(grepl("hg38", build)){
        build<-"hg38"
        Seq<-Seqinfo(genome=build)}}
    
    if(!is.null(Seq)){
    wigToBigWig(file.path(newDir, i), seqinfo = Seq, dest = paste(file_path_sans_ext(file.path(newDir, i), TRUE),"_", build, ".bw", sep=""))} else {
      print(paste("Could not create bigWig file from, ", i," because genome build is unknown."))}
  }
}
  
  
geo2bw(GSM, baseDir)
  
  