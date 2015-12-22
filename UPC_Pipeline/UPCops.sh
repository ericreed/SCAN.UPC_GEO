#!/bin/bash

######To Run
##chmod u+x UPCops.sh
##./UPCops.sh

#########################
####List of GSM#'s to run
# # # Space separated!
#########################
GEOacc=(GSM772724 GSM772968 GSM772974)
#########################

####Set e-mail address to alert when job is finished
email="reeder@bu.edu"

####Set Annotation File
annotFile="Promoters_EnsReg_hg19.txt"

####Specify the build of the annotation file
build="hg19"

####Set Directory for writing files
baseDir="/restricted/projectnb/combat/UPC_EGRM/data/UPCtest"

####Set file that has the genome build information
seqInfo="GenomBuilds.RData"


for i in ${GEOacc[@]}; do runMult="qsub -P combat -M $email -pe single_node 8 -v GSM=$i -v annot=$annotFile -v build=$build -v baseDir=$baseDir -v Seq=$seqInfo RunUPC.sh"; eval $runMult; done