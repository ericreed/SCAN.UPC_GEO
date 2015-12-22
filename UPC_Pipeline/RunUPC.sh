#!/bin/bash
#$ -S /bin/bash
#$ -N RunUPC
#$ -o logs/
#$ -e logs/
#$ -cwd
#$ -j y
#$ -V
#$ -m e
#$ -l h_rt=120:00:00

module load R/R-3.1.1
runIt="R CMD BATCH --no-save --no-restore '--args $GSM $annot $baseDir $Seq $build ' RunUPC.R Rout/RunUPC_$GSM.Rout"
echo $runIt
eval $runIt
