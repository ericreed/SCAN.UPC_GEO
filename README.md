# SCAN.UPC_GEO

This repository contains functions to download data files from GEO and implement the UPC_Generic function.

The current workflow is as follows:
  1.  Use geo2bw() function in getGEO_samples to write a bigWig file from .wig file in GEO.
    - Right now the .wig file needs documentation as to what genome build it is.  To do this I download the .soft file from GEO.         If no build is specified the bigWig file is not creates.  In this example, the GSM# file contains two .wig files, of            which only 1 has a specified build
  2.  Use UPC_bigWig to run UPC on bigWig file, specifying the desired bin size and number of base pairs to overlap on either side.
