# Script to run the snp_filtering function on all vcf files


# Source the R script

source("scripts/snp_filtering.R")

listvcf <- list.files(path = "data/vcf/" ,pattern = ".vcf")
thresh <- 0 # 0.05 0.001 10^(-6)

for(myfile in listvcf) {
  print(myfile)
  vcf_file <- paste("data/vcf/",myfile,sep="")
  #output <- paste("output/",strsplit(myfile,split = "\\.")[[1]][1],"_thresh=",thresh,"_out.txt",sep="")
  output <- paste("output/",strsplit(myfile,split = "\\.")[[1]][1],"nofilter_out.txt",sep="")
  genpop_estimate(FILE = vcf_file,REMOVEFULLHET = "no",OUT = output,THRESHOLD = thresh, CORR=T)
  }
