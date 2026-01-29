# Script to extract the number and length of contigs for each species

setwd("data/contigs")

list_sp <- list.files(pattern = "CONTINGS")
list_sp <- sapply(list_sp,function(x) strsplit(x,"_CONTINGS.txt")[[1]])

tab <- data.frame()
  
for(SP in list_sp){
  sp_name <- strsplit(SP,"_")[[1]][1:2]
  df1 <- read.table(paste(SP,"_CONTINGS.txt",sep=""),comment.char = "")
  df1$contig_id <- sapply(df1$V1,
         function(x) {
           strsplit(strsplit(x,"=")[[1]][[3]],",")[[1]][[1]]
         }
  )
  df1$contig_length <- sapply(df1$V1,
                          function(x) {
                            strsplit(strsplit(x,"=")[[1]][[4]],">")[[1]][[1]]
                          }
  )
  df1 <- df1[c(2,3)]
  df2 <- read.table(paste("white_list/",SP,".whitelist",sep=""))
  names(df2) <- c("contig_id","snp_position")
  df <- merge(df1,df2,by="contig_id")  
  Ntot <- length(df1$contig_id)
  Ltot <- sum(as.numeric(df1$contig_length))
  Nfiltered <- length(df$contig_id)
  Lfiltered <- sum(as.numeric(df$contig_length))
  tab <- rbind(tab, c(sp_name,Ntot,Ltot,Nfiltered,Lfiltered))
  print(SP)
}

names(tab) <- names(tab) <- c("Genus","Species","Ntot","Ltot","Nfiltered","Lfiltered")


write.table(tab,"info_species.txt",quote = F,row.names = F)
