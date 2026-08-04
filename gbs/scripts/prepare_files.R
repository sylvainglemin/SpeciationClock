# Script to combine different files


# List of Fis and selfing rates



df_summary <- data.frame()

for(thresh in c(0.05,0.001,1e-04,1e-06)){
  list_estimates <- list.files(path = "output/" ,pattern = paste(thresh,"_out.txt.txt",sep=""))
  for(myfile in list_estimates) {
    estimate_file <- paste("output/",myfile,sep="")
    df <- read.table(estimate_file,header=T)
    pop <- names(df)[grep(pattern = "Fis",names(df))]
    pop <- sub("Fis_","",pop)
    fis_pop <- df[,grep(pattern = "Fis",names(df))]
    self_pop <- df[,grep(pattern = "Self",names(df))]
    species <- rep(strsplit(myfile,"_")[[1]][[1]], length(pop))
    threshold <- rep(thresh, length(pop))
    df_summary <- rbind(df_summary,data.frame(species,pop,threshold,t(fis_pop),t(self_pop)))
  }
}


list_estimates <- list.files(path = "output/", pattern = "nofilter")
for(myfile in list_estimates) {
    estimate_file <- paste("output/",myfile,sep="")
    df <- read.table(estimate_file,header=T)
    pop <- names(df)[grep(pattern = "Fis",names(df))]
    pop <- sub("Fis_","",pop)
    fis_pop <- df[,grep(pattern = "Fis",names(df))]
    self_pop <- df[,grep(pattern = "Self",names(df))]
    species <- rep(strsplit(myfile,"_")[[1]][[1]], length(pop))
    threshold <- rep(0, length(pop))
    df_summary <- rbind(df_summary,data.frame(species,pop,threshold,t(fis_pop),t(self_pop)))
  }


names(df_summary) <- c("Species","Pop","Thresh","Fis","Self")

write.table(df_summary,"output/summary_selfingrate.txt",row.names = F,quote = F)


plot(df_summary[df_summary$Thresh==0.05,]$Self,df_summary[df_summary$Thresh==0,]$Self)
