# Script to analyses SpeciationClock GBS data
# Filter potentially erroneous SNPs and compute Fstatistics and selfing rate
# Also return the list of SNPs with their probability of erroneous mapping / based on heterozygosity excess

# Required R packages:
# vcfR: to read vcf
# adegenet: to manipulate files
# reshape2: to manipulate files
# hierfstat: to compute Fis (and other Fstat)
# inbreedR: to compute idendity disequilibrium statistics


packages <- c("vcfR","adegenet","reshape2","hierfstat","ggplot2")
A <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


# Functions to estimate selfing rates ####

# From the Fis
self_from_fis <- function(f) {
  2*f/(1+f)
}
# From the identity disequilibrium statistics (g2)
self_from_g2 <- function(g) {
  (1+5*g -sqrt(1+10*g+9*g^2))/(2*g)
}
# To compute Dxy from Hs and Ht
ff <- function(x) x*(x-1)/2
dxy <- function(n1,n2,H1,H2,Ht){
  n <- n1+n2
  d <- (ff(n)*Ht - ff(n1)*H1 -ff(n2)*H2) / (n1*n2)
  return(d)
}


# MAIN FUNCTION ####

# REMOVEFULLHET: to chose how to handle fully heterozygote sites: yes, no, or automatic: yes if there is an excess of the last category, no otherwise
# NMIN: minimum number of SNPs to do fdr correction, otherwise bonferonni
genpop_estimate <- function(FILE,OUT="output.txt",REMOVEFULLHET="automatic",THRESHOLD=0.05,CORR=T,NMIN = 1000) {
  
  print("Preparing files")
  # Read vcf file
  myvcf <- read.vcfR(FILE,verbose = F)
  # Extract genotype information
  gt <- extract.gt(myvcf)
  # extract popname
  # get the name of each individual, split the name according to "_", remove the two first part (individual name), fuse the terms whcih are left (corresponding to pop name)
  popname <- sapply(names(gt[1,]), function(x) apply(as.matrix(strsplit(x,split = "_")[[1]][-c(1,2)]),2,paste,collapse="_") )
  # Format datset for hierFstat
  genotype <- vcfR2genind(myvcf,return.alleles = T) # Convert into the genind format (adegent package)
  pop(genotype) <- popname # Add the pop names
  pop_label <- levels(genotype$pop)
  
  
  # Statistics on populations
  nb_pop <- length(pop_label) # Number of populations
  nb_ind <- summary(genotype$pop) # Number of individuals per pop
  sum_pop <- basic.stats(genotype)
  Hobs <- Ho(genotype)
  Hexp <- Hs(genotype)
  Fis <- 1 - Hobs/Hexp
  if(nb_pop==1){
    Fis <- Fis[1] # To avoid "dumpop" when there is only one pop
    names(Fis) <- pop_label
  }
 
  print("Checking Heterozygote distribution")
  # Plot of the raw Fis distributions for checking and remove fully heterozygote sites if required
  p_value <- rep(1,length(genotype$loc.n.all)) # Initialization of p_value
  for(p in pop_label) {
    geno <- genotype[which(popname==p),]
    geno_sum <- summary(geno)
    nind <- geno_sum$n.by.pop
    npol <- length(which(geno_sum$Hobs>=0))
    # Print histogram of % of Het obs
    #G <- ggplot(as.data.frame(geno_sum$Hobs),aes(x=geno_sum$Hobs)) + geom_histogram()
    #G <- G + scale_y_log10() + xlab("% of heterozygote individuals") + xlim(-0.1,1.1) + ggtitle(paste("Raw data - Population: ",p,sep=""))
    #pdf(paste("Raw_Hobs_distribution_",p,".pdf",sep=""),width = 12,height = 8)
    #  print(G)
    #dev.off()
 
    # Checking the distribution of heterozygote if the option of removing them is automatic
    histHo <- hist(geno_sum$Hobs,plot = F,breaks = 30)$counts
    #
    if(REMOVEFULLHET=="automatic") {
      histHo <- hist(geno_sum$Hobs,plot = F,breaks= nind*2)$counts
      n <- length(histHo)
      if(histHo[n]>2*histHo[n-1]) remove <- T else remove <- F
    }
    if(REMOVEFULLHET=="yes") remove <- T
    if(REMOVEFULLHET=="no") remove <- F
    
    if(remove){
      tokeep <- which(geno_sum$Hobs<1)
      newgeno <- geno[,sort(c(2*tokeep-1,2*tokeep))]
      if(Hs(newgeno)[1]==0) Fis[p] <- NA
      else Fis[p] <- 1 - Ho(newgeno)[1]/Hs(newgeno)[1]
      # P_value=0 for fully heterozygote sites
      p_value[which(geno_sum$Hobs==1)] <- 0
    }
  }

  print("Filtering SNPs")
  #Initialization
  STOP <- F
  # Loop
  while(STOP==F) {
    # Computation of p_values
    for(p in pop_label) {
      geno <- genotype[which(popname==p),]
      geno_sum <- summary(geno)
      npol <- length(which(geno_sum$Hobs>=0))
      n <- apply(geno$tab,2,function(x) sum(!is.na(x)))[c(T,F)] # number of non missing individuals, should take only one value per locus as there is  a value for each allele
      n_ho <- n*geno_sum$Hobs # number of observed homozygotes
      he <- geno_sum$Hexp*(1-Fis[p])
      if(CORR) {
        if(npol > NMIN) p_value <- p_value*p.adjust(dbinom(n_ho,n,he),method = "fdr")
        else p_value <- p_value*p.adjust(dbinom(n_ho,n,he),method = "bon")
      } else p_value <- p_value*dbinom(n_ho,n,he)
    }
    select <- which(p_value>=THRESHOLD)
    if(length(select)>0) {
      newgenotype <- genotype[,sort(c(2*select-1,2*select))]
      newFis <- 1 - Ho(newgenotype)/Hs(newgenotype)
      if(nb_pop==1){
        newFis <- newFis[1] # To avoid "dumpop" when there is only one pop
        names(newFis) <- pop_label
      }
      print(newFis)
      if(sum(newFis-Fis)==0) STOP <- T
      Fis <- newFis
    } else STOP <- T
  }
  
  # Analysis of the new dataset without out-filtered SNPs
  
  print("Estimation")
  # Re-estimation of parameters on the filtered dataset
  Hobs <- Ho(newgenotype)
  Hexp <- Hs(newgenotype)
  Fis <- 1 - Hobs/Hexp
  if(nb_pop==1){
    Hobs <- Hobs[1]
    names(Hobs) <- pop_label
    Hexp <- Hexp[1]
    names(Hexp) <- pop_label
    Fis <- Fis[1] # To avoid "dumpop" when there is only one pop
    names(Fis) <- pop_label
  }
  Self <- self_from_fis(Fis)
  if(nb_pop>1) {
    Fst <- as.vector(genet.dist(newgenotype,method = "WC84"))
    # Home made script to compute Dxy on all pairs
    Dxy <- c()
    for(i in 1:(nb_pop-1)) {
      for(j in (i+1):nb_pop){
        geno <- newgenotype[which(popname==pop_label[i] | popname==pop_label[j]),]
        bs <- basic.stats(geno)$overall
        Ht <- bs["Ht"]
        He <- Hs(geno)
        nb <- nb_ind[c(i,j)]
        Dxy <- c(Dxy,dxy(nb[1],nb[2],He[1],He[2],Ht))
      }
    }
  } else {
    Fst <- c()
    Dxy <- c()
  }
  
  

  print("Writing outputs")
  # Write results
  fst_name <- c()
  dxy_name <- c()
  if(nb_pop > 1) {
    for(i in 1:(nb_pop-1)) {
      for(j in (i+1):nb_pop){
        fst_name <- c(fst_name,paste("Fst",pop_label[i],pop_label[j],sep ="_"))
        dxy_name <- c(dxy_name,paste("Dxy",pop_label[i],pop_label[j],sep ="_"))
      }
    }
  }
  
  
  var_names <- c("dataset",
                 paste("Ho",pop_label,sep ="_"),
                 paste("He",pop_label,sep ="_"),
                 paste("Fis",pop_label,sep ="_"),
                 paste("Self",pop_label,sep ="_"),
                 fst_name,
                 dxy_name
                 )
  results <- c(FILE,Hobs,Hexp,Fis,Self,Fst,Dxy)
  
  write.table(rbind(var_names,results),file = paste(OUT,".txt",sep=""),sep = "\t",quote=F,row.names = F,col.names = F)
  write.table(cbind(names(gt[,1]),p_value[c(T,F)]),file = paste(OUT,"SNP_pvalues.txt",sep="_"),sep = "\t",quote=F,row.names = F,col.names = c("snp","p_value"))
  
  # To check by printing histograms
  # for(p in pop_label) {
  #   geno <- newgenotype[which(popname==p),]
  #   geno_sum <- summary(geno)
  #   # Print histogram of het obs
  #   G <- ggplot(as.data.frame(geno_sum$Hobs),aes(x=geno_sum$Hobs)) + geom_histogram()
  #   G <- G + scale_y_log10() + xlab("% of heterozygote individuals") + xlim(-0.1,1.1) + ggtitle(paste("Filtered data - Population: ",p,sep=""))
  #   pdf(paste("Filtered_Hobs_distribution_",p,".pdf",sep=""))
  #     print(G)
  #   dev.off()
  # }
  
  print("DONE!")
  
}


tryCatch(1/q,error = function(e) NA)
