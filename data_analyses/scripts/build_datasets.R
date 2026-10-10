# Script supporting the data analysis presented in Gizaw et al.
# Author: Lucas Marie-Orleach / modified Sylvain Glémin
# Date: April 2024 / September 2024

## PACKAGES ####
library(readr) #for read.delim
library(dplyr) #for left_join

## LOADINGS ####


# Correspondence table for the different ways of naming species
cortable <- cbind(
  c("AAU","AT","Aal","BE","CB","CG","CH","CM","DN","EV","HI","HP","LC","LH","MS","MT","RP","RR","SA","SAC","SH","SU"),
  c("Arabis_auriculata","Arabis_turrita","Arabis_alpina","Bunias_erucago","Cardamine_bellidifolia","Cochlearia_groenlandica","Cardamine_hirsuta","Cakile_maritima","Draba_nivalis","Eruca_vesicaria","Hirschfeldia_incana","Hornungia_petraea","Lepidium_campestre","Lepidium_hirtum","Matthiola_sinuata","Matthiola_tricuspidata","Ranunculus_pygmaeus","Rapistrum_rugosum","Sinapis_arvensis","Silene_acaulis","Saxifraga_hyperborea","Silene_uralensis"),
  c("Med","Med","Med","Med","Arc","Arc","Med","Med","Arc","Med","Med","Med","Med","Med","Med","Med","Arc","Med","Med","Arc","Arc","Arc")
)


# Parental Dataset
Med_ParData <- read_delim("data/experiment_data/Med_ParentalData.txt","\t")
Med_ParData$Species <- gsub("[()]","",Med_ParData$Species)
Med_ParData$Species <- gsub(" autonomous","",Med_ParData$Species)
Med_ParData$Species <- gsub(" ","_",Med_ParData$Species)
Med_ParData$Species <- gsub("Cakile_maritimaa愼㸰", "Cakile_maritima", Med_ParData$Species)
Med_ParData$Species <- gsub("Arabis_auricculata", "Arabis_auriculata", Med_ParData$Species)
Med_ParData$Species <- gsub("Sinapis_avensis", "Sinapis_arvensis", Med_ParData$Species)
colnames(Med_ParData)[4] <- "pop"
colnames(Med_ParData)[7] <- "TotalSeed"
colnames(Med_ParData)[8] <- "GerminatedSeed"
colnames(Med_ParData)[12] <- "SelfRate"
Med_ParData$PollenFert <- Med_ParData$FertilePollen/Med_ParData$TotalPollen
Med_ParData$SeedGerm <- Med_ParData$GerminatedSeed/Med_ParData$TotalSeed

Arc_ParData <- subset(read_delim("data/experiment_data/Arc_Data.txt","\t"),Cross_type=="PP" | Cross_type=="WPC")
Arc_ParData$pop <- substr(Arc_ParData$Maternal_pop,6,6)
Arc_ParData$Species <- gsub("C_bellidifolia", "Cardamine_bellidifolia", Arc_ParData$Species)
Arc_ParData$Species <- gsub("C_groenlandica", "Cochlearia_groenlandica", Arc_ParData$Species)
Arc_ParData$Species <- gsub("R_pygmaeus", "Ranunculus_pygmaeus", Arc_ParData$Species)
Arc_ParData$Species <- gsub("S_acaulis", "Silene_acaulis", Arc_ParData$Species)
Arc_ParData$Species <- gsub("S_hyperborea", "Saxifraga_hyperborea", Arc_ParData$Species)
Arc_ParData$Species <- gsub("S_uralensis", "Silene_uralensis", Arc_ParData$Species)
colnames(Arc_ParData)[10] <- "GerminatedSeed"
colnames(Arc_ParData)[11] <- "NonGerminatedSeed"
colnames(Arc_ParData)[13] <- "FertilePollen"
colnames(Arc_ParData)[14] <- "UnFertilePollen"
Arc_ParData$TotalPollen <- Arc_ParData$FertilePollen+Arc_ParData$UnFertilePollen
Arc_ParData$TotalSeed <- Arc_ParData$GerminatedSeed+Arc_ParData$NonGerminatedSeed
Arc_ParData$PollenFert <- Arc_ParData$FertilePollen/Arc_ParData$TotalPollen
Arc_ParData$SeedGerm <- Arc_ParData$GerminatedSeed/Arc_ParData$TotalSeed
Arc_GBSData_pop <- read_delim("outputs/gbs_results/Arc_GBSData_pop.txt","\t")
Arc_GBSData_pop$Species <- gsub("C_bellidifolia", "Cardamine_bellidifolia", Arc_GBSData_pop$Species)
Arc_GBSData_pop$Species <- gsub("C_groenlandica", "Cochlearia_groenlandica", Arc_GBSData_pop$Species)
Arc_GBSData_pop$Species <- gsub("R_pygmaeus", "Ranunculus_pygmaeus", Arc_GBSData_pop$Species)
Arc_GBSData_pop$Species <- gsub("S_acaulis", "Silene_acaulis", Arc_GBSData_pop$Species)
Arc_GBSData_pop$Species <- gsub("S_hyperborea", "Saxifraga_hyperborea", Arc_GBSData_pop$Species)
Arc_GBSData_pop$Species <- gsub("S_uralensis", "Silene_uralensis", Arc_GBSData_pop$Species)
Arc_GBSData_pop$pop <- substr(Arc_GBSData_pop$Pop,1,1)
for (i in 1:nrow(Arc_GBSData_pop)) {if(Arc_GBSData_pop$pop[i]=="O") {Arc_GBSData_pop$pop[i]<-"N"}}
Arc_GBSData_pop <- Arc_GBSData_pop[-c(1,5,12),] # REMOVE DUPLICATES IN THE GBS DATASET
Arc_ParData <- left_join(Arc_ParData, Arc_GBSData_pop[,c(1,3,4,5)], by=c("Species","pop"))

Dra_RawParData <- read_csv("data/experiment_data/Dra_ParentalData.txt")
Dra_RawHybData <- read_delim("data/experiment_data/Dra_HybridData.txt","\t")
Dra_RawHybData[, c(15:18)] <- apply(Dra_RawHybData[, c(15:18)], 2, function(x) as.numeric(as.character(x)))
Species <- as.data.frame(c(Dra_RawParData$Species, subset(Dra_RawHybData,Cross_type=="WPC")$Species))
FertilePollen <- as.data.frame(c(Dra_RawParData$FertilePollen, subset(Dra_RawHybData,Cross_type=="WPC")$F1_fertile_pollen))
TotalPollen <- as.data.frame(c(Dra_RawParData$FertilePollen+Dra_RawParData$InfertilePollen, subset(Dra_RawHybData,Cross_type=="WPC")$F1_fertile_pollen+subset(Dra_RawHybData,Cross_type=="WPC")$F1_unfertile_pollen))
GerminatedSeed <- as.data.frame(c(rep("NA",nrow(Dra_RawParData)), subset(Dra_RawHybData,Cross_type=="WPC")$F1_total_no_of_seeds-subset(Dra_RawHybData,Cross_type=="WPC")$F1_no_of_undeveloped_ovules))
TotalSeed <- as.data.frame(c(rep("NA",nrow(Dra_RawParData)), subset(Dra_RawHybData,Cross_type=="WPC")$F1_total_no_of_seeds))

Dra_ParData <- data.frame(matrix(nrow = nrow(Species), ncol=5))
Dra_ParData <- cbind(Species, FertilePollen, TotalPollen, GerminatedSeed, TotalSeed)
Dra_ParData[, c(2:5)] <- apply(Dra_ParData[, c(2:5)], 2, function(x) as.numeric(as.character(x)))
colnames(Dra_ParData) <- c("Species","FertilePollen","TotalPollen","GerminatedSeed","TotalSeed")
Dra_ParData$PollenFert <- Dra_ParData$FertilePollen/Dra_ParData$TotalPollen
Dra_ParData$SeedGerm <- Dra_ParData$GerminatedSeed/Dra_ParData$TotalSeed
Dra_GBSData_pop <- read_delim("outputs/gbs_results/Dra_GBSData_pop.txt","\t")[c(1,2,4),] ## NEED TO CHECK WITH Ch, Ab, Lo
Dra_ParData <- left_join(Dra_ParData, Dra_GBSData_pop[,c(1,3,4)], by=c("Species"))
Dra_ParData$Species <- gsub("D_fladnizensis", "Draba_fladnizensis", Dra_ParData$Species)
Dra_ParData$Species <- gsub("D_nivalis", "Draba_nivalis", Dra_ParData$Species)
Dra_ParData$Species <- gsub("D_subcapitata", "Draba_subcapitata", Dra_ParData$Species)
Dra_ParData$pop <- NA

ParData <- rbind(
  Med_ParData[,c("Species","pop","Fis","SelfRate","FertilePollen","TotalPollen","PollenFert","GerminatedSeed","TotalSeed","SeedGerm")], 
  Arc_ParData[,c("Species","pop","Fis","SelfRate","FertilePollen","TotalPollen","PollenFert","GerminatedSeed","TotalSeed","SeedGerm")],
  Dra_ParData[,c("Species","pop","Fis","SelfRate","FertilePollen","TotalPollen","PollenFert","GerminatedSeed","TotalSeed","SeedGerm")]
)
ParData$Dataset <- c(rep("Med",times=nrow(Med_ParData)),rep("Arc",times=nrow(Arc_ParData)),rep("Dra",times=nrow(Dra_ParData)))

SummaryPar_pop <- ParData %>% dplyr::group_by(Species,pop) %>% 
  dplyr::summarise(
    Fis=mean(Fis),
    SelfRate=mean(SelfRate),
    PollenFert=mean(FertilePollen/TotalPollen, na.rm=T),
    SeedGerm=mean(GerminatedSeed/TotalSeed, na.rm=T)
  )

SummaryPar_sp <- SummaryPar_pop %>% 
  dplyr::group_by(Species) %>% 
  dplyr::summarise(
    Fis=mean(Fis),
    SelfRate=mean(SelfRate),
    PollenFert=mean(PollenFert, na.rm=T),
    SeedGerm=mean(SeedGerm, na.rm=T)
  )

##################### #
## HYBRID DATASET ####
################### #

Med_HybData <- read_delim("data/experiment_data/Med_HybridData.txt","\t")
Med_HybData$Species <- gsub(" ","_",Med_HybData$Species)
Med_HybData$Species <- gsub("\\.","",Med_HybData$Species)
Med_HybData$Species <- gsub("\\Cakile\\w+","Cakile_maritima",Med_HybData$Species)
Med_HybData <- subset(Med_HybData,Species!="Lepidium_campestre" & Species!="Lepidium_hirtum")
colnames(Med_HybData)[3] <- "MatPop"
colnames(Med_HybData)[4] <- "PatPop"
colnames(Med_HybData)[5] <- "MatPop_PollenFert"
colnames(Med_HybData)[6] <- "PatPop_PollenFert"
colnames(Med_HybData)[7] <- "MatPop_SeedGerm"
colnames(Med_HybData)[8] <- "PatPop_SeedGerm"
colnames(Med_HybData)[9] <- "Fst"
colnames(Med_HybData)[15] <- "MatPop_SelfRate"
colnames(Med_HybData)[19] <- "PatPop_SelfRate"
colnames(Med_HybData)[20] <- "FertilePollen"
colnames(Med_HybData)[22] <- "TotalSeed" ##GerminatedSeed and TotalSownSeed should be swapped
colnames(Med_HybData)[23] <- "GerminatedSeed"
Med_HybData$UnFertilePollen <- Med_HybData$TotalPollen - Med_HybData$FertilePollen
Med_HybData$NonGerminatedSeed <- Med_HybData$TotalSeed - Med_HybData$GerminatedSeed
Med_HybData$MatPop_PollenFert <- pull(left_join(Med_HybData, SummaryPar_pop[,c(1,2,5)],by=c("Species","MatPop"="pop"))[,ncol(Med_HybData)+1])
Med_HybData$PatPop_PollenFert <- pull(left_join(Med_HybData, SummaryPar_pop[,c(1,2,5)],by=c("Species","PatPop"="pop"))[,ncol(Med_HybData)+1])
Med_HybData$MatPop_SeedGerm <- pull(left_join(Med_HybData, SummaryPar_pop[,c(1,2,6)],by=c("Species","MatPop"="pop"))[,ncol(Med_HybData)+1])
Med_HybData$PatPop_SeedGerm <- pull(left_join(Med_HybData, SummaryPar_pop[,c(1,2,6)],by=c("Species","PatPop"="pop"))[,ncol(Med_HybData)+1])

Arc_HybData <- subset(read_delim("data/experiment_data/Arc_Data.txt","\t"), Cross_type=="BRC")
Arc_HybData$Species <- gsub("C_bellidifolia", "Cardamine_bellidifolia", Arc_HybData$Species)
Arc_HybData$Species <- gsub("C_groenlandica", "Cochlearia_groenlandica", Arc_HybData$Species)
Arc_HybData$Species <- gsub("R_pygmaeus", "Ranunculus_pygmaeus", Arc_HybData$Species)
Arc_HybData$Species <- gsub("S_acaulis", "Silene_acaulis", Arc_HybData$Species)
Arc_HybData$Species <- gsub("S_hyperborea", "Saxifraga_hyperborea", Arc_HybData$Species)
Arc_HybData$Species <- gsub("S_uralensis", "Silene_uralensis", Arc_HybData$Species)
Arc_HybData$MatPop <- substr(Arc_HybData$Maternal_pop, 6, 6)
Arc_HybData$PatPop <- substr(Arc_HybData$Paternal_pop, 6, 6)
SummaryPar_pop$MatPop <- SummaryPar_pop$pop
SummaryPar_pop$PatPop <- SummaryPar_pop$pop
Arc_HybData$MatPop_PollenFert <- pull(left_join(Arc_HybData, SummaryPar_pop[,c("Species","MatPop","PollenFert")], by=c("Species","MatPop"))[,ncol(Arc_HybData)+1])
Arc_HybData$PatPop_PollenFert <- pull(left_join(Arc_HybData, SummaryPar_pop[,c("Species","PatPop","PollenFert")], by=c("Species","PatPop"))[,ncol(Arc_HybData)+1])
Arc_HybData$MatPop_SeedGerm   <- pull(left_join(Arc_HybData, SummaryPar_pop[,c("Species","MatPop","SeedGerm")], by=c("Species","MatPop"))[,ncol(Arc_HybData)+1])
Arc_HybData$PatPop_SeedGerm   <- pull(left_join(Arc_HybData, SummaryPar_pop[,c("Species","PatPop","SeedGerm")], by=c("Species","PatPop"))[,ncol(Arc_HybData)+1])
Arc_HybData$MatPop_SelfRate   <- pull(left_join(Arc_HybData, SummaryPar_pop[,c("Species","MatPop","SelfRate")], by=c("Species","MatPop"))[,ncol(Arc_HybData)+1])
Arc_HybData$PatPop_SelfRate   <- pull(left_join(Arc_HybData, SummaryPar_pop[,c("Species","PatPop","SelfRate")], by=c("Species","PatPop"))[,ncol(Arc_HybData)+1])
Arc_GBSData_PopPair <- read_delim("outputs/gbs_results/Arc_GBSData_PopPair.txt")
Arc_GBSData_PopPair <- Arc_GBSData_PopPair[-c(9,10,54,55,56,57,63,64),-c(6:9)] # manually removing redundant genomic data
Arc_GBSData_PopPair$Species <- gsub("C_bellidifolia", "Cardamine_bellidifolia", Arc_GBSData_PopPair$Species)
Arc_GBSData_PopPair$Species <- gsub("C_groenlandica", "Cochlearia_groenlandica", Arc_GBSData_PopPair$Species)
Arc_GBSData_PopPair$Species <- gsub("R_pygmaeus", "Ranunculus_pygmaeus", Arc_GBSData_PopPair$Species)
Arc_GBSData_PopPair$Species <- gsub("S_acaulis", "Silene_acaulis", Arc_GBSData_PopPair$Species)
Arc_GBSData_PopPair$Species <- gsub("S_hyperborea", "Saxifraga_hyperborea", Arc_GBSData_PopPair$Species)
Arc_GBSData_PopPair$Species <- gsub("S_uralensis", "Silene_uralensis", Arc_GBSData_PopPair$Species)
Arc_HybData$Cross <- paste(Arc_HybData$MatPop, "x", Arc_HybData$PatPop, sep="")
Arc_HybData$New_Mean_Dxy <- as.numeric(unlist(left_join(Arc_HybData, Arc_GBSData_PopPair[,c("Species","Cross","New_Mean_Dxy")], by=c("Species","Cross"))[,ncol(Arc_HybData)+1]))
Arc_HybData$Fst <- as.numeric(unlist(left_join(Arc_HybData, Arc_GBSData_PopPair[,c("Species","Cross","Fst")], by=c("Species","Cross"))[,ncol(Arc_HybData)+1]))
colnames(Arc_HybData)[13] <- "FertilePollen"
colnames(Arc_HybData)[14] <- "UnFertilePollen"
Arc_HybData$TotalPollen <- Arc_HybData$FertilePollen + Arc_HybData$UnFertilePollen
colnames(Arc_HybData)[10] <- "GerminatedSeed"
colnames(Arc_HybData)[11] <- "NonGerminatedSeed"
Arc_HybData$TotalSeed <- Arc_HybData$GerminatedSeed + Arc_HybData$NonGerminatedSeed

Dra_HybData <- subset(read_delim("data/experiment_data/Dra_HybridData.txt","\t"), Cross_type!="WPC" & FST!=0)
Dra_HybData[, c(15:18)] <- apply(Dra_HybData[, c(15:18)], 2, function(x) as.numeric(as.character(x)))
Dra_HybData$Species <- gsub("D_fladnizensis", "Draba_fladnizensis", Dra_HybData$Species)
Dra_HybData$Species <- gsub("D_nivalis", "Draba_nivalis", Dra_HybData$Species)
Dra_HybData$Species <- gsub("D_subcapitata", "Draba_subcapitata", Dra_HybData$Species)
colnames(Dra_HybData)[14] <- "Fst"
colnames(Dra_HybData)[15] <- "FertilePollen"
colnames(Dra_HybData)[16] <- "UnFertilePollen"
Dra_HybData$TotalPollen <- Dra_HybData$FertilePollen + Dra_HybData$UnFertilePollen
colnames(Dra_HybData)[17] <- "GerminatedSeed"
colnames(Dra_HybData)[18] <- "NonGerminatedSeed"
Dra_HybData$TotalSeed <- Dra_HybData$GerminatedSeed + Dra_HybData$NonGerminatedSeed
Dra_HybData <- left_join(Dra_HybData, SummaryPar_sp[,c(1,3,4,5)], by=c("Species")) #add parental SelfRate, PollenFert, SeedGerm
Dra_HybData <- Dra_HybData %>% rename("ParentalFertility"="PollenFert", "ParentalGermination"="SeedGerm")

HybData <- as.data.frame(matrix(ncol=16, nrow=nrow(Med_HybData)+nrow(Arc_HybData)+nrow(Dra_HybData)))
colnames(HybData) = c("Species","MatPop","PatPop","MatPop_SelfRate","PatPop_SelfRate","MatPop_PollenFert","PatPop_PollenFert","MatPop_SeedGerm","PatPop_SeedGerm","Fst", "New_Mean_Dxy","FertilePollen","UnFertilePollen","GerminatedSeed","NonGerminatedSeed","Dataset")
HybData$Species <- c(Med_HybData$Species, Arc_HybData$Species, Dra_HybData$Species)
HybData$MatPop <- c(Med_HybData$MatPop, Arc_HybData$MatPop, rep("NA", nrow(Dra_HybData)))
HybData$PatPop <- c(Med_HybData$PatPop, Arc_HybData$PatPop, rep("NA", nrow(Dra_HybData)))
HybData$MatPop_SelfRate <- c(Med_HybData$MatPop_SelfRate, Arc_HybData$MatPop_SelfRate, rep("NA", nrow(Dra_HybData)))
HybData$PatPop_SelfRate <- c(Med_HybData$PatPop_SelfRate, Arc_HybData$PatPop_SelfRate, rep("NA", nrow(Dra_HybData)))
HybData$MatPop_PollenFert <- c(Med_HybData$MatPop_PollenFert, Arc_HybData$MatPop_PollenFert, rep("NA", nrow(Dra_HybData)))
HybData$PatPop_PollenFert <- c(Med_HybData$PatPop_PollenFert, Arc_HybData$PatPop_PollenFert, rep("NA", nrow(Dra_HybData)))
HybData$MatPop_SeedGerm <- c(Med_HybData$MatPop_SeedGerm, Arc_HybData$MatPop_SeedGerm, rep("NA", nrow(Dra_HybData)))
HybData$PatPop_SeedGerm <- c(Med_HybData$PatPop_SeedGerm, Arc_HybData$PatPop_SeedGerm, rep("NA", nrow(Dra_HybData)))
HybData$Fst <- c(Med_HybData$Fst, Arc_HybData$Fst, Dra_HybData$Fst)
HybData$Mean_Dxy <- c(Med_HybData$New_Mean_Dxy, Arc_HybData$New_Mean_Dxy, rep("NA", nrow(Dra_HybData)))
HybData$FertilePollen <- c(Med_HybData$FertilePollen, Arc_HybData$FertilePollen, Dra_HybData$FertilePollen)
HybData$UnFertilePollen <- c(Med_HybData$UnFertilePollen, Arc_HybData$UnFertilePollen, Dra_HybData$UnFertilePollen)
HybData$GerminatedSeed <- c(Med_HybData$GerminatedSeed, Arc_HybData$GerminatedSeed, Dra_HybData$GerminatedSeed)
HybData$NonGerminatedSeed <- c(Med_HybData$NonGerminatedSeed, Arc_HybData$NonGerminatedSeed, Dra_HybData$NonGerminatedSeed)
HybData$Dataset <- c(rep("Med",times=nrow(Med_HybData)),rep("Arc",times=nrow(Arc_HybData)),rep("Dra",times=nrow(Dra_HybData)))
HybData[,c(4:15,17)]<-lapply(HybData[,c(4:15,17)],as.numeric)
HybData$PollenFert <- HybData$FertilePollen / (HybData$FertilePollen + HybData$UnFertilePollen)
HybData$SeedGerm <- HybData$GerminatedSeed / (HybData$GerminatedSeed + HybData$NonGerminatedSeed)
HybData$MeanSelfRate <- rowMeans(cbind(as.numeric(unlist(HybData$MatPop_SelfRate)), as.numeric(unlist(HybData$PatPop_SelfRate))), na.rm=T)
HybData[c((1+nrow(HybData)-nrow(Dra_HybData)):nrow(HybData)),]$MeanSelfRate <- Dra_HybData$SelfRate
HybData$MeanParentalFertility <- rowMeans(cbind(as.numeric(unlist(HybData$MatPop_PollenFert)), as.numeric(unlist(HybData$PatPop_PollenFert))), na.rm=T)
HybData[c((1+nrow(HybData)-nrow(Dra_HybData)):nrow(HybData)),]$MeanParentalFertility <- Dra_HybData$ParentalFertility
HybData$MeanParentalGermination <- rowMeans(cbind(as.numeric(unlist(HybData$MatPop_SeedGerm)), as.numeric(unlist(HybData$PatPop_SeedGerm))), na.rm=T)
HybData[c((1+nrow(HybData)-nrow(Dra_HybData)):nrow(HybData)),]$MeanParentalGermination <- Dra_HybData$ParentalGermination
HybData$STDMeanSelfRate = (HybData$MeanSelfRate - mean(HybData$MeanSelfRate, na.rm=T)) / sd(HybData$MeanSelfRate, na.rm=T)
HybData$STDFst = (HybData$Fst - mean(HybData$Fst, na.rm=T)) / sd(HybData$Fst, na.rm=T)
HybData$STDMean_Dxy = (HybData$Mean_Dxy - mean(HybData$Mean_Dxy, na.rm=T)) / sd(HybData$Mean_Dxy, na.rm=T)
HybData$STDMeanParentalFertility = (HybData$MeanParentalFertility - mean(HybData$MeanParentalFertility, na.rm=T)) / sd(HybData$MeanParentalFertility, na.rm=T)
HybData$STDMeanParentalGermination = (HybData$MeanParentalGermination - mean(HybData$MeanParentalGermination, na.rm=T)) / sd(HybData$MeanParentalGermination, na.rm=T)
Summary_HybData <- HybData %>% 
  dplyr::group_by(Dataset, Species, Fst) %>% 
  dplyr::summarise(
    PollenFert=mean(PollenFert, na.rm=T),
    SeedGerm=mean(SeedGerm, na.rm=T),
    SelfRate=mean(MeanSelfRate, na.rm=T),
    Dxy=mean(Mean_Dxy, na.rm=T),
    MeanParentalFertility=mean(MeanParentalFertility, na.rm=T),
    MeanParentalGermination=mean(MeanParentalGermination, na.rm=T),  
    STDMeanSelfRate=mean(STDMeanSelfRate, na.rm=T),
    STDFst=mean(STDFst, na.rm=T),
    STDMeanDxy=mean(STDMean_Dxy, na.rm=T),
    STDMeanParentalFertility=mean(STDMeanParentalFertility, na.rm=T),
    STDMeanParentalGermination=mean(STDMeanParentalGermination, na.rm=T),      
  )

HybData <- subset(HybData, !is.na(FertilePollen) & !is.na(STDFst) & !is.na(STDMeanSelfRate) & !is.na(STDMeanParentalFertility))
HybData$Species_Fst <-  paste(HybData$Species,HybData$Fst,sep = "_")

Summary_HybData <- left_join(Summary_HybData, 
                             subset(HybData, Fst!="NA" & MeanSelfRate!="NA" & PollenFert!="NA") %>% dplyr::group_by(Dataset, Species, Fst) %>% dplyr::summarise(nPollenFert=n()), 
                             by=c("Species","Fst"))

HybData$nPollenFert <- left_join(HybData, Summary_HybData, by=c("Species","Fst"))[,"nPollenFert"]

# Save datasets
write.table(HybData, file="outputs/synthetized_data/SpArc_HybData.txt", row.names = FALSE)
write.table(Summary_HybData, file="outputs/synthetized_data/SpArc_Summary_HybData.txt", row.names = FALSE)


#################### #
## DADI DATASETS ####
################## #


# Standard ∂a∂i analysis - no inbreeding ####

list.results <- list.files(path = "outputs/dadi_results/standard/", pattern = ".out$")

dadi <- data.frame()

for (input in list.results) {
  df <- read.table(paste("outputs/dadi_results/standard/",input,sep=""),header = T)
  df$dof <- c(0,1,2,2,3) 
  df$aic <- 2*df$dof-2*df$lnL
  df$Daic <- df$aic - min(df$aic)
  df$p_value <- 1- pchisq(2*(df$lnL-df$lnL[1]),df$dof)
  df$weight <- exp(-df$Daic/2)
  df$weight <- df$weight/sum(df$weight)
  df$Ttotal <- ifelse(is.na(df$Tinter),df$Tsplit,df$Tsplit+df$Tinter)
  df$Mtotal <- ifelse(is.na(df$Tinter),df$Mrec,(df$Tsplit*df$Manc+df$Tinter*df$Mrec)/(df$Tsplit+df$Ttotal))
  df$sTinter[c(1,2)] <- rep(sqrt(1/12),2)
  av <- apply(df[,-c(1:7)],2,function(x) sum(x*df$weight))
  df <- rbind(df,c(NA,unlist((df[1,c(2:6)])),NA,as.numeric(av)))
  # Here NA is used to keep the vector as numeric.
  # Then the name of the dataset and of the model are added
  df$dataset[6] <- df$dataset[1]
  df$model[6] <- "average"
  dadi <- rbind(dadi,df)
  rm(df)
}

# dataset with average parameters (model averaging procedure based on AIC)
# including also the pure allopatry model as a control
dadi_average <- dadi[dadi$model=="average" | dadi$model=="split_mig",]


# Add the variable species to the dadi datasets
for (i in 1:nrow(dadi_average)) (dadi_average$Species[i] <- cortable[match(strsplit(dadi_average$dataset[i],"_")[[1]][1],cortable[,1]),2])
for (i in 1:nrow(dadi_average)) (dadi_average$Dataset[i] <- cortable[match(strsplit(dadi_average$dataset[i],"_")[[1]][1],cortable[,1]),3])

# Add theta and TotalLength (ie, nb of total SNP) from the CorrectingFactorGBS.txt file, and compute K, the correcting factor
CorrectingFactorGBS <- read.table("outputs/gbs_results/CorrectingFactorGBS.txt",header = T)
dadi_average <- inner_join(dadi_average,CorrectingFactorGBS,by="Species")
colnames(dadi_average)[colnames(dadi_average) == 'Total_length'] <- 'TotalLength'
mu <- 0.59*10^(-8) # From Ossowski et al. 2010 in Arabidopsis thaliana
#dadi_average$K <- dadi_average$theta / (4*dadi_average$TotalLength*mu)
p_snp <- dadi_average$Nsnp/dadi_average$Ntot
dadi_average$meanPi <- (dadi_average$pi0+dadi_average$pi1)/(2*dadi_average$Ntot*dadi_average$Lcontig)
dadi_average$Nref <- dadi_average$theta/(4*mu*dadi_average$Ntot*dadi_average$Lcontig)
dadi_average$Nrefalt <- dadi_average$theta* (-log(1-p_snp)/p_snp)/(4*mu*dadi_average$Ntot*dadi_average$Lcontig)
# Remove pairs of populations not used in the rest of the analysis
dadi_average <- dadi_average[dadi_average$dataset!="CB_Oppdal-Norway",]
dadi_average <- dadi_average[dadi_average$dataset!="CB_Oppdal-Svalbard",]
dadi_average <- dadi_average[dadi_average$dataset!="RP_Svalbard-Norway",]
dadi_average <- dadi_average[dadi_average$dataset!="SAC_Svalbard_Lovisa-Alaska",]
dadi_average <- dadi_average[dadi_average$dataset!="SAC_Svalbard_Lovisa-Oppdal",]
dadi_average <- dadi_average[dadi_average$dataset!="SAC_Svalbard_Lovisa-Svalbard_Pernille",]
dadi_average <- dadi_average[dadi_average$dataset!="SU_Svalbard_Lovisa-Svalbard_Pernille",]
dadi_average <- dadi_average[dadi_average$dataset!="DN_Dovre-Svalbard_Pernille",]

# Create new species*pop index for RP and SAC where two and three pairs are used
dadi_average$Species <- ifelse(dadi_average$dataset=="RP_Alaska-Norway","Ranunculus_pygmaeus_AN",dadi_average$Species)
dadi_average$Species <- ifelse(dadi_average$dataset=="RP_Svalbard-Alaska","Ranunculus_pygmaeus_AS",dadi_average$Species)
dadi_average$Species <- ifelse(dadi_average$dataset=="SAC_Alaska-Oppdal","Silene_acaulis_AN",dadi_average$Species)
dadi_average$Species <- ifelse(dadi_average$dataset=="SAC_Alaska-Svalbard_Pernille","Silene_acaulis_AS",dadi_average$Species)
dadi_average$Species <- ifelse(dadi_average$dataset=="SAC_Oppdal-Svalbard_Pernille","Silene_acaulis_NS",dadi_average$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="A" & HybData$PatPop=="N","Ranunculus_pygmaeus_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="N" & HybData$PatPop=="A","Ranunculus_pygmaeus_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="A" & HybData$PatPop=="S","Ranunculus_pygmaeus_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="S" & HybData$PatPop=="S","Ranunculus_pygmaeus_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="A" & HybData$PatPop=="N","Silene_acaulis_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="N" & HybData$PatPop=="A","Silene_acaulis_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="A" & HybData$PatPop=="S","Silene_acaulis_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="S" & HybData$PatPop=="A","Silene_acaulis_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="N" & HybData$PatPop=="S","Silene_acaulis_NS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="S" & HybData$PatPop=="N","Silene_acaulis_NS",HybData$Species)

HybData_dadi <- merge(HybData,dadi_average,by="Species")
HybData_dadi <- subset(HybData_dadi,!is.na(Fst)&!is.na(PollenFert)&!is.na(MeanSelfRate)) #remove rows with missing data

Summary_HybData_dadi <- HybData_dadi %>% 
  dplyr::group_by(Dataset.y,Species,Fst,model) %>% 
  dplyr::summarise(
    PollenFert=mean(PollenFert, na.rm=T),
    nPollenFert=mean(nPollenFert, na.rm=T),
    MeanParentalFertility=mean(MeanParentalFertility, na.rm=T),
    MeanSelfRate=mean(MeanSelfRate, na.rm=T),
    SeedGerm=mean(SeedGerm,na.rm=T),
    MeanParentalGermination=mean(MeanParentalGermination,na.rm=T),
    STDMeanParentalGermination=mean(STDMeanParentalGermination,na.rm=T),
    meanPi=mean(meanPi),
    Nref=mean(Nref),
    Nrefalt=mean(Nrefalt),
    Mtotal=mean(Mtotal),
    Tsplit=mean(Tsplit),
    N0=mean(N0),
    N1=mean(N1),
    Mean_Dxy=mean(Mean_Dxy)
  )

# Save dadi datasets
write.table(HybData_dadi, file="outputs/synthetized_data/SpArc_HybData_dadi.txt", row.names = FALSE)
write.table(Summary_HybData_dadi, file="outputs/synthetized_data/SpArc_Summary_HybData_dadi.txt", row.names = FALSE)




# ∂a∂i analysis with fixed inbreeding ####

# We extract only the split_mig model

list.results <- list.files(path = "outputs/dadi_results/fixed_inbreeding/", pattern = ".out$")

# List with one fixed F coefficient for the two populations
list_fix1 <- list.results[grep("single",list.results)]
# List with two fixed F coefficients, one per population
list_fix2 <- list.results[grep("single",list.results,invert = T)]

# One F for the two populations ####
dadi_fix1 <- data.frame()
for (input in list_fix1) {
  df <- read.table(paste("outputs/dadi_results/fixed_inbreeding/",input,sep=""),header = T)
  dadi_fix1 <- rbind(dadi_fix1,df[df$model=="split_mig",])
  rm(df)
}
names(dadi_fix1)[which(names(dadi_fix1)=="N1")] <- "N0"
names(dadi_fix1)[which(names(dadi_fix1)=="N2")] <- "N1"
names(dadi_fix1)[which(names(dadi_fix1)=="sN1")] <- "sN0"
names(dadi_fix1)[which(names(dadi_fix1)=="sN2")] <- "sN1"
# Add the variable species to the dadi datasets
for (i in 1:nrow(dadi_fix1)) (dadi_fix1$Species[i] <- cortable[match(strsplit(dadi_fix1$dataset[i],"_")[[1]][1],cortable[,1]),2])
for (i in 1:nrow(dadi_fix1)) (dadi_fix1$Dataset[i] <- cortable[match(strsplit(dadi_fix1$dataset[i],"_")[[1]][1],cortable[,1]),3])

# Add theta and TotalLength (ie, nb of total SNP) from the CorrectingFactorGBS.txt file, and compute K, the correcting factor
CorrectingFactorGBS <- read.table("outputs/gbs_results/CorrectingFactorGBS.txt",header = T)
dadi_fix1 <- inner_join(dadi_fix1,CorrectingFactorGBS,by="Species")
colnames(dadi_fix1)[colnames(dadi_fix1) == 'Total_length'] <- 'TotalLength'
mu <- 0.59*10^(-8) # From Ossowski et al. 2010 in Arabidopsis thaliana
#dadi_free_average$K <- dadi_free_average$theta / (4*dadi_free_average$TotalLength*mu)
p_snp <- dadi_fix1$Nsnp/dadi_fix1$Ntot
dadi_fix1$Nref <- dadi_fix1$theta/(4*mu*dadi_fix1$Ntot*dadi_fix1$Lcontig)
dadi_fix1$Nrefalt <- dadi_fix1$theta* (-log(1-p_snp)/p_snp)/(4*mu*dadi_fix1$Ntot*dadi_fix1$Lcontig)
# Remove pairs of populations not used in the rest of the analysis
dadi_fix1 <- dadi_fix1[dadi_fix1$dataset!="CB_Oppdal-Norway",]
dadi_fix1 <- dadi_fix1[dadi_fix1$dataset!="CB_Oppdal-Svalbard",]
dadi_fix1 <- dadi_fix1[dadi_fix1$dataset!="RP_Svalbard-Norway",]
dadi_fix1 <- dadi_fix1[dadi_fix1$dataset!="SAC_Svalbard_Lovisa-Alaska",]
dadi_fix1 <- dadi_fix1[dadi_fix1$dataset!="SAC_Svalbard_Lovisa-Oppdal",]
dadi_fix1 <- dadi_fix1[dadi_fix1$dataset!="SAC_Svalbard_Lovisa-Svalbard_Pernille",]
dadi_fix1 <- dadi_fix1[dadi_fix1$dataset!="SU_Svalbard_Lovisa-Svalbard_Pernille",]
dadi_fix1 <- dadi_fix1[dadi_fix1$dataset!="DN_Dovre-Svalbard_Pernille",]

# Create new species*pop index for RP and SAC where two and three pairs are used
dadi_fix1$Species <- ifelse(dadi_fix1$dataset=="RP_Alaska-Norway","Ranunculus_pygmaeus_AN",dadi_fix1$Species)
dadi_fix1$Species <- ifelse(dadi_fix1$dataset=="RP_Svalbard-Alaska","Ranunculus_pygmaeus_AS",dadi_fix1$Species)
dadi_fix1$Species <- ifelse(dadi_fix1$dataset=="SAC_Alaska-Oppdal","Silene_acaulis_AN",dadi_fix1$Species)
dadi_fix1$Species <- ifelse(dadi_fix1$dataset=="SAC_Alaska-Svalbard_Pernille","Silene_acaulis_AS",dadi_fix1$Species)
dadi_fix1$Species <- ifelse(dadi_fix1$dataset=="SAC_Oppdal-Svalbard_Pernille","Silene_acaulis_NS",dadi_fix1$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="A" & HybData$PatPop=="N","Ranunculus_pygmaeus_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="N" & HybData$PatPop=="A","Ranunculus_pygmaeus_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="A" & HybData$PatPop=="S","Ranunculus_pygmaeus_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="S" & HybData$PatPop=="S","Ranunculus_pygmaeus_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="A" & HybData$PatPop=="N","Silene_acaulis_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="N" & HybData$PatPop=="A","Silene_acaulis_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="A" & HybData$PatPop=="S","Silene_acaulis_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="S" & HybData$PatPop=="A","Silene_acaulis_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="N" & HybData$PatPop=="S","Silene_acaulis_NS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="S" & HybData$PatPop=="N","Silene_acaulis_NS",HybData$Species)

HybData_dadi_fix1 <- merge(HybData,dadi_fix1,by="Species")
HybData_dadi_fix1 <- subset(HybData_dadi_fix1,!is.na(Fst)&!is.na(PollenFert)&!is.na(MeanSelfRate)) #remove rows with missing data

Summary_HybData_dadi_fix1 <- HybData_dadi_fix1 %>% 
  dplyr::group_by(Dataset.y,Species,Fst,model) %>% 
  dplyr::summarise(
    PollenFert=mean(PollenFert, na.rm=T),
    nPollenFert=mean(nPollenFert, na.rm=T),
    MeanParentalFertility=mean(MeanParentalFertility, na.rm=T),
    MeanSelfRate=mean(MeanSelfRate, na.rm=T),
    SeedGerm=mean(SeedGerm,na.rm=T),
    MeanParentalGermination=mean(MeanParentalGermination,na.rm=T),
    STDMeanParentalGermination=mean(STDMeanParentalGermination,na.rm=T),
    Mean_Dxy=mean(Mean_Dxy),
    Fst=mean(Fst),
    F1=mean(F1),
    F2=mean(F2),
    Self1=mean(2*F1/(1+F1)),
    Self2=mean(2*F2/(1+F2)),
    Self=mean(F1/(1+F1)+F2/(1+F2)),
    Nref=mean(Nref),
    Nrefalt=mean(Nrefalt),
    Mtotal=mean(Manc),  # Only split_mig model so Manc = Mrec = Mtot
    Tsplit=mean(Tsplit),
    N0=mean(N0),
    N1=mean(N1),

  )

# Save dadi datasets
write.table(HybData_dadi_fix1, file="outputs/synthetized_data/SpArc_HybData_dadi_fix1.txt", row.names = FALSE)
write.table(Summary_HybData_dadi_fix1, file="outputs/synthetized_data/SpArc_Summary_HybData_dadi_fix1.txt", row.names = FALSE)


# One F per population ####
dadi_fix2 <- data.frame()
for (input in list_fix2) {
  df <- read.table(paste("outputs/dadi_results/fixed_inbreeding/",input,sep=""),header = T)
  dadi_fix2 <- rbind(dadi_fix2,df[df$model=="split_mig",])
  rm(df)
}
names(dadi_fix2)[which(names(dadi_fix2)=="N1")] <- "N0"
names(dadi_fix2)[which(names(dadi_fix2)=="N2")] <- "N1"
names(dadi_fix2)[which(names(dadi_fix2)=="sN1")] <- "sN0"
names(dadi_fix2)[which(names(dadi_fix2)=="sN2")] <- "sN1"
# Add the variable species to the dadi datasets
for (i in 1:nrow(dadi_fix2)) (dadi_fix2$Species[i] <- cortable[match(strsplit(dadi_fix2$dataset[i],"_")[[1]][1],cortable[,1]),2])
for (i in 1:nrow(dadi_fix2)) (dadi_fix2$Dataset[i] <- cortable[match(strsplit(dadi_fix2$dataset[i],"_")[[1]][1],cortable[,1]),3])

# Add theta and TotalLength (ie, nb of total SNP) from the CorrectingFactorGBS.txt file, and compute K, the correcting factor
CorrectingFactorGBS <- read.table("outputs/gbs_results/CorrectingFactorGBS.txt",header = T)
dadi_fix2 <- inner_join(dadi_fix2,CorrectingFactorGBS,by="Species")
colnames(dadi_fix2)[colnames(dadi_fix2) == 'Total_length'] <- 'TotalLength'
mu <- 0.59*10^(-8) # From Ossowski et al. 2010 in Arabidopsis thaliana
#dadi_free_average$K <- dadi_free_average$theta / (4*dadi_free_average$TotalLength*mu)
p_snp <- dadi_fix2$Nsnp/dadi_fix2$Ntot
dadi_fix2$Nref <- dadi_fix2$theta/(4*mu*dadi_fix2$Ntot*dadi_fix2$Lcontig)
dadi_fix2$Nrefalt <- dadi_fix2$theta* (-log(1-p_snp)/p_snp)/(4*mu*dadi_fix2$Ntot*dadi_fix2$Lcontig)
# Remove pairs of populations not used in the rest of the analysis
dadi_fix2 <- dadi_fix2[dadi_fix2$dataset!="CB_Oppdal-Norway",]
dadi_fix2 <- dadi_fix2[dadi_fix2$dataset!="CB_Oppdal-Svalbard",]
dadi_fix2 <- dadi_fix2[dadi_fix2$dataset!="RP_Svalbard-Norway",]
dadi_fix2 <- dadi_fix2[dadi_fix2$dataset!="SAC_Svalbard_Lovisa-Alaska",]
dadi_fix2 <- dadi_fix2[dadi_fix2$dataset!="SAC_Svalbard_Lovisa-Oppdal",]
dadi_fix2 <- dadi_fix2[dadi_fix2$dataset!="SAC_Svalbard_Lovisa-Svalbard_Pernille",]
dadi_fix2 <- dadi_fix2[dadi_fix2$dataset!="SU_Svalbard_Lovisa-Svalbard_Pernille",]
dadi_fix2 <- dadi_fix2[dadi_fix2$dataset!="DN_Dovre-Svalbard_Pernille",]

# Create new species*pop index for RP and SAC where two and three pairs are used
dadi_fix2$Species <- ifelse(dadi_fix2$dataset=="RP_Alaska-Norway","Ranunculus_pygmaeus_AN",dadi_fix2$Species)
dadi_fix2$Species <- ifelse(dadi_fix2$dataset=="RP_Svalbard-Alaska","Ranunculus_pygmaeus_AS",dadi_fix2$Species)
dadi_fix2$Species <- ifelse(dadi_fix2$dataset=="SAC_Alaska-Oppdal","Silene_acaulis_AN",dadi_fix2$Species)
dadi_fix2$Species <- ifelse(dadi_fix2$dataset=="SAC_Alaska-Svalbard_Pernille","Silene_acaulis_AS",dadi_fix2$Species)
dadi_fix2$Species <- ifelse(dadi_fix2$dataset=="SAC_Oppdal-Svalbard_Pernille","Silene_acaulis_NS",dadi_fix2$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="A" & HybData$PatPop=="N","Ranunculus_pygmaeus_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="N" & HybData$PatPop=="A","Ranunculus_pygmaeus_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="A" & HybData$PatPop=="S","Ranunculus_pygmaeus_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="S" & HybData$PatPop=="S","Ranunculus_pygmaeus_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="A" & HybData$PatPop=="N","Silene_acaulis_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="N" & HybData$PatPop=="A","Silene_acaulis_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="A" & HybData$PatPop=="S","Silene_acaulis_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="S" & HybData$PatPop=="A","Silene_acaulis_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="N" & HybData$PatPop=="S","Silene_acaulis_NS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="S" & HybData$PatPop=="N","Silene_acaulis_NS",HybData$Species)

HybData_dadi_fix2 <- merge(HybData,dadi_fix2,by="Species")
HybData_dadi_fix2 <- subset(HybData_dadi_fix2,!is.na(Fst)&!is.na(PollenFert)&!is.na(MeanSelfRate)) #remove rows with missing data

Summary_HybData_dadi_fix2 <- HybData_dadi_fix2 %>% 
  dplyr::group_by(Dataset.y,Species,Fst,model) %>% 
  dplyr::summarise(
    PollenFert=mean(PollenFert, na.rm=T),
    nPollenFert=mean(nPollenFert, na.rm=T),
    MeanParentalFertility=mean(MeanParentalFertility, na.rm=T),
    MeanSelfRate=mean(MeanSelfRate, na.rm=T),
    SeedGerm=mean(SeedGerm,na.rm=T),
    MeanParentalGermination=mean(MeanParentalGermination,na.rm=T),
    STDMeanParentalGermination=mean(STDMeanParentalGermination,na.rm=T),
    Mean_Dxy=mean(Mean_Dxy),
    Fst=mean(Fst),
    F1=mean(F1),
    F2=mean(F2),
    Self1=mean(2*F1/(1+F1)),
    Self2=mean(2*F2/(1+F2)),
    Self=mean(F1/(1+F1)+F2/(1+F2)),
    Nref=mean(Nref),
    Nrefalt=mean(Nrefalt),
    Mtotal=mean(Manc),  # Only split_mig model so Manc = Mrec = Mtot
    Tsplit=mean(Tsplit),
    N0=mean(N0),
    N1=mean(N1),
    
  )

# Save dadi datasets
write.table(HybData_dadi_fix2, file="outputs/synthetized_data/SpArc_HybData_dadi_fix2.txt", row.names = FALSE)
write.table(Summary_HybData_dadi_fix2, file="outputs/synthetized_data/SpArc_Summary_HybData_dadi_fix2.txt", row.names = FALSE)





# ∂a∂i analysis with free inbreeding ####

# We extract only the split_mig model

list.results <- list.files(path = "outputs/dadi_results/free_inbreeding/", pattern = ".out$")

# List with one free F coefficient for the two populations
list_free1 <- list.results[grep("single",list.results)]
# List with two free F coefficients, one per population
list_free2 <- list.results[grep("single",list.results,invert = T)]

# One F for the two populations ####
dadi_free1 <- data.frame()
for (input in list_free1) {
  df <- read.table(paste("outputs/dadi_results/free_inbreeding/",input,sep=""),header = T)
  dadi_free1 <- rbind(dadi_free1,df[df$model=="split_mig",])
  rm(df)
}
names(dadi_free1)[which(names(dadi_free1)=="N1")] <- "N0"
names(dadi_free1)[which(names(dadi_free1)=="N2")] <- "N1"
names(dadi_free1)[which(names(dadi_free1)=="sN1")] <- "sN0"
names(dadi_free1)[which(names(dadi_free1)=="sN2")] <- "sN1"
# Add the variable species to the dadi datasets
for (i in 1:nrow(dadi_free1)) (dadi_free1$Species[i] <- cortable[match(strsplit(dadi_free1$dataset[i],"_")[[1]][1],cortable[,1]),2])
for (i in 1:nrow(dadi_free1)) (dadi_free1$Dataset[i] <- cortable[match(strsplit(dadi_free1$dataset[i],"_")[[1]][1],cortable[,1]),3])

# Add theta and TotalLength (ie, nb of total SNP) from the CorrectingFactorGBS.txt file, and compute K, the correcting factor
CorrectingFactorGBS <- read.table("outputs/gbs_results/CorrectingFactorGBS.txt",header = T)
dadi_free1 <- inner_join(dadi_free1,CorrectingFactorGBS,by="Species")
colnames(dadi_free1)[colnames(dadi_free1) == 'Total_length'] <- 'TotalLength'
mu <- 0.59*10^(-8) # From Ossowski et al. 2010 in Arabidopsis thaliana
#dadi_free_average$K <- dadi_free_average$theta / (4*dadi_free_average$TotalLength*mu)
p_snp <- dadi_free1$Nsnp/dadi_free1$Ntot
dadi_free1$Nref <- dadi_free1$theta/(4*mu*dadi_free1$Ntot*dadi_free1$Lcontig)
dadi_free1$Nrefalt <- dadi_free1$theta* (-log(1-p_snp)/p_snp)/(4*mu*dadi_free1$Ntot*dadi_free1$Lcontig)
# Remove pairs of populations not used in the rest of the analysis
dadi_free1 <- dadi_free1[dadi_free1$dataset!="CB_Oppdal-Norway",]
dadi_free1 <- dadi_free1[dadi_free1$dataset!="CB_Oppdal-Svalbard",]
dadi_free1 <- dadi_free1[dadi_free1$dataset!="RP_Svalbard-Norway",]
dadi_free1 <- dadi_free1[dadi_free1$dataset!="SAC_Svalbard_Lovisa-Alaska",]
dadi_free1 <- dadi_free1[dadi_free1$dataset!="SAC_Svalbard_Lovisa-Oppdal",]
dadi_free1 <- dadi_free1[dadi_free1$dataset!="SAC_Svalbard_Lovisa-Svalbard_Pernille",]
dadi_free1 <- dadi_free1[dadi_free1$dataset!="SU_Svalbard_Lovisa-Svalbard_Pernille",]
dadi_free1 <- dadi_free1[dadi_free1$dataset!="DN_Dovre-Svalbard_Pernille",]

# Create new species*pop index for RP and SAC where two and three pairs are used
dadi_free1$Species <- ifelse(dadi_free1$dataset=="RP_Alaska-Norway","Ranunculus_pygmaeus_AN",dadi_free1$Species)
dadi_free1$Species <- ifelse(dadi_free1$dataset=="RP_Svalbard-Alaska","Ranunculus_pygmaeus_AS",dadi_free1$Species)
dadi_free1$Species <- ifelse(dadi_free1$dataset=="SAC_Alaska-Oppdal","Silene_acaulis_AN",dadi_free1$Species)
dadi_free1$Species <- ifelse(dadi_free1$dataset=="SAC_Alaska-Svalbard_Pernille","Silene_acaulis_AS",dadi_free1$Species)
dadi_free1$Species <- ifelse(dadi_free1$dataset=="SAC_Oppdal-Svalbard_Pernille","Silene_acaulis_NS",dadi_free1$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="A" & HybData$PatPop=="N","Ranunculus_pygmaeus_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="N" & HybData$PatPop=="A","Ranunculus_pygmaeus_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="A" & HybData$PatPop=="S","Ranunculus_pygmaeus_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="S" & HybData$PatPop=="S","Ranunculus_pygmaeus_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="A" & HybData$PatPop=="N","Silene_acaulis_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="N" & HybData$PatPop=="A","Silene_acaulis_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="A" & HybData$PatPop=="S","Silene_acaulis_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="S" & HybData$PatPop=="A","Silene_acaulis_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="N" & HybData$PatPop=="S","Silene_acaulis_NS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="S" & HybData$PatPop=="N","Silene_acaulis_NS",HybData$Species)

HybData_dadi_free1 <- merge(HybData,dadi_free1,by="Species")
HybData_dadi_free1 <- subset(HybData_dadi_free1,!is.na(Fst)&!is.na(PollenFert)&!is.na(MeanSelfRate)) #remove rows with missing data

Summary_HybData_dadi_free1 <- HybData_dadi_free1 %>% 
  dplyr::group_by(Dataset.y,Species,Fst,model) %>% 
  dplyr::summarise(
    PollenFert=mean(PollenFert, na.rm=T),
    nPollenFert=mean(nPollenFert, na.rm=T),
    MeanParentalFertility=mean(MeanParentalFertility, na.rm=T),
    MeanSelfRate=mean(MeanSelfRate, na.rm=T),
    SeedGerm=mean(SeedGerm,na.rm=T),
    MeanParentalGermination=mean(MeanParentalGermination,na.rm=T),
    STDMeanParentalGermination=mean(STDMeanParentalGermination,na.rm=T),
    Mean_Dxy=mean(Mean_Dxy),
    Fst=mean(Fst),
    F1=mean(F1),
    F2=mean(F2),
    Self1=mean(2*F1/(1+F1)),
    Self2=mean(2*F2/(1+F2)),
    Self=mean(F1/(1+F1)+F2/(1+F2)),
    Nref=mean(Nref),
    Nrefalt=mean(Nrefalt),
    Mtotal=mean(Manc),  # Only split_mig model so Manc = Mrec = Mtot
    Tsplit=mean(Tsplit),
    N0=mean(N0),
    N1=mean(N1),
    
  )

# Save dadi datasets
write.table(HybData_dadi_free1, file="outputs/synthetized_data/SpArc_HybData_dadi_free1.txt", row.names = FALSE)
write.table(Summary_HybData_dadi_free1, file="outputs/synthetized_data/SpArc_Summary_HybData_dadi_free1.txt", row.names = FALSE)


# One F per population ####
dadi_free2 <- data.frame()
for (input in list_free2) {
  df <- read.table(paste("outputs/dadi_results/free_inbreeding/",input,sep=""),header = T)
  dadi_free2 <- rbind(dadi_free2,df[df$model=="split_mig",])
  rm(df)
}
names(dadi_free2)[which(names(dadi_free2)=="N1")] <- "N0"
names(dadi_free2)[which(names(dadi_free2)=="N2")] <- "N1"
names(dadi_free2)[which(names(dadi_free2)=="sN1")] <- "sN0"
names(dadi_free2)[which(names(dadi_free2)=="sN2")] <- "sN1"
# Add the variable species to the dadi datasets
for (i in 1:nrow(dadi_free2)) (dadi_free2$Species[i] <- cortable[match(strsplit(dadi_free2$dataset[i],"_")[[1]][1],cortable[,1]),2])
for (i in 1:nrow(dadi_free2)) (dadi_free2$Dataset[i] <- cortable[match(strsplit(dadi_free2$dataset[i],"_")[[1]][1],cortable[,1]),3])

# Add theta and TotalLength (ie, nb of total SNP) from the CorrectingFactorGBS.txt file, and compute K, the correcting factor
CorrectingFactorGBS <- read.table("outputs/gbs_results/CorrectingFactorGBS.txt",header = T)
dadi_free2 <- inner_join(dadi_free2,CorrectingFactorGBS,by="Species")
colnames(dadi_free2)[colnames(dadi_free2) == 'Total_length'] <- 'TotalLength'
mu <- 0.59*10^(-8) # From Ossowski et al. 2010 in Arabidopsis thaliana
#dadi_free_average$K <- dadi_free_average$theta / (4*dadi_free_average$TotalLength*mu)
p_snp <- dadi_free2$Nsnp/dadi_free2$Ntot
dadi_free2$Nref <- dadi_free2$theta/(4*mu*dadi_free2$Ntot*dadi_free2$Lcontig)
dadi_free2$Nrefalt <- dadi_free2$theta* (-log(1-p_snp)/p_snp)/(4*mu*dadi_free2$Ntot*dadi_free2$Lcontig)
# Remove pairs of populations not used in the rest of the analysis
dadi_free2 <- dadi_free2[dadi_free2$dataset!="CB_Oppdal-Norway",]
dadi_free2 <- dadi_free2[dadi_free2$dataset!="CB_Oppdal-Svalbard",]
dadi_free2 <- dadi_free2[dadi_free2$dataset!="RP_Svalbard-Norway",]
dadi_free2 <- dadi_free2[dadi_free2$dataset!="SAC_Svalbard_Lovisa-Alaska",]
dadi_free2 <- dadi_free2[dadi_free2$dataset!="SAC_Svalbard_Lovisa-Oppdal",]
dadi_free2 <- dadi_free2[dadi_free2$dataset!="SAC_Svalbard_Lovisa-Svalbard_Pernille",]
dadi_free2 <- dadi_free2[dadi_free2$dataset!="SU_Svalbard_Lovisa-Svalbard_Pernille",]
dadi_free2 <- dadi_free2[dadi_free2$dataset!="DN_Dovre-Svalbard_Pernille",]

# Create new species*pop index for RP and SAC where two and three pairs are used
dadi_free2$Species <- ifelse(dadi_free2$dataset=="RP_Alaska-Norway","Ranunculus_pygmaeus_AN",dadi_free2$Species)
dadi_free2$Species <- ifelse(dadi_free2$dataset=="RP_Svalbard-Alaska","Ranunculus_pygmaeus_AS",dadi_free2$Species)
dadi_free2$Species <- ifelse(dadi_free2$dataset=="SAC_Alaska-Oppdal","Silene_acaulis_AN",dadi_free2$Species)
dadi_free2$Species <- ifelse(dadi_free2$dataset=="SAC_Alaska-Svalbard_Pernille","Silene_acaulis_AS",dadi_free2$Species)
dadi_free2$Species <- ifelse(dadi_free2$dataset=="SAC_Oppdal-Svalbard_Pernille","Silene_acaulis_NS",dadi_free2$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="A" & HybData$PatPop=="N","Ranunculus_pygmaeus_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="N" & HybData$PatPop=="A","Ranunculus_pygmaeus_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="A" & HybData$PatPop=="S","Ranunculus_pygmaeus_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Ranunculus_pygmaeus" & HybData$MatPop=="S" & HybData$PatPop=="S","Ranunculus_pygmaeus_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="A" & HybData$PatPop=="N","Silene_acaulis_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="N" & HybData$PatPop=="A","Silene_acaulis_AN",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="A" & HybData$PatPop=="S","Silene_acaulis_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="S" & HybData$PatPop=="A","Silene_acaulis_AS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="N" & HybData$PatPop=="S","Silene_acaulis_NS",HybData$Species)
HybData$Species <- ifelse(HybData$Species=="Silene_acaulis" & HybData$MatPop=="S" & HybData$PatPop=="N","Silene_acaulis_NS",HybData$Species)

HybData_dadi_free2 <- merge(HybData,dadi_free2,by="Species")
HybData_dadi_free2 <- subset(HybData_dadi_free2,!is.na(Fst)&!is.na(PollenFert)&!is.na(MeanSelfRate)) #remove rows with missing data

Summary_HybData_dadi_free2 <- HybData_dadi_free2 %>% 
  dplyr::group_by(Dataset.y,Species,Fst,model) %>% 
  dplyr::summarise(
    PollenFert=mean(PollenFert, na.rm=T),
    nPollenFert=mean(nPollenFert, na.rm=T),
    MeanParentalFertility=mean(MeanParentalFertility, na.rm=T),
    MeanSelfRate=mean(MeanSelfRate, na.rm=T),
    SeedGerm=mean(SeedGerm,na.rm=T),
    MeanParentalGermination=mean(MeanParentalGermination,na.rm=T),
    STDMeanParentalGermination=mean(STDMeanParentalGermination,na.rm=T),
    Mean_Dxy=mean(Mean_Dxy),
    Fst=mean(Fst),
    F1=mean(F1),
    F2=mean(F2),
    Self1=mean(2*F1/(1+F1)),
    Self2=mean(2*F2/(1+F2)),
    Self=mean(F1/(1+F1)+F2/(1+F2)),
    Nref=mean(Nref),
    Nrefalt=mean(Nrefalt),
    Mtotal=mean(Manc),  # Only split_mig model so Manc = Mrec = Mtot
    Tsplit=mean(Tsplit),
    N0=mean(N0),
    N1=mean(N1),
    
  )

# Save dadi datasets
write.table(HybData_dadi_free2, file="outputs/synthetized_data/SpArc_HybData_dadi_free2.txt", row.names = FALSE)
write.table(Summary_HybData_dadi_free2, file="outputs/synthetized_data/SpArc_Summary_HybData_dadi_free2.txt", row.names = FALSE)
