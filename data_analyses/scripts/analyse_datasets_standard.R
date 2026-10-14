# Script supporting the data analysis presented in Gizaw et al.
# Author: Lucas Marie-Orleach / modified by Sylvain Glémin
# Date: April 2024 / September 2024 / November 2025



# LOADING PACKAGES --------------------------------------------------------

library(readr) # for read.delim
library(lme4) # for GLMM
library(car) # for Anova
library(piecewiseSEM) # for psem
library(dplyr) # for %>%
library(glmmTMB) # Alternative for GLMM
library(Hmisc) # Miscellaneous functions for plotting, data manipulation...
library(ggplot2) # for plotting
library(gridExtra) # for plotting
library(gridGraphics) # for plotting
library(lemon) # for plotting
library(corrplot)
library(ggcorrplot) # For pairwise correlation plots
library(ggrepel) # To avoid overlapping labels on graphics
library(MuMIn) # For model choice
library(DiagrammeR)


# LOADING AND PREPARING ----------------------------------------------------


# 1) Full dataset

HybData <- read_delim("data/experiment_data/SpArc_HybData.txt")
HybData$Species <- as.factor(HybData$Species)
HybData$MS <- ifelse(HybData$MeanSelfRate>0.5,"self","out")
HybData$Region <- ifelse(HybData$Dataset=="Med","Mediterranean","Arctic")
HybData$nSeed <- HybData$GerminatedSeed + HybData$NonGerminatedSeed
Data_lifespan <- read.table("data/experiment_data/Species_lifespan.txt",header = T,sep = "\t")
Data_ms_classification <- read.table("data/MS-species.csv",header = T, sep = ";")
Data_ms_classification <- Data_ms_classification[,c("Species","Classification")]
HybData <- merge(HybData,Data_lifespan,by = "Species")
HybData <- merge(HybData,Data_ms_classification,by = "Species")

Summary_HybData <- HybData %>% dplyr::group_by(Dataset,Species,STDFst) %>% 
  dplyr::summarise(
    NbPollen = sum(!is.na(FertilePollen)),
    NbSeed = sum(!is.na(SeedGerm)),
    Fst=mean(Fst, na.rm=T),
    PollenFert=mean(PollenFert, na.rm=T),
    FertilePollen=mean(FertilePollen,na.rm=T),
    UnFertilePollen=mean(UnFertilePollen,na.rm=T),
    SeedGerm=mean(SeedGerm, na.rm=T),
    nSeed=mean(nSeed, na.rm=T),
    nPollenFert=mean(nPollenFert, na.rm=T),
    STDMeanSelfRate=mean(STDMeanSelfRate, na.rm=T),
    MeanSelfRate=mean(MeanSelfRate, na.rm=T),
    STDMeanParentalFertility=mean(STDMeanParentalFertility, na.rm=T),
    STDMeanParentalGermination=mean(STDMeanParentalGermination, na.rm=T),
    MeanParentalGermination=mean(MeanParentalGermination, na.rm=T),
    MeanParentalFertility=mean(MeanParentalFertility, na.rm=T)
  )
Summary_HybData$Species <- as.factor(Summary_HybData$Species)
Summary_HybData$Region <- ifelse(Summary_HybData$Dataset=="Med","Mediterranean","Arctic")
Summary_HybData <- merge(Summary_HybData,Data_lifespan,by = "Species")
Summary_HybData <- merge(Summary_HybData,Data_ms_classification,by = "Species")

# Getting basic information on the dataset
nrow(HybData) # N individual
nrow(unique(HybData[c("Species","Fst")])) # N crosses
length(unique(HybData$Species)) # N species
# Number of crosses per species
Summary_HybData[,c("Species","NbPollen","NbSeed","Fst")]
sum(Summary_HybData$NbPollen) # With pollen data
summary(HybData[HybData$nSeed>0,]$Dataset) # With seed data


# 2) Dataset with ∂a∂i information

HybData_dadi <- read_delim("outputs/synthetized_data/SpArc_HybData_dadi.txt")
HybData_dadi$RelativeRecentSize <- (HybData_dadi$N0+HybData_dadi$N1)/2
HybData_dadi$div <- as.numeric(scale(HybData_dadi$Mean_Dxy/HybData_dadi$Nref))
HybData_dadi$div <- ifelse(HybData_dadi$div>3,NA,HybData_dadi$div) # Suppression of an outlier
HybData_dadi$AbsoluteRecentSize <- HybData_dadi$Nref*HybData_dadi$RelativeRecentSize
HybData_dadi$TsplitAbs <- 2*HybData_dadi$Nref*HybData_dadi$Tsplit
HybData_dadi$MtotalAbs <- HybData_dadi$Mtotal / (4*HybData_dadi$AbsoluteRecentSize)
# Standard deviation for Pop gene parameters for re-sampling
HybData_dadi$s_AbsoluteRecentSize <- HybData_dadi$Nref*(HybData_dadi$sN0+HybData_dadi$sN1)/2
HybData_dadi$s_TsplitAbs <- 2*HybData_dadi$Nref*HybData_dadi$sTsplit
HybData_dadi$s_Mtotal <- ifelse(!is.na(HybData_dadi$Tinter),
                                   (HybData_dadi$Tsplit*HybData_dadi$sManc+HybData_dadi$Tinter*HybData_dadi$sMrec)/(HybData_dadi$Tsplit+HybData_dadi$Ttotal),
                                   HybData_dadi$sManc
                                )

HybData_dadi$MS <- ifelse(HybData_dadi$MeanSelfRate>0.5,"self","out")
HybData_dadi$Species <- as.factor(HybData_dadi$Species)
HybData_dadi$Dataset <- HybData_dadi$Dataset.x
HybData_dadi$Region <- as.factor(ifelse(HybData_dadi$Dataset=="Med","Mediterranean","Arctic"))
HybData_dadi <- merge(HybData_dadi,Data_ms_classification,by = "Species")

Summary_HybData_dadi <- read_delim("outputs/synthetized_data/SpArc_Summary_HybData_dadi.txt",delim = " ")
Summary_HybData_dadi$RelativeRecentSize <- (Summary_HybData_dadi$N0+Summary_HybData_dadi$N1)/2
Summary_HybData_dadi$div <- as.numeric(scale(Summary_HybData_dadi$Mean_Dxy/Summary_HybData_dadi$Nref))
Summary_HybData_dadi$div <- ifelse(Summary_HybData_dadi$div>3,NA,Summary_HybData_dadi$div) # Suppression of an outlier
Summary_HybData_dadi$AbsoluteRecentSize <- Summary_HybData_dadi$Nref*Summary_HybData_dadi$RelativeRecentSize
Summary_HybData_dadi$TsplitAbs <- 2*Summary_HybData_dadi$Nref*Summary_HybData_dadi$Tsplit
Summary_HybData_dadi$MtotalAbs <- Summary_HybData_dadi$Mtotal / (4*Summary_HybData_dadi$AbsoluteRecentSize)
Summary_HybData_dadi$Species <- as.factor(Summary_HybData_dadi$Species)
Summary_HybData_dadi$Dataset <- Summary_HybData_dadi$Dataset.y
Summary_HybData_dadi$Region <- as.factor(ifelse(Summary_HybData_dadi$Dataset=="Med","Mediterranean","Arctic"))

# Adding generation time for plotting in years instead of generations
Summary_HybData_dadi <- merge(Summary_HybData_dadi,Data_lifespan,by = "Species")
# Adding classification of MS for analyses with only three discrete categories
Summary_HybData_dadi <- merge(Summary_HybData_dadi,Data_ms_classification,by = "Species")

# Getting basic information on the dataset
nrow(HybData_dadi) # N individual
unique(HybData_dadi[c("Species","Fst")]) # N crosses (check names /!\)
length(unique(HybData_dadi$Species)) # N species




# ANALYSING DATA - FULL DATASET -------------------------------------------


# 1) Seed Germination
m1 <- lm(asin(sqrt(SeedGerm)) ~ STDMeanParentalGermination + Region + STDFst*STDMeanSelfRate, data=Summary_HybData)
summary(m1)
Anova(m1, type=3)
sink("results/tab_seed_germ_total.txt")
summary(m1)
Anova(m1, type=3)
sink()

# 2) Pollen Fertility
m2 <- glmer(cbind(FertilePollen, UnFertilePollen) ~ (1|Species_Fst) + Region + STDFst*STDMeanSelfRate + STDMeanParentalFertility,
            data=HybData,
            family=binomial,
            control = glmerControl(optimizer = "bobyqa",autoscale = TRUE)
            )
m2fix <- glm(cbind(FertilePollen, UnFertilePollen) ~  Region + (STDFst * STDMeanSelfRate) + STDMeanParentalGermination, data=HybData, family=binomial)
summary(m2)
Anova(m2, type=3)
plot(m2)
# Random effect is highly significant
summary(m2)$AICtab[[1]]
summary(m2fix)$aic
sink("results/tab_pollen_fert_total.txt")
summary(m2)
Anova(m2, type=3)
sink()

              # Equivalent analyses using three MS categories 

# 3) Seed Germination
m1cat <- lm(asin(sqrt(SeedGerm)) ~  STDMeanParentalGermination + Region + STDFst*Classification, data=Summary_HybData)
summary(m1cat)
Anova(m1cat, type=3)
plot(m1cat)
sink("results/tab_seed_germ_total_catMS.txt")
summary(m1cat)
Anova(m1cat, type=3)
sink()

# 4) Pollen Fertility 3 MS categories
m2cat <- glmer(cbind(FertilePollen, UnFertilePollen) ~ (1|Species_Fst) + STDMeanParentalFertility + Region + STDFst*Classification, data=HybData, family=binomial)
summary(m2cat)
Anova(m2cat, type=3)
plot(m2cat)
sink("results/tab_pollen_fert_total_catMS.txt")
summary(m2cat)
Anova(m2cat, type=3)
sink()


          


# ANALYSINGF DATA - SUBSET WITH DADI DATA - DIVERSITY STATISTICS ----------

        # Only pollen fertility is analysed in the following steps
        # In this part the explanatory variables are the diversity statistics


        # We extract only the average model to avoid doubling the data
HybData_dadi_average <- HybData_dadi[HybData_dadi$model=="average",]
Summary_HybData_dadi_average <- Summary_HybData_dadi[Summary_HybData_dadi$model=="average",]

          # The explanatory variables are scaled
HybData_dadi_average$meanPi <- scale(HybData_dadi_average$meanPi)[,1]
HybData_dadi_average$Fst <- scale(HybData_dadi_average$Fst)[,1]
HybData_dadi_average$Mean_Dxy <- scale(HybData_dadi_average$Mean_Dxy)[,1]
Summary_HybData_dadi_average$meanPi <- scale(Summary_HybData_dadi_average$meanPi)[,1]
Summary_HybData_dadi_average$Fst <- scale(Summary_HybData_dadi_average$Fst)[,1]
Summary_HybData_dadi_average$Mean_Dxy <- scale(Summary_HybData_dadi_average$Mean_Dxy)[,1]

# 1) Correlation between pop gen variables

mat <- as.matrix(Summary_HybData_dadi_average[,c("MeanSelfRate","meanPi","Fst","Mean_Dxy")])
colnames(mat) <- c("Selfing rate","Pi","Fst","Dxy")
corr_mat <- rcorr(mat)
corrplot(corr_mat$r, type = "upper", p.mat = corr_mat$P,sig.level = 0.05, diag = F,insig = "label_sig",
         tl.col = "black", tl.srt = 45,)
          # The variables are highly correlated, especially Fst and Dxy
sink("results/correlation_diversity_statistics.txt")
corr_mat
sink()

# 2) Global analysis 
m3 <- glmer(
  cbind(FertilePollen, UnFertilePollen) ~ (1|Species) + STDMeanParentalFertility + Region + (meanPi + Fst + Mean_Dxy) * STDMeanSelfRate,
  data=HybData_dadi_average,
  family=binomial,
  na.action = "na.fail",
  control = glmerControl(optimizer = "bobyqa",autoscale = TRUE)
)
summary(m3)
Anova(m3, type=3)
plot(m3)
m3fix <- glm(cbind(FertilePollen, UnFertilePollen) ~  STDMeanParentalFertility + Region + (meanPi + Fst + Mean_Dxy) * STDMeanSelfRate,
                data=HybData_dadi_average,
                family=binomial)
# Random effect is highly significant
summary(m3)$AICtab[[1]]
summary(m3fix)$aic

modelselection3 <- dredge(
  m3,
  rank="AICc"
)
bestmodels3 <- get.models(modelselection3,subset = delta < 2)
length(bestmodels3)
summary(bestmodels3[[1]])

sink("results/tab_pollen_fert_dadi_divstat.txt")
summary(m3)
Anova(m3, type=3)
summary(bestmodels3[[1]])
Anova(bestmodels3[[1]], type=3)
sink()

# 3) Same analysis with three MS categories

HybData_dadi_average$Classification2 <- ifelse(HybData_dadi_average$MeanSelfRate<0.5,"out","self")

m3cat <- glmer(
  cbind(FertilePollen, UnFertilePollen) ~ (1|Species) + STDMeanParentalFertility + Region + (meanPi + Fst + Mean_Dxy) * Classification, 
  data=HybData_dadi_average,
  family=binomial,
  na.action = "na.fail",
  control = glmerControl(optimizer = "bobyqa",autoscale = TRUE)
)
summary(m3cat)
Anova(m3cat, type=3)
plot(m3cat)
m3catfix <- glm(cbind(FertilePollen, UnFertilePollen) ~  STDMeanParentalFertility + Region + (meanPi + Fst + Mean_Dxy) * Classification,
             data=HybData_dadi_average,
             family=binomial)
# Random effect is highly significant
summary(m3cat)$AICtab[[1]]
summary(m3catfix)$aic

modelselection3cat <- dredge(
  m3cat,
  rank="AICc"
)
bestmodels3cat <- get.models(modelselection3cat,subset = delta < 1)
length(bestmodels3cat)
summary(bestmodels3cat[[1]])

sink("results/tab_pollen_fert_dadi_divstat_catMS.txt")
summary(m3cat)
Anova(m3cat, type=3)
summary(bestmodels3cat[[1]])
Anova(bestmodels3cat[[1]], type=3)
sink()






############################################## #
# ANALYSING DATA - SUBSET WITH DADI DATA ####
############################################ #


# Only pollen fertility is analysed in the following steps


## Main analysis using the averaged parameters  ####

HybData_dadi_average <- HybData_dadi[HybData_dadi$model=="average",]
Summary_HybData_dadi_average <- Summary_HybData_dadi[Summary_HybData_dadi$model=="average",]

# Correlation between pop gen variables
mat <- as.matrix(Summary_HybData_dadi_average[,c("MeanSelfRate","AbsoluteRecentSize","Tsplit","Mtotal")])
colnames(mat) <- c("Selfing rate","Mean effective size","Divergence time","Migration rate")
corr_mat <- rcorr(mat)

corrplot(corr_mat$r, type = "upper", p.mat = corr_mat$P,sig.level = 0.05, diag = F,insig = "label_sig",
         tl.col = "black", tl.srt = 45,)

sink("results/correlation_genpop_param.txt")
corr_mat
sink()


#Pollen Fertility - Analysis using pop gene parameters as explanatory variables
m4 <- glmer(cbind(FertilePollen, UnFertilePollen) ~ (1|Species) +  STDMeanParentalFertility + Region + (log10(AbsoluteRecentSize) + log10(TsplitAbs) + log10(Mtotal) )* STDMeanSelfRate , 
            data=HybData_dadi_average,
            family=binomial,
            na.action = "na.fail",
            control = glmerControl(optimizer = "bobyqa",autoscale = TRUE)
            )
summary(m4)
Anova(m4, type=3)
plot(m4)

modelselection4 <- dredge(
  m4,
  rank="AIC"
)
bestmodels4 <- get.models(modelselection4,subset = delta < 1)
length(bestmodels4)
summary(bestmodels4[[1]])

sink("results/tab_pollen_fert_dadi.txt")
summary(m4)
Anova(m4, type=3)
summary(bestmodels4[[1]])
Anova(bestmodels4[[1]], type=3)
sink()


# Equivalent analysis using three MS categories but without interactions (matrix rank deficient no convergence)
m4cat <- glmer(cbind(FertilePollen, UnFertilePollen) ~ (1|Species) + STDMeanParentalFertility + Region + (log10(AbsoluteRecentSize) + log10(TsplitAbs) + log10(Mtotal)) * Classification, 
            data=HybData_dadi_average,
            family=binomial,
            na.action = "na.fail",
            control = glmerControl(optimizer = "bobyqa",autoscale = TRUE)
)
summary(m4cat)
Anova(m4cat, type=3)

modelselection4cat <- dredge(
  m4cat,
  rank="AICc"
)
bestmodels4cat <- get.models(modelselection4cat,subset = delta < 1)
length(bestmodels4cat)
summary(bestmodels4cat[[1]])

sink("results/tab_pollen_fert_dadi_catMS.txt")
summary(m4cat)
Anova(m4cat, type=3)
summary(bestmodels4cat[[1]])
Anova(bestmodels4cat[[1]], type=3)
sink()


# Path Analysis

# For path analysis we only use quantitative variables
# The effect of region is not directly tested but thre separated analyses are performed:
# Whole dataset, Arctic species only, Mediterranean species only

# Formating the dataset
HybData_dadi_Bino <- HybData_dadi_average[rep(seq_len(nrow(HybData_dadi_average)), HybData_dadi_average$FertilePollen + HybData_dadi_average$UnFertilePollen), ]
HybData_dadi_Bino$Replicate <- rep(seq_len(nrow(HybData_dadi_average)), HybData_dadi_average$FertilePollen + HybData_dadi_average$UnFertilePollen)
HybData_dadi_Bino$PollenFertBino <- rep( rep(c(0,1),nrow(HybData_dadi_average)), c(rbind(HybData_dadi_average$UnFertilePollen, HybData_dadi_average$FertilePollen)))
HybData_dadi_Bino$Species <- as.factor(HybData_dadi_Bino$Species)
sp <- levels(HybData_dadi_Bino$Species)
HybData_dadi_Bino$weight  <- unlist(sapply(sp,function(x){
  c(1,rep(0,length(which(HybData_dadi_Bino$Species==x))-1))
}))
Summary_HybData_dadi$STDMeanSelfRate = (Summary_HybData_dadi$MeanSelfRate - mean(Summary_HybData_dadi$MeanSelfRate)) / sd(Summary_HybData_dadi$MeanSelfRate)
Summary_HybData_dadi$STDMeanParentalFertility = (Summary_HybData_dadi$MeanParentalFertility - mean(Summary_HybData_dadi$MeanParentalFertility)) / sd(Summary_HybData_dadi$MeanParentalFertility)
Summary_HybData_dadi_average$STDMeanSelfRate = (Summary_HybData_dadi_average$MeanSelfRate - mean(Summary_HybData_dadi_average$MeanSelfRate)) / sd(Summary_HybData_dadi_average$MeanSelfRate)
Summary_HybData_dadi_average$STDMeanParentalFertility = (Summary_HybData_dadi_average$MeanParentalFertility - mean(Summary_HybData_dadi_average$MeanParentalFertility)) / sd(Summary_HybData_dadi_average$MeanParentalFertility)


# Good model to keep
# Check Tests of directed separation: what is suppose to be independent: it's OK, nothing is correlated
# Then predicted path analysis
# To use the same dataset but not inflate the number of parameters for species level analyses:
# weight 1 for the first replicate and 0 for all others: so just one point per species
# glmmTMB is faster than glmer and seems to be more robust


m5 <- psem(
  glmmTMB(
    PollenFertBino ~ log10(Mtotal) +  log10(AbsoluteRecentSize) + log10(TsplitAbs) + STDMeanSelfRate +
      STDMeanParentalFertility + (1|Species/Replicate), data=HybData_dadi_Bino,
    family=binomial(link = "logit")),
  lm(log10(Mtotal) ~ STDMeanSelfRate,weights = weight,data=HybData_dadi_Bino),
  lm(log10(AbsoluteRecentSize) ~ STDMeanSelfRate ,weights = weight,data=HybData_dadi_Bino),
  lm(log10(TsplitAbs) ~ STDMeanSelfRate,weights = weight,data=HybData_dadi_Bino)
)
m5_sum <- summary(m5)

# Same analysis for Mediterranean species only
m5Med <- psem(
  glmmTMB(
    PollenFertBino ~ log10(Mtotal) +  log10(AbsoluteRecentSize) + log10(TsplitAbs) + STDMeanSelfRate +
      STDMeanParentalFertility + (1|Species/Replicate), data=subset(HybData_dadi_Bino,Dataset.x=="Med"),
    family=binomial(link = "logit")),
  lm(log10(Mtotal) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Med")),
  lm(log10(AbsoluteRecentSize) ~ STDMeanSelfRate ,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Med")),
  lm(log10(TsplitAbs) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Med"))
)
m5Med_sum <- summary(m5Med)

# Same analysis for Arctic species only
m5Arc <- psem(
  glmmTMB(
    PollenFertBino ~ log10(Mtotal) +  log10(AbsoluteRecentSize) + log10(TsplitAbs) + STDMeanSelfRate +
      STDMeanParentalFertility + (1|Species/Replicate), data=subset(HybData_dadi_Bino,Dataset.x=="Arc"),
    family=binomial(link = "logit")),
  lm(log10(Mtotal) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Arc")),
  lm(log10(AbsoluteRecentSize) ~ STDMeanSelfRate ,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Arc")),
  lm(log10(TsplitAbs) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Arc"))
)
m5Arc_sum <- summary(m5Arc)

# to retreive paths with coefficients
plot(m5)
plot(m5Med)
plot(m5Arc)

sink("results/tab_pollen_fert_dadi_SEM.txt")
m5_sum
m5Med_sum
m5Arc_sum
sink()



## Analysis with the parameter estimated in the split_mig model ####

HybData_dadi_split_mig <- HybData_dadi[HybData_dadi$model=="split_mig",]
Summary_HybData_dadi_split_mig <- Summary_HybData_dadi[Summary_HybData_dadi$model=="split_mig",]

# Correlation between pop gen variables
mat <- as.matrix(Summary_HybData_dadi_split_mig[,c("MeanSelfRate","AbsoluteRecentSize","Tsplit","Mtotal")])
colnames(mat) <- c("Selfing rate","Mean effective size","Divergence time","Migration rate")
corr_mat <- rcorr(mat)

corrplot(corr_mat$r, type = "upper", p.mat = corr_mat$P,sig.level = 0.05, diag = F,insig = "n",
         tl.col = "black", tl.srt = 45,)



#Pollen Fertility - Analysis using pop gene parameters as explanatory variables
m6 <- glmer(cbind(FertilePollen, UnFertilePollen) ~ (1|Species)  + STDMeanParentalFertility + Region + ( log10(AbsoluteRecentSize) + log10(TsplitAbs) + log10(Mtotal))*STDMeanSelfRate , 
            data=HybData_dadi_split_mig,
            family=binomial,
            na.action = "na.fail",
            control = glmerControl(optimizer = "bobyqa",autoscale = TRUE))
summary(m6)
Anova(m6, type=3)

modelselection6 <- dredge(
  m6,
  rank="AIC"
)
bestmodels6 <- get.models(modelselection6,subset = delta < 2)
length(bestmodels6)
summary(bestmodels6[[1]])


sink("results/tab_pollen_fert_dadi_split_mig.txt")
summary(m6)
Anova(m6, type=3)
summary(bestmodels6[[1]])
Anova(bestmodels6[[1]], type=3)
sink()


# Equivalent analysis using three MS categories
m6cat <- glmer(cbind(FertilePollen, UnFertilePollen) ~ (1|Species) + STDMeanParentalFertility + Region + (log10(AbsoluteRecentSize) + log10(TsplitAbs) + log10(Mtotal)) * Classification,
               data=HybData_dadi_split_mig,
               family=binomial,
               na.action = "na.fail",
               control = glmerControl(optimizer = "bobyqa",autoscale = TRUE)
               )
summary(m6cat)
Anova(m6cat, type=3)

modelselection6cat <- dredge(
  m6cat,
  rank="AIC"
)
bestmodels6cat <- get.models(modelselection6cat,subset = delta < 2)
length(bestmodels6cat)
summary(bestmodels6cat[[1]])


sink("results/tab_pollen_fert_dadi_catMS_split_mig.txt")
summary(m6cat)
Anova(m6cat, type=3)
summary(bestmodels6cat[[1]])
Anova(bestmodels6cat[[1]], type=3)
sink()


# Path Analysis

# For path analysis we only use quantitative variables
# The effect of region is not directly tested but thre separated analyses are performed:
# Whole dataset, Arctic species only, Mediterranean species only

# Formating the dataset
HybData_dadi_Bino <- HybData_dadi_split_mig[rep(seq_len(nrow(HybData_dadi_split_mig)), HybData_dadi_split_mig$FertilePollen + HybData_dadi_split_mig$UnFertilePollen), ]
HybData_dadi_Bino$Replicate <- rep(seq_len(nrow(HybData_dadi_split_mig)), HybData_dadi_split_mig$FertilePollen + HybData_dadi_split_mig$UnFertilePollen)
HybData_dadi_Bino$PollenFertBino <- rep( rep(c(0,1),nrow(HybData_dadi_split_mig)), c(rbind(HybData_dadi_split_mig$UnFertilePollen, HybData_dadi_split_mig$FertilePollen)))
HybData_dadi_Bino$Species <- as.factor(HybData_dadi_Bino$Species)
sp <- levels(HybData_dadi_Bino$Species)
HybData_dadi_Bino$weight  <- unlist(sapply(sp,function(x){
  c(1,rep(0,length(which(HybData_dadi_Bino$Species==x))-1))
}))
Summary_HybData_dadi$STDMeanSelfRate = (Summary_HybData_dadi$MeanSelfRate - mean(Summary_HybData_dadi$MeanSelfRate)) / sd(Summary_HybData_dadi$MeanSelfRate)
Summary_HybData_dadi$STDMeanParentalFertility = (Summary_HybData_dadi$MeanParentalFertility - mean(Summary_HybData_dadi$MeanParentalFertility)) / sd(Summary_HybData_dadi$MeanParentalFertility)
Summary_HybData_dadi_split_mig$STDMeanSelfRate = (Summary_HybData_dadi_split_mig$MeanSelfRate - mean(Summary_HybData_dadi_split_mig$MeanSelfRate)) / sd(Summary_HybData_dadi_split_mig$MeanSelfRate)
Summary_HybData_dadi_split_mig$STDMeanParentalFertility = (Summary_HybData_dadi_split_mig$MeanParentalFertility - mean(Summary_HybData_dadi_split_mig$MeanParentalFertility)) / sd(Summary_HybData_dadi_split_mig$MeanParentalFertility)


# Good model to keep
# Check Tests of directed separation: what is suppose to be independent: it's OK, nothing is correlated
# Then predicted path analysis
# To use the same dataset but not inflate the number of parameters for species level analyses:
# weight 1 for the first replicate and 0 for all others: so just one point per species
# glmmTMB is faster than glmer and seems to be more robust


m7 <- psem(
  glmmTMB(
    PollenFertBino ~ log10(Mtotal) +  log10(AbsoluteRecentSize) + log10(TsplitAbs) + STDMeanSelfRate +
      STDMeanParentalFertility + (1|Species/Replicate), data=HybData_dadi_Bino,
    family=binomial(link = "logit")),
  lm(log10(Mtotal) ~ STDMeanSelfRate,weights = weight,data=HybData_dadi_Bino),
  lm(log10(AbsoluteRecentSize) ~ STDMeanSelfRate ,weights = weight,data=HybData_dadi_Bino),
  lm(log10(TsplitAbs) ~ STDMeanSelfRate,weights = weight,data=HybData_dadi_Bino)
)
m7_sum <- summary(m7)


# Same analysis for Mediterranean species only
m7Med <- psem(
  glmmTMB(
    PollenFertBino ~ log10(Mtotal) +  log10(AbsoluteRecentSize) + log10(TsplitAbs) + STDMeanSelfRate +
      STDMeanParentalFertility + (1|Species/Replicate), data=subset(HybData_dadi_Bino,Dataset.x=="Med"),
    family=binomial(link = "logit")),
  lm(log10(Mtotal) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Med")),
  lm(log10(AbsoluteRecentSize) ~ STDMeanSelfRate ,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Med")),
  lm(log10(TsplitAbs) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Med"))
)
m7Med_sum <- summary(m7Med)

# Same analysis for Arctic species only
m7Arc <- psem(
  glmmTMB(
    PollenFertBino ~ log10(Mtotal) +  log10(AbsoluteRecentSize) + log10(TsplitAbs) + STDMeanSelfRate +
      STDMeanParentalFertility + (1|Species/Replicate), data=subset(HybData_dadi_Bino,Dataset.x=="Arc"),
    family=binomial(link = "logit")),
  lm(log10(Mtotal) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Arc")),
  lm(log10(AbsoluteRecentSize) ~ STDMeanSelfRate ,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Arc")),
  lm(log10(TsplitAbs) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Arc"))
)
m7Arc_sum <- summary(m7Arc)

# to retreive paths with coefficients
plot(m4)
plot(m4Med)
plot(m4Arc)

sink("results/tab_pollen_fert_dadi_SEM_sup_splitmig.txt")
m7_sum
m7Med_sum
m7Arc_sum
sink()




############### #
## GRAPHICS ####
############# #


# FIGURE 1

# 11: seed
fig1a <- ggplot(Summary_HybData, aes(x=Fst, y=SeedGerm/MeanParentalGermination, colour=MeanSelfRate,shape=Region,size=nSeed)) + geom_point() +
  labs(x="Fst", y="Seed germination rate", size="n", colour="Selfing rate") +
  ylim(0,1.3) + xlim(0,1) + ggtitle("A") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme_classic(base_size=12) + 
  theme(aspect.ratio=1)

# 1B: pollen
fig1b <- ggplot(Summary_HybData, aes(x=Fst, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate)) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert,shape=Region)) +
  #stat_summary(fun = mean, geom="text", size=3, aes(label=Species, vjust=1.5)) +
  labs(x="Fst", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + xlim(0,1) + ggtitle("B") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme_classic(base_size=12) +
  theme(aspect.ratio=1)

fig1 <- grid.arrange(fig1a,fig1b,ncol=2)
ggsave("figures/Figure1.pdf",fig1,width = 14,height = 6)

# FIGURE 1 - SAME BUT WITH SPECIES NAMES

# 11: seed
fig1a <- ggplot(Summary_HybData, aes(x=Fst, y=SeedGerm/MeanParentalGermination, colour=MeanSelfRate,shape=Region,size=nSeed)) + geom_point() +
  labs(x="Fst", y="Seed germination rate", colour="Selfing rate") +
  geom_text_repel(aes(label=Short_name,size=150),col="black",fontface="italic",max.overlaps=20) +
  ylim(0,1.3) + xlim(0,1) + ggtitle("A") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme_classic(base_size=12) + 
  theme(aspect.ratio=1)

# 1B: pollen
fig1b <- ggplot(Summary_HybData, aes(x=Fst, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,size=nPollenFert,shape=Region)) + geom_point() +
  geom_text_repel(aes(label=Short_name,size=15),col="black",fontface="italic",max.overlaps=20) +
  labs(x="Fst", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + xlim(0,1) + ggtitle("B") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme_classic(base_size=12) +
  theme(aspect.ratio=1)

fig1 <- grid.arrange(fig1a,fig1b,ncol=2)
ggsave("figures/Figure1_WithSpeciesName.pdf",fig1,width = 14,height = 6)

# FIGURE 2

# Correlation between pop gen variables
mat <- as.matrix(Summary_HybData_dadi_average[,c("MeanSelfRate","AbsoluteRecentSize","Tsplit","Mtotal")])
colnames(mat) <- c("Selfing rate","Mean effective size","Divergence time","Migration rate")
corr_mat <- rcorr(mat)

#corrplot(corr_mat$r, type = "upper", p.mat = corr_mat$P,sig.level = 0.05, diag = F,insig = "label_sig",
#           tl.col = "black", tl.srt = 45,mar = c(1,1,1,1),title = "A                                 ",)
#grid.echo()
#fig2a <- grid.grab()

significance <- ifelse(corr_mat$P<0.05,1,0)

fig2a <- ggcorrplot(corr = corr_mat$r,type = "upper",ggtheme = theme_classic(),method = "square",
          outline.color ="black",lab = T,lab_size = 4,tl.cex = 12,
          colors = c("darkblue","white","darkred"),title = "A")


fig2b <- ggplot(Summary_HybData_dadi_average, aes(x=AbsoluteRecentSize, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
  theme_classic(base_size=17) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
  #stat_summary(fun = mean, geom="text", size=3, aes(label=Species, vjust=1.5)) +
  labs(x="Mean effective size", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + scale_x_log10() + ggtitle("B") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme(aspect.ratio=1)

fig2c <- ggplot(Summary_HybData_dadi_average, aes(x=TsplitAbs, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
  theme_classic(base_size=17) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
  #stat_summary(fun = mean, geom="text", size=3, aes(label=Species, vjust=1.5)) +
  labs(x="Divergent time (in generations)", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + scale_x_log10() +  ggtitle("C") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme(aspect.ratio=1)

fig2d <- ggplot(Summary_HybData_dadi_average, aes(x=Mtotal, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
  theme_classic(base_size=17) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
  #stat_summary(fun = mean, geom="text", size=3, aes(label=Species, vjust=1.5)) +
  labs(x="Migration rate", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + scale_x_log10() + ggtitle("D") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme(aspect.ratio=1)

fig2 <- grid.arrange(fig2a,fig2b,fig2c,fig2d,ncol=2)

ggsave("figures/Figure2.pdf",fig2,width = 14,height = 10)


# Additional figures
figsup1 <- ggplot(Summary_HybData_dadi_average, aes(x=TsplitAbs*Lifespan1, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
  theme_classic(base_size=17) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
  #stat_summary(fun = mean, geom="text", size=3, aes(label=Species, vjust=1.5)) +
  labs(x="Divergent time (in years)", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + scale_x_log10() +  ggtitle("A") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme(aspect.ratio=1)

figsup2 <- ggplot(Summary_HybData_dadi_average, aes(x=TsplitAbs*Lifespan2, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
  theme_classic(base_size=17) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
  #stat_summary(fun = mean, geom="text", size=3, aes(label=Species, vjust=1.5)) +
  labs(x="Divergent time (in years)", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + scale_x_log10() +  ggtitle("B") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme(aspect.ratio=1)

figSup <- grid.arrange(figsup1,figsup2,ncol=2)
ggsave("figures/FigureSup.pdf",figSup,width = 14,height = 8)



# FIGURE 2 - SAME BUT WITH SPECIES NAMES

# Correlation between pop gen variables
mat <- as.matrix(Summary_HybData_dadi_average[,c("MeanSelfRate","AbsoluteRecentSize","Tsplit","Mtotal")])
colnames(mat) <- c("Selfing rate","Mean effective size","Divergence time","Migration rate")
corr_mat <- rcorr(mat)

#corrplot(corr_mat$r, type = "upper", p.mat = corr_mat$P,sig.level = 0.05, diag = F,insig = "label_sig",
#           tl.col = "black", tl.srt = 45,mar = c(1,1,1,1),title = "A                                 ",)
#grid.echo()
#fig2a <- grid.grab()

significance <- ifelse(corr_mat$P<0.05,1,0)

fig2a <- ggcorrplot(corr = corr_mat$r,type = "upper",ggtheme = theme_classic(),method = "square",
                    outline.color = "black",lab = T,lab_size = 4,tl.cex = 12,
                    colors = c("darkblue","white","darkred"),title = "A")

fig2b <- ggplot(Summary_HybData_dadi_average, aes(x=AbsoluteRecentSize, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
  theme_classic(base_size=17) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
  geom_text_repel(aes(label=Short_name),col="black",fontface="italic",max.overlaps=20) +
  labs(x="Mean effective size", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + scale_x_log10() + ggtitle("B") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme(aspect.ratio=1) +
  geom_abline(intercept = -0.31,slope = 0.15)


fig2c <- ggplot(Summary_HybData_dadi_average, aes(x=TsplitAbs, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
  theme_classic(base_size=17) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
  geom_text_repel(aes(label=Short_name),col="black",fontface="italic") +
  labs(x="Divergent time (in generations)", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + scale_x_log10() +  ggtitle("C") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme(aspect.ratio=1)

fig2d <- ggplot(Summary_HybData_dadi_average, aes(x=Mtotal, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
  theme_classic(base_size=17) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
  geom_text_repel(aes(label=Short_name),col="black",fontface="italic") +
  labs(x="Migration rate", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + scale_x_log10() + ggtitle("D") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme(aspect.ratio=1)

fig2 <- grid.arrange(fig2a,fig2b,fig2c,fig2d,ncol=2)

ggsave("figures/Figure2_WithSpeciesName.pdf",fig2,width = 14,height = 10)


# Additional figures
figsup1 <- ggplot(Summary_HybData_dadi, aes(x=TsplitAbs*Lifespan1, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
  theme_classic(base_size=17) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
  geom_text_repel(aes(label=Short_name),col="black",fontface="italic") +
  labs(x="Divergent time (in years)", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + scale_x_log10() +  ggtitle("A") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme(aspect.ratio=1)

figsup2 <- ggplot(Summary_HybData_dadi, aes(x=TsplitAbs*Lifespan2, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
  theme_classic(base_size=17) +
  stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
  geom_text_repel(aes(label=Short_name),col="black",fontface="italic") +
  labs(x="Divergent time (in years)", y="Pollen fertility", size="n", colour="Selfing rate") +
  ylim(0,1.3) + scale_x_log10() +  ggtitle("B") +
  scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
  theme(aspect.ratio=1)

figSup <- grid.arrange(figsup1,figsup2,ncol=2)
ggsave("figures/FigureSup_WithSpeciesName.pdf",figSup,width = 14,height = 8)





