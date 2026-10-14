# Script supporting the data analysis presented in Gizaw et al.
# Author: Lucas Marie-Orleach / modified by Sylvain Glémin
# Date: April 2024 / September 2024 / November 2025 / September 2026

#################### #
## LOAD PACKAGES ####
################## #
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



####################### #
## LOADING DATASETS ####
##################### #


# Full dataset #

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


################################ #

# Dataset with ∂a∂i information

HybData_dadi <- read_delim("outputs/synthetized_data/SpArc_HybData_dadi_free1.txt")
HybData_dadi$Mtotal <- HybData_dadi$Manc #Manc = Mrec = Mtotal
HybData_dadi$RelativeRecentSize <- (HybData_dadi$N0+HybData_dadi$N1)/2
HybData_dadi$div <- as.numeric(scale(HybData_dadi$Mean_Dxy/HybData_dadi$Nref))
HybData_dadi$div <- ifelse(HybData_dadi$div>3,NA,HybData_dadi$div) # Suppression of an outlier
HybData_dadi$AbsoluteRecentSize <- HybData_dadi$Nref*HybData_dadi$RelativeRecentSize
HybData_dadi$TsplitAbs <- 2*HybData_dadi$Nref*HybData_dadi$Tsplit
HybData_dadi$MtotalAbs <- HybData_dadi$Mtotal / (4*HybData_dadi$AbsoluteRecentSize)
# Standard deviation for Pop gene parameters for re-sampling
HybData_dadi$s_AbsoluteRecentSize <- HybData_dadi$Nref*(HybData_dadi$sN0+HybData_dadi$sN1)/2
HybData_dadi$s_TsplitAbs <- 2*HybData_dadi$Nref*HybData_dadi$sTsplit
HybData_dadi$s_MtotalAbs <- ifelse(!is.na(HybData_dadi$Tinter),
                                   (HybData_dadi$Tsplit*HybData_dadi$sManc+HybData_dadi$Tinter*HybData_dadi$sMrec)/((HybData_dadi$Tsplit+HybData_dadi$Ttotal)*(4*HybData_dadi$AbsoluteRecentSize)),
                                   HybData_dadi$sManc/(4*HybData_dadi$AbsoluteRecentSize)
)
HybData_dadi$MS <- ifelse(HybData_dadi$MeanSelfRate>0.5,"self","out")
HybData_dadi$Species <- as.factor(HybData_dadi$Species)
HybData_dadi$Dataset <- HybData_dadi$Dataset.x
HybData_dadi$Region <- as.factor(ifelse(HybData_dadi$Dataset=="Med","Mediterranean","Arctic"))
HybData_dadi <- merge(HybData_dadi,Data_ms_classification,by = "Species")

Summary_HybData_dadi <- read_delim("outputs/synthetized_data/SpArc_Summary_HybData_dadi_free1.txt")
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





############################################## #
## ANALYSINGF DATA - SUBSET WITH DADI DATA ####
############################################ #

# Only pollen fertility is analysed in the following steps

# Correlation between pop gen variables
mat <- as.matrix(Summary_HybData_dadi[,c("MeanSelfRate","AbsoluteRecentSize","Tsplit","Mtotal")])
colnames(mat) <- c("Selfing rate","Mean effective size","Divergence time","Migration rate")
corr_mat <- rcorr(mat)


#Pollen Fertility - Analysis using pop gene parameters as explanatory variables
m3 <- glmer(cbind(FertilePollen, UnFertilePollen) ~ (1|Species) + STDMeanParentalFertility + Region + (log10(AbsoluteRecentSize) + log10(TsplitAbs) + log10(MtotalAbs)) * STDMeanSelfRate , 
            data=HybData_dadi,
            family=binomial,
            na.action = "na.fail",
            control = glmerControl(optimizer = "bobyqa",autoscale = TRUE)
)
summary(m3)
Anova(m3, type=3)

modelselection3 <- dredge(
  m3,
  rank="AICc"
)
bestmodels3 <- get.models(modelselection3,subset = delta < 1)
length(bestmodels3)
summary(bestmodels3[[1]])

sink("results/tab_pollen_fert_dadi_supfree1.txt")
summary(m3)
Anova(m3, type=3)
summary(bestmodels3[[1]])
Anova(bestmodels3[[1]], type=3)
sink()

# Equivalent analysis using three MS categories 
m3cat <- glmer(cbind(FertilePollen, UnFertilePollen) ~ (1|Species) + STDMeanParentalFertility + Region + (log10(MtotalAbs) + log10(TsplitAbs) + log10(AbsoluteRecentSize)) * Classification, 
            data=HybData_dadi,
            family=binomial,
            na.action = "na.fail",
            control = glmerControl(optimizer = "bobyqa",autoscale = TRUE)
)
summary(m3cat)
Anova(m3cat, type=3)

modelselection3cat <- dredge(
  m3cat,
  rank="AICc"
)
bestmodels3cat <- get.models(modelselection3cat,subset = delta < 1)
length(bestmodels3cat)
summary(bestmodels3cat[[1]])

sink("results/tab_pollen_fert_dadi_catMS_supfree1.txt")
summary(m3cat)
Anova(m3cat, type=3)
summary(bestmodels3cat[[1]])
Anova(bestmodels3cat[[1]], type=3)
sink()



# Path Analysis

# For path analysis we only use quantitative variables
# The effect of region is not directly tested but thre separated analyses are performed:
# Whole dataset, Arctic species only, Mediterranean species only

# Formating the dataset
HybData_dadi_Bino <- HybData_dadi[rep(seq_len(nrow(HybData_dadi)), HybData_dadi$FertilePollen + HybData_dadi$UnFertilePollen), ]
HybData_dadi_Bino$Replicate <- rep(seq_len(nrow(HybData_dadi)), HybData_dadi$FertilePollen + HybData_dadi$UnFertilePollen)
HybData_dadi_Bino$PollenFertBino <- rep( rep(c(0,1),nrow(HybData_dadi)), c(rbind(HybData_dadi$UnFertilePollen, HybData_dadi$FertilePollen)))
HybData_dadi_Bino$Species <- as.factor(HybData_dadi_Bino$Species)
sp <- levels(HybData_dadi_Bino$Species)
HybData_dadi_Bino$weight  <- unlist(sapply(sp,function(x){
  c(1,rep(0,length(which(HybData_dadi_Bino$Species==x))-1))
}))
Summary_HybData_dadi$STDMeanSelfRate = (Summary_HybData_dadi$MeanSelfRate - mean(Summary_HybData_dadi$MeanSelfRate)) / sd(Summary_HybData_dadi$MeanSelfRate)
Summary_HybData_dadi$STDMeanParentalFertility = (Summary_HybData_dadi$MeanParentalFertility - mean(Summary_HybData_dadi$MeanParentalFertility)) / sd(Summary_HybData_dadi$MeanParentalFertility)


# Good model to keep
# Check Tests of directed separation: what is suppose to be independent: it's OK, nothing is correlated
# Then predicted path analysis
# To use the same dataset but not inflate the number of parameters for species level analyses:
# weight 1 for the first replicate and 0 for all others: so just one point per species
# glmmTMB is faster than glmer and seems to be more robust


m4 <- psem(
  glmmTMB(
    PollenFertBino ~ log10(Mtotal) +  log10(AbsoluteRecentSize) + log10(TsplitAbs) + STDMeanSelfRate +
      STDMeanParentalFertility + (1|Species/Replicate), data=HybData_dadi_Bino,
    family=binomial(link = "logit")),
  lm(log10(Mtotal) ~ STDMeanSelfRate,weights = weight,data=HybData_dadi_Bino),
  lm(log10(AbsoluteRecentSize) ~ STDMeanSelfRate ,weights = weight,data=HybData_dadi_Bino),
  lm(log10(TsplitAbs) ~ STDMeanSelfRate,weights = weight,data=HybData_dadi_Bino)
)
m4_sum <- summary(m4)
plot(m4)

# Same analysis for Mediterranean species only
m4Med <- psem(
  glmmTMB(
    PollenFertBino ~ log10(Mtotal) +  log10(AbsoluteRecentSize) + log10(TsplitAbs) + STDMeanSelfRate +
      STDMeanParentalFertility + (1|Species/Replicate), data=subset(HybData_dadi_Bino,Dataset.x=="Med"),
    family=binomial(link = "logit")),
  lm(log10(Mtotal) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Med")),
  lm(log10(AbsoluteRecentSize) ~ STDMeanSelfRate ,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Med")),
  lm(log10(TsplitAbs) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Med"))
)
m4Med_sum <- summary(m4Med)

# Same analysis for Arctic species only
m4Arc <- psem(
  glmmTMB(
    PollenFertBino ~ log10(Mtotal) +  log10(AbsoluteRecentSize) + log10(TsplitAbs) + STDMeanSelfRate +
      STDMeanParentalFertility + (1|Species/Replicate), data=subset(HybData_dadi_Bino,Dataset.x=="Arc"),
    family=binomial(link = "logit")),
  lm(log10(Mtotal) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Arc")),
  lm(log10(AbsoluteRecentSize) ~ STDMeanSelfRate ,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Arc")),
  lm(log10(TsplitAbs) ~ STDMeanSelfRate,weights = weight,data=subset(HybData_dadi_Bino,Dataset.x=="Arc"))
)
m4Arc_sum <- summary(m4Arc)


sink("results/tab_pollen_fert_dadi_SEM_supfree1.txt")
m4_sum
m4Med_sum
m4Arc_sum
sink()


# ############### #
# ## GRAPHICS ####
# ############# #
# 
# 
# # FIGURE 2
# 
# # Correlation between pop gen variables
# mat <- as.matrix(Summary_HybData_dadi[,c("MeanSelfRate","AbsoluteRecentSize","Tsplit","Mtotal")])
# colnames(mat) <- c("Selfing rate","Mean effective size","Divergence time","Migration rate")
# corr_mat <- rcorr(mat)
# 
# #corrplot(corr_mat$r, type = "upper", p.mat = corr_mat$P,sig.level = 0.05, diag = F,insig = "label_sig",
# #           tl.col = "black", tl.srt = 45,mar = c(1,1,1,1),title = "A                                 ",)
# #grid.echo()
# #fig2a <- grid.grab()
# 
# significance <- ifelse(corr_mat$P<0.05,1,0)
# 
# fig2a <- ggcorrplot(corr = corr_mat$r,type = "upper",ggtheme = theme_classic(),method = "square",
#           outline.color ="black",lab = T,lab_size = 4,tl.cex = 12,
#           colors = c("darkblue","white","darkred"),title = "A")
# 
# fig2b <- ggplot(Summary_HybData_dadi, aes(x=AbsoluteRecentSize, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
#   theme_classic(base_size=17) +
#   stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
#   #stat_summary(fun = mean, geom="text", size=3, aes(label=Species, vjust=1.5)) +
#   labs(x="Mean effective size", y="Pollen fertility", size="n", colour="Selfing rate") +
#   ylim(0,1.3) + scale_x_log10() + ggtitle("B") +
#   scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
#   theme(aspect.ratio=1)
# 
# fig2c <- ggplot(Summary_HybData_dadi, aes(x=TsplitAbs, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
#   theme_classic(base_size=17) +
#   stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
#   #stat_summary(fun = mean, geom="text", size=3, aes(label=Species, vjust=1.5)) +
#   labs(x="Divergent time (in generations)", y="Pollen fertility", size="n", colour="Selfing rate") +
#   ylim(0,1.3) + scale_x_log10() +  ggtitle("C") +
#   scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
#   theme(aspect.ratio=1)
# 
# fig2d <- ggplot(Summary_HybData_dadi, aes(x=Mtotal, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
#   theme_classic(base_size=17) +
#   stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
#   #stat_summary(fun = mean, geom="text", size=3, aes(label=Species, vjust=1.5)) +
#   labs(x="Migration rate", y="Pollen fertility", size="n", colour="Selfing rate") +
#   ylim(0,1.3) + scale_x_log10() + ggtitle("D") +
#   scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
#   theme(aspect.ratio=1)
# 
# fig2 <- grid.arrange(fig2a,fig2b,fig2c,fig2d,ncol=2)
# 
# ggsave("figures/Figure2_supfree1.pdf",fig2,width = 14,height = 10)
# 
# 
# 
# # FIGURE 2 - SAME BUT WITH SPECIES NAMES
# 
# # Correlation between pop gen variables
# mat <- as.matrix(Summary_HybData_dadi[,c("MeanSelfRate","AbsoluteRecentSize","Tsplit","Mtotal")])
# colnames(mat) <- c("Selfing rate","Mean effective size","Divergence time","Migration rate")
# corr_mat <- rcorr(mat)
# 
# #corrplot(corr_mat$r, type = "upper", p.mat = corr_mat$P,sig.level = 0.05, diag = F,insig = "label_sig",
# #           tl.col = "black", tl.srt = 45,mar = c(1,1,1,1),title = "A                                 ",)
# #grid.echo()
# #fig2a <- grid.grab()
# 
# significance <- ifelse(corr_mat$P<0.05,1,0)
# 
# fig2a <- ggcorrplot(corr = corr_mat$r,type = "upper",ggtheme = theme_classic(),method = "square",
#                     outline.color = "black",lab = T,lab_size = 4,tl.cex = 12,
#                     colors = c("darkblue","white","darkred"),title = "A")
# 
# fig2b <- ggplot(Summary_HybData_dadi, aes(x=AbsoluteRecentSize, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
#   theme_classic(base_size=17) +
#   stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
#   geom_text_repel(aes(label=Short_name),col="black",fontface="italic",max.overlaps=20) +
#   labs(x="Mean effective size", y="Pollen fertility", size="n", colour="Selfing rate") +
#   ylim(0,1.3) + scale_x_log10() + ggtitle("B") +
#   scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
#   theme(aspect.ratio=1)
# 
# fig2c <- ggplot(Summary_HybData_dadi, aes(x=TsplitAbs, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
#   theme_classic(base_size=17) +
#   stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
#   geom_text_repel(aes(label=Short_name),col="black",fontface="italic") +
#   labs(x="Divergent time (in generations)", y="Pollen fertility", size="n", colour="Selfing rate") +
#   ylim(0,1.3) + scale_x_log10() +  ggtitle("C") +
#   scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
#   theme(aspect.ratio=1)
# 
# fig2d <- ggplot(Summary_HybData_dadi, aes(x=Mtotal, y=PollenFert/MeanParentalFertility, colour=MeanSelfRate,shape=Region)) +
#   theme_classic(base_size=17) +
#   stat_summary(fun = mean, geom="point", aes(size=nPollenFert)) +
#   geom_text_repel(aes(label=Short_name),col="black",fontface="italic") +
#   labs(x="Migration rate", y="Pollen fertility", size="n", colour="Selfing rate") +
#   ylim(0,1.3) + scale_x_log10() + ggtitle("D") +
#   scale_colour_gradient(low="#1965B0", high="#DC050C", limits=c(0,1)) +
#   theme(aspect.ratio=1)
# 
# fig2 <- grid.arrange(fig2a,fig2b,fig2c,fig2d,ncol=2)
# 
# ggsave("figures/Figure2_WithSpeciesName_supfree1.pdf",fig2,width = 14,height = 10)
# 




