# Script supporting the data analysis presented in Gizaw et al.
# Author: Lucas Marie-Orleach / modified by Sylvain Glémin
# Date: September 2026

library(ggplot2)

ref <- read.table("outputs/synthetized_data/SpArc_Summary_HybData_dadi.txt",header = T)
refav <- ref[ref$model=="average",]
refsp <- ref[ref$model=="split_mig",]
fix1 <- read.table("outputs/synthetized_data/SpArc_Summary_HybData_dadi_fix1.txt",header = T)
fix2 <- read.table("outputs/synthetized_data/SpArc_Summary_HybData_dadi_fix2.txt",header = T)
free1 <- read.table("outputs/synthetized_data/SpArc_Summary_HybData_dadi_free1.txt",header = T)
free2 <- read.table("outputs/synthetized_data/SpArc_Summary_HybData_dadi_free2.txt",header = T)

# Correlation between models

matT <- matrix(c(refav$Tsplit,refsp$Tsplit,fix1$Tsplit,fix2$Tsplit,free1$Tsplit,free2$Tsplit),20,6)
colnames(matT) <- c("Reference","Fis = 0","One fixed F","Two fixed F","One free F","Two free F")
corr_matT <- rcorr(matT)
gT <- ggcorrplot(corr = corr_matT$r,type = "upper",ggtheme = theme_classic(),method = "square",
                 outline.color = "black",lab = T,lab_size = 4,tl.cex = 12,
                 colors = c("lightblue","white","pink"),title = "Time split")

matM <- matrix(c(refav$Mtotal,refsp$Mtotal,fix1$Mtotal,fix2$Mtotal,free1$Mtotal,free2$Mtotal),20,6)
colnames(matM) <- c("Reference","Fis = 0","One fixed F","Two fixed F","One free F","Two free F")
corr_matM <- rcorr(matM)
gM <- ggcorrplot(corr = corr_matM$r,type = "upper",ggtheme = theme_classic(),method = "square",
                 outline.color = "black",lab = T,lab_size = 4,tl.cex = 12,
                 colors = c("lightblue","white","pink"),title = "Migration rate")

matN0 <- matrix(c(refav$N0,refsp$N0,fix1$N0,fix2$N0,free1$N0,free2$N0),20,6)
colnames(matN0) <- c("Reference","Fis = 0","One fixed F","Two fixed F","One free F","Two free F")
corr_matN0 <- rcorr(matN0)
gN0 <- ggcorrplot(corr = corr_matN0$r,type = "upper",ggtheme = theme_classic(),method = "square",
                 outline.color = "black",lab = T,lab_size = 4,tl.cex = 12,
                 colors = c("lightblue","white","pink"),title = "Effective size - pop 1")

matN1 <- matrix(c(refav$N1,refsp$N1,fix1$N1,fix2$N1,free1$N1,free2$N1),20,6)
colnames(matN1) <- c("Reference","Fis = 0","One fixed F","Two fixed F","One free F","Two free F")
corr_matN1 <- rcorr(matN1)
gN1 <- ggcorrplot(corr = corr_matN1$r,type = "upper",ggtheme = theme_classic(),method = "square",
                  outline.color = "black",lab = T,lab_size = 4,tl.cex = 12,
                  colors = c("lightblue","white","pink"),title = "Effective size - pop 2")


fig <- grid.arrange(gT,gM,gN0,gN1,ncol=2)

ggsave("figures/Figure_sup_model_cor.pdf",fig,width = 14,height = 10)
