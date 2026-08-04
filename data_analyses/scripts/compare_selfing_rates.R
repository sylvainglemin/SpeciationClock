# Comparison of the selfing rates estimates with different thresholds

library(ggplot2)
library(gridExtra) # for plotting
library(gridGraphics) # for plotting

SelfData <- read.table("data//summary_selfingrate.txt",header=T)
df1 <- SelfData[SelfData$Thresh==0.05,]
df2 <- SelfData[SelfData$Thresh==0.001,]
df3 <- SelfData[SelfData$Thresh==0.0001,]
df4 <- SelfData[SelfData$Thresh==0.000001,]
df5 <- SelfData[SelfData$Thresh==0,]

G1 <- ggplot(data = df5,aes(x = Fis)) + geom_histogram() +
  xlab("Fis") + ylab("Number of populations") + theme_classic(base_size=12) +
  geom_vline(xintercept = 0) + xlim(-1,1)

df <- rbind(
  merge(df1,df2,by=c("Species","Pop")),
  merge(df1,df4,by=c("Species","Pop"))
)
df$Threshold <- as.factor(df$Thresh.y)

G2 <- ggplot(data = df,aes(x = Fis.x, y = Fis.y,col=Threshold)) + 
  geom_point() + geom_abline(slope = 1,intercept = 0) + theme_classic(base_size=12) +
  xlab("Selfing rates in main analyses") + ylab("Alternative threshold") 

fig <- grid.arrange(G1,G2,ncol=2)

ggsave("figures/FigureSup_compare_selfingrates.pdf",fig,width = 10,height = 4)
