####Figure 1A

library(ggpubr)
axisSetting <- theme(  # remove grid line
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.border = element_rect(color = "NA", fill = "transparent"),
  axis.line.y = element_line(colour = "black"),
  #axis.line.x = element_line(colour = "black"),
  axis.line.x = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=10,color='black',face='bold'),
  axis.title.y = element_text(size=10,color='black',face='bold'),
  axis.text.x = element_text(size=10,color='black',face='bold'),
  axis.text.y = element_text(size=10,color='black',face='bold'),
  #axis.ticks.x = element_blank(),
  #legend.position ="top",
  # legend.position ="none",
  legend.position = c(0.80, 0.90),
  legend.title=element_blank(),
  legend.text=element_text(size = 10,color = 'black',face='bold'),
  legend.background=element_rect(size=1),
  legend.key.size=unit(0.2,'cm'),
  legend.key.width=unit(0.6,'cm'),
)

load("pd0_689.Rdata")
pd0_689$loc<-0
pd0_689$loc[pd0_689$PC1<0]<-1#N
pd0_689$loc[pd0_689$PC1>0]<-2#S
pd0_689$loc <- factor(pd0_689$loc, levels = c("1", "2", "0"), labels = c("North", "South", "Unknow"))
newdata<-pd0_689
newdata$Loc[which(is.na(newdata$Loc)==1)] <- c("Unknow")
newdata$PC1_N <- newdata$PC1
newdata$PC1_N[which(newdata$loc!="North")] <- NA
newdata$PC2_N <- newdata$PC2
newdata$PC2_N[which(newdata$loc!="North")] <- NA
newdata$PC1_S <- newdata$PC1
newdata$PC1_S[which(newdata$loc!="South")] <- NA
newdata$PC2_S <- newdata$PC2
newdata$PC2_S[which(newdata$loc!="South")] <- NA

#plot
plot <- ggplot(newdata, aes(x = PC1, y = PC2, col = factor(loc))) + geom_point(size = 1)
plot <- plot + axisSetting
plot
plot <- plot + scale_color_manual(values = c("red", "blue","gray60"))
plot
plot <- plot + geom_point(aes(x = PC1, y = PC2), size = 1, col = "gray60")
plot
plot <- plot + geom_point(aes(x = PC1_N, y = PC2_N), size = 1, col = "red")
plot
plot <- plot + scale_x_continuous(expand = c(0, 0), name = "PC1", limits = c(min(newdata$PC1) - 0.005, 0.09)) + 
  scale_y_continuous(expand = c(0, 0), name = "PC2", 
                     limits = c(min(newdata$PC2) - 0.005, 0.06))
plot <- plot + theme(axis.line.x = element_line(color = "black")) + 
  theme(axis.text.x = element_text(size = 10, color = 'black', face = 'bold'))
plot <- plot + geom_point(aes(x = PC1_S, y = PC2_S), size = 1, col = "blue")

plot <- plot + geom_vline(xintercept = 0, linetype = "dashed") 
plot

