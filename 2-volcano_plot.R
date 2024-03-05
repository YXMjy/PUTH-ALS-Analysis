####################### volcano plot ##########################################################
rm(list = ls())

#install.package
#install.packages("ggpubr")
#install.packages("ggthemes")
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggrepel)

#
axisSetting <- theme(  # remove grid line
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.border = element_blank(),
  axis.line.y = element_line(colour = "black"),
  axis.line.x = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size = 12, color = 'black', face = 'bold', family = "serif"),
  axis.title.y = element_text(size = 12, color = 'black', face = 'bold', family = "serif"),
  axis.text.x = element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
  axis.text.y = element_text(size = 10, color = 'black', face = 'bold', family = "serif"),
  legend.position ="none",
  legend.title=element_blank(),
  plot.title = element_text(hjust = 0.5,color = 'black',face = 'bold', family = "serif")
)

# load.data
load("D:/Desktop/ALS_202306/ALS_analysis_202306/Fig2_DMP/DMP_all.Rdata")
vol_data<-myDMP_all$control_to_case
vol_data$p.chisq<-pchisq(vol_data$t^2/2.134,df=1,lower.tail = FALSE)
vol_data$adj.p.chisq<-vol_data$p.chisq*302052
head(vol_data)
class(vol_data)


# group
vol_data$Group="not-significant"
vol_data$Group[which((vol_data$adj.p.chisq<0.05)&(vol_data$deltaBeta>0))]="hypermethylation"
vol_data$Group[which((vol_data$adj.p.chisq<0.05)&(vol_data$deltaBeta<0))]="hypomethylation"
table(vol_data$Group)

# label key genes
vol_data$label=""
vol_data$label[vol_data$gene=="ANKLE2"&vol_data$p.chisq==4.593269e-28]="ANKLE2"
vol_data$label[vol_data$gene=="SSH2"&vol_data$p.chisq==2.414414e-19]="SSH2"
vol_data$label[vol_data$gene=="CDC42BPB"&vol_data$p.chisq==3.173190e-14]="CDC42BPB"
vol_data$label[vol_data$gene=="ELAVL3"&vol_data$p.chisq==2.547125e-11]="ELAVL3"
vol_data$label[vol_data$gene=="CLEC14A"&vol_data$p.chisq==3.476843e-11]="CLEC14A"
vol_data$label[vol_data$gene=="USP53"&vol_data$p.chisq==3.700688e-09]="USP53"
vol_data$label[vol_data$gene=="KDM5A"&vol_data$p.chisq==5.220846e-09]="KDM5A"

vol_data$logP<-(-log10(vol_data$p.chisq))
ggscatter(vol_data,x="deltaBeta",y="logP",
          ggtheme =  theme_bw())

p<-ggscatter(vol_data,x="deltaBeta",y="logP",
          color = "Group",
          palette = c("#cc0000","#2f5688","#BBBBBB"),
          size = 2,
          #label = vol_data$label,
          #font.label = 12,
          #repel = T,
          #label.rectangle = T,
          #nudge_x = 0.1, nudge_y = 0.1,
          xlab="delta beta",
          ylab="-log10(adjust p-value)")+
  theme_base()+
  geom_hline(yintercept=(-log10(1.65e-7)),linetype="dashed")+
  geom_vline(xintercept = 0,linetype="dashed")+
  axisSetting
p
ggsave("D:/Desktop/ALS_202306/ALS_analysis_202306/Fig2_DMP/volcano_plot-new.tiff", plot = p, width = 4.8, height = 6, units = "in", dpi = 300)
