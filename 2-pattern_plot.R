####Figure S3A-C

library(readxl)
library(ggplot2)
library(ggthemes)
library("dplyr")
library("scales")

dat=read_excel("DMP_Ha_2710.xlsx",sheet = "pattern") 

axisSetting <- theme(  # remove grid line
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.border = element_blank(),
  axis.line.y = element_line(colour = "black"),
  axis.line.x = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
  axis.title.y = element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
  axis.text.x = element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
  axis.text.y = element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
  legend.position ="none",
  legend.title=element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
  plot.title = element_text(hjust = 0.5,color = 'black',face = 'bold', family = "serif")
)

#### chr_barplot
chr<-table(dat$Chromosomes,dat$subgroup)
chr<-data.frame(chr)
chr$Freq.sum<-NA
chr$Freq.sum[chr$Var2=="hypermethylated"]<-2416
chr$Freq.sum[chr$Var2=="hypomethylated"]<-294
chr<- within(chr,{rate <- round(Freq/Freq.sum,digits=4)*100})
chr<- subset(chr,select=c(Var1,Var2,rate))
write.csv(chr,"D:/Desktop/ALS_202308/Fig_Tab/pattern/pattern_chr.csv")

ggplot(data=chr, aes(x=Var1,y=rate))+
  geom_bar(stat="identity",aes(fill=Var2),width=0.5,position='stack')+
  scale_fill_manual(values=c('#FF8800','#5599FF'))+
  labs(x="Chromosomes",y="Distribution percentage of the 2710 DMPs (%)") +
  coord_flip()+
  axisSetting

#### gene region_barplot 
gene_reg<-table(dat$`Gene region`,dat$subgroup)
gene_reg<-data.frame(gene_reg)
gene_reg$Freq.sum<-NA
gene_reg$Freq.sum[gene_reg$Var2=="hypermethylated"]<-2416
gene_reg$Freq.sum[gene_reg$Var2=="hypomethylated"]<-294
gene_reg<- within(gene_reg,{rate <- round(Freq/Freq.sum,digits=4)*100})
gene_reg<- subset(gene_reg,select=c(Var1,Var2,rate))

ggplot(data=gene_reg, aes(x=Var1,y=rate))+
  geom_bar(stat="identity",aes(fill=Var2),width=0.5,position='stack')+
  scale_fill_manual(values=c('#FF8800','#5599FF'))+
  labs(x=NULL,y="Distribution percentage of the 2710 DMPs (%)") +
  coord_flip()+
  axisSetting

#### CpG island region_barplot 
cpg_reg<-table(dat$`CpG island region`,dat$subgroup)
cpg_reg<-data.frame(cpg_reg)
cpg_reg$Freq.sum<-NA
cpg_reg$Freq.sum[cpg_reg$Var2=="hypermethylated"]<-2416
cpg_reg$Freq.sum[cpg_reg$Var2=="hypomethylated"]<-294
cpg_reg<- within(cpg_reg,{rate <- round(Freq/Freq.sum,digits=4)*100})
cpg_reg<- subset(cpg_reg,select=c(Var1,Var2,rate))
write.csv(cpg_reg,"D:/Desktop/ALS_202308/Fig_Tab/pattern/pattern_cpg.reg.csv")

ggplot(data=cpg_reg, aes(x=Var1,y=rate))+
  geom_bar(stat="identity",aes(fill=Var2),width=0.5,position='stack')+
  scale_fill_manual(values=c('#FF8800','#5599FF'))+
  labs(x=NULL,y="Distribution percentage of the 2710 DMPs (%)") +
  coord_flip()+
  axisSetting
