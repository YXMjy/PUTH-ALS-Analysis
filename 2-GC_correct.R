####Figure 2/Table S1

library(ggplot2)
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
  axis.text.y = element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
  #legend.position ="none",
  legend.title=element_blank(),
  plot.title = element_text(hjust = 0.5,color = 'black',face = 'bold', family = "serif")
)

######################## Hannum--Figure2 ######################################################
##### load data
load("DMP_all.Rdata")
data<-myDMP_all$control_to_case
data$Probe<-rownames(data)

##### manhattan plot
data$CHR <- trimws(data$CHR) #remove additional level
data$ID<-rownames(data)
head(data)

mh.data<-data[, c("ID","CHR", "MAPINFO","P.Value")]
mh.data$CHR<-as.integer(mh.data$CHR)
head(mh.data)

library(qqman)
manhattan(
  mh.data,
  chr = "CHR",
  bp = "MAPINFO",
  p = "P.Value",
  snp = "ID",
  col = c("#4169E1", "#87CEFA"),
  chrlabs = NULL,
  suggestiveline = -log10(1.65e-07),
  genomewideline = -log10(1.65e-07),
  highlight = NULL,
  annotatePval = 1.49e-16,
  annotateTop = F)


##QQ plot
qq(data$P.Value)

#calculate inflation factor
pmed=median(data$P.Value)
gc=qchisq(pmed,df=1,lower.tail = F)/qchisq(0.5, df=1)
gc

##Frequency of p-value
summary(data$P.Value)
hist(data$P.Value,
     breaks = seq(0, 1, 0.05),
     xlab = "P-Value", 
     ylab = "Frequency", 
     col = "skyblue", 
     border = "skyblue")


#######correction for genomic inflation##################
data$p.chisq<-pchisq(data$t^2/2.134,df=1,lower.tail = FALSE)
head(data)
sum(data$p.chisq<1.65e-07)#51

newdat_ha <- data[data$p.chisq<=1.65e-07, ]# select 51 probe
save(newdat_ha,file="DMP_ha_51.Rdata")
write.csv(newdat_ha,"DMP_ha_51.csv")

##QQ plot
qq(data$p.chisq,col="#4169e1")

#calculate inflation factor
pmed=median(data$p.chisq)
gc=qchisq(pmed,df=1,lower.tail = F)/qchisq(0.5, df=1)
gc

##Frequency of p-value
summary(data$p.chisq)
hist(data$p.chisq,
     breaks = seq(0, 1, 0.05),
     xlab = "P-Value", 
     ylab = "Frequency", 
     col = "skyblue", 
     border = "skyblue")

##### manhattan plot
data$CHR <- trimws(data$CHR) #remove additional level
data$ID<-rownames(data)
head(data)

mh.data<-data[, c("ID","CHR", "MAPINFO","p.chisq")]
mh.data$CHR<-as.integer(mh.data$CHR)
head(mh.data)


manhattan(
  mh.data,
  chr = "CHR",
  bp = "MAPINFO",
  p = "p.chisq",
  snp = "ID",
  col = c("#4169E1", "#87CEFA"),
  chrlabs = NULL,
  suggestiveline = -log10(1.65e-07),
  genomewideline = -log10(1.65e-07),
  highlight = NULL,
  #annotatePval = 7.0e-9,
  annotateTop = F)

######################## Horvath ######################################################
##### load data
load("DMP_all_Ho.Rdata")
data<-myDMP_all_Ho$control_to_case

##### manhattan plot
data$CHR <- trimws(data$CHR) #remove additional level
data$ID<-rownames(data)
head(data)

mh.data<-data[, c("ID","CHR", "MAPINFO","P.Value")]
mh.data$CHR<-as.integer(mh.data$CHR)
head(mh.data)


library(qqman)
manhattan(
  mh.data,
  chr = "CHR",
  bp = "MAPINFO",
  p = "P.Value",
  snp = "ID",
  col = c("#4169E1", "#87CEFA"),
  chrlabs = NULL,
  suggestiveline = -log10(1.65e-07),
  genomewideline = -log10(1.65e-07),
  highlight = NULL,
  annotatePval = 1.49e-16,
  annotateTop = F)


##QQ plot
qq(data$P.Value)
#calculate inflation factor
pmed=median(data$P.Value)
gc=qchisq(pmed,df=1,lower.tail = F)/qchisq(0.5, df=1)
gc

##Frequency of p-value
summary(data$P.Value)
hist(data$P.Value,
     breaks = seq(0, 1, 0.05),
     xlab = "P-Value", 
     ylab = "Frequency", 
     col = "skyblue", 
     border = "skyblue")


#######correction for genomic inflation############################
data$p.chisq<-pchisq(data$t^2/2.619,df=1,lower.tail = FALSE)
head(data)
summary(data$p.chisq)
sum(data$p.chisq<1.65e-07)#12

newdat_ho <- data[data$p.chisq<=1.65e-07, ]# select 12 probe
save(newdat_ho,file="DMP_ho_12.Rdata")
write.csv(newdat_ho,"DMP_ho_12.csv")


##QQ plot
qq(data$p.chisq,col="#ffcc22")
#calculate inflation factor
pmed=median(data$p.chisq)
gc=qchisq(pmed,df=1,lower.tail = F)/qchisq(0.5, df=1)
gc

##Frequency of p-value
summary(data$p.chisq)
hist(data$p.chisq,
     breaks = seq(0, 1, 0.025),
     xlab = "P-Value", 
     ylab = "Frequency", 
     col = "#FFDD55", 
     border = "#FFDD55")

##### manhattan plot
data$CHR <- trimws(data$CHR) #remove additional level
data$ID<-rownames(data)
head(data)

mh.data<-data[, c("ID","CHR", "MAPINFO","p.chisq")]
mh.data$CHR<-as.integer(mh.data$CHR)
head(mh.data)


manhattan(
  mh.data,
  chr = "CHR",
  bp = "MAPINFO",
  p = "p.chisq",
  snp = "ID",
  col = c("#FFCC22", "#FFDD55"),
  chrlabs = NULL,
  suggestiveline = -log10(1.65e-07),
  genomewideline = -log10(1.65e-07),
  highlight = NULL)#,
  #annotatePval = 1.65e-7,
  #annotateTop = F)


######################## Horvath2 ######################################################
##### load data
load("DMP_all_Ho2.Rdata")
data<-myDMP_all_Ho2$control_to_case

##### manhattan plot
data$CHR <- trimws(data$CHR) #remove additional level
data$ID<-rownames(data)
head(data)

mh.data<-data[, c("ID","CHR", "MAPINFO","P.Value")]
mh.data$CHR<-as.integer(mh.data$CHR)
head(mh.data)


library(qqman)
manhattan(
  mh.data,
  chr = "CHR",
  bp = "MAPINFO",
  p = "P.Value",
  snp = "ID",
  col = c("#4169E1", "#87CEFA"),
  chrlabs = NULL,
  suggestiveline = -log10(1.65e-07),
  genomewideline = -log10(1.65e-07),
  highlight = NULL,
  annotatePval = 1.49e-16,
  annotateTop = F)


##QQ plot
qq(data$P.Value)
#calculate inflation factor
pmed=median(data$P.Value)
gc=qchisq(pmed,df=1,lower.tail = F)/qchisq(0.5, df=1)
gc #1.718525

##Frequency of p-value
summary(data$P.Value)
hist(data$P.Value,
     breaks = seq(0, 1, 0.05),
     xlab = "P-Value", 
     ylab = "Frequency", 
     col = "skyblue", 
     border = "skyblue")


####### correction for genomic inflation ##############################
data$p.chisq<-pchisq(data$t^2/1.719,df=1,lower.tail = FALSE)
head(data)
summary(data$p.chisq)
sum(data$p.chisq<1.65e-07)#36

newdat_ho2 <- data[data$p.chisq<=1.65e-07, ]# select 36 probes
save(newdat_ho2,file="DMP_ho2_36.Rdata")
write.csv(newdat_ho2,"DMP_ho2_36.csv")


##QQ plot
qq(data$p.chisq,col="#3cb371")
#calculate inflation factor
pmed=median(data$p.chisq)
gc=qchisq(pmed,df=1,lower.tail = F)/qchisq(0.5, df=1)
gc #1.001023

##Frequency of p-value
summary(data$p.chisq)
hist(data$p.chisq,
     breaks = seq(0, 1, 0.025),
     xlab = "P-Value", 
     ylab = "Frequency", 
     col = "#66CDAA", 
     border = "#66CDAA")

##### manhattan plot
data$CHR <- trimws(data$CHR) #remove additional level
data$ID<-rownames(data)
head(data)

mh.data<-data[, c("ID","CHR", "MAPINFO","p.chisq")]
mh.data$CHR<-as.integer(mh.data$CHR)
head(mh.data)


manhattan(
  mh.data,
  chr = "CHR",
  bp = "MAPINFO",
  p = "p.chisq",
  snp = "ID",
  col = c("#3cb371", "#66cdaa"),
  chrlabs = NULL,
  suggestiveline = -log10(1.65e-07),
  genomewideline = -log10(1.65e-07),
  highlight = NULL#,
  #annotatePval = 1.5e-8,
  #annotateTop = F
  )


######################## Levine ######################################################
##### load data
load("DMP_all_L.Rdata")
data<-myDMP_all_L$control_to_case

##### manhattan plot
data$CHR <- trimws(data$CHR) #remove additional level
data$ID<-rownames(data)
head(data)

mh.data<-data[, c("ID","CHR", "MAPINFO","P.Value")]
mh.data$CHR<-as.integer(mh.data$CHR)
head(mh.data)


library(qqman)
manhattan(
  mh.data,
  chr = "CHR",
  bp = "MAPINFO",
  p = "P.Value",
  snp = "ID",
  col = c("#4169E1", "#87CEFA"),
  chrlabs = NULL,
  suggestiveline = -log10(1.65e-07),
  genomewideline = -log10(1.65e-07),
  highlight = NULL,
  annotatePval = 1.49e-16,
  annotateTop = F)


##QQ plot
qq(data$P.Value)
#calculate inflation factor
pmed=median(data$P.Value)
gc=qchisq(pmed,df=1,lower.tail = F)/qchisq(0.5, df=1)
gc #1.87871

##Frequency of p-value
summary(data$P.Value)
hist(data$P.Value,
     breaks = seq(0, 1, 0.05),
     xlab = "P-Value", 
     ylab = "Frequency", 
     col = "skyblue", 
     border = "skyblue")


#######correction for genomic inflation ####################################
data$p.chisq<-pchisq(data$t^2/1.879,df=1,lower.tail = FALSE)
head(data)
summary(data$p.chisq)
sum(data$p.chisq<1.65e-07)

newdat_l <- data[data$p.chisq<=1.65e-07, ]# select 564 probes
save(newdat_l,file="DMP_l_564.Rdata")
write.csv(newdat_l,"DMP_l_564.csv")


##QQ plot
qq(data$p.chisq,col="#7B68EE")
#calculate inflation factor
pmed=median(data$p.chisq)
gc=qchisq(pmed,df=1,lower.tail = F)/qchisq(0.5, df=1)
gc #1.001023

##Frequency of p-value
summary(data$p.chisq)
hist(data$p.chisq,
     breaks = seq(0, 1, 0.025),
     xlab = "P-Value", 
     ylab = "Frequency", 
     col = "#9370DB", 
     border = "#9370DB")

##### manhattan plot
data$CHR <- trimws(data$CHR) #remove additional level
data$ID<-rownames(data)
head(data)

mh.data<-data[, c("ID","CHR", "MAPINFO","p.chisq")]
mh.data$CHR<-as.integer(mh.data$CHR)
head(mh.data)


manhattan(
  mh.data,
  chr = "CHR",
  bp = "MAPINFO",
  p = "p.chisq",
  snp = "ID",
  col = c("#7B68EE", "#9370DB"),
  chrlabs = NULL,
  suggestiveline = -log10(1.65e-07),
  genomewideline = -log10(1.65e-07),
  highlight = NULL#,
  #annotatePval = 1.15e-10,
  #annotateTop = F
  )




######### impgenes-pattern--Figure3F-G ##################################################################

pd<-read.table("PD_187.txt",sep="\t",header=T)

#####hannum###########
load("tmplm_Ha.Rdata")
head(tmp_lm_Ha)
ha1<-t(tmp_lm_Ha)
ha2<-data.frame(ha1[,c("cg26515084","cg16209303","cg16125874","cg08296288","cg14338936")])
ha2$X<-rownames(ha2)
ha3<-merge(ha2,pd,by="X",all=T)

## plot
library(patchwork)
library(ggplot2)
library(tidyverse)

data<-ha3
#ANKLE2-cg26515084
p1<-ggplot(data=data,aes(x=cg26515084,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p1

#SSH2-cg14338936
p2<-ggplot(data=data,aes(x=cg14338936,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p2

#CDC42BPB-cg16209303
p3<-ggplot(data=data,aes(x=cg16209303,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p3


#CLEC14A-cg16125874
p4<-ggplot(data=data,aes(x=cg16125874,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p4


#ELAVL3-cg08296288

p5<-ggplot(data=data,aes(x=cg08296288,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p5





#####horvath##########
load("tmplm_Ho.Rdata")

head(tmp_lm_Ho)
ho1<-t(tmp_lm_Ho)
ho2<-data.frame(ho1[,c("cg26515084","cg16209303","cg16125874","cg08296288","cg14338936")])
ho2$X<-rownames(ho2)
ho3<-merge(ho2,pd,by="X",all=T)

## plot
data<-ho3

#ANKLE2-cg26515084
p1<-ggplot(data=data,aes(x=cg26515084,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p1

#SSH2-cg14338936
p2<-ggplot(data=data,aes(x=cg14338936,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p2

#CDC42BPB-cg16209303
p3<-ggplot(data=data,aes(x=cg16209303,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p3


#CLEC14A-cg16125874
p4<-ggplot(data=data,aes(x=cg16125874,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p4


#ELAVL3-cg08296288
p5<-ggplot(data=data,aes(x=cg08296288,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p5




#####horvath2##########
load("tmplm_Ho2.Rdata")

head(tmp_lm_Ho2)
ho21<-t(tmp_lm_Ho2)
ho22<-data.frame(ho21[,c("cg26515084","cg16209303","cg16125874","cg08296288","cg14338936")])
ho22$X<-rownames(ho22)
ho23<-merge(ho22,pd,by="X",all=T)

## plot
data<-ho23

#ANKLE2-cg26515084
p1<-ggplot(data=data,aes(x=cg26515084,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p1

#SSH2-cg14338936
p2<-ggplot(data=data,aes(x=cg14338936,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p2

#CDC42BPB-cg16209303
p3<-ggplot(data=data,aes(x=cg16209303,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p3


#CLEC14A-cg16125874
p4<-ggplot(data=data,aes(x=cg16125874,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p4


#ELAVL3-cg08296288
p5<-ggplot(data=data,aes(x=cg08296288,fill=Group))+
  geom_histogram(binwidth = 0.02,color='white',mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p5




#####levine##########
load("tmplm_L.Rdata")

head(tmp_lm_L)
l1<-t(tmp_lm_L)
l2<-data.frame(l1[,c("cg26515084","cg16209303","cg16125874","cg08296288","cg14338936")])
l2$X<-rownames(l2)
l3<-merge(l2,pd,by="X",all=T)

## plot
data<-l3

#ANKLE2-cg26515084
p1<-ggplot(data=data,aes(x=cg26515084,fill=Group))+
  geom_histogram(binwidth = 0.02,color="white",mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p1

#SSH2-cg14338936
p2<-ggplot(data=data,aes(x=cg14338936,fill=Group))+
  geom_histogram(binwidth = 0.02,color="white",mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p2

#CDC42BPB-cg16209303
p3<-ggplot(data=data,aes(x=cg16209303,fill=Group))+
  geom_histogram(binwidth = 0.02,color="white",mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p3


#CLEC14A-cg16125874
p4<-ggplot(data=data,aes(x=cg16125874,fill=Group))+
  geom_histogram(binwidth = 0.02,color="white",mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p4


#ELAVL3-cg08296288
p5<-ggplot(data=data,aes(x=cg08296288,fill=Group))+
  geom_histogram(binwidth = 0.02,color="white",mapping=aes(y=..density..))+
  scale_x_continuous( limits=c(0.3, 1.0))+
  scale_y_continuous( limits=c(0, 20))+
  labs(x=NULL)+
  axisSetting
p5

