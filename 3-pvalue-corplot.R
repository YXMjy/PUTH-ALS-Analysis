rm(list = ls())
library(readxl)
library(dplyr)
library("ggplot2")

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
  #axis.text.x = NULL,
  axis.text.y = element_text(size = 10, color = 'black', face = 'bold', family = "serif"),
  #legend.position ="none",
  legend.title=element_blank(),
  plot.title = element_text(hjust = 0.5,color = 'black',face = 'bold', family = "serif")
)

#### corplot--Figure3 C-E ###########
load("DMP_Ha_2710.Rdata")
indata<-myDMP_Ha$control_to_case
indata<-vol_data[,c("gene","logP")]
outdata<-read.csv("impgene_out.csv")
outdata<-rename(outdata,gene=symbol)
outdata<-outdata[,c("gene","log_p")]
data<-merge(indata,outdata,by='gene')
data$log_p<-as.numeric(data$log_p)
data<-na.omit(data)
length(unique(data$gene))

data$label=""
data<-data[order(-data$logP),]
top.gene<-head(data$gene[which(data$logP>=0.05)],9)
top10_gene<-c(as.character(top.gene))
data$label[match(top10_gene,data$gene)]<-top10_gene

data$Group="not-significant"
data$Group[which(data$logP>6.78)]="better"
data$Group[which(data$log_p>5.77)]="better"
data$Group[which((data$logP>6.78)&(data$log_p>5.77))]="best"


newdt <- data[!duplicated(data$gene), ]
table(newdt$Group)

plot <- ggplot(data = newdt, aes(logP, log_p)) + geom_point(shape = 20, size = 3.6, aes(col = Group)) +
  scale_color_manual(values = c("#034081","#b0c4de","darkred","grey"))
plot <- plot + axisSetting
plot <- plot + geom_abline(slope = 1, intercept = 0)
plot <- plot + xlab("-log10(p-value)(PUTH cohort)") + 
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 25), expand = c(0, 0)) +
  ylab("-log10(p-value)(NYGC cohort)")+
  geom_text(aes(label = label,col="darkred"), hjust = -0.5)+
  geom_hline(yintercept = 5.77, linetype = "dashed", color = "gray") +  
  geom_vline(xintercept = 6.78, linetype = "dashed", color = "gray")   
plot



##############correlation between expression level(Methy-Gene) ###########################
########### Hannum #############################

## meth data
load("DMP_all.Rdata")
ha1<-myDMP_all$control_to_case
ha1$p.chisq<-pchisq(ha1$t^2/2.134,df=1,lower.tail = FALSE)

ha2<-ha1[,c("t","P.Value","adj.P.Val","control_AVG","case_AVG","deltaBeta","gene","p.chisq")]
ha2$group.pu<-NA
ha2$group.pu[which(ha2$deltaBeta>0)]="hyper"
ha2$group.pu[which(ha2$deltaBeta<0)]="hypo"

## gene data
outdata<-read_excel("impgene_out.xlsx",sheet="Sheet2")
outdata2<-na.omit(outdata)
outdata2$deltaExp<-outdata2$case_mean-outdata2$control_mean
outdata2$group.ny<-NA
outdata2$group.ny[which(outdata2$deltaExp>0)]="high"
outdata2$group.ny[which(outdata2$deltaExp<0)]="low"

## merge data
cor.ha<-merge(ha2,outdata2,by="gene",all.ha2=T)
length(unique(cor.ha$gene))

## chisq.test
t<-table(cor.ha$group.pu,cor.ha$group.ny)
a<-chisq.test(t) #X-squared = 340.79
a$p.value #4.293162e-76

## both 0.05
cor.ha1<-subset(cor.ha,cor.ha$adj.P.Val<0.05&cor.ha$pvalue<0.05)
length(unique(cor.ha1$gene))

## chisq.test
t1<-table(cor.ha1$group.pu,cor.ha1$group.ny)
a1<-chisq.test(t1) #X-squared = 7.3087, df = 1, p-value = 0.006862
a1$p.value

## both adjust 0.05
cor.ha2<-subset(cor.ha,cor.ha$p.chisq<1.65e-7&cor.ha$adjp<0.05)
table(cor.ha2$Group)

t2<-table(cor.ha2$group.pu,cor.ha2$group.ny)
a2<-fisher.test(t2)
a2$p.value
  
## GC-0.1;NYGC-adj0.05
cor.ha3<-subset(cor.ha,cor.ha$p.chisq<3.31e-7&cor.ha$adjp<0.1)
table(cor.ha3$Group) #number:6

cor.ha4<-subset(cor.ha,cor.ha$p.chisq<1.65e-6&cor.ha$adjp<0.5)
table(cor.ha4$Group) #number:9

t4<-table(cor.ha4$group.pu,cor.ha4$group.ny)
a4<-fisher.test(t4)
a4$p.value

########### Horvath #############################

## meth data
load("DMP_all_Ho.Rdata")
ho1<-myDMP_all_Ho$control_to_case
ho1$p.chisq<-pchisq(ho1$t^2/2.619,df=1,lower.tail = FALSE)

ho2<-ho1[,c("t","P.Value","adj.P.Val","control_AVG","case_AVG","deltaBeta","gene","p.chisq")]
ho2$group.pu<-NA
ho2$group.pu[which(ho2$deltaBeta>0)]="hyper"
ho2$group.pu[which(ho2$deltaBeta<0)]="hypo"
remove(myDMP_all_Ho)


## merge data
cor.ho<-merge(ho2,outdata2,by="gene",all.ho2=T)
length(unique(cor.ho$gene))

## chisq.test
t<-table(cor.ho$group.pu,cor.ho$group.ny)
t
a<-chisq.test(t) #X-squared = 488.06
a$p.value #3.768444e-108

cor.ha$Group<-NA
cor.ha$Group[which((cor.ha$deltaBeta>0)&(cor.ha$deltaExp>0))]="++"
cor.ha$Group[which((cor.ha$deltaBeta>0)&(cor.ha$deltaExp<0))]="+-"
cor.ha$Group[which((cor.ha$deltaBeta<0)&(cor.ha$deltaExp>0))]="-+"
cor.ha$Group[which((cor.ha$deltaBeta<0)&(cor.ha$deltaExp<0))]="--"


## both 0.05
cor.ho1<-subset(cor.ho,cor.ho$adj.P.Val<0.05&cor.ho$pvalue<0.05)
length(unique(cor.ho1$gene))

## chisq.test
t1<-table(cor.ho1$group.pu,cor.ho1$group.ny)
t1
a1<-chisq.test(t1) #X-squared = 3.2175, df = 1, p-value = 0.07285482
a1$p.value


## both adjust 0.05
cor.ho2<-subset(cor.ho,cor.ho$p.chisq<1.65e-7&cor.ho$adjp<0.05)
table(cor.ho2$Group)

## chisq.test
t2<-table(cor.ho2$group.pu,cor.ho2$group.ny)
t2
a2<-fisher.test(t2)
a2$p.value

## GC-0.1;NYGC-adj0.05
cor.ho3<-subset(cor.ho,cor.ho$p.chisq<3.31e-7&cor.ho$adjp<0.1)
table(cor.ho3$Group) #number:4

cor.ho4<-subset(cor.ho,cor.ho$p.chisq<1.65e-6&cor.ho$adjp<0.5)
table(cor.ho4$Group) #number:7

t4<-table(cor.ho4$group.pu,cor.ho4$group.ny)
a4<-fisher.test(t4)
a4$p.value


########### Horvath2 #############################

## meth data
load("DMP_all_Ho2.Rdata")
ho2.1<-myDMP_all_Ho2$control_to_case
ho2.1$p.chisq<-pchisq(ho2.1$t^2/1.719,df=1,lower.tail = FALSE)

ho2.2<-ho2.1[,c("t","P.Value","adj.P.Val","control_AVG","case_AVG","deltaBeta","gene","p.chisq")]
ho2.2$group.pu<-NA
ho2.2$group.pu[which(ho2.2$deltaBeta>0)]="hyper"
ho2.2$group.pu[which(ho2.2$deltaBeta<0)]="hypo"
remove(myDMP_all_Ho2)


## merge data
cor.ho2.0<-merge(ho2.2,outdata2,by="gene",all.ho2.2=T)
length(unique(cor.ho2.0$gene))

## chisq.test
t<-table(cor.ho2.0$group.pu,cor.ho2.0$group.ny)
t
a<-chisq.test(t) #X-squared = 241.411
a$p.value #1.93629e-54


## both 0.05
cor.ho2.1<-subset(cor.ho2.0,cor.ho2.0$adj.P.Val<0.05&cor.ho2.0$pvalue<0.05)
length(unique(cor.ho2.1$gene))

## chisq.test
t1<-table(cor.ho2.1$group.pu,cor.ho2.1$group.ny)
t1
a1<-chisq.test(t1) #X-squared = 4.3265, df = 1, p-value = 0.03752302
a1$p.value


## both adjust 0.05
cor.ho2.2<-subset(cor.ho2.0,cor.ho2.0$p.chisq<1.65e-7&cor.ho2.0$adjp<0.05)
table(cor.ho2.2$Group)

## chisq.test
t2<-table(cor.ho2.2$group.pu,cor.ho2.2$group.ny)
t2
a2<-fisher.test(t2)
a2$p.value

## GC-0.1;NYGC-adj0.05
cor.ho2.3<-subset(cor.ho2.0,cor.ho2.0$p.chisq<3.31e-7&cor.ho2.0$adjp<0.1)
table(cor.ho2.3$Group) #number:7

cor.ho2.4<-subset(cor.ho2.0,cor.ho2.0$p.chisq<1.65e-6&cor.ho2.0$adjp<0.5)
table(cor.ho2.4$Group) #number:9

t4<-table(cor.ho2.4$group.pu,cor.ho2.4$group.ny)
a4<-fisher.test(t4)
a4$p.value

########### Levine2 #############################

## meth data
load("DMP_all_L.Rdata")
l1<-myDMP_all_L$control_to_case
l1$p.chisq<-pchisq(l1$t^2/1.879,df=1,lower.tail = FALSE)

l2<-l1[,c("t","P.Value","adj.P.Val","control_AVG","case_AVG","deltaBeta","gene","p.chisq")]
l2$group.pu<-NA
l2$group.pu[which(l2$deltaBeta>0)]="hyper"
l2$group.pu[which(l2$deltaBeta<0)]="hypo"
remove(myDMP_all_L)


## merge data
cor.l<-merge(l2,outdata2,by="gene",all.l2=T)
length(unique(cor.l$gene))

## chisq.test
t<-table(cor.l$group.pu,cor.l$group.ny)
t
a<-chisq.test(t) #X-squared = 195.435
a$p.value #2.070616e-44


## both 0.05
cor.l1<-subset(cor.l,cor.l$adj.P.Val<0.05&cor.l$pvalue<0.05)
length(unique(cor.l1$gene))

## chisq.test
t1<-table(cor.l1$group.pu,cor.l1$group.ny)
t1
a1<-chisq.test(t1) #X-squared = 11.693, df = 1, p-value = 0.0006271471
a1$p.value


## both adjust 0.05
cor.l2<-subset(cor.l,cor.l$p.chisq<1.65e-7&cor.l$adjp<0.05)
table(cor.l$Group)

## chisq.test
t2<-table(cor.l2$group.pu,cor.l2$group.ny)
t2
a2<-fisher.test(t2)
a2$p.value

## GC-0.1;NYGC-adj0.05
cor.l3<-subset(cor.l,cor.l$p.chisq<3.31e-7&cor.l$adjp<0.1)
table(cor.l3$Group) #number:17

cor.l4<-subset(cor.l,cor.l$p.chisq<1.65e-6&cor.l$adjp<0.5)
table(cor.l4$Group) #number:48

t4<-table(cor.l4$group.pu,cor.l4$group.ny)
a4<-fisher.test(t4)
a4$p.value #0.67

######## barplot--Figure2 M-N ################################################################

dat1=read_excel("bar_exp_M&G_ha.xlsx",sheet = "Sheet1") 

#total
barplot1<-ggplot(data=dat1,aes(x=status,y=total,fill=Group))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  #scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  #scale_fill_manual(values = c("#0044BB","#228B22","orange","#B22222"
  #),
  scale_fill_manual(values = c("#87CEFA", "#ffcc22","#66CDAA","#7B68EE"
  ),
  name="")+
  labs(y="Number")+axisSetting
barplot1

#both 0.05/adjust 0.05
barplot3<-ggplot(data=dat1,aes(x=status,y=`adjust p<0.05`,fill=Group))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  #scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  #scale_fill_manual(values = c("#0044BB","#228B22","orange","#B22222"
  #),
  scale_fill_manual(values = c("#87CEFA", "#ffcc22","#66CDAA","#7B68EE"
  ),
  name="")+
  labs(y="Number")+axisSetting
barplot3

