
rm(list = ls())
setwd("D:/Desktop/ALS_202306/ALS_analysis_202306/Fig2_DMP")
################################ install.package_load.data #######################

#### step1.1-install.package ####
options()$repos 
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror


# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("minfi",ask = F,update = F)
BiocManager::install("ChAMP",ask = F,update = F)
BiocManager::install("methylationArrayAnalysis",ask = F,update = F) 
BiocManager::install("wateRmelon",ask = F,update = F)
BiocManager::install("clusterProfiler",ask = F,update = F)


#### step1.2-load.data ####
rm(list = ls())
options(stringsAsFactors = F)
require(GEOquery)
require(Biobase)
library("impute")
library(data.table)
info=fread("D:/Students/2021yaoxm/20230215/PD687_13.txt",data.table = F )
b=info
rownames(b)=b[,1]
b=b[,-1]
head(b)

c=fread("D:/Students/2021yaoxm/20230215/data303.txt",data.table = F )
c=c[,-1]
c[1:4,1:4]
rownames(c)=c[,1]
c=c[,-1]
beta=as.matrix(c)
beta=impute.knn(beta)##deal with NA
betaData=beta$data
betaData=betaData+0.00001
c=betaData
c[1:4,1:4]
identical(colnames(c),rownames(b))

library(ChAMP)
myLoad=champ.filter(beta = c,pd = info,filterSNPs = F)
myLoad
save(myLoad,file = 'step1-output.Rdata')




############################# normalization_covariate adjusted ###################

#### step2.1-quality check ####
rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-myLoad.Rdata')
#QC.GUI(beta = myLoad$beta,pheno =myLoad$pd$Group,arraytype = "450k")

##density plot
control<-which(myLoad$pd$Group=="Control")
case<-which(myLoad$pd$Group=="Case")

i=case[1]
plot(density(myLoad$beta[,i],from=0,to=1),ylim=c(0,4.5),type="n",axes = FALSE,xlab="beta")
for (i in case) {
  lines(density(myLoad$beta[,i],from=0,to=1),col=2)
}
axis(1) 
axis(2)
for (i in control) {
  lines(density(myLoad$beta[,i],from=0,to=1),col=3,axes = FALSE,xlab="beta")
}


#### step2.2-normalization #####
load(file = 'step1-myLoad.Rdata')
load(file="pd0_689.Rdata")
load(file="pd0_687.Rdata")
if(F){ }
myNorm<-champ.norm(beta =myLoad$beta, 
                   arraytype = "450k",
                   cores = 5)
dim(myNorm) 
pD=myLoad$pd
save(myNorm,pD,file = 'step2-champ_myNorm.Rdata')


#### step2.3-SVD(check effect of confounding factors) ####
setwd("D:/Desktop/ALS_202306/ALS_analysis_202306/Fig2_DMP")
load(file = 'step2-champ_myNorm.Rdata')

pd0_689$region[pd0_689$PC1<0]<-'1'
pd0_689$region[pd0_689$PC1>0]<-'2'
pd_svd<-pd0_689[,c('Group','Sex','Hannum','region','Sentrix_ID','Sentrix_Position','Bcell',
                          'CD4T','CD8T','Eos','Mono','Neu','NK')]
#save(pd_svd,file='pd_svd.Rdata')
#load(file='pd_svd.Rdata')
pd_svd<-pd_svd[order(rownames(pd_svd))]
data<-t(norm.beta)
newdata<-data[order(rownames(data)),]
newdata<-t(newdata)
identical(rownames(pd_svd),colnames(newdata))

champ.SVD(beta = newdata%>% as.data.frame(),
          pd=pd_svd,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages-0/")
champ.SVD(beta = myNorm%>% as.data.frame(),
          pd=pd_svd,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages-1/")


#### step2.4-batch effect ####
myCombat<-champ.runCombat(beta = myNorm,
                          pd=pd_svd,
                          variablename="Group",
                          batchname = c("Sentrix_ID","Sentrix_Position"),
                          logitTrans=T)
save(myCombat,pd_svd,file = 'step3-combat_ID&P.Rdata')

#### step2.5-cellcounts ####
load(file = 'step3-combat_ID&P.Rdata')
mycell<-champ.refbase(beta = myCombat,
                      arraytype="450K")
save(mycell,file = 'step4-mycell_ID&P.Rdata')


#### step2.6-adjusting age&gender ####
PD<-read.table("PD687_11.txt",sep="\t",header=T)
load(file = 'step4-mycell_ID&P.Rdata')
mycell8=mycell$CorrectedBeta[,-232]
mycell7=mycell8[,-233]
colnames(mycell$CorrectedBeta)
colnames(mycell7)
identical(PD$ID,colnames(mycell7))
save(mycell7,file = 'step4-mycell7.Rdata')

#hannum
load('step4-mycell7.Rdata')
lm_model2<- lm(t(mycell7) ~ as.factor(PD$Sex)+ PD$Hannum+PD$loc,)
tmp_lm_Ha <- t(lm_model2$res)+rowMeans(mycell7)
save(tmp_lm_Ha,file="tmplm_Ha.Rdata")
#horvath
lm_model3<- lm(t(mycell7) ~ as.factor(PD$Sex)+PD$Horvath+PD$loc,)
tmp_lm_Ho <- t(lm_model3$res)+rowMeans(mycell7)
save(tmp_lm_Ho,file="tmplm_Ho.Rdata")
#horbath2
lm_model4<- lm(t(mycell7) ~ as.factor(PD$Sex)+PD$Horvath2+PD$loc,)
tmp_lm_Ho2 <- t(lm_model4$res)+rowMeans(mycell7)
save(tmp_lm_Ho2,file="tmplm_Ho2.Rdata")
#levine
lm_model5<- lm(t(mycell7) ~ as.factor(PD$Sex)+PD$Levine+PD$loc,)
tmp_lm_L <- t(lm_model5$res)+rowMeans(mycell7)
save(tmp_lm_L,file="tmplm_L.Rdata")


#### step2.7-SVD again(evaluate the result of QC) ####
#myCombat%>% as.data.frame()
PD3<-read.table("PD687_11.txt",sep="\t",header=T)
load(file="pd0_687.Rdata")
pd0_687$Group[pd0_687$Group=="Control"]<-0
pd0_687$Group[pd0_687$Group=="Case"]<-1
pd0_687$Sex[pd0_687$Sex=="M"]<-1
pd0_687$Sex[pd0_687$Sex=="F"]<-2
pd0_687$region[pd0_687$PC1<0]<-1
pd0_687$region[pd0_687$PC1>0]<-2

pd_svd2<-pd0_687[,c('Group','Sex','Hannum',"region",'Sentrix_ID','Sentrix_Position',
                    'Bcell','CD4T','CD8T','Eos','Mono','Neu','NK')]
pd_svd2<-pd_svd2[order(rownames(pd_svd2)),]

load(file="tmplm_Ha.Rdata")
data<-t(tmp_lm_Ha)
newdata<-data[order(rownames(data)),]
newdata<-t(newdata)
identical(rownames(pd_svd2),colnames(newdata))


champ.SVD(beta =newdata%>% as.data.frame(),
          pd=pd_svd2,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages-22/")


######################### step3-DMP ##############################################

rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load("tmplm_Ha.Rdata")
PD<-read.table("PD_187.txt",sep="\t",header=T)
group_list=PD$Group
table(group_list)


##Bonferroni-0.05
#hannum
identical(colnames(tmp_lm_Ha),PD$X)
myDMP_Ha<- champ.DMP(beta = tmp_lm_Ha,pheno=group_list,adjPVal = 0.05,adjust.method = "bonferroni",
                     arraytype = "450K")
save(myDMP_Ha,file = 'DMP_Ha_2710.Rdata')
write.csv(myDMP_Ha,"DMP_Ha_2710.csv")

#horvath
load(file = "tmplm_Ho.Rdata")
identical(colnames(tmp_lm_Ho),PD$X)
myDMP_Ho<- champ.DMP(beta = tmp_lm_Ho,pheno=group_list,adjPVal = 0.05,adjust.method = "bonferroni",
                     arraytype = "450K")
save(myDMP_Ho,file = 'DMP_Ho_924.Rdata')
write.csv(myDMP_Ho,"DMP_Ho_924.csv")

#horvath2
load(file = "tmplm_Ho2.Rdata")
identical(colnames(tmp_lm_Ho2),PD$X)
myDMP_Ho2<- champ.DMP(beta = tmp_lm_Ho2,pheno=group_list,adjPVal = 0.05,adjust.method = "bonferroni",
                      arraytype = "450K")
save(myDMP_Ho2,file = 'DMP_Ho2_866.Rdata')
write.csv(myDMP_Ho2,"DMP_Ho2_866.csv")

#levine
load(file="tmplm_L.Rdata")
identical(colnames(tmp_lm_L),PD$X)
myDMP_L<- champ.DMP(beta = tmp_lm_L,pheno=group_list,adjPVal = 0.05,adjust.method = "bonferroni",
                    arraytype = "450K")
save(myDMP_L,file = 'DMP_Le_4055.Rdata')
write.csv(myDMP_L,"DMP_Le_4055.csv")
