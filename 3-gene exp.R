rm(list = ls())
setwd("D:/Desktop/gene_exp")

#############differential expression analysis ##################################

## load data
# This cohort is from ALS Consortium of New York Genomic Center.
# ID:GSE153960; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153960

dat<-read.table("GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt",header = T)
data<-separate(dat,EnsemblID,into = c("EnsemblID","drop"),sep = "[.]") %>%
  select(-drop)


## ENSG convert Symbol
library("AnnotationDbi")
library("org.Hs.eg.db")

ens<-data$EnsemblID
data$symbol <- mapIds(org.Hs.eg.db,
                    
                    keys=ens,
                    
                    column="SYMBOL",
                    
                    keytype="ENSEMBL",
                    
                    multiVals="first")
save(data,file="gene_exp.Rdata")
write.csv(data,"gene_exp.csv")


### differial expression analysis
load("gene_exp.Rdata")

library(dplyr)
data <- data %>% distinct(symbol, .keep_all = T)
rownames(data)<-data[,1]
data<-data[,-1]
data<-data[,-1660]
newdata<-as.data.frame(t(data))
newdata$ID<-rownames(newdata)

library(readxl)
pd<-read_excel("gene_pd.xlsx",sheet = "Sheet3")
merdata<-merge(pd,newdata,by="ID",all.pd=T)
merdata$group<-as.factor(merdata$group)
save(merdata,file="impgenes_data.Rdata")

#t test
load("impgenes_data.Rdata")
result=data.frame(ID=NA,t=NA,pvalue=NA,adjp=NA)
for(i in 4:35666){
  tt<-t.test(merdata[,i]~merdata$group)
  result[i,1]=colnames(merdata)[i]
  result[i,2]=tt$statistic
  result[i,3]=tt$p.value
}
table(result$pvalue<0.05)

result$adjp<-p.adjust(result$pvalue,
                      method = "bonferroni")
result1<-result[-(1:3),]
write.csv(result1,"impgene_out.csv")



############### key genes expression plot--Figure3 K-O ###############################

library(ggplot2)
library(ggthemes)

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
  #legend.position ="none",
  legend.title=element_text(size = 8, color = 'black', face = 'bold', family = "serif"),
  plot.title = element_text(hjust = 0.5,color = 'black',face = 'bold', family = "serif")
)

###ANKLE2
df1 <- data.frame(
  x = c(0.922,0.859),
  y = c(1126.63,1370.15),
  error.x = c(0.035/sqrt(480),0.067/sqrt(207)),
  error.y = c(574.70/sqrt(809),695.04/sqrt(190)),
  group = rep(c("Case", "Control"), each = 1)  
  
)

ggplot(df1, aes(x = x, y = y)) +
  geom_point(aes(color = group), size = 3) +  
  geom_errorbar(aes(ymin = y - error.y, ymax = y+error.y), width = 0.001) + 
  geom_errorbar(aes(xmin = x - error.x, xmax = x + error.x), width = 0.001) +
  labs( title ="ANKLE2",x = "methylation expression (PUTH-ALS)", y = "gene expression (NYGC-ALS)")+
  axisSetting



###SSH2
df2 <- data.frame(
  x = c(0.920,0.857),
  y = c(1588.49,1782.40),
  error.x = c(0.041/sqrt(480),0.086/sqrt(207)),
  error.y = c(875.13/sqrt(809),1071.66/sqrt(190)),
  group = rep(c("Case", "Control"), each = 1)  
  
)

ggplot(df2, aes(x = x, y = y)) +
  geom_point(aes(color = group), size = 3) +  
  geom_errorbar(aes(ymin = y - error.y, ymax = y+error.y), width = 0.001) + 
  geom_errorbar(aes(xmin = x - error.x, xmax = x + error.x), width = 0.001) +
  labs(title ="SSH2", x = "methylation expression (PUTH-ALS)", y = "gene expression (NYGC-ALS)")+
  axisSetting



###CDC42BPB
df3 <- data.frame(
  x = c(0.895,0.878),
  y = c(4931.64,6717.49),
  error.x = c(0.017/sqrt(480),0.021/sqrt(207)),
  error.y = c(2634.49/sqrt(809),3610.69/sqrt(190)),
  group = rep(c("Case", "Control"), each = 1)  
  
)

ggplot(df3, aes(x = x, y = y)) +
  geom_point(aes(color = group), size = 3) +
  geom_errorbar(aes(ymin = y - error.y, ymax = y+error.y), width = 0.0005) + 
  geom_errorbar(aes(xmin = x - error.x, xmax = x + error.x), width = 0.001) +
  labs( title ="CDC42BPB",x = "methylation expression (PUTH-ALS)", y = "gene expression (NYGC-ALS)")+
  axisSetting


###ELAVL3
df4 <- data.frame(
  x = c(0.637,0.596),
  y = c(4073.67,5946.90),
  error.x = c(0.048/sqrt(480),0.056/sqrt(207)),
  error.y = c(2420.19/sqrt(809),3500.20/sqrt(190)),
  group = rep(c("Case", "Control"), each = 1)  
  
)

ggplot(df4, aes(x = x, y = y)) +
  geom_point(aes(color = group), size = 3) +  
  geom_errorbar(aes(ymin = y - error.y, ymax = y+error.y), width = 0.001) + 
  geom_errorbar(aes(xmin = x - error.x, xmax = x + error.x), width = 0.001) +
  labs( title ="ELAVL3",x = "methylation expression (PUTH-ALS)", y = "gene expression (NYGC-ALS)")+
  axisSetting



###CLEC14A
df5 <- data.frame(
  x = c(0.560,0.485),
  y = c(154.84,245.86),
  error.x = c(0.091/sqrt(480),0.096/sqrt(207)),
  error.y = c(136.42/sqrt(809),199.07/sqrt(190)),
  group = rep(c("Case", "Control"), each = 1)  
  
)

ggplot(df5, aes(x = x, y = y)) +
  geom_point(aes(color = group), size = 3) +  
  geom_errorbar(aes(ymin = y - error.y, ymax = y+error.y), width = 0.001) + 
  geom_errorbar(aes(xmin = x - error.x, xmax = x + error.x), width = 0.001) +
  labs(title ="CLEC14A", x = "methylation expression (PUTH-ALS)", y = "gene expression (NYGC-ALS)")+
  axisSetting

