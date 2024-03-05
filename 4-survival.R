rm(list = ls())

## load data
load("D:/Desktop/ALS_analysis_202306/Fig3_survival/sur_data.Rdata")

## scale
sur_data_scale<-as.data.frame(scale(sur_data[,(3:2719)]))


#### cox-unimodel ####
library(survival)

newdata=data.frame(ID=NA,HR=NA,L95CI=NA,H95CI=NA,pvalue=NA)
for(i in 1:2717){
  scox1<-coxph(Surv(month)~sur_data_scale[,i],data=sur_data_scale)
  scox2<-summary(scox1)
  newdata[i,1]=colnames(sur_data_scale)[i]
  newdata[i,2]=scox2$conf.int[,"exp(coef)"]
  newdata[i,3]=scox2$conf.int[,"lower .95"]
  newdata[i,4]=scox2$conf.int[,"upper .95"]
  newdata[i,5]=scox2$coefficients[,"Pr(>|z|)"]
}
table(newdata$pvalue<0.05)
write.csv(newdata,"D:/Desktop/ALS_analysis_202306/Fig3_survival/unicox_2710.csv")

#### cox-mulmodel ####

sdmultiCox=coxph(Surv(month) ~ Sex+age+cg06102777+cg09429700+cg27443844+cg15554338+
                           cg09664186+cg13116071+cg19716038+cg15685223+cg16412370+cg01371004, data = sur_data_scale) #这里有个“.”，代表分析td数据中所有变量（列名）
sur_result<-broom::tidy(sdmultiCox,exponentiate=T,conf.int=T)
tdmultiCoxSum=summary(sdmultiCox)
outResult=data.frame()
outResult=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CIH=tdmultiCoxSum$conf.int[,"upper .95"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"])
outResult=cbind(id=row.names(outResult),outResult)
write.csv(outResult,"D:/Desktop/ALS_analysis_202306/Fig3_survival/mulcox_2710.csv")


#### kaplan-meier--Figure4A-D ####
#install.packages("survminer")
library("survminer")
library(foreign)

fit <- survfit(Surv(month) ~age_code, 
               data = sur_data) 

ggsurvplot(fit, data = sur_data,
           conf.int = TRUE, 
           risk.table = TRUE,
           pval = T,
           surv.median.line = "hv",
           legend.labs = c("1", "2") ,
           xlab = "Time(months)")

##### forest_plot--Figure4E-F #####
library(readxl)
library(ggplot2)
dat=read_excel("D:/Desktop/ALS_analysis_202306/Fig3_survival/forest_plot.xlsx",sheet="mul")

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
p <- ggplot(dat, aes(HR, id))

p + geom_point(size=3.6,col="lightpink") +
  
  geom_errorbarh(aes(xmax =H95CIH, xmin = L95CI,col="lightpink"), height = 0.4) +
  
  scale_x_continuous(limits= c(0.1, 2.6), breaks= seq(0, 2.5, 0.5)) +
  
  geom_vline(aes(xintercept = 1)) +
  
  xlab('Odds Ratio ') + ylab(' ')+axisSetting
p