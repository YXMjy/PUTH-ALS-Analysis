
####load packages###############################################
library(r2redux)
library(dplyr)
#install.packages("fmsb")
library(fmsb)
#install.packages('tidyverse')
library(tidyverse)
library(readxl)
library(ggplot2)
#install.packages('ggthemes')
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

############ PRS_liner ###################################################
pp<-read.csv("pp.csv")

#all data
ppp<-pp[30:39]
result<-data.frame(p=NA,R2=NA,R2l=NA,R2O=NA)
for (i in 1:9) {
  lm<-lm(outcome~ppp[[i]],data=ppp)
  lm_summary <- summary(lm)
  pvalue<- lm_summary$coef[2,"Pr(>|t|)"]
  rsq <- lm_summary$r.squared
  se <- sqrt((1 - rsq) / (672 - 2))
  cc <- cc_trf(rsq, se, 1/100000, 474/672)
  R2l = cc$R2l
  R2O = cc$R2O
  result[i,1] <- pvalue
  result[i,2] <- rsq
  result[i,3] <- R2l
  result[i,4] <- R2O
}
write.csv(result,"R2_prs_all.csv")

#North_liner
North<-subset(pp,PC1<0,select = s1:outcome)
result.n<-data.frame(p=NA,R2=NA,R2l=NA,R2O=NA)
for (i in 1:9) {
  lm<-lm(outcome~North[[i]],data=North)
  lm_summary <- summary(lm)
  pvalue<- lm_summary$coef[2,"Pr(>|t|)"]
  rsq <- lm_summary$r.squared
  se <- sqrt((1 - rsq) / (467 - 2))
  cc <- cc_trf(rsq, se, 1/100000, 354/467)
  R2l = cc$R2l
  R2O = cc$R2O
  result.n[i,1] <- pvalue
  result.n[i,2] <- rsq
  result.n[i,3] <- R2l
  result.n[i,4] <- R2O
}
write.csv(result.n,"R2_prs_n.csv")

#South_liner
South<-subset(pp,PC1>0,select = s1:outcome)
result.s<-data.frame(p=NA,R2=NA,R2l=NA,R2O=NA)
for (i in 1:9) {
  lm<-lm(outcome~South[[i]],data=South)
  lm_summary <- summary(lm)
  pvalue<- lm_summary$coef[2,"Pr(>|t|)"]
  rsq <- lm_summary$r.squared
  se <- sqrt((1 - rsq) / (205- 2))
  cc <- cc_trf(rsq, se, 1/100000, 120/205)
  R2l = cc$R2l
  R2O = cc$R2O
  result.s[i,1] <- pvalue
  result.s[i,2] <- rsq
  result.s[i,3] <- R2l
  result.s[i,4] <- R2O
}
write.csv(result.s,"R2_prs_s.csv")


#### MS&PRS liner model ############################################################

## all data
pm<-read.csv("lasso_pms_data.csv")
pm<-read.csv("xgb_pms_data.csv")
pms<-pm[c("ID","Group","s5","outcome","y_test","pla.value")]
pms<-rename(pms,y_pred=pla.value)

#MS
lm2<- lm(outcome ~ y_pred, data = pms)
summary(lm2)
rsq2<-summary(lm2)$r.squared 
se2 <- sqrt((1 - rsq2) / (672 - 2))
cc_trf(rsq2, se2, 1/100000, 474/672)

#MS+PRS
lm3<- lm(outcome ~ s5+y_pred, data = pms)
summary(lm3)
rsq3<-summary(lm3)$r.squared 
se3 <- sqrt((1 - rsq3) / (672 - 3))
cc_trf(rsq3, se3, 1/100000, 474/672)

#MS+PRS+MS*PRS
lm4<- lm(outcome ~ s5+y_pred+s5*y_pred, data = pms)
summary(lm4)
rsq4<-summary(lm4)$r.squared 
se4 <- sqrt((1 - rsq4) / (672 - 4))
cc_trf(rsq4, se4, 1/100000, 474/672)



#North##################
North<-subset(pm,PC1<0,select = c("ID","Group","s5","outcome","y_test","pla.value"))
North<-rename(North,y_pred=pla.value)
table(North$Group)

#MS
lm2<- lm(outcome ~ y_pred, data = North)
summary(lm2)
rsq2<-summary(lm2)$r.squared 
se2 <- sqrt((1 - rsq2) / (467 - 2))
cc_trf(rsq2, se2, 1/100000, 354/467)

#MS+PRS
lm3<- lm(outcome ~ s5+y_pred, data = North)
summary(lm3)
rsq3<-summary(lm3)$r.squared 
se3 <- sqrt((1 - rsq3) / (467 - 3))
cc_trf(rsq3, se3, 1/100000, 354/467)

#MS+PRS+MS*PRS
lm4<- lm(outcome ~ s5+y_pred+s5*y_pred, data = North)
summary(lm4)
rsq4<-summary(lm4)
se4 <- sqrt((1 - rsq4) / (467 - 4))
cc_trf(rsq4, se4, 1/100000, 354/467)



#South###############
South<-subset(pm,PC1>0,select = c("ID","Group","s5","outcome","y_test","pla.value"))
South<-rename(South,y_pred=pla.value)


#MS
lm2<- lm(outcome ~ y_pred, data = South)
summary(lm2)
rsq2<-summary(lm2)$r.squared
se2 <- sqrt((1 - rsq2) / (205 - 2))
cc_trf(rsq2, se2, 1/100000, 120/205)

#MS+PRS
lm3<- lm(outcome ~ s5+y_pred, data = South)
summary(lm3)
rsq3<-summary(lm3)$r.squared 
se3 <- sqrt((1 - rsq3) / (205 - 3))
cc_trf(rsq3, se3, 1/100000, 120/205)

#MS+PRS+MS*PRS
lm4<- lm(outcome ~ s5+y_pred+s5*y_pred, data = South)
summary(lm4)
rsq4<-summary(lm4)$r.squared 
se4 <- sqrt((1 - rsq4) / (205 - 4))
cc_trf(rsq4, se4, 1/100000, 120/205)





################################### plot ########################################

##### barplot1--Figure5F ########
dat1=read_excel("R2_plot.xlsx",sheet = "R2_prs") 

barplot1<-ggplot(data=dat1,aes(x=Group,y=R2l,fill=pT))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("#FFF0F5","#FFE4E1","#FFC0CB","#CD5C5C",
                               "#B22222","#A52A2A","#800000","#500000","#300000"
  ),
  name="")+
  labs(x="Cohort",y="Variance explained (adjusted R2)")+axisSetting
barplot1




##### barplot2--Figure5G ######
dat2=read_excel("R2_plot.xlsx",sheet = "R2_pms") 
barplot2<-ggplot(data=dat2,aes(x=Group,
                               y=R2l.xgb.calib,
                               fill=model))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9),color = "transparent")+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("#B22222",
                               "#f8ac8c",
                               "#9ac9db",
                               "#2878b5"
                               
  ),
  name="")+
  labs(x="Cohort",y="Variance explained (adjusted R2)")+axisSetting
barplot2




#### density plot--Figure5C-D #######
cv.xgb<-read.csv("xgb_pms_data.csv")
cv.lasso<-read.csv("lasso_pms_data.csv")

#LASSO
cv.lasso$Group[cv.lasso$y_test==1]<-"case"
cv.lasso$Group[cv.lasso$y_test==0]<-"control"
plot.lasso.pla<-ggplot(cv.lasso,aes(x=pla.value))+
  geom_density(aes(color = Group),show.legend = T)+
  geom_rug(aes(color = Group))+
  labs(x="risk score")+
  scale_color_manual(values = c("#B22222","#2878b5"))+
  axisSetting

plot.lasso.pla

#XGBoost
cv.xgb$Group[cv.xgb$y_test==1]<-"case"
cv.xgb$Group[cv.xgb$y_test==0]<-"control"
plot.xgb.pla<-ggplot(cv.xgb,aes(x=pla.value))+
  geom_density(aes(color = Group),show.legend = T)+
  geom_rug(aes(color = Group))+
  labs(x="risk score")+
  scale_color_manual(values = c("#B22222","#2878b5"))+
  axisSetting+
  scale_x_continuous( limits=c(0.0, 1.0))
plot.xgb.pla



##### decile plot--Figure5E ###########
dat3<-read_excel("quantile_plot_calib.xlsx",sheet = "plot")
dat3$decile<-as.factor(dat3$decile)
barplot3<-ggplot(data=dat3,aes(x=decile,y=`Odds`,fill=model))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9),color = "transparent")+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("#c82423",
                               #"#ff8884",
                               #"#f8ac8c",
                               #"#9ac9db"
                               "#2878b5"
  ))+
  labs(x="decile",y="odds")+axisSetting
barplot3
