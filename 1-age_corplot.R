rm(list = ls())

########### age_DNAm age ###########################################################################
library(readxl)
library(plyr)
library(gridExtra)
library(ggplot2)
library(ggpubr)
cols <- c("blue", "#490C65", "#FFA300", "#0081C9", "#46A040", "#F6313E")
count_tmp  <- 1
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
  legend.position ="none",
  legend.title=element_blank(),
  plot.title = element_text(hjust = 0.5,color = 'black',face = 'bold', family = "serif")
)
age_data<-read_excel("D:/Desktop/ALS_202308/Fig_Tab/age/pd_689.xlsx",sheet="case_429")
cor_ave <- cor(age_data$ave_age, age_data$age)
cor_ave_p <- cor.test(age_data$ave_age, age_data$age)[3]
cor1 <- cor(age_data$Hannum, age_data$age)
cor_1_p <- cor.test(age_data$Hannum, age_data$age)[3]
cor2 <- cor(age_data$Horvath, age_data$age)
cor_2_p <- cor.test(age_data$Horvath, age_data$age)[3]
cor3 <- cor(age_data$Horvath2, age_data$age)
cor_3_p <- cor.test(age_data$Horvath2, age_data$age)[3]
cor4 <- cor(age_data$Levine, age_data$age)
cor_4_p <- cor.test(age_data$Levine, age_data$age)[3]
cor5 <- cor(age_data$Zhangqixin, age_data$age)
cor_5_p <- cor.test(age_data$Zhangqixin, age_data$age)[3]
plot <- ggplot(age_data, aes(age, Hannum)) + geom_point(shape = 20, size = 2, col = "blue") 
plot <- plot + axisSetting
plot <- plot + geom_abline(slope = 1, intercept = 0)
plot <- plot + xlab("chronological age") + scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  ylab("DNAm Age(Hannum)")+stat_cor(method = "pearson",
                                    label.x = 50, label.y = 30)

plot

dev.off()

########## DNAm_age&DNAm_age ####################################################################################
#install.packages("GGally")
library(GGally)
library(ggplot2)

cor1 <- cor(age_data$Hannum, age_data$Horvath)#0.94
cor_1_p <- cor.test(age_data$Hannum, age_data$Horvath)[3]#4.01e-196
cor2 <- cor(age_data$Hannum, age_data$Horvath2)#0.96
cor_2_p <- cor.test(age_data$Hannum, age_data$Horvath2)[3]#1.25e-234
cor3 <- cor(age_data$Hannum, age_data$Levine)#0.91
cor_3_p <- cor.test(age_data$Hannum, age_data$Levine)[3]#3.6e-165
cor4 <- cor(age_data$Horvath, age_data$Horvath2)#0.94
cor_4_p <- cor.test(age_data$Horvath, age_data$Horvath2)[3]#2.78e-204
cor5 <- cor(age_data$Horvath, age_data$Levine)#0.88
cor_5_p <- cor.test(age_data$Horvath, age_data$Levine)[3]#4.59e-139
cor6 <- cor(age_data$Horvath2, age_data$Levine)#0.91
cor_6_p <- cor.test(age_data$Horvath2, age_data$Levine)[3]#8.66e-163

df<-data
ggscatter <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  df <- data.frame(x = x, y = y)
  sp1 <- ggplot(df, aes(x=x, y=y)) +
    geom_point(size=0.5,color="#F6313E") +
    geom_abline(intercept = 0, slope = 1, col = 'black')
  return(sp1)
}

ggdehist <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  df <- data.frame(x = x)
  dh1 <- ggplot(df, aes(x=x)) +
    geom_histogram(aes(y=..density..), bins = 50, fill = "#46A040", color='transparent', alpha=.4) +
    geom_density(aes(y=..density..)) + 
    theme_minimal()
  return(dh1)
}

ggpairs(df,
        lower = list(continuous = wrap(ggscatter)),
        diag = list(continuous = wrap(ggdehist)),
        upper = "blank") + 
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text =  element_text(color='black'))+axisSetting 
  







