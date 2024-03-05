rm(list = ls())

library(pROC)
library(caret)
library(ROCR)
library(dplyr)

##lasso

# input data
load('D:/Students/2021yaoxm/validation/seed123_300th/lasso_new100/pred.test.Rdata')

# 
for (i in 1:length(pred_te)) {
  data_df <- pred_te[[i]]

  predictions <- as.numeric(data_df[, 1])   
  actuals <- data_df[, 2]
  
  # calculate AUC value
  pred<-ROCR::prediction(predictions,actuals)
  perf<-ROCR::performance(pred,"tpr","fpr")
 
  # plot ROC curve
  plot(perf,lwd=1,lty=1,avg="vertical", add=TRUE,
       col= "pink")
}

# mean ROC curve
#roc_average <- 0.788
pred_te100<-as.data.frame(pred_te[[58]])
pred100<-ROCR::prediction(as.numeric((pred_te100$s0)),as.numeric(pred_te100$y_test))
perf100<-ROCR::performance(pred100,"tpr","fpr")

plot(perf100,lwd=1,lty=1,avg="vertical", add=TRUE,col= "red")
legend("bottomright", legend = c("Individual ROC", "Average ROC"), 
       col = c(rainbow(length(pred_te))[1], "red"), 
       lty = c(1, 1), 
       lwd = c(1, 2))


#xgb
# input data
load('D:/Students/2021yaoxm/validation/seed123_300th/xgb_new100/pred.test.Rdata')

#
for (i in 1:length(pred_te)) {
  data_df <- pred_te[[i]]
  
  predictions <- as.numeric(data_df[, 2])   
  actuals <- data_df[, 3]
  
  pred<-ROCR::prediction(predictions,actuals)
  perf<-ROCR::performance(pred,"tpr","fpr")
  
  plot(perf,lwd=1,lty=1,avg="vertical", add=TRUE,
       col= "lightblue")
}

# mean ROC curve
#roc_average <- 0.838
pred_te100<-as.data.frame(pred_te[[79]])
pred100<-ROCR::prediction(as.numeric((pred_te100$te_pred)),as.numeric(pred_te100$y_test))
perf100<-ROCR::performance(pred100,"tpr","fpr")

plot(perf100,lwd=1,lty=1,avg="vertical", add=TRUE,col= "blue")
legend("bottomright", legend = c("Individual ROC", "Average ROC"), 
       col = c(rainbow(length(pred_te))[1], "blue"), 
       lty = c(1, 1), 
       lwd = c(1, 2))
