rm(list = ls())
library(glmnet)
library(pROC)
library(xgboost)
library(caret)
library(ROCR)

set.seed(123)
load("D:/Students/2021yaoxm/validation/seed123_300th/data.Rdata")

#################################lasso###############################################
# run 100 times
genes<-list()
pred_tr<-list()
pred_te<-list()
pred_all<-list()
aucs=data.frame(train.lower=NA,train.auc=NA,train.upper=NA,
                test.lower=NA,test.auc=NA,test.upper=NA
                #all.lower=NA,all.auc=NA,all.upper=NA
                )
for (i in 1:10) {
  train_index <- sample(1:nrow(data), nrow(data) * 2 / 3)
  train_data <- data[train_index, ]
  
  x_train <- as.matrix(train_data[,-302052])
  y_train <- as.numeric(train_data[,302052])
  
  #model
  lasso_model <- cv.glmnet(x_train, y_train, family = "binomial", type.measure = "auc",nfolds = 5)
  lasso_model_min <- glmnet(x_train, y_train, alpha = 1, lambda=lasso_model$lambda.min)
  
  #select key genes
  choose_gene_min=rownames(lasso_model_min$beta)[as.numeric(lasso_model_min$beta)!=0]
  genes_matrix<-as.matrix(cbind(choose_gene_min,lasso_model_min$beta@x))
  genes[[i]]<-genes_matrix
  
  #training
  tr_pred<-predict(lasso_model_min, newx = x_train)
  pred_tr[[i]]<-tr_pred
  
  #testing
  test_index <- setdiff(1:nrow(data), train_index)
  test_data <- data[test_index, ]
  x_test <- as.matrix(test_data[,-302052])
  y_test <- as.numeric(test_data[,302052])
  te_pred <- predict(lasso_model_min, newx = x_test)
  pred_te[[i]]<-cbind(te_pred,y_test)
  
  #all sample
  x <- as.matrix(data[,-302053])
  y <- as.numeric(data[,302053])
  all_pred <- predict(lasso_model_min, newx = x)
  pred_all[[i]]<-all_pred
  
  # AUC value
  ci.tr<-ci.auc(as.numeric(y_train), as.numeric(tr_pred))
  ci.te<-ci.auc(as.numeric(y_test), as.numeric(te_pred))
  ci.all<-ci.auc(as.numeric(y), as.numeric(all_pred))
  aucs[i,1:3] <- ci.tr
  aucs[i,4:6] <- ci.te
  aucs[i,7:9] <- ci.all
}

####save result
 write.csv(aucs,"D:/Students/2021yaoxm/validation/seed123_300th/lasso_no.ankle2/aucs123.csv")
 save(pred_tr,file="D:/Students/2021yaoxm/validation/seed123_300th/lasso_new100/pred.train.Rdata")
 save(pred_te,file="D:/Students/2021yaoxm/validation/seed123_300th/lasso_no.ankle2/pred.test.Rdata")
 save(pred_all,file="D:/Students/2021yaoxm/validation/seed123_300th/lasso_new100/pred.all.Rdata")
 save(genes,file ="D:/Students/2021yaoxm/validation/seed123_300th/lasso_no.ankle2/genes.Rdata" )
 
 
 
 
 
 
 ###############################xgboost##########################################
 imp<-list()
 pred_tr<-list()
 pred_te<-list()
 pred_all<-list()
 aucs=data.frame(train.lower=NA,train.auc=NA,train.upper=NA,
                 test.lower=NA,test.auc=NA,test.upper=NA,
                 all.lower=NA,all.auc=NA,all.upper=NA)
for (i in 1:50) {
  train_index <- sample(1:nrow(data), nrow(data) * 2 / 3)
  train_data <- data[train_index, ]
  
  x_train <- as.matrix(train_data[,-302053])
  y_train <- as.numeric(train_data[,302053])
  
  #model
  dtrain <- xgb.DMatrix(data=x_train , label=as.matrix(y_train))
  bstDMatrix <- xgboost(data = dtrain, max.depth = 3,
                        eta = 0.01, nthread = 2, nrounds = 10,
                        objective = "binary:logistic")
  
  names <- bstDMatrix$feature_names
  # select key genes
  importance_matrix <- xgb.importance(names, model = bstDMatrix)
  imp[[i]]<-importance_matrix
  
  #training
  tr_pred<-predict(bstDMatrix, data.matrix(x_train))
  pred_tr[[i]]<-tr_pred
  
  #testing
  test_index <- setdiff(1:nrow(data), train_index)
  test_data <- data[test_index, ]
  x_test <- as.matrix(test_data[,-302053])
  y_test <- as.numeric(test_data[,302053])
  te_pred <-predict(bstDMatrix, data.matrix(x_test))
  pred_te[[i]]<-cbind(rownames(test_data),te_pred,y_test)
  
  #all sample
  x <- as.matrix(data[,-302053])
  y <- as.numeric(data[,302053])
  all_pred <- predict(bstDMatrix, data.matrix(x))
  pred_all[[i]]<-all_pred
  
  # AUC value
  ci.tr<-ci.auc(as.numeric(y_train), as.numeric(tr_pred))
  ci.te<-ci.auc(as.numeric(y_test), as.numeric(te_pred))
  ci.all<-ci.auc(as.numeric(y), as.numeric(all_pred))
  aucs[i,1:3] <- ci.tr
  aucs[i,4:6] <- ci.te
  aucs[i,7:9] <- ci.all
}

#save result
write.csv(aucs,"D:/Students/2021yaoxm/validation/seed123_300th/xgb_new50/aucs123.csv")
save(pred_tr,file="D:/Students/2021yaoxm/validation/seed123_300th/xgb_new50/pred.train.Rdata")
save(pred_te,file="D:/Students/2021yaoxm/validation/seed123_300th/xgb_new50/pred.test.Rdata")
save(pred_all,file="D:/Students/2021yaoxm/validation/seed123_300th/xgb_new50/pred.all.Rdata")
save(imp,file="D:/Students/2021yaoxm/validation/seed123_300th/xgb_new50/genes.Rdata")
