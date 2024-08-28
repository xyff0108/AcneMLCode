palette <- c("#313695", "#D73027")

windowsFonts(TNM = windowsFont("Times New Roman"))

library(doParallel) 
cl <- makeCluster(2)  
registerDoParallel(cl) 

library(openxlsx)
library(caret)
library(skimr)
library(DataExplorer)
library(tictoc)
library('plyr')
library(ROSE)
library(psych)

change2Level <- function(y){
  y <- factor(y, levels=c(0, 1, "No", "Yes"))
  y[y == '0'] <- 'No'
  y[y == '1'] <- 'Yes'
  y <- factor(y, levels=c("No", "Yes"))
}

############################### read.csv data#####################################
acne_num_train<-readRDS('data/acne_num_train.rds')
acne_num_test<-readRDS('data/acne_num_test.rds')

################################ tuning parameters############################################### 
## build adaptive resampling tuning 
## max running times = nboot x ntunelength
nboot = 10  #bootstrap times 
ntunelength = 20 
tunetrl <- trainControl(method = "adaptive_boot", number = nboot, # bootstrap n times
                        adaptive = list(min = 5, alpha = 0.05, #futility loss 0.05
                        method = "gls", #can use BT for alternative
                        complete = TRUE),
                        classProbs = TRUE,
                        search = "random",
                        #sampling = "rose",
                        summaryFunction = twoClassSummary)

# caret models used to compare
mymodels <-  c('rf','gbm','avNNet', 'svmLinear','knn')
modellist <- list()


####### running the adaptive resampling method for each model######### #############################################
i.Acne <- match("Acne",colnames(acne_num_train))
for (i in c(1:5)) { 
  currentmodel <- mymodels[i]
  print(paste0("current running model is ",mymodels[i],sep=""))
  
  tic()
  cmodel <- train(acne_num_train[,-i.Acne], 
                 acne_num_train$Acne,
                  metric = "ROC",
                  method = currentmodel, 
                  trControl = tunetrl, #use adaptive resampling method
                  #verbose = FALSE, 
                  tuneLength = ntunelength) #each bootstrap tuning length # Maximum number of hyperparameter combinations
  runningtime <- toc() #save runningtime
  cname <- paste("TuningModel_num_acne.",currentmodel,".rds",sep = "")
  saveRDS(cmodel,cname)
  mytime<-toc()
  saveRDS(mytime,paste("runningtime_num_acne.",currentmodel,".rds",sep = ""))
  modellist[[i]]<-cmodel #save to the big list
  }


################ Compare ML models performance #######################################################
RF<-readRDS('TuningModel_num_acne.rf.rds')
GBDT<-readRDS('TuningModel_num_acne.gbm.rds')
KNN<-readRDS('TuningModel_num_acne.knn.rds')
SVM<-readRDS('TuningModel_num_acne.svmLinear.rds')
NN<-readRDS('TuningModel_num_acne.avNNet.rds')



mymodels<-c('RF', 'GBDT','KNN', 'SVM', 'NN')
xrow<-nrow(RF$results)+
  nrow(GBDT$results)+
  nrow(KNN$results)+
  nrow(SVM$results)+
  nrow(NN$results)

performance_Acne <- data.frame(matrix(NA, nrow = xrow, ncol = 4))
colnames(performance_Acne) <- c("Model","ROC","Sens","Spec")
performance_Acne$Model <- c(rep(mymodels[1],nrow(RF$results)),
                    rep(mymodels[2],nrow(GBDT$results)),
                    rep(mymodels[3],nrow(KNN$results)),
                    rep(mymodels[4],nrow(SVM$results)),
                    rep(mymodels[5],nrow(NN$results)))
performance_Acne$ROC <-c(rep(RF$results$ROC),
                 rep(GBDT$results$ROC),
                 rep(KNN$results$ROC),
                 rep(SVM$results$ROC),
                 rep(NN$results$ROC))
performance_Acne$Sens <-c(rep(RF$results$Sens),
                  rep(GBDT$results$Sens),
                  rep(KNN$results$Sens),
                  rep(SVM$results$Sens),
                  rep(NN$results$Sens))
performance_Acne$Spec <-c(rep(RF$results$Spec),
                  rep(GBDT$results$Spec),
                  rep(KNN$results$Spec),
                  rep(SVM$results$Spec),
                  rep(NN$results$Spec))

dir.create(paste0("MLresults"))
write.csv(performance_Acne,"MLresults/performance_Acne.csv")
