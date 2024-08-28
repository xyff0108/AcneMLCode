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

############################### function to change label to Y/N#####################################
acnerank_num_train<-readRDS('data/acnerank_num_train.rds')
acnerank_num_test<-readRDS('data/acnerank_num_test.rds')


acnerank_num_train_rose<-ROSE(Acnerank ~., data =acnerank_num_train[,-1], seed = 1)$data
acnerank_num_train_rose$Acnerank<-change2Level(acnerank_num_train_rose$Acnerank)
skim(acnerank_num_train_rose)

acnerank_num_test$Acnerank<-change2Level(acnerank_num_test$Acnerank)
acnerank_num_test$Acne_rank<-NULL
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
i.Acnerank <- 60
for (i in c(1:5)) { 
  currentmodel <- mymodels[i]
  print(paste0("current running model is ",mymodels[i],sep=""))
  tic()
  cmodel <- train(acnerank_num_train_rose[,-i.Acnerank ], 
                 acnerank_num_train_rose$Acnerank,
                  metric = "ROC",
                  method = currentmodel, 
                  trControl = tunetrl, #use adaptive resampling method
                  #verbose = FALSE, 
                  tuneLength = ntunelength) #each bootstrap tuning length # Maximum number of hyperparameter combinations
  runningtime <- toc() #save runningtime
  cname <- paste("TuningModel_num_acnerank_rose.",currentmodel,".rds",sep = "")
  saveRDS(cmodel,cname)
  mytime<-toc()
  saveRDS(mytime,paste("runningtime_num_acnerank_rose.",currentmodel,".rds",sep = ""))
  modellist[[i]]<-cmodel #save to the big list
  }

#######################################################

################ Compare ML models performance #######################################################
RF<-readRDS('TuningModel_num_acnerank_rose.rf.rds')
GBDT<-readRDS('TuningModel_num_acnerank_rose.gbm.rds')
KNN<-readRDS('TuningModel_num_acnerank_rose.knn.rds')
SVM<-readRDS('TuningModel_num_acnerank_rose.svmLinear.rds')
NN<-readRDS('TuningModel_num_acnerank_rose.avNNet.rds')


mymodels<-c('RF', 'GBDT','KNN', 'SVM', 'NN')
xrow<-nrow(RF$results)+
  nrow(GBDT$results)+
  nrow(KNN$results)+
  nrow(SVM$results)+
  nrow(NN$results)

performance_Acnerank <- data.frame(matrix(NA, nrow = xrow, ncol = 4))
colnames(performance_Acnerank) <- c("Model","ROC","Sens","Spec")
performance_Acnerank$Model <- c(rep(mymodels[1],nrow(RF$results)),
                    rep(mymodels[2],nrow(GBDT$results)),
                    rep(mymodels[3],nrow(KNN$results)),
                    rep(mymodels[4],nrow(SVM$results)),
                    rep(mymodels[5],nrow(NN$results)))
performance_Acnerank$ROC <-c(rep(RF$results$ROC),
                 rep(GBDT$results$ROC),
                 rep(KNN$results$ROC),
                 rep(SVM$results$ROC),
                 rep(NN$results$ROC))
performance_Acnerank$Sens <-c(rep(RF$results$Sens),
                  rep(GBDT$results$Sens),
                  rep(KNN$results$Sens),
                  rep(SVM$results$Sens),
                  rep(NN$results$Sens))
performance_Acnerank$Spec <-c(rep(RF$results$Spec),
                  rep(GBDT$results$Spec),
                  rep(KNN$results$Spec),
                  rep(SVM$results$Spec),
                  rep(NN$results$Spec))

dir.create("MLresults")
write.csv(performance_Acnerank,"MLresults/performance_Acnerank.csv")

