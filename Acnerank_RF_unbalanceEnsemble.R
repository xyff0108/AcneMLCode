windowsFonts(TNM = windowsFont("Times New Roman"))
palette <- c("#313695", "#D73027")

library(doParallel) 
cl <- makeCluster(2)  
registerDoParallel(cl) 

library(openxlsx)
library(skimr)
library(DataExplorer)
library(caret)
library(randomForest)
library(pROC)
library(dplyr)
library(ROSE)
library(stats)
############################### function to change label to Y/N#####################################
change2Level <- function(y){y <- factor(y, levels=c(0, 1, "No", "Yes"))
y[y == '0'] <- 'No'
y[y == '1'] <- 'Yes'
y <- factor(y, levels=c("No", "Yes"))}

#########################read input data ###########################################

acnerank_num_train<-readRDS('data/acnerank_num_train.rds')
acnerank_num_test<-readRDS('data/acnerank_num_test.rds')
acnerank_num_test$Acne_rank<-NULL

acnerank_num_train$Acnerank<-change2Level(acnerank_num_train$Acnerank)
acnerank_num_test$Acnerank<-change2Level(acnerank_num_test$Acnerank)


#########################process unbalanced data##############################################
dir.create("balanceRF")
acnerank_num_train_balanced_both <- ovun.sample(Acnerank~., data =acnerank_num_train[,-1],
                            N = nrow(acnerank_num_train),p=0.5,seed=1,method = "both")$data

acnerank_num_train_balanced_under <- ovun.sample(Acnerank~., data =acnerank_num_train[,-1],
                                                seed=1,method = "under")$data

acnerank_num_train_balanced_over <- ovun.sample(Acnerank~., data =acnerank_num_train[,-1],
                                                 seed=1,method = "over")$data

acnerank_num_train_balanced_rose <- ROSE(Acnerank~., data =acnerank_num_train[,-1], seed = 1)$data


trainlist <- list(acnerank_num_train_balanced_both,acnerank_num_train_balanced_under,
                  acnerank_num_train_balanced_over,acnerank_num_train_balanced_rose)
names(trainlist) <- c("both","under","over","rose")

saveRDS(trainlist,file = "balanceRF/acnerank_num_train_balancedlist.rds")


#################################build rf model##############################################
i.Acnerank <- 60
x_names<-colnames(acnerank_num_train_balanced_both[,-i.Acnerank])
form_cls<-as.formula(paste0("Acnerank~",paste(x_names,collapse = "+")))
form_cls
WT <- c(1,1,1,1) # weights of different balancing method
CUTOFF <- 0.5 # probability threshold

hyper_grid <- expand.grid(
  ntree =(2:10)*250) 

#predlist.ensemble <- list()
hyper_grid <-c(2:10)*250

for (i in c(1:9)) {
  cntree <- hyper_grid[i]
  currentmodel <- i
  print(paste0("current running model is ",i,sep="")) 
  print(cntree)
  
  set.seed(100*i)
  cmodel.both <-randomForest(
    form_cls,
    data=acnerank_num_train_balanced_both,
    ntree=cntree, 
    mtry=5,         #depending on the final values in adaptive boosting process
    importance=T   
  )
  
  set.seed(100)
  cmodel.under <-randomForest(
    form_cls,
    data=acnerank_num_train_balanced_under,
    ntree=cntree, 
    mtry=5,         #depending on the final values in adaptive boosting process
    importance=T   
  )
  
  set.seed(100)
  cmodel.over <-randomForest(
    form_cls,
    data=acnerank_num_train_balanced_over,
    ntree=cntree, 
    mtry=5,         #depending on the final values in adaptive boosting process
    importance=T   
  )
  
  set.seed(100)
  cmodel.rose <- randomForest(
    form_cls,
    data=acnerank_num_train_balanced_rose,
    ntree=cntree, 
    mtry=5,         #depending on the final values in adaptive boosting process
    importance=T   
  )
  

  cmodellist <- list(cmodel.both,cmodel.under,cmodel.over,cmodel.rose)
  names(cmodellist) <- c("both","under","over","rose")
  
  saveRDS(cmodellist,file = paste0("balanceRF/SevereAcne_tune_balanceRFlist.",i,".rds"))
  
  
  rm(cmodel.both)
  rm(cmodel.under)
  rm(cmodel.over)
  rm(cmodel.rose)

}


#################################predict testset & perfromance analysis##############################################
for(i in c(1:9)){
  setwd(resultpath)
  print(i)
  cmodellist <- readRDS(file = paste0("balanceRF/SevereAcne_tune_balanceRFlist.",i,".rds"))
  cmodel.both <- cmodellist$both
  cmodel.under <- cmodellist$under
  cmodel.over <- cmodellist$over
  cmodel.rose <- cmodellist$rose
  
  pred.both.prob <- predict(cmodel.both, acnerank_num_test,type = 'prob')
  pred.under.prob <- predict(cmodel.under, acnerank_num_test,type = 'prob')
  pred.over.prob <- predict(cmodel.over, acnerank_num_test,type = 'prob')
  pred.rose.prob <- predict(cmodel.rose, acnerank_num_test,type = 'prob')
  
  pred.ensemble.Nomat <- cbind(pred.both.prob[,1],pred.under.prob[,1],pred.over.prob[,1],pred.rose.prob[,1])
  pred.ensemble.Yesmat <- cbind(pred.both.prob[,2],pred.under.prob[,2],pred.over.prob[,2],pred.rose.prob[,2])
  
  pred.ensemble.No <- apply(pred.ensemble.Nomat[,1:4],1,function(x) weighted.mean(x,WT))
  pred.ensemble.Yes <- apply(pred.ensemble.Yesmat[,1:4],1,function(x) weighted.mean(x,WT))
  pred.ensemble.prob <- data.frame(cbind(pred.ensemble.No,pred.ensemble.Yes))
  colnames(pred.ensemble.prob) <- colnames(pred.both.prob)
  
  pred.ensemble <- ifelse(pred.ensemble.prob$No>0.5,"No","Yes")
  pred.ensemble <- factor(pred.ensemble,levels = c("No","Yes"))
  
  pred.both <- predict(cmodel.both, acnerank_num_test,type = 'response')
  pred.under <- predict(cmodel.under, acnerank_num_test,type = 'response')
  pred.over <- predict(cmodel.over, acnerank_num_test,type = 'response')
  pred.rose <- predict(cmodel.rose, acnerank_num_test,type = 'response')
  
  confusion.both <- confusionMatrix(acnerank_num_test$Acnerank,pred.both,positive='Yes', mode = 'everything')
  confusion.under <-confusionMatrix(acnerank_num_test$Acnerank,pred.under,positive='Yes', mode = 'everything')
  confusion.over <-confusionMatrix(acnerank_num_test$Acnerank,pred.over,positive='Yes', mode = 'everything')
  confusion.rose <-confusionMatrix(acnerank_num_test$Acnerank,pred.rose,positive='Yes', mode = 'everything')
  
  confusion.ensemble <- confusionMatrix(acnerank_num_test$Acnerank,pred.ensemble,positive='Yes', mode = 'everything')
  
  
  
  roc.both <- roc(acnerank_num_test$Acnerank,pred.both.prob[,2],aur=TRUE,ci=TRUE)
  roc.under <- roc(acnerank_num_test$Acnerank,pred.under.prob[,2],aur=TRUE,ci=TRUE)
  roc.over <- roc(acnerank_num_test$Acnerank,pred.over.prob[,2],aur=TRUE,ci=TRUE)
  roc.rose <-  roc(acnerank_num_test$Acnerank,pred.rose.prob[,2],aur=TRUE,ci=TRUE)
  roc.ensemble <- roc(acnerank_num_test$Acnerank,pred.ensemble.prob[,2],aur=TRUE,ci=TRUE)
  
  
  
  overall.merge <- cbind(confusion.ensemble$overall,
                         confusion.both$overall,confusion.under$overall,
                         confusion.over$overall,confusion.rose$overall)
  colnames(overall.merge) <- c("Ensemble","Both","Under","Over","ROSE")
  
  
  byClass.merge <- cbind(confusion.ensemble$byClass,
                         confusion.both$byClass,confusion.under$byClass,
                         confusion.over$byClass,confusion.rose$byClass)
  colnames(byClass.merge) <- c("Ensemble","Both","Under","Over","ROSE")
  
  
  aucci.merge <- cbind(roc.ensemble$ci,
                       roc.both$ci,roc.under$ci,
                       roc.over$ci,roc.rose$ci)
  colnames(aucci.merge) <- c("Ensemble","Both","Under","Over","ROSE")
  rownames(aucci.merge) <- c("aucciLower","auccimedian","aucciupper")
  
  
  auc.merge <- cbind(roc.ensemble$auc,
                     roc.both$auc,roc.under$auc,
                     roc.over$auc,roc.rose$auc)
  rownames(auc.merge) <- c("auc")
  colnames(auc.merge) <- c("Ensemble","Both","Under","Over","ROSE")
  
  
  perRF.merge <- rbind.data.frame(overall.merge,byClass.merge,auc.merge,aucci.merge)
  perRF.merge <- round(perRF.merge,5)
  setwd(outpath)
  write.csv(perRF.merge,file = paste0("perRF.",i,".csv"))
  #saveRDS(perRF.merge,file = paste0("perRF.",i,".rds"))

}

## merge performance from 1:9 RF
setwd(outpath)
perRFlist <- list()
for (i in c(1:9)) {
  perRF <- read.csv(file = paste0("perRF.",i,".csv"),row.names = 1)
  perRFlist[[i]] <- perRF
}
names(perRFlist) <- paste0("RF",c(1:9))
perRFlist.ensemble <- lapply(perRFlist,function(x) x$Ensemble)
perRFmat <- data.frame(do.call(rbind,perRFlist.ensemble))
colnames(perRFmat) <- rownames(perRFlist$RF1)
permat <- round(perRFmat,5)
permat$ntree <- (2:10)*250
permat$RF <- rownames(permat)
permat$Models <- paste0(permat$RF,"(ntree=",permat$ntree,")")
permat$'AUC (95%CI)' <- paste0(sprintf("%.5f",perRFmat$auc),"(",sprintf("%.5f",perRFmat$aucciLower),"-",sprintf("%.5f",perRFmat$aucciupper),")")
permat$'Accuracy (95%CI)' <- paste0(sprintf("%.5f",perRFmat$Accuracy),"(",sprintf("%.5f",perRFmat$AccuracyLower),"-",sprintf("%.5f",perRFmat$AccuracyUpper),")")
permat.final <- permat[,c(25:27,8:9,14,2)]
write.csv(permat.final,file = "permat.final.csv")

############################# plot optimal model ROC curve ########################
library(RColorBrewer)
mycols <- brewer.pal(5,"Set1")
setwd(resultpath)
i = 3
WT = c(1,1,1,1)
bestrflist <- readRDS(file = paste0("SevereAcne_tune_balanceRFlist.",i,".rds"))
cmodel.both <- bestrflist$both
cmodel.under <- bestrflist$under
cmodel.over <- bestrflist$over
cmodel.rose <- bestrflist$rose

pred.both.prob <- predict(cmodel.both, acnerank_num_test,type = 'prob')
pred.under.prob <- predict(cmodel.under, acnerank_num_test,type = 'prob')
pred.over.prob <- predict(cmodel.over, acnerank_num_test,type = 'prob')
pred.rose.prob <- predict(cmodel.rose, acnerank_num_test,type = 'prob')

pred.ensemble.Nomat <- cbind(pred.both.prob[,1],pred.under.prob[,1],pred.over.prob[,1],pred.rose.prob[,1])
pred.ensemble.Yesmat <- cbind(pred.both.prob[,2],pred.under.prob[,2],pred.over.prob[,2],pred.rose.prob[,2])

pred.ensemble.No <- apply(pred.ensemble.Nomat[,1:4],1,function(x) weighted.mean(x,WT))
pred.ensemble.Yes <- apply(pred.ensemble.Yesmat[,1:4],1,function(x) weighted.mean(x,WT))
pred.ensemble.prob <- data.frame(cbind(pred.ensemble.No,pred.ensemble.Yes))
colnames(pred.ensemble.prob) <- colnames(pred.both.prob)

pred.ensemble <- ifelse(pred.ensemble.prob$No>0.5,"No","Yes")
pred.ensemble <- factor(pred.ensemble,levels = c("No","Yes"))

roc.both <- roc(acnerank_num_test$Acnerank,pred.both.prob[,2],aur=TRUE,ci=TRUE)
roc.under <- roc(acnerank_num_test$Acnerank,pred.under.prob[,2],aur=TRUE,ci=TRUE)
roc.over <- roc(acnerank_num_test$Acnerank,pred.over.prob[,2],aur=TRUE,ci=TRUE)
roc.rose <-  roc(acnerank_num_test$Acnerank,pred.rose.prob[,2],aur=TRUE,ci=TRUE)
roc.ensemble <- roc(acnerank_num_test$Acnerank,pred.ensemble.prob[,2],aur=TRUE,ci=TRUE)

aucs <- c(roc.ensemble$auc[1],roc.both$auc[1],roc.under$auc[1],roc.over$auc[1],roc.rose$auc[1])
aucs <- round(aucs,5)
legendtext <- paste0(c("Ensemble","Both","Under","Over","ROSE")," AUC=",aucs)
setwd(outpath)
tiff("bestRF.AUCcurve.tiff",width = 6,height = 5,units = "in",res = 300)
plot(roc.ensemble,col=mycols[1],print.auc=F,
     print.thres=F,
     grid=c(0.2,0.2),grid.col=c("grey","grey"))

plot(roc.both,add=TRUE,col=mycols[2],print.auc=F)
plot(roc.under,add=TRUE,col=mycols[3],print.auc=F)
plot(roc.over,add=TRUE,col=mycols[4],print.auc=F)
plot(roc.rose,add=TRUE,col=mycols[5],print.auc=F)
legend("bottomright", legend=legendtext,
       col=mycols,lty=1,cex=0.7)
title("RF3(ntree=1000)",line = 3)
dev.off()