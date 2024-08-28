palette <- c("#313695", "#D73027")


windowsFonts(TNM = windowsFont("Times New Roman"))

library(doParallel) 
cl <- makeCluster(2)  
registerDoParallel(cl) 

#stopCluster(cl)
library(openxlsx)
library(skimr)
library(DataExplorer)
library(caret)
library(randomForest)
library(pROC)
library(dplyr)
library(ROSE)
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
acnerank_num_train_rose<-ROSE(Acnerank~., data =acnerank_num_train[,-1], seed = 1)$data

#################################build rf model##############################################
i.Acnerank <- 60
x_names<-colnames(acnerank_num_train_rose[,-i.Acnerank])
form_cls<-as.formula(paste0("Acnerank~",paste(x_names,collapse = "+")))
form_cls
 
###########################grid research#####################
hyper_grid <- expand.grid(
  ntree =(2:10)*250) 

dir.create("RF_optimize")
for (i in 1:nrow( hyper_grid)){
  set.seed(100)
  currentmodel <- i
  print(paste0("current running model is ",i,sep="")) 

  cmodel<-randomForest(
  form_cls,
  data=acnerank_num_train_rose,
  ntree=hyper_grid$ntree[i], 
  mtry=5,         #depending on the final values in adaptive boosting process
  importance=T   
  )
  cname <- paste("RF_optimize/SevereAcne_tuneRF.",i,".rds",sep = "")
  saveRDS(cmodel,cname)
}


##################modelling performance compare#######################################
perRFlist <- list()
for (i in c(1:9)) {
  rf_acnerank<-readRDS(paste("RF_optimize/SevereAcne_tuneRF.",i,".rds",sep = ""))  
  
  testpredprop_acnerank<-predict(rf_acnerank,newdata = acnerank_num_test,type = 'prob')
  rfPred_acnerank<-predict(rf_acnerank,acnerank_num_test,type="response")
  rfperf_acnerank<-confusionMatrix(acnerank_num_test$Acnerank,rfPred_acnerank,positive='Yes', mode = 'everything')
  rocrf_acnerank<-roc(acnerank_num_test$Acnerank,testpredprop_acnerank[,2],aur=TRUE,ci=TRUE)
  
  overall <- data.frame(rfperf_acnerank$overall)
  overall_1$names<-rownames(overall_1) 
  overall_1$performance<-overall_1 $rfperf1_acnerank.overall
  overall_1 $rfperf1_acnerank.overall<-NULL
  
  overall <- data.frame(rfperf_acnerank$overall)
  overall$names<-rownames(overall) 
  overall$performance<-overall$rfperf_acnerank.overall
  overall$rfperf_acnerank.overall<-NULL
  
  
  byClass <-data.frame(rfperf_acnerank$byClass)
  byClass$names<-rownames(byClass)
  byClass$performance<-byClass$rfperf_acnerank.byClass
  byClass$rfperf_acnerank.byClass<-NULL
  
  
  aucci<-data.frame(rocrf_acnerank$ci)
  aucci$names<-c("aucciLower","auccimedian","aucciupper")
  aucci$performance<-aucci$rocrf_acnerank.ci
  aucci$rocrf_acnerank.ci<-NULL
  
  auc<-data.frame(rocrf_acnerank$auc)
  auc$names<-c("auc")
  auc$performance<-auc$rocrf_acnerank.auc
  auc$rocrf_acnerank.auc<-NULL
  
  perRF<-rbind.data.frame(overall,byClass,auc,aucci)
  perRF$performance<-round(perRF$performance,5)
  
  perRFlist[[i]] <- perRF
}

per_acnerank<-do.call(cbind,perRFlist)
per_acnerank<-per_acnerank[,c(1,2,4,6,8,10,12,14,16,18)]
per_acnerank<-data.frame(t(per_acnerank))
per_acnerank<-per_acnerank[2:10,]
per_acnerank$AUC<- per_acnerank$X1
per_acnerank$AUCLower<- per_acnerank$X11
per_acnerank$AUCUpper<- per_acnerank$X3
per_acnerank$X1<-NULL
per_acnerank$X11<-NULL
per_acnerank$X2<-NULL
per_acnerank$X3<-NULL
per_acnerank1<-as.data.frame(lapply(per_acnerank,as.numeric))
per_acnerank1$Model<-c("RF1","RF2","RF3","RF4","RF5","RF6","RF7","RF8","RF9")

optimalRF_acnerank<-per_acnerank1$Model[which(per_acnerank1$F1==max(per_acnerank1$F1))] # selecting the optimal RF using largest AUC
optimalRF_acnerank

write.csv(per_acnerank1, 'RF_optimize/RFperformance_Acnerank.csv')


################################# feature importance unsing the optimal RF  ##############################################
besti <- 3
rfbest_acnerank<-readRDS(paste("RF_optimize/SevereAcne_tuneRF.",besti,".rds",sep = ""))  

acnerankrf_importance_scale<-data.frame(importance(rfbest_acnerank,scale = TRUE), check.names = FALSE)
acnerankrf_importance<-data.frame(importance(rfbest_acnerank), check.names = FALSE)
acnerankrf_importance_sort<- data.frame(acnerankrf_importance[order(acnerankrf_importance$MeanDecreaseAccuracy,decreasing = TRUE),])
acnerankrf_importance_sort$features=rownames(acnerankrf_importance_sort)
acnerankrf_importance_sort$relativeinf<-round(acnerankrf_importance_sort$MeanDecreaseAccuracy/sum(acnerankrf_importance_sort$MeanDecreaseAccuracy)*100,2)

acnerankrf_importance_sort$group<-ifelse(acnerankrf_importance_sort$features%in%c("Gender","Age","Ethnicity","Household_income","Grade"),"DSCs",
                                         ifelse(acnerankrf_importance_sort$features%in%c("Skin_type","Sensitive_skin","Genetic","Menstrual_blood_volume",
                                                                                         "Menstrual_colic","Menstrual_cycle","","PHQ9","BMI"),"BMAs",
                                                ifelse(acnerankrf_importance_sort$features%in%c("Den_barbecue_shop","Den_busstop","Den_cafeterias_shop","Den_fruit_shop",
                                                                                                "Den_hotpot_restaurant","Den_intersection","Den_KFC_McDonald.s",
                                                                                                "Den_milktea_shop","Den_roadlength","Denpop","NDVI"),"BEs",
                                                       ifelse(acnerankrf_importance_sort$features%in%c("Daylight","Temperature","Relative_humidity",
                                                                                                       "CO","NO2","O3","PM2.5","SO2"),"NEs","LFs"))))
write.csv(acnerankrf_importance_sort,'RF_optimize/acnerankrf_importance.csv')
