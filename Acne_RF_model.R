palette <- c("#313695", "#D73027")


windowsFonts(TNM = windowsFont("Times New Roman"))

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
############################### function to change label to Y/N#####################################
change2Level <- function(y){
  y <- factor(y, levels=c(0, 1, "No", "Yes"))
  y[y == '0'] <- 'No'
  y[y == '1'] <- 'Yes'
  y <- factor(y, levels=c("No", "Yes"))
}
#########################read input data ###########################################
acne_num_train<-readRDS(paste0("data/acne_num_train.rds"))
acne_num_test<-readRDS(paste0("data/acne_num_test.rds"))
#################################build rf model##############################################
i.Acne <- 2
x_names<-colnames(acne_num_train[,-i.Acne])
form_cls<-as.formula(paste0("Acne~",paste(x_names,collapse = "+")))
form_cls

###########################grid search#####################
hyper_grid <- expand.grid(
  ntree =(2:10)*250) 

dir.create("RF_optimize")

for (i in 1:nrow( hyper_grid)){
  set.seed(100)
  currentmodel <- i
  print(paste0("current running model is ",i,sep="")) 

  cmodel<-randomForest(
  form_cls,
  data=acne_num_train,
  ntree=hyper_grid$ntree[i],    
  mtry=9,         #depending on the final values in adaptive boosting process
  importance=T   
  )
  cname <- paste("RF_optimize/acne_tuneRF.",i,".rds",sep = "")
  saveRDS(cmodel,cname, )
}



###########################RF performance compare#####################
perRFlist <- list()
for (i in c(1:9)) {
  fit_rf_cls <- readRDS(paste("RF_optimize/acne_tuneRF.",i,".rds",sep = ""))
  testpredprop<-predict(fit_rf_cls,newdata = acne_num_test,type = 'prob')
  
  rfPred<-predict(fit_rf_cls, acne_num_test,type="response")
  rfperf<-confusionMatrix(acne_num_test$Acne,rfPred,positive='Yes', mode="everything")
  rocrf<-roc(acne_num_test$Acne,testpredprop[,2],aur=TRUE,ci=TRUE)
  
  
  overall <- data.frame(rfperf$overall)
  overall$names<-rownames(overall) 
  overall$performance<-overall$rfperf.overall
  overall$rfperf.overall<-NULL
  
  byClass<-data.frame(rfperf$byClass)
  byClass$names<-rownames(byClass)
  byClass$performance<-byClass$rfperf.byClass
  byClass$rfperf.byClass<-NULL
  
  aucci<-data.frame(rocrf$ci)
  aucci$names<-c("aucciLower","auccimedian","aucciupper")
  aucci$performance<-aucci$rocrf.ci
  aucci$rocrf.ci<-NULL
  
  auc<-data.frame(rocrf$auc)
  auc$names<-c("auc")
  auc$performance<-auc$rocrf.auc
  auc$rocrf.auc<-NULL
  
  perRF<-rbind.data.frame(overall,byClass,auc,aucci)
  perRF$performance<-round(perRF$performance,5)
  
  perRFlist[[i]] <- perRF
  
}

per_acne <- do.call(cbind,perRFlist)
per_acne<-per_acne[,c(1,2,4,6,8,10,12,14,16,18)]
per_acne<-data.frame(t(per_acne))
per_acne<-per_acne[2:10,]
per_acne$AUC<- per_acne$X1
per_acne$AUCLower<- per_acne$X11
per_acne$AUCUpper<- per_acne$X3
per_acne$X1<-NULL
per_acne$X11<-NULL
per_acne$X2<-NULL
per_acne$X3<-NULL
per_acne1<-as.data.frame(lapply(per_acne,as.numeric))
per_acne1$Model<-c("RF1","RF2","RF3","RF4","RF5","RF6","RF7","RF8","RF9")

optimalRF<-per_acne1$Model[which(per_acne1$AUC==max(per_acne1$AUC))] # selecting the optimal RF using largest AUC
optimalRF

write.csv(per_acne1, 'RF_optimize/RFperformance_Acne.csv')


################################# feature importance  ##############################################
optimalRF
besti <- 7
rfbest <- readRDS(paste("RF_optimize/acne_tuneRF.",besti,".rds",sep = ""))

acnerf_importance_scale<-data.frame(importance(rfbest,scale = TRUE), check.names = FALSE)
acnerf_importance<-data.frame(importance(rfbest), check.names = FALSE)
acnerf_importance_sort<- data.frame(acnerf_importance[order(acnerf_importance$MeanDecreaseAccuracy,decreasing = TRUE),])
acnerf_importance_sort$features=rownames(acnerf_importance_sort)
acnerf_importance_sort$relativeinf<-sprintf("%0.2f", acnerf_importance_sort$MeanDecreaseAccuracy/sum(acnerf_importance_sort$MeanDecreaseAccuracy)*100)

acnerf_importance_sort$group<-ifelse(acnerf_importance_sort$features%in%c("Gender","Age","Ethnicity","Household_income","Grade"),"DSCs",
                                     ifelse(acnerf_importance_sort$features%in%c("Skin_type","Sensitive_skin","Genetic","Menstrual_blood_volume",
                                                                                 "Menstrual_colic","Menstrual_cycle","","PHQ9","BMI"),"BMAs",
                                            ifelse(acnerf_importance_sort$features%in%c("Den_barbecue_shop","Den_busstop","Den_cafeterias_shop","Den_fruit_shop",
                                                                                        "Den_hotpot_restaurant","Den_intersection","Den_KFC_McDold.s",
                                                                                        "Den_milktea_shop","Den_roadlength","Denpop","NDVI"),"BEs",
                                                   ifelse(acnerf_importance_sort$features%in%c("Daylight","Temperature","Relative_humidity",
                                                                                               "CO","NO2","O3","PM2.5","SO2"),"NEs","LFs"))))

write.csv(acnerf_importance_sort,'RF_optimize/acnerf_importance.csv')