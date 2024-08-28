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
####################################################################
acnerank <-read.csv("data/Dataacnerank.csv")
acnerank<-acnerank[,6:75]

acnerank$Having_children<-NULL
acnerank$Marital_status <-NULL
acnerank$Height<-NULL
acnerank$Weight<-NULL
skim(acnerank)

acnerank1<-acnerank
acnerank1$Waist<-NULL
acnerank1$Hip<-NULL
acnerank1$hukou<-NULL
acnerank1$Acnerank<-ifelse(acnerank1$Acne_rank<=3,0,1)

############################### function to change label to Y/N#####################################
acnerank2<-acnerank1
set.seed(100)
index1=createDataPartition(acnerank2$Acnerank, p = 0.8, list = FALSE)
acnerank_train=acnerank2[index1,]
acnerank_test=acnerank2[-index1,]
library(car)
lm<-lm(acnerank_train$Acnerank~., data=acnerank_train)
summary(vif(lm))
acnerank_train$PM10<-NULL
acnerank_test$PM10<-NULL 

############################### Lasso#####################################
acnerank3<-acnerank2
acnerank3$PM10<-NULL  ##according to VIF over 10
acnerank3$Living_expenses <-NULL   ##according lasso result 
acnerank3$IR_saltyfood<-NULL ##according lasso result 


############################### function to change label to Y/N#####################################

acnerank3$Gender<-change2Level(acnerank3$Gender)
acnerank3$Sensitive_skin<-change2Level(acnerank3$Sensitive_skin)
acnerank3$Ethnicity<-change2Level(acnerank3$Ethnicity)
acnerank3$Caring_sunscreen<-change2Level(acnerank3$Caring_sunscreen)
acnerank3$Using_SCP_cleaning<-change2Level(acnerank3$Using_SCP_cleaning)
acnerank3$Using_SCP_moisturizin<-change2Level(acnerank3$Using_SCP_moisturizin)
acnerank3$Using_SCP_sunscreen<-change2Level(acnerank3$Using_SCP_sunscreen)
acnerank3$Using_SCP_whitening<-change2Level(acnerank3$Using_SCP_whitening)
acnerank3$Using_SCP_oil_control<-change2Level(acnerank3$Using_SCP_oil_control)
acnerank3$Using_SCP_shumin<-change2Level(acnerank3$Using_SCP_shumin)
acnerank3$Using_SCP_anti.aging<-change2Level(acnerank3$Using_SCP_anti.aging)
acnerank3$Using_SCP_cosmetics<-change2Level(acnerank3$Using_SCP_cosmetics)
acnerank3$Using_SCP<-change2Level(acnerank3$Using_SCP)
acnerank3$Menstrual_colic<-change2Level(acnerank3$Menstrual_colic)
acnerank3$Acnerank<-change2Level(acnerank2$Acnerank)



acnerank_num_train=acnerank3[index1,]
acnerank_num_test=acnerank3[-index1,]
saveRDS(acnerank_num_train,'data/acnerank_num_train.rds')
saveRDS(acnerank_num_test,'data/acnerank_num_test.rds')