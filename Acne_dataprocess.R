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
datacne <-read.csv("data/Dataacne.csv")
datacne<-datacne[,6:75]

datacne$Having_children<-NULL
datacne$Marital_status <-NULL
datacne$Height<-NULL
datacne$Weight<-NULL

datacne1<-datacne
datacne1$Waist<-NULL
datacne1$Hip<-NULL
datacne1$hukou<-NULL

############################### function to change label to Y/N#####################################
datacne2<-datacne1
set.seed(100)
index1=createDataPartition(datacne2$Acne, p = 0.8, list = FALSE)
acne_train=datacne2[index1,]
acne_test=datacne2[-index1,]

library(car)
lm<-lm(acne_train$Acne~., data=acne_train)
summary(vif(lm))
acne_train$PM10<-NULL
acne_test$PM10<-NULL 

############################### Lasso#####################################
if(T){
  library(glmnet)
  library(plotmo)
  tmp_y <- acne_train$Acne
  tmp_y <- factor(tmp_y, levels=c('0', '1'))
  cox_model <- cv.glmnet(as.matrix(acne_train[,-1]), tmp_y, family = 'binomial')   
  if(T){
    dir.create(paste0("Lasso"))
    pdf('Lasso/lasso_coef1.pdf')
    plot_glmnet(cox_model$glmnet.fit, s=cox_model$lambda.min) #, col= brewer.pal(dim(train_SF_X)[2], 'Dark2'))
    library(ggpubr)
    tmp <- data.frame(coef(cox_model, s='lambda.min')[-1,])
    colnames(tmp) <- c('Coefficient')
    tmp$Feature <- rownames(tmp)
    tmp <- tmp[order(tmp$Coefficient), ]
    tmp$color <- 'redundant'
    tmp$color[tmp$Coefficient > 0] <- 'positive'
    tmp$color[tmp$Coefficient < 0] <- 'negative'
    tmp$Coefficient <- round(tmp$Coefficient, 4)
    f <- ggbarplot(tmp, 'Feature', 'Coefficient', color='color', fill='color', palette = c('steelblue', 'red', 'gray'),
                   label=T, #repel=T,
                   ylab=paste0('Coefficient (s=',round(cox_model$lambda.min,3),')')) + coord_flip()
    print(f)
    dev.off()
  }
  write.csv(coef(cox_model, s='lambda.min')[-1,], 'Lasso/lasso_coef.csv')
  full_features <- colnames(acne_train[,-1])
  features <- coef(cox_model, s='lambda.min')[-1,]
  features <- names(features[abs(features) > 0])
  features
} 

datacne.lasso <-cbind.data.frame(datacne2[,features],Acne=datacne2$Acne)

############################### function to change label to Y/N#####################################
datacne3<-datacne2
datacne3$Gender<-change2Level(datacne3$Gender)
datacne3$Sensitive_skin<-change2Level(datacne3$Sensitive_skin)
datacne3$Ethnicity<-change2Level(datacne3$Ethnicity)
datacne3$Caring_sunscreen<-change2Level(datacne3$Caring_sunscreen)
datacne3$Using_SCP_cleaning<-change2Level(datacne3$Using_SCP_cleaning)
datacne3$Using_SCP_moisturizin<-change2Level(datacne3$Using_SCP_moisturizin)
datacne3$Using_SCP_sunscreen<-change2Level(datacne3$Using_SCP_sunscreen)
datacne3$Using_SCP_whitening<-change2Level(datacne3$Using_SCP_whitening)
datacne3$Using_SCP_oil_control<-change2Level(datacne3$Using_SCP_oil_control)
datacne3$Using_SCP_shumin<-change2Level(datacne3$Using_SCP_shumin)
datacne3$Using_SCP_anti.aging<-change2Level(datacne3$Using_SCP_anti.aging)
datacne3$Using_SCP_cosmetics<-change2Level(datacne3$Using_SCP_cosmetics)
datacne3$Using_SCP<-change2Level(datacne3$Using_SCP)
datacne3$Menstrual_colic<-change2Level(datacne3$Menstrual_colic)
datacne3$Acne<-change2Level(datacne2$Acne)

acne_num_train=datacne3[index1,]
acne_num_test=datacne3[-index1,]
saveRDS(acne_num_train,'data/acne_num_train.rds')
saveRDS(acne_num_test,'data/acne_num_test.rds')