library(data.table)
library(nnet)
library(caret)
library(randomForest)
library(broom)

features<-fread("../LynchModel/data/A_thal_genes_PDS5_enrich.csv")
qtls<-fread("data/qtls_10cMwindow_fixed_overlap_5.5cm_march18.txt")

qtldt <- qtls[, .(qtl = if (any(type == "cis") & any(type == "trans")) "both" 
                    else if (any(type == "cis")) "cis" 
                    else "trans"), 
                by = .(txname)]
qtldt[,txname :=toupper(txname)]

features<-features[gene %in% qtldt$txname][,c(1,8:20)]
setkey(qtldt,txname)
setkey(features,gene)
dt<-qtldt[features, nomatch=0]
dt$qtl<-as.factor(dt$qtl)
dt<-na.omit(dt)

model<-multinom(qtl ~ ., data = dt)
summary(model)

tmodel<-tidy(model)
ggplot(tmodel,aes(x = term, y = estimate))+
  geom_point()+
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
  theme_minimal() +
  labs(title = "Coefficient Plot for Multinomial Logistic Regression",
       y = "Coefficient Estimate", x = "") +
  coord_flip()




###RANDOM FOREST###
set.seed(123)
trainIndex <-createDataPartition(dt$qtl, p = 0.8, list = FALSE)
trainData<- dt[trainIndex,]
testData<- dt[-trainIndex,]
rfModel<-randomForest(qtl ~ .,data = trainData, ntree = 500, mtry = 3, importance = TRUE)
##predict
predictions<-predict(rfModel,testData)

# Calculate confusion matrix and accuracy
confMatrix <- table(testData$qtl, predictions)
accuracy <- sum(diag(confMatrix)) / sum(confMatrix)
print(confMatrix)
print(accuracy)

# Calculate other performance metrics
performance <- confusionMatrix(predictions, testData$qtl)
print(performance)

