#' @export

boostci <- function(filename, ntrees=1000, nfolds=5)
{
  
  data <- read.csv(filename ,header = TRUE,sep = ",")
  
  nrows <- nrow(data)
  
  auc <- c()
  brier <- c()
  
  for(i in 1:10)
  {
    train.indices <- sample(1:nrows, .8*nrows)
    
    data.train <- data[train.indices,]
    data.test <- data[-train.indices,]
    
    gbm1 <- gbm(Surv(data.train$time, data.train$event) ~ .,
                data=data.train,
                distribution="coxph",
                n.trees=ntrees,         
                shrinkage=0.001,        
                interaction.depth=3,    
                bag.fraction = 0.5,     
                train.fraction = 0.5,   
                cv.folds = nfolds,       
                n.minobsinnode = 10,    
                keep.data = TRUE)
    
    best.iter <- gbm.perf(gbm1,method="cv")
    
	# Prediction and evaluation metrics
    event_prediction <- predict.gbm(gbm1,data.test,best.iter, type="response")
    
    surv <- Surv(data.test$time, data.test$event)
    
    auc[i] <- survConcordance(surv ~ event_prediction)$concordance
    
    brier[i] <- sum((data.test$event-event_prediction)^2)/nrow(data.test)
     
  }
  
  return(c(auc=mean(auc), auc.sd=sd(auc), brier.score=mean(brier)))

}