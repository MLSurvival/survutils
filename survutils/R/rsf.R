#' @export

rsf <- function(filename, nfolds=5){
  
  data <- read.csv(filename, sep=",", head=TRUE)
  
  nr <- nrow(data)
  test.size <- round(nr*0.2)
  ind.cens <- sample(which(data$event==0),(test.size/2))
  ind.uncens <- sample(which(data$event==1),(test.size/2))
  
  # Create test set with equal no of censored and uncensored items
  data.test <- data[c(ind.cens, ind.uncens),]
  data.train <- data[-c(ind.cens, ind.uncens),]
  
  foldid <- coxsplit(data.train, nfolds)
  
  # N-fold cross validation creates n models
  error.rates <- sapply(c(1:nfolds), function(n){
    
    data.fold.train <- data.train[which(foldid!=n),]
    data.fold.test <- data.train[which(foldid==n),]
    model.obj <- rfsrc(Surv(time, event) ~ ., data=data.fold.train, nsplit = 10, tree.err=FALSE)
    pred <- predict.rfsrc(model.obj, data.fold.test)
    return(mean(pred$err.rate, na.rm=T))
    
  })
  
  # Find the best model based on least error rate
  bestmodel.id <- which.min(error.rates)
  data.fold <- data.train[which(foldid!=bestmodel.id),]
  bestmodel <- rfsrc(Surv(time, event) ~ ., data=data.fold, nsplit = 10, tree.err=FALSE)
  
  # Test model on test data
  
  time <- as.vector(data.test$time)
  status <- as.vector(data.test$event)
  
  preds <- predictSurvProb(bestmodel, data.test, times=c(1000))
  
  event_prediction <- 1-preds
  
  surv <- Surv(time, status)
  
  auc <- survConcordance(surv ~ event_prediction)$concordance
  
  brier_score = sum((status-event_prediction)^2)/nrow(data.test)
  
  return(c(auc=auc, brier.score=brier_score))
  
}
