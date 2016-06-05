#' @export

spc <- function(filename, nfolds=5)
{
  
  data <- read.csv(filename ,header = TRUE,sep = ",")
  ncols <- ncol(data)
  nrows <- nrow(data)
  
  auc <- c()
  brier <- c()
  
  for(i in 1:10)
  {
    train.indices <- sample(1:nrows, .8*nrows)

    train.X <- as.matrix(t(data[train.indices,3:ncols]))
    train.time <- as.vector(data[train.indices,]$time)
    train.censor <- as.vector(data[train.indices,]$event)
    rownames(train.X) <- NULL
    
    test.X <- as.matrix(t(data[-train.indices,3:ncols]))
    test.time <- as.vector(data[-train.indices,]$time)
    test.censor <- as.vector(data[-train.indices,]$event)
    rownames(test.X) <- NULL
    
    feature_names <- colnames(data)
    
    datatrain <- list(x=train.X,y=train.time,censoring.status=train.censor,featurenames=feature_names)
    datatest <- list(x=test.X,y=test.time,censoring.status=test.censor,featurenames=feature_names)
    a <- superpc.train(datatrain,type="survival")
    fit <- superpc.predict(a,datatrain, datatest, threshold=1.0, n.components=1, prediction.type="continuous")
    
	# Prediction and evaluation metrics
    event_prediction <- 1-fit$v.pred
    time <- test.time
    status <- test.censor
    surv <- Surv(time, status)
    auc[i] <- survConcordance(surv ~ event_prediction)$concordance
    brier[i] <- sum((status-event_prediction)^2)/ncol(test.X)
  }
  
  return(c(auc=mean(auc), auc.sd=sd(auc), brier.score=mean(brier)))

}