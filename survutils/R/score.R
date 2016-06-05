#' @export

score <- function(filename, numtimes=0, optimizer="GLMNET", nfeatures=0)
{
  
  data <- read.csv(filename, sep=",")
  
  if(nfeatures != 0)
	cols <- nfeatures
  else 
    cols <- ncol(data)

  if((cols-2) > 1000 & optimizer=="MRCE")
  {
	warning("MRCE cannot run for more than 1000 features. Result obtained is using GLMNET")
	optimizer <- "GLMNET"
  }
  
  features <- data[,3:cols]
  features <- as.matrix(features)
  
  if(numtimes !=0 )
  {
    cutoffs <- sort(unique(data[which(data$event==1),1]))[1:numtimes]
    output <- pseudosurv(time=data$time,event=data$event, tmax=cutoffs)
  }
  else
  {
    output <- pseudosurv(time=data$time,event=data$event)
  }
  
  multi_output_survival <- as.matrix(output$pseudo)
  
  success <- 0
  if(optimizer == "MRCE")
  {
  # Error handling in MRCE
  tryCatch(
	{
		mrcefit <- mrce(Y=multi_output_survival, X=features, lam1=1e-3, lam2=1e-3, method="single", cov.tol=0.1, tol.out=1e-10)
		finalbeta <- as.matrix(mrcefit$Bhat)
		predictions <- features %*% finalbeta
		predictions <- as.matrix(predictions)
		success <- 1
	},
	error = function(e)
   {
    print(e$message)
	print("Proceeding with GLMNET")
   },
   finally = {
	if(success == 0)
		optimizer = "GLMNET"
   }
  )
  }
	
  if(optimizer == "GLMNET")
  {
  mfit <- cv.glmnet(features, multi_output_survival, family = "mgaussian", alpha=.5)
  predictions <- predict(mfit, type="response", newx = features, s = 0.01)
  predictions <- predictions[,,1]
  }
  
  time <- data$time
  status <- data$event
  surv <- Surv(time, status)
  
  auc <- sapply(1:ncol(predictions), function(i){
    event_prediction <- 1-predictions[,i]
    survConcordance(surv ~ event_prediction)$concordance
  })
  
  brier <- sapply(1:ncol(predictions), function(i){
    event_prediction <- 1-predictions[,i]
    sum((status-event_prediction)^2)/nrow(data)
  })

  return(c(auc=mean(auc), auc.sd=sd(auc), brier.score=mean(brier)))
  
}