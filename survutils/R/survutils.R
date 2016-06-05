#' Predict on survival data
#' @param filename Path to survival data file. Data should have the features in the format (time, status and other feature columns)
#' @param type Type of function to be applied. Values are {"RSF", "SUPERPC", "BOOSTCI", "SCORE", "TREC", "REC"}
#' @param rhor Row regularization parameter. Default value is 1. Used only in type "TREC"
#' @param rhoc Column regularization parameter. Default value is 1. Used only in type "TREC"
#' @param qr Norm on the row. When 1 it is L_1 norm and when 2 it is L_2 norm. This value should be kept <=2 for reducing computational cost. Used only in type "TREC"
#' @param qc Norm on the column. When 1 it is L_1 norm and when 2 it is L_2 norm. This value should be kept <=2 for reducing computational cost. Used only in type "TREC"
#' @param tol Convergence tolerance. Default value is 1e-3. Used only in type "TREC"
#' @param ntrees Number of trees for BOOSTCI model. Default is 1000
#' @param nfolds Number of folds for BOOSTCI cross validation. Default is 5
#' @param ntimes Number of distinct times of events to be considered for SCORE. Computation gets slow with large values. If no value is specified all times are taken
#' @param optimizer Type of optimization to apply for SCORE. Values are {"GLMNET", "MRCE"}. Default is "GLMNET". "MRCE" cannot run for more than 1000 features.
#' @param nfeatures Number of features for SCORE. If no value is specified all features are taken
#' @export

survutils <- function(filename, type, rhor=1,rhoc=1, qr=2, qc=2, tol = 1e-03, ntrees=1000, nfolds=5, ntimes=0, optimizer="GLMNET", nfeatures=0)
{
  if(type == "RSF")
   return(rsf(filename))
  
  else if(type == "SUPERPC")
    return(spc(filename))
  
  else if(type == "BOOSTCI")
    return(boostci(filename, ntrees, nfolds))
  
  else if(type == "SCORE")
    return(score(filename, ntimes, optimizer, nfeatures))

  else{
    
  data <- read.csv(filename, header=TRUE, sep=",")
  data[data$event==0, 1] <- NA
  res <- NULL
  
  if(type == "TREC")
  {
    ans1 <- mean.transpose(data[,-2], tol)
    xc <- ans1$xcen
    ans2 <- cov.transpose(xc,rhor,rhoc,qr,qc)
    res <- TREC(data[,-2], ans2$Sigmahat,ans2$Deltahat,ans2$Sigmaihat,ans2$Deltaihat,ans1$M)
  }
  
  else if(type == "REC")
    res <- REC(data[,-2])
  
  else
    print("Unknown type")
  
  if(!is.null(res)){
  imputed_matrix <- res$xhat
  imputed_matrix <- data.frame(time=round(imputed_matrix[,1],2), data[,2:ncol(data)])
  return(imputed_matrix)
  }
  }
  
}