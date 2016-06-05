library(survutils)
library(RegCox)

filename = "synthetic1-survival.csv"
#filename = "breast.csv"

rhor = 1
rhoc = 1
qr = 2
qc = 2

train_data <- survutils(filename, "TREC", rhor, rhoc, qr, qc)
#train_data <- survutils(filename, "REC")

train_time <- train_data$time
train_event <-train_data$event
x <- as.matrix(train_data[,3:ncol(train_data)])
y <- as.matrix(cbind(time=train_time,status=train_event))
t <- mean(train_data$time)

fit <- regcox(x, y, t, "Enet",nfolds=5 , alpha=0.1, stdbeta=F, standardize=F)
print("Enet Results")
print(fit)

