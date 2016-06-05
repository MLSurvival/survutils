library(survutils)

filename = "synthetic1-survival.csv"
#filename = "breast.csv"

res <- survutils(filename, "BOOSTCI", ntrees=100, nfolds=3)
print(res)