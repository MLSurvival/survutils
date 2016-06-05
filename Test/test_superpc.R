library(survutils)

filename = "synthetic1-survival.csv"
#filename = "breast.csv"

res <- survutils(filename, "SUPERPC")
print(res)