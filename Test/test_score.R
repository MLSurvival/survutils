library(survutils)

filename = "synthetic1-survival.csv"
#filename = "breast.csv"

res <- survutils(filename, "SCORE")
#res <- survutils(filename, "SCORE", optimizer="MRCE")

print(res)