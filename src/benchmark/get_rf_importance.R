#!/usr/bin/Rscript --vanilla

library(randomForest)

args <- commandArgs(trailingOnly=TRUE)

for (f in args) {
  write(f, stderr())
  load(f)
  ## By total mean decrease in accuracy
  ##worst <- sort(importance(rf)[,3])
  i <- importance(rf)
  worst <- sort(apply(i[,1:2], 1, max))
  writeLines(names(worst)[1])
  rm(rf)
}