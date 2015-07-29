require(compiler)
enableJIT(3)

library(foreach)
library(parallel)
library(iterators)
library(doParallel)
library(nnet)
library(zoo)

source("Variables.R")
source("Graphs.R")
source("generateDistro.R")
source("ORD.R")

ptm <- proc.time()
prepareSimulationData()
runSimulation()
resultFileConn <- file("results.txt")
results <- do.call(rbind,features)
write.table(results, file = resultFileConn, sep = '\t', col.names = FALSE, 
            row.names = FALSE)
proc.time() - ptm