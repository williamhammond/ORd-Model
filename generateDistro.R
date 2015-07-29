#Author: William Hammond
getDistribution <- function(x){
  GNa  <- 75     #1
  gNal <- 0.0075 #2
  Gto  <- 0.02   #3
  PCa  <- 0.0001 #4
  GKr  <- 0.046  #5
  GKs  <- 0.0034 #6
  GK1  <- 0.1908 #7
  Gncx <- 0.0008 #8
  Pnak <- 30     #9
  bt   <-4.7    #10
  input_variables <- c(GNa,gNal,Gto,PCa,GKr,GKs,GK1,Gncx,Pnak,bt)
  distributions <- 
    matrix( rnorm(x*10,mean=input_variables,sd=input_variables*.3), 10, x) 
}