#Author: William Hammond
#function to return the index and value of numeric vector x
getMaximums <- function(x){
  library(zoo)
  xz <- as.zoo(x)
  rxz <- rollapply(xz, 3, function(x) which.max(x)==2)
  index(rxz)[coredata(rxz)]
}

getMinimums <- function(x){
  library(zoo)
  xz <- as.zoo(x)
  rxz <- rollapply(xz, 3, function(x) which.min(x)==2)
  index(rxz)[coredata(rxz)]
}
extractAPFeatures <- function(voltage)
{
  peaks <- getMaximums(voltage)
  T2_MV <- peaks[1]
  v_peak <- voltage[T2_MV]
  
  valleys <- getMinimums(voltage)
  v_notch <- voltage[valleys[1]] 
  T_notch <- (peaks[1]+valleys[1])
  if (v_notch < 25 || v_notch > 15)
  {
    T_notch <- -1
    v_notch <- 1
  }
  RMP <- min(voltage)
  
  v_amp <- v_peak - RMP
  
  v90 <- RMP*0.9 + v_peak*0.1;
  v50 <- RMP*0.5 + v_peak*0.5;
  
  V_triangulation <- abs(v90-v50);
  
  #Determine the value and location of maximum rate of rise
  diff_v <- diff(voltage)
  #change 0.02 to value of dt 
  V_maxvelocity  <- max(diff_v)/0.02;
  #loc_maxv <- which(diff_v==max(diff_v), arr.ind=TRUE)
  #V_maxvelocity <- voltage[max_v]/(loc_maxv*0.02);
  
  minv90 <- 180;
  T1_v90 <- 0;
  
  for (i in 1:peaks[1]){
    if(abs(v90-voltage[i] < minv90)){
      minv90 <- abs(v90-voltage[i])
      T1_v90 <- time_array[i]  
    }
  }
  
  T2_v90 <- 0
  
  min_v50 <- 180
  T2_v50 <- 0
  for (i in peaks[1]:length(voltage)){
    if(abs(v90-voltage[i]) < minv90)
    {
      minv90 <- abs(v90-voltage[i])
      T2_v90 <- time_array[i]
    }
    if(abs(v50-voltage[i] < min_v50)){
      min_v50 <- abs(v90-voltage[i]);
      T2_v50 <- time_array[i]
    }
  }
  
  v_apd <- (T2_v90-T1_v90);
  
  T_Triangluation <- abs(T2_v50-T2_v90);
  
  d <- c(v_apd,RMP,v_peak,v_amp,v_notch,V_maxvelocity,T_Triangluation,V_triangulation,v90,v50);  
}

extractIonFeatures <- function(ion)
{
  peaks <- getMaximums(ion)
  T2_MI <- peaks[1]
  I_peak <- ion[T2_MI]
  
  diastolic <- min(ion)
  
  I_amp <- I_peak - diastolic
  
  I90 <- diastolic*0.9 + I_peak*0.1;
  I50 <- diastolic*0.5 + I_peak*0.5;
  
  
  minI90 <- 180;
  T1_I90 <- 0;
  
  for (i in 1:peaks[1]){
    if(abs(I90-ion[i] < minI90)){
      minI90 <- abs(I90-ion[i])
      T1_I90 <- time_array[i]  
    }
  }
  d <- c(diastolic,I_peak,I_amp,I90,I50,T1_I90);  
}