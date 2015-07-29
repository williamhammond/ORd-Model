# Author: William Hammond; whammondse@gmail.com
# File: ORD.R
# Description: The ORd model is described in the article "Simulation of the Undiseased
#   Human Cardia#Ventricular Action Potential: Model Formulation and
#   Experimental Validation"
#   by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy
#   This implementation is based off of work done by Mohamed Elshrif


library(zoo)
library(compiler)
library(parallel)
library(iterators)
library(doParallel)

source('generateDistro.R')
source('Graphs.R')
source('Variables.R')

prepareSimulationData <- function(){
  features <- list()
  #Initial Conditions for State Variables
  for (i in 1:nx)
    #for(i=1:nx, .combine=rbind)
  {
    v[i]          <<- -85.68
    vt[i]         <<- -85.68
    Nai[i]        <<- 7.753
    Nass[i]       <<- 7.754
    Ki[i]         <<- 143.0
    Kss[i]        <<- 143.0
    Cai[i]        <<- 0.0001177
    Cass[i]       <<- 0.0001207
    Cansr[i]      <<- 2.526
    Cajsr[i]      <<- 1.791
    xm[i]         <<- 0.007779
    xhf[i]        <<- 0.6779
    xhs[i]        <<- 0.6778
    xj[i]         <<- 0.6748
    hCaMK_slow[i] <<- 0.4311
    jCaMK[i]      <<- 0.6647
    xmL[i]        <<- 0.0002099
    xhL[i]        <<- 0.3122
    xhLCaMK[i]    <<- 0.1317
    xa[i]         <<- 0.001040
    i_fast[i]     <<- 0.9995
    i_slow[i]     <<- 0.2120
    a_CaMK[i]     <<- 0.0005301
    iCaMK_fast[i] <<- 0.9995
    iCaMK_slow[i] <<- 0.2335
    xd[i]         <<- 2.680e-9
    xff[i]        <<- 1.0
    xfs[i]        <<- 0.6514
    xfcaf[i]      <<- 1.0
    xfcas[i]      <<- 0.9463
    xjca[i]       <<- 0.9666
    xnca[i]       <<- 0.01084
    xffp[i]       <<- 1.0
    xfcafp[i]     <<- 1.0
    xrf[i]        <<- 0.01395
    xrs[i]        <<- 0.8062
    xs1[i]        <<- 0.5877
    xs2[i]        <<- 0.0005728
    xk1[i]        <<- 0.9969
    Jrel_NP[i]    <<- 1.479e-6
    Jrel_CaMK[i]  <<- 1.841e-6
    CaMK_trap[i]  <<- 0.09942
    JnakNa[i]     <<- 0.1279
    JNakK[i]      <<- -0.008524
  }
  
  #Parameters
  # Initalize parameters with the proper distribution size
  #Numerical Parameters
  #CL  <<-  1000 passing length in ms
  dt      <<- 0.02
  endtime <<- 1000.0
  #      nend <<- round(endtime / dt)
  
  icl <<- 1
  cyclelength <<- cls[icl]
  icycle <<- round(cyclelength / dt)
  #   do whole cycles, don't stop in the middle of one
  ncycles <<- round (endtime  /  cyclelength)
  
  if(ncycles < (endtime  /  cyclelength))
  {
    ncycles <<- ncycles + 1
  }
  nend<<-icycle  *  ncycles
  
  #write(6, * ) 'write n cycles of length x', ncycles,cyclelength
  
  outputevery <<- 0.5
  iout <<- round(outputevery / dt)
  
  FoRT <<- F / (R * T)
  RToF <<- 1 / FoRT
  
  #Make this an enum, endo  <-  0, epi  <-  1, M  <-  2?
  #Check all these in the same if? 
  
}

runSimulation <- function(distrubtion){
  GNa  <- distrubtion [1]
  gNal <- distrubtion [2]
  Gto  <- distrubtion [3]
  PCa  <- distrubtion [4]
  GKr  <- distrubtion [5]
  GKs  <- distrubtion [6]
  GK1  <- distrubtion [7]
  Gncx <- distrubtion [8]
  Pnak <- distrubtion [9]
  bt   <- distrubtion [10]
  
  icelltype  <-  1 
  if (icelltype == 1)
  {
    gNal <- gNal * 0.6
  }
  if (icelltype !=  0 )
  {
    Gto <- Gto  *  4.0 
  }
  Aff <- 0.6
  Afs <- 1 - Aff
  
  if (icelltype  ==  1)
  {
    PCa <- PCa  *  1.2
  }else if (icelltype  ==  1){
    PCa <- PCa  *  2.5 }
  
  PCap   <- 1.1       *  PCa
  PCaNa  <- 0.00125   *  PCa
  PCaK   <- 3.574e-4  *  PCa
  PCaNap <- 0.00125   *  PCap
  PCaKp  <-3.574e-4   *  PCap
  
  if(icelltype  ==  1)
  {
    GKr <- GKr  *  1.3  
  } else if (icelltype  ==  2){
    GKr <- GKr  *  0.8  }
  
  if(icelltype  ==  1)
  {
    GKs <- GKs  *  1.4  
  }
  
  if(icelltype  ==  1)
  {
    GK1 <- GK1  *  1.2 
  } else if(icelltype  ==  2){
    GK1 <- GK1  *  1.3  }
  
  if(icelltype  ==  1)
  {
    Gncx <- Gncx  *  1.1  
  } else if(icelltype  ==  2)  {
    Gncx <- Gncx * 1.4  }
  k3p <- 1899.0
  k4p <- 639.0
  
  if(icelltype  ==  1)
  {
    Pnak <- Pnak * 0.9  
  } else if( icelltype  ==  2 ){
    Pnak <- Pnak * 0.7 }
  
  GKb <- 0.003
  if(icelltype  ==  1)
  {
    GKb <- GKb * 0.6  
  }
  
  #Calcium Buffer Constant
  cmdnmax <-  0.05
  if(icelltype  ==  1)
  {
    cmdnmax <- cmdnmax  *  1.3  
  }
  # tables setup
  # nvt was being converted to a float 
  dvt <- (vhi-vlo) / nvt
  
  #Check the indexing on this 
  #ii was being changed to a float
  for (ii in 1:nvt)
  {
    vv <- vlo + ii * dvt
    if(abs(vv) < 1e-10)
    {
      vv <- 0.01
      vffrt[ii] <- vv * F * FoRT
      vfrt[ii]  <- vv * FoRT  
    }else{
      vffrt[ii] <- vv * F * FoRT
      vfrt[ii]  <- vv * FoRT   }
    mss[ii] <- 1.0 / (1.0 + exp((-(vv + 39.57)) / 9.871))
    #Equation (4)
    tm[ii] <- 1.0 / (6.765 * exp((vv + 11.64) / 34.77)
                     + 8.552 * exp(-(vv + 77.42) / 5.955))
    exptaumt[ii] <- exp(-dt / tm[ii])
    # Equation (5)
    hss[ii] <- 1.0 / (1 + exp((vv + 82.90) / 6.086))
    # Equation (7)
    thf[ii] <- 1.0 / (1.432e-5 * exp(-(vv + 1.196) / 6.285)
                      + 6.149 * exp((vv + 0.5096) / 20.27))
    exptauht_fast[ii] <- exp(-dt / thf[ii])
    # Equation (8)
    ths[ii] <- 1.0 / (0.009794 * exp(-(vv + 17.95) / 28.05)
                      + 0.3343 * exp((vv + 5.730) / 56.66))
    exptauht_slow[ii] <- exp(-dt / ths[ii])
    # Equation (9)
    jss[ii] <- hss[ii]
    # Equation (13)
    tj[ii] <- 2.038 + 1.0 / (0.02136 * exp(-(vv + 100.6) / 8.281)
                             + 0.3052 * exp((vv + 0.9941) / 38.45))
    exptaujt[ii] <- exp(-dt / tj[ii])
    # Equation (14)
    # Check this 
    hssp[ii] <- 1.0 / (1 + exp((vv + 89.1) / 6.086))
    #[ii] <- 1.0 / (1 + exp((vv + 89.1) / 6.086)) #hssp-114
    # Equation (16)
    thsp[ii] <- 3.0 * ths[ii]
    #            tauhCaMK_slow <- 3.0 * tauht_slow !thsp-115
    exptauhCaMK_slow[ii] <- exp(-dt / thsp[ii])
    # Equation (17)
    #            jCaMKinft[ii] <- 1.0 / (1 + exp((vv + 82.90) / 6.086))
    # Equation (21)
    tjp[ii] <- 1.46 * tj[ii]
    #            taujCaMK <- 1.46 * taujt !tjp-118
    exptaujCaMK[ii] <- exp(-dt / tjp[ii])
    # Equation (22)
    mLss[ii] <- 1.0 / (1.0 + exp((-(vv + 42.85)) / 5.264))
    # Equation (26)
    tmL[ii] <- 1.0 / (6.765 * exp((vv + 11.64) / 34.77)
                      + 8.552 * exp(-(vv + 77.42) / 5.955))
    exptaumLt[ii] <- exp(-dt / tmL[ii])
    # Equation (27)
    hLinft[ii] <- 1.0 / (1.0 + exp((vv + 87.61) / 7.488))
    # Equation (29)
    #            thL <- 200.0
    exptauhLt[ii] <- exp(-dt / thL)
    # Equation (c2)
    hLCaMKinft[ii] <- 1.0 / (1.0 + exp((vv + 93.81) / 7.488)) #hLssp-131
    # Equation (31)
    tauhLCaMKt <- 3.0 * thL #thLp-132
    exptauhLCaMK[ii] <- exp(-dt / tauhLCaMKt)
    # Equation (32)
    ainft[ii] <- 1.0 / (1.0 + exp((-(vv-14.34)) / 14.82))
    # Equation (39)
    tauat <- 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(vv-18.4099) / 29.3814)))
                       + 3.5 / (1.0 + exp((vv + 100.0) / 29.3814)))
    exptauat[ii] <- exp(-dt / tauat)
    # Equation (40)
    iinft[ii] <- 1.0 / (1.0 + exp((vv + 43.94) / 5.711))
    # Equation (42)
    
    if (icelltype  <-  1)
    {
      delta_epi <- 1.0-(0.95 / (1.0 + exp((vv + 70.0) / 5.0))) 
    }else{
      delta_epi <- 1.0 }
    
    tauit_fast <- 4.562 + 1 / (0.3933 * exp((-(vv + 100.0)) / 100.0)
                               + 0.08004 * exp((vv + 50.0) / 16.59))
    tauit_fast <- tauit_fast * delta_epi
    exptauit_fast[ii] <- exp(-dt / tauit_fast) #tiF-151
    # Equation (43)
    
    tauit_slow <- 23.62 + 1 / (0.001416 * exp((-(vv + 96.52)) / 59.05)
                               + 1.780e-8 * exp((vv + 114.1) / 8.079)) #tiS-152
    tauit_slow <- tauit_slow * delta_epi
    exptauit_slow[ii] <- exp(-dt / tauit_slow)
    # Equation (44)
    Ai_fast[ii] <- 1.0 / (1.0 + exp((vv-213.6) / 151.2))
    # Equation (45)
    Ai_slow[ii] <- 1.0-Ai_fast[ii]
    # Equation (46)
    aCaMKinft[ii] <- 1.0 / (1.0 + exp((-(vv-24.34)) / 14.82)) #assp-160
    # Equation (50)
    tauaCaMKt <- 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(vv-18.4099) / 29.3814)))
                           + 3.5 / (1.0 + exp((vv + 100.0) / 29.3814)))
    exptauaCaMKt[ii] <- exp(-dt / tauaCaMKt)
    # Equation (51)
    iCaMKinft[ii] <- 1.0 / (1.0 + exp((vv + 43.94) / 5.711))
    # Equation (53)
    DeltaCaMK_develop[ii] <- 1.354 + 1.0e-4 / (exp((vv-167.4) / 15.89)
                                               + exp(-(vv-12.23) / 0.2154))
    # Equation (54)
    DeltaCaMK_recover[ii] <- 1.0-0.5 / (1.0 + exp((vv + 70.0) / 20.0))
    # Equation (55)
    tauiCaMK_fast <- tauit_fast * DeltaCaMK_develop[ii] * DeltaCaMK_recover[ii]
    exptauiCaMK_fast[ii] <- exp(-dt / tauiCaMK_fast)
    # Equation (56)
    tauiCaMK_slow <- tauit_slow * DeltaCaMK_develop[ii]* DeltaCaMK_recover[ii]
    exptauiCaMK_slow[ii] <- exp(-dt / tauiCaMK_slow)
    # Equation (57)
    AiCaMK_fast[ii] <- Ai_fast[ii]
    # Equation (58)
    AiCaMK_slow[ii] <- Ai_slow[ii]
    # Equation (59)
    dinft[ii] <- 1.0 / (1.0 + exp((-(vv + 3.940)) / 4.230))
    # Equation (65)
    taudt <- 0.6 + 1.0 / (exp(-0.05 * (vv + 6.0)) + exp(0.09 * (vv + 14.0)))
    exptaudt[ii] <- exp(-dt / taudt)
    # Equation (66)
    finft[ii] <- 1.0 / (1.0 + exp((vv + 19.58) / 3.696))
    # Equation (68)
    tauft_fast <- 7.0 + 1.0 / (0.0045 * exp(-(vv + 20.0) / 10.0)
                               + 0.0045 * exp((vv + 20.0) / 10.0))
    exptauft_fast[ii] <- exp(-dt / tauft_fast)
    # Equation (69)
    tauft_slow <- 1000.0 + 1.0 / (0.000035 * exp(-(vv + 5.0) / 4.0)
                                  + 0.000035 * exp((vv + 5.0) / 6.0))
    exptauft_slow[ii] <- exp(-dt / tauft_slow)
    # Equation (70)
    fCainft[ii] <- finft[ii]
    # Equation (74)
    taufCat_fast <- 7.0 + 1.0 / (0.04 * exp(-(vv-4.0) / 7.0)
                                 + 0.04 * exp((vv-4.0) / 7.0))
    exptaufCat_fast[ii] <- exp(-dt / taufCat_fast)
    # Equation (75)
    taufCat_slow <- 100.0 + 1.0 / (0.00012 * exp(-vv / 3.0)
                                   + 0.00012 * exp(vv / 7.0))
    exptaufCat_slow[ii] <- exp(-dt / taufCat_slow)
    # Equation (76)
    AfCa_fast[ii] <- 0.3 + 0.6 / (1.0 + exp((vv-10.0) / 10.0))
    # Equation (77)
    AfCa_slow[ii] <- 1.0-AfCa_fast[ii]
    # Equation (78)
    jCainft[ii] <- finft[ii]
    # Equation (82)
    taujCat <- 75.0
    exptaujCat[ii] <- exp(-dt / taujCat)
    # Equation (82-83)
    fCaMKinft[ii] <- finft[ii]
    # Equation (84)
    taufCaMKt_fast <- 2.5 * tauft_fast
    exptaufCaMKt_fast[ii] <- exp(-dt / taufCaMKt_fast)
    # Equation (85)
    fCa_CaMKinft[ii] <- finft[ii]
    # Equation (89)
    taufCa_CaMKt_fast <- 2.5 * taufCat_fast
    exptaufCa_CaMKt_fast[ii] <- exp(-dt / taufCa_CaMKt_fast)
    # Equation (90)
    AfCa_CaMK_fast[ii] <- AfCa_fast[ii]
    # Equation (91)
    AfCa_CaMK_slow[ii] <- 1.0-AfCa_CaMK_fast[ii]
    # Equation (92)
    Xrinft[ii] <- 1.0 / (1.0 + exp((-(vv + 8.337)) / 6.789))
    # Equation (111)
    tauXrt_fast <- 12.98 + 1.0 / (0.3652 * exp((vv-31.66) / 3.869)
                                  + 4.123e-5 * exp((-(vv-47.78)) / 20.38))
    exptauXrt_fast[ii] <- exp(-dt / tauXrt_fast)
    # Equation (112)
    tauXrt_slow <- 1.865 + 1.0 / (0.06629 * exp((vv-34.70) / 7.355)
                                  + 1.128e-5 * exp((-(vv-29.74)) / 25.94))
    exptauXrt_slow[ii] <- exp(-dt / tauXrt_slow)
    # Equation (113)
    AXrt_fast[ii] <- 1.0 / (1.0 + exp((vv + 54.81) / 38.21))
    # Equation (114)
    AXrt_slow[ii] <- 1.0 - AXrt_fast[ii]
    # Equation (115)
    RKr[ii] <- 1.0 / (1.0 + exp((vv + 55.0) / 75.0)) * 1.0 /
      (1.0 + exp((vv-10.0) / 30.0))
    # Equation (119)
    Xs1inft[ii] <- 1.0 / (1.0 + exp((-(vv + 11.60)) / 8.932))
    # Equation (121)
    tauXs1t <- 817.3 + 1.0 / (2.326e-4 * exp((vv + 48.28) / 17.80)
                              +  0.001292 * exp((-(vv + 210.0)) / 230.0))
    exptauXs1t[ii] <- exp(-dt / tauXs1t)
    # Equation (122)
    Xs2inft[ii] <- Xs1inft[ii]
    # Equation (124)
    tauXs2t <- 1.0 / (0.01 * exp((vv-50.0) / 20.0)
                      +  0.0193 * exp((-(vv + 66.54)) / 31.0))
    exptauXs2t[ii] <- exp(-dt / tauXs2t)
    # Equation (125)
    Xk1inft[ii] <- 1.0 / (1.0 + exp(-(vv + 2.5538 * ko + 144.59)
                                    / (1.5692 * ko + 3.8115)))
    # Equation (128)
    tauXk1t <- 122.2 / (exp((-(vv + 127.2)) / 20.36)
                        +  exp((vv + 236.8) / 69.33))
    exptauXk1t[ii] <- exp(-dt / tauXk1t)
    # Equation (129)
    RK1[ii] <- 1.0 / (1.0 + exp((vv + 105.8-2.6 * ko) / 9.493))
    # Equation (131)
    hCa[ii] <- exp(qca * vv * FoRT)
    # Equation (133)
    hNa[ii] <- exp(qna * vv * FoRT)
    # Equation (134)
    KNai[ii] <- KNai0 * exp(delta * vv * FoRT / 3)
    # Equation (173)
    Knao[ii] <- Knao0 * exp((1.0-delta) * vv * FoRT / 3)
    # Equation (174)
    X_Kb[ii] <- 1.0 / (1.0 + exp(-(vv-14.48) / 18.34))
    # Equation (197) 
    
  }
  caiC     <- 1
  count    <- 1
  initoutput   <- "fort."
  initvoltfile <- "Voltage."
  #      current2 <- 'Incx.xxxx'
  CaMKa = CaMK_trap[1]
  # file.create(voltfile)
  # voltFileConn <- file(voltfile)
  # ============ TIME LOOP ============ 
  #Start main loop  
  
  for (icl in 1:1) #ncl
  {
    intcl <- round(cls[icl])
    #write(6,*) 'intcl <-  ' ,intcl
    print (intcl)
    idig1 <- intcl / 1000
    istp  <- intcl %% 1000
    idig2 <- istp / 100
    istp  <- istp %% 100
    idig3 <- istp / 10
    idig4 <- istp %% 10
    #     voltFile <- paste("Voltage.",dCount, sep='')
    #     voltFileConn <- file(voltFile)
    #     output   <- paste("Fort.",dCount, sep='')
    #     outputFileConn <- file(output)
    
    cyclelength <- cls[icl]
    icycle      <- round(cyclelength / dt)
    ncycles     <- round(endtime / cyclelength)
    if (ncycles > endtime / cyclelength)
    {
      ncycles <- ncycles+1 
    }
    
    nend <- icycle*ncycles
    #nend may need to be multiplied by NCL later
    time_array <<-numeric(nend)
    v_array <<-numeric(nend)
    cai_array <<-numeric(nend)
    # nend <- round(endtime/dt)
    cat("write n cycles of length x\n",ncycles,cyclelength)
    
    nups1   <- 0
    ndowns1 <- 0
    ntime   <- 0
    time    <- 0.0
    
    #start time loop
    for (ntime in 1:nend)
    {
      ntime <- ntime + 1
      time  <- ntime * dt
      
      if(ntime %% 25000==1) print(time)
      
      #Start Cell Loop
      for (i in 1:nx)
      {
        
        #compute indices into tables
        zindexv <- (v[i]-vlo) / dvt
        iv1 <- round(zindexv)
        if(iv1 < 0 || iv1 > nvt)
        {
          cat("Voltage index outside voltage table: \n",ntime,ntime*dt,v[i],iv1)
          stop  
        }
        
        #reversal potentials
        ENa <- RToF*log(nao/Nai[i])
        EK  <- RToF*log(ko/Ki[i])
        EKs <- RToF*log((ko+PKNa*nao)/(Ki[i]+PKNa*Nai[i]))
        
        #calculate INa
        
        xm[i]   <-  mss[iv1]-(mss[iv1]-xm[i])*exptaumt[iv1]
        
        xhf[i]  <-  hss[iv1]-(hss[iv1]-xhf[i])*exptauht_fast[iv1]
        xhs[i]  <-  hss[iv1]-(hss[iv1]-xhs[i])*exptauht_slow[iv1]
        xh <- Ahf*xhf[i]+Ahs*xhs[i]
        
        xj[i]  <-  jss[iv1]-(jss[iv1]-xj[i])*exptaujt[iv1]
        
        hCaMK_slow[i]  <-  hssp[iv1]-(hssp[iv1]-hCaMK_slow[i]) * 
          exptauhCaMK_slow[iv1]
        hp <- Ahf*xhf[i]+Ahs*hCaMK_slow[i]
        
        jCaMK[i]  <-  jss[iv1]-(jss[iv1]-jCaMK[i])*exptaujCaMK[iv1]
        fINap <- (1.0/(1.0+KmCaMK/CaMKa))
        INa <- GNa*(v[i]-ENa)*xm[i]*xm[i]*xm[i]*((1.0-fINap) *
                                                   xh*xj[i]+fINap*hp*jCaMK[i])
        
        #calculate INaL
        
        xmL[i]  <-  mLss[iv1]-(mLss[iv1]-xmL[i])*exptaumLt[iv1]
        
        xhL[i]  <-  hLinft[iv1]-(hLinft[iv1]-xhL[i])*exptauhLt[iv1]
        
        xhLCaMK[i]  <-  hLCaMKinft[iv1]-(hLCaMKinft[iv1]-xhLCaMK[i]) * 
          exptauhLCaMK[iv1]
        
        fINaLp <- (1.0/(1.0+KmCaMK/CaMKa))
        
        INaL <- gNal*(v[i]-ENa)*xmL[i]*((1.0-fINaLp)*xhL[i] + 
                                                  fINaLp*xhLCaMK[i])
        
        #calculate Ito
        
        xa[i]  <-  ainft[iv1]-(ainft[iv1]-xa[i])*exptauat[iv1]
        
        i_fast[i]  <-  iinft[iv1]-(iinft[iv1]-i_fast[i])*exptauit_fast[iv1]
        i_slow[i]  <-  iinft[iv1]-(iinft[iv1]-i_slow[i])*exptauit_slow[iv1]
        
        xi  <-  Ai_fast[iv1]*i_fast[i]+Ai_slow[iv1]*i_slow[i]
        
        a_CaMK[i] <- aCaMKinft[iv1]-(aCaMKinft[iv1]-a_CaMK[i])*
          exptauaCaMKt[iv1]
        
        iCaMK_fast[i]  <-  iCaMKinft[iv1]-(iCaMKinft[iv1]-iCaMK_fast[i])*
          exptauiCaMK_fast[iv1]
        iCaMK_slow[i]  <-  iCaMKinft[iv1]-(iCaMKinft[iv1]-iCaMK_slow[i])* 
          exptauiCaMK_slow[iv1]
        
        ip  <-  AiCaMK_fast[iv1]*iCaMK_fast[i]+AiCaMK_slow[iv1]*iCaMK_slow[i]
        
        fItop <- (1.0/(1.0+KmCaMK/CaMKa))
        Ito <- Gto*(v[i]-EK)*((1.0-fItop)*xa[i]*xi+fItop*a_CaMK[i]*ip)
        
        #calculate ICaL, ICaNa, ICaK
        
        xd[i] <- dinft[iv1]-(dinft[iv1]-xd[i])*exptaudt[iv1]
        
        xff[i] <- finft[iv1]-(finft[iv1]-xff[i])*exptauft_fast[iv1]
        xfs[i] <- finft[iv1]-(finft[iv1]-xfs[i])*exptauft_slow[iv1]
        
        xf <- Aff*xff[i]+Afs*xfs[i]
        
        xfcaf[i]  <-  fCainft[iv1]-(fCainft[iv1]-xfcaf[i]) * exptaufCat_fast[iv1]
        
        xfcas[i]  <-  fCainft[iv1]-(fCainft[iv1]-xfcas[i]) * exptaufCat_slow[iv1]
        
        fca <- AfCa_fast[iv1]*xfcaf[i]+AfCa_slow[iv1]*xfcas[i]
        
        xjca[i]  <-  jCainft[iv1]-(jCainft[iv1]-xjca[i]) * exptaujCat[iv1]
        
        xffp[i]  <-  fCaMKinft[iv1]-(fCaMKinft[iv1]-xffp[i]) * exptaufCaMKt_fast[iv1]
        
        fp <- Aff*xffp[i]+Afs*xfs[i]
        
        xfcafp[i]  <-  fCa_CaMKinft[iv1]-(fCa_CaMKinft[iv1]-xfcafp[i]) * 
          exptaufCa_CaMKt_fast[iv1]
        
        fcap <- AfCa_CaMK_fast[iv1]*xfcafp[i]+AfCa_CaMK_slow[iv1]*xfcas[i]
        
        km2n <- xjca[i]*1.0
        anca <- 1.0/(k2n/km2n+(1.0+Kmn/Cass[i])**4)
        dnca <- anca*k2n-xnca[i]*km2n
        xnca[i] <- xnca[i]+dt*dnca
        PhiCaL <- 4.0*vffrt[iv1]*(Cass[i]*exp(2.0*vfrt[iv1])-0.341*cao)/
          (exp(2.0*vfrt[iv1])-1.0)
        PhiCaNa <- 1.0*vffrt[iv1]*(0.75*Nass[i]*exp(1.0*vfrt[iv1])
                                   -0.75*nao)/(exp(1.0*vfrt[iv1])-1.0)
        PhiCaK <- 1.0*vffrt[iv1]*(0.75*Kss[i]*exp(1.0*vfrt[iv1])-0.75*ko)/
          (exp(1.0*vfrt[iv1])-1.0)
        
        
        fICaLp <- (1.0/(1.0+KmCaMK/CaMKa))
        
        ICaL <- (1.0-fICaLp)*PCa*PhiCaL*xd[i]*(xf*(1.0-xnca[i])+
                                                         xjca[i]*fca*xnca[i])
        +fICaLp*PCap*PhiCaL*xd[i]*(fp*(1.0-xnca[i])+
                                     xjca[i]*fcap*xnca[i])
        ICaNa <- (1.0-fICaLp)*PCaNa*PhiCaNa*xd[i]*(xf*(1.0-xnca[i])+
                                                     xjca[i]*fca*xnca[i])+
          fICaLp*PCaNap*PhiCaNa*xd[i]*(fp*(1.0-xnca[i])+
                                         xjca[i]*fcap*xnca[i])
        ICaK <- (1.0-fICaLp)*PCaK*PhiCaK*xd[i]*(xf*(1.0-xnca[i])+
                                                  xjca[i]*fca*xnca[i])+
          fICaLp*PCaKp*PhiCaK*xd[i]*(fp*(1.0-xnca[i])+
                                       xjca[i]*fcap*xnca[i])
        
        #calculate IKr
        
        xrf[i]  <-  Xrinft[iv1]-(Xrinft[iv1]-xrf[i])*exptauXrt_fast[iv1]
        xrs[i]  <-  Xrinft[iv1]-(Xrinft[iv1]-xrs[i])*exptauXrt_slow[iv1]
        
        xr <- AXrt_fast[iv1]*xrf[i]+AXrt_slow[iv1]*xrs[i]
        
        IKr <- GKr*sqrt(ko/5.4)*xr*RKr[iv1]*(v[i]-EK)
        
        #calculate IKs
        
        xs1[i]  <-  Xs1inft[iv1]-(Xs1inft[iv1]-xs1[i])*exptauXs1t[iv1]
        xs2[i]  <-  Xs2inft[iv1]-(Xs2inft[iv1]-xs2[i])*exptauXs2t[iv1]
        
        KsCa <- 1.0+0.6/(1.0+(3.8e-5/Cai[i])**1.4)
        
        IKs <- GKs*KsCa*xs1[i]*xs2[i]*(v[i]-EKs)
        
        # calculate IK1
        
        xk1[i]  <-  Xk1inft[iv1]-(Xk1inft[iv1]-xk1[i])*exptauXk1t[iv1]
        
        IK1 <- GK1*sqrt(ko)*RK1[iv1]*xk1[i]*(v[i]-EK)
        
        #calculate INaCa_i
        
        h1 <- 1+Nai[i]/kna3*(1+hNa[iv1])
        h2 <- (Nai[i]*hNa[iv1])/(kna3*h1)
        h3 <- 1.0/h1
        h4 <- 1.0+Nai[i]/kna1*(1+Nai[i]/kna2)
        h5 <- Nai[i]*Nai[i]/(h4*kna1*kna2)
        h6 <- 1.0/h4
        h7 <- 1.0+nao/kna3*(1.0+1.0/hNa[iv1])
        h8 <- nao/(kna3*hNa[iv1]*h7)
        h9 <- 1.0/h7
        h10 <- kasymm+1.0+nao/kna1*(1.0+nao/kna2)
        h11 <- nao*nao/(h10*kna1*kna2)
        h12 <- 1.0/h10
        k1 <- h12*cao*kcaon
        k2 <- kcaoff
        k3p <- h9*wca
        k3pp <- h8*wnaca
        k3 <- k3p+k3pp
        k4p <- h3*wca/hCa[iv1]
        k4pp <- h2*wnaca
        k4 <- k4p+k4pp
        k5 <- kcaoff
        k6 <- h6*Cai[i]*kcaon
        k7 <- h5*h2*wna
        k8 <- h8*h11*wna
        x1 <- k2*k4*(k7+k6)+k5*k7*(k2+k3)
        x2 <- k1*k7*(k4+k5)+k4*k6*(k1+k8)
        x3 <- k1*k3*(k7+k6)+k8*k6*(k2+k3)
        x4 <- k2*k8*(k4+k5)+k3*k5*(k1+k8)
        E1 <- x1/(x1+x2+x3+x4)
        E2 <- x2/(x1+x2+x3+x4)
        E3 <- x3/(x1+x2+x3+x4)
        E4 <- x4/(x1+x2+x3+x4)
        allo <- 1.0/(1.0+(KmCaAct/Cai[i])**2)
        JncxNa <- 3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp
        JncxCa <- E2*k2-E1*k1
        
        INaCa_i <- 0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa)
        
        #calculate INaCa_ss
        
        h1     <- 1+Nass[i]/kna3*(1+hNa[iv1])
        h2     <- (Nass[i]*hNa[iv1])/(kna3*h1)
        h3     <- 1.0/h1
        h4     <- 1.0+Nass[i]/kna1*(1+Nass[i]/kna2)
        h5     <- Nass[i]*Nass[i]/(h4*kna1*kna2)
        h6     <- 1.0/h4
        h7     <- 1.0+nao/kna3*(1.0+1.0/hNa[iv1])
        h8     <- nao/(kna3*hNa[iv1]*h7)
        h9     <- 1.0/h7
        h10    <- kasymm+1.0+nao/kna1*(1+nao/kna2)
        h11    <- nao*nao/(h10*kna1*kna2)
        h12    <- 1.0/h10
        k1     <- h12*cao*kcaon
        k2     <- kcaoff
        k3p    <- h9*wca
        k3pp   <- h8*wnaca
        k3     <- k3p+k3pp
        k4p    <- h3*wca/hCa[iv1]
        k4pp   <- h2*wnaca
        k4     <- k4p+k4pp
        k5     <- kcaoff
        k6     <- h6*Cass[i]*kcaon
        k7     <- h5*h2*wna
        k8     <- h8*h11*wna
        x1     <- k2*k4*(k7+k6)+k5*k7*(k2+k3)
        x2     <- k1*k7*(k4+k5)+k4*k6*(k1+k8)
        x3     <- k1*k3*(k7+k6)+k8*k6*(k2+k3)
        x4     <- k2*k8*(k4+k5)+k3*k5*(k1+k8)
        E1     <- x1/(x1+x2+x3+x4)
        E2     <- x2/(x1+x2+x3+x4)
        E3     <- x3/(x1+x2+x3+x4)
        E4     <- x4/(x1+x2+x3+x4)
        allo   <- 1.0/(1.0+(KmCaAct/Cass[i])**2)
        JncxNa <- 3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp
        JncxCa <- E2*k2-E1*k1
        
        INaCa_ss <- 0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa)
        
        #calculate INaK
        
        P <- eP/(1.0+H/Khp+Nai[i]/Knap+Ki[i]/Kxkur)
        a1 <- (k1p*(Nai[i]/KNai[iv1])**3)/
          ((1.0+Nai[i]/KNai[iv1])**3+(1.0+Ki[i]/KKi)**2-1.0)
        b1 <- k1m*MgADP
        a2 <- k2p
        b2 <- (k2m*(nao/Knao[iv1])**3)/((1.0+nao/Knao[iv1])**3+
                                          (1.0+ko/Kko)**2-1.0)
        a3 <- (k3p*(ko/Kko)**2)/((1.0+nao/Knao[iv1])**3+
                                   (1.0+ko/Kko)**2-1.0)
        b3 <- (k3m*P*H)/(1.0+MgATP/Kmgatp)
        a4 <- (k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp)
        b4 <- (k4m*(Ki[i]/KKi)**2)/((1.0+Nai[i]/KNai[iv1])**3+
                                      (1.0+Ki[i]/KKi)**2-1.0)
        x1 <- a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2
        x2 <- b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4
        x3 <- a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1
        x4 <- b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1
        E1 <- x1/(x1+x2+x3+x4)
        E2 <- x2/(x1+x2+x3+x4)
        E3 <- x3/(x1+x2+x3+x4)
        E4 <- x4/(x1+x2+x3+x4)
        JnakNa[i] <- 3.0*(E1*a3-E2*b3)
        JNakK[i] <- 2.0*(E4*b1-E3*a1)
        
        INaK <- Pnak*(zna*JnakNa[i]+zk*JNakK[i])
        
        #calculate IKb
        #        xkb <- 1.0/(1.0+exp(-(v[i]-14.48)/18.34))
        
        IKb <- GKb*X_Kb[iv1]*(v[i]-EK)
        
        #calculate INab
        INab <- PNab*vffrt[iv1]*(Nai[i]*exp(vfrt[iv1])-nao)/
          (exp(vfrt[iv1])-1.0)
        
        #calculate ICab
        ICab <- PCab*4.0*vffrt[iv1]*(Cai[i]*exp(2.0*vfrt[iv1])-0.341*cao)/
          (exp(2.0*vfrt[iv1])-1.0)
        
        #calculate IpCa
        IpCa <- GpCa*Cai[i]/(0.0005+Cai[i])
        
        #calculate the stimulus current, Istim
        if ( (ntime*dt) %% cyclelength <= duration) 
        {
          Istim <- amp 
        }else{
          Istim <- 0.0}
        
        
        
        #update the membrane voltage
        vt[i] <- v[i]-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+
                            INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+Istim)
        
        
        #add an inner time loop to update calcium more often
        #Start Calcium Loop
        for (iloop in 1:20)
        {
          #use a smaller dt for updates within this loop!
          dtinner <- dt*0.05
          
          #update CaMK
          CaMKb        <- CaMKo*(1.0-CaMK_trap[i])/(1.0+KmCaM/Cass[i])
          CaMKa        <- CaMKb+CaMK_trap[i]
          dCaMK_trap   <- aCaMK*CaMKb*(CaMKb+CaMK_trap[i])-bCaMK*CaMK_trap[i]
          CaMK_trap[i] <- CaMK_trap[i]+dtinner*dCaMK_trap
          
          #calculate diffusion fluxes
          Jdiff <- (Cass[i]-Cai[i])/0.2
          #            dJdiff[i] <- (Cass[i]-Cai[i])/0.2
          #            Jdiff[i] <- Jdiff[i]+dt*dJdiff[i]
          
          #calculate ryanodione receptor calcium induced calcium release from the jsr
          a_rel <- 0.5*bt
          Jrel_inf <- a_rel*(-ICaL)/(1.0+(1.5/Cajsr[i])**8)
          if(icelltype==2)
          {
            Jrel_inf <- Jrel_inf*1.7 
          }
          tau_rel <- bt/(1.0+0.0123/Cajsr[i])
          
          if (tau_rel > 0.001)
          {
            tau_rel <- 0.001
          }
          
          dJrel_NP <- (Jrel_inf-Jrel_NP[i])/tau_rel
          Jrel_NP[i] <- Jrel_NP[i]+dtinner*dJrel_NP
          btp <- 1.25*bt
          a_relp <- 0.5*btp
          Jrel_infp <- a_relp*(-ICaL)/(1.0+(1.5/Cajsr[i])**8)
          if(icelltype==2)
          {
            Jrel_infp <- Jrel_infp*1.7
          }
          tau_relp <- btp/(1.0+0.0123/Cajsr[i])
          
          if(tau_relp > 0.001)
          {
            tau_relp <- 0.001 
          }
          
          dJrel_CaMK   <- (Jrel_infp-Jrel_CaMK[i])/tau_relp
          Jrel_CaMK[i] <- Jrel_CaMK[i]+dtinner*dJrel_CaMK
          fJrel_CaMK   <- (1.0/(1.0+KmCaMK/CaMKa))
          Jrel         <- (1.0-fJrel_CaMK)*Jrel_NP[i]+fJrel_CaMK*Jrel_CaMK[i]
          
          # if(Jrel > -0.1) 
          # {
          #   Jrel <- -0.1 
          # }        
          
          #calculate serca pump, ca uptake flux
          Jupnp <- 0.004375*Cai[i]/(Cai[i]+0.00092)
          Jupp <- 2.75*0.004375*Cai[i]/(Cai[i]+0.00092-0.00017)
          if(icelltype==1) 
          {
            Jupnp <- Jupnp*1.3
            Jupp <- Jupp*1.3 
          }
          fJupp <- (1.0/(1.0+KmCaMK/CaMKa))
          Jleak <- 0.0039375*Cansr[i]/15.0
          Jup <- (1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak
          
          #calculate tranlocation flux
          Jtr <- (Cansr[i]-Cajsr[i])/100.0
          
          
          BCai <- 1.0/(1.0+cmdnmax*kmcmdn/(kmcmdn+Cai[i])**2+
                         trpnmax*kmtrpn/(kmtrpn+Cai[i])**2)
          
          #     Cai[i] <- Cai[i]+dt*(BCai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)
          #         -Jup*vnsr/vmyo+Jdiff*vss/vmyo))
          dCai <- BCai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-
                          ((1.0-fJupp)*Jupnp+fJupp*Jupp-(0.0039375*Cansr[i]/15.0))*
                          vnsr/vmyo+((Cass[i]-Cai[i])/0.2)*vss/vmyo)
          Cai[i] <- Cai[i]+dtinner*dCai
          
          BCass <- 1.0/(1.0+BSRmax*KmBSR/(KmBSR+Cass[i])**2
                        +BSLmax*KmBSL/(KmBSL+Cass[i])**2)
          dCass <- BCass*(-(ICaL-2.0*INaCa_ss)*Acap/
                            (2.0*F*vss)+Jrel*vjsr/vss-Jdiff)
          Cass[i] <- Cass[i]+dtinner*dCass
          
          dCansr <- Jup-Jtr*vjsr/vnsr
          Cansr[i] <- Cansr[i]+dtinner*dCansr
          
          BCajsr <- 1.0/(1.0+csqnmax*kmcsqn/(kmcsqn+Cajsr[i])**2)
          dCajsr <- BCajsr*(Jtr-Jrel)
          Cajsr[i] <- Cajsr[i]+dtinner*dCajsr
        } # end inner loop for calcium
        
        #finish by updating K and Na
        #calculate diffusion fluxes
        JdiffNa <- (Nass[i]-Nai[i])/2.0
        JdiffK <- (Kss[i]-Ki[i])/2.0
        
        #update intracellular concentrations, using buffers for Cai, Cass, Cajsr
        dNai <- -(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+
          JdiffNa*vss/vmyo 
        Nai[i] <- Nai[i]+dt*dNai
        
        dNass <- -(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa
        Nass[i] <- Nass[i]+dt*dNass
        
        dKi <- -(Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)+
          JdiffK*vss/vmyo
        Ki[i] <- Ki[i]+dt*dKi
        
        dKss <- -(ICaK)*Acap/(F*vss)-JdiffK
        Kss[i] <- Kss[i]+dt*dKss
        
      } # End of the Cell Loop 
      
      time_array[count] <<- time
      v_array[count] <<- v
      cai_array[count] <<- Cai
      count <- count + 1  
      
      voltfileDF <- data.frame(v_array,cai_array)
      
      #cat (c(paste(time,v,Cai)), file = voltFileConn,append=TRUE,sep="\n")
      
      #check for APD threshold crossings
      #APD90 threshold for 0d, 30s, CL <- 1000ms
      #Maximum of the Action Potential value  <-  40.35
      #Resting Potential value  <-  -87.94
      #%90 of the action potential will be (-87.94*0.9) + (40.35*0.1)
      #APD_90  <-  -79.146 + 4.035  <-  -75.111
      i1 <- 1
      if(v[i1] > vru && vt[i1] >= vr && nups1==ndowns1)
      {
        nups1 <- nups1+1
        ups1[nups1] <- ntime*dt+dt*(vr-v[i1])/(vt[i1]-v[i1]) 
      }else if(v[i1] >vr && vt[i1] <= vr && nups1==ndowns1+1){ 
        ndowns1 <- ndowns1+1
        downs1[ndowns1] <- ntime*dt+dt*(vr-v[i1])/(vt[i1]-v[i1]) }
      
      for (i in 1:nx)
      {
        v[i] <- vt[i]  
      }
      
    } # End Time Loop
    
    #writeLines(c(paste(time,di,apd,cl,prevdi,prevapd,prevcl)))
    
    #     cat (c(paste(v,Nai,Nass,Ki,Kss,Cai,Cass,Cansr,Cajsr,xm,xhf,
    #                  xhs,xj,hCaMK_slow,jCaMK,xmL,xhL,xhLCaMK,xa,
    #                  i_fast,i_slow,a_CaMK,iCaMK_fast,iCaMK_slow,xd,
    #                  xff,xfs,xfcaf,xfcas,xjca,xnca,xffp,xfcafp,xrf,
    #                  xrs,xs1,xs2,xk1,Jrel_NP,Jrel_CaMK,CaMK_trap,
    #                  JnakNa,JNakK)), file = output,append=TRUE,sep="\n")
  }# End Main Loop
  
  ionFeatures <<- extractIonFeatures(cai_array)
  voltageFeatures <<- extractAPFeatures(v_array)
  features <- c(ionFeatures,voltageFeatures)
  return (features)
}

executeSimulation <- function(cl,dSize){
  distros <- getDistribution(dSize)
  write.table(distros, file = "distrubtions.txt", col.names = FALSE, row.names = FALSE )
  clusterExport (cl=cl, list = c("vlo","vhi","nvt","FoRT","vffrt","vfrt","mss","tm"))
  prepareSimulationData()
  runSimulationc <- cmpfun(runSimulation)

  return (clusterApply(cl,distros,runSimulationc))
  #return (clusterApply(cl = cl, distros, runSimulationc))
  #return (parApply (cl,distros,2,runSimulationc))
}
analyzeResults <- function(){
  results <- read.table("results.txt")
  resultsHF <- read.table("results_hf.txt")
  distroHF <- read.table("normal_ions.txt")
  distro <- read.table("HF_ions.txt")
  
  features <- do.call(cbind, list(results,resultsHF))
  features[,1] <- NULL
  features[,51] <- NULL
  distroMatrix <- do.call(cbind, list(distro,distroHF))
}
errorReader <- function (dSize){
  tryCatch({
    runSimulation(dSize)
  },
  error=function(e){
    print(e)
    stop(e)
  })
} 
# cores <- detectCores()
# cl <- makePSOCKcluster(cores)
# registerDoParallel(cl)

ptm <- proc.time()
cl <- makeCluster(mc <- 1)
results <- executeSimulation(cl,4)
#results <- executeSimulation(cl,2)
#results <- executeSimulation(5)
#results <- executeSimulation(2)
proc.time() - ptm
# stopCluster(cl = cl)

#    The diffrence between implementing 0d and 1d for any model:
#       1- For the Membrane voltage and Variables for 0d is nx and for 1d 0:nx
#       2- For updating the membrane voltage of nodes the loop for 0d is 1,nx and for 1d is 0,nx


