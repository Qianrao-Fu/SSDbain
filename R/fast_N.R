# judge whether med1 and med2 reach the threshold
calflagH1_fast<-function(m11,m12,m21,m22,sd1,sd2,N,J,MedBF,N_SIM,Variances){

  med1<-calBF01medium1_fast(m11,m12,sd1,sd2,N,J,N_SIM,Variances)
  #calculate median bf01 when H1 is ture
  if(length(m21)>1){
    med2<-calBF01medium2_fast(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
  }
  else{
    med2<-calBF01medium1_fast(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
  }
  if(med1>MedBF&&med2>MedBF)
    return(1)
  else
    return(0)
}

calH0vsH1SSD_fast<-function(m11,m12,m21,m22,J,sd1,sd2,MedBF=5,Nmin=10,Nmax=1000,N_SIM,Variances){
  #set the initial value
  med1<-0
  med2<-0
  N_a<-Nmin
  N_b<-Nmax
  cat(paste('N=',as.character(Nmin),'\n'))
  flag_a<-calflagH1_fast(m11,m12,m21,m22,sd1,sd2,N_a,J,MedBF,N_SIM,Variances)
  cat(paste('N=',as.character(Nmax),'\n'))
  flag_b<-calflagH1_fast(m11,m12,m21,m22,sd1,sd2,N_b,J,MedBF,N_SIM,Variances)
  if(flag_a+flag_b==0){
    cat('Nmax is too small, please increase the value of Nmax, and assign the previous Nmax value to Nmin')
    N<-Nmax
  }
  else{
    if(flag_a==1)
      N<-N_a
    else{
      N_mid<-floor((N_a+N_b)/2)
      while(N_b>N_a) {
        N_mid<-floor((N_a+N_b)/2)
        if(N_mid==N_a)
          N_mid=N_a+1
        cat(paste('N=',as.character(N_mid),'\n'))
        flag_mid<-calflagH1_fast(m11,m12,m21,m22,sd1,sd2,N_mid,J,MedBF,N_SIM,Variances)
        if(flag_mid==1){
          if(N_mid==N_a||N_mid==N_b||N_mid==N_a+1||N_mid==N_b-1)
          {N<-N_mid
          break}
          N_b<-N_mid
        }
        else{
          N_a<-N_mid
        }
      }
    }
  }
  cat(' \n')
    med1<-calBF01medium1_fast(m11,m12,sd1,sd2,N,J,N_SIM, Variances)

    #calculate median bf01 when H1 is ture
    if(length(m21)>1){
      med2<-calBF01medium2_fast(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
    }
    else{
      med2<-calBF01medium1_fast(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
    }
    BF1<-med1
    BF2<-med2
  Results<-list(BF1,BF2,N)
  return(Results)
}



calflag_fast<-function(m11,m12,m21,m22,sd1,sd2,N,J,MedBF,N_SIM, Variances){
  # if(N_SIM==99)

  med1<-calBF02medium1_fast(m11,m12,sd1,sd2,N,J,N_SIM, Variances)
  #calculate median bf02 when H2 is ture
  if(length(m21)>1){
    med2<-calBF02medium2_fast(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
  }
  else{
    med2<-calBF02medium1_fast(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
  }
  if(med1>MedBF&&med2>MedBF)
    return(1)
  else
    return(0)
}
calH0vsH2SSD_fast<-function(m11,m12,m21,m22,J,sd1,sd2,MedBF=5,Nmin=10,Nmax=1000,N_SIM,Variances){
  #set the initial value
  med1<-numeric(N_SIM)
  med2<-numeric(N_SIM)

  N_a<-Nmin
  N_b<-Nmax
  cat(paste('N=',as.character(Nmin),'\n'))
  flag_a<-calflag_fast(m11,m12,m21,m22,sd1,sd2,N_a,J,MedBF,N_SIM, Variances)
  cat(paste('N=',as.character(Nmax),'\n'))

  flag_b<-calflag_fast(m11,m12,m21,m22,sd1,sd2,N_b,J,MedBF,N_SIM, Variances)
  if(flag_a+flag_b==0){
    cat('Nmax is too small, please increase the value of Nmax, and assign the previous Nmax value to Nmin')
    N<-Nmax
  }
  else{
    if(flag_a==1)
      N<-N_a
    else{
      N_mid<-floor((N_a+N_b)/2)
      while(N_b>N_a) {
        N_mid<-floor((N_a+N_b)/2)
        if(N_mid==N_a)
          N_mid=N_a+1
        cat(paste('N=',as.character(N_mid),'\n'))
        flag_mid<-calflag_fast(m11,m12,m21,m22,sd1,sd2,N_mid,J,MedBF,N_SIM, Variances)
        if(flag_mid==1){
          if(N_mid==N_a||N_mid==N_b||N_mid==N_a+1||N_mid==N_b-1)
          {N<-N_mid
          break}
          N_b<-N_mid
        }
        else{
          N_a<-N_mid
        }
      }
    }
  }
  # cat(paste('The second verification','\n'))
  #
  # for(i in N:Nmax)
  # {
  #   cat(paste('N=',as.character(i),'\n'))

  #   if(med1>MedBF&&med2>MedBF){
  #     N<-i
  #     break
  #   }
  #
  # }
  cat(' \n')
    med1<-calBF02medium1_fast(m11,m12,sd1,sd2,N,J,N_SIM, Variances)

    #calculate median bf01 when H1 is ture
    if(length(m21)>1){
      med2<-calBF02medium2_fast(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
    }
    else{
      med2<-calBF02medium1_fast(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
    }

    BF1<-med1
    BF2<-med2

  Results<-list(BF1,BF2,N)
  return(Results)
}
calSSD_fast<-function(m21,m22,var1,var2,Variances,J,MedBF,Nmin,Nmax,N_SIM,Hypothesis){
  #set the mean of H0:mu1=mu2
  m11<-m12<-0

    sd1<-sqrt(var1)
    sd2<-sqrt(var2)

  # when Hypothesis is H0 vs H1
  if(Hypothesis=='two-sided'){
    #calculate the results (including median p-value(H0 is ture), median p-value(H1 is ture),median BF(H0 is ture), median BF(H1 is ture) and sample size) when the Hypothesis is "H0 vs H1"
    Results=calH0vsH1SSD_fast(m11,m12,m21,m22,J,sd1,sd2,MedBF,Nmin,Nmax,N_SIM,Variances)
    #extract the sample size from results. The variable Results is in the form of list.

    Results<-list(Results)
  }
  if(Hypothesis=='one-sided'){
    Results=calH0vsH2SSD_fast(m11,m12,m21,m22,J,sd1,sd2,MedBF,Nmin,Nmax,N_SIM,Variances)

    Results<-list(Results)

  }
  return(Results)
}
