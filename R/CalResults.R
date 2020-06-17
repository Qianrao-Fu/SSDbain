# judge whether med1 and med2 reach the threshold
calflagH1<-function(m11,m12,m21,m22,sd1,sd2,N,J,threshold,N_SIM,Variances,power){

    threshold1<-calBF01medium1(m11,m12,sd1,sd2,N,J,threshold,N_SIM,Variances,power)
    #calculate median bf01 when H1 is ture
    if(length(m21)>1){
      threshold2<-calBF01medium2(m21,m22,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
    }
    else{
      threshold2<-calBF01medium1(m21,m22,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
    }
    if(threshold1[[2]]>power&&threshold2[[2]]>power)
      return(1)
    else
      return(0)
}

calH0vsH1SSD<-function(m11,m12,m21,m22,J,sd1,sd2,threshold,power,Nmin,Nmax,N_SIM,Variances){
#set the initial value
med1<-0
med2<-0
#dichotomy
N_a<-Nmin
N_b<-Nmax
cat(paste('N=',as.character(Nmin),'\n'))
flag_a<-calflagH1(m11,m12,m21,m22,sd1,sd2,N_a,J,threshold,N_SIM,Variances,power)
cat(paste('N=',as.character(Nmax),'\n'))
flag_b<-calflagH1(m11,m12,m21,m22,sd1,sd2,N_b,J,threshold,N_SIM,Variances,power)

while(flag_b==0){
  cat('Nmax is too small, please increase the value of Nmax, and assign the previous Nmax value to Nmin','\n')
  N_a<-N_b
  N_b<-N_b+100
  cat(paste('N=',as.character(N_b),'\n'))
  flag_b<-calflagH1(m11,m12,m21,m22,sd1,sd2,N_b,J,threshold,N_SIM,Variances,power)
}
N_min<-N_a
N_max<-N_b
  if(flag_a==1)
    N<-N_a
  else{
    N_mid<-floor((N_a+N_b)/2)
    while(N_b>N_a) {
      N_mid<-floor((N_a+N_b)/2)
      if(N_mid==N_a)
        N_mid=N_a+1
      cat(paste('N=',as.character(N_mid),'\n'))
      flag_mid<-calflagH1(m11,m12,m21,m22,sd1,sd2,N_mid,J,threshold,N_SIM,Variances,power)
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


cat(' \n')
#prevent the genaral judgement computing and printing, and save time
  threshold1<-calBF01medium1(m11,m12,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
  med1<-threshold1[[1]]
  Error_H0_BF<-threshold1[[2]]
 # CI_H0_BF<-threshold1[[3]]
  #calculate median bf01 when H1 is ture
  if(length(m21)>1){
    threshold2<-calBF01medium2(m21,m22,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
  }
  else{
    threshold2<-calBF01medium1(m21,m22,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
  }
  med2<-threshold2[[1]]
  Error_H1_BF<-threshold2[[2]]
 # CI_H1_BF<-threshold2[[3]]
#Error_H0_BF<-calBF01error1(m11,m12,sd1,sd2,N,J,N_SIM, Variances,threshold)
#CI_H0_BF<-calBF01CI1(m11,m12,sd1,sd2,N,J,N_SIM, Variances)
# if(length(m21)==1)
# {
#   Error_H1_BF<-calBF01error1(m21,m22,sd1,sd2,N,J,N_SIM, Variances,threshold)
#   CI_H1_BF<-calBF01CI1(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
# }
# #when hypothesis H1 is a normal distribution
# if(length(m21)>1)
# {
#   Error_H1_BF<-calBF01error2(m21,m22,sd1,sd2,N,J,N_SIM, Variances,threshold)
#   CI_H1_BF<-calBF01CI2(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
# }
# calculate median bf
BF1<-med1
BF2<-med2

#Results<-list(BF1,BF2,Error_H0_BF,CI_H0_BF,Error_H1_BF,CI_H1_BF,N)


Results<-list(BF1,BF2,Error_H0_BF,Error_H1_BF,N)
return(Results)
}



calflag<-function(m11,m12,m21,m22,sd1,sd2,N,J,threshold,N_SIM, Variances,power){
  # if(N_SIM==99)

    threshold1<-calBF02medium1(m11,m12,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
    med1<-threshold1[[1]]
  #calculate median bf02 when H2 is ture
  if(length(m21)>1){
    threshold2<-calBF02medium2(m21,m22,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
  }
  else{
    threshold2<-calBF02medium1(m21,m22,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
  }
    med2<-threshold2[[1]]
  if(threshold1[[2]]>power&&threshold2[[2]]>power)
    return(1)
  else
    return(0)

 #  else{
 #  if(flag_H0){
 #   med1<-calBF02medium1(m11,m12,sd1,sd2,N,J,N_SIM, Variances)
 #
 #   if(med1>MedBF)
 #     return(1)
 #   else
 #     return(0)
 # }
 # else
 # {
 #
 #   #calculate median bf02 when H2 is ture
 #   if(length(m21)>1){
 #     med2<-calBF02medium2(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
 #   }
 #   else{
 #     med2<-calBF02medium1(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
 #   }
 #   if(med2>MedBF)
 #     return(1)
 #   else
 #     return(0)
 # }
 #  }
}
calH0vsH2SSD<-function(m11,m12,m21,m22,J,sd1,sd2,threshold,power,Nmin=10,Nmax=1000,N_SIM,Variances){
#set the initial value
 med1<-numeric(N_SIM)
 med2<-numeric(N_SIM)
# med1<-calBF02medium1(m11,m12,sd1,sd2,100,J,N_SIM, Variances)
#calculate median bf02 when H2 is ture
# if(length(m21)>1){
#   med2<-calBF02medium2(m21,m22,sd1,sd2,100,J,N_SIM, Variances)
# }
# else{
#   med2<-calBF02medium1(m21,m22,sd1,sd2,100,J,N_SIM, Variances)
# }

# if(med1<med2)
#  flag_H0=1
# else
#   flag_H0=0

N_a<-Nmin
N_b<-Nmax
cat(paste('N=',as.character(Nmin),'\n'))
flag_a<-calflag(m11,m12,m21,m22,sd1,sd2,N_a,J,threshold,N_SIM, Variances,power)
cat(paste('N=',as.character(Nmax),'\n'))

flag_b<-calflag(m11,m12,m21,m22,sd1,sd2,N_b,J,threshold,N_SIM, Variances,power)
while(flag_b==0){
  N_a=N_b
  N_b=2*N_a
  flag_b<-calflag(m11,m12,m21,m22,sd1,sd2,N_b,J,threshold,N_SIM, Variances,power)
  cat(paste('N=',as.character(N_b),'\n'))
}
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
    flag_mid<-calflag(m11,m12,m21,m22,sd1,sd2,N_mid,J,threshold,N_SIM, Variances,power)
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

  threshold1<-calBF02medium1(m11,m12,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
  med1<-threshold1[[1]]
 Error_H0_BF<-threshold1[[2]]
  #CI_H0_BF<-threshold1[[3]]
  #calculate median bf01 when H1 is ture
  if(length(m21)>1){
    threshold2<-calBF02medium2(m21,m22,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
  }
  else{
    threshold2<-calBF02medium1(m21,m22,sd1,sd2,N,J,threshold,N_SIM, Variances,power)
  }
  med2<-threshold2[[1]]
  Error_H2_BF<-threshold2[[2]]
  #CI_H2_BF<-medBF2[[3]]
  # Error_H0_BF<-calBF02error1(m11,m12,sd1,sd2,N,J,N_SIM, Variances,MedBF)
  # CI_H0_BF<-calBF02CI1(m11,m12,sd1,sd2,N,J,N_SIM, Variances)
  # if(length(m21)==1)
  # {
  #   Error_H2_BF<-calBF02error1(m21,m22,sd1,sd2,N,J,N_SIM, Variances,MedBF)
  #   CI_H2_BF<-calBF02CI1(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
  # }
  # if(length(m21)>1)
  # {
  #   Error_H2_BF<-calBF02error2(m21,m22,sd1,sd2,N,J,N_SIM, Variances,MedBF)
  #   CI_H2_BF<-calBF02CI2(m21,m22,sd1,sd2,N,J,N_SIM, Variances)
  # }
  BF1<-med1
  BF2<-med2

Results<-list(BF1,BF2,Error_H0_BF,Error_H2_BF,N)

#Results<-list(BF1,BF2,Error_H0_BF,CI_H0_BF,Error_H2_BF,CI_H2_BF,N)
  return(Results)
}
calSSD<-function(m21,m22,var1,var2,Variances,J,threshold,power,Nmin,Nmax,N_SIM,Hypothesis){
  #set the mean of H0:mu1=mu2
   m11<-m12<-0

    sd1<-sqrt(var1)
    sd2<-sqrt(var2)

   # when Hypothesis is H0 vs H1
  if(Hypothesis=='two-sided'){
#calculate the results (including median p-value(H0 is ture), median p-value(H1 is ture),median BF(H0 is ture), median BF(H1 is ture) and sample size) when the Hypothesis is "H0 vs H1"
    Results=calH0vsH1SSD(m11,m12,m21,m22,J,sd1,sd2,threshold,power,Nmin,Nmax,N_SIM,Variances)
#extract the sample size from results. The variable Results is in the form of list.

    #Results<-list(Results)
  }
  if(Hypothesis=='one-sided'){
    Results=calH0vsH2SSD(m11,m12,m21,m22,J,sd1,sd2,threshold,power,Nmin,Nmax,N_SIM,Variances)

    #Results<-list(Results)

  }
  return(Results)
}
