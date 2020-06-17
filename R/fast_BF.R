#BF01 when H0 is true
calBF01medium1_fast<-function(m1,m2,sd1,sd2,N,J,N_SIM,Variances){
  N_SIM<-1
  N1<-N2<-N
  samph<-c(N1,N2)
  ngroup<-samph/J
  #simulate data

  ##caculate estimate and variance for group mean
  estimate<-list(c(m1,m2))
  variance<-list(c(sd1^2/N1,sd2^2/N2))

  ##caculate bf
  library(bain)
  bf<-c()
  if(Variances=='equal'){
    for(i in 1:N_SIM){
      names(estimate[[i]]) <- c("mu1","mu2")
      variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
      covlist<-list(matrix(variance_sp),matrix(variance_sp))
      res<-bain(estimate[[i]],n=ngroup, "mu1 = mu2",Sigma=covlist,group_parameters=1, joint_parameters=0)
      bf[i]<-res$fit$BF[1]
    }
  }
  else {
    for(i in 1:N_SIM){
      names(estimate[[i]]) <- c("mu1","mu2")
      covlist<-list(matrix(variance[[i]][1]),matrix(variance[[i]][2]))
      res<-bain(estimate[[i]],n=ngroup, "mu1 = mu2",Sigma=covlist,group_parameters=1, joint_parameters=0)
      bf[i]<-res$fit$BF[1]
    }
  }

  if(m1==m2)
    medianbf<- median(bf)
  else
    medianbf<- 1/median(bf)
  return(medianbf)
}

##BF10 H1 is ture
calBF01medium2_fast<-function(m1,m2,sd1,sd2,N,J,N_SIM,Variances){
  N1<-N2<-N
  samph<-c(N1,N2)
  ngroup<-samph/J
  set.seed(10)
  ##simulate data
  estimate<-list()
  variance<-list()
  for(i in 1:N_SIM){
    estimate[[i]]<-c(m1[i],m2[i])
    variance[[i]]<-c(sd1^2/N1,sd2^2/N2)

  }

  ##caculates estimate and variance for group mean

  ##caculate bf
  library(bain)
  bf<-c()
  if(Variances=='equal'){
    for(i in 1:N_SIM){
      names(estimate[[i]]) <- c("mu1","mu2")
      variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
      covlist<-list(matrix(variance_sp),matrix(variance_sp))
      res<-bain(estimate[[i]],n=ngroup,"mu1 = mu2", Sigma=covlist,group_parameters = 1, joint_parameters=0)
      bf[i]<-res$fit$BF[1]
    }
  }
  else {
    for(i in 1:N_SIM){
      names(estimate[[i]]) <- c("mu1","mu2")
      covlist<-list(matrix(variance[[i]][1]),matrix(variance[[i]][2]))
      res<-bain(estimate[[i]],n=ngroup,"mu1 = mu2", Sigma=covlist,group_parameters=1, joint_parameters=0)
      bf[i]<-res$fit$BF[1]
    }
  }
  medianbf<-1/median(bf)

  return(medianbf)
}

##define function BF02 H0 is ture H0 VS H2
calBF02medium1_fast<-function(m1,m2,sd1,sd2,N,J,N_SIM,Variances){
  N_SIM<-1
  N1<-N2<-N
  samph<-c(N1,N2)
  ngroup<-samph/J
  ##simulate data
  ##Repeat the same simulation results

  estimate<-list(c(m1,m2))
  variance<-list(c(sd1^2/N1,sd2^2/N2))
  ##caculate estimate and variance for group mean

  ##caculate bf
  library(bain)
  bf<-c()
  bf1<-c()

  library(bain)
  bf<-c()
  if(Variances=='equal'){
    for(i in 1:N_SIM){
      names(estimate[[i]]) <- c("mu1","mu2")
      variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
      covlist<-list(matrix(variance_sp),matrix(variance_sp))
      res<-bain(estimate[[i]],n=ngroup,"mu1=mu2;mu1>mu2",Sigma=covlist,group_parameters=1, joint_parameters=0)
      bf_fit<-res$fit$Fit
      bf_com<-res$fit$Com
      bf1<-bf_fit[[1]]/bf_com[[1]]
      bf2<-bf_fit[[2]]/bf_com[[2]]
      bf[i]<-bf1/bf2
    }
  }
  else {
    for(i in 1:N_SIM){
      names(estimate[[i]]) <- c("mu1","mu2")
      covlist<-list(matrix(variance[[i]][1]),matrix(variance[[i]][2]))
      res<-bain(estimate[[i]],n=ngroup,"mu1=mu2;mu1>mu2",Sigma=covlist,group_parameters=1, joint_parameters=0)
      bf_fit<-res$fit$Fit
      bf_com<-res$fit$Com
      bf1<-bf_fit[[1]]/bf_com[[1]]
      bf2<-bf_fit[[2]]/bf_com[[2]]
      bf[i]<-bf1/bf2
    }
  }

  if(m1==m2){
    medianbf<- median(bf)
  }
  else{
    medianbf<- 1/median(bf)
  }
  return(medianbf)
}

calBF02medium2_fast<-function(m1,m2,sd1,sd2,N,J,N_SIM,Variances){
  N1<-N2<-N
  samph<-c(N1,N2)
  n<-samph/J
  ##simulate data
  set.seed(10)
  # m1<-rnorm(1000,mean=0,sd=sqrt(2))
  # m2<-rnorm(1000,mean=0,sd=sqrt(2))
  ##give an initial value
  miu1<-miu2<-rnorm(N1,0,1)*0
  datalist<-list()
  mu1<-c()
  mu2<-c()
  mu10<-m1[m1>m2]
  mu20<-m2[m1>m2]
  mu1<-mu10[1:N_SIM]
  mu2<-mu20[1:N_SIM]
  estimate<-list()
  variance<-list()
  for(i in 1:N_SIM){
    estimate[[i]]<-c(mu1[i],mu2[i])
    variance[[i]]<-c(sd1^2/N1,sd2^2/N2)

  }
##caculate bf
library(bain)
bf<-c()

if(Variances=='equal'){
  for(i in 1:N_SIM){

    names(estimate[[i]]) <- c("mu1","mu2")
    variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
    covlist<-list(matrix(variance_sp),matrix(variance_sp))
    #variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
    res<-bain(estimate[[i]],n=ngroup, "mu1=mu2;mu1>mu2",Sigma=covlist,group_parameters=1, joint_parameters=0)
    bf_fit<-res$fit$Fit
    bf_com<-res$fit$Com
    bf1<-bf_fit[[1]]/bf_com[[1]]
    bf2<-bf_fit[[2]]/bf_com[[2]]
    bf[i]<-bf1/bf2

  }
}
else {
  for(i in 1:N_SIM){
    names(estimate[[i]]) <- c("mu1","mu2")
    covlist<-list(matrix(variance[[i]][1]),matrix(variance[[i]][2]))
    res<-bain(estimate[[i]],n=ngroup, "mu1=mu2;mu1>mu2",Sigma=covlist,group_parameters=1, joint_parameters=0)
    bf_fit<-res$fit$Fit
    bf_com<-res$fit$Com
    bf1<-bf_fit[[1]]/bf_com[[1]]
    bf2<-bf_fit[[2]]/bf_com[[2]]
    bf[i]<-bf1/bf2
  }
}
medianbf<-1/median(bf)
return(medianbf)
}
