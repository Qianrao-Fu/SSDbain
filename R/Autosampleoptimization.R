#' @import conflicted
#' @import doParallel
#' @import stringr
#' @import pkgmaker
#' @import foreach
#' @import doRNG

SSDttest<-function(type='equal',Population_mean=c(0.5,0),var=NULL,BFthresh=3,eta=0.8,Hypothesis='two-sided',T=10000){
  options(warn=-1)
  # necessary <- c('bain', 'conflicted','doParallel','stringr','pkgmaker','foreach','doRNG')
  # installed <- necessary %in% installed.packages()[, 'Package']
  # if (length(necessary[!installed]) >=1)
  #   install.packages(necessary[!installed])
  library(stringr)
  library(conflicted)
  library(pkgmaker, warn.conflicts = FALSE)
  library(foreach)
  library(doParallel)
  library(doRNG)
  no_cores <- detectCores() - 1
  ## Loading required package: iterators
  cl <- makeCluster(no_cores)

  registerDoParallel(cl)
 # clusterEvalQ(cl, .libPaths("D:/Fu/library"))
  clusterEvalQ(cl, .libPaths())

  if(length(var)==0){
  if(type=='equal')
    var1<-var2<-1
  else{
    var1<-1.33
    var2<-0.67
  }

}
  else{
  var1<-var[1]
  var2<-var[2]
}


  N_SIM=T
  Nmin<-10

  Nmax<-1000
  m11<-m12<-0
  J<-1
  cat('sensitivity analysis with fraction b')
  cat('\n')
  cat('\n')

  if(length(Population_mean)==2)
  {m21<-Population_mean[1]
  m22<-Population_mean[2]}
  else
  {    set.seed(5)
    m21<-rnorm(T*10,mean=0,sd=sqrt(2/J))
    m22<-rnorm(T*10,mean=0,sd=sqrt(2/J))
  }
  # if(length(m21)==1){
  #   cat(paste('initial step: fast estimate the sample size ','\n'))
  #
  #   #general sample size impression when simulation data set equals to 99
  #   t1<-system.time(Res<-calSSD_fast(m21,m22,var1,var2,type,J,BFthresh,Nmin,Nmax,N_SIM,Hypothesis))
  #   cat('\n')
  #   #cat(paste('Program running time of step 1 is ',as.character(t1[[3]]),'\n',sep=''))
  #   #cat(paste('second step:  estimate the sample size based on the true mean and variance','\n'))
  #   N<-Res[[1]][[3]]
  #
  #   Nmin<-max(N-100,10)
  #   Nmax<-N+100
  #   #final sample size when simulation data set equals 1000
  #   if(N-100<10)
  #     cat(paste('final step: simulation',as.character(T), 'data sets with Nmin=10,','Nmax=',as.character(N),'+100','\n'))
  #   else
  #     cat(paste('final step: simulation',as.character(T), 'data sets with Nmin=',as.character(N),'-100,','Nmax=',as.character(N),'+100','\n'))
  # }


  t2=system.time(Res1<-calSSD(m21,m22,var1,var2,type,J=1,BFthresh,eta,Nmin,Nmax,N_SIM,Hypothesis))
  #cat(paste('Program running time of step 2 (J=1) is ',as.character(t2[[3]]),'\n',sep=''))
  #cat('\n')
  N<-Res1[[5]]
  Nmin<-10#max(N-200,10)
  Nmax<-1000#N+200
  J<-2
  cat('sensitivity analysis with fraction 2b')
  cat('\n')
  if(length(Population_mean)==2)
  {m21<-Population_mean[1]
  m22<-Population_mean[2]}
  else
  {    set.seed(5)
    m21<-rnorm(T*10,mean=0,sd=sqrt(2/J))
    m22<-rnorm(T*10,mean=0,sd=sqrt(2/J))
  }
  # if(length(m21)==1){
  #   cat(paste('initial step: fast estimate the sample size ','\n'))
  #
  #   #general sample size impression when simulation data set equals to 99
  #   t1<-system.time(Res<-calSSD_fast(m21,m22,var1,var2,type,J,BFthresh,eta,Nmin,Nmax,N_SIM,Hypothesis))
  #   cat('\n')
  #   #cat(paste('Program running time of step 1 is ',as.character(t1[[3]]),'\n',sep=''))
  #   #cat(paste('second step: detailed estimate the sample size based on the true mean and variance','\n'))
  #
  #
  #   N<-Res[[1]][[3]]
  #
  #   Nmin<-max(N-100,10)
  #   Nmax<-N+100
  #   #final sample size when simulation data set equals 1000
  #   if(N-100<10)
  #     cat(paste('final step: simulation',as.character(T), 'data sets with Nmin=10,','Nmax=',as.character(N),'+100','\n'))
  #   else
  #     cat(paste('final step: simulation',as.character(T), 'data sets with Nmin=',as.character(N),'-100,','Nmax=',as.character(N),'+100','\n'))
  # }

  if(length(Population_mean)==2)
  {m21<-Population_mean[1]
  m22<-Population_mean[2]}
  else
  {    set.seed(5)
    m21<-rnorm(T*10,mean=0,sd=sqrt(2/J))
    m22<-rnorm(T*10,mean=0,sd=sqrt(2/J))
  }
  t3=system.time(Res2<-calSSD(m21,m22,var1,var2,type,J=2,BFthresh,eta,Nmin,Nmax,N_SIM,Hypothesis))
  cat('\n')

  J<-3
  cat('sensitivity analysis with fraction 3b')
  cat('\n')
  N<-Res2[[5]]
  Nmin<-10#max(N-200,10)
  Nmax<-1000#N+200
  if(length(Population_mean)==2)
  {m21<-Population_mean[1]
  m22<-Population_mean[2]}
  else
  {    set.seed(5)
    m21<-rnorm(T*10,mean=0,sd=sqrt(2/J))
    m22<-rnorm(T*10,mean=0,sd=sqrt(2/J))
  }
  # if(length(m21)==1){
  #   cat(paste('initial step: fast estimate the sample size ','\n'))
  #
  #     #general sample size impression when simulation data set equals to 99
  #     t1<-system.time(Res<-calSSD_fast(m21,m22,var1,var2,type,J,BFthresh,eta,Nmin,Nmax,N_SIM,Hypothesis))
  #   cat('\n')
  #   #cat(paste('Program running time of step 1 is ',as.character(t1[[3]]),'\n',sep=''))
  #   #cat(paste('second step: detailed estimate the sample size based on the true mean and variance','\n'))
  #
  #
  #   N<-Res[[1]][[3]]
  #
  #   Nmin<-max(N-100,10)
  #   Nmax<-N+100
  #   #final sample size when simulation data set equals 1000
  #   if(N-100<10)
  #     cat(paste('final step: simulation',as.character(T), 'data sets with Nmin=10,','Nmax=',as.character(N),'+100','\n'))
  #   else
  #     cat(paste('final step: simulation',as.character(T), 'data sets with Nmin=',as.character(N),'-100,','Nmax=',as.character(N),'+100','\n'))
  # }
  # # cat(paste('Program running time of step 2 (J=2) is ',as.character(t3[[3]]),'\n',sep=''))
  # # cat('\n')

  if(length(Population_mean)==2)
  {m21<-Population_mean[1]
  m22<-Population_mean[2]}
  else
  {    set.seed(5)
    m21<-rnorm(T*10,mean=0,sd=sqrt(2/J))
    m22<-rnorm(T*10,mean=0,sd=sqrt(2/J))
  }
  t4=system.time(Res3<-calSSD(m21,m22,var1,var2,type,J=3,BFthresh,eta,Nmin,Nmax,N_SIM,Hypothesis))
  cat('\n')
  #  cat(paste('Program running time of step 2 (J=3) is ',as.character(t4[[3]]),'\n',sep=''))
  cat('\n')
  Res_output(m21,m22,N,J=1,BFthresh,eta,type,N_SIM,Res1,Hypothesis)
  Res_output(m21,m22,N,J=2,BFthresh,eta,type,N_SIM,Res2,Hypothesis)
  Res_output(m21,m22,N,J=3,BFthresh,eta,type,N_SIM,Res3,Hypothesis)
  #   cat(paste('Total program running time is ',as.character(t1[[3]]+t2[[3]]+t3[[3]]+t4[[3]]),'\n',sep=''))
  stopCluster(cl)

  }
