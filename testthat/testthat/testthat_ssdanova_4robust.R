rm(list=ls())
library(SSDbain)
library(testthat)
library(stringr)
library(conflicted)
library(pkgmaker, warn.conflicts = FALSE)
library(foreach)
library(doParallel)
library(doRNG)
library(pkgmaker, warn.conflicts = FALSE)
no_cores <- detectCores() - 1
## Loading required package: iterators
cl <- makeCluster(no_cores)

registerDoParallel(cl)

clusterEvalQ(cl, .libPaths())
#========================================================================================")
#  using SSDanova to test the BF01 for the two-sided ANOVA if H0 is true
#========================================================================================")
N<-100
varnames<-c('mu1','mu2','mu3','mu4')
hyp1<-"mu1=mu2=mu3=mu4"
hyp2<-"Ha"
BFthresh<-3
type<-'unequal'
T=10000
f1<-0
f2<-0.25
skews=c(1.75,1.75,1.75,1.75)
kurts=c(5.89,5.89,5.89,5.89)
var<-c(4/3,2/3,4/3,2/3)
N_min<-10
N_max<-1000
sym<-unlist(strsplit(hyp1,split='[>=]'))
k<-length(sym)
if(length(var)==0){
  if(type=="equal"){
    var<-c()
    var[1:k]<-1
  }else{
    var<-numeric(k)
    if((k %% 2) == 0) {
      for (i in 1:k){
        if((i %% 2) == 0){
          var[i]<-1
        }else{
          var[i]<-2
        }
      }
    } else {
      for (i in 1:k){
        if((i %% 2) == 0){
          var[i]<-1
        }else{
          var[i]<-2
        }
      }
      var[k]<-1
    }
    var<-var/sum(var)*k
  }
}
varnames<-c()
b<-character()
for (i in 1:k){
  varnames[i]=paste('mu',as.character(i),sep='')
}
sd<-sqrt(var)
EI_matrix1<-matrix_trans_anova(varnames,hyp1)
ERr1<-EI_matrix1[[1]]
IRr1<-EI_matrix1[[2]]
EI_matrix2<-matrix_trans_anova(varnames,hyp2)
ERr2<-EI_matrix2[[1]]
IRr2<-EI_matrix2[[2]]
if(length(f1)==k){
  mu1<-f1
  mu2<-f2
  logic1<-Istrue_anova(varnames,hyp1,f1)
  if(!logic1[[1]]){
    stop("Please check the values of mean for Hypothesis 1", call. = FALSE)
  }
  logic2<-Istrue_anova(varnames,hyp2,f2)
  if(!logic2[[1]]){
    stop("Please check the values of mean for Hypothesis 2", call. = FALSE)
  }

}else if(length(f1)==1){
  mu<-cal_mu_anova(k,f1,f2,hyp1,hyp2,ERr1,ERr2)
  mu1<-mu[[1]]
  mu2<-mu[[2]]
}else{
  stop('Please input a correct group of effect size or mean value!')
}


m<-mu1
J<-1
samph<-rep(N,length=length(m))
sampm<-samph/J
bf<-c()
var1<-c()
Sigma<-list()
datalist_temp<-matrix(0,length(m),N)
L_temp<-matrix(0,length(m),N)

#library(PearsonDS)
#library(gk)
#library(SimMultiCorrData)
#simulate data
set.seed(seed=10,kind="Mersenne-Twister",normal.kind="Inversion")
library(WRS2)
g<-c()
h<-c()
library(EnvStats)
if(var(skews)==0&&var(kurts==0))
{    skews0<<-skews[1]
kurts0<<-kurts[1]
#calculate g and h
result<-optim(c(0.1, 0.1), f_cost)
x<-result$par

g[1:k]<-x[1]
h[1:k]<-x[2]

}else{
  for (i in 1:k){
    skews0<<-skews[i]
    kurts0<<-kurts[i]
    # GA main function
    result<-optim(c(0.1, 0.1), f_cost)
    x<-result$par
    # y<-gh2sk(x[1],x[2])
    g[i]<-x[1]
    h[i]<-x[2]
  }}

datalist<-list()
estimate<-c()
BF<- foreach(i = 1:T) %dorng% {

  library(bain)
  var1<-c()
  Sigma<-list()
  data<-list()
  for (i in 1:length(m)) {
    X<-rnorm(N,0,1)
    A<-m[i]
    B<-sd[i]
    if(g[i]==0){
      data_temp<-A+B*exp(h[i]/2*X^2)*X
    }else{
      data_temp<-A+B*exp(h[i]/2*X^2)*(exp(g[i]*X)-1)/g[i]
    }
    data[[i]]<-data_temp
  }
  estimate<-lapply(1:length(m),function (x) mean(data[[x]],tr=0.2))
  Sigma<-lapply(1:length(m),function (x) matrix(WRS2::trimse(data[[x]],tr=0.2)^2,nrow=1,ncol=1))
  estimate<-unlist(estimate)
  names(estimate)<-varnames


  library(bain)
  capture.output(res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
  if(hyp2=='Hc'){
    bf1u<-res$fit$BF[[1]]
  } else{
    bf1u<-res$fit$Fit[[1]]/res$fit$Com[[1]]
  }
  if(hyp2=='Ha'){
    bf2u<-1
  }else{
    capture.output(res<-bain::bain(estimate,hyp2,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
    if(hyp2=='Hc'){
      bf2u<-res$fit$BF[[1]]
    } else{
      bf2u<-res$fit$Fit[[1]]/res$fit$Com[[1]]
    }}
  bf12<-bf1u/bf2u

  return(bf12)
}
bf12<-unlist(BF)
T0<-length(bf12)
kk<-lapply(1:T0, function (x) bf12[x]>BFthresh )
eta1_test<-sum(unlist(kk)*1)/T0

m<-mu2
set.seed(seed=10,kind="Mersenne-Twister",normal.kind="Inversion")
BF<- foreach(i = 1:T) %dorng% {

  library(bain)
  var1<-c()
  Sigma<-list()
  data<-list()
  for (i in 1:length(m)) {
    X<-rnorm(N,0,1)
    A<-m[i]
    B<-sd[i]
    if(g[i]==0){
      data_temp<-A+B*exp(h[i]/2*X^2)*X
    }else{
      data_temp<-A+B*exp(h[i]/2*X^2)*(exp(g[i]*X)-1)/g[i]
    }
    data[[i]]<-data_temp
  }
  estimate<-lapply(1:length(m),function (x) mean(data[[x]],tr=0.2))
  Sigma<-lapply(1:length(m),function (x) matrix(WRS2::trimse(data[[x]],tr=0.2)^2,nrow=1,ncol=1))
  estimate<-unlist(estimate)
  names(estimate)<-varnames


  library(bain)
  capture.output(res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
  if(hyp2=='Hc'){
    bf1u<-res$fit$BF[[1]]
  } else{
    bf1u<-res$fit$Fit[[1]]/res$fit$Com[[1]]
  }
  if(hyp2=='Ha'){
    bf2u<-1
  }else{
    capture.output(res<-bain::bain(estimate,hyp2,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
    if(hyp2=='Hc'){
      bf2u<-res$fit$BF[[1]]
    } else{
      bf2u<-res$fit$Fit[[1]]/res$fit$Com[[1]]
    }}
  bf21<-bf2u/bf1u

  return(bf21)
}
bf21<-unlist(BF)
T0<-length(bf21)
kk<-lapply(1:T0, function (x) bf21[x]>BFthresh )
eta2_test<-sum(unlist(kk)*1)/T0




eta<-cal_medbf_robust_anova(N,mu1,mu2,sd,g,h,T,J=1,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta=0.8,flag_fast=0,type,seed=10)
eta1<-eta[[1]]
eta2<-eta[[2]]

eta_seed<-cal_medbf_robust_anova(N,mu1,mu2,sd,g,h,T,J=1,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta=0.8,flag_fast=0,type,seed=100)
eta_seed1<-eta_seed[[1]]
eta_seed2<-eta_seed[[2]]

test_that("SSDANOVA", {expect_equal(eta1_test,eta1,tolerance = .01)})
test_that("SSDANOVA", {expect_equal(eta2_test, eta2,tolerance = .01)})
test_that("SSDANOVA", {expect_equal(min(eta1_test,eta2_test),0.8,tolerance = .03)})

test_that("SSDANOVA", {expect_equal(eta1,eta_seed1 ,tolerance = .01)})
test_that("SSDANOVA", {expect_equal(eta2,eta_seed2,tolerance = .01)})

