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

Res<-SSDANOVA(hyp1<-"mu1=mu2=mu3=mu4",hyp2="Ha",f1=0,f2=0.25,var=NULL,T=10000,BFthresh=3,eta=0.8,type='unequal')
N<-Res[[1]][1]
eta1<-Res[[2]][[1]][[1]]
eta2<-Res[[2]][[1]][[2]]



#The results are
#  N=71 when b is used

# the probability that the Bayes factor is larger than 3 when hyp1 or hyp2 is true
# P(BF01>3|H0) P(BF10>3|H1)
# 0.971        0.805
varnames<-c('mu1','mu2','mu3','mu4')
hyp1<-"mu1=mu2=mu3=mu4"
hyp2<-"Ha"
BFthresh<-3
type<-'unequal'
T=10000
f1<-0
f2<-0.25
EI_matrix1<-list()
EI_matrix1<-matrix_trans_anova(varnames,hyp1)
ERr1<-EI_matrix1[[1]]
IRr1<-EI_matrix1[[2]]
EI_matrix2<-list()
EI_matrix2<-matrix_trans_anova(varnames,hyp2)
ERr2<-EI_matrix2[[1]]
IRr2<-EI_matrix2[[2]]
k<-4
#calculate mu
mu<-cal_mu_anova(k,f1,f2,hyp1,hyp2,ERr1,ERr2)
var<-c(4/3,2/3,4/3,2/3)
m1<-mu[[1]]*sqrt(mean(var))
m2<-mu[[2]]*sqrt(mean(var))
sd<-sqrt(var)
J<-1

##Calculate eta1 and eta2 in SSDANOVA with the final sample size N=71
m<-m1
samph<-rep(N,length=length(m))
sampm<-samph/J
bf<-c()
var1<-c()
Sigma<-list()
datalist_temp<-matrix(0,length(m),N)
L_temp<-matrix(0,length(m),N)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

registerDoParallel(cl)

clusterEvalQ(cl, .libPaths())
#simulate data

set.seed(seed=10,kind="Mersenne-Twister",normal.kind="Inversion")
datalist<-list()
estimate<-c()
BF<- foreach(i = 1:T) %dorng% {
  # out <- matrix()

  library(bain)
  data<-lapply(1:length(m),function(x){
    rnorm(N,m[x],sd[x])})
  estimate<-lapply(1:length(m),function(x) mean(data[[x]]))
  Sigma<-lapply(1:length(m),function(x) matrix(var(data[[x]])/samph[x],nrow=1,ncol=1))
  estimate<-unlist(estimate)
  names(estimate)<-varnames
  #var1[l]<-var(a);
  # Sigma[[l]]<-matrix(var1[l]/n[l],nrow=1,ncol=1)

  capture.output(res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
  if(hyp2=='Hc'){
    bf1u<-res$fit$BF[[1]]
  }
  else{
    bf1u<-res$fit$Fit[[1]]/res$fit$Com[[1]]
  }
  if(hyp2=='Ha'){
    bf2u<-1
  }else{
    capture.output(res<-bain::bain(estimate,hyp2,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
    if(hyp2=='Hc'){
      bf2u<-res$fit$BF[[1]]
    }else{
      bf2u<-res$fit$Fit[[1]]/res$fit$Com[[1]]
    }}
  bf12<-bf1u/bf2u
  return(list(bf12))
}




bf12<-unlist(BF)
T0<-length(bf12)
kk<-lapply(1:T0, function (x) bf12[x]>BFthresh )
eta1_test<-sum(unlist(kk)*1)/T0

#simulate data
m<-m2
set.seed(seed=10,kind="Mersenne-Twister",normal.kind="Inversion")
datalist<-list()
estimate<-c()
BF<- foreach(i = 1:T) %dorng% {
  # out <- matrix()

  library(bain)

  data<-lapply(1:length(m),function(x){
    rnorm(N,m[x],sd[x])})
  estimate<-lapply(1:length(m),function(x) mean(data[[x]]))
  Sigma<-lapply(1:length(m),function(x) matrix(var(data[[x]])/samph[x],nrow=1,ncol=1))
  estimate<-unlist(estimate)
  names(estimate)<-varnames
  #var1[l]<-var(a);
  # Sigma[[l]]<-matrix(var1[l]/n[l],nrow=1,ncol=1)
  if(hyp2=='Ha'){
    bf2u<-1
  }else{
    capture.output(res<-bain::bain(estimate,hyp2,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
    if(hyp2=='Hc'){
      bf2u<-res$fit$BF[[1]]
    }
    else{
      bf2u<-res$fit$Fit[[1]]/res$fit$Com[[1]]
    }}

  capture.output(res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
  if(hyp2=='Hc'){
    bf1u<-res$fit$BF[[1]]
  }
  else{
    bf1u<-res$fit$Fit[[1]]/res$fit$Com[[1]]
  }
  bf21<-bf2u/bf1u
  return(bf21)
}
bf21<-unlist(BF)
T0<-length(bf21)
kk<-lapply(1:T0, function (x) bf21[x]>BFthresh )
eta2_test<-sum(unlist(kk)*1)/T0



eta_seed<- cal_medbf_anova(N,m1,m2,sd,T=10000,J=1,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast=0,type,seed=100)
eta_seed1<-eta_seed[[1]]
eta_seed2<-eta_seed[[2]]


test_that("SSDANOVA", {expect_equal(eta1_test,eta1 ,tolerance = .01)})
test_that("SSDANOVA", {expect_equal(eta2_test, eta2,tolerance = .01)})
test_that("SSDANOVA", {expect_equal(min(eta1_test,eta2_test),0.8,tolerance = .02)})


test_that("SSDANOVA", {expect_equal(eta_seed1,eta1 ,tolerance = .02)})
test_that("SSDANOVA", {expect_equal(eta_seed2,eta2 ,tolerance = .02)})



