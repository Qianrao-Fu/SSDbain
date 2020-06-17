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
#To save time, omit the running below, but use the sample size directly.
# Res<-SSDANOVA(hyp1<-"mu1=mu2=mu3",hyp2="mu1>mu2>mu3",f1=0,f2=0.25,var=NULL,T=10000,BFthresh=3,eta=0.8,type='equal')
# N<-Res[[1]][1]
# eta1<-Res[[2]][[1]][[1]]
# eta2<-Res[[2]][[1]][[2]]
#The hypothesis can be H0 vs Ha, H0 vs H1, H1 vs Hc, H1 vs H2
varnames<-c('mu1','mu2','mu3')
hyp1<-"mu1=mu2>mu3"
hyp2<-"mu2=mu3>mu1"
BFthresh<-3
type<-'equal'
T=10000
f1<-0.1
f2<-0.1
EI_matrix1<-list()
EI_matrix1<-matrix_trans_anova(varnames,hyp1)
ERr1<-EI_matrix1[[1]]
IRr1<-EI_matrix1[[2]]
EI_matrix2<-list()
EI_matrix2<-matrix_trans_anova(varnames,hyp2)
ERr2<-EI_matrix2[[1]]
IRr2<-EI_matrix2[[2]]
k<-3
#calculate mu
mu<-cal_mu_anova(k,f1,f2,hyp1,hyp2,ERr1,ERr2)
m1<-mu[[1]]
m2<-mu[[2]]
sd<-1
J<-1

N<-103
##Calculate eta1 and eta2 in SSDANOVA software with the final sample size N=71
eta<- cal_medbf_anova(N,m1,m2,sd,T,J=1,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast=0,type,seed=10)
m<-m1
samph<-rep(N,length=length(m))
sampm<-samph/J
bf<-c()
var1<-c()
Sigma<-list()
datalist_temp<-matrix(0,length(m),N)
L_temp<-matrix(0,length(m),N)

#Calculate eta1,eta2 manually.

set.seed(seed=10,kind="Mersenne-Twister",normal.kind="Inversion")
datalist<-list()
estimate<-c()
BF<- foreach(i = 1:T) %dorng% {
  # out <- matrix()

  library(bain)


  datalist_temp<-unlist(lapply(1:length(m), function(x) rnorm(N,m[x],mean(sd))))
  L_temp<-unlist(lapply(1:length(m), function(x) rep(x,length=N)))

  datalist<-matrix(datalist_temp,nrow=length(m)*N, ncol=1)
  L<-matrix(L_temp,nrow=length(m)*N, ncol=1)
  dd<- data.frame(datalist,L)

  mu <- factor(dd$L)
  #sampm<-table(mu)
  y <- lm(dd$datalist~mu-1)
  estimate <- coef(y)
  n <- table(dd$L)

  var1 <- summary(y)$sigma**2
  Sigma<-lapply(1:length(m), function(x) matrix(var1/n[x],nrow=1,ncol=1))



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
    }
    else{
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

  datalist_temp<-unlist(lapply(1:length(m), function(x) rnorm(N,m[x],mean(sd))))
  L_temp<-unlist(lapply(1:length(m), function(x) rep(x,length=N)))

  datalist<-matrix(datalist_temp,nrow=length(m)*N, ncol=1)
  L<-matrix(L_temp,nrow=length(m)*N, ncol=1)
  dd<- data.frame(datalist,L)

  mu <- factor(dd$L)
  #sampm<-table(mu)
  y <- lm(dd$datalist~mu-1)
  estimate <- coef(y)
  n <- table(dd$L)

  var1 <- summary(y)$sigma**2
  Sigma<-lapply(1:length(m), function(x) matrix(var1/n[x],nrow=1,ncol=1))
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

eta1<-eta[[1]]
eta2<-eta[[2]]

eta_seed<- cal_medbf_anova(N,m1,m2,sd,T,J=1,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast=0,type,seed=100)
eta_seed1<-eta_seed[[1]]
eta_seed2<-eta_seed[[2]]

#test Tables in paper
test_that("SSDANOVA", {expect_equal(eta1_test,eta1 ,tolerance = .01)})
test_that("SSDANOVA", {expect_equal(eta2_test,eta2,tolerance = .01)})
#test min(eta1,eta2) vs 0.8
test_that("SSDANOVA", {expect_equal(min(eta1_test,eta2_test),0.8,tolerance = .02)})
#test set.seed=10 vs seet.seed=100
test_that("SSDANOVA", {expect_equal(eta1,eta_seed1 ,tolerance = .01)})
test_that("SSDANOVA", {expect_equal(eta2,eta_seed2,tolerance = .01)})




