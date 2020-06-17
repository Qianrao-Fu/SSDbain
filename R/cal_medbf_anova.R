cal_medbf_anova<-function(N,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed){
  if(hyp2=='Hc'){
    bf0j_H0<-cal_bf_anova(N,mu1,sd,T,J,varnames,hyp1,flag_output=0,flag_fast,type,hyp2,seed)
  }
  else{
    #BF0u(BFiu)
    bf0u_H0<-cal_bf_anova(N,mu1,sd,T,J,varnames,hyp1,flag_output=0,flag_fast,type,hyp2,seed)
    #BFju(BFju)
    if(length(ERr2)==0&&length(IRr2)==0){
      bfju_H0<-1
    }
    else{
      bfju_H0<-cal_bf_anova(N,mu1,sd,T,J,varnames,hyp2,flag_output=0,flag_fast,type,hyp2,seed)
    }
    #BF0j(BFij)
    #medbf0j_H0<-median(bf0u_H0/bfju_H0)
    bf0j_H0<-bf0u_H0/bfju_H0
  }
  # medbf0j_H0<-median(bf0j_H0)
  bf0j_H0<-na.omit(bf0j_H0)
  T0<-length(bf0j_H0)
  k<-lapply(1:T0, function (x) bf0j_H0[x]>BFthresh )
  E_H0<-sum(unlist(k)*1)/T0
  # k<-lapply(1:T, function (x) bf0j_H0[x]<1/BFthresh )
  # M_H0<-sum(unlist(k))/T
  # k<-lapply(1:T, function (x) bf0j_H0[x]>1/BFthresh&bf0j_H0[x]<BFthresh )
  # W_H0<-sum(unlist(k))/T
  # CI_H0<-round(quantile(bf0j_H0,c(0.05,0.1,0.2,0.3,0.4,0.5)),2)

  #cat(paste('Median BF0j = ',as.character(medbf0j_H0),'\n',sep=''))


  #Hj is true
  if(hyp2=='Hc'){
    bfj0_Hj<-cal_bf_anova(N,mu2,sd,T,J,varnames,hyp1,flag_output=0,flag_fast,type,hyp2,seed)
    bfj0_Hj<-1/bfj0_Hj
    }
  else{
  #BF0u(BFiu)
  bf0u_Hj<-cal_bf_anova(N,mu2,sd,T,J,varnames,hyp1,flag_output=0,flag_fast,type,hyp2,seed)
  #BFju(BFju)
  if(length(ERr2)==0&&length(IRr2)==0){
    bfju_Hj<-1}
  else{
    bfju_Hj<-cal_bf_anova(N,mu2,sd,T,J,varnames,hyp2,flag_output=0,flag_fast,type,hyp2,seed)
  }
  #BFj0(BFji)
  bfj0_Hj<-bfju_Hj/bf0u_Hj
  }
  bfj0_Hj<-na.omit(bfj0_Hj)
  T0<-length(bfj0_Hj)
  k<-lapply(1:T0, function (x) bfj0_Hj[x]>BFthresh )
  E_Hj<-sum(unlist(k)*1)/T0
  # medianbfj0_Hj<-median(bfj0_Hj)
  # k<-lapply(1:T, function (x) bfj0_Hj[x]<1 )
  # E_Hj<-sum(unlist(k))/T0
  # k<-lapply(1:T, function (x) bfj0_Hj[x]<1/BFthresh )
  # M_Hj<-sum(unlist(k))/T0
  # k<-lapply(1:T, function (x) bf0j_H0[x]>1/BFthresh&bf0j_H0[x]<BFthresh )
  # W_Hj<-sum(unlist(k))/T0
  # CI_Hj<-round(quantile(bfj0_Hj,c(0.05,0.1,0.2,0.3,0.4,0.5)),2)
  #cat(paste('Median BFj0 = ',as.character(medianbfj0_Hj),'\n',sep=''))

  return(list(E_H0,E_Hj))
}



cal_medbf_robust_anova<-function(N,mu1,mu2,sd,g,h,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed){
  if(hyp2=='Hc'){
    bf0j_H0<-cal_bf_robust_anova(N,mu1,sd,g,h,T,J,varnames,hyp1,flag_output=0,flag_fast,type,hyp2,seed)
  }
  else{
    #BF0u(BFiu)
    bf0u_H0<-cal_bf_robust_anova(N,mu1,sd,g,h,T,J,varnames,hyp1,flag_output=0,flag_fast,type,hyp2,seed)
    #BFju(BFju)
    if(length(ERr2)==0&&length(IRr2)==0){
      bfju_H0<-1
    }
    else{
      bfju_H0<-cal_bf_robust_anova(N,mu1,sd,g,h,T,J,varnames,hyp2,flag_output=0,flag_fast,type,hyp2,seed)
    }
    #BF0j(BFij)
    #medbf0j_H0<-median(bf0u_H0/bfju_H0)
    bf0j_H0<-bf0u_H0/bfju_H0
  }
  # medbf0j_H0<-median(bf0j_H0)
  bf0j_H0<-na.omit(bf0j_H0)
  T0<-length(bf0j_H0)
  k<-lapply(1:T0, function (x) bf0j_H0[x]>BFthresh )
  E_H0<-sum(unlist(k)*1)/T0
  # k<-lapply(1:T, function (x) bf0j_H0[x]<1/BFthresh )
  # M_H0<-sum(unlist(k))/T
  # k<-lapply(1:T, function (x) bf0j_H0[x]>1/BFthresh&bf0j_H0[x]<BFthresh )
  # W_H0<-sum(unlist(k))/T
  # CI_H0<-round(quantile(bf0j_H0,c(0.05,0.1,0.2,0.3,0.4,0.5)),2)

  #cat(paste('Median BF0j = ',as.character(medbf0j_H0),'\n',sep=''))


  #Hj is true
  if(hyp2=='Hc'){
    bfj0_Hj<-cal_bf_robust_anova(N,mu2,sd,g,h,T,J,varnames,hyp1,flag_output=0,flag_fast,type,hyp2,seed)
    bfj0_Hj<-1/bfj0_Hj
  }
  else{
    #BF0u(BFiu)
    bf0u_Hj<-cal_bf_robust_anova(N,mu2,sd,g,h,T,J,varnames,hyp1,flag_output=0,flag_fast,type,hyp2,seed)
    #BFju(BFju)
    if(length(ERr2)==0&&length(IRr2)==0){
      bfju_Hj<-1}
    else{
      bfju_Hj<-cal_bf_robust_anova(N,mu2,sd,g,h,T,J,varnames,hyp2,flag_output=0,flag_fast,type,hyp2,seed)
    }
    #BFj0(BFji)
    bfj0_Hj<-bfju_Hj/bf0u_Hj
  }
  bfj0_Hj<-na.omit(bfj0_Hj)
  T0<-length(bfj0_Hj)
  k<-lapply(1:T0, function (x) bfj0_Hj[x]>BFthresh )
  E_Hj<-sum(unlist(k)*1)/T0
  # medianbfj0_Hj<-median(bfj0_Hj)
  # k<-lapply(1:T, function (x) bfj0_Hj[x]<1 )
  # E_Hj<-sum(unlist(k))/T0
  # k<-lapply(1:T, function (x) bfj0_Hj[x]<1/BFthresh )
  # M_Hj<-sum(unlist(k))/T0
  # k<-lapply(1:T, function (x) bf0j_H0[x]>1/BFthresh&bf0j_H0[x]<BFthresh )
  # W_Hj<-sum(unlist(k))/T0
  # CI_Hj<-round(quantile(bfj0_Hj,c(0.05,0.1,0.2,0.3,0.4,0.5)),2)
  #cat(paste('Median BFj0 = ',as.character(medianbfj0_Hj),'\n',sep=''))

  return(list(E_H0,E_Hj))
}
