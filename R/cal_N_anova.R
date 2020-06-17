gh2sk<-function(g,h){
  if(g==0){
    skews<-0
    kurts<-3*(1-2*h)^3*(1/(1-4*h)^(5/2)+1/(2*h-1)^3)
  }else{
    skews<-((3*exp(g^2/(2-6*h))+exp(9*g^2/(2-6*h))-3*exp(2*g^2/(1-3*h))-1)/
              (1-3*h)^(1/2)-3*(1-2*exp(g^2/(2-4*h))+exp(2*g^2/(1-2*h)))*(exp(g^2/(2-2*h))-1)/
              ((1-2*h)^(1/2))+2*(exp(g^2/(2-2*h))-1)^3/(1-h)^(3/2))/
      (g^3*(((1-2*exp(g^2/(2-4*h))+exp(2*g^2/(1-2*h)))/(1-2*h)^(1/2)+
               (exp(g^2/(2-2*h))-1)^2/(h-1))/g^2)^(3/2))

    kurts<-(exp(8*g^2/(1-4*h))*(1+6*exp(6*g^2/(4*h-1))+exp(8*g^2/(4*h-1))-4*exp(7*g^2/(8*h-2))-
                                  4*exp(15*g^2/(8*h-2)))/(1-4*h)^(1/2)-4*(3*exp(g^2/(2-6*h))+
                                                                            exp(9*g^2/(2-6*h))-3*exp(2*g^2/(1-3*h))-1)*(exp(g^2/(2-2*h))-1)/
              ((1-3*h)^(1/2)*(1-h)^(1/2))-6*(exp(g^2/(2-2*h))-1)^4/(h-1)^2-12*(1-2*exp(g^2/(2-4*h))+
                                                                                 exp(2*g^2/(1-2*h)))*
              (exp(g^2/(2-2*h))-1)^2/((1-2*h)^(1/2)*(h-1))+3*(1-2*exp(g^2/(2-4*h))+
                                                                exp(2*g^2/(1-2*h)))^2/(2*h-1))/
      ((1-2*exp(g^2/(2-4*h))+exp(2*g^2/(1-2*h)))/(1-2*h)^(1/2)+(exp(g^2/(2-2*h))-1)^2/(h-1))^2
  }
  return(c(skews,kurts))
}

f_cost<-function(x){
  y<-gh2sk(x[1],x[2])
  return(c((skews0-y[1])^2+(kurts0-y[2])^2))
}


SSDANOVA<-function(hyp1="mu1=mu2=mu3",hyp2="mu1>mu2>mu3",type="equal",f1,f2,var=NULL,BFthresh=3,eta=0.8,T=10000,seed=10){
  if(BFthresh<1){
    stop("BFthresh should be larger than 1!")
  }
  if(eta<0){
    stop("eta should be larger than 0")
  }
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
  }else{
    if(length(var)!=k){
      stop("Please check the variance!",call. = TRUE)
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
    if(hyp2!='Hc'&&hyp2!='Ha'){
    logic2<-Istrue_anova(varnames,hyp2,f2)
    if(!logic2[[1]]){
      stop("Please check the values of mean for Hypothesis 2", call. = FALSE)
    }}

  }
  else if(length(f1)==1){
    sym_eq<- unlist(strsplit(hyp1,split='[>]'))
    if(length(sym_eq)==1){
      if(f1!=0){
        stop("Please check the effect size!",call. = TRUE)
      }
    }
    mu<-cal_mu_anova(k,f1,f2,hyp1,hyp2,ERr1,ERr2)
    mu1<-mu[[1]]*sqrt(mean(var))
    mu2<-mu[[2]]*sqrt(mean(var))
  }
  else{
    stop('Please input a correct group of effect size or mean value!')
  }
  N_J<-c()
  medBF_J<-list()
  N_min_J<-c(N_min,N_min,N_min)
  N_max_J<-c(N_max,N_max,N_max)
  options(warn=-1)

  # library(parallel)
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
 # cl <- parallel::makeCluster(46)

  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths())
  if(length(ERr1)|length(ERr1))
    Jmax<-3
  else
    Jmax<-1
  for (J in 1:Jmax){
    if(J>1){
      cat(paste('###Sensitive Analysis  ', as.character(J), '*', 'b','###','\n',sep=''))
    }
    ##step1 sample size for the population that the sample exactly equal to  the population mean and variance
    #increase sample size step by step till the BF reach the threshold
    N_min<-N_min_J[J]
    N_max<-N_max_J[J]
    # cat(paste('Initial Step: fast determine sample size with the real mean and variance.','\n',sep=''))
    # capture.output(Results_Step1<-Fast_SSD(N_min,N_max,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BF_thresh), file='NUL')
    # N_Step1<-Results_Step1[[1]]
    # medBF0j_H0_Step1<-Results_Step1[[2]]
    # medBFj0_Hj_Step1<-Results_Step1[[3]]
    # cat(paste('The initial sample size is ',as.character(N_Step1),'\n',sep=''))
    # cat(paste('Median BF0j is ',as.character(signif(medBF0j_H0_Step1,4)),'\n',sep=''))
    # cat(paste('Median BFj0 is ',as.character(signif(medBFj0_Hj_Step1,4)),'\n','\n',sep=''))
    #
    # #step2 substitute the N_Step1 to the function
    # cat(paste('Second Step: exactly determine sample size.','\n',sep=''))
    N_Step1<-100
    flag_fast<-0
    medBF<- cal_medbf_anova(N_Step1,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
    medBF0j_H0_Step2<-medBF[[1]]
    medBFj0_Hj_Step2<-medBF[[2]]
    cat(paste('The sample size is N=',as.character(N_Step1),'\n',sep=''))
    cat(paste('P(BF12>',as.character(BFthresh),') is ',as.character(signif(medBF0j_H0_Step2,4)),'\n',sep=''))
    cat(paste('P(BF21>',as.character(BFthresh),') is ',as.character(signif(medBFj0_Hj_Step2,4)),'\n','\n',sep=''))

    if(medBF[[1]]>eta&&medBF[[2]]>eta)
    {
      N_max<-N_Step1
      #cat(paste('The median BF is larger than the threshold. Therefore,','\n',sep=''))
      cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n','\n',sep=''))
      N_min<-max(floor(N_max/2),10);
      medBF<- cal_medbf_anova(N_min,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
      cat(paste('The lower boundary of sample size is N_min=',as.character(N_min),'\n',sep=''))
      N_temp<-c()
      while((medBF[[1]]>eta&&medBF[[2]]>eta)&&N_min>10){
        N_temp<-N_min
        N_min<-max(floor(N_min/2),10);
        medBF<- cal_medbf_anova(N_min,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
        cat(paste('The lower boundary of sample size is N_min=',as.character(N_min),'\n',sep=''))
      }
      if(N_min*2<(N_Step1-1))
      {
        N_max<-N_temp
      }
      else
        N_max<-N_Step1#min(N_min*2,N_Step1)
      cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n','\n',sep=''))
    }
    else{
      N_max<-N_Step1
      cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n','\n',sep=''))

      while(!(medBF[[1]]>eta&&medBF[[2]]>eta)){
        N_max<-2*N_max;
        medBF<- cal_medbf_anova(N_max,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
        cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n',sep=''))
      }
      if(N_max>N_Step1*2+1)
        N_min<-N_max/2
      else
        N_min<-N_Step1#max(N_max/2,N_Step1)
      cat(paste('The lower boundary of sample size is N_min=',as.character(N_min),'\n','\n',sep=''))
    }
    flag_end<-0
    if(N_min==10)
    {
      medBF<- cal_medbf_anova(N_min,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
      if(medBF[[1]]>eta&&medBF[[2]]>eta)
      {
        flag_end<-1
        N<-N_min
        medBF_N10<-medBF
        #cat(paste('N=',as.character(N),'\n',sep=''))
      }
    }


    N_med<-floor((N_min+N_max)/2)
    while((N_med-N_min>=1&&N_max-N_med>=1)&&!flag_end){
      cat(paste('The sample size in current loop is N=',as.character(N_med),'\n',sep=''))

      medBF<- cal_medbf_anova(N_med,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
      cat(paste('P(BF12>',as.character(BFthresh),') is ',as.character(signif(medBF[[1]],4)),'\n',sep=''))
      cat(paste('P(BF21>',as.character(BFthresh),') is ',as.character(signif(medBF[[2]],4)),'\n','\n',sep=''))

      if(medBF[[1]]>eta&&medBF[[2]]>eta){
        N_max<-N_med
      }
      else{
        N_min<-N_med
      }
      N_med<-floor((N_min+N_max)/2)
    }
    if(flag_end){
      cat(paste('The final sample size is N=',as.character(N),'\n',sep=''))
      medBF<-medBF_N10

    }
    else{
      N<-N_max
      cat(paste('The final sample size is N=',as.character(N),'\n',sep=''))
      medBF<- cal_medbf_anova(N,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)

      # medBF_Vec<-c(medBF[[1]],medBF[[2]])
      # E_Vec<-c(medBF[[3]],medBF[[4]])
      # WE<-(medBF[[9]]+medBF[[10]])/2
      # M_Vec<-c(medBF[[5]],medBF[[6]],WE)
      # CI_Vec<-c(medBF[[7]],medBF[[8]])
      # names(medBF_Vec)<-c('medBF01','medBF10')
      # names(E_Vec)<-c('Error I','Error II')
      # names(M_Vec)<-c('ME I ','ME II','WE')
      # print(medBF_Vec)
      # cat(paste('\n',sep=''))
      # print(E_Vec)
      # cat(paste('\n',sep=''))
      # print(M_Vec)
      # cat(paste('\n',sep=''))
      # cat(paste('CI1',sep=''))
      # cat(paste('\n',sep=''))
      # print(medBF[[7]])
      # cat(paste('CI2',sep=''))
      # cat(paste('\n',sep=''))
      # print(medBF[[8]])

    }

    # medBF<- cal_medbf_anova(N_med,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BF_thresh,flag_fast,type,seed)
    # if(medBF[[1]]>eta&&medBF[[2]]>eta){
    #   N<-N_med
    # }
    # else{
    #   N<-N_max
    #   cat(paste('N=',as.character(N),'\n',sep=''))
    #   medBF<- cal_medbf_anova(N,mu1,mu2,sd,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BF_thresh,flag_fast,type,seed)
    #
    # }
    N_J[J]<-N
    medBF_J[[J]]<-medBF
  }
  for (J in 1:Jmax){
    N<-N_J[J]
    medBF<-medBF_J[[J]]
    medBF_Vec<-c(medBF[[1]],medBF[[2]])
    # E_Vec<-c(medBF[[3]],medBF[[4]])
    # WE<-(medBF[[9]]+medBF[[10]])/2
    # M_Vec<-c(medBF[[5]],medBF[[6]],WE)
    # CI_Vec<-c(medBF[[7]],medBF[[8]])
    if(Jmax>1 && J>1)
    {cat(paste('#####Sensitive analysis ',as.character(J), '*', 'b','#####','\n',sep=''))
    }else{
      cat(paste('#####Sensitive analysis ','b','#####','\n',sep=''))
    }
    cat(paste('Hypothesis: H1: ',hyp1,' vs H2: ',hyp2,'\n',sep=''))
   # cat(paste('Effect size: f1=',as.character(f1),' vs f2=',as.character(f2),'\n',sep=''))
    if(type=='equal'){
      cat(paste('Type: ANOVA','\n',sep=''))
    }else if(type=='unequal'){
      cat(paste('Type: Welch','\'s ANOVA','\n',sep=''))
    }else if(type=='robust'){
      cat(paste('Type: robust ANOVA','\n',sep=''))
    }else{
      cat(paste('Please check your input about type!','\n',sep=''))
    }
    cat(paste('\n',sep=''))
    cat(paste('The sample size is N=',as.character(N),'\n',sep=''))

    cat(paste('P(BF12>',as.character(BFthresh),'|H1)=',as.character(medBF[[1]]),'\n',sep=''))
    cat(paste('P(BF21>',as.character(BFthresh),'|H2)=',as.character(medBF[[2]]),'\n',sep=''))
    cat(paste('\n',sep=''))

    # names(medBF_Vec)<-c('eta01','eta10')
    # names(E_Vec)<-c('Error I','Error II')
    # names(M_Vec)<-c('ME I ','ME II','WE')
    # print(medBF_Vec)
    # cat(paste('\n',sep=''))
    # print(E_Vec)
    # cat(paste('\n',sep=''))
    # print(M_Vec)
    # cat(paste('\n',sep=''))
    # cat(paste('quantile1',sep=''))
    # cat(paste('\n',sep=''))
    # print(medBF[[7]])
    # cat(paste('quantile2',sep=''))
    # cat(paste('\n',sep=''))
    # print(medBF[[8]])

  }


  # cat(paste('Effect size:','\n',sep=''))
  # cat(paste(as.character(f1),' vs ',as.character(f2),'\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('BFthreth:','\n',sep=''))
  # cat(paste(as.character(BFthresh),'\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('Variance:','\n',sep=''))
  # cat(paste(as.character(var),'\n',sep=''))
  #
  # cat(paste('\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('Type:','\n',sep=''))
  # cat(paste(type,'\n',sep=''))

  stopCluster(cl)

  return(list(N_J,medBF_J))

}



SSDANOVA_robust<-function(hyp1,hyp2,f1,f2,var=NULL,skews,kurts,T,BFthresh,eta,seed=10){
  if(BFthresh<1){
    stop("BFthresh should be larger than 1!")
  }
  if(eta<0){
    stop("eta should be larger than 0")
  }
   sym<- unlist(strsplit(hyp1,split='[>=]'))
   k<-length(sym)

   if(length(var)==0){
     var<-numeric(k)
     if((k %% 2) == 0) {
       for (i in 1:k){
         if((i %% 2) == 0){
           var[i]<-1
         }else{
           var[i]<-2
         }
       }
     }  else {
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



  N_min<-10
  N_max<-1000
 # library(mcga)
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
  }

}


  varnames<-c()
  type<-'robust'
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
  if(length(var)<k){
  var[1:k]<-1
  }
  if(length(f1)==k){
    mu1<-f1
    mu2<-f2
    logic1<-Istrue_anova(varnames,hyp1,f1)
    if(!logic1[[1]]){
      stop("Please check the values of mean for Hypothesis 1")
    }
    if(hyp2!='Hc'&&hyp2!='Ha'){
      logic2<-Istrue_anova(varnames,hyp2,f2)
      if(!logic2[[1]]){
        stop("Please check the values of mean for Hypothesis 2", call. = FALSE)
      }}

  }
  else if(length(f1)==1){
    sym_eq<- unlist(strsplit(hyp1,split='[>]'))
    if(length(sym_eq)==1){
      if(f1!=0){
        stop("Please check the effect size!",call. = TRUE)
      }
    }
    mu<-cal_mu_anova(k,f1,f2,hyp1,hyp2,ERr1,ERr2)
    mu1<-mu[[1]]*sqrt(mean(var))
    mu2<-mu[[2]]*sqrt(mean(var))
  }
  else{
    stop('Please input a correct group of effect size or mean value!')
  }
  N_J<-c()
  medBF_J<-list()
  N_min_J<-c(N_min,N_min,N_min)
  N_max_J<-c(N_max,N_max,N_max)
  options(warn=-1)

  # library(parallel)
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
  #cl <- parallel::makeCluster(46)

  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths())
  if(length(ERr1)|length(ERr1))
    Jmax<-3
  else
    Jmax<-1
  for (J in 1:Jmax){
    if(J>1){
      cat(paste('###Sensitive Analysis  ', as.character(J), '*', 'b','###','\n',sep=''))
    }
    ##step1 sample size for the population that the sample exactly equal to  the population mean and variance
    #increase sample size step by step till the BF reach the threshold
    N_min<-N_min_J[J]
    N_max<-N_max_J[J]
    # cat(paste('Initial Step: fast determine sample size with the real mean and variance.','\n',sep=''))
    # capture.output(Results_Step1<-Fast_SSD(N_min,N_max,mu1,mu2,sd,skews,kurts,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BF_thresh), file='NUL')
    # N_Step1<-Results_Step1[[1]]
    # medBF0j_H0_Step1<-Results_Step1[[2]]
    # medBFj0_Hj_Step1<-Results_Step1[[3]]
    # cat(paste('The initial sample size is ',as.character(N_Step1),'\n',sep=''))
    # cat(paste('Median BF0j is ',as.character(signif(medBF0j_H0_Step1,4)),'\n',sep=''))
    # cat(paste('Median BFj0 is ',as.character(signif(medBFj0_Hj_Step1,4)),'\n','\n',sep=''))
    #
    # #step2 substitute the N_Step1 to the function
    # cat(paste('Second Step: exactly determine sample size.','\n',sep=''))
    N_Step1<-100
    flag_fast<-0
    medBF<- cal_medbf_robust_anova(N_Step1,mu1,mu2,sd,g,h,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
    medBF0j_H0_Step2<-medBF[[1]]
    medBFj0_Hj_Step2<-medBF[[2]]
    cat(paste('The sample size is ',as.character(N_Step1),'\n',sep=''))
    cat(paste('P(BF12>',as.character(BFthresh),') is ',as.character(signif(medBF0j_H0_Step2,4)),'\n',sep=''))
    cat(paste('P(BF21>',as.character(BFthresh),') is ',as.character(signif(medBFj0_Hj_Step2,4)),'\n','\n',sep=''))

    if(medBF[[1]]>eta&&medBF[[2]]>eta)
    {
      N_max<-N_Step1
      #cat(paste('The median BF is larger than the threshold. Therefore,','\n',sep=''))
      cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n','\n',sep=''))
      N_min<-max(floor(N_max/2),10);
      medBF<- cal_medbf_robust_anova(N_min,mu1,mu2,sd,g,h,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
      cat(paste('The lower boundary of sample size is N_min=',as.character(N_min),'\n',sep=''))
      N_temp<-c()
      while((medBF[[1]]>eta&&medBF[[2]]>eta)&&N_min>10){
        N_temp<-N_min
        N_min<-max(floor(N_min/2),10);
        medBF<- cal_medbf_robust_anova(N_min,mu1,mu2,sd,g,h,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
        cat(paste('The lower boundary of sample size is N_min=',as.character(N_min),'\n',sep=''))
      }
      if(N_min*2<(N_Step1-1))
      {
        N_max<-N_temp
      }
      else
        N_max<-N_Step1#min(N_min*2,N_Step1)
      cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n','\n',sep=''))
    }
    else{
      N_max<-N_Step1
      cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n','\n',sep=''))

      while(!(medBF[[1]]>eta&&medBF[[2]]>eta)){
        N_max<-2*N_max;
        medBF<- cal_medbf_robust_anova(N_max,mu1,mu2,sd,g,h,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
        cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n',sep=''))
      }
      if(N_max>N_Step1*2+1)
        N_min<-N_max/2
      else
        N_min<-N_Step1#max(N_max/2,N_Step1)
      cat(paste('The lower boundary of sample size is N_min=',as.character(N_min),'\n','\n',sep=''))
    }
    flag_end<-0
    if(N_min==10)
    {
      medBF<- cal_medbf_robust_anova(N_min,mu1,mu2,sd,g,h,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
      if(medBF[[1]]>eta&&medBF[[2]]>eta)
      {
        flag_end<-1
        N<-N_min
        medBF_N10<-medBF
        #cat(paste('N=',as.character(N),'\n',sep=''))
      }
    }


    N_med<-floor((N_min+N_max)/2)
    while((N_med-N_min>=1&&N_max-N_med>=1)&&!flag_end){
      cat(paste('The sample size in current loop is N=',as.character(N_med),'\n',sep=''))

      medBF<- cal_medbf_robust_anova(N_med,mu1,mu2,sd,g,h,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)
      cat(paste('P(BF12>',as.character(BFthresh),') is ',as.character(signif(medBF[[1]],4)),'\n',sep=''))
      cat(paste('P(BF21>',as.character(BFthresh),') is ',as.character(signif(medBF[[2]],4)),'\n','\n',sep=''))

      if(medBF[[1]]>eta&&medBF[[2]]>eta){
        N_max<-N_med
      }
      else{
        N_min<-N_med
      }
      N_med<-floor((N_min+N_max)/2)
    }
    if(flag_end){
      cat(paste('The final sample size is N=',as.character(N),'\n',sep=''))
      medBF<-medBF_N10

    }
    else{
      N<-N_max
      cat(paste('The final sample size is N=',as.character(N),'\n',sep=''))
      medBF<- cal_medbf_robust_anova(N,mu1,mu2,sd,g,h,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BFthresh,eta,flag_fast,type,seed)

      # medBF_Vec<-c(medBF[[1]],medBF[[2]])
      # E_Vec<-c(medBF[[3]],medBF[[4]])
      # WE<-(medBF[[9]]+medBF[[10]])/2
      # M_Vec<-c(medBF[[5]],medBF[[6]],WE)
      # CI_Vec<-c(medBF[[7]],medBF[[8]])
      # names(medBF_Vec)<-c('medBF01','medBF10')
      # names(E_Vec)<-c('Error I','Error II')
      # names(M_Vec)<-c('ME I ','ME II','WE')
      # print(medBF_Vec)
      # cat(paste('\n',sep=''))
      # print(E_Vec)
      # cat(paste('\n',sep=''))
      # print(M_Vec)
      # cat(paste('\n',sep=''))
      # cat(paste('CI1',sep=''))
      # cat(paste('\n',sep=''))
      # print(medBF[[7]])
      # cat(paste('CI2',sep=''))
      # cat(paste('\n',sep=''))
      # print(medBF[[8]])

    }

    # medBF<- cal_medbf_robust_robust(N_med,mu1,mu2,sd,g,h,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BF_thresh,flag_fast,type,seed)
    # if(medBF[[1]]>eta&&medBF[[2]]>eta){
    #   N<-N_med
    # }
    # else{
    #   N<-N_max
    #   cat(paste('N=',as.character(N),'\n',sep=''))
    #   medBF<- cal_medbf_robust_anova(N,mu1,mu2,sd,g,h,T,J,varnames,hyp1,hyp2,ERr2,IRr2,BF_thresh,flag_fast,type,seed)
    #
    # }
    N_J[J]<-N
    medBF_J[[J]]<-medBF
  }
  for (J in 1:Jmax){
    N<-N_J[J]
    medBF<-medBF_J[[J]]
    medBF_Vec<-c(medBF[[1]],medBF[[2]])
    # E_Vec<-c(medBF[[3]],medBF[[4]])
    # WE<-(medBF[[9]]+medBF[[10]])/2
    # M_Vec<-c(medBF[[5]],medBF[[6]],WE)
    # CI_Vec<-c(medBF[[7]],medBF[[8]])
    if(Jmax>1 && J>1)
    {cat(paste('#####Sensitive analysis ',as.character(J), '*', 'b','#####','\n',sep=''))
    }else{
      cat(paste('#####Sensitive analysis ','b','#####','\n',sep=''))
    }
    # cat(paste('The final sample size is N=',as.character(N),'\n',sep=''))
    #
    # names(medBF_Vec)<-c('eta01','eta10')
    # names(E_Vec)<-c('Error I','Error II')
    # names(M_Vec)<-c('ME I ','ME II','WE')
    # print(medBF_Vec)
    # cat(paste('\n',sep=''))
    # print(E_Vec)
    # cat(paste('\n',sep=''))
    # print(M_Vec)
    # cat(paste('\n',sep=''))
    # cat(paste('quantile1',sep=''))
    # cat(paste('\n',sep=''))
    # print(medBF[[7]])
    # cat(paste('quantile2',sep=''))
    # cat(paste('\n',sep=''))
    # print(medBF[[8]])
    cat(paste('Hypothesis: H1: ',hyp1,' vs H2: ',hyp2,'\n',sep=''))
   # cat(paste('Effect size: f1=',as.character(f1),' vs f2=',as.character(f2),'\n',sep=''))
    if(type=='equal'){
      cat(paste('Type: ANOVA','\n',sep=''))
    }else if(type=='unequal'){
      cat(paste('Type: Welch','\'s ANOVA','\n',sep=''))
    }else if(type=='robust'){
      cat(paste('Type: robust ANOVA','\n',sep=''))
    }else{
      cat(paste('Please check your input about type!','\n',sep=''))
    }
    cat(paste('\n',sep=''))
    cat(paste('The sample size is N=',as.character(N),'\n',sep=''))

    cat(paste('P(BF12>',as.character(BFthresh),'|H1)=',as.character(medBF[[1]]),'\n',sep=''))
    cat(paste('P(BF21>',as.character(BFthresh),'|H2)=',as.character(medBF[[2]]),'\n',sep=''))
    cat(paste('\n',sep=''))

  }

  # cat(paste('\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('Hypothesis:','\n',sep=''))
  # cat(paste(hyp1,' vs ',hyp2,'\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('Effect size:','\n',sep=''))
  # cat(paste(as.character(f1),' vs ',as.character(f2),'\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('BFthreth:','\n',sep=''))
  # cat(paste(as.character(BFthresh),'\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('Variance:','\n',sep=''))
  # cat(paste(as.character(var),'\n',sep=''))
  #
  # cat(paste('\n',sep=''))
  # cat(paste('\n',sep=''))
  # cat(paste('Type:','\n',sep=''))
  # cat(paste(type,'\n',sep=''))

  stopCluster(cl)

  return(list(N_J,medBF_J))

}





