cal_bf_anova<-function(N,m,sd,T,J,varnames,hyp1,flag_output,flag_fast,type,hyp2,seed){
  #   N1<-N2<-N3<-N
  samph<-rep(N,length=length(m))
  sampm<-samph/J
  bf<-c()
  var1<-c()
  Sigma<-list()
  datalist_temp<-matrix(0,length(m),N)
  L_temp<-matrix(0,length(m),N)


  #simulate data
  set.seed(seed,kind="Mersenne-Twister",normal.kind="Inversion")

  datalist<-list()
  estimate<-c()
  #bf<-list()
  #flag_fast<-0
  if(!flag_fast)
  {
    # set.seed(10)
    #  seed <- doRNGseed()

    #    for (i in 1:T){
      BF<- foreach(i = 1:T) %dorng% {
       # out <- matrix()

        library(bain)
        if (type=='equal'){
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
        }
        else if(type=='unequal'){
          data<-lapply(1:length(m),function(x){
            rnorm(N,m[x],sd[x])})
          estimate<-lapply(1:length(m),function(x) mean(data[[x]]))
          Sigma<-lapply(1:length(m),function(x) matrix(var(data[[x]])/samph[x],nrow=1,ncol=1))
          estimate<-unlist(estimate)
          names(estimate)<-varnames
          #var1[l]<-var(a);
          # Sigma[[l]]<-matrix(var1[l]/n[l],nrow=1,ncol=1)

        }
        else{
          var1<-c()
          Sigma<-list()
          data<-lapply(1:length(m),function (x){
            rnorm(N,m[x],sd[x])})
          estimate<-lapply(1:length(m),function (x) mean(data[[x]],tr=0.2))
          Sigma<-lapply(1:length(m),function (x) matrix(WRS2::trimse(data[[x]],tr=0.2)^2,nrow=1,ncol=1))
          estimate<-unlist(estimate)
          names(estimate)<-varnames
        }

        library(bain)
        if(flag_output)
        {
          res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0)
        }
        else{
          capture.output(res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
        }
        if(hyp2=='Hc'){
          bf<-res$fit$BF[[1]]
        }
        else{
          bf<-res$fit$Fit[[1]]/res$fit$Com[[1]]
        }
        return(bf)
      }


    bf<-unlist(BF)
  }
  else
  {
    estimate<-c()
    Sigma<-list()
    for (x in 1:length(m)){
      estimate[x]<-m[x]
      Sigma[[x]]<-matrix(sd[x]^2/samph[x],nrow=1,ncol=1)
    }
    names(estimate)<-varnames#c("mu1","mu2","mu3","mu4")

    if(flag_output)
      res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0)
    else
      capture.output(res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
    if(hyp2=='Hc'){
      bf<-res$fit$BF[[1]]
    }
    else{
      bf<-res$fit$Fit[[1]]/res$fit$Com[[1]]
    }
    return(bf)

  }

  #medianbf<-median(bf)
  return(bf)
}


cal_bf_robust_anova<-function(N,m,sd,g,h,T,J,varnames,hyp1,flag_output,flag_fast,type,hyp2,seed){
  #   N1<-N2<-N3<-N
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
  set.seed(seed,kind="Mersenne-Twister",normal.kind="Inversion")
  library(WRS2)

    datalist<-list()
    estimate<-c()
    #bf<-list()
    #flag_fast<-0
    if(!flag_fast)
    {
      # set.seed(10)
      #  seed <- doRNGseed()
      #    for (i in 1:T){
      BF<- foreach(i = 1:T) %dorng% {
        # out <- matrix()

        library(bain)


          var1<-c()
          Sigma<-list()

          # data<-lapply(1:length(m),function(x){
          #   rgh(N,m[x],sd[x],g=skews[x],h=kurts[x])})


          # data<-list()
        # for (x in 1:length(m))
        # {
        #   seed_N<-round(runif(1,1,T*10))
        #
        #
        #   capture.output(data_res<-SimMultiCorrData::nonnormvar1(method = "Fleishman",
        #                                                          means = m[x], vars = sd[x]^2, skews = skews[x],
        #                                                          kurts = kurts[x], fifths = 0, sixths = 0,
        #                                                          Six = NULL, cstart = NULL, n = N,seed = seed_N), file='NUL')
        #   data[[x]]<-data_res$continuous_variable$V1
        # }

        # data=mnonr::unonr(N, mu=m, Sigma=diag(sd^2), skewness = skews, kurtosis = kurts, empirical = FALSE)
        #
        # data=split(data, rep(1:ncol(data), each = nrow(data)))
        data<-list()
        for (i in 1:length(m)) {
          #data_temp <-PearsonDS::rpearson(N, moments=c(mean=m[i],variance=sd[i]^2,skewness=skews[i],kurtosis=kurts[i]))
          #data_temp<-gk::rgh(N,m[i],sd[i],g[i],h[i])
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
      if(flag_output)
      {
        res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0)
      }
      else{
        capture.output(res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
      }
      if(hyp2=='Hc'){
        bf<-res$fit$BF[[1]]
      }
      else{
        bf<-res$fit$Fit[[1]]/res$fit$Com[[1]]
      }
      return(bf)
    }


    bf<-unlist(BF)
  }
  else
  {
    estimate<-c()
    Sigma<-list()
    for (x in 1:length(m)){
      estimate[x]<-m[x]
      Sigma[[x]]<-matrix(sd[x]^2/samph[x],nrow=1,ncol=1)
    }
    names(estimate)<-varnames#c("mu1","mu2","mu3","mu4")

    if(flag_output)
      res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0)
    else
      capture.output(res<-bain::bain(estimate,hyp1,n=sampm,Sigma=Sigma,group_parameters=1,joint_parameters = 0), file='NUL')
    if(hyp2=='Hc'){
      bf<-res$fit$BF[[1]]
    }
    else{
      bf<-res$fit$Fit[[1]]/res$fit$Com[[1]]
    }
    return(bf)

  }

  #medianbf<-median(bf)
  return(bf)
}




