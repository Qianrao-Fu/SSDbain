mu_Hc_anova<-function(k,mu){
  library(gtools)
  x <- c(1:k)
  #pick 2 balls from the urn with replacement
  #get all permutations
  database=permutations(n=k,r=k,v=x)
  num=nrow(database)

  database_temp<-list()
  database_temp[[1]]=database[1,]
  num_total<-sum(1:k-1)+1
  j<-1
  while(j<num_total){
    num_j<-length(database_temp[[j]])/k
    matrix_temp<-matrix(database_temp[[j]],num_j,k)
    vector_diff<-c()
    x<-1
    flag_row<-c()

    for (i in 1:num_j)
    {
      for (ii in 1:num)
      {
        vector_diff<-database[ii,]-matrix_temp[i,]
        index<-which(vector_diff!=0)
        if(length(index)==2)
        {
          if((index[2]-index[1])==1)
          {flag_row[x]<-ii
          x<-x+1}
        }
      }
    }
    database_temp[[j+1]]<-unique(database[flag_row,])
    j<-j+1
  }

  data_final<-rbind(database_temp[[1]],database_temp[[2]])
  for (i in 3:num_total)
  {
    data_final<-rbind(data_final,database_temp[[i]])
  }
  data_final<-unique(data_final)
  seq<-data_final[median(2:num),]
  return(mu[seq])}


Istrue_anova<-function(varnames,hyp,mu){
  for (i in 1:length(mu)){
    assign(varnames[i],mu[i])
  }
  temp <- unlist(gregexpr("[>=]", hyp))
  vari<-unlist(strsplit(hyp,split='[>=]'))
  logic<-numeric(length(vari)-1)
  for (i in 1:(length(vari)-1)){
    hyp_sub<- paste(vari[i],substr(hyp, temp[i],temp[i]),vari[i+1])
    if(substr(hyp, temp[i],temp[i])=='=')
    {
      hyp_sub<- paste(vari[i+1],'==',vari[i])
  #    eval(parse(text = hyp_sub))
    }
    ##hyp_sub<-gsub('=','==',hyp_sub)
    logic[i]<-(eval(parse(text = hyp_sub)))
  }
  for (i in 1:length(mu)){
    mu[i]<-get(varnames[i])
  }
  Logic<-sum(logic)==length(logic)
  return(list(Logic,mu))
}

cal_mu_anova<-function(k,f1,f2,hyp1,hyp2,ERr1,ERr2){

  library(stringr)

  if(f1==0)
  {
    mu1<-numeric(k)
  }
  else{
    hyp<-hyp1
    f<-f1
    mu<-numeric(k)
    mu[1]<-1;
    for (n in 1:(k-1)){
      eq_num<-0
      if(nrow(ERr1))
        for (l in 1:nrow(ERr1)){
          if(ERr1[l,n]==1){
            eq_num<-1
            break;
          }
        }
      else
        eq_num<-0;
      if(eq_num)
      {
        mu[n+1]<-mu[n];
      }
      else
      {
        mu[n+1]<-mu[n]-1;
      }
    }
    #mu[k]<-0;
    mu<-mu-mu[k]
    mu<-mu*f/(sqrt((k-1)/k)*sd(mu))
    mu1<-mu
    index<-unlist(str_extract_all(hyp1,"[0-9]"))
    index<-as.numeric(index)
    for (n in 1:k)
    {
      mu1[index[n]]<-mu[n]
    }

  }
  if(hyp2=='Hc')
  {
    #######mu2######
    hyp<-hyp2
    f<-f2
    #mu2<-mu_Hc_anova(k,mu1)
    #
    #order<-[3 1 2] [3 1 4 2] [5 1 2 3 4][ 1 5 6 4 3 2]
    if(k==3){
      order<-c(3, 1, 2)
    }else if(k==4){
      order<- c(3, 1, 4, 2)
    }else if(k==5){
     order<-c(5, 1, 2, 3, 4)
    }else if (k==6){
     order<-c( 1, 5, 6, 4, 3, 2)
    }else{
      stop('k should be between 3 and 6')
    }
    mu<-numeric(k)
    for (i in 1:k){
      mu[i]<-mu1[order[i]]
    }
    mu<-mu-min(mu)
    mu<-mu*f2/(sqrt((k-1)/k)*sd(mu))
    mu2<-mu
  }
  else if(hyp2=='Ha')
  {
    mu<-numeric(k)

    for (n in 1:(k-1)){

        mu[n+1]<-mu[n]-1;

    }
    #mu[k]<-0;
    mu<-mu-min(mu)
    mu<-mu*f2/(sqrt((k-1)/k)*sd(mu))
    mu2<-mu
  }
  else{
    hyp<-hyp2
    f<-f2
    mu<-numeric(k)
    mu[1]<-1;
    for (n in 1:(k-1)){
      eq_num<-0
      if(nrow(ERr2))
        for (l in 1:nrow(ERr2)){
          if(ERr1[l,n]==1){
            eq_num<-1
            break;
          }
        }
      else
        eq_num<-0;
      if(eq_num)
      {
        mu[n+1]<-mu[n];
      }
      else
      {
        mu[n+1]<-mu[n]-1;
      }
    }
    #mu[k]<-0;
    mu<-mu-min(mu)
    mu<-mu*f/(sqrt((k-1)/k)*sd(mu))
    mu2<-mu
    index<-unlist(str_extract_all(hyp2,"[0-9]"))
    index<-as.numeric(index)
    for (n in 1:k)
    {
      mu2[index[n]]<-mu[n]
    }
    }

  return(list(mu1,mu2))
}



matrix_trans_anova<-function(varnames,hyp1){
  Rr1<-create_matrices(varnames, hyp1)
  Rr1$ERr1[is.null(Rr1$ERr1)] <- 0
  nrow_ERr1<-nrow(Rr1$ERr1)
  ncol_ERr1<-ncol(Rr1$ERr1)
  nrow_ERr1[is.null(nrow_ERr1)] <- 0
  ncol_ERr1[is.null(ncol_ERr1)] <- 0

  Rr1$IRr1[is.null(Rr1$IRr1)] <- 0
  nrow_IRr1<-nrow(Rr1$IRr1)
  ncol_IRr1<-ncol(Rr1$IRr1)
  nrow_IRr1[is.null(nrow_IRr1)] <- 0
  ncol_IRr1[is.null(ncol_IRr1)] <- 0


  ERr1<-matrix(Rr1$ERr1,nrow_ERr1,ncol_ERr1)
  IRr1<-matrix(Rr1$IRr1,nrow_IRr1,ncol_IRr1)
  return(list(ERr1,IRr1))
}
