Res_output<-function(m21,m22,N,J,threshold,power,Variances,N_SIM,Res,Hypothesis){
  m11<-m12<-0
  if(Variances=='equal'){
    sd1<-sd2<-1
  }
  else
  {
    sd1<-sqrt(1.33)
    sd2<-sqrt(0.67)
  }
  if(Hypothesis=='two-sided'){
  BF1<-Res[[1]]
  BF1<-round(Res[[1]],2)
  BF2<-round(Res[[2]],2)
  BF<-c(BF1,BF2)
  names(BF)<-c('BF01','BF10')
   Error_H0_BF<-Res[[3]]
  # CI_H0_BF<-round(Res[[4]],2)
   Error_H1_BF<-Res[[4]]
   error<-c(Error_H0_BF,Error_H1_BF)
   names(error)<-c(paste('P(BF01>',as.character(threshold),')',sep=''),paste('P(BF10>',as.character(threshold),')',sep=''))

  # CI_H1_BF<-round(Res[[6]],2)
   N<-Res[[5]]
   if(J==1){

   cat(paste('populations used for simulating data',' (T= ',as.character(N_SIM), '  variances= ',as.character(Variances),')',sep=''),'\n')
   cat(' \n')
   cat(paste(' H0: mu1 is equal to mu2   H1: mu1 is not equal to mu2'),'\n')
   cat(' \n')

if(length(m21)>1)
  switch(J,

         cat(paste('      mu1=mu2=',as.character(m11),'                mu1~N(',as.character(m11),',','2)',' mu2~N(',as.character(m12),',','2)','\n',sep='')),


         cat(paste('      mu1=mu2=',as.character(m11),'           mu1~N(',as.character(m11),',','1)',' mu2~N(',as.character(m12),',','1)','\n',sep='')),

         cat(paste('      mu1=mu2=',as.character(m11),'           mu1~N(',as.character(m11),',','2/3)',' mu2~N(',as.character(m12),',','2/3)','\n',sep=''))
  )
else
  cat(paste('      mu1=mu2=',as.character(m11), '              mu1=',as.character(m21),'  mu2=', as.character(m22),'\n',sep=''))

cat(paste(' var1=',as.character(sd1^2),    '  var2=', as.character(sd2^2), '          var1=', as.character(sd1^2), '  var2=', as.character(sd2^2), '\n',sep=''     ) )

}
cat(' \n')
if(J==1){
  cat(paste('Using N=',as.character(N),' and b',sep=''),'\n')
}
else
{
  cat(paste('Using N=',as.character(N),' and ',as.character(J),'b',sep=''),'\n')
}
cat(' \n')

#cat(paste('P(BF> ',as.character(threshold),')=',as.character(power),' is obtained using N = ',as.character(N),sep=''),'\n')

#cat(paste('A median Bayes factor of ',as.character(MedBF),' is obtained using N = ',as.character(N),sep=''),'\n')
#cat(paste(as.character(100-power*100), '% quantile'   ,sep=''))

#print(BF)
cat(' \n')
print(error)

# cat('Quantile of BF01: \n')
# print(CI_H0_BF)
# cat(' \n')
#
# cat('Quantile of BF10: \n')
# print(CI_H1_BF)
# cat(' \n')

# cat('  ',paste(as.character(signif(BF1,4)),' (',
#                as.character(signif(CI_H0_BF[1],4)),',',as.character(signif(CI_H0_BF[2],4)),')      ',
#                as.character(signif(BF2,4)),' (',
#                as.character(signif(CI_H1_BF[1],4)),',',as.character(signif(CI_H1_BF[2],4)),')',
#                '\n',sep=''))
# cat(' \n')
# cat(paste('   P(BF01<1|H0) = ',as.character(signif(Error_H0_BF[[1]],4)),'    P(BF10<1|H1) = ',as.character(signif(Error_H1_BF[[1]],4)),'\n'))
# cat(' \n')
# cat(paste('   P(BF01<1/3|H0) = ',as.character(signif(Error_H0_BF[[2]],4)),'\n',sep=''))
# cat(paste('   P(BF10<1/3|H1) = ',as.character(signif(Error_H1_BF[[3]],4)),'\n',sep=''))
# cat(paste('   (P(1/3<BF01<3|H0)+P(1/3<BF10<3|H1))/2  = ' ,as.character(signif((Error_H0_BF[[4]]+Error_H1_BF[[4]])/2),4),'\n',sep=''))
cat(' \n')
}
else if (Hypothesis=='one-sided'){
  BF1<-round(Res[[1]],2)
  BF2<-round(Res[[2]],2)
  BF<-c(BF1,BF2)
  names(BF)<-c('BF02','BF20')
   Error_H0_BF<-Res[[3]]
  # CI_H0_BF<-round(Res[[4]],2)
   Error_H2_BF<-Res[[4]]
  # CI_H2_BF<-round(Res[[6]],2)
   error<-c(Error_H0_BF,Error_H2_BF)
   names(error)<-c(paste('P(BF02>',as.character(threshold),')',sep=''),paste('P(BF20>',as.character(threshold),')',sep=''))

  N<-Res[[5]]
 if(J==1) {
  cat(paste('populations used for simulating data',' (T= ',as.character(N_SIM), '  variances= ',as.character(Variances),')',sep=''),'\n')
  cat(' \n')
cat(paste(' H0: mu1 is equal to mu2       H2: mu1 is greater than mu2'),'\n')
cat(' \n')
if(length(m21)>1)
  switch(J,

         cat(paste('      mu1=mu2=',as.character(m11),'        mu1~N(',as.character(m11),',','2)',' mu2~N(',as.character(m12),',','2)','\n',sep='')),


         cat(paste('      mu1=mu2=',as.character(m11),'        mu1~N(',as.character(m11),',','1)',' mu2~N(',as.character(m12),',','1)','\n',sep='')),

         cat(paste('      mu1=mu2=',as.character(m11),'        mu1~N(',as.character(m11),',','2/3)',' mu2~N(',as.character(m12),',','2/3)','\n',sep=''))
  )
else
  cat(paste('      mu1=mu2=',as.character(m11),'                   mu1=',as.character(m21),'  mu2=', as.character(m22),'\n',sep=''))

cat(paste(' var1=',as.character(sd1^2),    '  var2=', as.character(sd2^2), '          var1=', as.character(sd1^2), '  var2=', as.character(sd2^2), '\n',sep=''     ) )
}
  cat(' \n')
  if(J==1){
    cat(paste('Using N=',as.character(N),' and b',sep=''),'\n')
  }
  else
  {
    cat(paste('Using N=',as.character(N),' and ',as.character(J),'b',sep=''),'\n')
  }
  cat(' \n')

#cat(paste('P(BF> ',as.character(threshold),')=',as.character(power),' is obtained using N = ',as.character(N),sep=''),'\n')
#cat(paste('A median Bayes factor of ',as.character(MedBF),' is obtained using N = ',as.character(N),sep=''),'\n')
#cat(paste(as.character(100-power*100), '% quantile'   ,sep=''))
#print(BF)
cat(' \n')
print(error)

# cat('Quantile of BF02: \n')
# print(CI_H0_BF)
# cat(' \n')
#
# cat('Quantile of BF20: \n')
# print(CI_H2_BF)
# cat(' \n')
#
cat(' \n')
# cat(paste('   P(BF02<1|H0) = ',as.character(signif(Error_H0_BF[[1]],4)),'        P(BF20<1|H2) = ',as.character(signif(Error_H2_BF[[1]],4)),'\n'))
# cat(' \n')
# cat(paste('   P(BF02<1/3|H0) =  ',as.character(signif(Error_H0_BF[[2]],4)),'\n',sep=''))
# cat(paste('   P(BF20<1/3|H2) =  ',as.character(signif(Error_H2_BF[[3]],4)),'\n',sep=''))
# cat(paste('   (P(1/3<BF02<3|H0)+P(1/3<BF20<3|H2))/2  = ',as.character(signif((Error_H0_BF[[4]]+Error_H2_BF[[4]])/2),4),'\n',sep=''))
# cat(' \n')
}
else {
  cat('Please input a correct Hypothesis type!',' \n')
}
}
