
expect_error(SSDANOVA_robust(hyp1="mu1=mu2=mu3",hyp2="mu1>mu2>mu3",f1=0,f2=0.25,skews=c(1.75,1.75,1.75),kurts=c(5.89,5.89,5.89),var=c(1.5,0.75,0.75),BFthresh=0.5,eta=0.8,T=10000,seed=10))
expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="Hc",type="unequal",f1=0.25,f2=0.25,var=c(1.5,0.75,0.75),BFthresh=1,eta=0.8,T=10000,seed=10))
