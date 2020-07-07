# SSDbain
SSDbain was built under R version 3.6.3.    
The function SSDttest is used to compute the sample size required per group for the Bayesian t-test and Bayesian Welch's test.
The function SSDANOVA and SSDANOVA_robust in the R package SSDbain computes the sample size for the Bayesian ANOVA, Welch's ANOVA, and robust ANOVA.

### install devtools
```
install.packages('devtools')
```

### Load devtools package for install_github()
```
library(devtools)

```
### get SSDbain from github
```
install_github("Qianrao-Fu/SSDbain",upgrade="never")

```
### Load SSDbain package
```
library(SSDbain)
```
