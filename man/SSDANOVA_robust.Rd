\name{SSDANOVA_robust}
\alias{SSDANOVA_robust}
\title{Sample size determination for the Bayesian robust ANOVA with bain}
\description{
The function \code{SSDANOVA_robust} in the \code{R} package \code{SSDbain} computes the sample size for the robust ANOVA using Bayes factors implemented in \code{bain} with informative hypothesis. \code{SSDbain} can be found at \url{https://github.com/Qianrao-Fu/SSDbain}. The details of the approach can be found in Fu, Moerbeek, & Hoijtink (unpublished). \bold{Users are required to install package \code{bain} that can be found at \url{https://informative-hypotheses.sites.uu.nl/software/bain/} before using \code{SSDbain} package}.}
\usage{
Res<-SSDANOVA_robust(hyp1,hyp2,f1,f2,skews,kurts,var,BFthresh,eta,T,seed)}
\arguments{
\item{hyp1}{strings that specify the former hypothesis of a pair of interested hypotheses. For example, if H0: mu1=mu2=mu3 vs H1: mu1>mu2>mu3, then \code{hyp1="mu1=mu2=mu3"}. }
\item{hyp2}{strings that specify the latter hypothesis. For example, as the previous example, then \code{hyp2="mu1>mu2>mu3"}. Note: when hyp2 is an unconstrained hypothesis, then \code{hyp2="Ha"}; when hyp2 is the complement hypothesis of hyp1, then \code{hyp2="Hc"}.}
\item{f1}{parameter used to specify the populations corresponding to hyp1. It can be specified as a number denoting the effect size or a vector denoting the means for each of group. \code{f1} is Cohen's f=sigma1/sigma, where sigma1 denotes the standard deviation of the means, and sigma denotes the common within-group standard deviation.}
\item{f2}{parameter used to specify the populations corresponding to hyp2. It can be specified as a number denoting the effect size or a vector denoting the means for each of group.}
\item{var}{a vector that specifies the within-group variances for each of the group. The user need to input Cohen's f or a vector of means and the variances for each group.}
\item{skews}{vector that specifies the skewness of each population. The default setting is \code{skews=c(0,0,0,...,0)}, which corrsponds to the normal distributions.}
\item{kurts}{vector that specifies the kurtosis of each population. The default setting is \code{kurts=c(0,0,0,...,0)}, which corresponds to the normal distributions.}
\item{\code{BFthresh}}{numeric value not less than 1 that specifies the required size of the Bayes factor. The default setting is \code{BFthresh=3}.}
\item{\code{eta}}{a numeric value that specifies the probability that the Bayes factor is larger than the \code{BFthresh} if either of the competing hypotheses is true. The default setting is \code{eta=0.80}.}
\item{\code{T}}{a positive integer that specifies the number of data sets sampled from the populations corresponding to the two hypotheses of interest. The default setting is \code{T=10000}.}
\item{\code{seed}}{a positive integer that specifies the seed of R's random number generator. The default setting is \code{seed=10}.}
}

\value{
The output resulting from analyses with \code{SSDANOVA_robust} contains:

1) required sample size N

2) the probability P(BFsv>BFthresh|Hs)

3) the probability P(BFvs>BFthresh|Hv)

}

\references{
Fu, Q., Mirjam, M., and Hoijtink, H. (unpublished). Sample size determination for the Bayesian ANOVAs with Informative hypotheses

Gu, X., Mulder, J., and Hoijtink, H. (2017). Approximated adjusted fractional Bayes factors: A general method for testing informative hypotheses. British Journal of Mathematical and Statistical Psychology. \emph{British Journal of Mathematical and Statistical Psychology}, 71(2), 229-261. doi:10.1111/bmsp.12110

Hoijtink, H., Gu, X., and Mulder, J. (2018). Bayesian evaluation of informative hypotheses for multiple populations. \emph{British Journal of Mathematical and Statistical Psychology}. doi:10.1111/bmsp.12145

}
\author{
Qianrao Fu
}


\examples{
#Example 1: if the data may include outliers, and using Cohen's f and within-group variances to specify the populations of interest for robust ANOVA

Res<-SSDANOVA_robust(hyp1="mu1=mu2=mu3",hyp2="mu1>mu2>mu3",f1=0,f2=0.25,skews=c(0,0,0),kurts=c(0,0,0),var=c(1.5,0.75,0.75),BFthresh=3,eta=0.8,T=10000,seed=10)

#Example 2: if the data is skewed or heavy tailed, and using Cohen's f and within-group variances to specify the populations of interest for robust ANOVA.

Res<-SSDANOVA_robust(hyp1="mu1=mu2=mu3",hyp2="mu1>mu2>mu3",f1=0,f2=0.25,skews=c(1.75,1.75,1.75),kurts=c(5.89,5.89,5.89),var=c(1.5,0.75,0.75),BFthresh=3,eta=0.8,T=10000,seed=10)

#Example 3: if the data is heavy tailed, and using means and within-group variances to specify the populations of interest for robust ANOVA.

Res<-SSDANOVA_robust(hyp1="mu1>mu2>mu3",hyp2="Hc",f1=c(7.33,6.13,5.00),f2=c(5.00,7.33,6.13),skews=c(0,0,0),kurts=c(6.94,6.94,6.94),var=c(2.330^2,2.875^2,2.059^2),BFthresh=3,eta=0.8,T=10000,seed=10)

}




