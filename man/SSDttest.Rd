\name{SSDttest}
\alias{SSDttest}
\title{Sample Size Determination for the Bayesian t-test and Bayesian Welch's test}
\description{The function \code{SSDttest} in the \code{R} package \code{SSDbain} computes the sample size for the Bayesian t-test and Welch's test using Bayes factors implemented in \code{bain}. \code{SSDbain} can be found at \url{https://github.com/Qianrao-Fu/SSDbain}. The details of the approach can be found in Fu et al. (unpublished). Users are required to install package \code{bain} that can be found at \url{https://informative-hypotheses.sites.uu.nl/software/bain/} from CRAN before using \code{SSDbain} package.}
\usage{
SSDttest(type='equal',Population_mean=c(0.5,0),var=NULL,BFthresh=3,eta=0.8,Hypothesis='two-sided',T=10000)}
\arguments{
\item{\code{type}}{a string that specifies the type of the test. If \code{type='equal'}, the t-test is used; if \code{type='unequal'}, Welch's test is used.}

\item{\code{Population_mean}}{vector of length 2 specifying the population means of the two groups under the alternative hypothesis H1 or H2. The default setting is \code{Population_mean=c(0.5,0)} when the effect size is d=0.5. Note that, if \code{var=NULL} and the population mean in Group 2 equals  0, the population mean in Group 1 is identical to Cohen's d.}

\item{\code{var}}{vector of length 2 giving the two within-group variances. If \code{type='equal'}, the default is \code{var=c(1,1)}; if \code{type='unequal'}, the default is \code{var=c(4/3,2/3)}. Of course, any values of the variances can be used as input for the argument \code{var}.}

\item{\code{BFthresh}}{a numeric value that specifies the magnitude of the Bayes factor, e.g. 1, 3, 5 or 10. The default setting is \code{BFthresh=3}.}
\item{\code{eta}}{a numeric value that specifies the probability that the Bayes factor is larger than the \code{BFthresh} if either the null hypothesis or the alternative hypothesis is true, e.g., 0.80, 0.90. That is, P(BF0i>BFthresh|H0) and P(BFi0>BFthresh|Hi) reach at least eta. The default setting is \code{eta=0.80}.}
\item{\code{Hypothesis}}{a string that specifies the hypothesis. \code{Hypothesis='two-sided'} when the competing hypotheses are H0:mu1=mu2, H1:mu1,mu2; \code{Hypothesis='one-sided'} when the competing hypotheses are H0:mu1=mu2, H2: mu1>mu2.}
\item{\code{T}}{a postive integer that indicating the number of data sets sampled from the null and alternative populations to determine the required sample size. The default setting is \code{T=10000}. Note: a minimum value of \bold{10000} is recommended.}
}
\value{
The output resulting from analyses with \code{SSDttest} contains:

1) required sample size N

2) the probability P(BF0i>BFthresh|H0)

3) the probability P(BFi0>BFthresh|Hi)

}

\references{
Fu, Q., Mulder, Hoijtink, H. and Mirjam, M. (unpublished). Sample size determination for the Bayesian t-test and Welch's test using the approximate adjusted fractional Bayes factor. doi:10.31234/osf.io/2v4js

Gu, X., Mulder, J., and Hoijtink, H. (2017). Approximated adjusted fractional Bayes factors: A general method for testing informative hypotheses. British Journal of Mathematical and Statistical Psychology. \emph{British Journal of Mathematical and Statistical Psychology}, 71(2), 229-261. doi:https://doi.org/10.1111/bmsp.12110

Hoijtink, H., Gu, X., and Mulder, J. (2018). Bayesian evaluation of informative hypotheses for multiple populations. \emph{British Journal of Mathematical and Statistical Psychology}, 72(2), 219-243. doi:https://doi.org/10.1111/bmsp.12145
}

\author{
Qianrao Fu
}

\examples{
#Example 1: for the two-sided t-test, effect size d=0.5, BFthresh=3, and eta=0.8.

SSDttest(type='equal',Population_mean=c(0.5,0),var=NULL,BFthresh=3,eta=0.8,Hypothesis='two-sided',T=10000)


#Example 2: for the one-sided t-test, effect size d=0.2, BFthresh=1, and eta=0.9.

SSDttest(type='equal',Population_mean=c(0.2,0),var=NULL,BFthresh=1, eta=0.9,Hypothesis='one-sided',T=10000)


#Example 3: for the two-sided Welch's test, effect size d=0.5, BFthresh=3,and eta=0.8.

SSDttest(type='unequal',Population_mean=c(0.5,0),var=NULL,BFthresh=3,eta=0.8,Hypothesis='two-sided',T=10000)

}
