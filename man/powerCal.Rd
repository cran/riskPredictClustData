\name{powerCal}
\alias{powerCal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the power for testing \eqn{\delta=0}{delta=0}
}
\description{
Calculate the power for testing \eqn{\delta=0}{delta=0}.
}
\usage{
powerCal(
  nSubj, 
  mu1, 
  triangle, 
  rho, 
  rho11, 
  rho22, 
  rho12, 
  p11, 
  p10, 
  p01, 
  alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nSubj}{
integer. number of subjects to be generated. Assume each subject has two observations.
}
  \item{mu1}{
\eqn{\mu_1=H(Y)-H(Y_c)}{mu_1=H(Y)-H(Y_c)} is the difference between probit transformation 
\eqn{H(Y)}{H(Y)} and probit-shift alternative \eqn{H(Y_c)}{H(Y_c)},
where \eqn{Y}{Y} is the prediction score of a randomly selected progressing
subunit, and \eqn{Y_c}{Y_c} is the counterfactual random variable
obtained if each subunit that had progressed actually had not progressed.
}
  \item{triangle}{
the difference of the expected value the the extended Mann-Whitney U statistics
between two prediction rules, i.e., \eqn{\triangle = \eta^{(1)}_c - \eta^{(2)}_c}{
triangle = eta^{(1)}_c - eta^{(2)}_c
}
}
  \item{rho}{
\eqn{\rho=corr\left(H\left(Z_{ij}\right), H\left(Z_{k\ell}\right)\right)}{
rho=corr(H(Z_{ij}), H(Z_{k ell}))
}, where \eqn{H=\Phi^{-1}}{H=Phi^{-1}} is the probit transformation.
}
  \item{rho11}{
\eqn{\rho_{11}=corr\left(H_{ij}^{(1)}, H_{i\ell}^{(1)}\right)}{
rho_{11}=corr(H_{ij}^{(1)}, H_{i ell}^{(1)})
}, where \eqn{H=\Phi^{-1}}{H=Phi^{-1}} is the probit transformation.

}
  \item{rho22}{
\eqn{\rho_{22}=corr\left(H_{ij}^{(2)}, H_{i\ell}^{(2)}\right)}{
rho_{22}=corr(H_{ij}^{(2)}, H_{i ell}^{(2)})
}, where \eqn{H=\Phi^{-1}}{H=Phi^{-1}} is the probit transformation.

}
  \item{rho12}{
\eqn{\rho_{12}=corr\left(H_{ij}^{(1)}, H_{i\ell}^{(2)}\right)}{
rho_{12}=corr(H_{ij}^{(1)}, H_{i ell}^{(2)})
}, where \eqn{H=\Phi^{-1}}{H=Phi^{-1}} is the probit transformation.

}
  \item{p11}{
\eqn{p_{11}=Pr(\delta_{i1}=1 \& \delta_{i2}=1)}{
p_{11}=Pr(delta_{i1}=1 & delta_{i2}=1)
}, where \eqn{\delta_{ij}=1}{delta_{ij}=1} if the \eqn{j}{j}-th subunit of the 
\eqn{i}{i}-th cluster has progressed.
}
  \item{p10}{
\eqn{p_{10}=Pr(\delta_{i1}=1 \& \delta_{i2}=0)}{
p_{11}=Pr(delta_{i1}=1 & delta_{i2}=0)
}, where \eqn{\delta_{ij}=1}{delta_{ij}=1} if the \eqn{j}{j}-th subunit of the 
\eqn{i}{i}-th cluster has progressed.

}
  \item{p01}{
\eqn{p_{01}=Pr(\delta_{i1}=0 \& \delta_{i2}=1)}{
p_{11}=Pr(delta_{i1}=0 & delta_{i2}=1)
}, where \eqn{\delta_{ij}=1}{delta_{ij}=1} if the \eqn{j}{j}-th subunit of the 
\eqn{i}{i}-th cluster has progressed.

}
  \item{alpha}{
type I error rate
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
the power
}
\references{
Rosner B, Qiu W, and Lee MLT.
Assessing Discrimination of Risk Prediction Rules in a Clustered Data Setting.
Lifetime Data Anal. 2013 Apr; 19(2): 242-256.
}
\author{
Bernard Rosner <stbar@channing.harvard.edu>, 
Weiliang Qiu <Weiliang.Qiu@gmail.com>, 
Meiling Ting Lee <MLTLEE@umd.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

 

set.seed(1234567)
mu1 = 0.8

power = powerCal(nSubj = 30, mu1 = mu1, 
  triangle = 0.05, rho = 0.93, rho11 = 0.59, rho22 = 0.56, rho12 = 0.52,
  p11 = 0.115, p10 = 0.142, p01 = 0.130, alpha = 0.05)

print(power)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
