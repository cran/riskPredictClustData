\name{powerCalData}
\alias{powerCalData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the power for testing \eqn{\delta=0}{delta=0} based on a dataset
}
\description{
Calculate the power for testing \eqn{\delta=0}{delta=0} based on a dataset.
}
\usage{
powerCalData(
  nSubj, 
  triangle, 
  frame,
  alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nSubj}{
integer. number of subjects to be generated. Assume each subject has two observations.
}
  \item{triangle}{
the difference of the expected value the the extended Mann-Whitney U statistics
between two prediction rules, i.e., \eqn{\triangle = \eta^{(1)}_c - \eta^{(2)}_c}{
triangle = eta^{(1)}_c - eta^{(2)}_c
}
}
  \item{frame}{
A data frame with 5 columns: cid, subuid, status, score1, and score2.
\code{cid} indicates cluster id; \code{subuid} indicates unit ID within a cluster;
\code{status=1} indicates an eye is progressed;
\code{status=0} indicates an eye is not progressed;
\code{score1} represents the score based on prediction rule 1.
\code{score2} represents the score based on prediction rule 2.
}
  \item{alpha}{
type I error rate
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with 11 elements. 
\item{power}{the esstimated power}
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
p_{10}=Pr(delta_{i1}=1 & delta_{i2}=0)
}, where \eqn{\delta_{ij}=1}{delta_{ij}=1} if the \eqn{j}{j}-th subunit of the 
\eqn{i}{i}-th cluster has progressed.

}
  \item{p01}{
\eqn{p_{01}=Pr(\delta_{i1}=0 \& \delta_{i2}=1)}{
p_{01}=Pr(delta_{i1}=0 & delta_{i2}=1)
}, where \eqn{\delta_{ij}=1}{delta_{ij}=1} if the \eqn{j}{j}-th subunit of the 
\eqn{i}{i}-th cluster has progressed.

}
  \item{p00}{
\eqn{p_{00}=Pr(\delta_{i1}=0 \& \delta_{i2}=0)}{
p_{00}=Pr(delta_{i1}=0 & delta_{i2}=0)
}, where \eqn{\delta_{ij}=1}{delta_{ij}=1} if the \eqn{j}{j}-th subunit of the 
\eqn{i}{i}-th cluster has progressed.

}

  \item{mu1}{
\eqn{\mu_1=H(Y)-H(Y_c)}{mu_1=H(Y)-H(Y_c)} is the difference between probit transformation 
\eqn{H(Y)}{H(Y)} and probit-shift alternative \eqn{H(Y_c)}{H(Y_c)} for the first prediction score,
where \eqn{Y}{Y} is the prediction score of a randomly selected progressing
subunit, and \eqn{Y_c}{Y_c} is the counterfactual random variable
obtained if each subunit that had progressed actually had not progressed.
}

  \item{mu2}{
\eqn{\mu_2=H(Y)-H(Y_c)}{mu_2=H(Y)-H(Y_c)} is the difference between probit transformation 
\eqn{H(Y)}{H(Y)} and probit-shift alternative \eqn{H(Y_c)}{H(Y_c)} for the second prediction score,
where \eqn{Y}{Y} is the prediction score of a randomly selected progressing
subunit, and \eqn{Y_c}{Y_c} is the counterfactual random variable
obtained if each subunit that had progressed actually had not progressed.
}

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

datFrame = genSimDataGLMEM(nSubj = 30, beta0 = -6, sd.beta0i = 1.58, 
                          beta1 = 1.58, beta2 = -3.95, beta3 = 3.15, beta4 = 2.06,
                          beta5 = 0.51, beta6 = 1.47, beta7 = 3.11, 
                          p.smkcur = 0.08, p.inieye31 = 0.44, p.inieye32 = 0.42,
                          p.inieye41 = 0.12, p.inieye42 = 0.11, sd.lncalorc = 0.33)

print(dim(datFrame))
print(datFrame[1:2,])

# prediction rule 1
tt1 = getScore(fmla = prog~smkcur+lncalorc+inieye3+inieye4+factor(rtotfat), 
  cidVar = "cid", subuidVar = "subuid", statusVar = "prog", 
  datFrame = datFrame, mycorstr = "exchangeable",
  verbose = FALSE)
myframe1=tt1$frame

print(dim(myframe1))
print(myframe1[1:3,])

####
# prediction rule 2
tt2 = getScore(fmla = prog~smkcur+lncalorc+inieye3+inieye4, 
  cidVar = "cid", subuidVar = "subuid", statusVar = "prog", 
  datFrame = datFrame, mycorstr = "exchangeable",
  verbose = FALSE)
myframe2=tt2$frame

print(dim(myframe2))
print(myframe2[1:3,])

# combine scores from two prediction rules
myframe12=myframe1[, c("cid", "subuid", "status")]
myframe12$score1=myframe1$score
myframe12$score2=myframe2$score
print(dim(myframe12))
print(myframe12[1:3,])


res = powerCalData(nSubj = 30, triangle = 0.05, frame=myframe12, alpha = 0.05)

print(res)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
