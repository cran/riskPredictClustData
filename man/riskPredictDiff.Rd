\name{riskPredictDiff}
\alias{riskPredictDiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Difference of two risk prediction rules for clustered data
}
\description{
Difference of two risk prediction rules for clustered data.
}
\usage{
riskPredictDiff(frame, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frame}{
A data frame with 5 columns: cid, subuid, status, score1, and score2.
\code{cid} indicates cluster id; \code{subuid} indicates unit ID within a cluster;
\code{status=1} indicates an eye is progressed;
\code{status=0} indicates an eye is not progressed;
\code{score1} represents the score based on prediction rule 1.
\code{score2} represents the score based on prediction rule 2.
}
  \item{alpha}{
numeric. The confidence level.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list of 7 elements:
  \item{diff}{the difference of test statistics 
\eqn{\hat{\eta}_c^{(1)}-\hat{\eta}_c^{(2)}}
{hat{eta}_c^{(1)}-hat{eta}_c^{(2)}} based on the 2 prediction rules.}
  \item{se.diff}{standard error of the difference under the null hypothesis.}
  \item{z}{z score \code{z=diff/se.diff}}
  \item{pval}{p-value of the test}
  \item{res1}{
output object of the function \code{riskPredict} for prediction rule 1.
}
  \item{res2}{
output object of the function \code{riskPredict} for prediction rule 2.
}
  \item{rhoVec}{
A vector of 4 correlations:
\eqn{\rho=cov(H_{ij}^{(1)}, H_{ij}^{(2)})}{
rho=cov(H_{ij}^{(1)}, H_{ij}^{(2)})
},
\eqn{\rho_{11}=cov(H_{ij}^{(1)}, H_{it}^{(1)})}{
rho_{11}=cov(H_{ij}^{(1)}, H_{it}^{(1)})
},
\eqn{\rho_{22}=cov(H_{ij}^{(2)}, H_{it}^{(2)})}{
rho_{22}=cov(H_{ij}^{(2)}, H_{it}^{(2)})
},
and \eqn{\rho_{12}=cov(H_{ij}^{(1)}, H_{it}^{(2)})}{
rho_{22}=cov(H_{ij}^{(1)}, H_{it}^{(2)})
}
}
  \item{E.diff.Ha}{
expectation of the difference under the alternative hypothesis.
}
  \item{se.diff.Ha}{
standard error of the difference under the alternative hypothesis.
}
  \item{CIlow.diff}{
Lower confidence limit.
}
  \item{CIup.diff}{
Upper confidence limit.
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

####
resDiff = riskPredictDiff(frame=myframe12)
print(names(resDiff))
print(resDiff)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")
