\name{riskPredict}
\alias{riskPredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assessing risk prediction performance for clustered data
}
\description{
Assessing risk prediction performance for clustered data.
}
\usage{
riskPredict(frame, alpha=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frame}{
A data frame with 4 columns: cid, subuid, status, and score.
\code{cid} indicates cluster id; \code{subuid} indicates unit ID within a cluster;
\code{status=1} indicates an eye is progressed;
\code{status=0} indicates an eye is not progressed;
\code{score} represents the risk score.
}
  \item{alpha}{
numeric. confidence level for \eqn{eta_c}.
}
}
\details{
To obtain 95\% confidence interval of \eqn{\eta_c}{eta_c},
we first obtain 95\% confidence interval \eqn{[c_1, c_2]}{[c_1, c_2]}
for 
\eqn{\Phi^{-1}(\eta_c)}{Phi^{-1}(eta_c)}, then
transform back:
\eqn{[\Phi(c_1), \Phi(c_2)]}{[Phi(c_1), Phi(c_2)]}.
}
\value{
A list of 6 elements:
  \item{stat}{the test statistics 
\eqn{\hat{\eta}_c^{(1)}}
{hat{eta}_c^{(1)}} based on the prediction rule.}
  \item{se.stat}{standard error of the test statistic under the 
null hypothesis.}
  \item{z}{z score \code{z=(stat - 0.5)/se.stat}}
  \item{pval}{p-value of the test}
  \item{rho}{
correlation between \eqn{H(Z_{ij}) and H(Z_{i \ell})}{
H(Z_{ij}) and H(Z_{i ell})
}
}
  \item{mu.hat}{
estimated \eqn{\mu}.
}
  \item{theta.hat}{
estimated \eqn{\theta}.
}
  \item{theta.c.hat}{
estimated \eqn{\theta_c}.
}
  \item{E.stat.Ha}{
expectation of \eqn{\hat{\eta}_c} under
the alternative hypothesis.
}
  \item{se.stat.Ha}{
standard error for \eqn{\hat{\eta}_c} under
the alternative hypothesis.
}
  \item{CIlow}{
lower confidence limit  for \eqn{\eta_c}.
}
  \item{CIupp}{
upper confidence limit  for \eqn{\eta_c}.
}
  \item{datHk}{
A nSubj by 2 matrix of probit transformed risk scores by using only
the first 2 observations of each subject.
}
  \item{ci}{
the vector of \eqn{c_i}{c_i}, the number of progressing subunits for the \eqn{i}{i}-th subject.
}
  \item{di}{
the vector of \eqn{d_i}{d_i}, the number of non-progressing subunits for the \eqn{i}{i}-th subject.
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

tt1 = getScore(fmla = prog~smkcur+lncalorc+inieye3+inieye4+factor(rtotfat), 
  cidVar = "cid", subuidVar = "subuid", statusVar = "prog", 
  datFrame = datFrame, mycorstr = "exchangeable",
  verbose = FALSE)
myframe1=tt1$frame

print(dim(myframe1))
print(myframe1[1:3,])

res1 = riskPredict(myframe1)
print(names(res1))
print(res1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")
