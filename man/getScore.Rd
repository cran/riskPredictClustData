\name{getScore}
\alias{getScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get data frame for the function riskPredict
}
\description{
Get data frame for the function riskPredict.
}
\usage{
getScore(fmla, cidVar, subuidVar, statusVar, datFrame, mycorstr = "exchangeable",
  verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fmla}{A formula object for the function \code{gee}}
  \item{cidVar}{character. Phenotype variable name for cluster id}
  \item{subuidVar}{character. Phenotype variable name for unit id}
  \item{statusVar}{character. Phenotype variable name for progression status}
  \item{datFrame}{
A data frame with at least 3 columns corresponding to \code{cid} (indicated by \code{cidVar}), 
\code{subuid} (indicated by \code{subuidVar}), 
\code{status} (indicated by \code{statusID}).
\code{cid} indicates cluster id; \code{subuid} indicates unit ID within a cluster;
\code{status=1} indicates an eye is progressed;
\code{status=0} indicates an eye is not progressed.
}
  \item{mycorstr}{character. indicates correlation structure. see the manual for the function \code{gee} in the R library \code{gee}}
  \item{verbose}{
logical. indicating if summary of gee results should be printed out.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with two elements: \code{frame} and \code{gee.obj}.
\code{frame} is a data frame with at least 4 columns: cid, subuid, status, and score.
\code{cid} indicates cluster id; \code{subuid} indicates unit ID within a cluster;
\code{status=1} indicates an eye is progressed;
\code{status=0} indicates an eye is not progressed;
\code{score} represents the risk score.

\code{gee.obj} is the object returned by \code{gee} function.
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

gee.obj=tt1$gee.obj
print(summary(gee.obj))

print(dim(myframe1))
print(myframe1[1:3,])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")
