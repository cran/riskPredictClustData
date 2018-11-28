\name{genSimDataGLMEM}
\alias{genSimDataGLMEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate simulated data from logistic mixed effects model based on the AMD data
}
\description{
Generate simulated data from logistic mixed effects model based on the AMD data.
}
\usage{
genSimDataGLMEM(
  nSubj = 131, 
  beta0 = -6, 
  sd.beta0i = 1.58, 
  beta1 = 1.58, 
  beta2 = -3.95, 
  beta3 = 3.15, 
  beta4 = 2.06, 
  beta5 = 0.51, 
  beta6 = 1.47, 
  beta7 = 3.11, 
  p.smkcur = 0.08, 
  p.inieye31 = 0.44, 
  p.inieye32 = 0.42, 
  p.inieye41 = 0.12, 
  p.inieye42 = 0.11, 
  sd.lncalorc = 0.33)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nSubj}{
integer. Number of subjects. Each subject would have data for 2 eyes.
}
  \item{beta0}{
mean of intercept \eqn{\beta_{0i}}{beta_{0i}}, 
which is assumed random and follows normal distribution
\eqn{N(\beta_0, \sigma^2_{\beta})}{
N(beta_0, sigma^2_{beta}).
}
}
  \item{sd.beta0i}{
standard deviation \eqn{\sigma^2_{\beta}}{sigma^2_{beta}} of the random intercept 
\eqn{\beta_{0i}}{beta_{0i}}.
}
  \item{beta1}{
slope for the binary covariate \eqn{cursmk}{cursmk} (current smoking status).
\eqn{cursmk=1}{cursmk=1} indicates current smokers.
\eqn{cursmk=0}{cursmk=0} indicates past smokers or never smokers.
}
  \item{beta2}{
slope for the continuous mean-centered covariate \eqn{lncalor_c}{lncalor_c}.
}
  \item{beta3}{
slope for the binary covariate \eqn{inieye3}{inieye3} indicating if an eye of a subject has initial grade equal to 3.
\eqn{inieye3=1}{inieye3=1} indicates the eye has initial grade equal to 3.
}
  \item{beta4}{
slope for the binary covariate \eqn{inieye4}{inieye4} indicating if an eye of a subject has initial grade equal to 4.
\eqn{inieye4=1}{inieye4=1} indicates the eye has initial grade equal to 4.
}
  \item{beta5}{
slope for the binary covariate \eqn{rtotfat_1}{rtotfat_1} indicating if the subject's
total fat intake is in the 2nd quartile of total fat intake.
\eqn{rtotfat_1=1}{rtotfat_1=1} indicates the subject is in the 2nd quartile.
}
  \item{beta6}{
slope for the binary covariate \eqn{rtotfat_2}{rtotfat_2} indicating if the subject's
total fat intake is in the 3rd quartile of total fat intake.
\eqn{rtotfat_2=1}{rtotfat_2=1} indicates the subject is in the 3rd quartile.
}
  \item{beta7}{
slope for the binary covariate \eqn{rtotfat_3}{rtotfat_3} indicating if the subject's
total fat intake is in the 4th quartile of total fat intake.
\eqn{rtotfat_3=1}{rtotfat_3=1} indicates the subject is in the 4th quartile.
}
  \item{p.smkcur}{
proportion of current smokers.
}
  \item{p.inieye31}{
proportion of left eye having inital grade equal to 3. 
}
  \item{p.inieye32}{
proportion of right eye having inital grade equal to 3. 
}
  \item{p.inieye41}{
proportion of left eye having inital grade equal to 4. 
}
  \item{p.inieye42}{
proportion of right eye having inital grade equal to 4. 
}
  \item{sd.lncalorc}{
standard deviation for \eqn{lncalor_c}{lncalor_c}.
}
}
\details{
We generate simulated data set from the following generalized linear mixed effects model:
\deqn{\log\left(\frac{p_{ij}}{(1-p_{ij})}\right)=\beta_{0i}+\beta_1 smkcur_i+
\beta_2 lncalor_{ci} + \beta_3 inieye3_{ij} + \beta_4 inieye4_{ij} 
+\beta_5 rtotfat_{1i} +\beta_6 rtotfat_{2i} + \beta_7 rtotfat_{3i},}{
logit(p_{ij})=beta_{0i}+beta_1 smkcur_i+
beta_2 lncalor_{ci} + beta_3 inieye3_{ij} + beta_4 inieye4_{ij} 
+\beta_5 rtotfat_{1i} +\beta_6 rtotfat_{2i} + \beta_7 rtotfat_{3i},
}
\eqn{i=1,\ldots, N, j=1, 2}{i=1, ..., N, j=1, 2},
\eqn{\beta_{0i}\sim N\left(\beta_0, \sigma^2_{\beta}\right).}{
beta_{0i}~ N(beta_0, sigma^2_{\beta}).
}
}
\value{
A data frame with 8 columns: cid, subuid, prog, smkcur, lncalorc, inieye3, inieye4, and rtotfat,
where cid is the subject id, subuid is the unit id, and prog is the progression status.
\eqn{prog=1}{prog=1} indicates the eye is progressed.
\eqn{prog=0}{prog=0} indicates the eye is not progressed.
There are \code{nSubj*2} rows. The first \code{nSubj} rows
are for the left eyes and the second \code{nSubj} rows are for the right eyes.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")
