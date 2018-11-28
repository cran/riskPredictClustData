# modified on June 13, 2018
#  (1) remove row names after create data frame in simulated data
#
# modified on March 7, 2018
#  (1) sort by cluster id
#  (2) change 'sid' to 'cid'
#  (3) change 'uid' to 'subuid'
#
# modified on Dec. 28, 2017
#  (1) set the default values: beta0=-6
#
# simulate data from logistic mixed effects model 
#
#  \log\left(\frac{p_{ij}}{(1-p_{ij})}\right)=&\beta_{0i}+\beta_1 smkcur_i+
# \beta_2 lncalor_{ci} + \beta_3 inieye3_{ij} + \beta_4 inieye4_{ij} \\
# &+\beta_5 rtotfat_{1i} +\beta_6 rtotfat_{2i} + \beta_7 rtotfat_{3i},
# i=1,\ldots, N, j=1, 2,\\
#  \beta_{0i}\sim & \N\left(\beta_0, \sigma^2_{\beta}\right),

genSimDataGLMEM=function(nSubj=131, beta0 = -6, sd.beta0i = 1.58, 
                          beta1=1.58, beta2=-3.95, beta3=3.15, beta4=2.06,
                          beta5=0.51, beta6=1.47, beta7=3.11, 
                          p.smkcur=0.08, p.inieye31=0.44, p.inieye32=0.42,
                          p.inieye41=0.12, p.inieye42=0.11, sd.lncalorc=0.33)
{
  # generate intercept
  beta0i=rnorm(nSubj, mean=beta0, sd=sd.beta0i)
  
  # generate current smoking status
  smkcuri=sample(c(1,0), size=nSubj, prob=c(p.smkcur, 1-p.smkcur), replace=TRUE)
  
  # generate lncalorc
  lncalorc = rnorm(nSubj, mean=0, sd=sd.lncalorc)
  
  # generate inieye3_1 (left eye)
  inieye31 = sample(c(1,0), size=nSubj, prob=c(p.inieye31, 1-p.inieye31),
                    replace = TRUE)
  
  # generate inieye3_2 (right eye)
  inieye32 = sample(c(1,0), size=nSubj, prob=c(p.inieye32, 1-p.inieye32),
                    replace = TRUE)
  
  # generate inieye4_1 (left eye)
  inieye41 = sample(c(1,0), size=nSubj, prob=c(p.inieye41, 1-p.inieye41),
                    replace = TRUE)
  
  # generate inieye4_2 (right eye)
  inieye42 = sample(c(1,0), size=nSubj, prob=c(p.inieye42, 1-p.inieye42),
                    replace = TRUE)
  
  # generate rtotfat quartiles
  rtotfat4=sample(c(1,2,3,4), size=nSubj, prob=c(1/4,1/4,1/4,1/4),
                  replace = TRUE)

  rtotfat42=as.numeric(rtotfat4==2)
  rtotfat43=as.numeric(rtotfat4==3)
  rtotfat44=as.numeric(rtotfat4==4)

  # generate outcome for left eye
  a1 = beta0i+beta1*smkcuri+beta2*lncalorc+beta3*inieye31+beta4*inieye41+
      beta5*rtotfat42+beta6*rtotfat43+beta7*rtotfat44
  ea1 = exp(a1)
  p1 = ea1/(1+ea1)
  y1 = unlist(lapply(1:nSubj, function(i) {
    tti=sample(c(1,0), size=1, prob=c(p1[i], 1-p1[i]), replace=TRUE)
    return(tti)
  }))

  a2 = beta0i+beta1*smkcuri+beta2*lncalorc+beta3*inieye32+beta4*inieye42+
      beta5*rtotfat42+beta6*rtotfat43+beta7*rtotfat44
  ea2 = exp(a2)
  p2 = ea2/(1+ea2)
  y2 = unlist(lapply(1:nSubj, function(i) {
    tti=sample(c(1,0), size=1, prob=c(p2[i], 1-p2[i]), replace=TRUE)
    return(tti)
  }))
 
  # construct data frame
  cid=c(1:nSubj, 1:nSubj)
  subuid=c(rep(1, nSubj), rep(2, nSubj))
  prog=c(y1, y2)
  smkcurVec=c(smkcuri, smkcuri)
  lncalorcVec=c(lncalorc, lncalorc)
  inieye3Vec=c(inieye31, inieye32)
  inieye4Vec=c(inieye41, inieye42)
  rtotfatVec=c(rtotfat4, rtotfat4)
  
  datFrame=data.frame(cid=cid, subuid=subuid, prog=prog, smkcur=smkcurVec, lncalorc=lncalorcVec,
                      inieye3=inieye3Vec, inieye4=inieye4Vec,
                      rtotfat=rtotfatVec)
  datFrame.s=datFrame[order(datFrame$cid, datFrame$subuid),]
  rownames(datFrame.s)=NULL
  # need to use print(dataFrame.s, row.names=FALSE)
  
  invisible(datFrame.s)
  
}

# test
#datFrame=genSimDataGLMEM(nSubj=131, beta0 = -6, sd.beta0i = 1.58, 
#                          beta1=1.58, beta2=-3.95, beta3=3.15, beta4=2.06,
#                          beta5=0.51, beta6=1.47, beta7=3.11, 
#                          p.smkcur=0.08, p.inieye31=0.44, p.inieye32=0.42,
#                          p.inieye41=0.12, p.inieye42=0.11, sd.lncalorc=0.33)
