# simulation study for clustered Wilcoxon test statistics
# created on June 13, 2018
#  (1) output regression coefficient of GEE in 'getScore' function
#  (2) write a power calculation function so that the user can input pilot data to obtain parameter estimates


# created on March 7, 2018 by Weiliang Qiu
#  (1) change 'value' to 'score' 
#  (2) change 'cID' to 'cid'
#  (3) change 'subuid' to 'subuid'
#
# v0.2.0 created on May 3, 2017 by Weiliang Qiu
#  (1) simplify the R code
#
# v0.1.9 created on April 25, 2017 by Weiliang Qiu
#  (1) add function print.riskPredict and print.riskPredictDiff
#
# v0.1.8 created on April 12, 2017 by Weiliang Qiu
#  (1) revised gee function to allow input variable number for cid
#  (2) add GEE wrapper function
#
# v0.1.7 created on Feb. 15, 2017 by Weiliang Qiu
#  (1) use general data format: 
#   cluster ID (cid), unitID (subuid), progression status (status), 
#   measurement (value)
#
# v0.1.6 created on Feb. 1, 2017 by Weiliang Qiu
#  (1) replace "rep=TRUE" by "replace=TRUE"
#
# v0.1.5 created on Oct. 18, 2011 by Weiliang Qiu
#  (1) added function powerCal
# v0.1.4 created on Jan. 3, 2010 by Weiliang Qiu
#  (1) revised function 'Tstat.func' to allow missing scores
#
# v0.1.3 created on Nov. 23, 2009 by Weiliang Qiu
#  (1) added function 'getMu2' to calculate mu2, based on 
#      mu1, rho11, rho22, triangle, and deltas
#
# v0.1.2 created on Nov. 11, 2009 by Weiliang Qiu
#  (1) added function 'genDelta2'
#
# v0.1.1 created on September 25, 2009 by Weiliang Qiu
#  (1) fixed a few bugs in formulas
#
# v0.1.0 created on July 9, 2009 by Weiliang Qiu
#   (1) fixed a bug in the function 'genDat1' (delta without definition)
#   (2) added function to generate data set containing
#       scores for two prediction rules
#
# v0.0.9 created on July 9, 2009 by Weiliang Qiu
#  (1) fixed some wrong formula in calculating variance
#
# v0.0.8 created on July 2, 2009 by Weiliang Qiu
#  (1) add function to calculate Var(theta^{(1)}_c-theta^{(2)}_c)
#
# v0.0.7 created on June 16, 2009 by Weiliang Qiu
#  (1) add function to generate simulated data set
#
# v0.0.6 created on June 11, 2009 by Weiliang Qiu
#  (1) the definition of test statistic changed, hence all formula
#      will change
#
# v0.0.5 created on March 4, 2009 by Weiliang Qiu
#  (1) added functions to obtain empirical CDF
#  (2) added the function 'calPval' to calculate two-sided 
#      p-value based on a Z score
#
# v0.0.4 created on Jan. 5, 2009 by Weiliang Qiu
#  (1) revise the parameter estimates 'paraEst.func'
#  (2) delete functions 'cal.err.func' and 'cal.err.func2'
#  (3) added functions 'cal.TypeIerr.func', 'cal.power.func' and
#      'cal.Thoe.power.func'
#  (4) added output 'paraEstVec' to 'cal.power.func'
#
# v0.0.3 created on Dec. 23, 2008 by Weiliang Qiu
#  (1) added function to calculate log-likelihood function
#
# v0.0.2 created on Dec. 22, 2008 by Weiliang Qiu
#   (1) added parameter estimate functions for rough estimate
#   (2) added simulation function 'cal.err.func' to check type I error and 
#       power
#
# created on Dec. 12, 2008, by Weiliang Qiu


# U function
# U(a)=1 if a<0
# U(a)=1/2 if a=0
# U(a)=0 if a>0
Ufunc<-function(x, y)
{
  a<-x-y
  if(a<0)
  { 
    return(1)
  } else if (a==0) {
    return(1/2)
  } else {
    return(0)
  }
}

# frame - data frame with columns: cid, subuid, status, score
estMuFunc=function(status, score)
{

  Zk<-as.numeric(score)
  HZk<-qnorm(eCDF.func(Zk))
  HZk2 = HZk[which(status==1)]
  mu = mean(HZk2, na.rm=TRUE)

  return(mu)
}


# datFrame - data frame with columns: cid, subuid, status, and other covariates
#  fmla - a formula object for GEE formula
# cidVar - variable indicates cluster id
getScore=function(fmla, cidVar, subuidVar, statusVar, datFrame, mycorstr="exchangeable", verbose=FALSE)
{

  datFrame.old=datFrame

  # gee() function  requires columns 'id' and 'status'
  datFrame$id=datFrame[, c(cidVar)]
  datFrame$cid=datFrame[, c(cidVar)]
  datFrame$subuid=datFrame[, c(subuidVar)]
  datFrame$status=datFrame[, c(statusVar)]
  datFrame.old$status=datFrame$status

  # suppress the output of gee
  mylog <- capture.output({
    gee.obj <- gee(formula=fmla, family=binomial,
      data=datFrame, corstr=mycorstr, silent=TRUE)
  })
  if(verbose)
  {
    print(summary(gee.obj))
  }
    
  fitted<-gee.obj$fitted.values

  datFrame.old$score=fitted

  res=list(frame=datFrame.old, gee.obj=gee.obj)
  invisible(res) 
}

# frame - data frame with columns: cid, subuid, status, score
riskPredict=function(frame, alpha=0.05)
{
  # obtain basic quantities
  # won't sort cid
  u.cid=unique(frame$cid)

  nSubj<-length(u.cid)
  # tapply will sort the cid
  nSubUnits<-tapply(frame$cid, frame$cid, length)

  pos=match(u.cid, names(nSubUnits))
  nSubUnits = nSubUnits[pos]

  # number of eyes that have progressed for the i-th subject
  ci<-tapply(frame$status, frame$cid, sum, na.rm=TRUE)

  # number of eyes that have *NOT* progressed for the i-th subject
  di<-tapply(frame$status, frame$cid, function(x) {sum(1-x, na.rm=TRUE)})

  # obtain test statistic
  stat = Tstat.func(frame = frame, u.cid=u.cid, nSubj=nSubj, nSubUnits=nSubUnits)

  # obtain variance of test statistic under H0

  # exclude subjects having only 1 observation
  pos1=which(nSubUnits <2)
  if(length(pos1))
  {
    u.cid2 = u.cid[-pos1]  
  } else {
    u.cid2 = u.cid
  }
  frame2=frame[which(frame$cid %in% u.cid2),]
  # no need to re-order it
  #frame2.s=frame2[order(frame2$cid, frame2$subuid),]

  # use only first 2 observations to calculate rho
  datk<-matrix(NA, nrow=nSubj, ncol=2, byrow=F)
  for(i in 1:nSubj)
  {
    posi=which(frame2$cid==u.cid2[i])
    framei=frame2[posi,]
    datk[i,1]=framei$score[1]
    datk[i,2]=framei$score[2]
  } 

  Zk<-as.numeric(c(datk))
  HZk<-qnorm(eCDF.func(Zk))
  datHk<-matrix(HZk, nrow=nSubj, ncol=2, byrow=F)
  rho11k<-stats::cor(x=datHk[,1], y=datHk[,2], use="complete.obs")
  #rho11k<-stats::cor(datHk[,1], datHk[,2], use="complete.obs")
  #rho11k<-stats::cov(datHk[,1], datHk[,2], use="complete.obs")

  # var under H0
  var.stat = var.Tstat.func(
    theta = 1/2, 
    theta.c = 1/2, 
    rho = rho11k,
    ci = ci, 
    di = di)

  se.stat = sqrt(var.stat)
  z = (stat - 0.5)/se.stat

  pval = calPval(Z = z)


  # var under Ha
  # theta = Phi(mu/sqrt(2))
  # theta.c = Phi(mu/sqrt(2*(1-rho)))
  # mu is the mean of Hij (Hij~ N(mu, 1) when delta_{ij}=1)

  # estimate mu
  mu.hat=estMuFunc(status=frame$status, score=frame$score)
  theta.hat = pnorm(mu.hat/sqrt(2))
  theta.c.hat = pnorm(mu.hat/sqrt(2*(1-rho11k)))
  
  E.stat.Ha = E.Tstat.func(theta=theta.hat, theta.c=theta.c.hat, ci=ci, di=di)

  var.stat.Ha = var.Tstat.func(
    theta = theta.hat, 
    theta.c = theta.c.hat, 
    rho = rho11k,
    ci = ci, 
    di = di)

  se.stat.Ha = sqrt(var.stat.Ha)
  za=qnorm(1-alpha/2)

  # 95% CI
  xi=qnorm(stat)
  se.xi=se.stat.Ha/abs(dnorm(qnorm(E.stat.Ha)))
  CIlow =  pnorm(xi - za*se.xi)
  CIupp =  pnorm(xi + za*se.xi)

  res=list(stat=stat, se.stat=se.stat, z=z, pval=pval,
    rho=rho11k, mu.hat=mu.hat,
    theta.hat = theta.hat, theta.c.hat=theta.c.hat,
    E.stat.Ha = E.stat.Ha, se.stat.Ha=se.stat.Ha,
    CIlow=CIlow, CIupp=CIupp, datHk=datHk,
    ci = ci, di = di)

  class(res)="class.riskPredict"

  invisible(res)
}

# frame - data frame with columns: cid, subuid, status, score1, score2
# score1 - score 1
# score2 - score 2
riskPredictDiff=function(frame, alpha=0.05)
{

  frame1=frame[, c("cid", "subuid", "status", "score1")]
  colnames(frame1)[4]="score"
  frame2=frame[, c("cid", "subuid", "status", "score2")]
  colnames(frame2)[4]="score"

  res1 = riskPredict(frame=frame1, alpha=alpha)
  res2 = riskPredict(frame=frame2, alpha=alpha)

  ci = res1$ci
  di = res1$di

  Tstat1k = res1$stat
  Tstat2k = res2$stat

  diff = Tstat1k - Tstat2k

  ##
  datHk1 = res1$datHk
  datHk2 = res2$datHk

  rho11 = res1$rho
  rho22 = res2$rho

  ttx=c(datHk1[,1], datHk1[,2])
  tty=c(datHk2[,1], datHk2[,2])
 
  rho = stats::cor(x=ttx, y=tty, use="complete.obs")

  ttx=c(datHk1[,1], datHk1[,2])
  tty=c(datHk2[,2], datHk2[,1])
 
  rho12 = stats::cor(x=ttx, y=tty, use="complete.obs")

  #rho = stats::cov(c(datHk1[,1], datHk1[,2]), c(datHk2[,1], datHk2[,2]),use="complete.obs")
  #rho12 = stats::cov(c(datHk1[,1], datHk1[,2]), c(datHk2[,2], datHk2[,1]), use="complete.obs")

  var.diff<-var.Tstat12.func(
    theta.1 = 1/2, 
    theta.c1 = 1/2, 
    theta.2 = 1/2, 
    theta.c2 = 1/2, 
    rho = rho, 
    rho11 = rho11, 
    rho22 = rho22, 
    rho12 = rho12,
    ci = ci,
    di = di)
  se.diff<-sqrt(var.diff)

  z = diff/se.diff

  pval = calPval(Z = z)

  rhoVec=c(rho, rho11, rho22, rho12)
  names(rhoVec)=c("rho", "rho11", "rho22", "rho12")

  # var under Ha
  # theta = Phi(mu/sqrt(2))
  # theta.c = Phi(mu/sqrt(2*(1-rho)))
  # mu is the mean of Hij (Hij~ N(mu, 1) when delta_{ij}=1)

  E.diff.Ha = res1$E.stat.Ha - res2$E.stat.Ha

  var.diff.Ha<-var.Tstat12.func(
    theta.1 = res1$theta.hat, 
    theta.c1 = res1$theta.c.hat, 
    theta.2 = res2$theta.hat, 
    theta.c2 = res2$theta.c.hat, 
    rho = rho, 
    rho11 = rho11, 
    rho22 = rho22, 
    rho12 = rho12,
    ci = ci,
    di = di)

  se.diff.Ha<-sqrt(var.diff.Ha)

  za=qnorm(1-alpha/2)

  # 95% CI
  CIlow.diff =  diff - za*se.diff.Ha
  CIupp.diff =  diff + za*se.diff.Ha

  res=list(diff=diff, se.diff=se.diff, z=z, pval=pval,
    res1=res1, res2=res2, rhoVec=rhoVec,
    E.diff.Ha = E.diff.Ha,
    se.diff.Ha=se.diff.Ha,
    CIlow.diff=CIlow.diff,
    CIupp.diff=CIupp.diff
  )

  class(res) = "class.riskPredictDiff"

  invisible(res)
}

# frame - data frame with columns: cid, subuid, status, score
Tstat.func<-function(frame, u.cid=u.cid, nSubj=nSubj, nSubUnits=nSubUnits)
{

  numer<-0
  denom<-0
  for(i in 1:nSubj)
  { 
    dati=frame[which(frame$cid == u.cid[i]),]
    scorei=dati$score
    statusi=dati$status
    for(j in 1:nSubUnits[i])
    { 
      score.ij=scorei[j]		   
      status.ij=statusi[j]		   
      for (k in 1:nSubj)
      { 
        datk=frame[which(frame$cid == u.cid[k]),]
        scorek=datk$score
        statusk=datk$status

        for (ell in 1:nSubUnits[k])
        {
          score.kell=scorek[ell]		   
          status.kell=statusk[ell]		   
          if(!is.na(score.ij) && !is.na(score.kell))
          { 
	    numer<-numer+Ufunc(score.ij, score.kell)*
              (1-status.ij)*status.kell
            denom<-denom+(1-status.ij)*status.kell
          }
        }
      }
    }
  }

  stat<-numer/denom

  return(stat)
}




# calculate E(Tstat)
# theta -- Pr(Zij< Zkell)
# theta.c -- Pr(Zij < Ziell)
# delta -- nSubj x nSubUnits indicator matrix
#         delta[i,j]=1 means the j-th eye of the i-th subject is progressed
#         delta[i,j]=0 means the j-th eye of the i-th subject is not progressed
# frame - data frame with columns: cid, subuid, status, score
E.Tstat.func<-function(theta, theta.c, ci, di)
{

  C<-sum(ci)
  D<-sum(di)
  CD<-C*D
  A<-sum(ci*di)
  B<-CD-A

  ETstat<-(theta.c*A+theta*B)/CD

  return(ETstat)

}

# calculate Var(Tstat)=(Var(A)+Var(B)+2*Cov(A, B))/(C*D)^2
#
# theta -- Pr(Zij< Zkell)
# theta.c -- Pr(Zij < Ziell)
# delta -- nSubj x nSubUnits indicator matrix
#         delta[i,j]=1 means the j-th eye of the i-th subject is progressed
#         delta[i,j]=0 means the j-th eye of the i-th subject is not progressed
# rho -- corr(H(Zij), H(Ziell))
# frame - data frame with columns: cid, subuid, status, score
var.Tstat.func<-function(theta, theta.c, rho, ci, di)
{
  C<-sum(ci)
  D<-sum(di)
  CD<-C*D

  sumcidi<-sum(ci*di)
  sumci2di2<-sum(ci^2*di^2)
  sumci2di<-sum(ci^2*di)
  sumcidi2<-sum(ci*di^2)
  sumci2<-sum(ci^2)
  sumdi2<-sum(di^2)

  i.theta<-qnorm(theta) 
  i.theta.c<-qnorm(theta.c) 


  #########
  # calculate varA
  #########
  part1<-theta.c*(1-theta.c)*sumcidi

  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-1/2
  corMat[2,1]<-corMat[1,2]
  Phi2.half<-as.numeric(pmvnorm(upper=c(i.theta.c, i.theta.c), corr=corMat))

  part2<-(Phi2.half-theta.c^2)*(sumci2di+sumcidi2-2*sumcidi)

  varA<-part1+part2


  #########
  # calculate varB
  #########

  ############
  # deltaB1
  ############
  part1<-theta*(1-theta)*(CD-sumcidi)

  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-(1+rho)/2
  corMat[2,1]<-corMat[1,2]
  Phi2.1p.rho.d2<-as.numeric(pmvnorm(upper=c(i.theta, i.theta), corr=corMat))
  tmp<-D*sumci2+C*sumdi2-2*CD-sumci2di-sumcidi2+2*sumcidi
  part2<-(Phi2.1p.rho.d2-theta^2)*tmp

  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-rho
  corMat[2,1]<-corMat[1,2]
  Phi2.rho<-as.numeric(pmvnorm(upper=c(i.theta, i.theta), corr=corMat))
  tmp<-(sumdi2-D)*(sumci2-C)-sumci2di2+sumcidi2+sumci2di-sumcidi
  part3<-(Phi2.rho-theta^2)*tmp

  deltaB1<-part1+part2+part3

  ############
  # deltaB2
  ############

  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-1/2
  corMat[2,1]<-corMat[1,2]
  tmp<-C^2*D-D*sumci2-2*(C*sumcidi-sumci2di)
  Phi2.half<-as.numeric(pmvnorm(upper=c(i.theta, i.theta), corr=corMat))

  part1<-(Phi2.half-theta^2)*tmp

  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-rho/2
  corMat[2,1]<-corMat[1,2]
  tmp<-(sumdi2-D)*(C^2-sumci2)-2*C*sumcidi2+2*C*sumcidi+2*sumci2di2-2*sumci2di
  Phi2.rho.d2<-as.numeric(pmvnorm(upper=c(i.theta, i.theta), corr=corMat))

  part2<-(Phi2.rho.d2-theta^2)*tmp

  deltaB2<-part1+part2

  ############
  # deltaB3
  ############

  tmp<-C*(D^2-sumdi2)-2*(D*sumcidi-sumcidi2)
  part1<-(Phi2.half-theta^2)*tmp

  tmp<-(sumci2-C)*(D^2-sumdi2)-2*(D*(sumci2di-sumcidi)-sumci2di2+sumcidi2)
  part2<-(Phi2.rho.d2-theta^2)*tmp

  deltaB3<-part1+part2

  ############
  # deltaB4
  ############
  deltaB4<-0


  ############
  # deltaB5
  ############
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<- -rho
  corMat[2,1]<-corMat[1,2]
  tmp<-sumcidi^2-sumci2di2
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta, i.theta), corr=corMat))

  deltaB5<-(Phi2-theta^2)*tmp

  ############
  # deltaB6
  ############
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<- -rho/2
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta, i.theta), corr=corMat))

  tmp<-C*(D*sumcidi-sumcidi2)-D*sumci2di+2*sumci2di2-sumcidi^2
  deltaB6<-(Phi2-theta^2)*tmp

  ############
  # deltaB7
  ############
  tmp<-D*(C*sumcidi-sumci2di)-C*sumcidi2+2*sumci2di2-sumcidi^2
  deltaB7<-(Phi2-theta^2)*tmp

  ############
  # deltaB8
  ############
  deltaB8<-0

  varB<-deltaB1+deltaB2+deltaB3+deltaB4+deltaB5+deltaB6+deltaB7+deltaB8

  #########
  # calculate covAB
  #########
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-sqrt(1-rho)/2
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta, i.theta.c), corr=corMat))

  tmp<-C*sumcidi-sumci2di
  deltaAB1<-(Phi2-theta^2)*tmp

  tmp<-D*sumcidi-sumcidi2
  deltaAB2<-(Phi2-theta^2)*tmp

  deltaAB3<-0

  covAB<-deltaAB1+deltaAB2+deltaAB3

  #######
  # var(hat{theta}_c)
  #######
  numer<-varA+varB+2*covAB
  varTstat<-numer/CD^2

  return(varTstat)
}




# empirical cumulative distribution function
eCDF.func<-function(y)
{
  n<-length(y)
  mat<-cbind(1:n, y)

  mat2<-mat[order(mat[,2]),]
  mat3<-cbind(mat2, (1:n)/(n+1))

  mat4<-mat3[order(mat3[,1]),]
  
  eCDF<-as.numeric(mat4[,3])

  return(eCDF)
}

# get two-sided pvalue based on z-value
calPval<-function(Z)
{
  # pvalue
  if(Z<=0)
  { 
    pval<-2*pnorm(Z)
  } else {
    pval<-2*(1-pnorm(Z))
  }
  return(pval)
}

# get intra-class correlation
# y.1 - replicate 1
# y.2 - replicate 2
intraCorr.func<-function(y.1, y.2)
{
  y<-c(y.1, y.2)
  nSubj<-length(y.1)
  memSubj<-rep(1:nSubj, 2)

  mat.aov<-aov(y~as.factor(memSubj))
  res<-summary(mat.aov)

  res2<-res[[1]]
  MSb<-res2[1,3] # between subject variation
  MSe<-res2[2,3] # within subject variation
  # intra-class correlation
  rho<-(MSb-MSe)/(MSb+MSe)

  return(rho)
  
}

# generate simulated data set for one prediction rule
# mu -- E(H(Z_{ij})|the j-th subunit progressed)
# delta -- nSubj x nSubUnits matrix
#          nSubj -- number of subjects
#          nSubUnits -- number of subunits per subject
# rho -- correlation between H(Z_{ij}) and H(Z_{i ell})
# pVec -- 'nSubUnits' x 1 vector. 
#      'pVec[j]'=probability of the $j$-th subunit progressed 
# output:
#    dat -- 'nSubj' x 'nSubUnits' matrix
genDat1<-function(nSubj, mu, pVec=c(0.2, 0.2), rho=0.8)
{
  nSubUnits = 2
  delta = genDelta(nSubj, nSubUnits=nSubUnits, pVec=pVec)

  dat<-matrix(0, nrow=nSubj, ncol=nSubUnits)

  covMat<-matrix(rho, nrow=nSubUnits, ncol=nSubUnits)
  diag(covMat)<-1
  for(i in 1:nSubj)
  {
    muVec<-rep(0, nSubUnits)
    muVec[which(delta[i,]==1)]<-mu
    dat[i,]<-MASS::mvrnorm(1, mu=muVec, Sigma=covMat)
  }

  score=c(dat)
  status=c(delta)

  cid=rep(1:nSubj,rep(nSubUnits, nSubj))
  subuid=rep(1:nSubUnits, nSubj)

  frame=data.frame(cid=cid, subuid=subuid, status=status, score=score)
  rownames(frame)=paste("subj", 1:nrow(frame), sep="")

  invisible(frame)
}

# generate data set for two prediction rules
# mu1 -- E(H(Zij^{(1)})) # mean for prediction rule 1
# mu2 -- E(H(Zij^{(2)})) # mean for prediction rule 2
# delta -- nSubj x nSubUnits matrix
#          nSubj -- number of subjects
#          nSubUnits -- number of subunits per subject#
# rho -- cov(Hij^{(1)}, Hij^{(2)})
# rho11 -- cov(Hij^{(1)}, Hit^{(1)})
# rho22 -- cov(Hij^{(2)}, Hit^{(2)})
# rho12 -- cov(Hij^{(1)}, Hit^{(2)})
#
# output
#  dat - 'nSubj' x '2*nSubUnits' matrix
#
#
#genDat12<-function(nSubj, mu1 = 0.8, triangle=0.05, 
#  probVec=c(0.115, 0.142, 0.130, 1-0.115-0.142-0.130),
#  rhoVec=c(0.93, 0.59, 0.56, 0.52))

# generate ci and di 
gencidi<-function(nSubj, mu1 = 0.8, triangle=0.05, 
  probVec=c(0.115, 0.142, 0.130, 1-0.115-0.142-0.130),
  rhoVec=c(0.93, 0.59, 0.56, 0.52))

{

  nSubUnits<-2
  delta = genDelta2(nSubj=nSubj, probVec=probVec)

  status=c(delta)

  cid=rep(1:nSubj,rep(nSubUnits, nSubj))

  # obtain basic quantities
  # won't sort cid
  u.cid=unique(cid)

  nSubj<-length(u.cid)
  # tapply will sort the cid
  nSubUnits<-tapply(cid, cid, length)

  pos=match(u.cid, names(nSubUnits))
  nSubUnits = nSubUnits[pos]

  # number of eyes that have progressed for the i-th subject
  ci<-tapply(status, cid, sum, na.rm=TRUE)

  # number of eyes that have *NOT* progressed for the i-th subject
  di<-tapply(status, cid, function(x) {sum(1-x, na.rm=TRUE)})

  res=list(ci=ci, di=di)

  invisible(res)
}


genDelta<-function(nSubj, nSubUnits=2, pVec=c(0.2, 0.2))
{
  delta<-matrix(0, nrow=nSubj, ncol=nSubUnits)

  for(j in 1:nSubUnits)
  { 
    delta[,j]<-sample(c(1,0), size=nSubj, prob=c(pVec[j], 1-pVec[j]), replace=TRUE)
  }

  colnames(delta)=paste("delta", 1:ncol(delta), sep="")
  rownames(delta)=paste("subj", 1:nrow(delta), sep="")
  invisible(delta)
}


# variance of the difference between two test statistics derived 
# from two prediction rules for the same set of subjects

# var(theta_c^{(1)}-theta_c^{(2)})
var.Tstat12.func<-function(theta.1, theta.c1, theta.2, theta.c2, 
rho, rho11, rho22, rho12, ci, di)
{
  varT1<- var.Tstat.func(theta.1, theta.c1, rho11, ci, di)
  varT2<- var.Tstat.func(theta.2, theta.c2, rho22, ci, di)
  covT1T2<-cov.Tstat12.func(theta.1, theta.c1, theta.2, theta.c2, 
rho, rho11, rho22, rho12, ci, di)  

  varT12<-varT1+varT2-2*covT1T2

  return(varT12)
}

# covariance for the test statistics derived from two prediction rules
# for the same set of subjects
#
# Cov(theta_c^{(1)}, theta_c^{(2)})
#
# theta.1 -- Pr(Z_{ij}^{(1)} < Z_{k ell}^{(1)})
# theta.c1 -- Pr(Z_{ij}^{(1)} < Z_{i ell}^{(1)})
# theta.2 -- Pr(Z_{ij}^{(2)} < Z_{k ell}^{(2)})
# theta.c2 -- Pr(Z_{ij}^{(2)} < Z_{i ell}^{(2)})
#
# rho=cov(H_{ij}^{(1)}, H_{ij}^{(2)})
# rho11=cov(H_{ij}^{(1)}, H_{i ell}^{(1)})
# rho22=cov(H_{ij}^{(2)}, H_{i ell}^{(2)})
# rho12=cov(H_{ij}^{(1)}, H_{i ell}^{(2)})
#
# delta -- nSubj x nSubUnits indicator matrix
#         delta[i,j]=1 means the j-th eye of the i-th subject is progressed
#         delta[i,j]=0 means the j-th eye of the i-th subject is not progressed
# frame - data frame with columns: cid, subuid, status, score
cov.Tstat12.func<-function(theta.1, theta.c1, theta.2, theta.c2, 
rho, rho11, rho22, rho12, ci, di)
{

  i.theta.1<-qnorm(theta.1) 
  i.theta.2<-qnorm(theta.2) 
  i.theta.c1<-qnorm(theta.c1) 
  i.theta.c2<-qnorm(theta.c2) 

  C<-sum(ci)
  D<-sum(di)
  CD<-C*D
  sumcidi<-sum(ci*di)
  sumci2di<-sum(ci^2*di)
  sumcidi2<-sum(ci*di^2)
  sumci2<-sum(ci^2)
  sumdi2<-sum(di^2)
  sumci2di2<-sum(ci^2*di^2)
  
  #######
  # calculate cov(A^{(1)}, A^{(2)})
  #######

  corMat<-matrix(1, 2, 2)
  tt<-(rho-rho12)/(sqrt((1-rho11)*(1-rho22)))
  if(tt>=1)
  { 
    tt=0.99
  }
  corMat[1,2]<-tt
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.c1, i.theta.c2), corr=corMat))

  part1<-(Phi2-theta.c1*theta.c2)*sumcidi

  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-(rho-rho12)/(2*sqrt((1-rho11)*(1-rho22)))
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.c1, i.theta.c2), corr=corMat))

  part2<-(Phi2-theta.c1*theta.c2)*(sumci2di+sumcidi2-2*sumcidi)

  covA1A2<-part1+part2

  #######
  # calculate cov(A^{(1)}, B^{(2)})
  #######
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-(rho-rho12)/(2*sqrt((1-rho11)))
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.c1, i.theta.2), corr=corMat))
  
  t1<-(C+D)*sumcidi-sumci2di-sumcidi2
  covA1B2<-(Phi2-theta.c1*theta.2)*t1
  
  #######
  # calculate cov(B^{(1)}, A^{(2)})
  #######
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-(rho-rho12)/(2*sqrt((1-rho22)))
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.1, i.theta.c2), corr=corMat))
  
  covB1A2<-(Phi2-theta.1*theta.c2)*t1


  #######
  # calculate cov(B^{(1)}, B^{(2)})
  #######
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-rho
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.1, i.theta.2), corr=corMat))
  
  t1<-C*D-sumcidi
  part1<-(Phi2-theta.1*theta.2)*t1

  ####
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-(rho+rho12)/2
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.1, i.theta.2), corr=corMat))
  
  t1<-2*sumcidi+C*sumdi2+D*sumci2-sumcidi2-sumci2di-2*C*D
  part2<-(Phi2-theta.1*theta.2)*t1

  ####
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-rho12
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.1, i.theta.2), corr=corMat))
  
  t1<-sumci2*sumdi2-C*sumdi2-D*sumci2+CD-sumci2di2+sumci2di+sumcidi2-sumcidi
 
  part3<-(Phi2-theta.1*theta.2)*t1

  ####
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-rho/2
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.1, i.theta.2), corr=corMat))
  
  t1<-C*(D^2-sumdi2)-2*(D*sumcidi-sumcidi^2)
 
  part4<-(Phi2-theta.1*theta.2)*t1

  ####
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-rho12/2
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.1, i.theta.2), corr=corMat))
  
  t1<-(sumci2-C)*(D^2-sumdi2)-2*(D*(sumci2di-sumcidi)-sumci2di2+sumcidi2)
 
  part5<-(Phi2-theta.1*theta.2)*t1

  ####
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-rho/2
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.1, i.theta.2), corr=corMat))
  
  t1<-D*(C^2-sumci2)-2*(C*sumcidi-sumci2di)
 
  part6<-(Phi2-theta.1*theta.2)*t1

  ####
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<-rho12/2
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.1, i.theta.2), corr=corMat))
  
  t1<-(sumdi2-D)*(C^2-sumci2)-2*(C*(sumcidi2-sumcidi)-sumci2di2+sumci2di)
 
  part7<-(Phi2-theta.1*theta.2)*t1

  ####
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<- -rho12
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.1, i.theta.2), corr=corMat))
  
  t1<-sumcidi^2-sumci2di2
 
  part8<-(Phi2-theta.1*theta.2)*t1

  ####
  corMat<-matrix(1, 2, 2)
  corMat[1,2]<- -rho12/2
  corMat[2,1]<-corMat[1,2]
  Phi2<-as.numeric(pmvnorm(upper=c(i.theta.1, i.theta.2), corr=corMat))
  
  t1<-CD*sumcidi-C*sumcidi2-D*sumci2di-sumcidi^2+2*sumci2di2
 
  part9<-2*(Phi2-theta.1*theta.2)*t1

  covB1B2<-part1+part2+part3+part4+part5+part6+part7+part8+part9

  covT1T2<-(covA1A2+covA1B2+covB1A2+covB1B2)/CD^2

  return(covT1T2)
}



genDelta2<-function(nSubj, probVec=c(0.115, 0.142, 0.130, 1-0.115-0.142-0.130))
{
  nSubj<-as.integer(nSubj)

  if(sum(probVec)>1)
  { stop("sum(probVec) should be less than or equal to one!\n") }
  if(sum(probVec<=0))
  { stop("all elements of 'probVec' should be positive!\n") }
  if(nSubj < 1)
  { stop("nSubj should be positive integer!\n") }
  
  probMat<-matrix(probVec, nrow=1, byrow=T)
  tmp<-rMultinom(probMat, nSubj)

  deltaMat<-matrix(0, nrow=nSubj, ncol=2)
  colnames(deltaMat)<-c("delta1", "delta2")
  for(i in 1:nSubj)
  {
    tmp.1i=tmp[1,i]	   
    if(tmp.1i==1)
    { deltaMat[i,1] <- 1
      deltaMat[i,2] <- 1
    } else if (tmp.1i==2) {
      deltaMat[i,1] <- 1
      deltaMat[i,2] <- 0
    } else if (tmp.1i==3) {
      deltaMat[i,1] <- 0
      deltaMat[i,2] <- 1
    }
  }

  colnames(deltaMat)=c("delta1", "delta2")  
  rownames(deltaMat)=paste("subj", 1:nrow(deltaMat), sep="")
  invisible(deltaMat) 
}

# get mu2 based on mu1, rho11, rho22, triangle, and deltas
getDiffTriangles<-function(mu2, mu1, rho11, rho22, triangle, sumcidi, CD)
{

  theta1<-pnorm(mu1/sqrt(2))
  theta1.c<-pnorm(mu1/sqrt(2*(1-rho11)))

  numer.c1<-theta1.c*sumcidi+theta1*(CD-sumcidi)
  denom.c1<-CD
  eta.c1<- numer.c1/denom.c1

  theta2<-pnorm(mu2/sqrt(2))
  theta2.c<-pnorm(mu2/sqrt(2*(1-rho22)))

  numer.c2<-theta2.c*sumcidi+theta2*(CD-sumcidi)
  denom.c2<-CD
  eta.c2<- numer.c2/denom.c2

  trianglenew <- eta.c1-eta.c2

  diff<-abs(trianglenew-triangle)

  return(diff)
}

getMu2<-function(mu2Ini, mu1, rho11, rho22, triangle, 
    ci, di, low.mu2=-Inf, upp.mu2=Inf, ...)
{
  C<-sum(ci)
  D<-sum(di)
  CD<-C*D

  sumcidi<-sum(ci*di)

  res.opt<-optim(par=mu2Ini, fn=getDiffTriangles, method="L-BFGS-B",
    lower= low.mu2, upper=upp.mu2, mu1=mu1, rho11=rho11,
    rho22=rho22, triangle=triangle, sumcidi=sumcidi, CD=CD, ...)

  mu2<-res.opt$par

  return(mu2)
}

# calculate power based on pilot data
# nSubj - number of subjects to be recruited in the experiment
# frame - data frame with columns: cid, subuid, status, score1, score2
# alpha - type I error rate
powerCalData<-function(nSubj, triangle, frame, alpha=0.05)
{
  # risk score 1
  frame1=frame[, c("cid", "subuid", "status", "score1")]
  colnames(frame1)[4]="score"

  frame2=frame[, c("cid", "subuid", "status", "score2")]
  colnames(frame2)[4]="score"

  res1 = riskPredict(frame=frame1, alpha=alpha)
  res2 = riskPredict(frame=frame2, alpha=alpha)

  ## calculate rho, rho11, rho12, rho22
  datHk1 = res1$datHk
  datHk2 = res2$datHk

  rho11 = res1$rho
  rho22 = res2$rho

  ttx=c(datHk1[,1], datHk1[,2])
  tty=c(datHk2[,1], datHk2[,2])
  rho = stats::cor(x=ttx, y=tty, use="complete.obs")

  ttx=c(datHk1[,1], datHk1[,2])
  tty=c(datHk2[,2], datHk2[,1])
  rho12 = stats::cor(x=ttx, y=tty, use="complete.obs")

  rhoVec = c(rho, rho11, rho22, rho12)

  ######
  # calculate p11, p10, p01, and p00
  #####
  u.subj=unique(frame$cid)
  nSubj.data=length(u.subj)
  n11=0
  n10=0
  n01=0
  n00=0
  for(i in 1:nSubj.data)
  {
    dati=frame[which(frame$cid==u.subj[i]),,drop=FALSE]
    dati.s=dati[order(dati$subuid),,drop=FALSE]
    if(dati.s$status[1]==1 & dati.s$status[2]==1)
    {
      n11=n11+1
    } else if(dati.s$status[1]==1 & dati.s$status[2]==0) {
      n10=n10+1	    
    } else if (dati.s$status[1]==0 & dati.s$status[2]==1) {
      n01=n01+1	    
    } else {
      n00=n00+1
    }
  }
  p11=n11/nSubj.data
  p10=n10/nSubj.data
  p01=n01/nSubj.data
  p00=n00/nSubj.data

  probVec = c(p11, p10, p01, p00)

  # calculate mu1
  mu1=estMuFunc(status=frame1$status, score=frame1$status)

  ###
  theta1=pnorm(mu1/sqrt(2))

  ttres = gencidi(
         nSubj = nSubj,
         mu1 = mu1,
         triangle = triangle,
         probVec = probVec,
         rhoVec = rhoVec)
  ci=ttres$ci
  di=ttres$di

  ###
  theta1.c<-pnorm(mu1/sqrt(2*(1-rho11)))

  mu2<-getMu2(mu2Ini=mu1, mu1=mu1, rho11=rho11, rho22=rho22, 
    triangle=triangle, ci = ci, di = di, low.mu2=-Inf, upp.mu2=Inf)
  theta2<-pnorm(mu2/sqrt(2))
  theta2.c<-pnorm(mu2/sqrt(2*(1-rho22)))

  # mu2.null = mu2 under H0: triangle = 0
  mu2.null<-getMu2(mu2Ini=mu1, mu1=mu1, rho11=rho11, rho22=rho22, 
    triangle=0, ci=ci, di=di, low.mu2=-Inf, upp.mu2=Inf)

  theta2.null<-pnorm(mu2.null/sqrt(2))
  theta2.c.null<-pnorm(mu2.null/sqrt(2*(1-rho22)))

  # true parameter
  E.Tstat1<-E.Tstat.func(theta1, theta1.c, ci, di)
  E.Tstat2<-E.Tstat.func(theta2, theta2.c, ci, di)
  E.diffT<-E.Tstat1-E.Tstat2

  # calculate the theoretical variance of the test statistic
  sd.Tstat12<-sqrt(var.Tstat12.func(theta1, theta1.c, theta2, theta2.c, 
    rho, rho11, rho22, rho12, ci, di))
  
  # variance under H0
  sd0.Tstat12<-sqrt(var.Tstat12.func(theta1, theta1.c, theta2.null, theta2.c.null, rho, rho11, rho22, rho12, ci, di))
  
  za<-qnorm(1-alpha/2)
  t1<-(sd0.Tstat12*za-E.diffT)/sd.Tstat12
  t2<-(-sd0.Tstat12*za-E.diffT)/sd.Tstat12
  power<-1-pnorm(t1)+pnorm(t2)

  res=list(power=power, rho=rho, rho11=rho11, rho22=rho22, rho12=rho12, p11=p11, p10=p10, p01=p01, p00=p00,
	   mu1=mu1,  mu2=mu2)
  return(res)
}

powerCal<-function(nSubj, mu1, triangle, rho, rho11, rho22, rho12,
  p11, p10, p01, alpha = 0.05)
{
  rhoVec = c(rho, rho11, rho22, rho12)
  probVec = c(p11, p10, p01, 1-p11-p10-p01)
  theta1=pnorm(mu1/sqrt(2))

  ttres = gencidi(
         nSubj = nSubj,
         mu1 = mu1,
         triangle = triangle,
         probVec = probVec,
         rhoVec = rhoVec)
  ci=ttres$ci
  di=ttres$di

  ###
  theta1.c<-pnorm(mu1/sqrt(2*(1-rho11)))

  mu2<-getMu2(mu2Ini=mu1, mu1=mu1, rho11=rho11, rho22=rho22, 
    triangle=triangle, ci = ci, di = di, low.mu2=-Inf, upp.mu2=Inf)
  theta2<-pnorm(mu2/sqrt(2))
  theta2.c<-pnorm(mu2/sqrt(2*(1-rho22)))

  # mu2.null = mu2 under H0: triangle = 0
  mu2.null<-getMu2(mu2Ini=mu1, mu1=mu1, rho11=rho11, rho22=rho22, 
    triangle=0, ci=ci, di=di, low.mu2=-Inf, upp.mu2=Inf)

  theta2.null<-pnorm(mu2.null/sqrt(2))
  theta2.c.null<-pnorm(mu2.null/sqrt(2*(1-rho22)))

  # true parameter
  E.Tstat1<-E.Tstat.func(theta1, theta1.c, ci, di)
  E.Tstat2<-E.Tstat.func(theta2, theta2.c, ci, di)
  E.diffT<-E.Tstat1-E.Tstat2

  # calculate the theoretical variance of the test statistic
  sd.Tstat12<-sqrt(var.Tstat12.func(theta1, theta1.c, theta2, theta2.c, 
    rho, rho11, rho22, rho12, ci, di))
  
  # variance under H0
  sd0.Tstat12<-sqrt(var.Tstat12.func(theta1, theta1.c, theta2.null, theta2.c.null, rho, rho11, rho22, rho12, ci, di))
  
  za<-qnorm(1-alpha/2)
  t1<-(sd0.Tstat12*za-E.diffT)/sd.Tstat12
  t2<-(-sd0.Tstat12*za-E.diffT)/sd.Tstat12
  power<-1-pnorm(t1)+pnorm(t2)

  return(power)
}

print.class.riskPredict = function(x, ...)
{
  summary.class.riskPredict(x)
  invisible(x)
}

summary.class.riskPredict = function(object, ...) 
{
  x=object
  stopifnot(inherits(x, "class.riskPredict"))
  cat("\t\neta.hat=", round(x$stat, 3), 
      ", se.eta|H0=", round(x$se.stat, 3), 
      "\nz-value=(eta.hat-0.5)/se.eta=", round(x$z, 3),
      ", p-value=", sprintf("%3.2e", x$pval), 
      "\nCI.eta.low=", round(x$CIlow, 3),
      ", CI.eta.upp=", round(x$CIupp, 3))
  cat("\n")

  invisible(object)
}

print.class.riskPredictDiff = function(x, ...)
{
  summary.class.riskPredictDiff(x)
  invisible(x)
}

summary.class.riskPredictDiff = function(object, ...) 
{
  x=object
  stopifnot(inherits(x, "class.riskPredictDiff"))
  x1=x$res1
  cat("\t\neta1.hat=", round(x1$stat, 3),
      ", se.eta1|H0=", round(x1$se.stat, 3),
      "\nz-value1=(eta1.hat-0.5)/se.eta1=", round(x1$z, 3),
      ", p-value1=", sprintf("%3.2e", x1$pval),
      "\nCI.eta1.low=", round(x1$CIlow, 3),
      ", CI.eta1.upp=", round(x1$CIupp, 3))
  cat("\n")

  x2=x$res2
  cat("\t\neta2.hat=", round(x2$stat, 3),
      ", se.eta2|H0=", round(x2$se.stat, 3),
      "\nz-value2=(eta2.hat-0.5)/se.eta2=", round(x2$z, 3),
      ", p-value2=", sprintf("%3.2e", x2$pval),
      "\nCI.eta2.low=", round(x2$CIlow, 3),
      ", CI.eta2.upp=", round(x2$CIupp, 3))
  cat("\n")

  cat("\t\neta1.hat-eta2.hat=", round(x$diff, 3),
      ", se(eta1.hat-eta2.hat)|H0=", round(x$se.diff, 3),
      "\nz-value=", round(x$z, 3),
      ", p-value=", sprintf("%3.2e", x$pval),
      "\nCI.diff.low=", round(x$CIlow.diff, 3),
      ", CI.diff.upp=", round(x$CIupp.diff, 3))
  cat("\n")

  invisible(object)
}

