#Functions
require(SuperLearner)
library(boot)
require(medicaldata)
require(dplyr)
require(dagitty)
require(dosearch)
require(randomForest)
require(gam,quietly = T)

#reads in the used data
data_read<-function(){
  #read in data
  dat=indo_rct
  #limit to the relevant sites
  dat=dat[(dat$site=="1_UM"|dat$site=="2_IU"),]
  #make it binary
  dat$site=ifelse(dat$site=="2_IU",1,0)
  dat$rx=ifelse(dat$rx=="1_indomethacin",1,0)
  #remove bleed
  dat=dat[,1:32]
  #recode factors as dummies
  for (k in 6:29){
    dat[,k]=ifelse(dat[,k]=="0_no",0,1)

  }
  dat$gender=ifelse(dat$gender=="1_female",1,0) #gender==1 is female

  dat$status=ifelse(dat$status=="0_inpatient",0,1)
  #make categorial variable into dummies
  dat$type1=ifelse(dat$type=="1_type 1",1,0)
  dat$type2=ifelse(dat$type=="2_type 2",1,0)
  dat$type3=ifelse(dat$type=="3_type 3",1,0)
  #remove the ID and the type variable
  dat=dat[,-31]
  dat=dat[,-1]
  #####
  #omit rows with missings
  data=na.omit(data.frame(subset(dat,select=-c(site,rx,outcome)),
                          Y=dat$outcome,Z=dat$rx,S=dat$site)
  )

  Y=data$Y
  S=data$S
  Z=data$Z
  X=subset(data,select = -c(Z,S,Y))
  return(list(S,Y,Z,X,data))
}

#implements the omnibus test
run_omnibus<-function(W.=X,A.=S,Y.=Y){
  mmd_test = function(R, S, D.R, D.S, sig.meth = 'eig', num.reps = 1e4,
                      return.cutoff = FALSE) {
    n = length(R)

    D.R.mat1 = matrix(rep(D.R,n),nrow=n)
    D.R.mat2 = matrix(rep(D.R,each=n),nrow=n)

    D.S.mat1 = matrix(rep(D.S,n),nrow=n)
    D.S.mat2 = matrix(rep(D.S,each=n),nrow=n)

    R.mat1 = matrix(rep(R,n),nrow=n)
    R.mat2 = matrix(rep(R,each=n),nrow=n)

    S.mat1 = matrix(rep(S,n),nrow=n)
    S.mat2 = matrix(rep(S,each=n),nrow=n)

    EE = ((2*(R.mat1-R.mat2)*(D.R.mat2-D.R.mat1) + 1 - (4*(R.mat1-R.mat2)^2-2)*D.R.mat1*D.R.mat2)*exp(-(R.mat1-R.mat2)^2)
          - ((2*(S.mat1-R.mat2)*(D.R.mat2-D.S.mat1) + 1 - (4*(S.mat1-R.mat2)^2-2)*D.S.mat1*D.R.mat2)*exp(-(S.mat1-R.mat2)^2))
          - ((2*(R.mat1-S.mat2)*(D.S.mat2-D.R.mat1) + 1 - (4*(R.mat1-S.mat2)^2-2)*D.R.mat1*D.S.mat2)*exp(-(R.mat1-S.mat2)^2))
          + (2*(S.mat1-S.mat2)*(D.S.mat2-D.S.mat1) + 1 - (4*(S.mat1-S.mat2)^2-2)*D.S.mat1*D.S.mat2)*exp(-(S.mat1-S.mat2)^2))

    # EE = exp(-(R.mat1-R.mat2)^2) - 2*exp(-(S.mat1-R.mat2)^2) + exp(-(S.mat1-S.mat2)^2)

    if(sig.meth=='eig'){
      line.means = rowMeans(EE)
      EE.ctrd = EE - matrix(rep(line.means,n),nrow=n) - matrix(rep(line.means,each=n),nrow=n) + matrix(rep(mean(line.means),n^2),nrow=n)
      num.eigs = min(200,n)
      # tmp = eigen(EE.ctrd)$values/n
      tmp = eigs_sym(EE.ctrd,num.eigs,which='LA')$values/n
      num.pos.eigs = num.eigs # sum(tmp>0)
      draws=c(matrix(rnorm(num.reps*num.pos.eigs)^2-1,nrow=num.reps,ncol=num.pos.eigs)%*%cbind(tmp[1:num.pos.eigs]))
    }

    # U-statistic
    diag(EE) = 0
    est = (rbind(rep(1/(n-1),n)) %*% EE %*% cbind(rep(1/n,n)))[1,1]
    # V-statistic
    # est = (rbind(rep(1/n,n)) %*% EE %*% cbind(rep(1/n,n)))[1,1]

    if(sig.meth=='eig'){
      pval = mean(draws>n*est)
    } else if(sig.meth=='var'){
      pval = pchisq(est/(2*var(D.R)/n)+1,df=1,lower.tail=FALSE)
    }

    return(if(!return.cutoff){
      c(est,pval)
    }else{
      c(est,pval,
        if(sig.meth=='eig'){
          quantile(draws,0.95)
        }else{
          2*var(D.R)*(qchisq(0.95,df=1)-1)})})
  }
  est_psi_prob_binom =
    function(W=W., A=A., Y=Y., W.train = NULL, A.train = NULL, Y.train = NULL,
             sig.meth = 'var', est.g = TRUE,
             g0 = NULL,
             SL.library = c('SL.glm', 'SL.step',   'SL.glm.interaction')) {
      n=length(A)

      if(is.null(W.train) | is.null(A.train) | is.null(Y.train)){
        W.train = W
        A.train = A
        Y.train = Y
      }

      # Estimate outcome regressions
      Qbar.est = SuperLearner(Y=Y.train,X=data.frame(W=W.train,A=A.train),newX=data.frame(W=rbind(W,W),A=rep(c(0,1),each=n)),SL.library=SL.library, family='binomial')
      Qbar.est.0 = Qbar.est$SL.predict[,1][1:n]
      Qbar.est.1 = Qbar.est$SL.predict[,1][(n+1):(2*n)]

      if(est.g){
        gg = SuperLearner(Y=A,X=W,SL.library=SL.library,family='binomial')$SL.predict[,1]
      } else {
        gg = g0(W)
      }

      # Plug-in estimate of blip
      R = Qbar.est.1 - Qbar.est.0
      S = rep(0,n)

      D.R = A/gg * (Y-Qbar.est.1) - (1-A)/(1-gg) * (Y-Qbar.est.0)
      D.S = rep(0,n)

      return(mmd_test(R,S,D.R,D.S,sig.meth=sig.meth))
    }

  test=est_psi_prob_binom(W=X,A=S,Y)
  res=data.frame(qestimate=test[1],pvalue=test[2])
  return(res)
}

#simulates the example used in the avoid confounding case
simulate_confounding_example_data<-function(){
  #Simulation code for the Avoid confounding case
  #S=0
  set.seed(300)
  n=5000
  t_conf1=rbinom(n,1,0.6)
  t_conf2=runif(n,-5,5)
  t_conf3=rbinom(n,1,0.2)

  t_treat=rbinom(n,1,t_conf3*0.2+t_conf2*0.01+0.3+t_conf1*0.05)
  t_rand=rbinom(n,1,0.5)
  t_outcome_treat=rbinom(n,1,t_conf3*0.3+t_conf2*0.02+0.3+t_conf1*0.1*t_treat)
  t_outcome_rand=rbinom(n,1,t_conf3*0.3+t_conf2*0.02+0.3+t_conf1*0.1*t_rand)
  ######
  #S=1
  n=5000
  s_conf1=rbinom(n,1,0.3)
  s_conf2=runif(n,-5,5)
  s_conf3=rbinom(n,1,0.6)
  s_rand=rbinom(n,1,0.2)
  s_outcome=rbinom(n,1,s_conf3*0.3+s_conf2*0.02+0.3+s_conf1*0.1*s_rand)
  conf1=c(s_conf1,t_conf1)
  conf2=c(s_conf2,t_conf2)
  conf3=c(s_conf3,t_conf3)
  out=c(s_outcome,t_outcome_treat)
  treat=c(s_rand,t_treat)

  data=data.frame(Y=out,S=c(rep(1,n),rep(0,n)),Z=treat,conf1,conf2)
  S=data$S
  Y=data$Y
  Z=data$Z
  X=subset(data,select = -c(Z,S,Y))

  estcol=numeric(0)
  cicol1=numeric(0)
  cicol2=numeric(0)

  return(list(S=S,Z=Z,Y=Y,X=X,estcol=estcol,cicol1=cicol1,cicol2=cicol2,
              t_conf1=t_conf1,t_conf2=t_conf2,t_conf3=t_conf3,t_treat=t_treat,t_rand=t_rand,
              t_outcome_treat=t_outcome_treat,t_outcome_rand=t_outcome_rand,s_outcome=s_outcome,
              s_rand=s_rand))
}

#this function calls the estimator and bootstraps over the sample
run_bootstrap<-function(R=10,purpose,data,method,lib=c("SL.glm","SL.glm.interaction","SL.gam")){


  validate_boot<-function(data,indices,method,lib){
    dat <- data[indices,]
    S=data$S
    X=subset(data,select = -c(Z,S,Y))
    Y=data$Y
    Z=data$Z
    res<-Validate_RCT_results(S=S,Z=Z,Y=Y,X=X,method=method,lib.=lib)
    res1=as.numeric(res[2])

    return(res1)
  }
  overcome_boot<-function(data,indices,method,lib){
    dat <- data[indices,]
    S=data$S
    X=subset(data,select = -c(Z,S,Y))
    Y=data$Y
    Z=data$Z
    res<-Avoid_confoundig(Y=Y,Z=Z,S=S,X=X)
    re1=res$est
    return(re1)
  }
  gennew_boot<-function(data,indices,method=parent.frame()$method,lib=parent.frame()$lib){
    dat <- data[indices,]
    S=data$S
    X=subset(data,select = -c(Z,S,Y))
    Y=data$Y
    Z=data$Z
    res<-Generate_new_Evidence()
    res1=as.numeric(res[2])
    return(res1)
  }
  Estimate_validate=NULL
  CI_validate=NULL
  Estimate_gennew=NULL
  CI_gennew=NULL
  Estimate_overcome=NULL
  CI_overcome=NULL


  if(purpose=="gennew"){
    bootstrapIntervalls_gennew=boot(data,gennew_boot,R=R,method=method,lib=lib)
    CI_gennew=c(bootstrapIntervalls_gennew$t0-1.96*sd(bootstrapIntervalls_gennew$t),bootstrapIntervalls_gennew$t0+1.96*sd(bootstrapIntervalls_gennew$t))
    Estimate_gennew=bootstrapIntervalls_gennew$t0
    #res<-Generate_new_Evidence()
    #ate=res[1]
    ret=list(
      #source_ATE=ate,
      Estimate_gennew=Estimate_gennew,CI_gennew=CI_gennew)
  }
  if(purpose=="validate"){
    bootstrapIntervalls=boot(data,validate_boot,R=R,method=method,lib=lib)
    CI_validate=c(bootstrapIntervalls$t0-1.96*sd(bootstrapIntervalls$t),bootstrapIntervalls$t0+1.96*sd(bootstrapIntervalls$t))
    Estimate_validate=bootstrapIntervalls$t0
    #res<-Validate_RCT_results()
    #ate=res[1]
    ret=list(
      #source_ate=ate,
      Estimate_validate=Estimate_validate,CI_validate=CI_validate)
  }
  if(purpose=="overcome"){
    bootstrapIntervalls_overcome=boot(data,overcome_boot,R=R,method=method,lib=lib)
    CI_overcome=c(mean(bootstrapIntervalls_overcome$t0)-1.96*sd(bootstrapIntervalls_overcome$t),mean(bootstrapIntervalls_overcome$t)+1.96*sd(bootstrapIntervalls_overcome$t))
    Estimate_overcome=bootstrapIntervalls_overcome$t0
    #res=Avoid_confoundig(Y=Y,Z=Z,S=S,X=X)
    #cate=res$CATE
    #sate=res$SATE
    ret=list(
      #source_ate=sate,confounded_ate=cate,
      Estimate_overcome=Estimate_overcome,CI_overcome=CI_overcome)
  }

  return(ret)
}

#the following functions are for the different Purposes of transportability they wwill be called by the run_bootstrap function
#starting with the generate new evidence case
Generate_new_Evidence<-function(S=parent.frame()$S,
                                Z=parent.frame()$Z,
                                Y=parent.frame()$Y,
                                X=parent.frame()$X,
                                lib.=parent.frame()$lib,
                                method=parent.frame()$method){
  ATE=mean(Y[S==1&Z==1])-mean(Y[S==1&Z==0])
  #TATE=mean(Y[S==0&Z==1])-mean(Y[S==0&Z==0])
  est=numeric(0)
  if (method=="TMLE"){
    #Outcome and treatment set to zero giving no information on these cases to the learner
    #Will not work with NA's as the Super Learner function cannot handle these
    Z[S==0]=0
    Y[S==0]=0
    n=transport_compare(z=Z,y=Y,site=S,w=X,lib=lib.)
    est=n$est
  }
  return(list(ATE,est))
}

#the function for the validateRCTresults case
Validate_RCT_results<-function(S,Z,Y,X,method=method,lib.=lib){
  ATE=mean(Y[S==1&Z==1])-mean(Y[S==1&Z==0])

  TATE=mean(Y[S==0&Z==1])-mean(Y[S==0&Z==0])
  est=numeric(0)
  if (method=="TMLE"){
    n=transport_compare(z=Z,y=Y,site=S,w=X,lib=lib.)
    est=n$est
  }
  return(list(ATE,est))
}

#the function when working on the overcome confounding case
Avoid_confoundig<-function(S,Z,Y,X,
                           lib.=parent.frame()$lib,
                           method=parent.frame()$method){
  CATE=mean(Y[S==0&Z==1])-mean(Y[S==0&Z==0])
  SATE=mean(Y[S==1&Z==1])-mean(Y[S==1&Z==0])
  if (method=="TMLE"){

    #TATE=mean(t_outcome_rand[t_rand==1])-mean(t_outcome_rand[t_rand==0])
    p=transport_compare(z=parent.frame()$Z,y=parent.frame()$Y,site=parent.frame()$S,w=parent.frame()$X,lib=lib.)
    est=p$est
  }
  Out=list(est=est,SATE=SATE,CATE=CATE)
  return(Out)
}
