#Functions used in the paper for eitehr creating outputs or plots
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

#function helping the ggplots generated to look the way in the publication
my_theme = function() {
  theme(
    # add border
    panel.border = element_blank(),
    axis.line.x = element_line(linetype = "solid", colour = "black",  linewidth = 1.5),
    axis.line.y = element_line(linetype = "solid", colour = "black",  linewidth = 1.5),
    # color background
    panel.background = element_blank(),
    # modify grid
    panel.grid.major.x = element_line(linetype = "dotted", color = "black"),
    panel.grid.major.y = element_line(linetype = "dotted", color = "black"),
    # modify text, axis, titles and colour
    plot.title = element_text(colour = "black", size = 30, hjust = 0.5),
    axis.text.y = element_text(colour = "black", face = "italic", size = 15),
    axis.text.x = element_text(colour = "black", face = "italic", size = 20),
    axis.title = element_text(colour = "black", size = 15),
    axis.ticks = element_blank(),
    # legend
    legend.position = "right",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 24, margin = margin(t=0.3, b =0.3, unit = "cm")),
    #legend.key = element_rect(fill="white"),
    legend.key.size=unit(2, 'cm'),
    #strip text size
    strip.text=element_text(colour = "black", size = 16)
  )
}

#overall estimation function for the bootstrapped confidence intervals

call_Tbility<-function(data,fun="DR1_est",
                       libraries=c("SL.glm","SL.glm.interaction","SL.gam"),
                       B=10){
  
  cores=parallel::detectCores()
  #require(doParallel)
  require(dplyr)
  cl=parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  foo <- parallel::clusterEvalQ(cl, require(dplyr))  
  require(foreach)
  if (fun== "TMLE"){
    foo <-  parallel::clusterEvalQ(cl,require(SuperLearner))
    foo <-  parallel::clusterEvalQ(cl, transport_compare <-function(z=parent.frame()$Z,
                                                                    y=parent.frame()$Y,
                                                                    site=parent.frame()$S ,
                                                                    w=parent.frame()$X,
                                                                    lib=c("SL.glm","SL.glm.interaction","SL.gam")) {
      datw <- w
      n.dat <- nrow(w)
      
      # Calculate components of clever covariate
      cps_sl<-SuperLearner(Y=site,X=data.frame(w),
                           family=binomial(),
                           SL.library =lib)
      
      cps=cps_sl$SL.predict
      
      
      ### Z~1
      nzmodel="z~1"
      glm_cpz = glm(formula = nzmodel, data = data.frame(cbind( z = z, datw)),
                    family = "binomial")
      cpz <- predict(glm_cpz, type = "response")
      
      # Calculate clever covariate.
      g0w <- ((1 - cpz) * cps) / (1 - cps)
      g1w <- (cpz * cps) / (1 - cps)
      h0w <- ((1 - z) * I(site == 1)) / g0w
      h1w <- (z * I(site == 1)) / g1w
      
      y_Z_0<-SuperLearner(Y=y[site==1],X=data.frame(w=w,z=z)[site==1,],
                          family=binomial(),
                          SL.library = lib)
      
      fit_0<-predict.SuperLearner(y_Z_0, newdata = dplyr::mutate(data.frame(w=w,z=z),z=0),type="response")$pred
      
      fit_1<-predict.SuperLearner(y_Z_0, newdata = dplyr::mutate(data.frame(w=w,z=z),z=1),type="response")$pred
      fit_y<-predict.SuperLearner(y_Z_0, newdata =data.frame(w=w,z=z),type="response")$pred
      
      q <- cbind(fit_y,fit_0,fit_1)
      
      
      epsilon <- coef(glm(y ~ -1 + offset(q[, 1]) + h0w + h1w, family = "binomial",
                          subset = site == 1))
      # Update initial prediction.
      q1 <- q + c((epsilon[1] * h0w + epsilon[2] * h1w),
                  epsilon[1] / g0w, epsilon[2] / g1w)
      
      # Get efficient influence curve values for everyone
      tmleest <- mean(plogis(q1[, 3][site == 0])) - mean(plogis(q1[, 2][site == 0]))
      ps0 <- mean(I(site == 0))
      eic <- (((z * h1w / ps0) - ((1 - z) * h0w / ps0)) * (y - plogis(q[, 1]))) +
        (I(site == 0) / ps0 * plogis(q1[, 3])) - (I(site == 0) / ps0 * plogis(q1[, 2]))- (tmleest / ps0)
      
      results = list("est" = tmleest,
                     "var" = var(eic) / n.dat,
                     "eic" = eic[site == 0])
      return(results)
    })
  }
  else if(fun=="OM_est"){
    foo <-  parallel::clusterEvalQ(cl,OM_est<-function(data){
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      S1data_A1<-dplyr::select(S1data_A1, -S,-A)
      OM1mod<-glm(formula=Y~., data=S1data_A1)
      p1<- predict(OM1mod,newdata=data, type="response") 
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      S1data_A0<-dplyr::select(S1data_A0, -S,-A)
      OM0mod<-glm(formula=Y~., data=S1data_A0)
      p0<- predict(OM0mod,newdata=data, type="response") 
      data$p1<-p1
      data$p0<-p0
      S0sub<-dplyr::filter(data, S==0)
      OM_1<-mean(S0sub$p1)
      OM_0<-mean(S0sub$p0)
      OM<-mean(S0sub$p1)-mean(S0sub$p0)
      list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
      return(list)
    })
  }else if(fun=="DR3_est") {
    foo <-  parallel::clusterEvalQ(cl,DR3_est<-function(data){
      S0data<-dplyr::filter(data, S==0)
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      S1data_A1<-dplyr::select(S1data_A1, -S,-A,-p1,-p0)
      DR1mod<<-glm(formula=Y~.-w, data=S1data_A1, weights=w)
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      S1data_A0<-dplyr::select(S1data_A0, -S,-A,-p1,-p0)
      DR0mod<-glm(formula=Y~.-w, data=S1data_A0, weights=w)
      p1<- predict(DR1mod,newdata=S0data, type="response") 
      p0<- predict(DR0mod,newdata=S0data, type="response") 
      DR3_1<- mean(p1)
      DR3_0<- mean(p0)
      DR3<-mean(p1)-mean(p0) 
      list<-list(DR3_1=DR3_1,DR3_0=DR3_0, DR3=DR3,DR1mod=DR1mod, DR0mod=DR0mod)
      return(list)
    })
    foo <-  parallel::clusterEvalQ(cl,OM_est<-function(data){
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      S1data_A1<-dplyr::select(S1data_A1, -S,-A)
      OM1mod<-glm(formula=Y~., data=S1data_A1)
      p1<- predict(OM1mod,newdata=data, type="response") 
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      S1data_A0<-dplyr::select(S1data_A0, -S,-A)
      OM0mod<-glm(formula=Y~., data=S1data_A0)
      p0<- predict(OM0mod,newdata=data, type="response") 
      data$p1<-p1
      data$p0<-p0
      S0sub<-dplyr::filter(data, S==0)
      OM_1<-mean(S0sub$p1)
      OM_0<-mean(S0sub$p0)
      OM<-mean(S0sub$p1)-mean(S0sub$p0)
      list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
      return(list)
    })
    foo <-  parallel::clusterEvalQ(cl,generate_weights<-function(Smod,Amod, data){
      S1data<-dplyr::filter(data, S==1)
      w_reg<-glm(Smod, family="binomial", data=data)
      ps<- predict(w_reg,newdata=data, type="response") 
      w_reg2<-glm(Amod, family="binomial", data=S1data)
      pa<- predict(w_reg2,newdata=data, type="response") 
      w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
      data$w<-w
      list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
      return(list)
    })
  }else if(fun=="DR2_est"){
    foo <-  parallel::clusterEvalQ(cl,OM_est<-function(data){
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      S1data_A1<-dplyr::select(S1data_A1, -S,-A)
      OM1mod<-glm(formula=Y~., data=S1data_A1)
      p1<- predict(OM1mod,newdata=data, type="response") 
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      S1data_A0<-dplyr::select(S1data_A0, -S,-A)
      OM0mod<-glm(formula=Y~., data=S1data_A0)
      p0<- predict(OM0mod,newdata=data, type="response") 
      data$p1<-p1
      data$p0<-p0
      S0sub<-dplyr::filter(data, S==0)
      OM_1<-mean(S0sub$p1)
      OM_0<-mean(S0sub$p0)
      OM<-mean(S0sub$p1)-mean(S0sub$p0)
      list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
      return(list)
    })
    foo <-  parallel::clusterEvalQ(cl,generate_weights<-function(Smod,Amod, data){
      S1data<-dplyr::filter(data, S==1)
      w_reg<-glm(Smod, family="binomial", data=data)
      ps<- predict(w_reg,newdata=data, type="response") 
      w_reg2<-glm(Amod, family="binomial", data=S1data)
      pa<- predict(w_reg2,newdata=data, type="response") 
      w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
      data$w<-w
      list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
      return(list)
    })
    foo <-  parallel::clusterEvalQ(cl,DR2_est<-function(data){
      A<-data$A
      S<-data$S
      Y<-data$Y
      p1<-data$p1
      p0<-data$p0
      w<-data$w
      sum1_DR2<-sum(S*A*w*(Y-p1)) 
      sum0_DR2<-sum(S*(1-A)*w*(Y-p0)) 
      norm1<-(sum(S*A*w))^-1
      norm0<-(sum(S*(1-A)*w))^-1
      DR2_1<-norm1*sum1_DR2 + (sum(1-S)^-1)*sum((1-S)*p1)
      DR2_0<-norm0*sum0_DR2 + (sum(1-S)^-1)*sum((1-S)*p0)
      DR2<-DR2_1-DR2_0
      list<-list(DR2_1=DR2_1,DR2_0=DR2_0, DR2=DR2)
      return(list)
    })
  }else if(fun=="DR1_est"){
    foo <-  parallel::clusterEvalQ(cl,OM_est<-function(data){
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      S1data_A1<-dplyr::select(S1data_A1, -S,-A)
      OM1mod<-glm(formula=Y~., data=S1data_A1)
      p1<- predict(OM1mod,newdata=data, type="response") 
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      S1data_A0<-dplyr::select(S1data_A0, -S,-A)
      OM0mod<-glm(formula=Y~., data=S1data_A0)
      p0<- predict(OM0mod,newdata=data, type="response") 
      data$p1<-p1
      data$p0<-p0
      S0sub<-dplyr::filter(data, S==0)
      OM_1<-mean(S0sub$p1)
      OM_0<-mean(S0sub$p0)
      OM<-mean(S0sub$p1)-mean(S0sub$p0)
      list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
      return(list)
    })
    foo <-  parallel::clusterEvalQ(cl,generate_weights<-function(Smod,Amod, data){
      S1data<-dplyr::filter(data, S==1)
      w_reg<-glm(Smod, family="binomial", data=data)
      ps<- predict(w_reg,newdata=data, type="response") 
      w_reg2<-glm(Amod, family="binomial", data=S1data)
      pa<- predict(w_reg2,newdata=data, type="response") 
      w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
      data$w<-w
      list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
      return(list)
    })
    foo <-  parallel::clusterEvalQ(cl, DR1_est<-function(data){  
      A<-data$A
      S<-data$S
      Y<-data$Y
      p1<-data$p1
      p0<-data$p0
      w<-data$w
      DR1_1<-(sum((1-S))^-1)* sum(S*A*w*(Y-p1) + (1-S)*p1)
      DR1_0<-(sum((1-S))^-1)* sum(S*(1-A)*w*(Y-p0) + (1-S)*p0)
      DR1<-DR1_1-DR1_0
      list<-list(DR1_1=DR1_1,DR1_0=DR1_0, DR1=DR1)
      return(list)
    })
    
  }else if(fun=="IOW1_est"){
    foo <-  parallel::clusterEvalQ(cl,generate_weights<-function(Smod,Amod, data){
      S1data<-dplyr::filter(data, S==1)
      w_reg<-glm(Smod, family="binomial", data=data)
      ps<- predict(w_reg,newdata=data, type="response") 
      w_reg2<-glm(Amod, family="binomial", data=S1data)
      pa<- predict(w_reg2,newdata=data, type="response") 
      w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
      data$w<-w
      list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
      return(list)
    })
    foo <-  parallel::clusterEvalQ(cl,IOW1_est<-function(data){ 
      A<-data$A
      S<-data$S
      w<-data$w
      Y<-data$Y
      IOW1_1 <-(sum((1-S))^-1)* sum(A*S*w*Y)
      IOW1_0 <-(sum((1-S))^-1)* sum((1-A)*S*w*Y)
      IOW1 = IOW1_1 - IOW1_0
      return(list(IOW1_1=IOW1_1,IOW1_0=IOW1_0, IOW1=IOW1))
    })
  }else if(fun=="IOW2_est"){
    foo <-  parallel::clusterEvalQ(cl,IOW2_est<-function(data){
      S0data<-dplyr::filter(data, S==0)  
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      IOW1mod<-glm(formula=Y~1, data=S1data_A1, weights=w)
      p1<- predict(IOW1mod,newdata=S0data, type="response") 
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      IOW0mod<-glm(formula=Y~1, data=S1data_A0, weights=w)
      p0<- predict(IOW0mod,newdata=S0data, type="response") 
      IOW2_1<-mean(p1)
      IOW2_0<-mean(p0)
      IOW2<-mean(p1)-mean(p0)
      list<-list(IOW2_1=IOW2_1,IOW2_0=IOW2_0, IOW2=IOW2,IOW1mod=IOW1mod,IOW0mod=IOW0mod)
      return(list)
    })
    foo <-  parallel::clusterEvalQ(cl,generate_weights<-function(Smod,Amod, data){
      S1data<-dplyr::filter(data, S==1)
      w_reg<-glm(Smod, family="binomial", data=data)
      ps<- predict(w_reg,newdata=data, type="response") 
      w_reg2<-glm(Amod, family="binomial", data=S1data)
      pa<- predict(w_reg2,newdata=data, type="response") 
      w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
      data$w<-w
      list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
      return(list)
    })
    
  }
  
  
  
  
  out <- numeric(B)
  target=data[data$S==0,]
  source=data[data$S==1,]
  out<-foreach(i=1:B) %dopar% {
    indices <- sample(rownames(source),size=nrow(source),replace=T)
    #sample_data<-source[1:dim(source)[1] %in% indices,]
    sample_data<-source[indices,]
    data_set=rbind(sample_data,target)
    
    # Calculate the mean of the sample
    if (fun== "TMLE"){
      out[i] <-transport_compare(z=data_set$Z,
                                 y=data_set$Y,
                                 site=data_set$S,
                                 w=subset(data_set,select=-c(Z,S,Y)),
                                 lib=libraries)$est
      
    }
    else if(fun=="OM_est"){
      colnames(data_set)[3]="A"
      estimate<-try(eval(parse(text=paste(fun,"(data_set)[[3]]",sep=""))))
    }else if(fun=="DR3_est" |fun=="DR2_est"|fun=="DR1_est"){
      colnames(data_set)[3]="A"
      OM<-OM_est(data=data_set)
      DRdata=generate_weights(Smod=S~.-A-Y,Amod=A~.-S-Y,data=data_set)$dat
      DRdata$p1<-OM$p1
      DRdata$p0<-OM$p0
      estimate<-try(eval(parse(text=paste(fun,"(data=DRdata)[[3]]",sep=""))))
    }else{
      colnames(data_set)[3]="A"
      estdat=generate_weights(Smod=S~.-A-Y,Amod=A~.-S-Y,data=data_set)$dat
      estimate<-try(eval(parse(text=paste(fun,"(data=estdat)[[3]]",sep=""))))
    }
    
  }
  data=rbind(source,target) 
  require(SuperLearner)
  
  if (fun== "TMLE"){
    transport_compare <-function(z=parent.frame()$Z,     y=parent.frame()$Y,
                                 site=parent.frame()$S ,
                                 w=parent.frame()$X,
                                 lib=c("SL.glm","SL.glm.interaction","SL.gam")) {
      datw <- w
      n.dat <- nrow(w)
      
      # Calculate components of clever covariate
      cps_sl<-SuperLearner(Y=site,X=data.frame(w),
                           family=binomial(),
                           SL.library =lib)
      
      cps=cps_sl$SL.predict
      
      
      ### Z~1
      nzmodel="z~1"
      glm_cpz = glm(formula = nzmodel, data = data.frame(cbind( z = z, datw)),
                    family = "binomial")
      cpz <- predict(glm_cpz, type = "response")
      
      # Calculate clever covariate.
      g0w <- ((1 - cpz) * cps) / (1 - cps)
      g1w <- (cpz * cps) / (1 - cps)
      h0w <- ((1 - z) * I(site == 1)) / g0w
      h1w <- (z * I(site == 1)) / g1w
      
      y_Z_0<-SuperLearner(Y=y[site==1],X=data.frame(w=w,z=z)[site==1,],
                          family=binomial(),
                          SL.library = lib)
      
      fit_0<-predict.SuperLearner(y_Z_0, newdata = dplyr::mutate(data.frame(w=w,z=z),z=0),type="response")$pred
      
      fit_1<-predict.SuperLearner(y_Z_0, newdata = dplyr::mutate(data.frame(w=w,z=z),z=1),type="response")$pred
      fit_y<-predict.SuperLearner(y_Z_0, newdata =data.frame(w=w,z=z),type="response")$pred
      
      q <- cbind(fit_y,fit_0,fit_1)
      
      
      epsilon <- coef(glm(y ~ -1 + offset(q[, 1]) + h0w + h1w, family = "binomial",
                          subset = site == 1))
      # Update initial prediction.
      q1 <- q + c((epsilon[1] * h0w + epsilon[2] * h1w),
                  epsilon[1] / g0w, epsilon[2] / g1w)
      
      # Get efficient influence curve values for everyone
      tmleest <- mean(plogis(q1[, 3][site == 0])) - mean(plogis(q1[, 2][site == 0]))
      ps0 <- mean(I(site == 0))
      eic <- (((z * h1w / ps0) - ((1 - z) * h0w / ps0)) * (y - plogis(q[, 1]))) +
        (I(site == 0) / ps0 * plogis(q1[, 3])) - (I(site == 0) / ps0 * plogis(q1[, 2]))- (tmleest / ps0)
      
      results = list("est" = tmleest,
                     "var" = var(eic) / n.dat,
                     "eic" = eic[site == 0])
      return(results)
    }
  }
  else if(fun=="OM_est"){
    OM_est<-function(data){
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      S1data_A1<-dplyr::select(S1data_A1, -S,-A)
      OM1mod<-glm(formula=Y~., data=S1data_A1)
      p1<- predict(OM1mod,newdata=data, type="response") 
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      S1data_A0<-dplyr::select(S1data_A0, -S,-A)
      OM0mod<-glm(formula=Y~., data=S1data_A0)
      p0<- predict(OM0mod,newdata=data, type="response") 
      data$p1<-p1
      data$p0<-p0
      S0sub<-dplyr::filter(data, S==0)
      OM_1<-mean(S0sub$p1)
      OM_0<-mean(S0sub$p0)
      OM<-mean(S0sub$p1)-mean(S0sub$p0)
      list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
      return(list)
    }
  }else if(fun=="DR3_est") {
    DR3_est<-function(data){
      S0data<-dplyr::filter(data, S==0)
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      S1data_A1<-dplyr::select(S1data_A1, -S,-A,-p1,-p0)
      DR1mod<<-glm(formula=Y~.-w, data=S1data_A1, weights=w)
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      S1data_A0<-dplyr::select(S1data_A0, -S,-A,-p1,-p0)
      DR0mod<-glm(formula=Y~.-w, data=S1data_A0, weights=w)
      p1<- predict(DR1mod,newdata=S0data, type="response") 
      p0<- predict(DR0mod,newdata=S0data, type="response") 
      DR3_1<- mean(p1)
      DR3_0<- mean(p0)
      DR3<-mean(p1)-mean(p0) 
      list<-list(DR3_1=DR3_1,DR3_0=DR3_0, DR3=DR3,DR1mod=DR1mod, DR0mod=DR0mod)
      return(list)
    }
    OM_est<-function(data){
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      S1data_A1<-dplyr::select(S1data_A1, -S,-A)
      OM1mod<-glm(formula=Y~., data=S1data_A1)
      p1<- predict(OM1mod,newdata=data, type="response") 
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      S1data_A0<-dplyr::select(S1data_A0, -S,-A)
      OM0mod<-glm(formula=Y~., data=S1data_A0)
      p0<- predict(OM0mod,newdata=data, type="response") 
      data$p1<-p1
      data$p0<-p0
      S0sub<-dplyr::filter(data, S==0)
      OM_1<-mean(S0sub$p1)
      OM_0<-mean(S0sub$p0)
      OM<-mean(S0sub$p1)-mean(S0sub$p0)
      list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
      return(list)
    }
    generate_weights<-function(Smod,Amod, data){
      S1data<-dplyr::filter(data, S==1)
      w_reg<-glm(Smod, family="binomial", data=data)
      ps<- predict(w_reg,newdata=data, type="response") 
      w_reg2<-glm(Amod, family="binomial", data=S1data)
      pa<- predict(w_reg2,newdata=data, type="response") 
      w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
      data$w<-w
      list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
      return(list)
    }
  }else if(fun=="DR2_est"){
    OM_est<-function(data){
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      S1data_A1<-dplyr::select(S1data_A1, -S,-A)
      OM1mod<-glm(formula=Y~., data=S1data_A1)
      p1<- predict(OM1mod,newdata=data, type="response") 
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      S1data_A0<-dplyr::select(S1data_A0, -S,-A)
      OM0mod<-glm(formula=Y~., data=S1data_A0)
      p0<- predict(OM0mod,newdata=data, type="response") 
      data$p1<-p1
      data$p0<-p0
      S0sub<-dplyr::filter(data, S==0)
      OM_1<-mean(S0sub$p1)
      OM_0<-mean(S0sub$p0)
      OM<-mean(S0sub$p1)-mean(S0sub$p0)
      list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
      return(list)
    }
    generate_weights<-function(Smod,Amod, data){
      S1data<-dplyr::filter(data, S==1)
      w_reg<-glm(Smod, family="binomial", data=data)
      ps<- predict(w_reg,newdata=data, type="response") 
      w_reg2<-glm(Amod, family="binomial", data=S1data)
      pa<- predict(w_reg2,newdata=data, type="response") 
      w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
      data$w<-w
      list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
      return(list)
    }
    DR2_est<-function(data){
      A<-data$A
      S<-data$S
      Y<-data$Y
      p1<-data$p1
      p0<-data$p0
      w<-data$w
      sum1_DR2<-sum(S*A*w*(Y-p1)) 
      sum0_DR2<-sum(S*(1-A)*w*(Y-p0)) 
      norm1<-(sum(S*A*w))^-1
      norm0<-(sum(S*(1-A)*w))^-1
      DR2_1<-norm1*sum1_DR2 + (sum(1-S)^-1)*sum((1-S)*p1)
      DR2_0<-norm0*sum0_DR2 + (sum(1-S)^-1)*sum((1-S)*p0)
      DR2<-DR2_1-DR2_0
      list<-list(DR2_1=DR2_1,DR2_0=DR2_0, DR2=DR2)
      return(list)
    }
  }else if(fun=="DR1_est"){
    OM_est<-function(data){
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      S1data_A1<-dplyr::select(S1data_A1, -S,-A)
      OM1mod<-glm(formula=Y~., data=S1data_A1)
      p1<- predict(OM1mod,newdata=data, type="response") 
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      S1data_A0<-dplyr::select(S1data_A0, -S,-A)
      OM0mod<-glm(formula=Y~., data=S1data_A0)
      p0<- predict(OM0mod,newdata=data, type="response") 
      data$p1<-p1
      data$p0<-p0
      S0sub<-dplyr::filter(data, S==0)
      OM_1<-mean(S0sub$p1)
      OM_0<-mean(S0sub$p0)
      OM<-mean(S0sub$p1)-mean(S0sub$p0)
      list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
      return(list)
    }
    generate_weights<-function(Smod,Amod, data){
      S1data<-dplyr::filter(data, S==1)
      w_reg<-glm(Smod, family="binomial", data=data)
      ps<- predict(w_reg,newdata=data, type="response") 
      w_reg2<-glm(Amod, family="binomial", data=S1data)
      pa<- predict(w_reg2,newdata=data, type="response") 
      w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
      data$w<-w
      list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
      return(list)
    }
    DR1_est<-function(data){  
      A<-data$A
      S<-data$S
      Y<-data$Y
      p1<-data$p1
      p0<-data$p0
      w<-data$w
      DR1_1<-(sum((1-S))^-1)* sum(S*A*w*(Y-p1) + (1-S)*p1)
      DR1_0<-(sum((1-S))^-1)* sum(S*(1-A)*w*(Y-p0) + (1-S)*p0)
      DR1<-DR1_1-DR1_0
      list<-list(DR1_1=DR1_1,DR1_0=DR1_0, DR1=DR1)
      return(list)
    }
  }else if(fun=="IOW1_est"){
    generate_weights<-function(Smod,Amod, data){
      S1data<-dplyr::filter(data, S==1)
      w_reg<-glm(Smod, family="binomial", data=data)
      ps<- predict(w_reg,newdata=data, type="response") 
      w_reg2<-glm(Amod, family="binomial", data=S1data)
      pa<- predict(w_reg2,newdata=data, type="response") 
      w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
      data$w<-w
      list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
      return(list)
    }
    IOW1_est<-function(data){ 
      A<-data$A
      S<-data$S
      w<-data$w
      Y<-data$Y
      IOW1_1 <-(sum((1-S))^-1)* sum(A*S*w*Y)
      IOW1_0 <-(sum((1-S))^-1)* sum((1-A)*S*w*Y)
      IOW1 = IOW1_1 - IOW1_0
      return(list(IOW1_1=IOW1_1,IOW1_0=IOW1_0, IOW1=IOW1))
    }
  }else if(fun=="IOW2_est"){
    IOW2_est<-function(data){
      S0data<-dplyr::filter(data, S==0)  
      S1data_A1<-dplyr::filter(data, S==1 & A==1)
      IOW1mod<-glm(formula=Y~1, data=S1data_A1, weights=w)
      p1<- predict(IOW1mod,newdata=S0data, type="response") 
      S1data_A0<-dplyr::filter(data, S==1 & A==0)
      IOW0mod<-glm(formula=Y~1, data=S1data_A0, weights=w)
      p0<- predict(IOW0mod,newdata=S0data, type="response") 
      IOW2_1<-mean(p1)
      IOW2_0<-mean(p0)
      IOW2<-mean(p1)-mean(p0)
      list<-list(IOW2_1=IOW2_1,IOW2_0=IOW2_0, IOW2=IOW2,IOW1mod=IOW1mod,IOW0mod=IOW0mod)
      return(list)
    }
    generate_weights<-function(Smod,Amod, data){
      S1data<-dplyr::filter(data, S==1)
      w_reg<-glm(Smod, family="binomial", data=data)
      ps<- predict(w_reg,newdata=data, type="response") 
      w_reg2<-glm(Amod, family="binomial", data=S1data)
      pa<- predict(w_reg2,newdata=data, type="response") 
      w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
      data$w<-w
      list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
      return(list)
    } 
  }
  if (fun== "TMLE"){
    estimate<-transport_compare(z=data$Z,
                                y=data$Y,
                                site=data$S,
                                w=subset(data,select=-c(Z,S,Y)),
                                lib=libraries)$est
  }else if(fun=="OM_est"){
    colnames(data)[3]="A"
    estimate<-try(eval(parse(text=paste(fun,"(data)[[3]]",sep=""))))
  }else if(fun=="DR3_est" |fun=="DR2_est"|fun=="DR1_est"){
    colnames(data)[3]="A"
    OM<-OM_est(data=data)
    DRdata=generate_weights(Smod=S~.-A-Y,Amod=A~.-S-Y,data=data)$dat
    DRdata$p1<-OM$p1
    DRdata$p0<-OM$p0
    estimate<-try(eval(parse(text=paste(fun,"(data=DRdata)[[3]]",sep=""))))
  }else{
    colnames(data)[3]="A"
    estdat=generate_weights(Smod=S~.-A-Y,Amod=A~.-S-Y,data=data)$dat
    estimate<-try(eval(parse(text=paste(fun,"(data=estdat)[[3]]",sep=""))))
  }
  se_boot <- sd(as.numeric(out),na.rm=T)
  return(list(est_boot = out,estimate=estimate, se_boot = se_boot))
}


#funtion to implement the TMLE for transportability following Rudolph et al. 2017

transport_compare <-function(z=parent.frame()$Z,
                             y=parent.frame()$Y,
                             site=parent.frame()$S ,
                             w=parent.frame()$X,
                             lib=c("SL.glm","SL.glm.interaction","SL.gam")) {
  datw <- w
  n.dat <- nrow(w)
  
  # Calculate components of clever covariate
  cps_sl<-SuperLearner(Y=site,X=data.frame(w),
                       family=binomial(),
                       SL.library =lib)
  
  cps=cps_sl$SL.predict
  
  
  ### Z~1
  nzmodel="z~1"
  glm_cpz = glm(formula = nzmodel, data = data.frame(cbind( z = z, datw)),
                family = "binomial")
  cpz <- predict(glm_cpz, type = "response")
  
  # Calculate clever covariate.
  g0w <- ((1 - cpz) * cps) / (1 - cps)
  g1w <- (cpz * cps) / (1 - cps)
  h0w <- ((1 - z) * I(site == 1)) / g0w
  h1w <- (z * I(site == 1)) / g1w
  
  y_Z_0<-SuperLearner(Y=y[site==1],X=data.frame(w=w,z=z)[site==1,],
                      family=binomial(),
                      SL.library = lib)
  
  fit_0<-predict.SuperLearner(y_Z_0, newdata = dplyr::mutate(data.frame(w=w,z=z),z=0),type="response")$pred
  
  fit_1<-predict.SuperLearner(y_Z_0, newdata = dplyr::mutate(data.frame(w=w,z=z),z=1),type="response")$pred
  fit_y<-predict.SuperLearner(y_Z_0, newdata =data.frame(w=w,z=z),type="response")$pred
  
  q <- cbind(fit_y,fit_0,fit_1)
  
  
  epsilon <- coef(glm(y ~ -1 + offset(q[, 1]) + h0w + h1w, family = "binomial",
                      subset = site == 1))
  # Update initial prediction.
  q1 <- q + c((epsilon[1] * h0w + epsilon[2] * h1w),
              epsilon[1] / g0w, epsilon[2] / g1w)
  
  # Get efficient influence curve values for everyone
  tmleest <- mean(plogis(q1[, 3][site == 0])) - mean(plogis(q1[, 2][site == 0]))
  ps0 <- mean(I(site == 0))
  eic <- (((z * h1w / ps0) - ((1 - z) * h0w / ps0)) * (y - plogis(q[, 1]))) +
    (I(site == 0) / ps0 * plogis(q1[, 3])) - (I(site == 0) / ps0 * plogis(q1[, 2]))- (tmleest / ps0)
  
  results = list("est" = tmleest,
                 "var" = var(eic) / n.dat,
                 "eic" = eic[site == 0])
  return(results)
}

#Estimation equation functions to run the transportability using sandwich variance  
{
  OM_EE <- function(data){
    A<-data$A
    S<- data$S
    Y <- data$Y
    X <- as.matrix(cbind(1, select(data,-S,-Y,-A))) 
    matA <- cbind(1, data$A) 
    A[is.na(A)] <- 0 
    Y[is.na(Y)] <- 0
    function(theta){ 
      #outcome model 
      beta<-theta[1:(length(coef(OM$OM1mod)))]
      alpha<-theta[(length(coef(OM$OM1mod))+1):(2*length(coef(OM$OM1mod)))]
      mu1<-theta[(2*length(coef(OM$OM1mod)))+1]
      mu0<-theta[(2*length(coef(OM$OM1mod)))+2]
      muate<-theta[(2*length(coef(OM$OM1mod)))+3]
      m_A1 <-X %*% beta
      m_A0<-X %*% alpha
      ols_A1 <-crossprod(X, (S*A)*(Y - m_A1))
      ols_A0 <-crossprod(X, (S*(1-A))*(Y - m_A0))
      #estimates
      mean1<-(1-S)*(m_A1-mu1) 
      mean0 <- (1-S)*(m_A0-mu0) 
      ate<-(1-S)*(m_A1-m_A0-muate) 
      c(ols_A1,ols_A0,mean1, mean0,ate)
    }
  }
  
  
  IOW1_EE <- function(data){ 
    A<-data$A
    S<- data$S
    Y <- data$Y  
    X <- as.matrix(cbind(1, select(data,-S,-Y,-A)))
    matA <- cbind(1, data$A) 
    A[is.na(A)] <- 0 
    Y[is.na(Y)] <- 0 
    function(theta){
      mu1<-theta[(2*length(coef(OM$OM1mod)))+1]
      mu0<-theta[(2*length(coef(OM$OM1mod)))+2]
      mu_ate<-theta[(2*length(coef(OM$OM1mod)))+3]
      #participation model
      lp  <- X %*% theta[1:(length(coef(OM$OM1mod)))]
      ps <- plogis(lp)
      score_eqns<-crossprod(X, S-ps)
      #treatment model
      lp2  <- X %*% theta[(length(coef(OM$OM1mod))+1):(2*length(coef(OM$OM1mod)))]
      pa<- plogis(lp2)
      score_eqns2<-crossprod(X,S*(A - pa) )
      w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
      #estimates
      summand1<-(w*A*S*Y)  
      summand0<-(w*(1-A)*S*Y)  
      mean1<-(summand1)- (1-S)*mu1
      mean0<-(summand0)- (1-S)*mu0
      ate<-(summand1)-(summand0)- (1-S)*mu_ate
      c(score_eqns,score_eqns2,mean1,mean0, ate)
    }
  }
  
  
  IOW2_EE <- function(data){
    A<-data$A
    S<- data$S
    Y <- data$Y
    X <- as.matrix(cbind(1, select(data,-S,-Y,-A)))
    matA <- cbind(1, data$A) 
    A[is.na(A)] <- 0 
    Y[is.na(Y)] <- 0 
    function(theta){
      #participation model
      lp  <- X %*% theta[1:(length(coef(OM$OM1mod)))]
      ps <- plogis(lp)
      score_eqns<-crossprod(X, S-ps)
      #treatment model
      lp2  <- X %*% theta[(length(coef(OM$OM1mod))+1):(2*length(coef(OM$OM1mod)))]
      pa<- plogis(lp2)
      score_eqns2<-crossprod(X,S*(A - pa) )
      w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
      #outcome model
      m_A1<-1 %*% theta[(2*length(coef(OM$OM1mod)))+1]
      m_A0<-1 %*% theta[(2*length(coef(OM$OM1mod)))+2]
      linear_eqns1<-crossprod(1, (S*A*w)*(Y -  m_A1) )
      linear_eqns0<-crossprod(1, (S*(1-A)*w)*(Y - 1 %*% theta[10]) )
      mu1<-theta[(2*length(coef(OM$OM1mod)))+3]
      mu0<-theta[(2*length(coef(OM$OM1mod)))+4]
      muate<-theta[(2*length(coef(OM$OM1mod)))+5]
      #estimates
      mean1 <- (1-S)*( m_A1 -mu1) 
      mean0 <- (1-S)*( m_A0 -mu0)
      ate <- (1-S)*(m_A1- m_A0 - muate)
      c(score_eqns,score_eqns2, linear_eqns1,linear_eqns0,  mean1,  mean0,ate)
    }
  }
  
  
  DR1_EE <- function(data){
    A<-data$A
    S<- data$S
    Y <- data$Y
    X <- as.matrix(cbind(1, select(data,-S,-Y,-A)))  
    matA <- cbind(1, data$A)  
    A[is.na(A)] <- 0 
    Y[is.na(Y)] <- 0 
    function(theta){
      #participation model
      lp  <- X %*% theta[1:(length(coef(OM$OM1mod)))]
      ps <- plogis(lp)
      score_eqns<-crossprod(X, S-ps)
      #treatment model
      lp2  <- X %*% theta[(length(coef(OM$OM1mod))+1):(2*length(coef(OM$OM1mod)))]
      pa<- plogis(lp2)
      score_eqns2<-crossprod(X,S*(A - pa) )
      w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
      #outcome model
      beta<-theta[(2*length(coef(OM$OM1mod))+1):(3*length(coef(OM$OM1mod)))]
      alpha<-theta[(3*length(coef(OM$OM1mod))+1):(4*length(coef(OM$OM1mod)))]
      mu1<-theta[(4*length(coef(OM$OM1mod)))+1]
      mu0<-theta[(4*length(coef(OM$OM1mod)))+2]
      mu<-theta[(4*length(coef(OM$OM1mod)))+3]
      m_A1 <-X %*% beta
      m_A0<-X %*% alpha
      ols_A1 <-crossprod(X, (S*A)*(Y - m_A1))
      ols_A0 <-crossprod(X, (S*(1-A))*(Y - m_A0))
      ey1<-w*S*A*(Y-m_A1) + (1-S)*m_A1
      ey0<-w*S*(1-A)*(Y-m_A0) + (1-S)*m_A0
      #estimates
      mean1<-ey1-(1-S)*mu1
      mean0<-ey0-(1-S)*mu0
      ate<-ey1-ey0- (1-S)*mu
      c(score_eqns,score_eqns2,ols_A1,ols_A0,mean1,mean0, ate)   
    }
  }
  
  
  DR2_EE <- function(data){
    A<-data$A
    S<- data$S
    Y <- data$Y
    X <- as.matrix(cbind(1,select(data,-S,-Y,-A)))  
    matA <- cbind(1, data$A) 
    A[is.na(A)] <- 0 
    Y[is.na(Y)] <- 0
    function(theta){
      #participation model
      lp  <- X %*% theta[1:(length(coef(OM$OM1mod)))]
      ps <- plogis(lp)
      score_eqns<-crossprod(X, S-ps)
      #treatment model
      lp2  <- X %*% theta[(length(coef(OM$OM1mod))+1):(2*length(coef(OM$OM1mod)))]
      pa<- plogis(lp2)
      score_eqns2<-crossprod(X,S*(A - pa) ) 
      w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
      #outcome model 
      beta<-theta[(2*length(coef(OM$OM1mod))+1):(3*length(coef(OM$OM1mod)))]
      alpha<-theta[(3*length(coef(OM$OM1mod))+1):(4*length(coef(OM$OM1mod)))]
      mu1<-theta[(4*length(coef(OM$OM1mod)))+4]
      mu0<-theta[(4*length(coef(OM$OM1mod)))+5]
      mu<-theta[(4*length(coef(OM$OM1mod)))+6]
      m_A1 <-X %*% beta
      m_A0<-X %*% alpha
      ols_A1 <-crossprod(X, (S*A)*(Y - m_A1))
      ols_A0 <-crossprod(X, (S*(1-A))*(Y - m_A0))
      mu_S<-theta[(4*length(coef(OM$OM1mod)))+1]
      propS1<-S-mu_S
      one_over<-(1/(1-mu_S))
      mu_norm1<-theta[(4*length(coef(OM$OM1mod)))+2]
      norm1eq<-(A*S*w)-mu_norm1
      norm1<-1/mu_norm1
      mu_norm0<-theta[(4*length(coef(OM$OM1mod)))+3]
      norm0eq<-((1-A)*S*w)-mu_norm0
      norm0<-1/mu_norm0
      ey1<-norm1*((w*S*A*(Y-m_A1))) + one_over*((1-S)*m_A1)
      ey0<-norm0*((w*S*(1-A)*(Y-m_A0))) + one_over*((1-S)*m_A0)
      #estimates
      mean1<-ey1-mu1
      mean0<-ey0-mu0
      ate<-ey1-ey0-mu
      c(score_eqns,score_eqns2,ols_A1,ols_A0,propS1, norm1eq, norm0eq,mean1,mean0,ate)
    }
  }
  
  
  DR3_EE <- function(data){
    A<-data$A
    S<- data$S
    Y <- data$Y
    X <- as.matrix(cbind(1, select(data,-S,-Y,-A)))  
    matA <- cbind(1, data$A) 
    A[is.na(A)] <- 0 
    Y[is.na(Y)] <- 0 
    function(theta){
      #participation model
      lp  <- X %*% theta[1:(length(coef(OM$OM1mod)))]
      ps <- plogis(lp)
      score_eqns<-crossprod(X, S-ps)
      #treatment model
      lp2  <- X %*% theta[(length(coef(OM$OM1mod))+1):(2*length(coef(OM$OM1mod)))]
      pa<- plogis(lp2)
      score_eqns2<-crossprod(X,S*(A - pa) )
      w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa))
      #outcome model
      beta<-theta[(2*length(coef(OM$OM1mod))+1):(3*length(coef(OM$OM1mod)))]
      alpha<-theta[(3*length(coef(OM$OM1mod))+1):(4*length(coef(OM$OM1mod)))]
      mu1<-theta[(4*length(coef(OM$OM1mod)))+1]
      mu0<-theta[(4*length(coef(OM$OM1mod)))+2]
      mu<-theta[(4*length(coef(OM$OM1mod)))+3]
      m_A1 <-X %*% beta
      m_A0<-X %*% alpha
      ols_A1 <-crossprod(X, (S*A*w)*(Y - m_A1))
      ols_A0 <-crossprod(X, (S*(1-A)*w)*(Y - m_A0))
      #estimates   
      mean1 <- (1-S)*(m_A1-mu1) 
      mean0 <- (1-S)*(m_A0-mu0) 
      mean <- (1-S)*(m_A1-m_A0-mu) 
      c(score_eqns,score_eqns2,ols_A1,ols_A0,mean1,mean0, mean)
    }
  }
  
  
  #Function to extract point estimate and SE from geex output
  extractEST<-function(geex_output=OM_mest, est_name="m1",param_start=param_start_OM){
    param_num_EST<-match(est_name,names(param_start))
    EST<-geex_output@estimates[param_num_EST]
    sandwich_se <- diag(geex_output@vcov)^0.5 
    SE<-sandwich_se[param_num_EST]
    return(c(EST, SE=SE))
  }
  
  
  generate_weights<-function(Smod,Amod, data){
    S1data<-dplyr::filter(data, S==1)
    w_reg<-glm(Smod, family="binomial", data=data)
    ps<- predict(w_reg,newdata=data, type="response") 
    w_reg2<-glm(Amod, family="binomial", data=S1data)
    pa<- predict(w_reg2,newdata=data, type="response") 
    w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
    data$w<-w
    list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
    return(list)
  }
  
  
  #functions that do not account for estimation of the working models 
  #most appropriate for bootstrap estimation of SEs
  
  OM_est<-function(data){
    S1data_A1<-dplyr::filter(data, S==1 & A==1)
    OM1mod<-glm(formula=Y~.-S-A-w, data=S1data_A1)
    p1<- predict(OM1mod,newdata=data, type="response") 
    data$p1<-p1
    S1data_A0<-dplyr::filter(data, S==1 & A==0)
    OM0mod<-glm(formula=Y~.-S-p1-A-w, data=S1data_A0)
    p0<- predict(OM0mod,newdata=data, type="response") 
    data$p0<-p0
    S0sub<-dplyr::filter(data, S==0)
    OM_1<-mean(S0sub$p1)
    OM_0<-mean(S0sub$p0)
    OM<-mean(S0sub$p1)-mean(S0sub$p0)
    list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
    return(list)
  }
  
  
  IOW1_est<-function(data){ 
    A<-data$A
    S<-data$S
    w<-data$w
    Y<-data$Y
    IOW1_1 <-(sum((1-S))^-1)* sum(A*S*w*Y)
    IOW1_0 <-(sum((1-S))^-1)* sum((1-A)*S*w*Y)
    IOW1 = IOW1_1 - IOW1_0
    return(list(IOW1_1=IOW1_1,IOW1_0=IOW1_0, IOW1=IOW1))
  }
  
  
  IOW2_est<-function(data){
    S0data<-dplyr::filter(data, S==0)  
    S1data_A1<-dplyr::filter(data, S==1 & A==1)
    IOW1mod<-glm(formula=Y~1, data=S1data_A1, weights=w)
    p1<- predict(IOW1mod,newdata=S0data, type="response") 
    S1data_A0<-dplyr::filter(data, S==1 & A==0)
    IOW0mod<-glm(formula=Y~1, data=S1data_A0, weights=w)
    p0<- predict(IOW0mod,newdata=S0data, type="response") 
    IOW2_1<-mean(p1)
    IOW2_0<-mean(p0)
    IOW2<-mean(p1)-mean(p0)
    list<-list(IOW2_1=IOW2_1,IOW2_0=IOW2_0, IOW2=IOW2,IOW1mod=IOW1mod,IOW0mod=IOW0mod)
    return(list)
  }
  
  
  DR1_est<-function(data){  
    A<-data$A
    S<-data$S
    Y<-data$Y
    p1<-data$p1
    p0<-data$p0
    w<-data$w
    DR1_1<-(sum((1-S))^-1)* sum(S*A*w*(Y-p1) + (1-S)*p1)
    DR1_0<-(sum((1-S))^-1)* sum(S*(1-A)*w*(Y-p0) + (1-S)*p0)
    DR1<-DR1_1-DR1_0
    list<-list(DR1_1=DR1_1,DR1_0=DR1_0, DR1=DR1)
    return(list)
  }
  
  
  DR2_est<-function(data){
    A<-data$A
    S<-data$S
    Y<-data$Y
    p1<-data$p1
    p0<-data$p0
    w<-data$w
    sum1_DR2<-sum(S*A*w*(Y-p1)) 
    sum0_DR2<-sum(S*(1-A)*w*(Y-p0)) 
    norm1<-(sum(S*A*w))^-1
    norm0<-(sum(S*(1-A)*w))^-1
    DR2_1<-norm1*sum1_DR2 + (sum(1-S)^-1)*sum((1-S)*p1)
    DR2_0<-norm0*sum0_DR2 + (sum(1-S)^-1)*sum((1-S)*p0)
    DR2<-DR2_1-DR2_0
    list<-list(DR2_1=DR2_1,DR2_0=DR2_0, DR2=DR2)
    return(list)
  }
  
  
  DR3_est<-function(data){
    S0data<-dplyr::filter(data, S==0)
    S1data_A1<-dplyr::filter(data, S==1 & A==1)
    DR1mod<-glm(formula=Y~.-A-S-w-p1-p0, data=S1data_A1, weights=w)
    p1<- predict(DR1mod,newdata=S0data, type="response") 
    S1data_A0<-dplyr::filter(data, S==1 & A==0)
    DR0mod<-glm(formula=Y~.-A-S-w-p1-p0, data=S1data_A0, weights=w)
    p0<- predict(DR0mod,newdata=S0data, type="response") 
    DR3_1<- mean(p1)
    DR3_0<- mean(p0)
    DR3<-mean(p1)-mean(p0) 
    list<-list(DR3_1=DR3_1,DR3_0=DR3_0, DR3=DR3,DR1mod=DR1mod, DR0mod=DR0mod)
    return(list)
  }
}


