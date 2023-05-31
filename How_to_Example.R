#example

source("howto_functions.R")
#include the Indomethacin trial dataset

sol=data_read()
#the Population Dummy
S=sol[[1]];
#the Output data
Y=sol[[2]];
#the treatment variable
Z=sol[[3]];
#the matrix of covariates
X=sol[[4]];
#the dataset as whole
data=sol[[5]];

#input the DAG that is assumed
SD_input=dagitty("dag{
  Riskfactors ->Y;
  Aspirin-> Bleed;
  X->Bleed;
  Bleed->Y;
  S-> Riskfactors;
  X -> Y;
  S-> Aspirin;
  Aspirin->Bleed;
  Riskfactors<->Aspirin;
  X [exposure]
  Y [outcome]
}")

#print the DAG for visual inspection
plot(SD_input)
data="p(Riskfactors,Y,Aspirin|do(X),S)
         p(Riskfactors,Aspirin)"

#generate and display the transport formula
res=dosearch(data,
             query= "p(Y|do(X))",graph=SD_input,transportability = "S")
res
#run the omnibustest (only applicable for the validate rct purpose)
run_omnibus()

#seed for comparability
set.seed(1)
data=sol[[5]];
lib=c("SL.glm","SL.glm.interaction","SL.gam")

#Bootstrap results for each purpose
purpose="validate"
g=run_bootstrap(R=100,purpose=purpose,lib=lib,method="TMLE",data=data)
print(g)

purpose="gennew"
g=run_bootstrap(R=100,purpose=purpose,lib=lib,method="TMLE",data=data)
print(g)

purpose="overcome"
set.seed(1)
simdat=simulate_confounding_example_data()
data=data.frame(cbind(X=simdat$X,Y=simdat$Y,S=simdat$S,Z=simdat$Z))
g=run_bootstrap(R=100,purpose=purpose,lib=lib,method="TMLE",data=data)

CATE=mean(data$Y[data$S==0&data$Z==1])-mean(data$Y[data$S==0&data$Z==0])
SATE=mean(data$Y[data$S==1&data$Z==1])-mean(data$Y[data$S==1&data$Z==0])
TATE=mean(simdat$t_outcome_rand[simdat$t_rand==1])-mean(simdat$t_outcome_rand[simdat$t_rand==0])
print(g);CATE;SATE;TATE;


