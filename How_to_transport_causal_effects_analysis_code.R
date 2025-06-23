

# change these paths to where the other scripts are stored
source("C:/Users/vb3/Desktop/whentohowto/howto_functions.R")
source("C:/Users/vb3/Desktop/T_functions.R")
#require(c(tidyr))
#this function reads in and preprocesses the data used in our example the returned dataset is cleaned
a=data_read()[[5]]

b=subset(a,select=c(Y,Z,S))
names(indo_rct)
# bleed and asa and train due to logic
dat=select(a,-asa81,-asa325,-asa)
#paninj and train due to no risk factors and type due to double info
dat=select(dat,-train,-paninj,-type1,-type2,-type3)
#check missingness, no missings
apply(is.na(dat),FUN=sum,MARGIN=2)
#check for rare cases see table cutoff <3% or less than 10
tableone::CreateTableOne(names(dat),strata = "S",data=dat,factorVars=names(dat)[-c(1,2)])
#removing violation: brush pneudil 
#removing <10 amp acinar chole pbmal status 
dat=select(dat,-brush,-pneudil,-amp,-acinar,-chole,-pbmal,-status)
#correlations drop pdstent
cor_dat=cor(select(dat,-Y,-Z,-S))

#plotting variables in a heatmap to see highly correlated variable pairs
data_melt=melt(round(cor_dat,1))

ggplot(data_melt, aes(Var1, Var2, fill=value)) +
  geom_tile() +geom_text(aes(label=round(value, 2)), color="black", size=3)+
  scale_fill_gradient2(low="blue", high="red", mid="white", 
                       midpoint=0, limit=c(-1,1), space="Lab", 
                       name="Correlation",guide = guide_colourbar(ticks.colour = "black", ticks.linewidth = 1, label.theme = element_text(color = "black")))+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(1.5, "cm"))  +
  coord_fixed()
dat=select(dat,-pdstent)

###code to construct cloudplots 
dat$S=ifelse(dat$S==1,0,1)
{data=dat
data$A=data$Z
data=select(data,-Z)
names=names(select(data,-A,-Y,-S))
names=sort(names)
B = 1000
#threshold 75% of B other threshold 80%
threshold = 0.75 * B
eval(parse(text=paste0("OR_y_mod = data.frame(intercept = rep(0, B),",paste0(names, collapse ="= 0," ),"=0)")))
eval(parse(text=paste0("OR_s_mod = data.frame(intercept = rep(0, B),",paste0(names, collapse ="= 0," ),"=0)")))
library(parameters)


#1000 bootstraps and models for outcome and sample model
{
  my_theme = function() {
    theme(
      # add border
      panel.border = element_rect(colour = "black", fill = NA),
      # axis.line.x = element_line(linetype = "solid", colour = "black"),
      # axis.line.y = element_line(linetype = "solid", colour = "black"),
      # color background
      panel.background = element_blank(),
      # modify grid
      # panel.grid.major.x = element_line(color = "lightgrey"),
      # panel.grid.major.y = element_line(color = "lightgrey"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      # modify text, axis, titles and colour
      plot.title = element_text(colour = "black", size = 20, hjust = 0.5),
      axis.text = element_text(colour = "black", size = 12),
      #axis.text.x = element_text(colour = "black", face = "italic", size = 12,
      #                           angle = 45, vjust = 1, hjust = 1),
      axis.title = element_text(colour = "black", size = 15),
      axis.ticks = element_blank(),
      # legend
      legend.position = "right",
      legend.title = element_text(size = 16),
      legend.title.align = 0.5,
      legend.text = element_text(size = 12), #, margin = margin(t=30, b =30, unit = "cm")),
      #legend.text = element_text(margin = margin(t = 10)))
      #legend.spacing.x = unit(4, 'cm'),
      #legend.key = element_rect(fill="white"),
      #legend.key.size=unit(2, 'cm'),
      #strip text size
      legend.background = element_blank(),
      legend.box.background = element_rect(fill = "white", color = "black"),
      strip.text=element_text(colour = "black", size = 16),
      #strip.background = element_rect(fill = "white", color = "black", vjust = 1.5),
      strip.background = element_blank()
    )
  }
  
}
for (n in (1:B)){
  set.seed(4827+n*3)
  data_boot = data[sample(nrow(data) , nrow(data) , replace = TRUE),]
  #outcome model
  y_mod = standardize_parameters(glm(Y ~ .,
                                     data = select(filter(data_boot, S == 1),-A, -S )),
                                 two_sd = TRUE, ci_method="wald")
  OR_y = data.frame(names = y_mod$Parameter, OR = y_mod$Std_Coefficient)
  OR_y_mod[n,] = pull(arrange(OR_y, names), OR)
  s_mod = standardize_parameters(glm(S ~ ., data = select(data_boot, -Y, -A)),
                                 two_sd = TRUE, ci_method="wald")
  OR_s = data.frame(names = s_mod$Parameter, OR =  s_mod$Std_Coefficient)
  OR_s_mod[n,] = pull(arrange(OR_s, names), OR)
}

OR_y_mod_long = pivot_longer(OR_y_mod, cols = everything(), names_to = "covariate", values_to = "aOR_y_mod")
OR_y_mod_long$id = rep(1:B, each = ncol(OR_y_mod))

OR_s_mod_long = pivot_longer(OR_s_mod, cols = everything(), names_to = "covariate", values_to = "aOR_s_mod")
OR_s_mod_long$id = rep(1:B, each = ncol(OR_s_mod))

OR_table = merge(OR_y_mod_long, OR_s_mod_long, by = c("covariate","id"))
OR_table = filter(OR_table, covariate != "intercept")

OR_table_means = OR_table %>% group_by(covariate) %>% summarize(aOR_y_mod = round(mean(aOR_y_mod),4),
                                                                aOR_s_mod = round(mean(aOR_s_mod),4))
OR_table_means = arrange(OR_table_means, covariate)

rbind(table(OR_table$aOR_y_mod >= 0, OR_table$covariate),
      mean = OR_table_means$aOR_y_mod)

rbind(table(OR_table$aOR_s_mod >= 0, OR_table$covariate),
      mean = OR_table_means$aOR_s_mod)

#significance level = 0.5 --> with a two-sided test this means 75% of points must lie on one side of an axis
threshold = (0.50 + (1-0.5)/2) * B

include_y_mod = data.frame(pivot_longer(
  rownames_to_column(data.frame(rbind(table(OR_table$aOR_y_mod >= 0, OR_table$covariate))),"y_mod_include") ,
  names_to = "covariate", cols = -1))
include_y_mod  = pivot_wider(include_y_mod , names_from = y_mod_include, values_from = value)
include_y_mod = mutate(include_y_mod, include_y_mod = case_when(`FALSE` > threshold | `TRUE` > threshold ~ TRUE, TRUE ~ FALSE))
OR_table_means = merge(OR_table_means, include_y_mod, by= "covariate")
OR_table_means = mutate(OR_table_means, conf_y_mod = case_when(`FALSE` >= `TRUE` ~ `FALSE`, TRUE ~ `TRUE`)/1000)
OR_table_means = select(OR_table_means, covariate, aOR_y_mod, aOR_s_mod, conf_y_mod)

include_s_mod = data.frame(pivot_longer(
  rownames_to_column(data.frame(rbind(table(OR_table$aOR_s_mod >= 0, OR_table$covariate))),"s_mod_include") ,
  names_to = "covariate", cols = -1))
include_s_mod  = pivot_wider(include_s_mod , names_from = s_mod_include, values_from = value)
include_s_mod = mutate(include_s_mod, include_s_mod = case_when(`FALSE` > threshold | `TRUE` > threshold ~ TRUE, TRUE ~ FALSE))
OR_table_means = merge(OR_table_means, include_s_mod, by= "covariate")
OR_table_means = mutate(OR_table_means, conf_s_mod = case_when(`FALSE` >= `TRUE` ~ `FALSE`, TRUE ~ `TRUE`)/1000)
OR_table_means = select(OR_table_means, covariate, aOR_y_mod, aOR_s_mod, conf_y_mod, conf_s_mod)

include = merge(select(include_y_mod, covariate, include_y_mod, "FALSE_y_mod" = `FALSE`,"TRUE_y_mod" =  `TRUE`),
                select(include_s_mod, covariate, include_s_mod, "FALSE_s_mod" =`FALSE`, "TRUE_s_mod" = `TRUE`), by = "covariate")

exclude = include[,"covariate"]

include = filter(include, include_y_mod == TRUE & include_s_mod == TRUE | FALSE_y_mod >= 0.975 * B | FALSE_y_mod >= 0.975 * B)$covariate
exclude = exclude[!exclude %in% include]

include;exclude


OR_table$group[OR_table$covariate %in% include] = "Include"
OR_table$group[OR_table$covariate %in% exclude] = "Exclude"
OR_table_means$group[OR_table_means$covariate %in% include] = "Include"
OR_table_means$group[OR_table_means$covariate %in% exclude] = "Exclude"
OR_table_means

OR_table = rbind(OR_table, list("\nInclude", 1, 0,0,"Include"),
                 list("\nExclude", 1, 0,0,"Exclude"))

OR_table$covariate = factor(OR_table$covariate, levels = c(include, exclude, "\nInclude", "\nExclude" ))


#colors for plot
brewer_colors = c(brewer.pal(n = 12, name = "Paired"), "grey80", "grey40")
colors = data.frame(covariate = unique(OR_table$covariate), colors = rep("black",length(unique(OR_table$covariate))))
for (i in (1:nrow(colors))){
  colors$colors[i] = case_when(colors$covariate[i] %in% c("\nInclude", "\nExclude") ~ "black",
                               TRUE ~ brewer_colors[i])
}
#colors$colors[is.na(colors$colors)] = "NULL"

colors$covariate = factor(colors$covariate, levels = c("\nInclude",include, "\nExclude",exclude))
colors = arrange(colors,covariate)

OR_table$covariate = factor(OR_table$covariate, levels = c("\nInclude",include, "\nExclude",exclude))
# OR_table_means$covariate = factor(OR_table_means$covariate, levels = c("Include",include, "Exclude",exclude))

# plot full

a = ggplot(OR_table) +
  geom_point(aes(x = aOR_s_mod, y = aOR_y_mod,
                 color = covariate, fill = covariate), shape = 21)+#, alpha = 0.1) +
  geom_point(data = OR_table_means, aes(x = aOR_s_mod, y = aOR_y_mod), color = "black")+
  geom_vline(xintercept = 0, color = "lightgrey")+
  geom_hline(yintercept = 0, color = "lightgrey")+
  # scale_x_log10(limits = c(min(OR_table$aOR_s_mod), max(OR_table$aOR_s_mod)), name = "aOR sampling model")+
  # scale_y_log10(limits = c(min(OR_table$aOR_y_mod), max(OR_table$aOR_y_mod)), name = "aOR outcome model")+
  scale_x_continuous(limits = c(-max(abs(c(OR_table$aOR_s_mod, OR_table$aOR_y_mod ))),
                                max(abs(c(OR_table$aOR_s_mod, OR_table$aOR_y_mod )))), name = "Association sampling model",
                     expand = c(0.01,0.01))+
  scale_y_continuous(limits = c(-max(abs(c(OR_table$aOR_s_mod, OR_table$aOR_y_mod ))),
                                max(abs(c(OR_table$aOR_s_mod, OR_table$aOR_y_mod )))), name = "Association outcome model",
                     expand = c(0.01,0.01))+
  
  #scale_y_continuous(limits = c(-0.2, 0.2), name = "OR outcome model")+
  #scale_color_brewer(palette = "Paired", name = "Covariate")+
  #scale_fill_brewer(palette = "Paired", name = "Covariate")+
  scale_color_manual(name = "Covariate", values = alpha(colors$colors, case_when(colors$covariate %in% c("\nInclude", "\nExclude") ~ 0,
                                                                                 TRUE ~ 1)))+
  scale_fill_manual(name = "Covariate", values = alpha(colors$colors, case_when(colors$covariate %in% c("\nInclude", "\nExclude") ~ 0,
                                                                                TRUE ~ 0.4)))+
  #scale_alpha_manual(name = "Covariate", values = alpha(c(colors[1:(ncol(OR_y_mod)-1)]), 0.05))+
  guides(color = guide_legend(override.aes = list(size = 3, stroke = 2))) +
  facet_wrap(~factor(group, levels=c("Include", "Exclude")))+
    geom_vline(xintercept = 0, color = "lightgrey")+
  geom_hline(yintercept = 0, color = "lightgrey")+
  my_theme()
a

ggsave("C:/Users/vb3/Desktop/Indomethacin.jpeg",
      units = "px", width = 2400, height = 1800)
}
my_theme = function() {
  theme(
    # add border
    panel.border = element_blank(),
    axis.line.x = element_line(linetype = "solid", colour = "black",  linewidth = 1.5),
    axis.line.y = element_line(linetype = "solid", colour = "black",  linewidth = 1.5),
    # color background
    panel.background = element_blank(),
    # modify grid
    panel.grid.major.x = element_line(color = "lightgrey"),
    panel.grid.major.y = element_line(color = "lightgrey"),
    # modify text, axis, titles and colour
    plot.title = element_text(colour = "black", size = 20, hjust = 0.5),
    axis.text = element_text(colour = "black", size = 12),
    #axis.text.x = element_text(colour = "black", face = "italic", size = 12,
    #                           angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(colour = "black", size = 15),
    axis.ticks = element_blank(),
    # legend
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.title.align = 0.5,
    legend.text = element_text(size = 12), #, margin = margin(t=30, b =30, unit = "cm")),
    #legend.text = element_text(margin = margin(t = 10)))
    #legend.spacing.x = unit(4, 'cm'),
    #legend.key = element_rect(fill="white"),
    #legend.key.size=unit(2, 'cm'),
    #strip text size
    legend.background = element_blank(),
    legend.box.background = element_rect(fill = "white", color = "black"),
    strip.text=element_text(colour = "black", size = 16)
  )
}
set=include
# now run transport 
# main analysis estimation equations '75%  95% @Ymodel

{
  DF=subset(data,select=c("Y","A","S",set))
  library(geex)
  
  DF$Y[is.na(DF$Y)] <- 0
  DF$A[is.na(DF$A)] <- 0
  DF$epsilon<-NULL
  DF$inter<-NULL
  DF$d_pop<-NULL
  DF$d_out<-NULL
  DF$Y_1<-NULL
  DF$Y_0<-NULL
  ####Smod Kovariaten modellieren S(sample dummy)
  ###Amod Kovariaten modellieren A(Treatment dummy)
  
  
  
  weights<-generate_weights(Smod=S~.-A-Y,Amod=A~.-S-Y, data=DF)
  DF2<-weights$dat
  OM<-OM_est(data=DF)
  #-epsilon-inter-d_pop-d_out-Y_1-Y_0
  DF2$p1<-OM$p1
  DF2$p0<-OM$p0
  
  
  
  
  IOW1<-IOW1_est(data=DF2)
  IOW2<-IOW2_est(data=DF2)
  
  DR1<-DR1_est(data=DF2)
  DR2<-DR2_est(data=DF2)
  DR3<-DR3_est(data=DF2)
  
  
  results_1<-rbind(OM_1=OM$OM_1, IOW1_1=IOW1$IOW1_1, IOW2_1=IOW2$IOW2_1,
                   DR1_1=DR1$DR1_1, DR2_1=DR2$DR2_1,  DR3_1=DR3$DR3_1)
  #print(results_1)
  
  #Estimates of mu(0)
  results_0<-rbind(OM_0=OM$OM_0, IOW1_0=IOW1$IOW1_0, IOW2_0=IOW2$IOW2_0,
                   DR1_0=DR1$DR1_0, DR2_0=DR2$DR2_0,  DR3_0=DR3$DR3_0)
  # print(results_0)
  
  #Estimates of ATE
  results<-rbind(OM=OM$OM, IOW1=IOW1$IOW1, IOW2=IOW2$IOW2,
                 DR1=DR1$DR1, DR2=DR2$DR2,  DR3=DR3$DR3)
  # print(results)
  
  
  param_start_OM<-c(coef(OM$OM1mod), coef(OM$OM0mod),
                    m1=OM$OM_1, m0=OM$OM_0,ate=OM$OM)
  
  OM_mest<-m_estimate(
    estFUN = OM_EE,
    data  = DF,
    root_control = setup_root_control(start = param_start_OM),
    compute_roots = T,
    compute_vcov = T
  )
  
  #return potential outcomes means and their difference (ate) and corresponding standard errors
  OM_m1<-extractEST(geex_output=OM_mest, est_name="m1",param_start=param_start_OM)
  OM_m0<-extractEST(geex_output=OM_mest, est_name="m0",param_start=param_start_OM)
  OM_ate<-extractEST(geex_output=OM_mest, est_name="ate",param_start=param_start_OM)
  


  (param_start_IOW1<-c(coef(weights$Smod) , coef(weights$Amod),
                       m1=IOW1$IOW1_1, m0=IOW1$IOW1_0, ate=IOW1$IOW1) )

  IOW1_mest <-m_estimate(
    estFUN = IOW1_EE,
    data  = DF,
    root_control = setup_root_control(start = param_start_IOW1),
    compute_roots = T,
    compute_vcov = T
  )

  #save variance + SE
  IOW2_m1<-extractEST(geex_output=IOW1_mest, est_name="m1",param_start=param_start_IOW1)
  IOW2_m0<-extractEST(geex_output=IOW1_mest, est_name="m0",param_start=param_start_IOW1)
  IOW2_ate<-extractEST(geex_output=IOW1_mest, est_name="ate",param_start=param_start_IOW1)


  
  #DR1
  
  (coef_DR1est<-c(coef(OM$OM1mod), coef(OM$OM0mod), m1=DR1$DR1_1, m0=DR1$DR1_0, ate=DR1$DR1))
  
  param_start_DR1<-c(coef(weights$Smod) , coef(weights$Amod), coef_DR1est)
  
  DR1_mest<-m_estimate(
    estFUN = DR1_EE,
    data  = DF,
    root_control = setup_root_control(start = param_start_DR1),
    compute_roots = T,
    compute_vcov = T
  )
  #save variance + SE
  DR1_m1<-extractEST(geex_output=DR1_mest, est_name="m1",param_start=param_start_DR1)
  DR1_m0<-extractEST(geex_output=DR1_mest, est_name="m0",param_start=param_start_DR1)
  DR1_ate<-extractEST(geex_output=DR1_mest, est_name="ate",param_start=param_start_DR1)
  #DR2
  
  param_start_DR2<-c(coef(weights$Smod), coef(weights$Amod), 0.5, 0.5, 0.5,
                     coef(OM$OM1mod), coef(OM$OM0mod),
                     m1=DR2$DR2_1, m0=DR2$DR2_0, ate=DR2$DR2)
  
  DR2_mest<-m_estimate(
    estFUN = DR2_EE,
    data  = DF,
    root_control = setup_root_control(start = param_start_DR2),
    compute_roots = T,
    compute_vcov = T
  ) 
  
  #save variance + SE
  DR2_m1<-extractEST(geex_output=DR2_mest, est_name="m1",param_start=param_start_DR2)
  DR2_m0<-extractEST(geex_output=DR2_mest, est_name="m0",param_start=param_start_DR2)
  DR2_ate<-extractEST(geex_output=DR2_mest, est_name="ate",param_start=param_start_DR2)
  
  
  #DR3
  
  param_start_DR3<-c(coef(weights$Smod) , coef(weights$Amod), 
                     coef(DR3$DR1mod), coef(DR3$DR0mod),
                     m1=DR3$DR3_1, m0=DR3$DR3_0, ate=DR3$DR3)
  
  DR3_mest<-m_estimate(
    estFUN = DR3_EE,
    data  = DF,
    root_control = setup_root_control(start = param_start_DR3),
    compute_roots = T,
    compute_vcov = T
  ) 
  DR3_m1<-extractEST(geex_output=DR3_mest, est_name="m1",param_start=param_start_DR3)
  DR3_m0<-extractEST(geex_output=DR3_mest, est_name="m0",param_start=param_start_DR3)
  DR3_ate<-extractEST(geex_output=DR3_mest, est_name="ate",param_start=param_start_DR3)
  
output=list(set=set,est=data.frame(rbind(OM_ate,IOW2_ate,DR1_ate,DR2_ate,DR3_ate)))
}
# calculation of the TATE and SATE using bootstrap intervals
{####TATE BOOTSTRAPPING
  B=10000
  TATE=numeric(B)
  for (i in 1:B){
    indices <- sample(rownames(data),size=nrow(data),replace=T)
    #sample_data<-source[1:dim(source)[1] %in% indices,]
    sample_data<-data[indices,]
    data_set=sample_data
    dat=data_set
    TATE[i]=mean(dat$Y[dat$A==1&dat$S==0])-mean(dat$Y[dat$A==0&dat$S==0])
  }
  ####SATE BOOTSTRAPPING
  B=10000
  SATE=numeric(B)
  for (i in 1:B){
    indices <- sample(rownames(data),size=nrow(data),replace=T)
    #sample_data<-source[1:dim(source)[1] %in% indices,]
    sample_data<-data[indices,]
    data_set=sample_data
    dat=data_set
    SATE[i]=mean(dat$Y[dat$A==1&dat$S==1])-mean(dat$Y[dat$A==0&dat$S==1])
  }
  
  SATE_point=mean(data$Y[data$A==1&data$S==1])-mean(data$Y[data$A==0&data$S==1])
  SATE_upper=mean(data$Y[data$A==1&data$S==1])-mean(data$Y[data$A==0&data$S==1])+sd(SATE)*1.96
  SATE_lower=mean(data$Y[data$A==1&data$S==1])-mean(data$Y[data$A==0&data$S==1])-sd(SATE)*1.96
  
  TATE_point=mean(data$Y[data$A==1&data$S==0])-mean(data$Y[data$A==0&data$S==0])
  TATE_upper=mean(data$Y[data$A==1&data$S==0])-mean(data$Y[data$A==0&data$S==0])+sd(TATE)*1.96
  TATE_lower=mean(data$Y[data$A==1&data$S==0])-mean(data$Y[data$A==0&data$S==0])-sd(TATE)*1.96
  
}

# plot the results
{
outdat=mutate(output$est,
              CI95_lb= ate-1.96*SE,CI95_ub=ate+1.96*SE,
              ATE=ate,
              Estimator=gsub("_ate","",rownames(output$est)))#,
              #Covariates=paste(set,sep=""))
sate_tate=rbind(c(0,0,TATE_lower,TATE_upper,TATE_point,"TATE"),
                c(0,0,SATE_lower,SATE_upper,SATE_point,"SATE"))
colnames(sate_tate)=names(outdat)
data_for_plot=rbind(outdat,sate_tate)
data4plot=data_for_plot[,-c(1,2,7)]
data4plot=mutate(data4plot,CI95_lb=round(as.numeric(CI95_lb),3),CI95_ub=round(as.numeric(CI95_ub),3),ATE=round(as.numeric(ATE),3))

library(ggplot2)

my_theme_forest = function() {
  theme(
    # add border
    panel.border = element_blank(),
    #axis.line.x = element_line(linetype = "solid", colour = "black",  linewidth = 1.5),
    axis.line.y = element_line(linetype = "solid", colour = "black",  linewidth = 1),
    # color background
    panel.background = element_blank(),
    # modify grid
    panel.grid.major.y = element_line(linetype = "solid", color = "lightgrey", linewidth = 0.7),
    panel.grid.major.x = element_blank(),
    # modify text, axis, titles and colour
    plot.title = element_text(colour = "black", size = 30, hjust = 0.5),
    axis.text.y = element_text(colour = "black", size = 12),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    axis.ticks = element_blank(),
    # legend
    legend.position = "right",
    #    legend.position = "inside",
    #    legend.position.inside = c(0.75,0.85),
    #    legend.direction = "horizontal",
    legend.title = element_text(size = 16, hjust = 0.5),
    #legend.title.position = "top",
    legend.text = element_text(size = 14, margin = margin(t=0.3, b =0.3, unit = "cm")),
    #legend.key = element_rect(fill="white", color = "white"),
    legend.background = element_blank(),
    legend.box.background = element_rect(fill = "white", color = "black"),
    legend.key.spacing.y = unit(0.2, 'cm'),
    legend.key.spacing.x = unit(0.5, 'cm'),
    #legend.spacing.y = unit(0.3, 'cm'),
    #strip text size
    strip.text=element_text(colour = "black", size = 16),
    strip.background = element_rect(fill = "white", color = "black")
  )
}

how_to_forest = data4plot

how_to_forest$Group = case_when(how_to_forest$Estimator == "SATE" ~ " Indiana\n original",
                                how_to_forest$Estimator == "TATE" ~ " Michigan\n original",
                                TRUE ~ " Transported")
how_to_forest$Group = factor(how_to_forest$Group, levels = c(" Indiana\n original", " Transported",
                                                             " Michigan\n original"))
how_to_forest$Width = case_when(how_to_forest$Estimator %in% c("SATE", "TATE") ~ 0.1,
                                TRUE ~ 0.5)

how_to_forest$Estimator = factor(how_to_forest$Estimator, levels = c("SATE", "OM", "IOW2", "DR1", "DR2", "DR3", "TATE"))


a = ggplot(how_to_forest) +
  #geom_hline(yintercept = 0, linetype = "solid") +
  geom_errorbar(aes(x=Estimator, y = ATE, ymax=CI95_ub, ymin=CI95_lb,
                    #shape = Group,
                    width = Width),
                size = 1, position=position_dodge(width = 0.9))+#, shape = estimator),  position=position_dodge2(width = 0.9, reverse = TRUE), size = 0.7)
  geom_point(aes(x = Estimator, y = ATE, shape = Group), size = 4, stroke = 1.6,
             position=position_dodge(width = 0.9))+
  xlab("")+
  ylab("Risk difference")+
  #guides(shape=guide_legend(title="Sample"))+
  scale_shape_manual(values = c(0,7,15), name = "Estimate") +
  guides(shape = "none") +
  # scale_y_continuous(limits = c(-max(abs(c(how_to_forest$CI95_lb, how_to_forest$CI95_ub ))),
  #                                max(abs(c(how_to_forest$CI95_lb, how_to_forest$CI95_ub )))))+
  scale_y_continuous(limits = c(-0.5,0.1))+
  
  #facet_grid(estimate~.)+
  facet_wrap(~Group, scales = "free_x")#+
  #coord_flip() +
 # my_theme_forest()
a

#ggsave("C:/Users/vb3/Desktop/Indomethacin_forest.jpeg",
 #      units = "px", width = 2400, height = 1800)
}


