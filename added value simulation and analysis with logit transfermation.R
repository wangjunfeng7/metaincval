memory.limit()

rm(list=ls())


library(MASS)
library(pROC)
library(rms)
library(metafor)
library(metamisc)


library(ggplot2)
library(reshape2)

library(clusterGeneration)
library(glmnet)

## Overall settings

set.seed(12345)

# Number of simulation
n_simulation <- 100

# Number of predictors (6 in the model, 1 for added value)
J  <- 7                    

# Totoal population
N_population_total <- 10^6
N_population <- 10^5

# Max number of studies
n_study <- c(5,10,20)
k_max <-  max(n_study)


# Sample size in each study

n200 <- round(rnorm(n = k_max, mean = 200, sd = 0.1*200))
n500 <- round(rnorm(n = k_max, mean = 500, sd = 0.1*500))
n1000 <- round(rnorm(n = k_max, mean = 1000, sd = 0.1*1000))







# performance[[1]]$result <- c(1,2,3,4,5)


## Data generation

# mean and correlation (rho = 0, 0.2, 0.4)

mu <- as.vector(matrix(0, 1, J)) 

sigma_0 <- matrix(0,nrow=J,ncol=J) 
diag(sigma_0) <- rep(1,J)


sigma_0.2 <- matrix(0.2,nrow=J,ncol=J) 
diag(sigma_0.2) <- rep(1,J)


sigma_0.4 <- matrix(0.4,nrow=J,ncol=J) 
diag(sigma_0.4) <- rep(1,J)


# coefficients

beta_x1 <- 0.4
beta_x2 <- 0.7
beta_x3 <- 0.9
beta_x4 <- 0.4
beta_x5 <- 0.7
beta_x6 <- 0.9


# coefficients (lower reference AUC)

beta_x1 <- 0.4
beta_x2 <- 0.4
beta_x3 <- 0.4
beta_x4 <- 0.4
beta_x5 <- 0.4
beta_x6 <- 0.4

# effect of added predictor (beta = 0, 0.4, 0.7, 0.9)
beta_x7_all <- c(0, 0.4, 0.7, 0.9)






####### Simulation starts here


# generate some lists to save results

performance <- vector(length(beta_x7_all),mode="list")
result_200 <- vector(n_simulation,mode="list")
# result[[1]] <- vector(length(beta_x7_all),mode="list")
result_500 <- vector(n_simulation,mode="list")
result_1000 <- vector(n_simulation,mode="list")



## get reference values




# only run one time in a large population as the reference performence


for (i in 1:1){
 #i <- 1
  
  # Simulation per correlation matrix
  
  data <- as.data.frame(mvrnorm(n=N_population_total, mu=mu, Sigma=sigma_0))
  
  names(data) <- c("X1", "X2","X3", "X4","X5", "X6","X7")
  
  # Dichotomize X4-X6
  
  data$X4_bin <-  ifelse(data$X4 < 0,0,1)
  data$X5_bin <-  ifelse(data$X5 < 0,0,1)
  data$X6_bin <-  ifelse(data$X6 < 0,0,1)
  
  # Dichotomize X7
  data$X7_bin <-  ifelse(data$X7 < 0,0,1)
  
  
  for (j in (length(beta_x7_all)):1){
    
    # define effect of X7
     # j <- 1
    beta_x7 <- beta_x7_all[j]
    
    
    # calculate lp and p
    
    data$lp_con <- -(beta_x4+beta_x5+beta_x6)*0.5 + beta_x1*data$X1 + beta_x2*data$X2 + beta_x3*data$X3 + beta_x4*data$X4_bin +
      beta_x5*data$X5_bin + beta_x6*data$X6_bin + beta_x7*data$X7
    
    data$lp_bin <- -(beta_x4+beta_x5+beta_x6)*0.5-0.5*beta_x7 + beta_x1*data$X1 + beta_x2*data$X2 + beta_x3*data$X3 + beta_x4*data$X4_bin +
      beta_x5*data$X5_bin + beta_x6*data$X6_bin + beta_x7*data$X7_bin
    
    
    
    
    data$p_con <- exp(data$lp_con)/(1+exp(data$lp_con))
    data$p_bin <- exp(data$lp_bin)/(1+exp(data$lp_bin))
    
    # generate outcome
    
    data$outcome_con <- rbinom(n = N_population_total, 1, data$p_con)
    data$outcome_bin <- rbinom(n = N_population_total, 1, data$p_bin)
    
    
    # calculate oracle performance (true model)
    
    performance[[j]]$auc_con_oracle[i] <- roc(outcome_con ~ p_con, data)$auc
    performance[[j]]$auc_bin_oracle[i] <- roc(outcome_bin ~ p_bin, data)$auc
    
    
    # refit model
    
    model_6_con <- glm(outcome_con~X1+X2+X3+X4_bin+X5_bin+X6_bin,data=data,family=binomial(link="logit"))
    model_6_bin <- glm(outcome_bin~X1+X2+X3+X4_bin+X5_bin+X6_bin,data=data,family=binomial(link="logit"))
    
    model_7_con <- glm(outcome_con~X1+X2+X3+X4_bin+X5_bin+X6_bin+X7,data=data,family=binomial(link="logit"))
    model_7_bin <- glm(outcome_bin~X1+X2+X3+X4_bin+X5_bin+X6_bin+X7_bin,data=data,family=binomial(link="logit"))
    
    
    # make prediction
    
    
    data$p_predicted_6_con <- predict(model_6_con, data, type="response")
    data$p_predicted_6_bin <- predict(model_6_bin, data, type="response")
    data$p_predicted_7_con <- predict(model_7_con, data, type="response")
    data$p_predicted_7_bin <- predict(model_7_bin, data, type="response")
    
    # calculate oracle performance (fitted model)
    
    performance[[j]]$auc_6_con[i] <- roc(outcome_con ~ p_predicted_6_con, data)$auc
    performance[[j]]$auc_6_bin[i] <- roc(outcome_bin ~ p_predicted_6_bin, data)$auc
    
    performance[[j]]$auc_7_con[i] <- roc(outcome_con ~ p_predicted_7_con, data)$auc
    performance[[j]]$auc_7_bin[i] <- roc(outcome_bin ~ p_predicted_7_bin, data)$auc
    
  }}
    












## simulate primary studies


set.seed(123)

for (i in 1:n_simulation){
  #i <- 1

# Simulation per correlation matrix

data <- as.data.frame(mvrnorm(n=N_population, mu=mu, Sigma=sigma_0))

names(data) <- c("X1", "X2","X3", "X4","X5", "X6","X7")

# Dichotomize X4-X6

data$X4_bin <-  ifelse(data$X4 < 0,0,1)
data$X5_bin <-  ifelse(data$X5 < 0,0,1)
data$X6_bin <-  ifelse(data$X6 < 0,0,1)

# Dichotomize X7
data$X7_bin <-  ifelse(data$X7 < 0,0,1)


for (j in (length(beta_x7_all)):1){

# define effect of X7
# j <- 2
beta_x7 <- beta_x7_all[j]


# calculate lp and p

data$lp_con <- -(beta_x4+beta_x5+beta_x6)*0.5 + beta_x1*data$X1 + beta_x2*data$X2 + beta_x3*data$X3 + beta_x4*data$X4_bin +
  beta_x5*data$X5_bin + beta_x6*data$X6_bin + beta_x7*data$X7

data$lp_bin <- -(beta_x4+beta_x5+beta_x6)*0.5-0.5*beta_x7 + beta_x1*data$X1 + beta_x2*data$X2 + beta_x3*data$X3 + beta_x4*data$X4_bin +
  beta_x5*data$X5_bin + beta_x6*data$X6_bin + beta_x7*data$X7_bin
  
  
  

data$p_con <- exp(data$lp_con)/(1+exp(data$lp_con))
data$p_bin <- exp(data$lp_bin)/(1+exp(data$lp_bin))

# generate outcome

data$outcome_con <- rbinom(n = N_population, 1, data$p_con)
data$outcome_bin <- rbinom(n = N_population, 1, data$p_bin)






# fit model in each primary study

for (k in 1:k_max){
  #k <- 1
  
# analysis for small sample size (n~200)  
  
  data_primary_200 <- data[((k-1)*N_population/k_max+1):((k-1)*N_population/k_max+n200[k]),c("X1","X2","X3", "X4_bin", "X5_bin", "X6_bin", "X7", "X7_bin" ,"outcome_con","outcome_bin")]
  # refit model
  
  model_6_con_primary_200 <- lrm(outcome_con~X1+X2+X3+X4_bin+X5_bin+X6_bin,data=data_primary_200)
  model_6_bin_primary_200 <- lrm(outcome_bin~X1+X2+X3+X4_bin+X5_bin+X6_bin,data=data_primary_200)
  
  model_7_con_primary_200 <- lrm(outcome_con~X1+X2+X3+X4_bin+X5_bin+X6_bin+X7,data=data_primary_200)
  model_7_bin_primary_200 <- lrm(outcome_bin~X1+X2+X3+X4_bin+X5_bin+X6_bin+X7_bin,data=data_primary_200)
  
  
  # make prediction
  
  
  data_primary_200$p_predicted_6_con <- predict(model_6_con_primary_200, data_primary_200, type="fitted")
  data_primary_200$p_predicted_6_bin <- predict(model_6_bin_primary_200, data_primary_200, type="fitted")
  data_primary_200$p_predicted_7_con <- predict(model_7_con_primary_200, data_primary_200, type="fitted")
  data_primary_200$p_predicted_7_bin <- predict(model_7_bin_primary_200, data_primary_200, type="fitted")
  
  # calculate oracle performance (fitted model)
  
  result_200[[i]][[j]]$auc_6_con[k] <- roc(outcome_con ~ p_predicted_6_con, data_primary_200)$auc
  result_200[[i]][[j]]$var_auc_6_con[k] <- var(roc(outcome_con ~ p_predicted_6_con, data_primary_200))
  
  result_200[[i]][[j]]$auc_6_bin[k] <- roc(outcome_bin ~ p_predicted_6_bin, data_primary_200)$auc
  result_200[[i]][[j]]$var_auc_6_bin[k] <- var(roc(outcome_bin ~ p_predicted_6_bin, data_primary_200))
  
  result_200[[i]][[j]]$auc_7_con[k] <- roc(outcome_con ~ p_predicted_7_con, data_primary_200)$auc
  result_200[[i]][[j]]$var_auc_7_con[k] <- var(roc(outcome_con ~ p_predicted_7_con, data_primary_200))
  
  result_200[[i]][[j]]$auc_7_bin[k] <- roc(outcome_bin ~ p_predicted_7_bin, data_primary_200)$auc
  result_200[[i]][[j]]$var_auc_7_bin[k] <- var(roc(outcome_bin ~ p_predicted_7_bin, data_primary_200))
  
  #roc.test(roc(outcome_con ~ p_predicted_6_con, data_primary_200), roc(outcome_con ~ p_predicted_7_con, data_primary_200))
  
  result_200[[i]][[j]]$cov_auc_con[k] <- result_200[[i]][[j]]$var_auc_7_con[k] + result_200[[i]][[j]]$var_auc_6_con[k] - 2*cov(roc(outcome_con ~ p_predicted_6_con, data_primary_200), roc(outcome_con ~ p_predicted_7_con, data_primary_200))
  result_200[[i]][[j]]$cov_auc_bin[k] <- result_200[[i]][[j]]$var_auc_7_bin[k] + result_200[[i]][[j]]$var_auc_6_bin[k] - 2*cov(roc(outcome_bin ~ p_predicted_6_bin, data_primary_200), roc(outcome_bin ~ p_predicted_7_bin, data_primary_200))
  
  
  
  
  # analysis for medium sample size (n~500)  
  
  data_primary_500 <- data[((k-1)*N_population/k_max+1):((k-1)*N_population/k_max+n500[k]),c("X1","X2","X3", "X4_bin", "X5_bin", "X6_bin", "X7", "X7_bin" ,"outcome_con","outcome_bin")]
  # refit model
  
  model_6_con_primary_500 <- lrm(outcome_con~X1+X2+X3+X4_bin+X5_bin+X6_bin,data=data_primary_500)
  model_6_bin_primary_500 <- lrm(outcome_bin~X1+X2+X3+X4_bin+X5_bin+X6_bin,data=data_primary_500)
  
  model_7_con_primary_500 <- lrm(outcome_con~X1+X2+X3+X4_bin+X5_bin+X6_bin+X7,data=data_primary_500)
  model_7_bin_primary_500 <- lrm(outcome_bin~X1+X2+X3+X4_bin+X5_bin+X6_bin+X7_bin,data=data_primary_500)
  
  
  # make prediction
  
  
  data_primary_500$p_predicted_6_con <- predict(model_6_con_primary_500, data_primary_500, type="fitted")
  data_primary_500$p_predicted_6_bin <- predict(model_6_bin_primary_500, data_primary_500, type="fitted")
  data_primary_500$p_predicted_7_con <- predict(model_7_con_primary_500, data_primary_500, type="fitted")
  data_primary_500$p_predicted_7_bin <- predict(model_7_bin_primary_500, data_primary_500, type="fitted")
  
  # calculate oracle performance (fitted model)
  
  result_500[[i]][[j]]$auc_6_con[k] <- roc(outcome_con ~ p_predicted_6_con, data_primary_500)$auc
  result_500[[i]][[j]]$var_auc_6_con[k] <- var(roc(outcome_con ~ p_predicted_6_con, data_primary_500))
  
  result_500[[i]][[j]]$auc_6_bin[k] <- roc(outcome_bin ~ p_predicted_6_bin, data_primary_500)$auc
  result_500[[i]][[j]]$var_auc_6_bin[k] <- var(roc(outcome_bin ~ p_predicted_6_bin, data_primary_500))
  
  result_500[[i]][[j]]$auc_7_con[k] <- roc(outcome_con ~ p_predicted_7_con, data_primary_500)$auc
  result_500[[i]][[j]]$var_auc_7_con[k] <- var(roc(outcome_con ~ p_predicted_7_con, data_primary_500))
  
  result_500[[i]][[j]]$auc_7_bin[k] <- roc(outcome_bin ~ p_predicted_7_bin, data_primary_500)$auc
  result_500[[i]][[j]]$var_auc_7_bin[k] <- var(roc(outcome_bin ~ p_predicted_7_bin, data_primary_500))
  
  #roc.test(roc(outcome_con ~ p_predicted_6_con, data_primary_500), roc(outcome_con ~ p_predicted_7_con, data_primary_500))
  
  result_500[[i]][[j]]$cov_auc_con[k] <- result_500[[i]][[j]]$var_auc_7_con[k] + result_500[[i]][[j]]$var_auc_6_con[k] - 2*cov(roc(outcome_con ~ p_predicted_6_con, data_primary_500), roc(outcome_con ~ p_predicted_7_con, data_primary_500))
  result_500[[i]][[j]]$cov_auc_bin[k] <- result_500[[i]][[j]]$var_auc_7_bin[k] + result_500[[i]][[j]]$var_auc_6_bin[k] - 2*cov(roc(outcome_bin ~ p_predicted_6_bin, data_primary_500), roc(outcome_bin ~ p_predicted_7_bin, data_primary_500))
  
  
  
  
  # analysis for large sample size (n~1000)  
  
  data_primary_1000 <- data[((k-1)*N_population/k_max+1):((k-1)*N_population/k_max+n1000[k]),c("X1","X2","X3", "X4_bin", "X5_bin", "X6_bin", "X7", "X7_bin" ,"outcome_con","outcome_bin")]
  # refit model
  
  model_6_con_primary_1000 <- lrm(outcome_con~X1+X2+X3+X4_bin+X5_bin+X6_bin,data=data_primary_1000)
  model_6_bin_primary_1000 <- lrm(outcome_bin~X1+X2+X3+X4_bin+X5_bin+X6_bin,data=data_primary_1000)
  
  model_7_con_primary_1000 <- lrm(outcome_con~X1+X2+X3+X4_bin+X5_bin+X6_bin+X7,data=data_primary_1000)
  model_7_bin_primary_1000 <- lrm(outcome_bin~X1+X2+X3+X4_bin+X5_bin+X6_bin+X7_bin,data=data_primary_1000)
  
  
  # make prediction
  
  
  data_primary_1000$p_predicted_6_con <- predict(model_6_con_primary_1000, data_primary_1000, type="fitted")
  data_primary_1000$p_predicted_6_bin <- predict(model_6_bin_primary_1000, data_primary_1000, type="fitted")
  data_primary_1000$p_predicted_7_con <- predict(model_7_con_primary_1000, data_primary_1000, type="fitted")
  data_primary_1000$p_predicted_7_bin <- predict(model_7_bin_primary_1000, data_primary_1000, type="fitted")
  
  # calculate oracle performance (fitted model)
  
  result_1000[[i]][[j]]$auc_6_con[k] <- roc(outcome_con ~ p_predicted_6_con, data_primary_1000)$auc
  result_1000[[i]][[j]]$var_auc_6_con[k] <- var(roc(outcome_con ~ p_predicted_6_con, data_primary_1000))
  
  result_1000[[i]][[j]]$auc_6_bin[k] <- roc(outcome_bin ~ p_predicted_6_bin, data_primary_1000)$auc
  result_1000[[i]][[j]]$var_auc_6_bin[k] <- var(roc(outcome_bin ~ p_predicted_6_bin, data_primary_1000))
  
  result_1000[[i]][[j]]$auc_7_con[k] <- roc(outcome_con ~ p_predicted_7_con, data_primary_1000)$auc
  result_1000[[i]][[j]]$var_auc_7_con[k] <- var(roc(outcome_con ~ p_predicted_7_con, data_primary_1000))
  
  result_1000[[i]][[j]]$auc_7_bin[k] <- roc(outcome_bin ~ p_predicted_7_bin, data_primary_1000)$auc
  result_1000[[i]][[j]]$var_auc_7_bin[k] <- var(roc(outcome_bin ~ p_predicted_7_bin, data_primary_1000))
  
  #roc.test(roc(outcome_con ~ p_predicted_6_con, data_primary_1000), roc(outcome_con ~ p_predicted_7_con, data_primary_1000))
  
  result_1000[[i]][[j]]$cov_auc_con[k] <- result_1000[[i]][[j]]$var_auc_7_con[k] + result_1000[[i]][[j]]$var_auc_6_con[k] - 2*cov(roc(outcome_con ~ p_predicted_6_con, data_primary_1000), roc(outcome_con ~ p_predicted_7_con, data_primary_1000))
  result_1000[[i]][[j]]$cov_auc_bin[k] <- result_1000[[i]][[j]]$var_auc_7_bin[k] + result_1000[[i]][[j]]$var_auc_6_bin[k] - 2*cov(roc(outcome_bin ~ p_predicted_6_bin, data_primary_1000), roc(outcome_bin ~ p_predicted_7_bin, data_primary_1000))
  
  
}}}



# meta analysis 



# generate some functions

f.logit <- function(c) {
  logit_c <- ifelse(c==1, 23, log(c/(1-c)))
  return(logit_c)
}



f.logit_var <- function(c,v) {
  logit_var <- ifelse(c==1, 10^10, v/((c*(1-c))^2))
  logit_var_floor <- ifelse(logit_var<10^(-10), 10^(-10), logit_var)
  return(logit_var_floor)
}


f.logit_rev <- function(logit) {
  c <- exp(logit)/(1+exp(logit))
  return(c)
}



# delta
con_analysis_200 <- vector(length(beta_x7_all),mode="list")
con_analysis_500 <- vector(length(beta_x7_all),mode="list")
con_analysis_1000 <- vector(length(beta_x7_all),mode="list")

# AUC1 vs AUC2
con_analysis2_200 <- vector(length(beta_x7_all),mode="list")
con_analysis2_500 <- vector(length(beta_x7_all),mode="list")
con_analysis2_1000 <- vector(length(beta_x7_all),mode="list")

# AUC1 vs delta
con_analysis3_200 <- vector(length(beta_x7_all),mode="list")
con_analysis3_500 <- vector(length(beta_x7_all),mode="list")
con_analysis3_1000 <- vector(length(beta_x7_all),mode="list")

# delta
bin_analysis_200 <- vector(length(beta_x7_all),mode="list")
bin_analysis_500 <- vector(length(beta_x7_all),mode="list")
bin_analysis_1000 <- vector(length(beta_x7_all),mode="list")


# AUC1 vs AUC2
bin_analysis2_200 <- vector(length(beta_x7_all),mode="list")
bin_analysis2_500 <- vector(length(beta_x7_all),mode="list")
bin_analysis2_1000 <- vector(length(beta_x7_all),mode="list")

# AUC1 vs delta
bin_analysis3_200 <- vector(length(beta_x7_all),mode="list")
bin_analysis3_500 <- vector(length(beta_x7_all),mode="list")
bin_analysis3_1000 <- vector(length(beta_x7_all),mode="list")

for (i in 1:n_simulation){
  #i <- 7

for (j in 1:(length(beta_x7_all))){
  
  # define effect of X7
 # j <- 4

  for (n in (length(n_study)):1){
    
    
     # meta analysis of delta AUC
    
    # continuous
    
     # sample size is 200
     con_fit_delta_200 <- uvmeta(r=(result_200[[i]][[j]]$auc_7_con-result_200[[i]][[j]]$auc_6_con)[1:n_study[n]], r.vi=ifelse(result_200[[i]][[j]]$cov_auc_con[1:n_study[n]]>rep(10^(-10),n_study[n]),result_200[[i]][[j]]$cov_auc_con[1:n_study[n]],rep(10^(-10),n_study[n])),method="REML",control=list(stepadj=0.5))
     con_analysis_200[[j]][[n]]$Est[i] <- con_fit_delta_200$est
     con_analysis_200[[j]][[n]]$ciL[i] <- con_fit_delta_200$ci.lb
     con_analysis_200[[j]][[n]]$ciU[i] <- con_fit_delta_200$ci.ub
     
     # sample size is 500
     con_fit_delta_500 <- uvmeta(r=(result_500[[i]][[j]]$auc_7_con-result_500[[i]][[j]]$auc_6_con)[1:n_study[n]], r.vi=ifelse(result_500[[i]][[j]]$cov_auc_con[1:n_study[n]]>rep(10^(-10),n_study[n]),result_500[[i]][[j]]$cov_auc_con[1:n_study[n]],rep(10^(-10),n_study[n])),method="REML",control=list(stepadj=0.5))
     con_analysis_500[[j]][[n]]$Est[i] <- con_fit_delta_500$est
     con_analysis_500[[j]][[n]]$ciL[i] <- con_fit_delta_500$ci.lb
     con_analysis_500[[j]][[n]]$ciU[i] <- con_fit_delta_500$ci.ub
     
     # sample size is 1000
     con_fit_delta_1000 <- uvmeta(r=(result_1000[[i]][[j]]$auc_7_con-result_1000[[i]][[j]]$auc_6_con)[1:n_study[n]], r.vi=ifelse(result_1000[[i]][[j]]$cov_auc_con[1:n_study[n]]>rep(10^(-10),n_study[n]),result_1000[[i]][[j]]$cov_auc_con[1:n_study[n]],rep(10^(-10),n_study[n])),method="REML",control=list(stepadj=0.5))
     con_analysis_1000[[j]][[n]]$Est[i] <- con_fit_delta_1000$est
     con_analysis_1000[[j]][[n]]$ciL[i] <- con_fit_delta_1000$ci.lb
     con_analysis_1000[[j]][[n]]$ciU[i] <- con_fit_delta_1000$ci.ub
     
     
     
     # meta analysis of AUC1 and AUC2
     
     # sample size is 200
     con_data_12_200 <- data.frame(f.logit(c=result_200[[i]][[j]]$auc_6_con),f.logit_var(c=result_200[[i]][[j]]$auc_6_con,v=result_200[[i]][[j]]$var_auc_6_con),f.logit(c=result_200[[i]][[j]]$auc_7_con),f.logit_var(c=result_200[[i]][[j]]$auc_7_con,v=result_200[[i]][[j]]$var_auc_7_con))
     
     # con_data_12_200 <- data.frame(result_200[[i]][[j]]$auc_6_con,result_200[[i]][[j]]$var_auc_6_con,result_200[[i]][[j]]$auc_7_con,result_200[[i]][[j]]$var_auc_7_con)
     names(con_data_12_200) <- c("Y1", "vars1","Y2", "vars2")
     
     con_fit_12_200 <- riley(con_data_12_200[1:n_study[n],])
     con_analysis2_200[[j]][[n]]$beta1[i] <- con_fit_12_200$coefficients[1]
     con_analysis2_200[[j]][[n]]$beta2[i] <- con_fit_12_200$coefficients[2]
     con_analysis2_200[[j]][[n]]$psi1[i] <- con_fit_12_200$coefficients[3]
     con_analysis2_200[[j]][[n]]$psi2[i] <- con_fit_12_200$coefficients[4]
     con_analysis2_200[[j]][[n]]$rhoT[i] <- con_fit_12_200$coefficients[5]
     
     # sample size is 500
     con_data_12_500 <- data.frame(f.logit(c=result_500[[i]][[j]]$auc_6_con),f.logit_var(c=result_500[[i]][[j]]$auc_6_con,v=result_500[[i]][[j]]$var_auc_6_con),f.logit(c=result_500[[i]][[j]]$auc_7_con),f.logit_var(c=result_500[[i]][[j]]$auc_7_con,v=result_500[[i]][[j]]$var_auc_7_con))
     
     # con_data_12_500 <- data.frame(result_500[[i]][[j]]$auc_6_con,result_500[[i]][[j]]$var_auc_6_con,result_500[[i]][[j]]$auc_7_con,result_500[[i]][[j]]$var_auc_7_con)
     names(con_data_12_500) <- c("Y1", "vars1","Y2", "vars2")
     
     con_fit_12_500 <- riley(con_data_12_500[1:n_study[n],])
     con_analysis2_500[[j]][[n]]$beta1[i] <- con_fit_12_500$coefficients[1]
     con_analysis2_500[[j]][[n]]$beta2[i] <- con_fit_12_500$coefficients[2]
     con_analysis2_500[[j]][[n]]$psi1[i] <- con_fit_12_500$coefficients[3]
     con_analysis2_500[[j]][[n]]$psi2[i] <- con_fit_12_500$coefficients[4]
     con_analysis2_500[[j]][[n]]$rhoT[i] <- con_fit_12_500$coefficients[5]
     
     
     # sample size is 1000
     con_data_12_1000 <- data.frame(f.logit(c=result_1000[[i]][[j]]$auc_6_con),f.logit_var(c=result_1000[[i]][[j]]$auc_6_con,v=result_1000[[i]][[j]]$var_auc_6_con),f.logit(c=result_1000[[i]][[j]]$auc_7_con),f.logit_var(c=result_1000[[i]][[j]]$auc_7_con,v=result_1000[[i]][[j]]$var_auc_7_con))
     
     #con_data_12_1000 <- data.frame(result_1000[[i]][[j]]$auc_6_con,result_1000[[i]][[j]]$var_auc_6_con,result_1000[[i]][[j]]$auc_7_con,result_1000[[i]][[j]]$var_auc_7_con)
     names(con_data_12_1000) <- c("Y1", "vars1","Y2", "vars2")
     
     con_fit_12_1000 <- riley(con_data_12_1000[1:n_study[n],])
     con_analysis2_1000[[j]][[n]]$beta1[i] <- con_fit_12_1000$coefficients[1]
     con_analysis2_1000[[j]][[n]]$beta2[i] <- con_fit_12_1000$coefficients[2]
     con_analysis2_1000[[j]][[n]]$psi1[i] <- con_fit_12_1000$coefficients[3]
     con_analysis2_1000[[j]][[n]]$psi2[i] <- con_fit_12_1000$coefficients[4]
     con_analysis2_1000[[j]][[n]]$rhoT[i] <- con_fit_12_1000$coefficients[5]
     
     
     
     # meta analysis of AUC1 and delta 
     
     # sample size is 200
     con_data_1d_200 <- data.frame(f.logit(c=result_200[[i]][[j]]$auc_6_con),f.logit_var(c=result_200[[i]][[j]]$auc_6_con,v=result_200[[i]][[j]]$var_auc_6_con),(result_200[[i]][[j]]$auc_7_con-result_200[[i]][[j]]$auc_6_con),result_200[[i]][[j]]$cov_auc_con)
     
     #con_data_1d_200 <- data.frame(result_200[[i]][[j]]$auc_6_con,result_200[[i]][[j]]$var_auc_6_con,(result_200[[i]][[j]]$auc_7_con-result_200[[i]][[j]]$auc_6_con),result_200[[i]][[j]]$cov_auc_con)
     names(con_data_1d_200) <- c("Y1", "vars1","Y2", "vars2")
     con_data_1d_200$vars2 <- ifelse(con_data_1d_200$vars2>rep(10^(-10),k_max),con_data_1d_200$vars2,rep(10^(-10),k_max))
     
     con_fit_1d_200 <- riley(con_data_1d_200[1:n_study[n],])
     con_analysis3_200[[j]][[n]]$beta1[i] <- con_fit_1d_200$coefficients[1]
     con_analysis3_200[[j]][[n]]$beta2[i] <- con_fit_1d_200$coefficients[2]
     con_analysis3_200[[j]][[n]]$psi1[i] <- con_fit_1d_200$coefficients[3]
     con_analysis3_200[[j]][[n]]$psi2[i] <- con_fit_1d_200$coefficients[4]
     con_analysis3_200[[j]][[n]]$rhoT[i] <- con_fit_1d_200$coefficients[5]
     
     
     # sample size is 500
     con_data_1d_500 <- data.frame(f.logit(c=result_500[[i]][[j]]$auc_6_con),f.logit_var(c=result_500[[i]][[j]]$auc_6_con,v=result_500[[i]][[j]]$var_auc_6_con),(result_500[[i]][[j]]$auc_7_con-result_500[[i]][[j]]$auc_6_con),result_500[[i]][[j]]$cov_auc_con)
     
    # con_data_1d_500 <- data.frame(result_500[[i]][[j]]$auc_6_con,result_500[[i]][[j]]$var_auc_6_con,(result_500[[i]][[j]]$auc_7_con-result_500[[i]][[j]]$auc_6_con),result_500[[i]][[j]]$cov_auc_con)
     names(con_data_1d_500) <- c("Y1", "vars1","Y2", "vars2")
     con_data_1d_500$vars2 <- ifelse(con_data_1d_500$vars2>rep(10^(-10),k_max),con_data_1d_500$vars2,rep(10^(-10),k_max))
     
     
     con_fit_1d_500 <- riley(con_data_1d_500[1:n_study[n],])
     con_analysis3_500[[j]][[n]]$beta1[i] <- con_fit_1d_500$coefficients[1]
     con_analysis3_500[[j]][[n]]$beta2[i] <- con_fit_1d_500$coefficients[2]
     con_analysis3_500[[j]][[n]]$psi1[i] <- con_fit_1d_500$coefficients[3]
     con_analysis3_500[[j]][[n]]$psi2[i] <- con_fit_1d_500$coefficients[4]
     con_analysis3_500[[j]][[n]]$rhoT[i] <- con_fit_1d_500$coefficients[5]
     
     # sample size is 1000
     con_data_1d_1000 <- data.frame(f.logit(c=result_1000[[i]][[j]]$auc_6_con),f.logit_var(c=result_1000[[i]][[j]]$auc_6_con,v=result_1000[[i]][[j]]$var_auc_6_con),(result_1000[[i]][[j]]$auc_7_con-result_1000[[i]][[j]]$auc_6_con),result_1000[[i]][[j]]$cov_auc_con)
     
     #con_data_1d_1000 <- data.frame(result_1000[[i]][[j]]$auc_6_con,result_1000[[i]][[j]]$var_auc_6_con,(result_1000[[i]][[j]]$auc_7_con-result_1000[[i]][[j]]$auc_6_con),result_1000[[i]][[j]]$cov_auc_con)
     names(con_data_1d_1000) <- c("Y1", "vars1","Y2", "vars2")
     con_data_1d_1000$vars2 <- ifelse(con_data_1d_1000$vars2>rep(10^(-10),k_max),con_data_1d_1000$vars2,rep(10^(-10),k_max))
     
     
     con_fit_1d_1000 <- riley(con_data_1d_1000[1:n_study[n],])
     con_analysis3_1000[[j]][[n]]$beta1[i] <- con_fit_1d_1000$coefficients[1]
     con_analysis3_1000[[j]][[n]]$beta2[i] <- con_fit_1d_1000$coefficients[2]
     con_analysis3_1000[[j]][[n]]$psi1[i] <- con_fit_1d_1000$coefficients[3]
     con_analysis3_1000[[j]][[n]]$psi2[i] <- con_fit_1d_1000$coefficients[4]
     con_analysis3_1000[[j]][[n]]$rhoT[i] <- con_fit_1d_1000$coefficients[5]
  
     
     # binary
     
     # sample size is 200
     bin_fit_delta_200 <- uvmeta(r=(result_200[[i]][[j]]$auc_7_bin-result_200[[i]][[j]]$auc_6_bin)[1:n_study[n]], r.vi=ifelse(result_200[[i]][[j]]$cov_auc_bin[1:n_study[n]]>rep(10^(-10),n_study[n]),result_200[[i]][[j]]$cov_auc_bin[1:n_study[n]],rep(10^(-10),n_study[n])),method="REML",control=list(stepadj=0.5))
     bin_analysis_200[[j]][[n]]$Est[i] <- bin_fit_delta_200$est
     bin_analysis_200[[j]][[n]]$ciL[i] <- bin_fit_delta_200$ci.lb
     bin_analysis_200[[j]][[n]]$ciU[i] <- bin_fit_delta_200$ci.ub
     
     # sample size is 500
     bin_fit_delta_500 <- uvmeta(r=(result_500[[i]][[j]]$auc_7_bin-result_500[[i]][[j]]$auc_6_bin)[1:n_study[n]], r.vi=ifelse(result_500[[i]][[j]]$cov_auc_bin[1:n_study[n]]>rep(10^(-10),n_study[n]),result_500[[i]][[j]]$cov_auc_bin[1:n_study[n]],rep(10^(-10),n_study[n])),method="REML",control=list(stepadj=0.5))
     bin_analysis_500[[j]][[n]]$Est[i] <- bin_fit_delta_500$est
     bin_analysis_500[[j]][[n]]$ciL[i] <- bin_fit_delta_500$ci.lb
     bin_analysis_500[[j]][[n]]$ciU[i] <- bin_fit_delta_500$ci.ub
     
     # sample size is 1000
     bin_fit_delta_1000 <- uvmeta(r=(result_1000[[i]][[j]]$auc_7_bin-result_1000[[i]][[j]]$auc_6_bin)[1:n_study[n]], r.vi=ifelse(result_1000[[i]][[j]]$cov_auc_bin[1:n_study[n]]>rep(10^(-10),n_study[n]),result_1000[[i]][[j]]$cov_auc_bin[1:n_study[n]],rep(10^(-10),n_study[n])),method="REML",control=list(stepadj=0.5))
     bin_analysis_1000[[j]][[n]]$Est[i] <- bin_fit_delta_1000$est
     bin_analysis_1000[[j]][[n]]$ciL[i] <- bin_fit_delta_1000$ci.lb
     bin_analysis_1000[[j]][[n]]$ciU[i] <- bin_fit_delta_1000$ci.ub
     
     
     
     # meta analysis of AUC1 and AUC2
     
     # sample size is 200
     bin_data_12_200 <- data.frame(f.logit(c=result_200[[i]][[j]]$auc_6_bin),f.logit_var(c=result_200[[i]][[j]]$auc_6_bin,v=result_200[[i]][[j]]$var_auc_6_bin),f.logit(c=result_200[[i]][[j]]$auc_7_bin),f.logit_var(c=result_200[[i]][[j]]$auc_7_bin,v=result_200[[i]][[j]]$var_auc_7_bin))
     
     #bin_data_12_200 <- data.frame(result_200[[i]][[j]]$auc_6_bin,result_200[[i]][[j]]$var_auc_6_bin,result_200[[i]][[j]]$auc_7_bin,result_200[[i]][[j]]$var_auc_7_bin)
     names(bin_data_12_200) <- c("Y1", "vars1","Y2", "vars2")
     
     bin_fit_12_200 <- riley(bin_data_12_200[1:n_study[n],])
     bin_analysis2_200[[j]][[n]]$beta1[i] <- bin_fit_12_200$coefficients[1]
     bin_analysis2_200[[j]][[n]]$beta2[i] <- bin_fit_12_200$coefficients[2]
     bin_analysis2_200[[j]][[n]]$psi1[i] <- bin_fit_12_200$coefficients[3]
     bin_analysis2_200[[j]][[n]]$psi2[i] <- bin_fit_12_200$coefficients[4]
     bin_analysis2_200[[j]][[n]]$rhoT[i] <- bin_fit_12_200$coefficients[5]
     
     # sample size is 500
     bin_data_12_500 <- data.frame(f.logit(c=result_500[[i]][[j]]$auc_6_bin),f.logit_var(c=result_500[[i]][[j]]$auc_6_bin,v=result_500[[i]][[j]]$var_auc_6_bin),f.logit(c=result_500[[i]][[j]]$auc_7_bin),f.logit_var(c=result_500[[i]][[j]]$auc_7_bin,v=result_500[[i]][[j]]$var_auc_7_bin))
     
     #bin_data_12_500 <- data.frame(result_500[[i]][[j]]$auc_6_bin,result_500[[i]][[j]]$var_auc_6_bin,result_500[[i]][[j]]$auc_7_bin,result_500[[i]][[j]]$var_auc_7_bin)
     names(bin_data_12_500) <- c("Y1", "vars1","Y2", "vars2")
     
     bin_fit_12_500 <- riley(bin_data_12_500[1:n_study[n],])
     bin_analysis2_500[[j]][[n]]$beta1[i] <- bin_fit_12_500$coefficients[1]
     bin_analysis2_500[[j]][[n]]$beta2[i] <- bin_fit_12_500$coefficients[2]
     bin_analysis2_500[[j]][[n]]$psi1[i] <- bin_fit_12_500$coefficients[3]
     bin_analysis2_500[[j]][[n]]$psi2[i] <- bin_fit_12_500$coefficients[4]
     bin_analysis2_500[[j]][[n]]$rhoT[i] <- bin_fit_12_500$coefficients[5]
     
     
     # sample size is 1000
     bin_data_12_1000 <- data.frame(f.logit(c=result_1000[[i]][[j]]$auc_6_bin),f.logit_var(c=result_1000[[i]][[j]]$auc_6_bin,v=result_1000[[i]][[j]]$var_auc_6_bin),f.logit(c=result_1000[[i]][[j]]$auc_7_bin),f.logit_var(c=result_1000[[i]][[j]]$auc_7_bin,v=result_1000[[i]][[j]]$var_auc_7_bin))
     
     #bin_data_12_1000 <- data.frame(result_1000[[i]][[j]]$auc_6_bin,result_1000[[i]][[j]]$var_auc_6_bin,result_1000[[i]][[j]]$auc_7_bin,result_1000[[i]][[j]]$var_auc_7_bin)
     names(bin_data_12_1000) <- c("Y1", "vars1","Y2", "vars2")
     
     bin_fit_12_1000 <- riley(bin_data_12_1000[1:n_study[n],])
     bin_analysis2_1000[[j]][[n]]$beta1[i] <- bin_fit_12_1000$coefficients[1]
     bin_analysis2_1000[[j]][[n]]$beta2[i] <- bin_fit_12_1000$coefficients[2]
     bin_analysis2_1000[[j]][[n]]$psi1[i] <- bin_fit_12_1000$coefficients[3]
     bin_analysis2_1000[[j]][[n]]$psi2[i] <- bin_fit_12_1000$coefficients[4]
     bin_analysis2_1000[[j]][[n]]$rhoT[i] <- bin_fit_12_1000$coefficients[5]
     
     
     
     # meta analysis of AUC1 and delta 
     
     # sample size is 200
     bin_data_1d_200 <- data.frame(f.logit(c=result_200[[i]][[j]]$auc_6_bin),f.logit_var(c=result_200[[i]][[j]]$auc_6_bin,v=result_200[[i]][[j]]$var_auc_6_bin),(result_200[[i]][[j]]$auc_7_bin-result_200[[i]][[j]]$auc_6_bin),result_200[[i]][[j]]$cov_auc_bin)
     
     #bin_data_1d_200 <- data.frame(result_200[[i]][[j]]$auc_6_bin,result_200[[i]][[j]]$var_auc_6_bin,(result_200[[i]][[j]]$auc_7_bin-result_200[[i]][[j]]$auc_6_bin),result_200[[i]][[j]]$cov_auc_bin)
     names(bin_data_1d_200) <- c("Y1", "vars1","Y2", "vars2")
     bin_data_1d_200$vars2 <- ifelse(bin_data_1d_200$vars2>rep(10^(-10),k_max),bin_data_1d_200$vars2,rep(10^(-10),k_max))
     
     
     bin_fit_1d_200 <- riley(bin_data_1d_200[1:n_study[n],])
     bin_analysis3_200[[j]][[n]]$beta1[i] <- bin_fit_1d_200$coefficients[1]
     bin_analysis3_200[[j]][[n]]$beta2[i] <- bin_fit_1d_200$coefficients[2]
     bin_analysis3_200[[j]][[n]]$psi1[i] <- bin_fit_1d_200$coefficients[3]
     bin_analysis3_200[[j]][[n]]$psi2[i] <- bin_fit_1d_200$coefficients[4]
     bin_analysis3_200[[j]][[n]]$rhoT[i] <- bin_fit_1d_200$coefficients[5]
     
     
     # sample size is 500
     bin_data_1d_500 <- data.frame(f.logit(c=result_500[[i]][[j]]$auc_6_bin),f.logit_var(c=result_500[[i]][[j]]$auc_6_bin,v=result_500[[i]][[j]]$var_auc_6_bin),(result_500[[i]][[j]]$auc_7_bin-result_500[[i]][[j]]$auc_6_bin),result_500[[i]][[j]]$cov_auc_bin)
     
     #bin_data_1d_500 <- data.frame(result_500[[i]][[j]]$auc_6_bin,result_500[[i]][[j]]$var_auc_6_bin,(result_500[[i]][[j]]$auc_7_bin-result_500[[i]][[j]]$auc_6_bin),result_500[[i]][[j]]$cov_auc_bin)
     names(bin_data_1d_500) <- c("Y1", "vars1","Y2", "vars2")
     bin_data_1d_500$vars2 <- ifelse(bin_data_1d_500$vars2>rep(10^(-10),k_max),bin_data_1d_500$vars2,rep(10^(-10),k_max))
     
     
     bin_fit_1d_500 <- riley(bin_data_1d_500[1:n_study[n],])
     bin_analysis3_500[[j]][[n]]$beta1[i] <- bin_fit_1d_500$coefficients[1]
     bin_analysis3_500[[j]][[n]]$beta2[i] <- bin_fit_1d_500$coefficients[2]
     bin_analysis3_500[[j]][[n]]$psi1[i] <- bin_fit_1d_500$coefficients[3]
     bin_analysis3_500[[j]][[n]]$psi2[i] <- bin_fit_1d_500$coefficients[4]
     bin_analysis3_500[[j]][[n]]$rhoT[i] <- bin_fit_1d_500$coefficients[5]
     
     # sample size is 1000
     bin_data_1d_1000 <- data.frame(f.logit(c=result_1000[[i]][[j]]$auc_6_bin),f.logit_var(c=result_1000[[i]][[j]]$auc_6_bin,v=result_1000[[i]][[j]]$var_auc_6_bin),(result_1000[[i]][[j]]$auc_7_bin-result_1000[[i]][[j]]$auc_6_bin),result_1000[[i]][[j]]$cov_auc_bin)
     
     #bin_data_1d_1000 <- data.frame(result_1000[[i]][[j]]$auc_6_bin,result_1000[[i]][[j]]$var_auc_6_bin,(result_1000[[i]][[j]]$auc_7_bin-result_1000[[i]][[j]]$auc_6_bin),result_1000[[i]][[j]]$cov_auc_bin)
     names(bin_data_1d_1000) <- c("Y1", "vars1","Y2", "vars2")
     bin_data_1d_1000$vars2 <- ifelse(bin_data_1d_1000$vars2>rep(10^(-10),k_max),bin_data_1d_1000$vars2,rep(10^(-10),k_max))
     
     
     bin_fit_1d_1000 <- riley(bin_data_1d_1000[1:n_study[n],])
     bin_analysis3_1000[[j]][[n]]$beta1[i] <- bin_fit_1d_1000$coefficients[1]
     bin_analysis3_1000[[j]][[n]]$beta2[i] <- bin_fit_1d_1000$coefficients[2]
     bin_analysis3_1000[[j]][[n]]$psi1[i] <- bin_fit_1d_1000$coefficients[3]
     bin_analysis3_1000[[j]][[n]]$psi2[i] <- bin_fit_1d_1000$coefficients[4]
     bin_analysis3_1000[[j]][[n]]$rhoT[i] <- bin_fit_1d_1000$coefficients[5]
  }}}




##### summarize results

# continuous marker

# sample size = 200
## strategy 1 delta

# sample size = 200, beta = 0

# con_est_0_0_200_5_delta
# correlation 0; beta 0; sample size 200; number of study 5; strategy 1 delta

con_est_0_0_200_5_delta <- data.frame(con_analysis_200[[1]][[1]]$Est,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(con_est_0_0_200_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_200_10_delta <- data.frame(con_analysis_200[[1]][[2]]$Est,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(con_est_0_0_200_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_200_20_delta <- data.frame(con_analysis_200[[1]][[3]]$Est,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(con_est_0_0_200_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.4

con_est_0_04_200_5_delta <- data.frame(con_analysis_200[[2]][[1]]$Est,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_200_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_200_10_delta <- data.frame(con_analysis_200[[2]][[2]]$Est,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_200_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_200_20_delta <- data.frame(con_analysis_200[[2]][[3]]$Est,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_200_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.7

con_est_0_07_200_5_delta <- data.frame(con_analysis_200[[3]][[1]]$Est,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_200_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_200_10_delta <- data.frame(con_analysis_200[[3]][[2]]$Est,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_200_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_200_20_delta <- data.frame(con_analysis_200[[3]][[3]]$Est,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_200_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 200, beta = 0.9

con_est_0_09_200_5_delta <- data.frame(con_analysis_200[[4]][[1]]$Est,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_200_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_200_10_delta <- data.frame(con_analysis_200[[4]][[2]]$Est,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_200_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_200_20_delta <- data.frame(con_analysis_200[[4]][[3]]$Est,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_200_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")




## strategy 2 AUC12

# sample size = 200, beta = 0

# con_est_0_0_200_5_AUC12
# correlation 0; beta 0; sample size 200; number of study 5; strategy 2 AUC12

con_est_0_0_200_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[1]][[1]]$beta2)-f.logit_rev(con_analysis2_200[[1]][[1]]$beta1),rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(con_est_0_0_200_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_200_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[1]][[2]]$beta2)-f.logit_rev(con_analysis2_200[[1]][[2]]$beta1),rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(con_est_0_0_200_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_200_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[1]][[3]]$beta2)-f.logit_rev(con_analysis2_200[[1]][[3]]$beta1),rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(con_est_0_0_200_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.4

con_est_0_04_200_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[2]][[1]]$beta2)-f.logit_rev(con_analysis2_200[[2]][[1]]$beta1),rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_200_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_200_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[2]][[2]]$beta2)-f.logit_rev(con_analysis2_200[[2]][[2]]$beta1),rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_200_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_200_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[2]][[3]]$beta2)-f.logit_rev(con_analysis2_200[[2]][[3]]$beta1),rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_200_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.7

con_est_0_07_200_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[3]][[1]]$beta2)-f.logit_rev(con_analysis2_200[[3]][[1]]$beta1),rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_200_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_200_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[3]][[2]]$beta2)-f.logit_rev(con_analysis2_200[[3]][[2]]$beta1),rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_200_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_200_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[3]][[3]]$beta2)-f.logit_rev(con_analysis2_200[[3]][[3]]$beta1),rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_200_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 200, beta = 0.9

con_est_0_09_200_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[4]][[1]]$beta2)-f.logit_rev(con_analysis2_200[[4]][[1]]$beta1),rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_200_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_200_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[4]][[2]]$beta2)-f.logit_rev(con_analysis2_200[[4]][[2]]$beta1),rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_200_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_200_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_200[[4]][[3]]$beta2)-f.logit_rev(con_analysis2_200[[4]][[3]]$beta1),rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_200_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


## strategy 3 AUC1d

# sample size = 200, beta = 0

# con_est_0_0_200_5_AUC1d
# correlation 0; beta 0; sample size 200; number of study 5; strategy 3 AUC1d

con_est_0_0_200_5_AUC1d <- data.frame(con_analysis3_200[[1]][[1]]$beta2,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(con_est_0_0_200_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_200_10_AUC1d <- data.frame(con_analysis3_200[[1]][[2]]$beta2,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(con_est_0_0_200_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_200_20_AUC1d <- data.frame(con_analysis3_200[[1]][[3]]$beta2,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(con_est_0_0_200_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.4

con_est_0_04_200_5_AUC1d <- data.frame(con_analysis3_200[[2]][[1]]$beta2,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_200_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_200_10_AUC1d <- data.frame(con_analysis3_200[[2]][[2]]$beta2,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_200_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_200_20_AUC1d <- data.frame(con_analysis3_200[[2]][[3]]$beta2,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_200_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.7

con_est_0_07_200_5_AUC1d <- data.frame(con_analysis3_200[[3]][[1]]$beta2,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_200_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_200_10_AUC1d <- data.frame(con_analysis3_200[[3]][[2]]$beta2,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_200_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_200_20_AUC1d <- data.frame(con_analysis3_200[[3]][[3]]$beta2,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_200_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 200, beta = 0.9

con_est_0_09_200_5_AUC1d <- data.frame(con_analysis3_200[[4]][[1]]$beta2,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_200_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_200_10_AUC1d <- data.frame(con_analysis3_200[[4]][[2]]$beta2,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_200_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_200_20_AUC1d <- data.frame(con_analysis3_200[[4]][[3]]$beta2,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_200_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")



# sample size = 500

## strategy 1 delta

# sample size = 500, beta = 0

# con_est_0_0_500_5_delta
# correlation 0; beta 0; sample size 500; number of study 5; strategy 1 delta

con_est_0_0_500_5_delta <- data.frame(con_analysis_500[[1]][[1]]$Est,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(con_est_0_0_500_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_500_10_delta <- data.frame(con_analysis_500[[1]][[2]]$Est,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(con_est_0_0_500_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_500_20_delta <- data.frame(con_analysis_500[[1]][[3]]$Est,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(con_est_0_0_500_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.4

con_est_0_04_500_5_delta <- data.frame(con_analysis_500[[2]][[1]]$Est,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_500_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_500_10_delta <- data.frame(con_analysis_500[[2]][[2]]$Est,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_500_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_500_20_delta <- data.frame(con_analysis_500[[2]][[3]]$Est,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_500_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.7

con_est_0_07_500_5_delta <- data.frame(con_analysis_500[[3]][[1]]$Est,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_500_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_500_10_delta <- data.frame(con_analysis_500[[3]][[2]]$Est,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_500_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_500_20_delta <- data.frame(con_analysis_500[[3]][[3]]$Est,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_500_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 500, beta = 0.9

con_est_0_09_500_5_delta <- data.frame(con_analysis_500[[4]][[1]]$Est,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_500_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_500_10_delta <- data.frame(con_analysis_500[[4]][[2]]$Est,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_500_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_500_20_delta <- data.frame(con_analysis_500[[4]][[3]]$Est,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_500_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")




## strategy 2 AUC12

# sample size = 500, beta = 0

# con_est_0_0_500_5_AUC12
# correlation 0; beta 0; sample size 500; number of study 5; strategy 2 AUC12

con_est_0_0_500_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[1]][[1]]$beta2)-f.logit_rev(con_analysis2_500[[1]][[1]]$beta1),rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(con_est_0_0_500_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_500_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[1]][[2]]$beta2)-f.logit_rev(con_analysis2_500[[1]][[2]]$beta1),rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(con_est_0_0_500_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_500_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[1]][[3]]$beta2)-f.logit_rev(con_analysis2_500[[1]][[3]]$beta1),rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(con_est_0_0_500_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.4

con_est_0_04_500_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[2]][[1]]$beta2)-f.logit_rev(con_analysis2_500[[2]][[1]]$beta1),rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_500_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_500_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[2]][[2]]$beta2)-f.logit_rev(con_analysis2_500[[2]][[2]]$beta1),rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_500_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_500_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[2]][[3]]$beta2)-f.logit_rev(con_analysis2_500[[2]][[3]]$beta1),rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_500_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.7

con_est_0_07_500_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[3]][[1]]$beta2)-f.logit_rev(con_analysis2_500[[3]][[1]]$beta1),rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_500_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_500_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[3]][[2]]$beta2)-f.logit_rev(con_analysis2_500[[3]][[2]]$beta1),rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_500_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_500_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[3]][[3]]$beta2)-f.logit_rev(con_analysis2_500[[3]][[3]]$beta1),rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_500_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 500, beta = 0.9

con_est_0_09_500_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[4]][[1]]$beta2)-f.logit_rev(con_analysis2_500[[4]][[1]]$beta1),rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_500_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_500_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[4]][[2]]$beta2)-f.logit_rev(con_analysis2_500[[4]][[2]]$beta1),rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_500_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_500_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_500[[4]][[3]]$beta2)-f.logit_rev(con_analysis2_500[[4]][[3]]$beta1),rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_500_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


## strategy 3 AUC1d

# sample size = 500, beta = 0

# con_est_0_0_500_5_AUC1d
# correlation 0; beta 0; sample size 500; number of study 5; strategy 3 AUC1d

con_est_0_0_500_5_AUC1d <- data.frame(con_analysis3_500[[1]][[1]]$beta2,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(con_est_0_0_500_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_500_10_AUC1d <- data.frame(con_analysis3_500[[1]][[2]]$beta2,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(con_est_0_0_500_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_500_20_AUC1d <- data.frame(con_analysis3_500[[1]][[3]]$beta2,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(con_est_0_0_500_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.4

con_est_0_04_500_5_AUC1d <- data.frame(con_analysis3_500[[2]][[1]]$beta2,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_500_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_500_10_AUC1d <- data.frame(con_analysis3_500[[2]][[2]]$beta2,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_500_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_500_20_AUC1d <- data.frame(con_analysis3_500[[2]][[3]]$beta2,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_500_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.7

con_est_0_07_500_5_AUC1d <- data.frame(con_analysis3_500[[3]][[1]]$beta2,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_500_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_500_10_AUC1d <- data.frame(con_analysis3_500[[3]][[2]]$beta2,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_500_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_500_20_AUC1d <- data.frame(con_analysis3_500[[3]][[3]]$beta2,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_500_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 500, beta = 0.9

con_est_0_09_500_5_AUC1d <- data.frame(con_analysis3_500[[4]][[1]]$beta2,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_500_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_500_10_AUC1d <- data.frame(con_analysis3_500[[4]][[2]]$beta2,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_500_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_500_20_AUC1d <- data.frame(con_analysis3_500[[4]][[3]]$beta2,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_500_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")




# sample size = 1000

## strategy 1 delta

# sample size = 1000, beta = 0

# con_est_0_0_1000_5_delta
# correlation 0; beta 0; sample size 1000; number of study 5; strategy 1 delta

con_est_0_0_1000_5_delta <- data.frame(con_analysis_1000[[1]][[1]]$Est,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(con_est_0_0_1000_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_1000_10_delta <- data.frame(con_analysis_1000[[1]][[2]]$Est,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(con_est_0_0_1000_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_1000_20_delta <- data.frame(con_analysis_1000[[1]][[3]]$Est,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(con_est_0_0_1000_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.4

con_est_0_04_1000_5_delta <- data.frame(con_analysis_1000[[2]][[1]]$Est,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_1000_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_1000_10_delta <- data.frame(con_analysis_1000[[2]][[2]]$Est,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_1000_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_1000_20_delta <- data.frame(con_analysis_1000[[2]][[3]]$Est,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_1000_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.7

con_est_0_07_1000_5_delta <- data.frame(con_analysis_1000[[3]][[1]]$Est,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_1000_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_1000_10_delta <- data.frame(con_analysis_1000[[3]][[2]]$Est,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_1000_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_1000_20_delta <- data.frame(con_analysis_1000[[3]][[3]]$Est,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_1000_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 1000, beta = 0.9

con_est_0_09_1000_5_delta <- data.frame(con_analysis_1000[[4]][[1]]$Est,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_1000_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_1000_10_delta <- data.frame(con_analysis_1000[[4]][[2]]$Est,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_1000_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_1000_20_delta <- data.frame(con_analysis_1000[[4]][[3]]$Est,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_1000_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")




## strategy 2 AUC12

# sample size = 1000, beta = 0

# con_est_0_0_1000_5_AUC12
# correlation 0; beta 0; sample size 1000; number of study 5; strategy 2 AUC12

con_est_0_0_1000_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[1]][[1]]$beta2)-f.logit_rev(con_analysis2_1000[[1]][[1]]$beta1),rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(con_est_0_0_1000_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_1000_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[1]][[2]]$beta2)-f.logit_rev(con_analysis2_1000[[1]][[2]]$beta1),rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(con_est_0_0_1000_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_1000_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[1]][[3]]$beta2)-f.logit_rev(con_analysis2_1000[[1]][[3]]$beta1),rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(con_est_0_0_1000_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.4

con_est_0_04_1000_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[2]][[1]]$beta2)-f.logit_rev(con_analysis2_1000[[2]][[1]]$beta1),rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_1000_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_1000_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[2]][[2]]$beta2)-f.logit_rev(con_analysis2_1000[[2]][[2]]$beta1),rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_1000_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_1000_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[2]][[3]]$beta2)-f.logit_rev(con_analysis2_1000[[2]][[3]]$beta1),rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_1000_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.7

con_est_0_07_1000_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[3]][[1]]$beta2)-f.logit_rev(con_analysis2_1000[[3]][[1]]$beta1),rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_1000_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_1000_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[3]][[2]]$beta2)-f.logit_rev(con_analysis2_1000[[3]][[2]]$beta1),rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_1000_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_1000_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[3]][[3]]$beta2)-f.logit_rev(con_analysis2_1000[[3]][[3]]$beta1),rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_1000_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 1000, beta = 0.9

con_est_0_09_1000_5_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[4]][[1]]$beta2)-f.logit_rev(con_analysis2_1000[[4]][[1]]$beta1),rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_1000_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_1000_10_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[4]][[2]]$beta2)-f.logit_rev(con_analysis2_1000[[4]][[2]]$beta1),rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_1000_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_1000_20_AUC12 <- data.frame(f.logit_rev(con_analysis2_1000[[4]][[3]]$beta2)-f.logit_rev(con_analysis2_1000[[4]][[3]]$beta1),rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_1000_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


## strategy 3 AUC1d

# sample size = 1000, beta = 0

# con_est_0_0_1000_5_AUC1d
# correlation 0; beta 0; sample size 1000; number of study 5; strategy 3 AUC1d

con_est_0_0_1000_5_AUC1d <- data.frame(con_analysis3_1000[[1]][[1]]$beta2,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(con_est_0_0_1000_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_1000_10_AUC1d <- data.frame(con_analysis3_1000[[1]][[2]]$beta2,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(con_est_0_0_1000_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_0_1000_20_AUC1d <- data.frame(con_analysis3_1000[[1]][[3]]$beta2,rep((performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(con_est_0_0_1000_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.4

con_est_0_04_1000_5_AUC1d <- data.frame(con_analysis3_1000[[2]][[1]]$beta2,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_1000_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_1000_10_AUC1d <- data.frame(con_analysis3_1000[[2]][[2]]$beta2,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_1000_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_04_1000_20_AUC1d <- data.frame(con_analysis3_1000[[2]][[3]]$beta2,rep((performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(con_est_0_04_1000_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.7

con_est_0_07_1000_5_AUC1d <- data.frame(con_analysis3_1000[[3]][[1]]$beta2,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_1000_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_1000_10_AUC1d <- data.frame(con_analysis3_1000[[3]][[2]]$beta2,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_1000_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_07_1000_20_AUC1d <- data.frame(con_analysis3_1000[[3]][[3]]$beta2,rep((performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(con_est_0_07_1000_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 1000, beta = 0.9

con_est_0_09_1000_5_AUC1d <- data.frame(con_analysis3_1000[[4]][[1]]$beta2,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_1000_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_1000_10_AUC1d <- data.frame(con_analysis3_1000[[4]][[2]]$beta2,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_1000_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

con_est_0_09_1000_20_AUC1d <- data.frame(con_analysis3_1000[[4]][[3]]$beta2,rep((performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(con_est_0_09_1000_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")




# binary marker

# sample size = 200
## strategy 1 delta

# sample size = 200, beta = 0

# bin_est_0_0_200_5_delta
# correlation 0; beta 0; sample size 200; number of study 5; strategy 1 delta

bin_est_0_0_200_5_delta <- data.frame(bin_analysis_200[[1]][[1]]$Est,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_200_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_200_10_delta <- data.frame(bin_analysis_200[[1]][[2]]$Est,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_200_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_200_20_delta <- data.frame(bin_analysis_200[[1]][[3]]$Est,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_200_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.4

bin_est_0_04_200_5_delta <- data.frame(bin_analysis_200[[2]][[1]]$Est,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_200_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_200_10_delta <- data.frame(bin_analysis_200[[2]][[2]]$Est,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_200_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_200_20_delta <- data.frame(bin_analysis_200[[2]][[3]]$Est,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_200_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.7

bin_est_0_07_200_5_delta <- data.frame(bin_analysis_200[[3]][[1]]$Est,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_200_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_200_10_delta <- data.frame(bin_analysis_200[[3]][[2]]$Est,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_200_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_200_20_delta <- data.frame(bin_analysis_200[[3]][[3]]$Est,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_200_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 200, beta = 0.9

bin_est_0_09_200_5_delta <- data.frame(bin_analysis_200[[4]][[1]]$Est,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_200_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_200_10_delta <- data.frame(bin_analysis_200[[4]][[2]]$Est,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_200_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_200_20_delta <- data.frame(bin_analysis_200[[4]][[3]]$Est,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_200_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")




## strategy 2 AUC12

# sample size = 200, beta = 0

# bin_est_0_0_200_5_AUC12
# correlation 0; beta 0; sample size 200; number of study 5; strategy 2 AUC12

bin_est_0_0_200_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[1]][[1]]$beta2)-f.logit_rev(bin_analysis2_200[[1]][[1]]$beta1),rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_200_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_200_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[1]][[2]]$beta2)-f.logit_rev(bin_analysis2_200[[1]][[2]]$beta1),rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_200_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_200_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[1]][[3]]$beta2)-f.logit_rev(bin_analysis2_200[[1]][[3]]$beta1),rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_200_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.4

bin_est_0_04_200_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[2]][[1]]$beta2)-f.logit_rev(bin_analysis2_200[[2]][[1]]$beta1),rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_200_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_200_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[2]][[2]]$beta2)-f.logit_rev(bin_analysis2_200[[2]][[2]]$beta1),rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_200_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_200_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[2]][[3]]$beta2)-f.logit_rev(bin_analysis2_200[[2]][[3]]$beta1),rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_200_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.7

bin_est_0_07_200_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[3]][[1]]$beta2)-f.logit_rev(bin_analysis2_200[[3]][[1]]$beta1),rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_200_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_200_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[3]][[2]]$beta2)-f.logit_rev(bin_analysis2_200[[3]][[2]]$beta1),rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_200_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_200_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[3]][[3]]$beta2)-f.logit_rev(bin_analysis2_200[[3]][[3]]$beta1),rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_200_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 200, beta = 0.9

bin_est_0_09_200_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[4]][[1]]$beta2)-f.logit_rev(bin_analysis2_200[[4]][[1]]$beta1),rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_200_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_200_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[4]][[2]]$beta2)-f.logit_rev(bin_analysis2_200[[4]][[2]]$beta1),rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_200_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_200_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_200[[4]][[3]]$beta2)-f.logit_rev(bin_analysis2_200[[4]][[3]]$beta1),rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_200_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


## strategy 3 AUC1d

# sample size = 200, beta = 0

# bin_est_0_0_200_5_AUC1d
# correlation 0; beta 0; sample size 200; number of study 5; strategy 3 AUC1d

bin_est_0_0_200_5_AUC1d <- data.frame(bin_analysis3_200[[1]][[1]]$beta2,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_200_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_200_10_AUC1d <- data.frame(bin_analysis3_200[[1]][[2]]$beta2,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_200_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_200_20_AUC1d <- data.frame(bin_analysis3_200[[1]][[3]]$beta2,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_200_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.4

bin_est_0_04_200_5_AUC1d <- data.frame(bin_analysis3_200[[2]][[1]]$beta2,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_200_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_200_10_AUC1d <- data.frame(bin_analysis3_200[[2]][[2]]$beta2,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_200_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_200_20_AUC1d <- data.frame(bin_analysis3_200[[2]][[3]]$beta2,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_200_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 200, beta = 0.7

bin_est_0_07_200_5_AUC1d <- data.frame(bin_analysis3_200[[3]][[1]]$beta2,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_200_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_200_10_AUC1d <- data.frame(bin_analysis3_200[[3]][[2]]$beta2,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_200_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_200_20_AUC1d <- data.frame(bin_analysis3_200[[3]][[3]]$beta2,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_200_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 200, beta = 0.9

bin_est_0_09_200_5_AUC1d <- data.frame(bin_analysis3_200[[4]][[1]]$beta2,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_200_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_200_10_AUC1d <- data.frame(bin_analysis3_200[[4]][[2]]$beta2,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_200_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_200_20_AUC1d <- data.frame(bin_analysis3_200[[4]][[3]]$beta2,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(200,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_200_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")



# sample size = 500

## strategy 1 delta

# sample size = 500, beta = 0

# bin_est_0_0_500_5_delta
# correlation 0; beta 0; sample size 500; number of study 5; strategy 1 delta

bin_est_0_0_500_5_delta <- data.frame(bin_analysis_500[[1]][[1]]$Est,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_500_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_500_10_delta <- data.frame(bin_analysis_500[[1]][[2]]$Est,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_500_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_500_20_delta <- data.frame(bin_analysis_500[[1]][[3]]$Est,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_500_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.4

bin_est_0_04_500_5_delta <- data.frame(bin_analysis_500[[2]][[1]]$Est,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_500_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_500_10_delta <- data.frame(bin_analysis_500[[2]][[2]]$Est,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_500_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_500_20_delta <- data.frame(bin_analysis_500[[2]][[3]]$Est,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_500_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.7

bin_est_0_07_500_5_delta <- data.frame(bin_analysis_500[[3]][[1]]$Est,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_500_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_500_10_delta <- data.frame(bin_analysis_500[[3]][[2]]$Est,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_500_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_500_20_delta <- data.frame(bin_analysis_500[[3]][[3]]$Est,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_500_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 500, beta = 0.9

bin_est_0_09_500_5_delta <- data.frame(bin_analysis_500[[4]][[1]]$Est,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_500_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_500_10_delta <- data.frame(bin_analysis_500[[4]][[2]]$Est,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_500_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_500_20_delta <- data.frame(bin_analysis_500[[4]][[3]]$Est,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_500_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")




## strategy 2 AUC12

# sample size = 500, beta = 0

# bin_est_0_0_500_5_AUC12
# correlation 0; beta 0; sample size 500; number of study 5; strategy 2 AUC12

bin_est_0_0_500_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[1]][[1]]$beta2)-f.logit_rev(bin_analysis2_500[[1]][[1]]$beta1),rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_500_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_500_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[1]][[2]]$beta2)-f.logit_rev(bin_analysis2_500[[1]][[2]]$beta1),rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_500_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_500_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[1]][[3]]$beta2)-f.logit_rev(bin_analysis2_500[[1]][[3]]$beta1),rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_500_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.4

bin_est_0_04_500_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[2]][[1]]$beta2)-f.logit_rev(bin_analysis2_500[[2]][[1]]$beta1),rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_500_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_500_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[2]][[2]]$beta2)-f.logit_rev(bin_analysis2_500[[2]][[2]]$beta1),rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_500_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_500_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[2]][[3]]$beta2)-f.logit_rev(bin_analysis2_500[[2]][[3]]$beta1),rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_500_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.7

bin_est_0_07_500_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[3]][[1]]$beta2)-f.logit_rev(bin_analysis2_500[[3]][[1]]$beta1),rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_500_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_500_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[3]][[2]]$beta2)-f.logit_rev(bin_analysis2_500[[3]][[2]]$beta1),rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_500_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_500_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[3]][[3]]$beta2)-f.logit_rev(bin_analysis2_500[[3]][[3]]$beta1),rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_500_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 500, beta = 0.9

bin_est_0_09_500_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[4]][[1]]$beta2)-f.logit_rev(bin_analysis2_500[[4]][[1]]$beta1),rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_500_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_500_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[4]][[2]]$beta2)-f.logit_rev(bin_analysis2_500[[4]][[2]]$beta1),rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_500_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_500_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_500[[4]][[3]]$beta2)-f.logit_rev(bin_analysis2_500[[4]][[3]]$beta1),rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_500_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


## strategy 3 AUC1d

# sample size = 500, beta = 0

# bin_est_0_0_500_5_AUC1d
# correlation 0; beta 0; sample size 500; number of study 5; strategy 3 AUC1d

bin_est_0_0_500_5_AUC1d <- data.frame(bin_analysis3_500[[1]][[1]]$beta2,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_500_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_500_10_AUC1d <- data.frame(bin_analysis3_500[[1]][[2]]$beta2,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_500_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_500_20_AUC1d <- data.frame(bin_analysis3_500[[1]][[3]]$beta2,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_500_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.4

bin_est_0_04_500_5_AUC1d <- data.frame(bin_analysis3_500[[2]][[1]]$beta2,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_500_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_500_10_AUC1d <- data.frame(bin_analysis3_500[[2]][[2]]$beta2,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_500_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_500_20_AUC1d <- data.frame(bin_analysis3_500[[2]][[3]]$beta2,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_500_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 500, beta = 0.7

bin_est_0_07_500_5_AUC1d <- data.frame(bin_analysis3_500[[3]][[1]]$beta2,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_500_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_500_10_AUC1d <- data.frame(bin_analysis3_500[[3]][[2]]$beta2,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_500_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_500_20_AUC1d <- data.frame(bin_analysis3_500[[3]][[3]]$beta2,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_500_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 500, beta = 0.9

bin_est_0_09_500_5_AUC1d <- data.frame(bin_analysis3_500[[4]][[1]]$beta2,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_500_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_500_10_AUC1d <- data.frame(bin_analysis3_500[[4]][[2]]$beta2,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_500_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_500_20_AUC1d <- data.frame(bin_analysis3_500[[4]][[3]]$beta2,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(500,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_500_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")




# sample size = 1000

## strategy 1 delta

# sample size = 1000, beta = 0

# bin_est_0_0_1000_5_delta
# correlation 0; beta 0; sample size 1000; number of study 5; strategy 1 delta

bin_est_0_0_1000_5_delta <- data.frame(bin_analysis_1000[[1]][[1]]$Est,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_1000_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_1000_10_delta <- data.frame(bin_analysis_1000[[1]][[2]]$Est,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_1000_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_1000_20_delta <- data.frame(bin_analysis_1000[[1]][[3]]$Est,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_1000_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.4

bin_est_0_04_1000_5_delta <- data.frame(bin_analysis_1000[[2]][[1]]$Est,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_1000_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_1000_10_delta <- data.frame(bin_analysis_1000[[2]][[2]]$Est,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_1000_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_1000_20_delta <- data.frame(bin_analysis_1000[[2]][[3]]$Est,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_1000_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.7

bin_est_0_07_1000_5_delta <- data.frame(bin_analysis_1000[[3]][[1]]$Est,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_1000_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_1000_10_delta <- data.frame(bin_analysis_1000[[3]][[2]]$Est,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_1000_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_1000_20_delta <- data.frame(bin_analysis_1000[[3]][[3]]$Est,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_1000_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 1000, beta = 0.9

bin_est_0_09_1000_5_delta <- data.frame(bin_analysis_1000[[4]][[1]]$Est,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_1000_5_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_1000_10_delta <- data.frame(bin_analysis_1000[[4]][[2]]$Est,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_1000_10_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_1000_20_delta <- data.frame(bin_analysis_1000[[4]][[3]]$Est,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(1,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_1000_20_delta) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")




## strategy 2 AUC12

# sample size = 1000, beta = 0

# bin_est_0_0_1000_5_AUC12
# correlation 0; beta 0; sample size 1000; number of study 5; strategy 2 AUC12

bin_est_0_0_1000_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[1]][[1]]$beta2)-f.logit_rev(bin_analysis2_1000[[1]][[1]]$beta1),rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_1000_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_1000_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[1]][[2]]$beta2)-f.logit_rev(bin_analysis2_1000[[1]][[2]]$beta1),rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_1000_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_1000_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[1]][[3]]$beta2)-f.logit_rev(bin_analysis2_1000[[1]][[3]]$beta1),rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_1000_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.4

bin_est_0_04_1000_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[2]][[1]]$beta2)-f.logit_rev(bin_analysis2_1000[[2]][[1]]$beta1),rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_1000_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_1000_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[2]][[2]]$beta2)-f.logit_rev(bin_analysis2_1000[[2]][[2]]$beta1),rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_1000_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_1000_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[2]][[3]]$beta2)-f.logit_rev(bin_analysis2_1000[[2]][[3]]$beta1),rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_1000_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.7

bin_est_0_07_1000_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[3]][[1]]$beta2)-f.logit_rev(bin_analysis2_1000[[3]][[1]]$beta1),rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_1000_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_1000_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[3]][[2]]$beta2)-f.logit_rev(bin_analysis2_1000[[3]][[2]]$beta1),rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_1000_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_1000_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[3]][[3]]$beta2)-f.logit_rev(bin_analysis2_1000[[3]][[3]]$beta1),rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_1000_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 1000, beta = 0.9

bin_est_0_09_1000_5_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[4]][[1]]$beta2)-f.logit_rev(bin_analysis2_1000[[4]][[1]]$beta1),rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_1000_5_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_1000_10_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[4]][[2]]$beta2)-f.logit_rev(bin_analysis2_1000[[4]][[2]]$beta1),rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_1000_10_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_1000_20_AUC12 <- data.frame(f.logit_rev(bin_analysis2_1000[[4]][[3]]$beta2)-f.logit_rev(bin_analysis2_1000[[4]][[3]]$beta1),rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(2,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_1000_20_AUC12) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


## strategy 3 AUC1d

# sample size = 1000, beta = 0

# bin_est_0_0_1000_5_AUC1d
# correlation 0; beta 0; sample size 1000; number of study 5; strategy 3 AUC1d

bin_est_0_0_1000_5_AUC1d <- data.frame(bin_analysis3_1000[[1]][[1]]$beta2,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_1000_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_1000_10_AUC1d <- data.frame(bin_analysis3_1000[[1]][[2]]$beta2,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_1000_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_0_1000_20_AUC1d <- data.frame(bin_analysis3_1000[[1]][[3]]$beta2,rep((performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0,n_simulation))
names(bin_est_0_0_1000_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.4

bin_est_0_04_1000_5_AUC1d <- data.frame(bin_analysis3_1000[[2]][[1]]$beta2,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_1000_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_1000_10_AUC1d <- data.frame(bin_analysis3_1000[[2]][[2]]$beta2,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_1000_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_04_1000_20_AUC1d <- data.frame(bin_analysis3_1000[[2]][[3]]$beta2,rep((performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.4,n_simulation))
names(bin_est_0_04_1000_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

# sample size = 1000, beta = 0.7

bin_est_0_07_1000_5_AUC1d <- data.frame(bin_analysis3_1000[[3]][[1]]$beta2,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_1000_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_1000_10_AUC1d <- data.frame(bin_analysis3_1000[[3]][[2]]$beta2,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_1000_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_07_1000_20_AUC1d <- data.frame(bin_analysis3_1000[[3]][[3]]$beta2,rep((performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.7,n_simulation))
names(bin_est_0_07_1000_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


# sample size = 1000, beta = 0.9

bin_est_0_09_1000_5_AUC1d <- data.frame(bin_analysis3_1000[[4]][[1]]$beta2,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(5,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_1000_5_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_1000_10_AUC1d <- data.frame(bin_analysis3_1000[[4]][[2]]$beta2,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(10,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_1000_10_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")

bin_est_0_09_1000_20_AUC1d <- data.frame(bin_analysis3_1000[[4]][[3]]$beta2,rep((performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]),n_simulation),rep(1000,n_simulation),rep(20,n_simulation),rep(3,n_simulation),rep(0.9,n_simulation))
names(bin_est_0_09_1000_20_AUC1d) <- c("Est", "Ref","N.Sample", "N.Study","Strategy","Beta")


















# present results


# continuous marker

# beta = 0
# sample size =200

con_est_0_0_200_all <- rbind(
con_est_0_0_200_5_delta, 
con_est_0_0_200_10_delta,
con_est_0_0_200_20_delta,  
con_est_0_0_200_5_AUC12, 
con_est_0_0_200_10_AUC12, 
con_est_0_0_200_20_AUC12, 
con_est_0_0_200_5_AUC1d, 
con_est_0_0_200_10_AUC1d, 
con_est_0_0_200_20_AUC1d)


con_est_0_0_200_all$N.Study <- as.factor(con_est_0_0_200_all$N.Study)
con_est_0_0_200_all$Strategy <- as.factor(con_est_0_0_200_all$Strategy)


ggplot(con_est_0_0_200_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot() +ylim(-0.002, 0.002)+
  geom_hline(yintercept=(performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]), linetype="dashed", color = "red")

# sample size = 500

con_est_0_0_500_all <- rbind(
  con_est_0_0_500_5_delta, 
  con_est_0_0_500_10_delta,
  con_est_0_0_500_20_delta,  
  con_est_0_0_500_5_AUC12, 
  con_est_0_0_500_10_AUC12, 
  con_est_0_0_500_20_AUC12, 
  con_est_0_0_500_5_AUC1d, 
  con_est_0_0_500_10_AUC1d, 
  con_est_0_0_500_20_AUC1d)


con_est_0_0_500_all$N.Study <- as.factor(con_est_0_0_500_all$N.Study)
con_est_0_0_500_all$Strategy <- as.factor(con_est_0_0_500_all$Strategy)


ggplot(con_est_0_0_500_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(-0.002, 0.002)+
  geom_hline(yintercept=(performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]), linetype="dashed", color = "red")


# sample size = 1000

con_est_0_0_1000_all <- rbind(
  con_est_0_0_1000_5_delta, 
  con_est_0_0_1000_10_delta,
  con_est_0_0_1000_20_delta,  
  con_est_0_0_1000_5_AUC12, 
  con_est_0_0_1000_10_AUC12, 
  con_est_0_0_1000_20_AUC12, 
  con_est_0_0_1000_5_AUC1d, 
  con_est_0_0_1000_10_AUC1d, 
  con_est_0_0_1000_20_AUC1d)


con_est_0_0_1000_all$N.Study <- as.factor(con_est_0_0_1000_all$N.Study)
con_est_0_0_1000_all$Strategy <- as.factor(con_est_0_0_1000_all$Strategy)


ggplot(con_est_0_0_1000_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(-0.002, 0.002)+
  geom_hline(yintercept=(performance[[1]]$auc_7_con[1]-performance[[1]]$auc_6_con[1]), linetype="dashed", color = "red")



# beta = 0.4
# sample size = 200


con_est_0_04_200_all <- rbind(
  con_est_0_04_200_5_delta, 
  con_est_0_04_200_10_delta,
  con_est_0_04_200_20_delta,  
  con_est_0_04_200_5_AUC12, 
  con_est_0_04_200_10_AUC12, 
  con_est_0_04_200_20_AUC12, 
  con_est_0_04_200_5_AUC1d, 
  con_est_0_04_200_10_AUC1d, 
  con_est_0_04_200_20_AUC1d)


con_est_0_04_200_all$N.Study <- as.factor(con_est_0_04_200_all$N.Study)
con_est_0_04_200_all$Strategy <- as.factor(con_est_0_04_200_all$Strategy)


ggplot(con_est_0_04_200_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.00, 0.03)+
  geom_hline(yintercept=(performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]), linetype="dashed", color = "red")


# sample size = 500


con_est_0_04_500_all <- rbind(
  con_est_0_04_500_5_delta, 
  con_est_0_04_500_10_delta,
  con_est_0_04_500_20_delta,  
  con_est_0_04_500_5_AUC12, 
  con_est_0_04_500_10_AUC12, 
  con_est_0_04_500_20_AUC12, 
  con_est_0_04_500_5_AUC1d, 
  con_est_0_04_500_10_AUC1d, 
  con_est_0_04_500_20_AUC1d)


con_est_0_04_500_all$N.Study <- as.factor(con_est_0_04_500_all$N.Study)
con_est_0_04_500_all$Strategy <- as.factor(con_est_0_04_500_all$Strategy)


ggplot(con_est_0_04_500_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.00, 0.03)+
  geom_hline(yintercept=(performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]), linetype="dashed", color = "red")

# sample size = 1000


con_est_0_04_1000_all <- rbind(
  con_est_0_04_1000_5_delta, 
  con_est_0_04_1000_10_delta,
  con_est_0_04_1000_20_delta,  
  con_est_0_04_1000_5_AUC12, 
  con_est_0_04_1000_10_AUC12, 
  con_est_0_04_1000_20_AUC12, 
  con_est_0_04_1000_5_AUC1d, 
  con_est_0_04_1000_10_AUC1d, 
  con_est_0_04_1000_20_AUC1d)


con_est_0_04_1000_all$N.Study <- as.factor(con_est_0_04_1000_all$N.Study)
con_est_0_04_1000_all$Strategy <- as.factor(con_est_0_04_1000_all$Strategy)


ggplot(con_est_0_04_1000_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.00, 0.03)+
  geom_hline(yintercept=(performance[[2]]$auc_7_con[1]-performance[[2]]$auc_6_con[1]), linetype="dashed", color = "red")


# beta = 0.7

# sample size = 200


con_est_0_07_200_all <- rbind(
  con_est_0_07_200_5_delta, 
  con_est_0_07_200_10_delta,
  con_est_0_07_200_20_delta,  
  con_est_0_07_200_5_AUC12, 
  con_est_0_07_200_10_AUC12, 
  con_est_0_07_200_20_AUC12, 
  con_est_0_07_200_5_AUC1d, 
  con_est_0_07_200_10_AUC1d, 
  con_est_0_07_200_20_AUC1d)


con_est_0_07_200_all$N.Study <- as.factor(con_est_0_07_200_all$N.Study)
con_est_0_07_200_all$Strategy <- as.factor(con_est_0_07_200_all$Strategy)


ggplot(con_est_0_07_200_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.01, 0.05)+
  geom_hline(yintercept=(performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]), linetype="dashed", color = "red")



# sample size = 500


con_est_0_07_500_all <- rbind(
  con_est_0_07_500_5_delta, 
  con_est_0_07_500_10_delta,
  con_est_0_07_500_20_delta,  
  con_est_0_07_500_5_AUC12, 
  con_est_0_07_500_10_AUC12, 
  con_est_0_07_500_20_AUC12, 
  con_est_0_07_500_5_AUC1d, 
  con_est_0_07_500_10_AUC1d, 
  con_est_0_07_500_20_AUC1d)


con_est_0_07_500_all$N.Study <- as.factor(con_est_0_07_500_all$N.Study)
con_est_0_07_500_all$Strategy <- as.factor(con_est_0_07_500_all$Strategy)


ggplot(con_est_0_07_500_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.01, 0.05)+
  geom_hline(yintercept=(performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]), linetype="dashed", color = "red")

# sample size = 1000


con_est_0_07_1000_all <- rbind(
  con_est_0_07_1000_5_delta, 
  con_est_0_07_1000_10_delta,
  con_est_0_07_1000_20_delta,  
  con_est_0_07_1000_5_AUC12, 
  con_est_0_07_1000_10_AUC12, 
  con_est_0_07_1000_20_AUC12, 
  con_est_0_07_1000_5_AUC1d, 
  con_est_0_07_1000_10_AUC1d, 
  con_est_0_07_1000_20_AUC1d)


con_est_0_07_1000_all$N.Study <- as.factor(con_est_0_07_1000_all$N.Study)
con_est_0_07_1000_all$Strategy <- as.factor(con_est_0_07_1000_all$Strategy)


ggplot(con_est_0_07_1000_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.01, 0.05)+
  geom_hline(yintercept=(performance[[3]]$auc_7_con[1]-performance[[3]]$auc_6_con[1]), linetype="dashed", color = "red")



# beta = 0.9
# sample size = 200

con_est_0_09_200_all <- rbind(
  con_est_0_09_200_5_delta, 
  con_est_0_09_200_10_delta,
  con_est_0_09_200_20_delta,  
  con_est_0_09_200_5_AUC12, 
  con_est_0_09_200_10_AUC12, 
  con_est_0_09_200_20_AUC12, 
  con_est_0_09_200_5_AUC1d, 
  con_est_0_09_200_10_AUC1d, 
  con_est_0_09_200_20_AUC1d)


con_est_0_09_200_all$N.Study <- as.factor(con_est_0_09_200_all$N.Study)
con_est_0_09_200_all$Strategy <- as.factor(con_est_0_09_200_all$Strategy)


ggplot(con_est_0_09_200_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.03, 0.09)+
  geom_hline(yintercept=(performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]), linetype="dashed", color = "red")


# sample size = 500

con_est_0_09_500_all <- rbind(
  con_est_0_09_500_5_delta, 
  con_est_0_09_500_10_delta,
  con_est_0_09_500_20_delta,  
  con_est_0_09_500_5_AUC12, 
  con_est_0_09_500_10_AUC12, 
  con_est_0_09_500_20_AUC12, 
  con_est_0_09_500_5_AUC1d, 
  con_est_0_09_500_10_AUC1d, 
  con_est_0_09_500_20_AUC1d)


con_est_0_09_500_all$N.Study <- as.factor(con_est_0_09_500_all$N.Study)
con_est_0_09_500_all$Strategy <- as.factor(con_est_0_09_500_all$Strategy)


ggplot(con_est_0_09_500_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.03, 0.09)+
  geom_hline(yintercept=(performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]), linetype="dashed", color = "red")



# sample size = 1000

con_est_0_09_1000_all <- rbind(
  con_est_0_09_1000_5_delta, 
  con_est_0_09_1000_10_delta,
  con_est_0_09_1000_20_delta,  
  con_est_0_09_1000_5_AUC12, 
  con_est_0_09_1000_10_AUC12, 
  con_est_0_09_1000_20_AUC12, 
  con_est_0_09_1000_5_AUC1d, 
  con_est_0_09_1000_10_AUC1d, 
  con_est_0_09_1000_20_AUC1d)


con_est_0_09_1000_all$N.Study <- as.factor(con_est_0_09_1000_all$N.Study)
con_est_0_09_1000_all$Strategy <- as.factor(con_est_0_09_1000_all$Strategy)


ggplot(con_est_0_09_1000_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.03, 0.09)+
  geom_hline(yintercept=(performance[[4]]$auc_7_con[1]-performance[[4]]$auc_6_con[1]), linetype="dashed", color = "red")





# binary marker

# beta = 0
# sample size =200

bin_est_0_0_200_all <- rbind(
  bin_est_0_0_200_5_delta, 
  bin_est_0_0_200_10_delta,
  bin_est_0_0_200_20_delta,  
  bin_est_0_0_200_5_AUC12, 
  bin_est_0_0_200_10_AUC12, 
  bin_est_0_0_200_20_AUC12, 
  bin_est_0_0_200_5_AUC1d, 
  bin_est_0_0_200_10_AUC1d, 
  bin_est_0_0_200_20_AUC1d)


bin_est_0_0_200_all$N.Study <- as.factor(bin_est_0_0_200_all$N.Study)
bin_est_0_0_200_all$Strategy <- as.factor(bin_est_0_0_200_all$Strategy)


ggplot(bin_est_0_0_200_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot() +ylim(-0.002, 0.002)+
  geom_hline(yintercept=(performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]), linetype="dashed", color = "red")

# sample size = 500

bin_est_0_0_500_all <- rbind(
  bin_est_0_0_500_5_delta, 
  bin_est_0_0_500_10_delta,
  bin_est_0_0_500_20_delta,  
  bin_est_0_0_500_5_AUC12, 
  bin_est_0_0_500_10_AUC12, 
  bin_est_0_0_500_20_AUC12, 
  bin_est_0_0_500_5_AUC1d, 
  bin_est_0_0_500_10_AUC1d, 
  bin_est_0_0_500_20_AUC1d)


bin_est_0_0_500_all$N.Study <- as.factor(bin_est_0_0_500_all$N.Study)
bin_est_0_0_500_all$Strategy <- as.factor(bin_est_0_0_500_all$Strategy)


ggplot(bin_est_0_0_500_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(-0.002, 0.002)+
  geom_hline(yintercept=(performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]), linetype="dashed", color = "red")


# sample size = 1000

bin_est_0_0_1000_all <- rbind(
  bin_est_0_0_1000_5_delta, 
  bin_est_0_0_1000_10_delta,
  bin_est_0_0_1000_20_delta,  
  bin_est_0_0_1000_5_AUC12, 
  bin_est_0_0_1000_10_AUC12, 
  bin_est_0_0_1000_20_AUC12, 
  bin_est_0_0_1000_5_AUC1d, 
  bin_est_0_0_1000_10_AUC1d, 
  bin_est_0_0_1000_20_AUC1d)


bin_est_0_0_1000_all$N.Study <- as.factor(bin_est_0_0_1000_all$N.Study)
bin_est_0_0_1000_all$Strategy <- as.factor(bin_est_0_0_1000_all$Strategy)


ggplot(bin_est_0_0_1000_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(-0.002, 0.002)+
  geom_hline(yintercept=(performance[[1]]$auc_7_bin[1]-performance[[1]]$auc_6_bin[1]), linetype="dashed", color = "red")



# beta = 0.4
# sample size = 200


bin_est_0_04_200_all <- rbind(
  bin_est_0_04_200_5_delta, 
  bin_est_0_04_200_10_delta,
  bin_est_0_04_200_20_delta,  
  bin_est_0_04_200_5_AUC12, 
  bin_est_0_04_200_10_AUC12, 
  bin_est_0_04_200_20_AUC12, 
  bin_est_0_04_200_5_AUC1d, 
  bin_est_0_04_200_10_AUC1d, 
  bin_est_0_04_200_20_AUC1d)


bin_est_0_04_200_all$N.Study <- as.factor(bin_est_0_04_200_all$N.Study)
bin_est_0_04_200_all$Strategy <- as.factor(bin_est_0_04_200_all$Strategy)


ggplot(bin_est_0_04_200_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.00, 0.03)+
  geom_hline(yintercept=(performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]), linetype="dashed", color = "red")


# sample size = 500


bin_est_0_04_500_all <- rbind(
  bin_est_0_04_500_5_delta, 
  bin_est_0_04_500_10_delta,
  bin_est_0_04_500_20_delta,  
  bin_est_0_04_500_5_AUC12, 
  bin_est_0_04_500_10_AUC12, 
  bin_est_0_04_500_20_AUC12, 
  bin_est_0_04_500_5_AUC1d, 
  bin_est_0_04_500_10_AUC1d, 
  bin_est_0_04_500_20_AUC1d)


bin_est_0_04_500_all$N.Study <- as.factor(bin_est_0_04_500_all$N.Study)
bin_est_0_04_500_all$Strategy <- as.factor(bin_est_0_04_500_all$Strategy)


ggplot(bin_est_0_04_500_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.00, 0.03)+
  geom_hline(yintercept=(performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]), linetype="dashed", color = "red")

# sample size = 1000


bin_est_0_04_1000_all <- rbind(
  bin_est_0_04_1000_5_delta, 
  bin_est_0_04_1000_10_delta,
  bin_est_0_04_1000_20_delta,  
  bin_est_0_04_1000_5_AUC12, 
  bin_est_0_04_1000_10_AUC12, 
  bin_est_0_04_1000_20_AUC12, 
  bin_est_0_04_1000_5_AUC1d, 
  bin_est_0_04_1000_10_AUC1d, 
  bin_est_0_04_1000_20_AUC1d)


bin_est_0_04_1000_all$N.Study <- as.factor(bin_est_0_04_1000_all$N.Study)
bin_est_0_04_1000_all$Strategy <- as.factor(bin_est_0_04_1000_all$Strategy)


ggplot(bin_est_0_04_1000_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.00, 0.03)+
  geom_hline(yintercept=(performance[[2]]$auc_7_bin[1]-performance[[2]]$auc_6_bin[1]), linetype="dashed", color = "red")


# beta = 0.7

# sample size = 200


bin_est_0_07_200_all <- rbind(
  bin_est_0_07_200_5_delta, 
  bin_est_0_07_200_10_delta,
  bin_est_0_07_200_20_delta,  
  bin_est_0_07_200_5_AUC12, 
  bin_est_0_07_200_10_AUC12, 
  bin_est_0_07_200_20_AUC12, 
  bin_est_0_07_200_5_AUC1d, 
  bin_est_0_07_200_10_AUC1d, 
  bin_est_0_07_200_20_AUC1d)


bin_est_0_07_200_all$N.Study <- as.factor(bin_est_0_07_200_all$N.Study)
bin_est_0_07_200_all$Strategy <- as.factor(bin_est_0_07_200_all$Strategy)


ggplot(bin_est_0_07_200_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.01, 0.05)+
  geom_hline(yintercept=(performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]), linetype="dashed", color = "red")



# sample size = 500


bin_est_0_07_500_all <- rbind(
  bin_est_0_07_500_5_delta, 
  bin_est_0_07_500_10_delta,
  bin_est_0_07_500_20_delta,  
  bin_est_0_07_500_5_AUC12, 
  bin_est_0_07_500_10_AUC12, 
  bin_est_0_07_500_20_AUC12, 
  bin_est_0_07_500_5_AUC1d, 
  bin_est_0_07_500_10_AUC1d, 
  bin_est_0_07_500_20_AUC1d)


bin_est_0_07_500_all$N.Study <- as.factor(bin_est_0_07_500_all$N.Study)
bin_est_0_07_500_all$Strategy <- as.factor(bin_est_0_07_500_all$Strategy)


ggplot(bin_est_0_07_500_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.01, 0.05)+
  geom_hline(yintercept=(performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]), linetype="dashed", color = "red")

# sample size = 1000


bin_est_0_07_1000_all <- rbind(
  bin_est_0_07_1000_5_delta, 
  bin_est_0_07_1000_10_delta,
  bin_est_0_07_1000_20_delta,  
  bin_est_0_07_1000_5_AUC12, 
  bin_est_0_07_1000_10_AUC12, 
  bin_est_0_07_1000_20_AUC12, 
  bin_est_0_07_1000_5_AUC1d, 
  bin_est_0_07_1000_10_AUC1d, 
  bin_est_0_07_1000_20_AUC1d)


bin_est_0_07_1000_all$N.Study <- as.factor(bin_est_0_07_1000_all$N.Study)
bin_est_0_07_1000_all$Strategy <- as.factor(bin_est_0_07_1000_all$Strategy)


ggplot(bin_est_0_07_1000_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.01, 0.05)+
  geom_hline(yintercept=(performance[[3]]$auc_7_bin[1]-performance[[3]]$auc_6_bin[1]), linetype="dashed", color = "red")



# beta = 0.9
# sample size = 200

bin_est_0_09_200_all <- rbind(
  bin_est_0_09_200_5_delta, 
  bin_est_0_09_200_10_delta,
  bin_est_0_09_200_20_delta,  
  bin_est_0_09_200_5_AUC12, 
  bin_est_0_09_200_10_AUC12, 
  bin_est_0_09_200_20_AUC12, 
  bin_est_0_09_200_5_AUC1d, 
  bin_est_0_09_200_10_AUC1d, 
  bin_est_0_09_200_20_AUC1d)


bin_est_0_09_200_all$N.Study <- as.factor(bin_est_0_09_200_all$N.Study)
bin_est_0_09_200_all$Strategy <- as.factor(bin_est_0_09_200_all$Strategy)


ggplot(bin_est_0_09_200_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.03, 0.09)+
  geom_hline(yintercept=(performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]), linetype="dashed", color = "red")


# sample size = 500

bin_est_0_09_500_all <- rbind(
  bin_est_0_09_500_5_delta, 
  bin_est_0_09_500_10_delta,
  bin_est_0_09_500_20_delta,  
  bin_est_0_09_500_5_AUC12, 
  bin_est_0_09_500_10_AUC12, 
  bin_est_0_09_500_20_AUC12, 
  bin_est_0_09_500_5_AUC1d, 
  bin_est_0_09_500_10_AUC1d, 
  bin_est_0_09_500_20_AUC1d)


bin_est_0_09_500_all$N.Study <- as.factor(bin_est_0_09_500_all$N.Study)
bin_est_0_09_500_all$Strategy <- as.factor(bin_est_0_09_500_all$Strategy)


ggplot(bin_est_0_09_500_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.03, 0.09)+
  geom_hline(yintercept=(performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]), linetype="dashed", color = "red")



# sample size = 1000

bin_est_0_09_1000_all <- rbind(
  bin_est_0_09_1000_5_delta, 
  bin_est_0_09_1000_10_delta,
  bin_est_0_09_1000_20_delta,  
  bin_est_0_09_1000_5_AUC12, 
  bin_est_0_09_1000_10_AUC12, 
  bin_est_0_09_1000_20_AUC12, 
  bin_est_0_09_1000_5_AUC1d, 
  bin_est_0_09_1000_10_AUC1d, 
  bin_est_0_09_1000_20_AUC1d)


bin_est_0_09_1000_all$N.Study <- as.factor(bin_est_0_09_1000_all$N.Study)
bin_est_0_09_1000_all$Strategy <- as.factor(bin_est_0_09_1000_all$Strategy)


ggplot(bin_est_0_09_1000_all, aes(x=N.Study, y=Est, fill=Strategy)) + 
  geom_boxplot()+ylim(0.03, 0.09)+
  geom_hline(yintercept=(performance[[4]]$auc_7_bin[1]-performance[[4]]$auc_6_bin[1]), linetype="dashed", color = "red")






 save.image("results_correlation0.RData")

# save.image("results_correlation02.RData")


 #save.image("results_correlation04.RData")





     