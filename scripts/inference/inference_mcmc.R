#### Setup #### 

library(purrr)
library(data.table)
library(Matrix)
library(HMM)
library(parallel)
library(mgcv)
library(dplyr)
library(splines)

ncl = detectCores()-1
cl = makeCluster(ncl, type = "FORK")
clusterSetRNGStream(cl)

source("scripts/inference/inference_functions.R")
clusterEvalQ(cl, {source("scripts/inference/inference_functions.R")})

dataTot = readRDS("data/clean/dataTot.rds")

#### Study 1 ENA #### 

configs = expand.grid("model_index"=1:5,"prop_of_observed"=1,"contr_simulated"=1,
                      "N"=5000, "obsburn"=4700, "size_simu"=1,"iter_metropolis" = 30,"ena_pop"=T)

groups = c("ena")
formulasList = list()

formulasList[[1]] = list("gamma0"= traj_lead ~ 1,
                         "gamma1"= I(1-traj_lead) ~ 1,
                         "lambda"= trace_obs ~ 1)

formulasList[[2]] = list("gamma0"= traj_lead ~ 1,
                         "gamma1"= I(1-traj_lead) ~ 1,
                         "lambda"= trace_obs ~ lagged1+lagged2+sqrt(lagged3)+I(lagged1==0)+as.double(lagged1==0):sqrt(lagged3))

# Run again 3
formulasList[[3]] = list("gamma0"= traj_lead ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T),
                         "gamma1"= I(1-traj_lead) ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T),
                         "lambda"= trace_obs ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+
                          lagged1+lagged2+sqrt(lagged3)+
                          as.double(lagged1==0):bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))

formulasList[[4]] = list("gamma0"= traj_lead ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+I(period-age_ena),
                         "gamma1"= I(1-traj_lead) ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+I(period-age_ena),
                         "lambda"= trace_obs ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+
                          lagged1+lagged2+sqrt(lagged3)+I(period-age_ena)+
                          as.double(lagged1==0):bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))

formulasList[[5]] = list("gamma0"= traj_lead ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)
                           + as.double(round(period*64) %in% c(9,10,23,24,33,34,43,44,53,54,63,64)),
                         "gamma1"= I(1-traj_lead) ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)
                           + as.double(round(period*64) %in% c(9,10,23,24,33,34,43,44,53,54,63,64)),
                         "lambda"= trace_obs ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+
                          lagged1+lagged2+sqrt(lagged3)+
                          as.double(round(period*64) %in% c(10,11,24,25,34,35,44,45,54,55,64,65))+
                          as.double(lagged1==0):bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))

i = 1
for (i in c(3)){
  print(paste0("Config ",i))
  
  config = as.list(configs[i,])[c("model_index","N","obsburn","size_simu","iter_metropolis","prop_of_observed","contr_simulated","ena_pop")]
  config[["formulas"]] = formulasList[[config[["model_index"]]]]
  config[["groups"]] = groups
  config[["dataTot"]] = dataTot
  
  result = do.call(inference_mod, config)
  result[["config"]] = config
  
  saveRDS(result[c("params","config","obs","given_data")], paste0("data/results/ena_config_",i,".rds"))
}

#### Study 2 Important organizations ####

configs = expand.grid("model_index"=1:2,"prop_of_observed"=1,"contr_simulated"=1,
                      "N"=5000, "obsburn"=4700, "size_simu"=1,"iter_metropolis" = 30,"ena_pop"=F, "remove_early"=F)
groupsTot =  c("cprefet","ccomptes","iga","ce","cdiplo","igf","dgfip","insee",
               "igas","dgtresor")

groupsList = list()

groupsList[[1]] = groupsTot
groupsList[[2]] = groupsTot

formulasList = list()

formulasList[[1]] = list("gamma0"= traj_lead ~ -1+group,
                         "gamma1"= I(1-traj_lead) ~ -1+group,
                         "lambda"= trace_obs ~ -1+group
                         +group:(lagged1+lagged2+sqrt(lagged3))+as.double(lagged1==0):bs(age_max,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))

formulasList[[2]] = list("gamma0"= traj_lead ~ -1+group+group:sex,
                         "gamma1"= I(1-traj_lead) ~ -1+group+group:sex,
                         "lambda"= trace_obs ~ -1+group+group:sex
                         +group:(lagged1+lagged2+sqrt(lagged3))+as.double(lagged1==0):bs(age_max,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))

for (i in c(2)){
  set.seed(42)
  print(paste0("Config ",i))
  
  config = as.list(configs[i,])[c("model_index","N","obsburn","size_simu","iter_metropolis","prop_of_observed","contr_simulated","ena_pop","remove_early")]
  config[["formulas"]] = formulasList[[config[["model_index"]]]]
  config[["groups"]] = groupsList[[config[["model_index"]]]]
  config[["dataTot"]] = dataTot


  result = do.call(inference_mod, config)
  result[["config"]] = config
  
  saveRDS(result[c("params","config","obs","given_data")], paste0("data/results/groups_config_",i,".rds"))
}


#### Descriptive statistics ####

