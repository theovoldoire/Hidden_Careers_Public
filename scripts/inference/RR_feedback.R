library(purrr)
library(data.table)
library(Matrix)
library(HMM)
library(parallel)
library(mgcv)
library(dplyr)
library(splines)
library(latex2exp)
library(kableExtra)

n1 = 5000
n2 = 4700
ncl = detectCores()-1
cl = makeCluster(ncl, type = "FORK")
clusterSetRNGStream(cl)

source("scripts/inference/inference_functions.R")
clusterEvalQ(cl, {source("scripts/inference/inference_functions.R")})

cleanCV = function(res_CV){
  
  res_CVb = merge(res_CV$given_data, 
                  res_CV$obs[,.(traj_pred=mean(as.double(traj))),.(ident,age_max)], 
                  by = c("ident","age_max"))[period!=0,c("ident","period","traj_obs","traj_pred","profile_observed","profile_included","group")]
  res_CVb = res_CVb[traj_obs != 98 & profile_included == FALSE,]
  
  res_CVb$traj_pred[res_CVb$traj_pred==0] = 0.01
  res_CVb$traj_pred[res_CVb$traj_pred==1] = 0.99
  res_CVb$loss = -(as.double(res_CVb$traj_obs) * log(res_CVb$traj_pred) + 
                     (1-as.double(res_CVb$traj_obs)) * log(1-res_CVb$traj_pred))
  
  return(res_CVb)
}

#### 1.1. Analysis with another prior ####

dataTot = readRDS("data/clean/dataTot.rds")

formulas4 = list("gamma0"= traj_lead ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+I(period-age_ena),
                 "gamma1"= I(1-traj_lead) ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+I(period-age_ena),
                 "lambda"= trace_obs ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+
                   lagged1+lagged2+sqrt(lagged3)+I(period-age_ena)+
                   as.double(lagged1==0):bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))

res_prior4 = inference_mod(formulas = formulas4, groups = ("ena"), dataTot = dataTot, model_index = "prior1",
                           N = 300, size_simu = 1, iter_metropolis = 30, obsburn = n2,
                           prop_of_observed = 1, contr_simulated = 1, remove_early = F, ena_pop = T)

saveRDS(res_prior4, "data/revise_resubmit/ena_prior4.rds")

res_prior4 = readRDS("data/revise_resubmit/ena_prior4.rds")

#### 1.2. Analysis with different weight ####

dataTot = readRDS("data/clean/dataTot.rds")

formulas4 = list("gamma0"= traj_lead ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+I(period-age_ena),
                 "gamma1"= I(1-traj_lead) ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+I(period-age_ena),
                 "lambda"= trace_obs ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+
                   lagged1+lagged2+sqrt(lagged3)+I(period-age_ena)+
                   as.double(lagged1==0):bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))

set.seed(42)
res_phi06 = inference_mod(formulas = formulas4, groups = ("ena"), dataTot = dataTot, model_index = "check_phi_06",
                          N = 300, size_simu = 1, iter_metropolis = 30, obsburn = 150, 
                          prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = T)
saveRDS(res_phi06, "data/revise_resubmit/phi_06.rds")

set.seed(42)
res_phi07 = inference_mod(formulas = formulas4, groups = ("ena"), dataTot = dataTot, model_index = "check_phi_07",
                           N = 300, size_simu = 1, iter_metropolis = 30, obsburn = 150, 
                           prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = T)

saveRDS(res_phi07, "data/revise_resubmit/phi_07.rds")

res_phi08 = inference_mod(formulas = formulas4, groups = ("ena"), dataTot = dataTot, model_index = "NA",
                          N = 300, size_simu = 1, iter_metropolis = 30, obsburn = 150, 
                          prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = T)
saveRDS(res_phi08, "data/revise_resubmit/phi_08.rds")

res_phi09 = inference_mod(formulas = formulas4, groups = ("ena"), dataTot = dataTot, model_index = "check_phi_09",
                          N = 300, size_simu = 1, iter_metropolis = 30, obsburn = 150, 
                          prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = T)
saveRDS(res_phi09, "data/revise_resubmit/phi_09.rds")

res_phi06 = readRDS("data/revise_resubmit/phi_06.rds")
res_phi07 = readRDS("data/revise_resubmit/phi_07.rds")
res_phi08 = readRDS("data/revise_resubmit/phi_08.rds")
res_phi09 = readRDS("data/revise_resubmit/phi_09.rds")

(res_phi06 |> cleanCV())$loss |> mean()
(res_phi07 |> cleanCV())$loss |> mean()
(res_phi08 |> cleanCV())$loss |> mean()
(res_phi09 |> cleanCV())$loss |> mean()

cor(res_phi06$given_data$lagged1[!is.na(res_phi06$given_data$lagged1)], res_phi06$given_data$lagged2[!is.na(res_phi06$given_data$lagged1)])

cor(res_phi09$given_data$lagged1[!is.na(res_phi09$given_data$lagged1)], res_phi09$given_data$lagged2[!is.na(res_phi09$given_data$lagged1)])


formulas4_all = list("gamma0"= traj_lead ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+I(period-age_ena),
                 "gamma1"= I(1-traj_lead) ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+I(period-age_ena),
                 "lambda"= trace_obs ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+
                   lagged1+lagged2+lagged2b+sqrt(lagged3)+I(period-age_ena)+
                   as.double(lagged1==0):bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))

res_phi_all = inference_mod(formulas = formulas4_all, groups = ("ena"), dataTot = dataTot, model_index = "check_phi_all",
                          N = 300, size_simu = 1, iter_metropolis = 30, obsburn = 150, 
                          prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = T)
(res_phi_all |> cleanCV())$loss |> mean()

saveRDS(res_phi_all, "data/revise_resubmit/phi_all.rds")

(res_phi_all |> cleanCV())$loss |> mean()
quantile(res_phi_all$params$gamma0[6,100:300], c(0.05, 0.95))

 #### 2. Simulation study on synthetic data ####

#### 2.1 Simulating the data ####

simulate_data = function(pars, prob_obs = 0.2, nind = 1000){

  sim_data = list()
  
  nper = 61
  len_traj = (nper-1)/2
  phi = 0.8; phiP = phi^-seq(1, len_traj)

  for (ident in 1:nind){
    
    cohort = runif(n = 1, min = 0, max = 0.5)
    age_max = (0:(len_traj-1))/(nper-1)
    
    lambda_ind = rnorm(n = 1, mean = 0, sd = 0.1)
    base_lambda = pars[1] + lambda_ind + pars[2] * age_max + pars[3] * cohort
    gamma0 = plogis(pars[4] + pars[5] * age_max + pars[6] * cohort)
    gamma1 = plogis(pars[7] + pars[8] * age_max + pars[9] * cohort)
      
    lambda = rep(NA,len_traj)
    traj = rep(NA,len_traj)
    traces = rep(NA,len_traj)
    traj[1] = 0
    traces[1] = 0
    
    for (t in 1:(len_traj-1)){
      
      if (traj[t] == 0) traj[t+1] = sample(x = c(0,1), size = 1, prob = c(1-gamma0[t],gamma0[t]))
      if (traj[t] == 1) traj[t+1] = sample(x = c(0,1), size = 1, prob = c(gamma1[t],1-gamma1[t]))
      
      lagged1 = lag(cumsum(traces)/(1:length(traces)))
      lagged2 = lag(cumsum(as.double(traces)*phiP)/cumsum(phiP))
      lambda[t+1] = plogis(base_lambda[t+1] + pars[10] * lagged1[t+1] + pars[11] * lagged2[t+1])
      
      if (traj[t+1] == 0){
        traces[t+1] = sample(x=c(0,1),size=1,prob=c(1-lambda[t+1],lambda[t+1]))
      } else {traces[t+1] = 0}
      
    }
    lagged1[1] = 0; lagged2[1] = 0
    
    profile_pres = as.double(ident <= (nind * prob_obs))
    if (profile_pres == 1){
      traj_obs = traj
    } else {
      traj_obs = rep(98,len_traj)
      #traj_obs[1] = 0
    }
    sim_data[[ident]] = data.table(ident = rep(ident,len_traj), group = rep("all",len_traj),
                                   traj_true = traj, trace_obs = traces, traj_obs = traj_obs, period = cohort+age_max, 
                                 lambda = lambda, gamma0 = gamma0, gamma1 = gamma1,
                                 lagged1 = lagged1, lagged2 = lagged2,
                                 age_ena = rep(NA,len_traj), age_max = age_max, profile_present = rep(profile_pres, len_traj))
  }
  names(sim_data) = as.character(1:nind)
  return(sim_data)
}

#### 2.2. Running the inference ####

#### 2.2.1. More statistical power ####

set.seed(42)
pars1 = c(-1.5, 0, 0, # base_lambda: intercept, age_max, cohort
          -4, 0, 1,  # gamma 0: intercept, age_max, cohort
          -4, 0, 0, # gamma 1: intercept, age_max, cohort
          0, 0) # lambda: A(1), A(0.8)

sim_data_small = simulate_data(pars = pars1, prob_obs = 1, nind = 200)
sim_data_medium = simulate_data(pars = pars1, prob_obs = 1, nind = 500)
sim_data_large = simulate_data(pars = pars1, prob_obs = 0.2, nind = 1000)

formulas_st = list("gamma0"= traj_lead ~ I(period-age_max),
                   "gamma1"= I(1-traj_lead) ~ 1,
                   "lambda"= trace_obs ~ 1)

res_small =  inference_mod(formulas = formulas_st, dataTot = sim_data_small,  groups = ("all"), 
                           N = 5000, size_simu = 1, iter_metropolis = 30, obsburn = Inf,
                           prop_of_observed = 1, contr_simulated = 1, remove_early = F, ena_pop = F)

res_medium = inference_mod(formulas = formulas_st, dataTot = sim_data_medium,  groups = ("all"), 
                           N = 5000, size_simu = 1, iter_metropolis = 30, obsburn = Inf,
                           prop_of_observed = 1, contr_simulated = 1, remove_early = F, ena_pop = F)

res_large =  inference_mod(formulas = formulas_st, dataTot = sim_data_large,  groups = ("all"), 
                           N = 5000, size_simu = 1, iter_metropolis = 30, obsburn = Inf,
                           prop_of_observed = 1, contr_simulated = 1, remove_early = F, ena_pop = F)

saveRDS(list("small"=res_small,"medium"=res_medium,"large"=res_large), "data/revise_resubmit/syn_power.rds")

#### 2.2.2. Issue of misspecification ####

set.seed(42)
pars1 = c(-2, 0, 0, # base_lambda: intercept, age_max, cohort
          -4, 1, 0,  # gamma 0: intercept, age_max, cohort
          -4, 0, 0, # gamma 1: intercept, age_max, cohort
          0, 2.5) # lambda: A(1), A(0.8)

sim_data1 = simulate_data(pars = pars1, nind = 2000, prob_obs = 0.3)

formulas_right = list("gamma0"= traj_lead ~ age_max,
                      "gamma1"= I(1-traj_lead) ~ 1,
                      "lambda"= trace_obs ~ lagged2)

formulas_mis = list("gamma0"= traj_lead ~ age_max,
                    "gamma1"= I(1-traj_lead) ~ 1,
                    "lambda"= trace_obs ~ 1)

res_right =  inference_mod(formulas = formulas_right,dataTot = sim_data1,  groups = ("all"), 
              N = n1, size_simu = 1, iter_metropolis = 30, obsburn = n2,
              prop_of_observed = 1, contr_simulated = 1, remove_early = F, ena_pop = F)

res_mis =  inference_mod(formulas = formulas_mis,dataTot = sim_data1,  groups = ("all"), 
                           N = n1, size_simu = 1, iter_metropolis = 30, obsburn = n2,
                           prop_of_observed = 1, contr_simulated = 1, remove_early = F, ena_pop = F)

saveRDS(list("well_specified"=res_right, "misspecified"=res_mis), file = "data/revise_resubmit/syn_miss.rds")


#### 2.2.3. Figures ####

# Figure 1
ll = readRDS("data/revise_resubmit/syn_power.rds")
res_small = ll$small
res_medium = ll$medium
res_large = ll$large

# Crafting the table

result = matrix(nrow = 6, ncol =  4)
colnames(result) = c("Case 1", "Case 2", "Case 3", "True value")
rownames(result) = c("n", "$n_\\text{obs}$", "$\\gamma^0$, Intercept", "$\\gamma^0$, Cohort",
                     "$\\gamma^1$", "$\\lambda$")
result[1,] = c("200", "500", "1000", "")
result[2,] = c("200", "500", "200", "")
Cinterval = c(0.05,0.95)
Dinterval = 201:2000

mm = function(x){
  return(paste0("[",
                paste0(round(quantile(x, Cinterval),2), collapse = ","),
                "]"))
}

result[3,] = c(mm(res_small$params$gamma0[1,Dinterval]),
               mm(res_medium$params$gamma0[1,Dinterval]),
               mm(res_large$params$gamma0[1,Dinterval]),
               "-4")
result[4,] = c(mm(res_small$params$gamma0[2,Dinterval]),
               mm(res_medium$params$gamma0[2,Dinterval]),
               mm(res_large$params$gamma0[2,Dinterval]),
               "1")
result[5,] = c(mm(res_small$params$gamma1[1,Dinterval]),
               mm(res_medium$params$gamma1[1,Dinterval]),
               mm(res_large$params$gamma1[1,Dinterval]),
               "-4")
result[6,] = c(mm(res_small$params$lambda[1,Dinterval]),
               mm(res_medium$params$lambda[1,Dinterval]),
               mm(res_large$params$lambda[1,Dinterval]),
               "-1.5")

kbl(result, format = "latex", booktabs=T, longtable = T, align = "r",
    row.names = T, escape = F, linesep = "", col.names = colnames(result)) |>
  kable_styling(font_size = 7) |>
  save_kable(file = "documents/figures/RR_synth1.tex") 

# Figure 2

res_right = readRDS("data/revise_resubmit/syn_miss.rds")[[1]]
res_mis = readRDS("data/revise_resubmit/syn_miss.rds")[[2]]

result = matrix(nrow = 5, ncol =  3)
colnames(result) = c("Well-specified", "Misspecified", "True value")
rownames(result) = c("$\\gamma^0$, Intercept", "$\\gamma^0$, Age",
                     "$\\gamma^1$, Intercept", "$\\lambda$, Intercept", "$\\lambda$, $A(0.8)$")
Cinterval = c(0.05,0.95)
Dinterval = 201:2000

result[1,]=c(mm(res_right$params$gamma0[1,Dinterval]),
             mm(res_mis$params$gamma0[1,Dinterval]),
             "-4")
result[2,]=c(mm(res_right$params$gamma0[2,Dinterval]),
             mm(res_mis$params$gamma0[2,Dinterval]),
             "1")
result[3,]=c(mm(res_right$params$gamma1[1,Dinterval]),
             mm(res_mis$params$gamma1[1,Dinterval]),
             "-4")
result[4,]=c(mm(res_right$params$lambda[1,100:200]),
             mm(res_mis$params$lambda[1,100:200]),
             "-2")
result[5,]=c(mm(res_right$params$lambda[2,100:200]),
             "",
             "2.5")

kbl(result, format = "latex", booktabs=T, longtable = T, align = "r",
    row.names = T, escape = F, linesep = "", col.names = colnames(result)) |>
  kable_styling(font_size = 7) |>
  save_kable(file = "documents/figures/RR_synth2.tex") 


#### 3. Cross validation on real data ####

#### 3.1. Running the inference with hold-out (for ENA) ####

dataTot = readRDS("data/clean/dataTot.rds")

formulasCV1 = list("gamma0"= traj_lead ~ 1,
                         "gamma1"= I(1-traj_lead) ~ 1,
                         "lambda"= trace_obs ~ 1)

formulasCV2 = list("gamma0"= traj_lead ~ 1,
                         "gamma1"= I(1-traj_lead) ~ 1,
                         "lambda"= trace_obs ~ lagged1+lagged2+sqrt(lagged3)+I(lagged1==0)+as.double(lagged1==0):sqrt(lagged3))

formulasCV3 = list("gamma0"= traj_lead ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T),
                   "gamma1"= I(1-traj_lead) ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T),
                   "lambda"= trace_obs ~ -1+bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T)+
                     lagged1+lagged2+sqrt(lagged3)+
                     as.double(lagged1==0):bs(age_ena,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))

set.seed(42)
res_CV1 = inference_mod(formulas = formulasCV1, groups = ("ena"), dataTot = dataTot, model_index = 3,
                        N = n1, size_simu = 1, iter_metropolis = 30, obsburn = n2,
                        prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = T)
saveRDS(res_CV1, "data/revise_resubmit/ena_CV1.rds")

res_CV1b = inference_mod(formulas = formulasCV1, groups = ("ena"), dataTot = dataTot, model_index = "CV_empty",
                       N = n1, size_simu = 1, iter_metropolis = 30, obsburn = n2,
                       prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = T)
saveRDS(res_CV1b, "data/revise_resubmit/ena_CV1b.rds")

set.seed(42)
res_CV2 = inference_mod(formulas = formulasCV2, groups = ("ena"), dataTot = dataTot, model_index = 3,
                         N = n1, size_simu = 1, iter_metropolis = 30, obsburn = n2,
                         prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = T)

saveRDS(res_CV2, "data/revise_resubmit/ena_CV2.rds")

set.seed(42)
res_CV3 = inference_mod(formulas = formulasCV3, groups = ("ena"), dataTot = dataTot, model_index = 3,
                        N = n1, size_simu = 1, iter_metropolis = 30, obsburn = n2,
                        prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = T)
saveRDS(res_CV3, "data/revise_resubmit/ena_CV3.rds")


# Which clipping


res_CV1 = readRDS("data/revise_resubmit/ena_CV1.rds") |> cleanCV()
res_CV1b = readRDS("data/revise_resubmit/ena_CV1b.rds") |> cleanCV()
res_CV2 = readRDS("data/revise_resubmit/ena_CV2.rds") |> cleanCV()
res_CV3 = readRDS("data/revise_resubmit/ena_CV3.rds") |> cleanCV()

mean(res_CV1$loss)
mean(res_CV1b$loss)
mean(res_CV2$loss)
mean(res_CV3$loss)


#### 3.2. Running the analysis with hold out (2nd population) ####

n1 = 2000
n2 = 1700

dataTot = readRDS("data/clean/dataTot.rds")
groupsTot =  c("cprefet","ccomptes","iga","ce","cdiplo","igf","insee","igas","dgtresor", "dgfip")

formulasCV1 = list("gamma0"= traj_lead ~ -1+group,
     "gamma1"= I(1-traj_lead) ~ -1+group,
     "lambda"= trace_obs ~ -1+group
     +group:(lagged1+lagged2+sqrt(lagged3))+as.double(lagged1==0):bs(age_max,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))

formulasCV2 = list("gamma0"= traj_lead ~ -1+group+group:sex,
                   "gamma1"= I(1-traj_lead) ~ -1+group+group:sex,
                   "lambda"= trace_obs ~ -1+group+group:sex
                   +group:(lagged1+lagged2+sqrt(lagged3))+as.double(lagged1==0):bs(age_max,df=5,degree=2,Boundary.knots=c(0,1),intercept=T))


set.seed(4)
res1 = inference_mod(formulas = formulasCV1, groups =groupsTot, dataTot = dataTot, model_index = 3,
              N = n1, size_simu = 1, iter_metropolis = 30, obsburn = n2,
              prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = F)
saveRDS(res1, "data/revise_resubmit/groups_CV1.rds")

set.seed(4)
res2 = inference_mod(formulas = formulasCV2, groups =groupsTot, dataTot = dataTot, model_index = 3,
                     N = n1, size_simu = 1, iter_metropolis = 30, obsburn = n2,
                     prop_of_observed = 0.8, contr_simulated = 1, remove_early = F, ena_pop = F)
saveRDS(res2, "data/revise_resubmit/groups_CV2.rds")

res_CV1 = readRDS("data/revise_resubmit/groups_CV1.rds") |> cleanCV()
res_CV2 = readRDS("data/revise_resubmit/groups_CV2.rds") |> cleanCV()

mean(res_CV1$loss)
mean(res_CV2$loss)


#### 4. Motivating data augmentation ####

# Study 1 
dataTot = readRDS("data/clean/dataTot.rds")
dataTot2 = parLapply(cl, dataTot, setting_up_data, NPER=65,
                     groups = c("ena"), ena_pop = T, remove_early = F,
                     prop_of_observed = 1, contr_simulated = 1)
dataPop = list()
for (x in names(dataTot)) if (nrow(dataTot2[[x]])>1) dataPop[[x]] = dataTot2[[x]]
data1 = rbindlist(dataPop)
table(data1$trace) |> prop.table()

# Study 2
dataTot = readRDS("data/clean/dataTot.rds")
dataTot2 = parLapply(cl, dataTot, setting_up_data, NPER=65,
                     groups = c("cprefet","ccomptes","iga","ce","cdiplo","igf","dgfip","insee",
                                "igas","dgtresor"), ena_pop = F, remove_early = F,
                     prop_of_observed = 1, contr_simulated = 1)
dataPop = list()
for (x in names(dataTot)) if (nrow(dataTot2[[x]])>1) dataPop[[x]] = dataTot2[[x]]
data2 = rbindlist(dataPop)
table(data2$trace) |> prop.table()


# Why not only work on X? 

# Motivating the data augmentation

# Proportion of missing values in total
round(prop.table(table(data$traj_obs))*100,1)

# Proportion of individuals with fully observed trajectories
1-mean(data[,.(mean(traj_obs==98)),.(ident)]$V1 == 0)

# Histogram of observed values 
hist(data[,.(mean(traj_obs==98)),.(ident)]$V1, 
     xlab = "Proportion of missing values",
     main = "Histogram of the proportion of missing values. One observation is corresponds to the proportion
     of missing values for one individual.")

#### 5. Convergence analysis ####

# Just a check with two chains for the main model
set.seed(1)
res1 = inference_mod(formulas = formulas4, groups = ("ena"), dataTot = dataTot, model_index = "",
                           N = 300, size_simu = 1, iter_metropolis = 30, obsburn = n2,
                           prop_of_observed = 1, contr_simulated = 1, remove_early = F, ena_pop = T)
set.seed(2)
res2 = inference_mod(formulas = formulas4, groups = ("ena"), dataTot = dataTot, model_index = "",
                     N = 300, size_simu = 1, iter_metropolis = 30, obsburn = n2,
                     prop_of_observed = 1, contr_simulated = 1, remove_early = F, ena_pop = T)

sims = cbind(res1$params$gamma0[6,100:300], res2$params$gamma0[6,100:300])
saveRDS(sims, "data/revise_resubmit/convergence_analysis.rds")
Rhat(sims)


observe_rhats = function(rr){
  conv = list()
  span = 201:5000
  conv[[1]] = apply(t(rr$params$gamma0[,span]), FUN = Rhat, MARGIN = 2)
  conv[[2]] = apply(t(rr$params$gamma1[,span]), FUN = Rhat, MARGIN = 2)
  conv[[3]] = apply(t(rr$params$lambda[,span]), FUN = Rhat, MARGIN = 2)
  return(unlist(conv))
}
rr = readRDS("data/results/ena_config_2.rds")
# Rhat(rr$params$gamma0[,span])
x3 = observe_rhats(readRDS("data/results/ena_config_3.rds"))
x4 = observe_rhats(readRDS("data/results/ena_config_4.rds"))
x5 = observe_rhats(readRDS("data/results/ena_config_5.rds"))

x6 = observe_rhats(readRDS("data/results/groups_config_1.rds"))
x7 = observe_rhats(readRDS("data/results/groups_config_2.rds"))

max(c(x3, x4, x5, x6, x7))

Rhat(rr2$params$gamma0[10,201:5000])
