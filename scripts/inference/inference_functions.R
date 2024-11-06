
  
  simulation_mod = function(dt, posterior_estim,it){
   
    #define the HMM
    states = c("0","1")
    symbols = c("Obs0","Obs1","0","1")
    transProbs = array(dim = c(2,2,nrow(dt)), dimnames = list(states, states))
    emissionProbs = array(dim = c(2,4,nrow(dt)), dimnames = list(states, symbols))
    
    for (t in 1:nrow(dt)){
      transProbs[,,t] = matrix(c(1-dt$gamma0[t], dt$gamma0[t], 
                                 dt$gamma1[t], 1-dt$gamma1[t]), 
                               nrow = 2, byrow = T)
      emissionProbs[,,t] = matrix(c(1, 0, 1-dt$lambda[t], dt$lambda[t], 
                                    0, 1, 1, 0), byrow = T, nrow = 2)
    }
    
    hmm = varying_initHMM(
      States = states,
      Symbols = symbols,
      startProbs = c(1, 0),
      transProbs = transProbs,
      emissionProbs = emissionProbs)
    
    # forward filter simulation and correction
    y = dt$trace
    y[1] = "Obs0"
    f = varying_forward(hmm,y)
    
    forwardSeq = f
    for (t in 1:ncol(forwardSeq)){
      if (forwardSeq["1",t] == -Inf){forwardSeq["1",t]=0; forwardSeq["0",t]=1 
      } else if (forwardSeq["0",t] == -Inf){forwardSeq["0",t]=0; forwardSeq["1",t]=1
      } else if (forwardSeq["0",t] ==  Inf){forwardSeq["0",t]=1; forwardSeq["1",t]=0
      } else if (forwardSeq["1",t] ==  Inf){forwardSeq["1",t]=1; forwardSeq["0",t]=0}
      else {
        forwardSeq[,t] = exp(forwardSeq[,t] - mean(forwardSeq[,t]))
        forwardSeq[,t] = forwardSeq[,t] / sum(forwardSeq[,t])}
    }
    
    # backward simulation
    prob_loc = rep(NA, ncol(forwardSeq))
    simuX = rep(NA, ncol(forwardSeq))
    endP = forwardSeq[,ncol(forwardSeq)] / sum(forwardSeq[,ncol(forwardSeq)])
    simuX[ncol(forwardSeq)] = states[which(rmultinom(1, 1, prob = endP)==1)]
    for (t in 1:(ncol(forwardSeq)-1)){
      state = simuX[ncol(forwardSeq)-t+1]
      tMat = t(transProbs[,,ncol(forwardSeq)-t]) %*% diag(forwardSeq[,ncol(forwardSeq)-t])
      tMat = tMat / rowSums(tMat)
      simuX[ncol(forwardSeq)-t] = states[which(rmultinom(1, 1, prob = tMat[state,])==1)]
      #prob_loc[ncol(forwardSeq)-t] = tMat[[state]]
    }
    
    dt$prob_loc = prob_loc
    dt$traj = simuX
    dt$traj_lead = as.double(lead(simuX))
    
    if (posterior_estim){
      
      b = varying_backward(hmm, y)
      posteriorSeq = varying_posterior(f,b,y,length(states))
      
      avg.nbjumps = matrix(nrow = ncol(forwardSeq)-1, ncol=2)
      for (t in 1:(ncol(forwardSeq)-1)){
        tMat = t(transProbs[,,t]) %*% diag(forwardSeq[,t])
        tMat = tMat / rowSums(tMat)
        pMat = matrix(rep(posteriorSeq[,t+1],2), ncol=2, byrow=F)
        diag(pMat) = 0
        rMat = tMat * pMat
        avg.nbjumps[t,] = colSums(rMat)
      }
      nrow(dt)

      dt$avgT0 = c(avg.nbjumps[,1],0)
      dt$avgT1 = c(avg.nbjumps[,2],0)
      dt$avgA0 = posteriorSeq["0",]
      dt$avgA1 = posteriorSeq["1",]

    }
    dt$it = it
    return(dt)
  }
  

#### 1. Sampler ####
  
  
  logitll=function(beta,y,X,w){
    if (is.matrix(beta)==F) beta=as.matrix(t(beta)) 
    n=dim(beta)[1]
    pll=rep(0,n)
    for (i in 1:n){
      K = X %*% matrix(beta[i,])
      lF1 = plogis(K,log=T)
      lF2 = plogis(-K,log=T)
      pll[i]=sum(y*lF1*w+(1-y)*lF2*w)
    }
    return(pll) 
  }
  
hmflatlogit=function(iter_metropolis,y,X,scale,w,prior){
  p=dim(X)[2] 
  mod = summary(glm(y~-1+X,family=quasibinomial(link = "logit"), weights = w))
  beta=matrix(0,iter_metropolis,p) 
  beta[1,]=as.vector(mod$coeff[,1]) 
  Sigma2=as.matrix(mod$cov.unscaled)
  for (i in 2:iter_metropolis){ 
    tildebeta=rmvn(1,beta[i-1,],scale*Sigma2)
    
    priorM = rep(0,length(beta[i-1,]))
    priorSD = c(10,rep(4, length(beta[i-1])-1))
    
    if (prior==F){
      logpriorR = 0
    } else if (prior=="norm1"){
      priorM = rep(0,length(tildebeta))
      priorSD = c(5, rep(5, length(tildebeta)-1))
      logpriorR = sum(dnorm(tildebeta, mean=priorM, sd=priorSD, log = TRUE)) -
        sum(dnorm(beta[i-1,], mean=priorM, sd=priorSD, log = TRUE))
    }

    loglikR=logitll(tildebeta,y,X,w)-logitll(beta[i-1,],y,X,w)
    
    Mratio = logpriorR+loglikR
    
    if (log(runif(1))<=Mratio) {beta[i,]=tildebeta
    } else { beta[i,]=beta[i-1,] }
  } 
  return(beta)
}

hmflatlogit_wrapper = function(ww, iter_metropolis, prior=F){
  return(hmflatlogit(iter_metropolis,
                     y = ww$outcome,
                     X = model.matrix(ww$formulas, ww$dataPopMatRest),
                     scale = 1,
                     w = ww$dataPopMatRest$weight,
                     prior=prior)[iter_metropolis,])
}

#### 2. Clean and update functions ####

setting_up_data = function(x, groups,prop_of_observed,contr_simulated,
                  ena_pop,remove_early,NPER=65,model_index=NA){
  
  # Restricting the the field 
  y = rle(x$group)
  w = "NA"
  for (i in 1:length(y$lengths)){
    if (!(y$values[i] %in% groups)) y$values[i] = w
    w = y$values[i]
  }
  x$group = inverse.rle(y)
  x = x[group != "NA",]
  
  if (ena_pop) {
    x = x[!is.na(age_ena),]
    if (sum(x[age_ena >= 0.03125 & age_ena <= 0.09375,trace_obs])==0) return(x[NA,])
    if (model_index == "CV_empty"){x$trace_obs=0}
    x = x[age_ena >= 0.0625,]
    x = x[((period-age_ena)*(NPER-1))<=58,]
    x$age_ena = x$age_ena - 0.0625
    x$trace[x$age_ena==0] = "Obs0"
  }
  
  if (nrow(x)==0) return(x)
  
  # Crafting the covariates matrix
  x$age_min = x$period - min(x$period)
  trace_obs_l = as.double(x$trace_obs)
  trace_obs_l[1] = 0
  phi = 0.8
  if (model_index == "check_phi_09"){
    phi = 0.9
  }
  if (model_index == "check_phi_07"){
    phi = 0.7
  }
  if (model_index == "check_phi_06"){
    phi = 0.6
  }
  
  phiP = phi^-seq(1, length(trace_obs_l))
  x$lagged1 = lag(cumsum(trace_obs_l)/(1:length(trace_obs_l)))
  x$lagged2 = lag(cumsum(trace_obs_l*phiP)/cumsum(phiP))
  if (model_index == "check_phi_all"){
  x$lagged2b = lag(cumsum(trace_obs_l*0.6)/cumsum(0.6))
  x$lagged2c = lag(cumsum(trace_obs_l*0.7)/cumsum(0.7))
  x$lagged2d = lag(cumsum(trace_obs_l*0.9)/cumsum(0.9))
  }
  
  lagged3 = c() 
  for (l in transpose(rle(c(x$trace_obs)))){
    if (l[[2]]==0) lagged3 = append(lagged3,seq(1,l[[1]]))
    else lagged3 = append(lagged3, rep(0,l[[1]]))} 
  x$lagged3 = lag(lagged3, default = NA)
  
  # Parameters
  x$gamma0 = 0.05
  x$gamma1 = 0.05
  x$lambda = 0.1
  x$weight = 1
  x$trace = as.character(x$trace_obs)
  x[age_min==0]$trace_obs = 1
  x$lambda[x$age_min==0]=1
  
  # correcting depending on setting
  test_included = runif(1)<=prop_of_observed
  if (test_included){
    x$trace = case_when(x$traj_obs == "1"~"Obs1",
                        x$traj_obs == "0"~"Obs0",
                        T ~ x$trace)
  }
  
  x$profile_observed = any(x[!is.na(traj_obs)]$traj_obs=="1" | x[!is.na(traj_obs)]$traj_obs == "0" & x[!is.na(traj_obs)]$age_max > 0)
  x$profile_included = test_included & x$profile_observed
  
  x[traj_obs == 98,"weight"] = x[traj_obs == 98,"weight"]*contr_simulated
  
  if (remove_early) {
    x = x[((x$period - x$age_max)*64)>5,]
  }
  
  x$sex = factor(x$sex, levels = c("M","F"))
  x$prob_loc = rep(NA,nrow(x))
  
  # test gender
  return(x)
}

compute_transformed_params = function(x, mods){
  for (par in c("gamma0","gamma1","lambda")){
    x[[par]] = predict(mods[[par]], x, type = "response")
  }
  x$lambda[x$age_max==0]=1
  return(x)
}

#### 3. Varying HMM ####
  
varying_initHMM = function (States, Symbols, startProbs = NULL, transProbs = NULL, 
                              emissionProbs = NULL) 
{
    nStates = length(States)
    nSymbols = length(Symbols)
    S = rep(1/nStates, nStates)
    T = transProbs
    E = emissionProbs
    names(S) = States
    if (!is.null(startProbs)) {
      S[] = startProbs[]
    }
    return(list(States = States, Symbols = Symbols, startProbs = S, 
                transProbs = T, emissionProbs = E))
  }
  
library(HMM)
forward

varying_forward = function (hmm, observation) 
{
  
  hmm$transProbs[is.na(hmm$transProbs)] = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations = length(observation)
  nStates = length(hmm$States)
  f = array(NA, c(nStates, nObservations))
  dimnames(f) = list(states = hmm$States, index = 1:nObservations)
  for (state in hmm$States) {
    f[state, 1] = log(hmm$startProbs[state] * hmm$emissionProbs[state, 
                                                                observation[1],1])
  }
  for (k in 2:nObservations) {
    for (state in hmm$States) {
      logsum = -Inf
      for (previousState in hmm$States) {
        temp = f[previousState, k - 1] + log(hmm$transProbs[previousState, 
                                                            state, k - 1])
        if (temp > -Inf) {
          logsum = temp + log(1 + exp(logsum - temp))
        }
      }
      f[state, k] = log(hmm$emissionProbs[state, observation[k], k]) + 
        logsum
    }
  }
  return(f)
}

varying_backward = function (hmm, observation) 
{
  hmm$transProbs[is.na(hmm$transProbs)] = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations = length(observation)
  nStates = length(hmm$States)
  b = array(NA, c(nStates, nObservations))
  dimnames(b) = list(states = hmm$States, index = 1:nObservations)
  for (state in hmm$States) {
    b[state, nObservations] = log(1)
  }
  for (k in (nObservations - 1):1) {
    for (state in hmm$States) {
      logsum = -Inf
      for (nextState in hmm$States) {
        temp = b[nextState, k + 1] + log(hmm$transProbs[state, nextState, k+1] * hmm$emissionProbs[nextState, observation[k + 1], k+1])
        if (temp > -Inf) {
          logsum = temp + log(1 + exp(logsum - temp))
        }
      }
      b[state, k] = logsum
    }
  }
  return(b)
}
  
  
varying_posterior = function (f, b, observation, nstates) 
{
  probObservations = f[1, length(observation)]
  for (i in 2:nstates) {
    j = f[i, length(observation)]
    if (j > -Inf) {
      probObservations = j + log(1 + exp(probObservations - 
                                           j))
    }
  }
  posteriorProb = exp((f + b) - probObservations)
  posteriorProb2 = posteriorProb
  
  posteriorProb2["0",] = posteriorProb2["0",]/colSums(posteriorProb)
  posteriorProb2["1",] = posteriorProb2["1",]/colSums(posteriorProb)
  return(posteriorProb2)
}
  
  #### 4. Wrapper inference functions ####
  
inference_mod = function(formulas, groups, dataTot, model_index = NA,
                         N = 1000, size_simu = 1, iter_metropolis = 30, obsburn = 900, 
                         prop_of_observed = 1, contr_simulated = 1,
                         ena_pop = F, remove_early = F){

  # Data selection and preparation
  dataTot2 = parLapply(cl, dataTot, setting_up_data, NPER=65,
                       groups = groups, ena_pop = ena_pop, remove_early = remove_early,
                       prop_of_observed = prop_of_observed, contr_simulated = contr_simulated, 
                       model_index = model_index)
  
  
  dataPop = list()
  for (x in names(dataTot)) if (nrow(dataTot2[[x]])>1) dataPop[[x]] = dataTot2[[x]]
  
  # quantities
  mods = list()
  save_params = list()
  save_obs = list()

  # first pass to initialize
  dataPop = parLapply(cl, dataPop, simulation_mod, posterior_estim = T, it = 0)
  dataPopMat = rbindlist(dataPop)
  #list_var = c("ident","group","traj_obs","lagged1","lagged2","lagged3","period","age_min","age_ena",
  #             "age_max","profile_present","trace","trace_obs","traj","traj_lead")
  given_data = dataPopMat
  
  vec_field = list("gamma0"=dataPopMat$traj==0 & !is.na(dataPopMat$traj_lead),
                   "gamma1"=dataPopMat$traj==1 & !is.na(dataPopMat$traj_lead),
                   "lambda"=dataPopMat$traj==0 & !is.na(dataPopMat$lagged1))

  for (par in c("gamma0","gamma1","lambda")){
    # Estimating so to have a model object with the right structure
    mods[[par]] = glm(formula = formulas[[par]], 
                      data = dataPopMat[vec_field[[par]],], 
                      family = binomial(link = "logit"))
    
    # Creating an array of the right format for saving parameters
    save_params[[par]] = array(dim = list(length(mods[[par]]$coefficients),N),
                               dimnames = list(names(mods[[par]]$coefficients)))
    
  }
  
  print(0)
  t = 1
  
  prior = F
  if (model_index == "prior1"){
    prior = "norm1"
  } else if (model_index == "prior2"){
    prior = "norm2"
  }
  print(prior)
  
  for (t in 1:N){
    if (t%%10 ==0) print(t)
    
    ## Gibbs step
    
    # Computing transformed parameters 
    dataPopMat = rbindlist(dataPop)
    for (par in c("gamma0","gamma1","lambda")) dataPopMat[[par]] = predict(mods[[par]], dataPopMat, type = "response")
    dataPopMat$lambda[dataPopMat$age_min==0]=1
    dataPop = split(dataPopMat,by="ident")
 
    # Simulating data
    subsel = sample(names(dataPop), round(size_simu*length(dataPop)), replace=F)
    dataPop[subsel] = parLapply(cl, dataPop[subsel], simulation_mod, 
                                posterior_estim = t>obsburn, it = t)
     
    if (t>obsburn) save_obs[[t]] = rbindlist(dataPop[subsel])[,c("ident","it","age_max","traj","avgA0","prob_loc")]
    
    ## Metropolis-Hastings step
    
    # Matrices for computation and data organization
    dataPopMat = rbindlist(dataPop,fill=T)
    vec_field = list("gamma0"=dataPopMat$traj==0 & !is.na(dataPopMat$traj_lead),
                     "gamma1"=dataPopMat$traj==1 & !is.na(dataPopMat$traj_lead),
                     "lambda"=dataPopMat$traj==0 & !is.na(dataPopMat$lagged1))
    
    metropolis_data = list()
    metropolis_data$formulas = formulas
    metropolis_data$outcome = list("gamma0"=dataPopMat[vec_field[["gamma0"]],traj_lead],
                   "gamma1"=1-dataPopMat[vec_field[["gamma1"]],traj_lead],
                   "lambda"=dataPopMat[vec_field[["lambda"]],trace_obs])
    
    metropolis_data$dataPopMatRest = list()
    for (par in c("gamma0","gamma1","lambda")) metropolis_data$dataPopMatRest[[par]] = dataPopMat[vec_field[[par]]]
    
    # Posterior sampling 
    
    metropolis_result = parLapply(cl, purrr::transpose(metropolis_data), hmflatlogit_wrapper,
                                  iter_metropolis = iter_metropolis, prior=prior)
  
    # Updating parameters within models
    for (par in c("gamma0","gamma1","lambda")){
      mods[[par]]$coefficients = metropolis_result[[par]]
      save_params[[par]][,t] = metropolis_result[[par]]
    }

    }
  
    return(list("params"=save_params, "obs"=rbindlist(save_obs), "given_data"=given_data))
}
  
#### 5. Miscellaneous #### 


mkdate = function(y,m) {(((y - 1990)*12 + (m-1)) %/% 6) +1}
