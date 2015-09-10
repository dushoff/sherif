
### Wrapping function for 'synlik'
###
sherif_wrap <- function(param, nsim, extraArgs, ...){
  
  # Unpack
  paramsModel <- extraArgs[["paramsModel"]]
  paramsSimul <- extraArgs[["paramsSimul"]]
  
  # NOTE: 'param' = model parameters to be calibrated
  # so must override model parameter given in 'paramsModel'
  for(nn in names(param)) paramsModel[nn] <- param[nn]
  
  # Overrides mc_iter by nsim
  paramsSimul[["mc_iter"]] <- nsim
  
  # Overrides the seed for random number generator
  if(!is.null(seed)) paramsSimul[["seed"]] <- as.integer(seed)
  
  # Run the simulation with appropriate parameters
  sim <- rcpp_sherif(paramsSimul=paramsSimul,
                     paramsModel=paramsModel)
  
  # Store the 'nsim' simulations into a list of data frames:
  
  # data frame of time series:
  ts <- list()
  # Generation intervals
  GIbck <- list()
  
  for(i in 1:nsim){
    df <- data.frame(time= sim$time[[i]],
                     cumInc= sim$cumIncidence[[i]],
                     cumInc_hcw= sim$cumIncidence_hcw[[i]],
                     deathsCum= sim$deathsCum[[i]],
                     deathsCum_hcw= sim$deathsCum_hcw[[i]])
    ts[[i]] <- df
    GIbck[[i]] <- sim$GIbck_gi[[i]]
  }
  return(list(ts=ts, GIbck=GIbck))
}




