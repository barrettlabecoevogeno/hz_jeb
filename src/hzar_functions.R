# HZAR FUNCTIONS -----------------------------------------------
# Wrappers around HZAR functions, for simulation.


fitNoneHZAR <- function(transect_df, c, w) {
  # This function takes in a data frame containing the information about
  # a transect, and returns a an hzar list object
  # with all the results of running hzar with the free.none model
  
  # It mostly follows along with the HZAR example script, with some modifications
  
  # Prep HZAR results object ------------------------------------------------
  locus <- list();
  ## Space to hold the observed data
  locus$obs <- list();
  ## Space to hold the models to fit
  locus$models <- list();
  ## Space to hold the compiled fit requests
  locus$fitRs <- list();
  ## Space to hold the output data chains
  locus$runs <- list();
  ## Space to hold the analysed data
  locus$analysis <- list();
  
  
  # Add observed data, from the transect_df data frame ----------------------
  locus$obs <- hzar.doMolecularData1DPops(transect_df$transect.dist, 
                                          transect_df$Cr.HYD, 
                                          transect_df$Cr.nSamples)
  
  print("HZAR results object prepped")
  # Specify the model to test ----------------------------------------------
  locus$models[["free.none"]] <- hzar.makeCline1DFreq(locus$obs,  scaling = "free", tails = "none")
  
  
  # Compile model ----------------------------------------------------------
  ## Make the first hzar.fitRequest object for each of the models
  locus$fitRs$init <- sapply(locus$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=locus$obs,
                             verbose=FALSE,
                             simplify=FALSE)
  locus$fitRs$init$free.none$modelParam$init$center <- rnorm(n = 1, mean = c, sd = 20)
  locus$fitRs$init$free.none$modelParam$init$width <- abs(rnorm(n = 1, mean = w, sd = 15))
  locus$fitRs$init$free.none$modelParam$init$pMin <- runif(n = 1, 0, 0.2)
  locus$fitRs$init$free.none$modelParam$init$pMax <- runif(n = 1, 0.8, 1)
  print("Models compiled")
  # Initialize chains -------------------------------------------------------
  locus$runs$init <- list() 
  # This takes a while, ~30 minutes, 
  # right now it runs each of the chains in serial
  # might be a way to speed this up with parallelization.
  for (model in names(locus$fitRs$init)) {
    model.id <- noquote(model)
    locus$runs$init[[model.id]] <- hzar.doFit(locus$fitRs$init[[model.id]])
  }
  print("Chains initialized")
  # Compile new fit requests from the initial chains ------------------------
  locus$fitRs$chains <-
    lapply(locus$runs$init, hzar.next.fitRequest)
  
  ## Replicate each fit request 3 times, keeping the original
  ## seeds while switching to a new seed channel.
  locus$fitRs$chains <-
    hzar.multiFitRequest(locus$fitRs$chains,
                         each=3,
                         baseSeed=NULL)
  print("New fit requests compiled")
  # Randomize seeds  --------------------------------------------------------
  # Randomize centers
  # Center is in all models, so in all 3 chains
  center.seeds <- rnorm(n = 3, mean = c, sd = 20)
  for (chain in 1:3) {
    locus$fitRs$chains[[chain]]$modelParam$init["center"] <- center.seeds[chain]
  }
  # Randomize widths
  # also in all 15 chains
  width.seeds <- abs(rnorm(n = 3, mean = w, sd = 15))
  for (chain in 1:3) {
    locus$fitRs$chains[[chain]]$modelParam$init["width"] <- width.seeds[chain]
  }
  #Randomize pMin and pMax, only in models where interval is free
  pmin.seeds <- runif(n = 3, 0, 0.2)
  pmax.seeds <- runif(n = 3, 0.8, 1)
  for (chain in 1:3) {
    locus$fitRs$chains[[chain]]$modelParam$init["pMin"] <- pmin.seeds[chain]
    locus$fitRs$chains[[chain]]$modelParam$init["pMax"] <- pmax.seeds[chain]
  }
  
  print("Seeds randomized")
  # Run the main analysis chains --------------------------------------------
  ## This is definitely gonna take a while, even in parallel. 
  locus$runs$chains <-  hzar.doChain.multi(locus$fitRs$chains,
                                           doPar=T,
                                           inOrder=FALSE,
                                           count=3)
  print("Analysis chains finished")
  # Compile chains into results objects -------------------------------------
  
  ## Create a model data group for the null model (expected allele
  ## frequency independent of distance along cline) to include in
  ## analysis.
  locus$analysis$initDGs <- list(
    nullModel =  hzar.dataGroup.null(locus$obs))
  
  ## Create a model data group (hzar.dataGroup object) for each
  ## model from the initial runs.
  for (model in names(locus$runs$init)) {
    model.id <- noquote(model)
    locus$analysis$initDGs[[model.id]] <- hzar.dataGroup.add(locus$runs$init[[model.id]])
  }
  
  
  ## Create a hzar.obsDataGroup object from the four hzar.dataGroup
  ## just created, copying the naming scheme.
  locus$analysis$oDG <-
    hzar.make.obsDataGroup(locus$analysis$initDGs)
  locus$analysis$oDG <-
    hzar.copyModelLabels(locus$analysis$initDGs,
                         locus$analysis$oDG)
  
  ## Convert all runs to hzar.dataGroup objects, adding them to
  ## the hzar.obsDataGroup object.
  locus$analysis$oDG <-
    hzar.make.obsDataGroup(lapply(locus$runs$chains,
                                  hzar.dataGroup.add),
                           locus$analysis$oDG)
  print("Results written to object")
  return(locus)
}
# HZAR model selection ----------------------------------------------------
modelSelect <- function(cline_list) {
  # This function takes in a list of cline results from HZAR,
  # Adds an AICc table to the results,
  # and extracts the hzar.dataGroup object for the model with the lowest AICc
  
  
  # Add AICc table to results
  cline_list$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(cline_list$analysis$oDG)
  # Print AICc table
  print(cline_list$analysis$AICcTable)
  
  ## Print out the model with the minimum AICc score
  print(cline_list$analysis$model.name <-
          rownames(cline_list$analysis$AICcTable
          )[[ which.min(cline_list$analysis$AICcTable$AICc )]])
  
  ## Extract the hzar.dataGroup object for the selected model
  cline_list$analysis$model.selected <-
    cline_list$analysis$oDG$data.groups[[cline_list$analysis$model.name]]
  
  return(cline_list)
}

# Extract model parameters from HZAR --------------------------------------
extractModelParams <- function(cline_list) {
  # This function takes in an HZAR object containing the results
  # of the cline fitting, and returns a data frame of all the parameter
  # estimates from each model, and their upper and lower 2LL limits
  
  # Get the models from the cline_list object, dropping the first one (the NULL model)
  models <- names(cline_list$analysis$oDG$data.groups)[-1]
  
  # start the results data frame, which will hold everything. Not the most efficient,
  # but fine for our purposes. 
  res <- data.frame(model = "filter", parameter = "filter", estimate = 0,
                    lower = 0, upper = 0, row.names = F)
  
  # for each model
  for (model in models) {
    #extract the names of the estimates parameters in the model
    name.Ps <- names(cline_list$analysis$oDG$data.groups[[model]]$ML.cline$param.free)
    # extract the values of those parameters
    Ps <- t(cline_list$analysis$oDG$data.groups[[model]]$ML.cline$param.free)
    # extract the upper and lower limits of the parameters
    limits <- t(hzar.getLLCutParam(cline_list$analysis$oDG$data.groups[[model]], name.Ps))
    
    # create a mini results data frame of the results of 1 model
    mod.res <- data.frame(rep(model, times = length(Ps)),
                          name.Ps, estimate = Ps,
                          limits[which(1:length(limits) %% 2 == 1)], # lower limits are the odd ones
                          limits[which(1:length(limits) %% 2 != 1)], row.names = NULL) # upper limits even
    names(mod.res) <- c("model", "parameter", "estimate", "lower", "upper")
    res <- rbind(res, mod.res) # add the individual model's results to the overall result df
  }
  res <- filter(res, model != "filter") # get rid of the extra row we used to initialize
  
  # To facilitate comparison later, add extract AIC information from the cline list
  aictable <- cbind(model = row.names(cline_list$analysis$AICcTable), cline_list$analysis$AICcTable) %>%
    arrange(AICc) %>%
    mutate(deltaAICc = AICc - min(AICc)) %>%
    filter(model != "nullModel")
  
  # Add the AICc results to the parameter estimates
  res <- left_join(res, aictable, by = "model")
  return(res)
}


# Cline equation functions ------------------------------------------------

# Functions implementing the cline equations from the models in Stan
# As arguments, they take in the distance along the transect and 
# the necessary parameters.
# They return a value of expected allele frequency at that location
# along the transect. 

# Using the R Vectorize() function, I've made vectorized versions of all of these,
# to ease plotting, the first argument can be a vector of distances along the transect,
# and it will return a vector of expected allele frequencies.

# All equations were copied from the Stan model code and 
# modified to work in R.

# No tail functions -------------------------------------------------------
none.eqn <- function(transectDist, center, width, pmin, pmax) {
  p <- pmin + (pmax - pmin) * (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
  return(p)
}

none.eqn.vec <- Vectorize(none.eqn)


# Left tail functions -----------------------------------------------------
left.eqn <- function(transectDist, center, width, pmin, pmax, deltaL, tauL) {
  if (transectDist < center - deltaL) {
    p <-  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist - center + deltaL)/width)/(1 + exp(-4*deltaL/width)))
  } 
  else {
    p <- pmin + (pmax - pmin) * (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
  }
  return(p)
}
left.eqn.vec <- Vectorize(left.eqn)

# Right tail functions ----------------------------------------------------
right.eqn <- function(transectDist, center, width, pmin, pmax, deltaR, tauR) {
  if (transectDist >= center + deltaR) {
    p <-  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist - center - deltaR)/width)/(1 + exp(-4*deltaR/width))))
  } 
  else {
    p <- pmin + (pmax - pmin) * (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
  }
  return(p)
}

right.eqn.vec <- Vectorize(right.eqn)

# Mirror tail functions ---------------------------------------------------

mirror.eqn <- function(transectDist, center, width, pmin, pmax, deltaM, tauM) {
  if (transectDist <= center - deltaM) { # left
    p <-  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaM/width)))*exp((4*tauM*(transectDist - center + deltaM)/width)/(1 + exp(-4*deltaM/width)))
  }
  else if(transectDist >= center + deltaM) { #right
    p <-  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaM/width)))*exp((-4*tauM*(transectDist - center - deltaM)/width)/(1 + exp(-4*deltaM/width))))
  } 
  else {
    p <- pmin + (pmax - pmin) * (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
  }
  return(p)
}

mirror.eqn.vec <- Vectorize(mirror.eqn)

# Independent tail functions ----------------------------------------------
ind.eqn <- function(transectDist, center, width, pmin, pmax, deltaR, tauR, deltaL, tauL) {
  if (transectDist <= center - deltaL) { # left
    p <-  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist - center + deltaL)/width)/(1 + exp(-4*deltaL/width)))
  }
  else if(transectDist >= center + deltaR) { #right
    p <-  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist - center - deltaR)/width)/(1 + exp(-4*deltaR/width))))
  } 
  else {
    p <- pmin + (pmax - pmin) * (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
  }
  return(p)
}

ind.eqn.vec <- Vectorize(ind.eqn)


