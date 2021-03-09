library(spmoran)
library(doParallel)
library(maptools)
library(raster)
library(rgeos)
library(vegan)
library(gstat)
library(igraph)
library(phylin)
library(reldist)
library(vegan)
library(betapart)
library(bipartite)
library(spdep)
library(spacemakeR)
library(adespatial)
library(reshape2)
library(zoo)

#precise whether parallelisation is used
isParallelised <- TRUE

#set working directory
setwd("~/Desktop/PhD/16_03_20_workFromHome/CleanUp_17_03_2020/Model/Model_Current/Model/newWorkingDirectory")


##### MODEL FUNCTIONS ######
#General view of functions structure:

#1) "createLandscape" prepares the baseLandscape object and launches createSpecies.
#2) "populateLandscapeWithSpecies" prepares and integrates species data into a baseLandscape object.
#3) "launchModels" first launches createLandscapeVariations,
# "createLandscapeVariations" creates tables of landscape variations for all times steps, 
# for all combinations of heterogeneity, frequency and magnitude levels.
#4) "launchModels" then launches the  "model" function several times, to fulfill all possible combinations of heterogeneity, frequency and magnitude levels.
#5) the model function simulates a meta-community over a specified timescale for a specific combination of heterogeneity, frequency and magnitude levels.

#6) several functions calculate the niche-neutral metrics for the models' results with : 
#Variation Partitioning, Abundance Beta-Null Deviation, Species Variation Across Replicates (SVAR), RAD shape, and clustering


#NNC MODEL FUNCTION
#This is the main function. The simulation of a species meta-community following a modification of Gravel's original stochastic niche model (Gravel et al. 2006).
# Simulates the dynamics of multiple species in a patchy landscape, with each patch having an environmental condition E {0,100}
# requires the following input parameters: 
# [1] z = exponential distance decay parameter that denotes the shape of the dispersal kernel
# [2] J = total number of sites across the landscape (across all patches combined) [default 1000]
# [3] MaxPatches = total number of patches across the landscape [default 50]
# [4] Nspecies = initial number of species [default 30]
# [5] Tmax = number of time steps [default 1000] { Notice that such a long period leads to the extinction of most species in the neutral and continuum cases}
# [6] d = mortality rate (fraction of individuals lost in each time step) [default 0.25]
# [7] filt = a string denoting the type of filter for community generation ('niche' [default], 'neutral', or 'continuum')
# [8] land = landscape list
# [-] storageReduction : IGNORED 
# [9] Niterations = total iterations
#hidden
# [10] iteration = not to be precised, program does it automatically as a recursive function
# [11] results = not to be precised, program passes results on through recursive levels

#ignored in this project
#abundances, spRichness, temporalList, abundanceStart, temporalSampling
model = function(z, 
                 J = J, 
                 MaxPatches = MaxPatches, 
                 fitnessVariation, 
                 Tmax = Tmax, 
                 d, 
                 filt = "niche", 
                 land, 
                 storageReduction,
                 Niterations = 1, 
                 iteration = 0, 
                 results = NULL,
                 
                 abundances = NULL,
                 spRichness = NULL,
                 temporalList = NULL,
                 abundanceStart = NULL,
                 temporalSampling = 0
)
{
  
  #initialise iterations
  if(!iteration){iteration <- Niterations}
  
  # initialize simulation variables
  Erichness = c()
  Srichness = c()
  
  #initialise dataframe for  results
  newResults <- as.data.frame(array(data = 0, dim = c(1, 8)))
  
  # if a landscape list is included, import landscape characteristics from it
  print("gettingLand")
  if(!is.null(land)){
    Npatches = land$Npatches
    Jpatches = land$Jpatches
    Spatches = land$Spatches
    EMpatches = land$EMpatches
    D = land$D
    
    #use another abundanceStart if it is provided, otherwise use landscape one
    #abundanceStart is IGNORED in this project
    if(!is.null(abundanceStart)){
      #if abundanceStart has iterations, use current one
      if(length(abundanceStart) == Niterations){
        #associate each abundance start to the relevant model
        abundance = abundanceStart[[iteration]]
      }else{
        abundance = abundanceStart
      }
    }else{
      abundance = land$speciesInfo$abundanceStart
    }
    abundanceStart0 = abundance
    storage = land$speciesInfo$storageStart
  }
  
  # Calculate heterogeneity metrics (Ignored in this project)
  tempHetero <- 0 	
  
  #create list to host species richness and evenness through time
  spRich <- c()
  spEvenness <- c()
  
  #Ti is equal to 1 or t, depending on whether patches vary through time (dynamic model) or not (static model.
  Ti <- 1 
  startShift <- 1
  
  #if vector given for tmax, use current iteration for Tmax
  if(length(Tmax) > 1){
    #prolong all iterations by 5 timesteps, but landscapevariations are randomly shifted
    startShift <- Tmax[[iteration]]
    finalTmax <- startShift + 5
  }else{
    finalTmax <- Tmax
  }
  # Simulate meta-community dynamics
  for (t in startShift:finalTmax)
  {
    
    abundance0 <- abundance # keep the abundance levels from the beginning of the current time step
    
    #storage IGNORED in this projected
    storage <- floor(storage * storageReduction) # apply storage death rate
    
    #if there are temporal variations, Ti == t. If there are no variations (static model), Ti remains 1.
    #Ti will only reach 200 and repeat itself afterwards (limits RAM used for environmental variation)
    if(length(fitnessVariation) > 1 ){Ti <- t%%200 }
    #Ti goes from 1 to 200, so if %% returns 0, make it 200.
    if(Ti == 0){Ti <- 200}
    
    if(Ti %% 100 == 0){
      print(t)
      print(Ti)
      print("---")
    }
    
    for (p in 1:Npatches){
      #apply death rate and get number of available sites
      survivors = floor((1 - d) * abundance0[p, ])	# abundance of surviving individuals (post mortality; first right-hand term in manuscript equation [1])
      available = Jpatches[p] - sum(survivors)		# number of available sites
      
      #debug
      if(available < 1){browser()}
      
      
      
      
      #main model, decomposed into new indivduals from storage, non-storage, total storage and total non-storage
      nonStor = 0
      nonStorTot = 0
      stor = 0
      storTot = 0
      
      if (filt == 'niche'){ 
        
        if(storageReduction != 0){
          #storageReduction ignored in this project
          stor = fitnessVariation[[Ti]][p, ]*storage[p,]
          storTot = sum(fitnessVariation[[Ti]][p, ]*storage[p,])
        }
        
        nonStor = fitnessVariation[[Ti]][p, ]* colSums(D[p, ] * (abundance0 > 0))
        nonStorTot = sum(fitnessVariation[[Ti]][p, ]* colSums(D[p, ] * (abundance0 > 0)))
        
      }
      
      if (filt == 'continuum') {
        
        #storageReduction ignored in this project
        if(storageReduction != 0){
          stor = fitnessVariation[[Ti]][p, ]*storage[p,]
          storTot = sum(fitnessVariation[[Ti]][p, ]*storage[p,])
        }
        
        nonStor = fitnessVariation[[Ti]][p, ] * colSums(D[p, ] * (abundance0))
        nonStorTot = sum(fitnessVariation[[Ti]][p, ] * colSums(D[p, ] * (abundance0)))
        
      }
      
      if (filt == 'neutral'){
        
        #storageReduction ignored in this project
        if(storageReduction != 0){
          stor = storage[p,]
          storTot = sum(storage[p,])
        }
        
        nonStor = colSums(D[p, ] * (abundance0) ) 
        nonStorTot =  sum(colSums( D[p, ]*(abundance0) ) )
        
      }
      
      ## Calculate number of establishing individuals of each species
      
      #intialise recruits
      recruits = 0
      recruitsFromStorage = 0
      
      #calculate probability of settling of each species
      #debug
      if(is.nan(((nonStor + stor)/(nonStorTot + storTot))[1])){browser()}
      # probability of establishment per species
      # storage is generally IGNORED, thus stor and storTot remain at 0
      # m gives immigration from implicit metacommunity
      m <- 0.1
      # P is the relative abundance for every species in the metacommunity, we use a uniform distribution
      P <- rep(1/Nspecies, Nspecies)
      R <- ( (1-m) * ( (nonStor + stor)/(nonStorTot + storTot) ) ) + (m * P)
      
      if (max(R) > 0){
        recruits =  rmultinom(1, available, R)
        recruitsFromStorage = floor(recruits * ((stor/(storTot+nonStorTot)) / ((stor+nonStor)/(storTot+nonStorTot))) ) #get amount of recruits that are from storage
      }	
      
      #apply decay and recruitement to storage
      #IGNORED as storage kept at 0
      storSurvivors = floor((1-storageReduction)*storage[p, ])
      storAvailable <- (Spatches[p] - sum(storage[p,]) ) * storageReduction #number of available storage sites in p
      storAvailable[storAvailable <0] <- 0
      #calculate probability of adding to storage of each species
      S <- ( colSums(D[p, ] * (abundance0) ) ) / sum(colSums( D[p, ]*(abundance0) ) )
      #debug
      if (!is.finite(S)){
        browser()
      }
      if (max(S) > 0){
        storRecruits <- rmultinom(1, storAvailable, S)
      }
      
      #abundances updated 
      abundance[p, ] = survivors + recruits		# update species abundance distribution per patch
      
      #IGNORED due to storage remaining at 0
      storage[p, ] <- storage[p, ] - recruitsFromStorage #remove recruited individuals from storage
      storage[p, ] = storSurvivors + storRecruits #storage recruits added to abundance
      
      #debug
      if (sum(abundance[p,]) > Jpatches[p]){browser()}
    }
    
    #update species richness
    spRich <- c( spRich, sum(colSums(abundance)>0) )
    
    #update community evenness score
    sortedAbundance <- sort(colSums(abundance), decreasing = TRUE)
    sortedAbundance <- sortedAbundance[sortedAbundance>0]
    
    #use evenness formula
    #S = number of species
    S <- length(sortedAbundance)
    #p = proportions of each species
    p <- sortedAbundance/sum(sortedAbundance)
    evennessScore <- 1 - ( ( (S*sum(p^2) - 1) / (S-1) ) )^(1/2)
    spEvenness <- c(spEvenness, evennessScore)
    
    
    
    #temporalSampling IGNORED in this project
    if(temporalSampling >0){
      if ( t %% temporalSampling == 0){
        #add current abundances to list
        if(!exists("temporalAbundances")){
          temporalAbundances <- list(abundance)
        }else{
          temporalAbundances <- c(temporalAbundances, list(abundance))
        }
      }
    }
    
    #end of time cycle
  }
  
  #create results at end of time cycle
  if(temporalSampling == 0){
    Srichness = c(Srichness, sum(colSums(abundance) > 0, na.rm = TRUE))	# calculate the number of species at the landscape scale
    newResults <- c(iteration, Npatches, Nspecies, Srichness, z, tempHetero)
    
    #create results if first run
    if(is.null(results)){
      results <- newResults
    }else{
      #otherwise appends to results passed through recursion
      results <- rbind(results, newResults) 
      colnames(results) <- c("iteration", "Npatches", "Nspecies", "Srichness", "z", "Shannon")
    }
    
    #create abundances if first run
    if(is.null(abundances)){
      abundances <- list(abundance)
    }else{
      #or append to existing abundance list
      abundances <- c(abundances, list(abundance))
    }
    
    #create species richness list if first run
    if(is.null(spRichness)){
      spRichness <- list(spRich)
    }else{
      #or append to existing species richnesses
      spRichness <- c(spRichness, list(spRich))
    }
    
  }else{
    #Temporal sampling IGNORED in this project
    #if there is temporal sampling, make iterative list of temporal list (gets confusing...)
    if(is.null(temporalList)){
      temporalList <- list(temporalAbundances)
    }else{
      temporalList <- c(temporalList, list(temporalAbundances))
    }
  }
  
  #replication
  #replicates are handled through a recursion
  #results are passed to each recursion
  iteration <- iteration - 1
  if(iteration < 1){
    if(temporalSampling >0){
      return(list(list(temporalList)))
    }else{
      return(list(list(res = results, abu = abundances, rich = spRichness, evenness = spEvenness) ) )
    }
  }else{
    #continue recursion through all iterations
    return(model(z = z,
                 J = J, 
                 MaxPatches =  MaxPatches,
                 fitness =  fitnessVariation, 
                 Tmax =  Tmax, 
                 d = d, 
                 filt =  filt, 
                 land = land, 
                 Niterations = Niterations,
                 storageReduction = storageReduction,
                 # nicheAbuDependant = nicheAbuDependant,
                 # nicheDispersalLimitation = nicheDispersalLimitation,
                 iteration =  iteration, 
                 results = results, 
                 abundances = abundances,
                 spRichness = spRichness,
                 temporalList = temporalList,
                 abundanceStart = abundanceStart0,
                 temporalSampling = temporalSampling))
  }
}


#nullAnalysis compares null assembly models with results of different patches
#CREATE LANDSCAPE FUNCTION
#creates basic landscape with no environmental variation
#takes: Npatches, J, z 
#returns: baseLandscape (a list containing: Npatches, Jpatches, EMpatches, xy, D, Spatches)
createLandscape <- function(Npatches, 
                            J, 
                            z){
  #N = patches
  #J = cells within patches
  #S = storage cells for each patch
  
  Jpatches = rmultinom(1, J, rep(1, Npatches))  # randomize number of sites per patch (all will sum up to J)
  
  XY = matrix(runif(Npatches * 2, 1, size), nrow = Npatches, ncol = 2) # generate random patch locations
  Ddist = as.matrix(dist(XY, diag = T, upper = T))	                  # calculate inter-patch distances
  D = exp(-z * Ddist)		                                       # calculate propagule arrival ratios based on an exponential dispersal kernel. z is the distance-decay parameter
  
  #### Generate autocorrelated grid of E values
  #### Extract E values based on N patch positions.
  xy = expand.grid(1:size, 1:size)
  names(xy) = c("x","y")
  gridded(xy) = ~x+y
  g.dummy = gstat(formula = z~1, dummy = TRUE, beta = 50, model = vgm(100,"Exp",10), nmax = 10)
  yy = predict(g.dummy, xy, nsim = 1)
  ry = raster(yy)
  EMpatches = extract(ry, XY) # mean patch E
  
  
  #save landscape as a list
  landscape = list(Npatches = Npatches, Jpatches = Jpatches, EMpatches = EMpatches, xy = XY, D = D, z = z)
  
  return(landscape)
}

#Storage IGNORED in this project
addStorageToLandscape <- function( baseLandscape, storageMultiple){
  
  Spatches <- baseLandscape$Jpatches * storageMultiple
  baseLandscape[[length(baseLandscape)+1]] <- Spatches
  listNames <- names(baseLandscape)
  listNames[[length(listNames)]] <- "Spatches"
  names(baseLandscape) <- listNames
  
  return(baseLandscape)
}

#POPULATE LANDSCAPE FUNCTION
#generates starting species communities in landscape
#takes: base landscape and Nspecies 
#returns: landscape with species data
populateLandscapeWithSpecies <- function(landscapeWithStorage, Nspecies){
  
  
  # GENERATE SPECIES
  #EMspecies = sample(EMpatches, Nspecies, replace = TRUE)  	
  
  # Assign a niche center for each species (niche center evenly distributed along niche axis), Env1
  Const <- 40
  #distribute niches equally
  EMspecies <- ((Const)/2) + (1:Nspecies *  ((100-Const) / (Nspecies)) )
  #ditribute species randomly
  #EMspecies <- rnorm(Nspecies, 50, 10)
  
  Npatches <- landscapeWithStorage$Npatches
  EMpatches <- landscapeWithStorage$EMpatches
  #generate generalist/specialist species, by varying their sigmaSpecies (niche breadth)
  sigmaSpecies = runif(Nspecies, specialist, generalist) # randomise specialist-generalist
  
  abundanceStart = array(data = 0, dim = c(Npatches, Nspecies))	# initialize the species abundance table (patch by species array)
  abundance.s = abundanceStart#useful?
  storageStart <- array(data = 0, dim = c(Npatches, Nspecies))
  
  EMfitness = array(data = 0.0, dim = c(Npatches, Nspecies))	# initialize the array of habitat suitability per species per patch
  
  #determine initial conditions (abundances and starting (mean) fitness)
  for (p in 1:Npatches){
    for (S in 1:Nspecies){
      EMfitness[p,S] = exp(-(EMspecies[S] - EMpatches[p])^2 / (2 * sigmaSpecies[S]^2))*(specialist/sigmaSpecies[S])	# calculate habitat suitability per species per patch, based on manuscript equation [2]
    }
    
    #distribute abundances equally between patches
    #abundanceStart[p, ] = distributeSum( landscapeWithStorage$Jpatches[[p]], Nspecies )#fitness[p, ])		# determine initial species abundances in patches
    #storageStart[p, ] = distributeSum( landscapeWithStorage$Spatches[[p]], Nspecies) #determine initial storage
    
    #distribute abundances randomly
    abundanceStart[p, ] = randomDistr( landscapeWithStorage$Jpatches[[p]], Nspecies )#fitness[p, ])		# determine initial species abundances in patches
    
    #R <- EMfitness[p, ]
    #abundanceStart[p, ] <- rmultinom(1, landscapeWithStorage$Jpatches[[p]] , R )
    
    storageStart[p, ] = randomDistr( landscapeWithStorage$Spatches[[p]], Nspecies) #determine initial storage
    #storageStart[p, ] <- rmultinom(1, landscapeWithStorage$Spatches[[p]] , R )
    
  }
  
  
  #get total abundances and correct for inequality
  goal <- J/Nspecies
  for(n in 1:Nspecies){
    diff <- goal - sum(abundanceStart[, n])
    if(diff != 0){
      #remove or add difference from throughout patches randomly
      distrDiff <- rmultinom(1, abs(diff), rep(abs(diff)/Npatches, Npatches))
      abundanceStart[, n] <- abundanceStart[, n] + (sign(diff)*distrDiff)
      for (p in 1:Npatches){
        if (abundanceStart[[p, n]] < 0){
          abundanceStart[[which(abundanceStart[, n] > -abundanceStart[[p, n]])[[1]], n]] <- abundanceStart[[which(abundanceStart[, n] > -abundanceStart[[p, n]])[[1]], n]] + abundanceStart[[p, n]]
          abundanceStart[[p, n]] <- 0
        }
      }
    }
    
    diff <- goal - sum(storageStart[, n])
    if (diff != 0){
      #remove or add difference from throughout patches randomly
      distrDiff <- rmultinom(1, abs(diff), rep(abs(diff)/Npatches, Npatches) )
      storageStart[, n] <- storageStart[, n] + (sign(diff)*distrDiff)
    }
  }
  
  #insert speciesInfo 
  
  speciesInfo <- list(EMspecies = EMspecies, EMfitness = EMfitness, abundanceStart = abundanceStart, storageStart = storageStart, Nspecies = Nspecies, sigmaSpecies = sigmaSpecies)
  
  landscapeWithStorage$speciesInfo<- speciesInfo
  
  return(landscapeWithStorage)
}

#distribute sum so that it is equal among all N. Leftover is evenly split as much as possible.
distributeSum <- function(sum, N){
  leftover <- sum %% N
  result <- rep( ((sum-leftover) / N) , N)
  while (leftover > 0){
    result[[leftover]] <- result[[leftover]] + 1
    leftover <- leftover -1
  }
  return(result)
}
randomDistr<- function(sum, N){
  result <- rmultinom(1, sum, rep(1/N, N))
  
  return(result)
}


#CREATE LANDSCAPE VARIATIONS FUNCTION
#takes: landscape, percentVariations, frequencies, amplitudes
#returns: varyingLandscape (a list containing: landscapeVariation, fitnessVariation, variationListList)
createLandscapeVariations <- function(landscape,
                                      percentVariations, 
                                      frequencies, 
                                      amplitudes, 
                                      Tmax = Tmax){
  
  varyingLandscapeArray_perc_freq_amp <- array(rep(list(), length(percentVariations)*length(frequencies)*length(amplitudes)), c(length(percentVariations), length(frequencies), length(amplitudes)))
  
  for(perc in 1:length(percentVariations)){
    for(freq in 1:length(frequencies)){
      for(amp in 1:length(amplitudes)){
        
        percentVariation = percentVariations[perc]
        frequency = frequencies[freq]
        amplitude = amplitudes[amp]
        
        Npatches = landscape$Npatches
        EMpatches = landscape$EMpatches
        EMspecies = landscape$speciesInfo$EMspecies
        speciesInfo = landscape$speciesInfo
        sigmaSpecies = landscape$speciesInfo$sigmaSpecies
        XY = landscape$xy
        fitness = array(data = 0.0, dim = c(Npatches, Nspecies))	# initialize the array of habitat suitability per species per patch
        
        Epatches <- EMpatches
        
        #create landscape variation object, to be used in all models
        landscapeVariation <- list()
        fitnessVariation <- list()
        variationListList <- list()
        
        variationList <- c()
        
        for(t in 1:(Tmax)){
          #create changes every X frame based on frequency
          if(t%%frequency == 0 || t ==1){
            #create NON-FIX list of Patches that will vary
            
            NpatchVarying <- ceiling(percentVariation/100*Npatches)
            if(NpatchVarying < 2 && NpatchVarying > 0){
              NpatchVarying <- 2
            }
            if(NpatchVarying > 0){
              variationList <- sample(Npatches, NpatchVarying, replace = FALSE )
              #create list of all other patches (dependent)
              otherList <- c(1:Npatches)[-variationList]
            }else if (NpatchVarying == 0){
              variationList <- 1:Npatches
            }
            
            
            #if there is heterogeneity
            if(NpatchVarying > 0){ 
              # Vary each patch in variationList with a Normal distr around EM.
              for (p in variationList){
                Epatches[p] = rnorm(1, EMpatches[p], amplitude) # get a new E following normal distribution around EM. Amplitude drives sd of normal distribution.
                if(Epatches[p] > 100){Epatches[p] <- 100}
                else if(Epatches[p] < 0){Epatches[p] <- 0}
              }
              
              # Change other patches based on varied patches (inverse distance weighted)
              df.XYE <- data.frame(cbind(array(XY[variationList, ], c(length(variationList), 2)), array(Epatches[variationList], c(length(Epatches[variationList]),1)) ) )
              colnames(df.XYE) = c("x", "y", "E")
              
              #get interpolations (if needed)
              if(percentVariation < 100){
                idw.otherList <-idw(df.XYE[,3], df.XYE[,1:2], XY[otherList,], p = 5)
                #replace dependant Epatches values with interpolated values
                for (i in 1:nrow(idw.otherList)){
                  Epatches[[otherList[i]]] <- idw.otherList[i,]
                }
              }
            }else{
              #NO heterogeneity (whole space varies uniformly, like TNTB)
              Epatches[1:length(Epatches)] <- rnorm(1, EMpatches[[1]], amplitude)
            }
            
            
            variationListList[[length(variationListList)+1]] <- variationList
            landscapeVariation[[length(landscapeVariation)+1]] <- Epatches
          } else {
            #copy precedent frame
            variationListList[[length(variationListList)+1]] <- variationListList[[length(variationListList)]]
            landscapeVariation[[length(landscapeVariation)+1]] <- landscapeVariation[[length(landscapeVariation)]]
          }
        }
        
        varyingLandscape <- list(landscapeVariation = landscapeVariation, fitnessVariation = fitnessVariation, variationListList = variationListList)
        
        varyingLandscapeArray_perc_freq_amp[[perc, freq, amp]] <- varyingLandscape
        
      }
    }
  }
  
  
  return(varyingLandscapeArray_perc_freq_amp)
}

#ADD FITNESS FUNCTION
#Generates fitness lists based on the variations of the environment and species affinities and appends it to varyingLandscape
#takes: populatedLandscape, varyingLandscape
#returns: updated varyingLandscape
addVaryingFitness <- function(populatedLandscape, varyingLandscape){
  NSpecies <- populatedLandscape$NSpecies
  Npatches <- populatedLandscape$Npatches
  EMspecies <- populatedLandscape$speciesInfo$EMspecies
  sigmaSpecies <- populatedLandscape$speciesInfo$sigmaSpecies
  
  fitness <- matrix(NA, Npatches, Nspecies)
  for(perc in 1:length(varyingLandscape[, 1, 1])){
    for(freq in 1:length(varyingLandscape[1, , 1])){
      for(amp in 1:length(varyingLandscape[1, 1, ])){
        
        fitnessVariation <- list()
        for(t in 1:Tmax){
          #create fitnesses based on patches and species
          for (r in 1:Npatches){
            for (S in 1:Nspecies)
            {
              fitness[[r,S]] <- exp(-(EMspecies[S] - varyingLandscape[[perc, freq, amp]]$landscapeVariation[[t]][[r]])^2 / (2 * sigmaSpecies[S]^2))	# calculate habitat suitability per species per patch, based on manuscript equation [2]
            }
          }
          
          fitnessVariation[[length(fitnessVariation)+1]] <- fitness
        }
        varyingLandscape[[perc, freq, amp]]$fitnessVariation <- fitnessVariation
        
      }
    }
  }
  
  return(varyingLandscape)
}


#GET TRUE MEAN VARIATION FUNCTION
#Calculates the true mean of each patch's environmental condition through time, then appends it to varyingLandscape 
#takes: populatedLandscape, finalVaryingLandscape
#returns: updated varyingLandscape
getTrueMeanVariation <- function(populatedLandscape, finalVaryingLandscape){
  Npatches <- populatedLandscape$Npatches
  Nspecies <- populatedLandscape$speciesInfo$Nspecies
  EMspecies <- populatedLandscape$speciesInfo$EMspecies
  sigmaSpecies <- populatedLandscape$speciesInfo$sigmaSpecies
  
  for(perc in 1:length(finalVaryingLandscape[, 1, 1])){
    for(freq in 1:length(finalVaryingLandscape[1, , 1])){
      for(amp in 1:length(finalVaryingLandscape[1, 1, ])){
        variationList <- finalVaryingLandscape[[perc, freq, amp]]$landscapeVariation
        variationFlattened <- sapply(variationList, unlist)
        EMpatches <- apply(variationFlattened,1,  mean )
        SDpatches <- apply(variationFlattened, 1, sd)
        TrueEMfitness <- matrix(NA, Npatches, Nspecies)
        for (p in 1:Npatches){
          for (S in 1:Nspecies){
            TrueEMfitness[p,S] = exp(-(EMspecies[S] - EMpatches[p])^2 / (2 * sigmaSpecies[S]^2))	# calculate habitat suitability per species per patch, based on manuscript equation [2]
          }
        }
        finalVaryingLandscape[[perc, freq, amp]]$TrueEMfitness <- TrueEMfitness
        finalVaryingLandscape[[perc, freq, amp]]$TrueEMpatches <- EMpatches
        finalVaryingLandscape[[perc, freq, amp]]$SDpatches <- SDpatches
      }
    }
  }
  return(finalVaryingLandscape)
}


#NOT USED IN THIS PROJECT
# prepareAbundanceStart <- function(percentVariations, baseLandscape, landscapeVariations){
#     
#   #initialise array of result lists
#   abundanceStartModels <- array(rep(list(), length(percentVariations)), c(length(percentVariations)), dimnames = list(percNames))
# 
#   
#   abundanceStartModelsFlattened<- foreach (perc = 1:length(percentVariations), .combine = "c", .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape",   "Nspecies", "z", "J") )%dopar%{
#   #for(perc in 1:length(percentVariations)){
#      #execute model and store results in array
#   #abundanceStartModels[[perc]] <- 
#     model(z = z, 
#         J = J, 
#         Tmax = 100, 
#         d = d,
#         fitnessVariation = list(landscapeVariations[[perc, 1, 1]]$TrueEMfitness), 
#         filt = "niche", 
#         land = baseLandscape,
#         nicheDispersalLimitation = TRUE,
#         nicheAbuDependant = FALSE,
#         storageReduction = 0,
#         Niterations = 1,
#         abundanceStart = NULL)
#   }
#   
#   abundanceStartModels <- array(abundanceStartModelsFlattened, c(length(percentVariations)) )
#   
#   return(abundanceStartModels)
#   
# }


#LAUNCH MODELS FUNCTION
#This launches simulations for every combination of parameters chosen. Function is split into two (variation and static models). This was to avoid computational limits
#takes: percentVariations, frequencies, amplitues, storageReductions, modelTypes, baseLandscape, landscapeVariations
#returns: a list of results containing : varModels , staticModels

#Dynamic models
launchVarModels <- function(percentVariations,
                            frequencies,
                            amplitudes,
                            storageReductions,
                            modelTypes, 
                            populatedLandscape,
                            landscapeVariations,
                            nicheDispersalLimitation,
                            nicheAbuDependant,
                            abundanceStartModels = NULL,
                            temporalSampling = 0){
  
  #initialise array of result lists
  varModels <- array(rep(list(), length(percentVariations)*length(frequencies)*length(amplitudes)*length(storageReductions)*length(modelTypes)), c(length(percentVariations), length(frequencies), length(amplitudes), length(storageReductions), length(modelTypes)), dimnames = varNames)
  
  if(isParallelised == TRUE){
    #execute variation models
    varModelsFlattened <- foreach (mdlTyp = 1:length(modelTypes), .combine = "c", .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape", "nicheDispersalLimitation", "nicheAbuDependant", "storageReductions",   "Nspecies", "z", "J", "abundanceStartModels", "populatedLandscape", "graph_from_data_frame") )%:%
      foreach (st = 1:length(storageReductions), .combine = "c" , .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape", "nicheDispersalLimitation", "nicheAbuDependant", "storageReductions",   "Nspecies", "z", "J", "abundanceStartModels", "populatedLandscape", "graph_from_data_frame"))%:%
      foreach (amp = 1:length(amplitudes), .combine = "c", .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape", "nicheDispersalLimitation", "nicheAbuDependant", "storageReductions",   "Nspecies", "z", "J", "abundanceStartModels", "populatedLandscape", "graph_from_data_frame") )%:%
      #explore different storage decay timelengths
      foreach (freq = 1:length(frequencies), .combine = "c", .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape", "nicheDispersalLimitation", "nicheAbuDependant", "storageReductions",   "Nspecies", "z", "J", "abundanceStartModels", "populatedLandscape", "graph_from_data_frame") )%:%
      #do all 3 modelseLa
      foreach (perc = 1:length(percentVariations), .combine = "c", .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape", "nicheDispersalLimitation", "nicheAbuDependant", "storageReductions",   "Nspecies", "z", "J", "abundanceStartModels", "populatedLandscape", "graph_from_data_frame") )%dopar%{
        #execute model and store results in array
        model(z = z, 
              J = J, 
              Tmax = Tmax, 
              d = d,
              fitnessVariation = landscapeVariations[[perc, freq, amp]]$fitnessVariation, 
              filt = modelTypes[mdlTyp], 
              land = populatedLandscape,
              storageReduction = storageReductions[st],
              Niterations = Niterations,
              temporalSampling = temporalSampling,
              abundanceStart = abundanceStartModels[[perc, freq, amp, 1, mdlTyp]]$abu)
      }
    varModels <- array(varModelsFlattened, c(length(percentVariations), length(frequencies), length(amplitudes), length(storageReductions), length(modelTypes)) )
  }else{
    
    for (mdlTyp in 1:length(modelTypes)){
      for (st in 1:length(storageReductions)){
        for (amp in 1:length(amplitudes)){
          #explore different storage decay timelengths
          for (freq in 1:length(frequencies)){
            #do all 3 modelseLa
            for (perc in 1:length(percentVariations)){
              #execute model and store results in array
              varModels[[perc, freq, amp, st, mdlTyp]] <- model(z = z, 
                                                                J = J, 
                                                                Tmax = Tmax, 
                                                                d = d,
                                                                fitnessVariation = landscapeVariations[[perc, freq, amp]]$fitnessVariation, 
                                                                filt = modelTypes[mdlTyp], 
                                                                land = populatedLandscape,
                                                                storageReduction = storageReductions[st],
                                                                Niterations = Niterations,
                                                                temporalSampling = temporalSampling,
                                                                abundanceStart = abundanceStartModels[[perc, freq, amp, 1, mdlTyp]]$abu)
            }
          }
        }
      }
    }
  }
  print("TEST")
  
  
  return(list( varModels = varModels, info = list(J = J, StorageMultiple = storageMultiple, landscape = baseLandscape, frequencies = frequencies, percentVariations = percentVariations, amplitudes = amplitudes, storageReductions = storageReductions)))
  
  
}

#Static models
launchStaticModels <- function(percentVariations,
                               frequencies,
                               amplitudes,
                               storageReductions,
                               modelTypes, 
                               populatedLandscape,
                               landscapeVariations,
                               nicheDispersalLimitation,
                               nicheAbuDependant,
                               abundanceStartModels = NULL,
                               temporalSampling = 0){
  
  #initialise array of result lists
  staticModels <- array(rep(list(), length(percentVariations)*length(frequencies)*length(amplitudes)*length(storageReductions)*length(modelTypes)), c(length(percentVariations), length(frequencies), length(amplitudes), length(storageReductions), length(modelTypes)), dimnames = varNames)
  
  if(isParallelised == TRUE){
    #execute static models
    staticModelsFlattened <- foreach (mdlTyp = 1:length(modelTypes) , .combine = "c", .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape", "nicheDispersalLimitation", "nicheAbuDependant", "storageReductions",   "Nspecies", "z", "J", "abundanceStartModels", "populatedLandscape", "graph_from_data_frame"))%:%
      foreach (st = 1:length(storageReductions), .combine = "c", .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape", "nicheDispersalLimitation", "nicheAbuDependant", "storageReductions",   "Nspecies", "z", "J", "abundanceStartModels", "populatedLandscape", "graph_from_data_frame"))%:%
      foreach (amp = 1:length(amplitudes) , .combine = "c", .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape", "nicheDispersalLimitation", "nicheAbuDependant", "storageReductions",   "Nspecies", "z", "J", "abundanceStartModels", "populatedLandscape", "graph_from_data_frame"))%:%
      #explore different storage decay timelengths
      foreach (freq = 1:length(frequencies), .combine = "c", .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape", "nicheDispersalLimitation", "nicheAbuDependant", "storageReductions",   "Nspecies", "z", "J", "abundanceStartModels", "populatedLandscape", "graph_from_data_frame") )%:%
      #do all 3 models
      foreach (perc = 1:length(percentVariations), .combine = "c", .export = c("model", "Niterations", "Tmax", "d", "landscapeVariations", "modelTypes", "baseLandscape", "nicheDispersalLimitation", "nicheAbuDependant", "storageReductions",   "Nspecies", "z", "J", "abundanceStartModels", "populatedLandscape", "graph_from_data_frame")) %dopar%{
        model(z = z, 
              J = J, 
              Tmax = Tmax, 
              d = d,
              fitnessVariation = list(landscapeVariations[[perc, freq, amp]]$TrueEMfitness), 
              filt = modelTypes[mdlTyp], 
              land = populatedLandscape, 
              storageReduction = storageReductions[st],
              Niterations = Niterations,
              temporalSampling = temporalSampling,
              abundanceStart = abundanceStartModels[[perc, freq, amp, 1, mdlTyp]]$abu)
      }
    staticModels <- array(staticModelsFlattened, c(length(percentVariations), length(frequencies), length(amplitudes), length(storageReductions), length(modelTypes)) )
    
  }else{
    for (mdlTyp in 1:length(modelTypes)){
      for (st in 1:length(storageReductions)){
        for (amp in 1:length(amplitudes)){
          #explore different storage decay timelengths
          for (freq in 1:length(frequencies)){
            #do all 3 modelseLa
            for (perc in 1:length(percentVariations)){
              #execute model and store results in array
              staticModels[[perc, freq, amp, st, mdlTyp]] <- model(z = z, 
                                                                   J = J, 
                                                                   Tmax = Tmax, 
                                                                   d = d,
                                                                   fitnessVariation = list(landscapeVariations[[perc, freq, amp]]$TrueEMfitness), 
                                                                   filt = modelTypes[mdlTyp], 
                                                                   land = populatedLandscape, 
                                                                   storageReduction = storageReductions[st],
                                                                   Niterations = Niterations,
                                                                   temporalSampling = temporalSampling,
                                                                   abundanceStart = abundanceStartModels[[perc, freq, amp, 1, mdlTyp]]$abu)
            }
          }
        }
      }
    }
  }
  
  return(list( staticModels = staticModels, info = list(J = J, StorageMultiple = storageMultiple, landscape = baseLandscape, frequencies = frequencies, percentVariations = percentVariations, amplitudes = amplitudes, storageReductions = storageReductions)))
  
}


#CREATE RESULTS ARRAY
#create an empty array to store information order by static, var, amp, freq, perc, stor, mdlTyp
createResultsArray <- function(percentVariations, frequencies, amplitudes, storageReductions = NULL, modelTypes = NULL, varNm = varNames){
  
  if(!is.null(storageReductions) && !is.null(modelTypes)){
    varModels <- array(rep(list(), length(percentVariations)*length(frequencies)*length(amplitudes)*length(storageReductions)*length(modelTypes)), c(length(percentVariations), length(frequencies), length(amplitudes), length(storageReductions), length(modelTypes)), dimnames = varNames)
    staticModels <- array(rep(list(), length(percentVariations)*length(frequencies)*length(amplitudes)*length(storageReductions)*length(modelTypes)), c(length(percentVariations), length(frequencies), length(amplitudes), length(storageReductions), length(modelTypes)), dimnames = varNames)
    
  }else{
    varModels <- array(rep(list(), length(percentVariations)*length(frequencies)*length(amplitudes)), c(length(percentVariations), length(frequencies), length(amplitudes)), dimnames = varNm[1:3])
    staticModels <- array(rep(list(), length(percentVariations)*length(frequencies)*length(amplitudes)), c(length(percentVariations), length(frequencies), length(amplitudes)), dimnames = varNm[1:3])
    
  }
  resultsArrayList <- list(staticModels = staticModels, varModels = varModels)
  
  return(resultsArrayList)
}

##### NICHE-NEUTRAL METRIC CALCULATION FUNCTIONS #####
####____________________________________________####
####ABUNDANCE BETA-NULL DEVIATION (Tucker et al. 2016)####
abuBetaDev <- function(modelResults){
  # Calculate the number of cores
  
  #no_cores set to 1 to avoid automatically saturating computer resources
  #set to detectCores()-1 to increase speed
  
  #no_cores <- detectCores()-1
  no_cores <- 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  mdlTypeList <- c("neutral", "continuum")
  
  #resultArrays <- createResultsArray
  
  staticAbuBetaDevFlattened <- foreach (mdlTyp = 1:length(mdlTypeList) , .combine = "c", .export = c( "calculateAbuBetaDev", "modelTypes", "storageReductions"), .packages = c("vegan"))%:%
    foreach (st = 1:length(modelResults$info$storageReductions), .combine = "c", .export = c("calculateAbuBetaDev", "mdlTypeList", "storageReductions"), .packages = c("vegan"))%:%
    foreach (amp = 1:length(modelResults$info$amplitudes) , .combine = "c", .export = c("calculateAbuBetaDev", "mdlTypeList", "storageReductions"), .packages = c("vegan"))%:%
    #explore different storage decay timelengths
    foreach (freq = 1:length(modelResults$info$frequencies), .combine = "c", .export = c("calculateAbuBetaDev", "mdlTypeList", "storageReductions") , .packages = c("vegan"))%:%
    #do all 3 models
    foreach (perc = 1:length(modelResults$info$percentVariations), .combine = "c", .export = c("calculateAbuBetaDev", "mdlTypeList", "storageReductions"), .packages = c("vegan")) %dopar%{
      
      calculateAbuBetaDev(modelResults$staticModels[[perc, freq, amp, st, mdlTyp]]$abu)[[1]]$mean
      
    }
  varAbuBetaDevFlattened <- foreach (mdlTyp = 1:length(mdlTypeList) , .combine = "c", .export = c("calculateAbuBetaDev", "modelTypes", "storageReductions"), .packages = c("vegan"))%:%
    foreach (st = 1:length(modelResults$info$storageReductions), .combine = "c", .export = c("calculateAbuBetaDev", "mdlTypeList", "storageReductions"), .packages = c("vegan"))%:%
    foreach (amp = 1:length(modelResults$info$amplitudes) , .combine = "c", .export = c("calculateAbuBetaDev", "mdlTypeList", "storageReductions"), .packages = c("vegan"))%:%
    #explore different storage decay timelengths
    foreach (freq = 1:length(modelResults$info$frequencies), .combine = "c", .export = c("calculateAbuBetaDev", "mdlTypeList", "storageReductions") , .packages = c("vegan"))%:%
    #do all 3 models
    foreach (perc = 1:length(modelResults$info$percentVariations), .combine = "c", .export = c("calculateAbuBetaDev", "mdlTypeList", "storageReductions"), .packages = c("vegan")) %dopar%{
      calculateAbuBetaDev(modelResults$varModels[[perc, freq, amp, st, mdlTyp]]$abu)[[1]]$mean
    }
  
  stopImplicitCluster()
  
  # staticArray <- array(staticAbuBetaDevFlattened, c(length(percentVariations), length(frequencies), length(amplitudes), length(storageReductions), length(modelTypes)) )
  # varArray <- array(varAbuBetaDevFlattened, c(length(percentVariations), length(frequencies), length(amplitudes), length(storageReductions), length(modelTypes)) )
  # 
  #results <- list(static = staticAbuBetaDevFlattened, var = varAbuBetaDevFlattened)
  
  staticArray <- array(staticAbuBetaDevFlattened, c( length(modelResults$info$percentVariations), length(modelResults$info$frequencies), length(modelResults$info$amplitudes), length(modelResults$info$storageReductions), length(mdlTypeList)) )
  varArray <- array(varAbuBetaDevFlattened, c( length(modelResults$info$percentVariations), length(modelResults$info$frequencies), length(modelResults$info$amplitudes), length(modelResults$info$storageReductions), length(mdlTypeList)) )
  
  
  return(results <- list(staticModels = staticArray, varModels = varArray))
}
calculateAbuBetaDev <- function(abundances){
  #vector of all mean null deviations for replicates
  abund_null_vector <- c()
  
  #loop through all replicates
  for (iter in 1:length(abundances) ){ 
    
    abundance <- abundances[[iter]]
    #null model parameterization
    dune <- abundance
    #rand determines null model replicates (can take a very long time if at 999)
    rand <- 1
    patches <- 100
    
    ### Prepare and calculate abundance beta-null deviation metric
    ## Adjusted from Stegen et al 2012 GEB
    bbs.sp.site <- dune
    null.alphas <- matrix(NA, ncol(dune), rand)
    expected_beta <- matrix(NA, 1, rand)
    null.gamma <- matrix(NA, 1, rand)
    null.alpha.comp <- numeric()
    bucket_bray_res <- matrix(NA, patches, rand)
    
    bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)) #tranform into presence absence?
    mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site) #mean.alpha
    gamma <- ncol(bbs.sp.site) #gamma
    obs_beta <- 1-mean.alpha/gamma
    obs_beta_all <- 1-rowSums(bbs.sp.site)/gamma
    
    ##Generate null patches
    for (randomize in 1:rand) {  
      null.dist = dune
      #randomly reorganise every species column (keeping total species abundances)
      for (species in 1:ncol(null.dist)) {
        tot.abund = sum(null.dist[,species])
        null.dist[,species] = 0
        for (individual in 1:tot.abund) {
          sampled.site = sample(c(1:nrow(bbs.sp.site)), 1)
          null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
        }
      }
      
      ##Calculate null deviation for null patches and store
      null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
      null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
      expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
      null.alpha <- mean(null.alphas[,randomize])
      null.alpha.comp <- c(null.alpha.comp, null.alpha)
      
      bucket_bray <- as.matrix(vegdist(null.dist, "bray"))
      diag(bucket_bray) <- NA
      bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE)
    } ## end randomize loop
    
    ## Calculate beta-diversity for obs metacommunity
    beta_comm_abund <- vegdist(dune, "bray")
    res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
    diag(res_beta_comm_abund) <- NA
    # output beta diversity (Bray)
    beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)
    
    # output abundance beta-null deviation
    abund_null_dev <- mean(beta_div_abund_stoch - mean(bucket_bray_res) )
    abund_null_vector <- c(abund_null_vector, abund_null_dev)
  }
  abund_null_dev_mean <- mean(abund_null_vector)
  abund_null_dev_mean10 <- mean(abund_null_vector[1:10])
  abund_null_dev_mean5 <- mean(abund_null_vector[1:5])
  
  abund_null_dev_SD <- sd(abund_null_vector)
  
  return( list(list(mean = abund_null_dev_mean, mean5 = abund_null_dev_mean5, mean10 = abund_null_dev_mean10, SD = abund_null_dev_SD)) )
}
#CONVERT ABUBETADEV
convertAbuBetaDev <- function(abubetadev, speciesN, z, model = 3){
  convertedResults <- createResultsArray(percentVariations, 
                                         frequencies,
                                         amplitudes 
  )
  
  finalTable <- c()
  #convert data
  for(freq in 1:length(frequencies)){
    for(perc in 1:length(percentVariations)){
      for(amp in 1:length(amplitudes)){
        
        #get mean of all iterations
        convertedResults$staticModels[[perc, freq, amp]] <- abubetadev$staticModels[[perc, freq, amp, 1, model]]$mean
        convertedResults$varModels[[perc, freq, amp]] <- abubetadev$varModels[[perc, freq, amp, 1, model]]$mean
        
      }
    }
    #melt into tables
    staticResults <- melt(matrix(unlist(convertedResults$staticModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedResults$static[ , freq, ])))
    varResults <- melt(matrix(unlist(convertedResults$varModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedResults$static[ , freq, ])))
    #add Species column
    staticResults$Sp <- rep(speciesN, nrow(staticResults))
    varResults$Sp <- rep(speciesN, nrow(staticResults))
    #add z column
    staticResults$z <- rep(z, nrow(staticResults))
    varResults$z <- rep(z, nrow(staticResults))
    #add period column
    staticResults$Freq <- rep(1/frequencies[[freq]], nrow(staticResults))
    varResults$Freq <- rep(1/frequencies[[freq]], nrow(staticResults))
    #add model column
    staticResults$Model <- rep("static", nrow(staticResults))
    varResults$Model <- rep("variable", nrow(staticResults))
    finalTable <- rbind(finalTable, 
                        staticResults,
                        varResults)
  }
  return(finalTable)
}

#### VARIATION PARTITIONING ####
convertVariationPartitioning <- function(modelresults, speciesN, z,baseLandscape, finalVaryingLandscape, MEM.select = NULL, model =3){
  
  convertedResults <- createResultsArray(modelresults$info$percentVariations, 
                                         modelresults$info$frequencies,
                                         modelresults$info$amplitudes 
  )
  convertedEandS <- createResultsArray(modelresults$info$percentVariations, 
                                       modelresults$info$frequencies,
                                       modelresults$info$amplitudes 
  )
  convertedResiduals <- createResultsArray(modelresults$info$percentVariations, 
                                           modelresults$info$frequencies,
                                           modelresults$info$amplitudes 
  )
  
  #z <- modelresults$info$landscape$z
  varTypes <-  c("staticModels", "varModels")
  
  finalTable <- c()
  finalTable_EandS <- c()
  finalTable_Residuals <- c()
  
  
  modelForVar <- model-1
  modelForStatic <- model-1
  for(freq in 1:length(modelresults$info$frequencies)){
    for(perc in 1:length(modelresults$info$percentVariations)){
      #TO CHECK potentially re-forward select for every perc
      
      for(amp in 1:length(modelresults$info$amplitudes)){
        #generate container for MEM (for both static and variable models)
        MEM.select <- list()
        for(Niter in 1:length(modelresults$varModels[[1,1,1,1,1]]$abu)){
          
          
          partitionResults <-rep(list( rep(NA, 4)) , 2 )
          
          for (var in 1:2){
            #correct model for missing niche model in static models: neutral,usually 2, becomes 1
            if (var == 1){model <- modelForStatic}else{model <- modelForVar}
            #get data
            Y <- modelresults[[varTypes[[var]]]][[ perc, freq, amp, 1, model]]$abu[[Niter]]
            if (varTypes[[var]] == "staticModels"){
              #static
              X <- data.frame(finalVaryingLandscape[[perc, freq, amp]]$TrueEMpatches )
            }else{
              #variable
              #last time frame ( >> very little explained by environment or both)
              #with Tmax = 2500, last time frame is at 100 of landscape variations. As there are only 200 timesteps of variation that are repeated.
              
              #get last timeframe based on repetitions
              # lastTimeFrame <- floor(Tmax%%(200*frequencies[freq])/frequencies[freq] )
              # if(lastTimeFrame == 0){ lastTimeFrame <- 200}
              X <- data.frame(finalVaryingLandscape[[perc, freq, amp]]$landscapeVariation[[ 200 ]])
              #average
              #X <- data.frame(finalVaryingLandscape[[perc, freq, amp]]$TrueEMpatches )
            }
            
            names(X) <- "X"
            #helinger transform data
            #Ycorr <- sqrt(Y/outer(apply(Y, 1,sum), rep(1, ncol(Y)), "*"))
            #Yfinal <- resid(lm(as.matrix(Ycorr) ~ as.matrix(modelresults$info$landscape$xy)))
            Y<-decostand(Y, "hellinger")
            
            #define MEM only on first 
            #OLD - is.null(MEM.select) && var == 1 && model == 1
            #calculate MEM for first replicate only
            if(Niter == 1){
              #define connectivity
              neighbours <- dnearneigh(modelresults$info$landscape$xy, 0, 100)
              distances <- nbdists(neighbours, baseLandscape$xy)
              
              #define distance decay
              fdist <- lapply(distances, function(x) 1/exp(z * x))
              #fdist <- lapply(distances, function(x) 1 - x/max(dist(modelresults$info$landscape$xy)))
              #fdist <- modelresults$info$landscape$D
              
              #get spatial information
              listw.final <- nb2listw(neighbours, glist = fdist, style = "B")
              #pcnmResults <- pcnm(dist(modelresults$info$landscape$xy))
              MEMresults <- scores.listw(listw.final, MEM.autocor = "positive")
              #MEMresults <- MEMresults$vectors[ ,MEMresults$values > 0]
              
              # #assume true after initial tests
              # #if (anova.cca(rda(Y, MEMresults), permutations = 9999)$Pr[1] <= 0.05) {
              # Global adjusted R-squared of the model
              R2adj <- RsquareAdj(rda(Y, MEMresults))$adj.r.squared
              # FWD with two stopping criteria
              fsel <- forward.sel(Y, MEMresults, adjR2thresh = R2adj, nperm = 999)
              # We order the selected MEM by decreasing eigenvalue
              sorted_sel <- fsel$order
              # Object containing the selected MEM
              
              MEM.select[[var]] <- as.data.frame(MEMresults)[, c(sorted_sel)]
              print("DONE")
              # #} else print("No significant spatial autocorrelation was detected in the response")
              # 
            }
            
            #calculate and plot VariationPartitioning
            mod <- varpart(Y, ~X, MEM.select[[var]], data = X )
            #showvarparts(2, Xnames = c("Env", "Space") , bg = 2:5 ) 
            #plot(mod,bg=2:5)
            
            #store results [environment, both, space, residuals]
            partitionResults[[var]] <- mod$part$indfract$Adj.R.squared
            
          }
          #convert preliminary results to index
          varPartStatic <- calculateVarPartIndex(partitionResults[[2]])
          varPartVar <- calculateVarPartIndex(partitionResults[[1]])                 
          
          convertedResults$staticModels[[perc, freq, amp]][[Niter]] <- varPartStatic[[1]]
          convertedResults$varModels[[perc, freq, amp]][[Niter]] <- varPartVar[[1]]
          
          convertedEandS$staticModels[[perc, freq, amp]][[Niter]] <- varPartStatic[[2]]
          convertedEandS$varModels[[perc, freq, amp]][[Niter]] <- varPartVar[[2]]
          
          convertedResiduals$staticModels[[perc, freq, amp]][[Niter]] <- varPartStatic[[3]]
          convertedResiduals$varModels[[perc, freq, amp]][[Niter]] <- varPartVar[[3]]
        }
        convertedResults$staticModels[[perc, freq, amp]] <- mean(convertedResults$staticModels[[perc, freq, amp]])
        convertedResults$varModels[[perc, freq, amp]] <- mean(convertedResults$varModels[[perc, freq, amp]])
        
        convertedEandS$staticModels[[perc, freq, amp]] <- mean(convertedEandS$staticModels[[perc, freq, amp]])
        convertedEandS$varModels[[perc, freq, amp]] <- mean(convertedEandS$varModels[[perc, freq, amp]])
        
        convertedResiduals$staticModels[[perc, freq, amp]] <- mean(convertedResiduals$staticModels[[perc, freq, amp]])
        convertedResiduals$varModels[[perc, freq, amp]] <- mean(convertedResiduals$varModels[[perc, freq, amp]])
      }
    }
    #results
    #melt into tables
    staticResults <- melt(matrix(unlist(convertedResults$staticModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedResults$static[ , freq, ])))
    varResults <- melt(matrix(unlist(convertedResults$varModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedResults$varModels[ , freq, ])))
    #add Species column
    staticResults$Sp <- rep(speciesN, nrow(staticResults))
    varResults$Sp <- rep(speciesN, nrow(varResults))
    #add z column
    staticResults$z <- rep(z, nrow(staticResults))
    varResults$z <- rep(z, nrow(varResults))
    
    #add period column
    staticResults$Freq <- rep(1/modelresults$info$frequencies[[freq]], nrow(staticResults))
    varResults$Freq <- rep(1/modelresults$info$frequencies[[freq]], nrow(varResults))
    #add model column
    staticResults$Model <- rep("static", nrow(staticResults))
    varResults$Model <- rep("variable", nrow(varResults))
    
    finalTable <- rbind(finalTable,
                        staticResults,
                        varResults)
    
    #EandS
    #melt into tables
    staticResults <- melt(matrix(unlist(convertedEandS$staticModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedEandS$static[ , freq, ])))
    varResults <- melt(matrix(unlist(convertedEandS$varModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedEandS$varModels[ , freq, ])))
    #add Species column
    staticResults$Sp <- rep(speciesN, nrow(staticResults))
    varResults$Sp <- rep(speciesN, nrow(varResults))
    #add z column
    staticResults$z <- rep(z, nrow(staticResults))
    varResults$z <- rep(z, nrow(varResults))
    
    #add period column
    staticResults$Freq <- rep(1/modelresults$info$frequencies[[freq]], nrow(staticResults))
    varResults$Freq <- rep(1/modelresults$info$frequencies[[freq]], nrow(varResults))
    #add model column
    staticResults$Model <- rep("static", nrow(staticResults))
    varResults$Model <- rep("variable", nrow(varResults))
    
    finalTable_EandS <- rbind(finalTable_EandS,
                              staticResults,
                              varResults)
    #residuals
    #melt into tables
    staticResults <- melt(matrix(unlist(convertedResiduals$staticModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedResiduals$static[ , freq, ])))
    varResults <- melt(matrix(unlist(convertedResiduals$varModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedResiduals$varModels[ , freq, ])))
    #add Species column
    staticResults$Sp <- rep(speciesN, nrow(staticResults))
    varResults$Sp <- rep(speciesN, nrow(varResults))
    #add z column
    staticResults$z <- rep(z, nrow(staticResults))
    varResults$z <- rep(z, nrow(varResults))
    
    #add period column
    staticResults$Freq <- rep(1/modelresults$info$frequencies[[freq]], nrow(staticResults))
    varResults$Freq <- rep(1/modelresults$info$frequencies[[freq]], nrow(varResults))
    #add model column
    staticResults$Model <- rep("static", nrow(staticResults))
    varResults$Model <- rep("variable", nrow(varResults))
    
    finalTable_Residuals <- rbind(finalTable_Residuals,
                                  staticResults,
                                  varResults)
  }
  
  return(list(list(finalTable, finalTable_EandS, finalTable_Residuals)))
}
calculateVarPartIndex <- function(results){
  #results in form : [environment, both, space, residuals]
  env <- results[[1]]
  spc <- results[[3]]
  EandS <- results[[2]]
  residuals <- results[[4]]
  
  #-1 = env, 0 = both, +1 = spc
  calculatedIndex <- (spc-env)
  # return(list(calculatedIndex, residuals)
  return(list(calculatedIndex, EandS, residuals))
}



#### SPECIES VARIATION ACROSS REPLICATES (SVAR) ####
convertSimpson <- function(modelresults, speciesN, z, model = 3){
  convertedResults <- createResultsArray(percentVariations, 
                                         frequencies,
                                         amplitudes 
  )
  
  finalTable <- c()
  #convert data
  for(freq in 1:length(frequencies)){
    for(perc in 1:length(percentVariations)){
      for(amp in 1:length(amplitudes)){
        
        #get mean of all iterations
        convertedResults$staticModels[[perc, freq, amp]] <- calculateSimpson( modelresults$staticModels[[perc, freq, amp, 1, model-1]]$abu)[[1]]
        convertedResults$varModels[[perc, freq, amp]] <- calculateSimpson(modelresults$varModels[[perc, freq, amp, 1, model-1]]$abu)[[1]]
        
      }
    }
    #melt into tables
    staticResults <- melt(matrix(unlist(convertedResults$staticModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedResults$static[ , freq, ])))
    varResults <- melt(matrix(unlist(convertedResults$varModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedResults$static[ , freq, ])))
    #add Species column
    staticResults$Sp <- rep(speciesN, nrow(staticResults))
    varResults$Sp <- rep(speciesN, nrow(staticResults))
    #add z column
    staticResults$z <- rep(z, nrow(staticResults))
    varResults$z <- rep(z, nrow(staticResults))
    
    #add period column
    staticResults$Freq <- rep(1/frequencies[[freq]], nrow(staticResults))
    varResults$Freq <- rep(1/frequencies[[freq]], nrow(staticResults))
    #add model column
    staticResults$Model <- rep("static", nrow(staticResults))
    varResults$Model <- rep("variable", nrow(staticResults))
    
    finalTable <- rbind(finalTable, 
                        staticResults,
                        varResults)
  }
  return(finalTable)
}

#Calculate measure adapted from Avi Bar-Massada et al. (2014) to work with patch network rather than grid (see equation X in manuscript)
calculateSimpson <- function(modelresultsWithReplicates){
  
  
  indexList <- c()
  for (patch in 1:100){
    patchAbundances <- c()
    #go through all replicates
    repNum <- length(modelresultsWithReplicates)
    for (rep in 1:repNum){
      #make a series of columns ,each column is a replicate of abundances
      patchAbundances <- cbind(patchAbundances, modelresultsWithReplicates[[rep]][patch,])
    }
    #get highest abundance in patch for all replicates
    mostCommon<- max.col(t(patchAbundances))
    mostCommonSpecies <- unique(mostCommon)                       
    mostCommonSpeciesCount <- length(mostCommonSpecies )
    
    p <- c()
    for (sp in mostCommonSpecies){
      p <- c(p, ( length(mostCommon[mostCommon == sp])) )
    }
    #Simpson index
    index <- 1- ( sum(p*(p-1)) / (repNum * repNum-1) )
    
    indexList <- c(indexList, index)
  }
  return(list(mean = mean(indexList), sd = sd(indexList), map = indexList))
}

#### EVENNESS ####
convertSADs <- function(modelresults, speciesN, z, model = 3, perPatch = FALSE, includeExtinctions = TRUE){
  convertedResults <- createResultsArray(modelresults$info$percentVariations, 
                                         modelresults$info$frequencies,
                                         modelresults$info$amplitudes 
  )
  finalTable <- c()
  #convert data
  for(freq in 1:length(modelresults$info$frequencies)){
    for(perc in 1:length(modelresults$info$percentVariations)){
      for(amp in 1:length(modelresults$info$amplitudes)){
        for(Niter in 1:length(modelresults$varModels[[1,1,1,1,1]]$abu)){
          
          staticContinuum <- calculate_E9_RADindex(modelresults$staticModels[[perc, freq, amp, 1, model-1]]$abu[[Niter]])
          
          varContinuum <- calculate_E9_RADindex(modelresults$varModels[[perc, freq, amp, 1, model-1]]$abu[[Niter]])
          
          convertedResults$staticModels[[perc, freq, amp]][[Niter]] <- staticContinuum
          convertedResults$varModels[[perc, freq, amp]][[Niter]] <- varContinuum
        }
        #get mean of all iterations
        convertedResults$staticModels[[perc, freq, amp]] <- mean(convertedResults$staticModels[[perc, freq, amp]])
        convertedResults$varModels[[perc, freq, amp]] <- mean(convertedResults$varModels[[perc, freq, amp]])
      }
    }
    #melt into tables
    staticResults <- melt(matrix(unlist(convertedResults$staticModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedResults$static[ , freq, ])))
    varResults <- melt(matrix(unlist(convertedResults$varModels[ , freq, ]), nrow= 4, dimnames = dimnames(convertedResults$static[ , freq, ])))
    #add Species column
    staticResults$Sp <- rep(speciesN, nrow(staticResults))
    varResults$Sp <- rep(speciesN, nrow(staticResults))
    #add z column
    staticResults$z <- rep(z, nrow(staticResults))
    varResults$z <- rep(z, nrow(staticResults))
    
    #add period column
    staticResults$Freq <- rep(1/frequencies[[freq]], nrow(staticResults))
    varResults$Freq <- rep(1/frequencies[[freq]], nrow(staticResults))
    #add model column
    staticResults$Model <- rep("static", nrow(staticResults))
    varResults$Model <- rep("variable", nrow(staticResults))
    
    finalTable <- rbind(finalTable, 
                        staticResults,
                        varResults)
  }
  
  return(finalTable)
}
calculate_E9_RADindex <- function(abundances, perPatch = FALSE, includeExtinctions = FALSE){
  
  
  
  if(perPatch == FALSE){
    #sort abundances into SADs
    sortedAbundance <- sort(colSums(abundances), decreasing = TRUE)
    
    if(includeExtinctions == FALSE){
      sortedAbundance <- sortedAbundance[sortedAbundance>0]
    }
    #use mathematical formula
    #S = number of species
    S <- length(sortedAbundance)
    #p = proportions of each species
    p <- sortedAbundance/sum(sortedAbundance)
    index <- 1 - ( ( (S*sum(p^2) - 1) / (S-1) ) )^(1/2)
  }
  return(index)
}



#### FUNCTION FOR COMBINING, ORGANIZING AND LABELING NICHE-NEUTRAL METRIC CALCULATIONS####
#Function to organize and correctly label niche, neutral and continuum data into a single table, 
#as well as standardize or inverse if necessary.
#This is mostly to prepare the data for plotting.
correctTableTerms <- function( nicheTable = NULL, neutralTable, continuumTable, standardise = FALSE,  inverse = FALSE){
  #function to prepare data for analysis and plotting
  if(standardise == TRUE){
    neutralTableValue <- (neutralTable$value - neutralTable$value)
    # /( nicheTable$value -  neutralTable$value)
    nicheTableValue <- (nicheTable$value -  neutralTable$value)
    # /(nicheTable$value -  neutralTable$value )
    
    continuumTableValue <- (continuumTable$value - neutralTable$value )
    # /(nicheTable$value -  neutralTable$value )
    
    neutralTable$value <- neutralTableValue
    nicheTable$value <- nicheTableValue
    continuumTable$value <- continuumTableValue
    
    if(inverse == TRUE){
      continuumTable$value <- (continuumTable$value-1)*-1
      nicheTable$value <- (nicheTable$value-1)*-1
      neutralTable$value <- (neutralTable$value-1)*-1
      
    }
  }
  
  
  if(is.null(nicheTable)){
    tableList <- list( neutralTable, continuumTable)
    modelType <- list("neutral", "continuum")
  }else{
    tableList <- list( nicheTable, neutralTable, continuumTable)
    modelType <- list("niche","neutral", "continuum")
  }
  
  finalTable <- list()
  
  for(model in 1:length(modelType)){
    
    tableList[[model]]$Var1 <- as.character(tableList[[model]]$Var1)
    tableList[[model]]$Var2 <- as.character(tableList[[model]]$Var2)
    
    tableList[[model]]$Var1[ tableList[[model]]$Var1 == "percVar_5"] <- 5
    tableList[[model]]$Var1[ tableList[[model]]$Var1 == "percVar_10"] <- 10
    tableList[[model]]$Var1[ tableList[[model]]$Var1 == "percVar_25"] <- 25
    tableList[[model]]$Var1[ tableList[[model]]$Var1 == "percVar_50"] <- 50
    
    tableList[[model]]$Var1 <- factor(tableList[[model]]$Var1,  levels = c(5, 10, 25, 50))
    
    tableList[[model]]$Var2[ tableList[[model]]$Var2 == "amp_1"] <- 1
    tableList[[model]]$Var2[ tableList[[model]]$Var2 == "amp_5"] <- 5
    tableList[[model]]$Var2[ tableList[[model]]$Var2 == "amp_10"] <- 10
    tableList[[model]]$Var2[ tableList[[model]]$Var2 == "amp_15"] <- 15
    tableList[[model]]$Var2 <- factor(tableList[[model]]$Var2,  levels = c(1, 5, 10, 15))
    
    
    names(tableList[[model]])[[6]] <- "Freq"
    
    #add a modelType category
    tableList[[model]]$modelType <- rep(modelType[[model]], nrow(tableList[[model]]))
    
    #fuse
    finalTable <- rbind(finalTable, tableList[[model]])
  }
  
  colnames(finalTable)[1] <- "Pattern"
  colnames(finalTable)[2] <- "Amplitude"
  
  #round and factorise frequency
  finalTable$Freq <- as.factor(round(finalTable$Freq, 1))
  
  
  print("DONE")
  
  return(finalTable)
}

#interupts Sourcing 
browser()

#### MODEL RUNS ####
####___________####
#Model parameters
#number of times each model is repeated
#set at 2 for minimal rapid results (SVAR requires a minimum of 2 replicates). 
#The study was done with 50 replicates but this is time consuming.
Niterations <- 2

#death rate occuring at every timestep
d <- 0.25 

#types of models
modelTypes = c("neutral", "continuum")

#environmental stochasticity parameters to combine 
#(amplitudes = magnitudes in manuscript)
#(percentVariations = SATEF in manuscript)
percentVariations = c( 5, 10 , 25, 50)
frequencies = c(1, 3,  5)
amplitudes = c( 1, 5, 10, 15)

#storage not used for this project (ignored in simulations)
storageMultiple <- 1 #multiple of Jpatches that serves as storage (stored individual for every living individual)
storageReductions = c(0)

# give names to parameters
percNames <- paste("percVar_", percentVariations, sep = "")
freqNames <- paste("freq_", frequencies, sep = "")
ampNames <- paste("amp_", amplitudes, sep = "")
storRedNames <- paste("storRed_", storageReductions, sep = "" )

# var model names
varNames <- list(percNames, freqNames, ampNames, storRedNames, modelTypes)
# static model names
staticNames <- list(storRedNames, modelTypes)

Nspecies <- 25 # number of species
# number of patches
N <- 100 
#total sites (to be divided among all patches)
J <- 20000 
#amount of total J patches is slightly altered to be divisible by amount of species
J <- Nspecies * round(J/Nspecies) 
#size of landscape (number of horizontal and vertical units)
size <- 100

#Niche widths - specialists and generalist are identical for this project
specialist <- 4 #sigma value for specialists
generalist <- 4 #sigma value for generalists

z = 0.3 #dispersal decay

#RUN MODELS

#generate base landscape
baseLandscape <<- createLandscape(N, J, z)

#storage has no effect in this project (storage effect is ignored)
#however it is still generated to avoid missing objects in the simulation
landscapeWithStorage <<- addStorageToLandscape(baseLandscape, storageMultiple)

Tmax <- 200 # timesteps for environmental stochasticity
#create environmental stochasticity table (limited to 200 timesteps which repeat indefinately)
varyingLandscapesArray <<- createLandscapeVariations(baseLandscape, percentVariations, frequencies, amplitudes, Tmax = Tmax)

#### Model with 25 SPECIES and z = 0.3 ####
Nspecies = 25 #species richness

#add species information
populatedLandscape_Nsp25 <<- populateLandscapeWithSpecies(landscapeWithStorage, Nspecies)
#calculate and append fitness tables for every species
finalVaryingLandscape_Nsp25 <<- addVaryingFitness(populatedLandscape_Nsp25, varyingLandscapesArray)
#calculate and append true mean E value of each patch
finalVaryingLandscape_Nsp25 <<- getTrueMeanVariation(populatedLandscape_Nsp25, finalVaryingLandscape_Nsp25)


Tmax <- 2000 #total timesteps

if(isParallelised == TRUE){
  # Calculate the number of cores
  no_cores <- detectCores()-1
  # Initiate cluster
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
}

#run static and dynamic models
VarModelresults <- launchVarModels( percentVariations, frequencies, amplitudes, storageReductions, modelTypes, populatedLandscape = populatedLandscape_Nsp25, landscapeVariations = finalVaryingLandscape_Nsp25, nicheDispersalLimitation = TRUE, nicheAbuDependant = FALSE)
StaticModelresults <- launchStaticModels(percentVariations, frequencies, amplitudes, storageReductions, modelTypes, populatedLandscape = populatedLandscape_Nsp25, landscapeVariations = finalVaryingLandscape_Nsp25, nicheDispersalLimitation = TRUE, nicheAbuDependant = FALSE)

stopImplicitCluster()

#store results
modelresults_Nsp25 <- list()
modelresults_Nsp25$varModels <- VarModelresults$varModels
modelresults_Nsp25$staticModels <- StaticModelresults$staticModels
modelresults_Nsp25$info <- VarModelresults$info

save(modelresults_Nsp25, file = "modelresults_Nsp25.RData")

#### Model with 100 SPECIES and z = 0.3 ####
Nspecies = 100 #species richness

#add species information
populatedLandscape_Nsp100 <<- populateLandscapeWithSpecies(landscapeWithStorage, Nspecies)

Tmax <- 200
#calculate and append fitness tables for every species
finalVaryingLandscape_Nsp100 <<- addVaryingFitness(populatedLandscape_Nsp100, varyingLandscapesArray)
#calculate and append true mean E value of each patch
finalVaryingLandscape_Nsp100 <<- getTrueMeanVariation(populatedLandscape_Nsp100, finalVaryingLandscape_Nsp100)

Tmax <- 2000

if(isParallelised == TRUE){
  # Calculate the number of cores
  no_cores <- detectCores()-1
  #no_cores <- 1
  # Initiate cluster
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
}

#run static and dynamic models
VarModelresults <- launchVarModels(percentVariations, frequencies, amplitudes, storageReductions, modelTypes, populatedLandscape_Nsp100, finalVaryingLandscape_Nsp100, nicheDispersalLimitation = TRUE, nicheAbuDependant = FALSE)
StaticModelresults <- launchStaticModels(percentVariations, frequencies, amplitudes, storageReductions, modelTypes, populatedLandscape_Nsp100, finalVaryingLandscape_Nsp100, nicheDispersalLimitation = TRUE, nicheAbuDependant = FALSE)

stopImplicitCluster()

#store results
modelresults_Nsp100 <- list()
modelresults_Nsp100$varModels <- VarModelresults$varModels
modelresults_Nsp100$staticModels <- StaticModelresults$staticModels
modelresults_Nsp100$info <- VarModelresults$info

save(modelresults_Nsp100, file = "modelresults_Nsp100.RData")


##### NICHE-NEUTRAL METRIC CALCULATIONS ####
####___________________________________####
#requires relevant modelResults from launchStaticModels() and launchVarModels()

####SPECIES VARIATION ACROSS REPLICATES (SVAR) ####
finalSimpsonTable <- rbind(convertSimpson(modelresults_Nsp25, speciesN = 25, z=0.3),
                           convertSimpson(modelresults_Nsp100, speciesN = 100, z=0.3)
)
finalSimpsonTable_neutral <- rbind(convertSimpson(modelresults_Nsp25, speciesN = 25, z=0.3, model = 2),
                                   convertSimpson(modelresults_Nsp100, speciesN = 100, z=0.3, model = 2)
)

neutralTable <- finalSimpsonTable_neutral
continuumTable <- finalSimpsonTable
SVARTables <- correctTableTerms(neutralTable =  neutralTable, continuumTable = continuumTable,standardise = FALSE, inverse = FALSE)

####EVENNESS####
finalSADTable <- rbind(
                       convertSADs(modelresults_Nsp25, speciesN = 25, z=0.3),
                       convertSADs(modelresults_Nsp100, speciesN = 100, z=0.3)
)
finalSADTable_neutral <- rbind(
                               convertSADs(modelresults_Nsp25, speciesN = 25, z=0.3, model = 2),
                               convertSADs(modelresults_Nsp100, speciesN = 100, z=0.3, model = 2)
)

neutralTable <- finalSADTable_neutral
continuumTable <- finalSADTable
SpeciesAbundanceDistributionTables <- correctTableTerms(neutralTable = neutralTable, continuumTable = continuumTable,standardise = FALSE, inverse = FALSE)


####ABUNDANCE BETA-NULL DEVIATION####
#calculate beta-null deviations
#(can take a very long time)
abuBetaDev_Nsp25 <- abuBetaDev(modelresults_Nsp25)
abuBetaDev_Nsp100 <- abuBetaDev(modelresults_Nsp100)

#convert into tables
finalAbuBetaTable_continuum <- rbind(
                                     convertAbuBetaDev(abuBetaDev_z03_Nsp25, speciesN = 25, z=0.3),
                                     convertAbuBetaDev(abuBetaDev_z03_Nsp100, speciesN = 100, z=0.3)
)
finalAbuBetaTable_neutral <- rbind(
                                   convertAbuBetaDev(abuBetaDev_z03_Nsp25, speciesN = 25, z=0.3, model = 2),
                                   convertAbuBetaDev(abuBetaDev_z03_Nsp100, speciesN = 100, z=0.3, model = 2)
)

neutralTable <- finalAbuBetaTable_neutral
continuumTable <- finalAbuBetaTable_continuum
abundanceBetaNullDeviationTable <- correctTableTerms(neutralTable = neutralTable,continuumTable = continuumTable,standardise = FALSE, inverse = FALSE)

#### VARIATION PARTITIONING ####
#first conversion gives table and MEM.select to use in all consequent conversions
VP_cont_Nsp100 <- convertVariationPartitioning(modelresults_Nsp100, speciesN = 100, z=0.3,baseLandscape = baseLandscape, finalVaryingLandscape = finalVaryingLandscape_Nsp100)
VP_cont_Nsp25 <- convertVariationPartitioning(modelresults_Nsp25, speciesN = 25, z=0.3,baseLandscape = baseLandscape, finalVaryingLandscape = finalVaryingLandscape_Nsp25)

VP_neutral_Nsp100 <- convertVariationPartitioning(modelresults_Nsp100, speciesN = 100, z=0.3,baseLandscape = baseLandscape, finalVaryingLandscape = finalVaryingLandscape_Nsp100, model = 2)
VP_neutral_Nsp25 <- convertVariationPartitioning(modelresults_Nsp25, speciesN = 25, z=0.3,baseLandscape = baseLandscape, finalVaryingLandscape = finalVaryingLandscape_Nsp25, model = 2)


#get index values (spc - env)
finalVarPartTable <- rbind(VP_cont_Nsp100[[1]][[1]],
                           VP_cont_Nsp25[[1]][[1]]
)
finalVarPartTable_neutral <- rbind(VP_neutral_Nsp100[[1]][[1]],
                                   VP_neutral_Nsp25[[1]][[1]]
)

neutralTable <- finalVarPartTable_neutral
continuumTable <- finalVarPartTable
variationPartitioningTable <- correctTableTerms(neutralTable = neutralTable,continuumTable = continuumTable,standardise = FALSE, inverse = FALSE)
variationPartitioningTable$value <- (variationPartitioningTable$value - 1 )/ -2 # correct index to go from -1:1 to 0:1. -1 becomes 1, 0 becomes 0.5

