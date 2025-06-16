#######################################################################################
## SPATIAL CLOSE-KIN MARK-RECAPTURE CODE FOR MOSQUITOES:                             ##
##                                                                                   ##
## This code calculates the kinship probabilities and likelihood of observed close-  ##
## kin pairs in a spatially-structured population considering parent-offspring and   ##
## full-sibling pairs where each individual may be a larva or adult. The code then   ##
## uses an optimization algorithm to find the dispersal parameters - here, the mean  ##
## daily dispersal distance (1 / lambda_d), and the strength of a barrier to         ##
## movement (delta) - that maximize the likelihood given the data.                   ##
##                                                                                   ##
## Code written by: John Marshall: john.marshall@berkeley.edu                        ##
##                  Shuyi Yang: shuyiyang@berkeley.edu                               ##
## Date: June 14th, 2025                                                             ##
## Reference: Marshall JM, Yang S, Bennett JB, Filipovic I, Rasic G (2025) Spatial   ##
## close-kin mark-recapture methods to estimate dispersal parameters and barrier     ##
## strength for mosquitoes. bioRxiv doi: https://doi.org/10.1101/2025.05.11.653364   ##
#######################################################################################

# Load required libraries:
library(optimx)
library(dplyr)
library(tidyr)
library(doParallel)

# Clear stored variables from the environment:
rm(list=ls())

# Directory where the data is:
setwd("...")

# Initialize the number of sampled individuals of each life stage:
numSampledLarvae <- 0
numSampledAdultMales <- 0
numSampledAdultFemales <- 0
numSampledAdults <- 0

# helper function used in vectorization
calc_numFullSib <- function(a, b) {
  if (length(a)>0) {
    return (a[b])
  } else {
    return (0)
  }
}

# Load data corresponding to sampled individuals of each life stage:
if (file.exists("cut_L.csv")) { 
  sampledLarvae <- read.csv("cut_L.csv") 
  numSampledLarvae <- length(sampledLarvae[,1])
  sampledPatchesLarvae <- sort(unique(sampledLarvae$Patch))
  numSampledPatchesLarvae <- length(sampledPatchesLarvae)
}

if (file.exists("cut_M.csv")) { 
  sampledAdultMales <- read.csv("cut_M.csv") 
  numSampledAdultMales <- length(sampledAdultMales[,1])
  sampledPatchesAdultMales <- sort(unique(sampledAdultMales$Patch))
  numSampledPatchesAdultMales <- length(sampledPatchesAdultMales)
}

if (file.exists("cut_F.csv")) { 
  sampledAdultFemales <- read.csv("cut_F.csv") 
  numSampledAdultFemales <- length(sampledAdultFemales[,1])
  sampledPatchesAdultFemales <- sort(unique(sampledAdultFemales$Patch))
  numSampledPatchesAdultFemales <- length(sampledPatchesAdultFemales)
}

numSampledAdults <- numSampledAdultMales + numSampledAdultFemales

# Concatenate adult female & adult male data into one adult data set:
if (numSampledAdultMales > 0) {
  sampledAdultMales$Mate <- NA
}
if ((numSampledAdultMales > 0) && (numSampledAdultFemales > 0)) {
  sampledAdults <- rbind(sampledAdultFemales,sampledAdultMales)
}
if ((numSampledAdultMales > 0) && (numSampledAdultFemales == 0)) {
  sampledAdults <- sampledAdultMales
}
if ((numSampledAdultMales == 0) && (numSampledAdultFemales > 0)) {
  sampledAdults <- sampledAdultFemales
}

if (numSampledAdults > 0) {
  sampledPatchesAdults <- sort(unique(sampledAdults$Patch))
  numSampledPatchesAdults <- length(sampledPatchesAdults)
}

# Record start & stop time of sampling of mosquitoes of any life stage:
if ((numSampledAdults > 0) && (numSampledLarvae == 0)) {
  tSamplingStart <- min(sampledAdults$Time)
  tSamplingEnd <- max(sampledAdults$Time)
  sampledPatches <- sort(unique(sampledPatchesAdults))
  numSampledPatches <- length(sampledPatches)
}
if ((numSampledAdults == 0) && (numSampledLarvae > 0)) {
  tSamplingStart <- min(sampledLarvae$Time)
  tSamplingEnd <- max(sampledLarvae$Time)
  sampledPatches <- sort(unique(sampledPatchesLarvae))
  numSampledPatches <- length(sampledPatches)
}
if ((numSampledAdults > 0) && (numSampledLarvae > 0)) {
  tSamplingStart <- min(min(sampledAdults$Time), min(sampledLarvae$Time))
  tSamplingEnd <- max(max(sampledAdults$Time), max(sampledLarvae$Time))
  sampledPatches <- sort(unique(c(sampledPatchesLarvae,
                                  sampledPatchesAdults)))
  numSampledPatches <- length(sampledPatches)
}

# Time period over which likelihood calculations will be performed:
samplingDays <- seq(tSamplingStart, tSamplingEnd, by=1)

# Total number of days over which sampling occurred:
numSamplingDays <- tSamplingEnd - tSamplingStart + 1

#######################################################################################
## MOTHER-OFFSPRING (LARVA) KINSHIP PROBABILITIES & LIKELIHOOD:                      ##
#######################################################################################

logLike_MOL <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                        beta, mu_E, mu_L, mu_P, 
                        lambda, barrier_strength, numPatches, distMat) {
  
  #####################################################################################
  ## CALCULATE TRANSITION PROBABILITIES GIVEN DISPERSAL PARAMETERS:                  ##
  #####################################################################################
  
  # First, calculate the daily transition matrix, which describes the probability that,
  # given a mosquito is at node i on one day, it is at node j on the next. Assuming a
  # mean daily dispersal distance, 1/lambda, and a barrier to movement of strength,
  # delta, we have:

  # Initialize an empty square matrix M (size: numPatches x numPatches) 
  # This will store movement probabilities between patches.
  M <- matrix(rep(0, numPatches^2), nrow = numPatches, byrow = TRUE)
  
  # Loop over each origin patch i.
  for (i in 1:numPatches) {
    
    # Initialize the denominator (normalization constant) for row i.
    Denominator <- 0
    
    # First loop: calculate the denominator as the sum of movement weights from patch i to all j.
    for (j in 1:numPatches) {
      
      # Check if patches i and j are on the same side of the barrier.
      # If both are on the same side (≤ barrier_node or > barrier_node):
      if (((i <= barrier_node) && (j <= barrier_node)) || 
          ((i > barrier_node) && (j > barrier_node))) {
        
        # Movement weight without barrier penalty.
        Denominator <- Denominator + exp(-lambda * distMat[i, j])
        
      } else {
        # If crossing the barrier, apply the barrier_strength penalty.
        Denominator <- Denominator + (1 - barrier_strength) * exp(-lambda * distMat[i, j])
      }
    }
    
    # Second loop: assign normalized movement probabilities from i to j.
    for (j in 1:numPatches) {
      
      if (((i <= barrier_node) && (j <= barrier_node)) || 
          ((i > barrier_node) && (j > barrier_node))) {
        
        # Same-side movement probability: weight divided by the total.
        M[i, j] <- exp(-lambda * distMat[i, j]) / Denominator
        
      } else {
        # Cross-barrier movement probability: penalized weight divided by the total.
        M[i, j] <- (1 - barrier_strength) * exp(-lambda * distMat[i, j]) / Denominator
      }
    }
  }
  
  # Next, calculate the transition matrices for 1, ..., T_A days, which describe the 
  # probability that, given a mosquito is at node i on one day, it is at node j 1, ...,
  # T_A days later. Given the daily transition matrix, M, defined earlier, we have:
  
  rho <- replicate(T_A,M)
  for (i in 2:T_A) {
    # Multiply matrices recursively to get transition probabilities over multiple days
    rho[,,i] <- rho[,,(i-1)] %*% M
  }
  
  #####################################################################################
  ## MOTHER-OFFSPRING (LARVA) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################
  
  # Given a larva sampled at location x2 and time t2, this is the probability that an
  # adult female sampled at location x1 and time t1 is their mother.
  
  # First, calculate the denominator:
  # This is the expected number of surviving larvae at location x2 time t2 from adult
  # females at any consistent time (assuming a constant population size, this is 
  # independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - (T_L-1)):(t2 - T_E)) {
    # Accumulate expected larvae contributions over time windows from female parents
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * 
                                    ((1 - mu_L)^(t2 - y2 - T_E)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of larvae at location x2 and day t2 from an adult 
  # female sampled at location x1 and time t1.
  # * By default, let t1 = 0, as the same equation will apply at all times.
  # * Earliest possible t2 is (-(T_A-1) + T_E), if the mother was caught at the end 
  #   of her life & gave birth to the offspring soon after emergence.
  # * Latest possible t2 is (T_E + (T_L-1)), if the the mother gave birth at the time 
  #   of sampling & the larva was caught at the end of its life.
  # * So we will explore (-(T_A-1) + T_E) <= t2 <= (T_E + (T_L-1))
  
  # Probability of adult surviving from 0 to (T_A-1) days:
  AdultSurvivalProbability <- rep(0, T_A)
  AdultAge <- rep(0, T_A)
  for (i in 1:T_A){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and (T_A-1) days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-(T_A-1) + T_E), (T_E + (T_L-1)), by=1)
  
  Numerator <- replicate(length(t2), matrix(rep(0, numSampledPatches^2), 
                                            nrow = numSampledPatches, byrow = TRUE))
  
  # Loop over all pairs of sampled patches to calculate the expected number 
  # of individuals from the two locations that are of mother-larval offspring kinship:
  for (x1 in 1:numSampledPatches) {
    for (x2 in 1:numSampledPatches) {
      for (i in 1:length(t2)) { 
        for (y2 in (t2[i] - T_E - (T_L-1)):(t2[i] - T_E)) {
          if ((y2 >= (t1 - (T_A-1))) && (y2 <= t1)) {
            Numerator[x1, x2, i] <- Numerator[x1, x2, i] + ((1 - mu_A)^(t1 - y2)
                                                            * ((rho[sampledPatches[x2], sampledPatches[x1], (1 + t1 - y2)]) 
                                                               / (sum(rho[, sampledPatches[x1], (1 + t1 - y2)])))
                                                            * beta * ((1 - mu_E)^T_E)
                                                            * ((1 - mu_L)^(t2[i] - y2 - T_E)))
          }
        }
      }
    }
  }
  
  # Mother-larval offspring probability:
  # Given a larva sampled at location x2 and time t2, the probability that an adult
  # female sampled at location x1 and time t1 is their mother is given by:
  
  P_MOL <- Numerator/Denominator
  
  # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
  # in general, the indices relate to times t2-t1, i.e.:
  
  t2Minust1 <- t2
  
  #####################################################################################
  ## FIND MOTHER-OFFSPRING (LARVA) PAIRS IN THE DATA:                                ##
  #####################################################################################
  
  larvaID <- 0
  motherID <- 0
  larvaSamplingTime <- 0
  motherSamplingTime <- 0
  larvaSamplingPatch <- 0
  motherSamplingPatch <- 0
  j <- 0
  
  # Loop over sampled larvae and match to mothers by ID:
  for (i in 1:numSampledLarvae) {
    if (any(sampledAdultFemales$myID == sampledLarvae[i, "momID"])) {
      j <- j + 1
      larvaID[j] <- sampledLarvae[i, "myID"] # Larva ID
      motherID[j] <- sampledLarvae[i, "momID"] # Mother ID
      larvaSamplingTime[j] <- sampledLarvae[i, "Time"] # Day larva sampled
      larvaSamplingPatch[j] <- sampledLarvae[i, "Patch"] # Patch larva sampled from
      # Day mother sampled and patch mother sampled from:
      motherSamplingTime[j] <- sampledAdultFemales$Time[sampledAdultFemales$myID==motherID[j]]
      motherSamplingPatch[j] <- sampledAdultFemales$Patch[sampledAdultFemales$myID==motherID[j]]
    }
  }
  
  # Array of mother-larval offspring pairs keeping track of: i) larva ID, ii) mother ID,
  # iii) day larva sampled, iv) patch larva sampled from, v) day mother sampled & 
  # vi) patch mother sampled from:
  MOL_Pairs_Data <- cbind(larvaID, motherID, larvaSamplingTime, larvaSamplingPatch, 
                          motherSamplingTime, motherSamplingPatch)
  
  # Record number of adult females sampled on each day at each sampling patch:
  dailySampledAdultFemales <- matrix(rep(0, numSampledPatches*numSamplingDays), 
                                     nrow = numSampledPatches, byrow = TRUE)
  
  for (x1 in 1:numSampledPatches) {
    for (i in 1:numSamplingDays) {
      # Count number of adult females sampled at patch x1 on day i
      dailySampledAdultFemales[x1,i] <- sum(sampledAdultFemales[
        which(sampledAdultFemales$Patch==sampledPatches[x1]),]$Time==samplingDays[i])
    }
  }
  
  # Record number of larvae sampled on each day at each sampling patch:
  dailySampledLarvae <- matrix(rep(0, numSampledPatches*numSamplingDays), 
                               nrow = numSampledPatches, byrow = TRUE)
  
  # Count number of larvae sampled at patch x1 on day i
  for (x1 in 1:numSampledPatches) {
    for (i in 1:numSamplingDays) {
      dailySampledLarvae[x1,i] <- sum(sampledLarvae[
        which(sampledLarvae$Patch==sampledPatches[x1]),]$Time==samplingDays[i])
    }
  }
  
  #####################################################################################
  ## CALCULATE LIKELIHOOD OF MOTHER-OFFSPRING (LARVA) PAIR DATA:                     ##
  #####################################################################################
  
  logLike <- 0
  
  for (x1 in 1:numSampledPatches) {
    for (x2 in 1:numSampledPatches) {
      
      # Calculate the loglikelihood that, for x sampled larvae at location x2 on day t2,
      # w have mothers among the y sampled adult females at location x1 on day t1:
      
      for (i in 1:length(samplingDays)) {
        
        larvaSamplingTimeI <- samplingDays[i]
        x <- dailySampledLarvae[x2,i] # Number of sampled larvae on day t2 at location x2
        
        for (j in 1:length(t2Minust1)) {
          motherSamplingTimeJ <- larvaSamplingTimeI - t2Minust1[j]
          
          if ((motherSamplingTimeJ >= tSamplingStart) & (motherSamplingTimeJ <= tSamplingEnd)) {
            # Number of sampled adult females on day t1 at location x1:
            y <- dailySampledAdultFemales[x1, which(samplingDays==motherSamplingTimeJ)]
            
            # Number of sampled larvae at location x2 on day t2 that have mothers among the 
            # sampled adult females at location x1 on day t1:
            w <- sum(MOL_Pairs_Data[,"motherSamplingTime"]
                     [((MOL_Pairs_Data[,"larvaSamplingTime"]==larvaSamplingTimeI) & 
                         ((MOL_Pairs_Data[,"larvaSamplingPatch"]==sampledPatches[x2]) &
                            (MOL_Pairs_Data[,"motherSamplingPatch"]==sampledPatches[x1])))]==motherSamplingTimeJ)
            
            if ((y > 0) & (P_MOL[x1,x2,j] > 0)) {
              # Probability that a given sampled larva on day t2 at location x2 has a mother 
              # among the y sampled adult females on day t1 at location x1:
              z <- 1 - ((1 - P_MOL[x1,x2,j])^y)
              
              # Now calculate the log likelihood that w sampled larvae on day t2 at location
              # x2 have mothers among the sampled adult females on day t1 at location x1:
              logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
            }
          }
        }
      }
    }
    print(x1/numSampledPatches)
  }
  
  -logLike
}

#######################################################################################
## MOTHER-OFFSPRING (ADULT) KINSHIP PROBABILITIES & LIKELIHOOD:                      ##
#######################################################################################

logLike_MOA <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                        beta, mu_E, mu_L, mu_P, 
                        lambda, barrier_strength, numPatches, distMat) {
  
  #####################################################################################
  ## CALCULATE TRANSITION PROBABILITIES GIVEN DISPERSAL PARAMETERS:                  ##
  #####################################################################################
  
  # First, calculate the daily transition matrix, which describes the probability that,
  # given a mosquito is at node i on one day, it is at node j on the next. Assuming a
  # mean daily dispersal distance, 1/lambda, and a barrier to movement of strength,
  # delta, we have:
  
  # Initialize an empty square matrix M (size: numPatches x numPatches) 
  # This will store movement probabilities between patches.
  M <- matrix(rep(0, numPatches^2), nrow = numPatches, byrow = TRUE)
  
  # Loop over each origin patch i.
  for (i in 1:numPatches) {
    
    # Initialize the denominator (normalization constant) for row i.
    Denominator <- 0
    
    # First loop: calculate the denominator as the sum of movement weights from patch i to all j.
    for (j in 1:numPatches) {
      
      # Check if patches i and j are on the same side of the barrier.
      # If both are on the same side (≤ barrier_node or > barrier_node):
      if (((i <= barrier_node) && (j <= barrier_node)) || 
          ((i > barrier_node) && (j > barrier_node))) {
        
        # Movement weight without barrier penalty.
        Denominator <- Denominator + exp(-lambda * distMat[i, j])
        
      } else {
        # If crossing the barrier, apply the barrier_strength penalty.
        Denominator <- Denominator + (1 - barrier_strength) * exp(-lambda * distMat[i, j])
      }
    }
    
    # Second loop: assign normalized movement probabilities from i to j.
    for (j in 1:numPatches) {
      
      if (((i <= barrier_node) && (j <= barrier_node)) || 
          ((i > barrier_node) && (j > barrier_node))) {
        
        # Same-side movement probability: weight divided by the total.
        M[i, j] <- exp(-lambda * distMat[i, j]) / Denominator
        
      } else {
        # Cross-barrier movement probability: penalized weight divided by the total.
        M[i, j] <- (1 - barrier_strength) * exp(-lambda * distMat[i, j]) / Denominator
      }
    }
  }
  
  # Next, calculate the transition matrices for 1, ..., T_A days, which describe the 
  # probability that, given a mosquito is at node i on one day, it is at node j 1, ...,
  # T_A days later. Given the daily transition matrix, M, we have:
  
  # Initialize a 3D array rho: [origin, destination, time lag]
  rho <- replicate(T_A, M)
  for (i in 2:T_A) {
    # Multiply matrices recursively to get transition probabilities over multiple days
    rho[,,i] <- rho[,,(i-1)] %*% M
  }
  
  #####################################################################################
  ## MOTHER-OFFSPRING (ADULT) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################
  
  # Given an adult sampled at location x2 and time t2, this is the probability that an
  # adult female sampled at location x1 and time t1 is their mother.
  
  # First, calculate the denominator:
  # This is the expected number of surviving adults at location x2 and time t2 from
  # adult females at any consistent time (assuming a constant population size, this is
  # independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - T_L - T_P - (T_A-1)):(t2 - T_E - T_L - T_P)) {
    # Accumulate expected adults contributions over time windows from females
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                  * ((1 - mu_P)^T_P) 
                                  * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of adults at location x2 on day t2 from an adult female
  # sampled at location x1 on day t1.
  # * By default, let t1 = 0, as the same equation will apply at all times (it's the
  #   difference between t1 & t2 that matters).
  # * Earliest possible t2 is (-(T_A-1) + T_E + T_L + T_P), if the mother was caught 
  #   at the end of her life & gave birth to the offspring soon after emergence.
  # * Latest possible t2 is (T_E + T_L + + T_P + (T_A-1)), if the the mother gave birth 
  #   at the time of sampling & the adult offspring was caught at the end of its life.
  # * So we will explore (-(T_A-1) + T_E + T_L + T_P) <= t2 <= (T_E + T_L + T_P + (T_A-1)).
  
  # Probability of adult surviving from 0 to (T_A-1) days:
  AdultSurvivalProbability <- rep(0, T_A)
  AdultAge <- rep(0, T_A)
  for (i in 1:T_A){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and (T_A-1) days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-(T_A-1) + T_E + T_L + T_P), (T_E + T_L + T_P + (T_A-1)), by=1)
  
  Numerator <- replicate(length(t2), matrix(rep(0, numSampledPatches^2), 
                                            nrow = numSampledPatches, byrow = TRUE))
  
  # Loop over all pairs of sampled patches to calculate the expected number of
  # individuals from the two locations are of mother-adult offspring kinship
  for (x1 in 1:numSampledPatches) {
    for (x2 in 1:numSampledPatches) {
      for (i in 1:length(t2)) { 
        for (y2 in (t2[i] - T_E - T_L - T_P - (T_A-1)):(t2[i] - T_E - T_L - T_P)) {
          if ((y2 >= (t1 - (T_A-1))) && (y2 <= t1)) {
            # First calculate the probability that, given a mother sampled at location
            # x1, the adult offspring is sampled at location x2 (rho_x1_x2):
            rho_xi_x1 <- rho[, sampledPatches[x1], (1 + t1 - y2)]
            sum_rho_xi_x1 <- sum(rho_xi_x1)
            rho_xi_x2 <- rho[, sampledPatches[x2], (1 + t2[i] - (y2 + T_E + T_L + T_P))]
            rho_x1_x2 <- as.numeric((rho_xi_x1 / sum_rho_xi_x1) %*% rho_xi_x2)
            
            # Then calculate the numerator:
            Numerator[x1, x2, i] <- Numerator[x1, x2, i] + ((1 - mu_A)^(t1 - y2) 
                                                            * rho_x1_x2 * beta * ((1 - mu_E)^T_E)
                                                            * ((1 - mu_L)^T_L) * ((1 - mu_P)^T_P)
                                                            * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
          }
        }
      }
    }
  }
  
  # Mother-adult offspring probability:
  # Given an adult sampled at time t2 & location x2, the probability that an adult 
  # female sampled at time t1 & location x1 is their mother is given by:
  
  P_MOA <- Numerator/Denominator
  
  # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
  # in general, the indices relate to times t2-t1, i.e.:
  
  t2Minust1 <- t2
  
  #####################################################################################
  ## FIND MOTHER-OFFSPRING (ADULT) PAIRS IN THE DATA:                                ##
  #####################################################################################
  
  adultID <- 0
  motherID <- 0
  adultSamplingTime <- 0
  motherSamplingTime <- 0
  adultSamplingPatch <- 0
  motherSamplingPatch <- 0
  j <- 0
  
  # Loop over sampled adults and match to mothers by ID
  for (i in 1:numSampledAdults) {
    if (any(sampledAdultFemales$myID == sampledAdults[i, "momID"])) {
      j <- j + 1
      adultID[j] <- sampledAdults[i, "myID"] # Adult ID
      motherID[j] <- sampledAdults[i, "momID"] # Mother ID
      adultSamplingTime[j] <- sampledAdults[i, "Time"] # Day adult sampled
      adultSamplingPatch[j] <- sampledAdults[i, "Patch"] # Patch adult sampled from
      # Day mother sampled & patch mother sampled from:
      motherSamplingTime[j] <- sampledAdultFemales$Time[sampledAdultFemales$myID==motherID[j]]
      motherSamplingPatch[j] <- sampledAdultFemales$Patch[sampledAdultFemales$myID==motherID[j]]
    }
  }
  
  # Array of mother-adult offspring pairs keeping track of: i) adult ID, ii) mother ID,
  # iii) day adult offspring sampled, iv) patch adult offspring sampled from, v) day 
  # mother sampled & vi) patch mother sampled from:
  MOA_Pairs_Data <- cbind(adultID, motherID, adultSamplingTime, adultSamplingPatch, 
                          motherSamplingTime, motherSamplingPatch)
  
  # Record number of adult females sampled on each day at each sampling patch:
  dailySampledAdultFemales <- matrix(rep(0, numSampledPatches*numSamplingDays), 
                                     nrow = numSampledPatches, byrow = TRUE)
  
  for (x1 in 1:numSampledPatches) {
    for (i in 1:numSamplingDays) {
      dailySampledAdultFemales[x1,i] <- sum(sampledAdultFemales[
        which(sampledAdultFemales$Patch==sampledPatches[x1]),]$Time==samplingDays[i])
    }
  }
  
  # Record number of adults sampled on each day at each sampling patch:
  dailySampledAdults <- matrix(rep(0, numSampledPatches*numSamplingDays), 
                               nrow = numSampledPatches, byrow = TRUE)
  
  # Count number of adults sampled at patch x1 on day i:
  for (x1 in 1:numSampledPatches) {
    for (i in 1:numSamplingDays) {
      dailySampledAdults[x1,i] <- sum(sampledAdults[
        which(sampledAdults$Patch==sampledPatches[x1]),]$Time==samplingDays[i])
    }
  }
  
  #####################################################################################
  ## CALCULATE LIKELIHOOD OF MOTHER-OFFSPRING (ADULT) PAIR DATA:                     ##
  #####################################################################################
  
  logLike <- 0
  
  for (x1 in 1:numSampledPatches) {
    for (x2 in 1:numSampledPatches) {
      
      # Calculate the loglikelihood that, for x sampled adults on day t2 at location x2,
      # w have mothers among the y sampled adult females on day t1 at location x1:
      
      for (i in 1:length(samplingDays)) {
        
        adultSamplingTimeI <- samplingDays[i]
        x <- dailySampledAdults[x2,i] # Number of sampled adults on day t2 at location x2
        
        for (j in 1:length(t2Minust1)) {
          motherSamplingTimeJ <- adultSamplingTimeI - t2Minust1[j]
          
          if ((motherSamplingTimeJ >= tSamplingStart) & (motherSamplingTimeJ <= tSamplingEnd)) {
            # Number of sampled adult females on day t1 at location x1:
            y <- dailySampledAdultFemales[x1, which(samplingDays==motherSamplingTimeJ)]
            if ((j == which(t2Minust1==0)) && (x1 == x2)) { y <- (y - 1) }
            
            # Number of sampled adults at location x2 on day t2 that have mothers among
            # the sampled adult females at x1 on day t1:
            w <- sum(MOA_Pairs_Data[,"motherSamplingTime"]
                     [((MOA_Pairs_Data[,"adultSamplingTime"]==adultSamplingTimeI) & 
                         ((MOA_Pairs_Data[,"adultSamplingPatch"]==sampledPatches[x2]) &
                            (MOA_Pairs_Data[,"motherSamplingPatch"]==sampledPatches[x1])))]==motherSamplingTimeJ)
            
            if ((y > 0) & (P_MOA[x1,x2,j] > 0)) {
              # Probability that a given sampled adult on day t2 at location x2 has a 
              # mother among the y sampled adult females on day t1 at location x1:
              z <- 1 - ((1 - P_MOA[x1,x2,j])^y)
              
              # Now calculate the log likelihood that w sampled adults on day t2 at 
              # location x2 have mothers among the sampled adult females at location x1 
              # on day t1:
              logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
            }
          }
        }
      }
    }
    print(x1/numSampledPatches)
  }
  
  -logLike
}

#######################################################################################
## FATHER-OFFSPRING (LARVA) KINSHIP PROBABILITIES & LIKELIHOOD:                      ##
#######################################################################################

logLike_FOL <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                        beta, mu_E, mu_L, mu_P, 
                        lambda, barrier_strength, numPatches, distMat) {
  
  #####################################################################################
  ## CALCULATE TRANSITION PROBABILITIES GIVEN DISPERSAL PARAMETERS:                  ##
  #####################################################################################
  
  # First, calculate the daily transition matrix, which describes the probability that,
  # given a mosquito is at node i on one day, it is at node j on the next. Assuming a
  # mean daily dispersal distance, 1/lambda, and a barrier to movement of strength,
  # delta, we have:
  
  # Initialize an empty square matrix M (size: numPatches x numPatches) 
  # This will store movement probabilities between patches.
  M <- matrix(rep(0, numPatches^2), nrow = numPatches, byrow = TRUE)
  
  # Loop over each origin patch i.
  for (i in 1:numPatches) {
    
    # Initialize the denominator (normalization constant) for row i.
    Denominator <- 0
    
    # First loop: calculate the denominator as the sum of movement weights from patch i to all j.
    for (j in 1:numPatches) {
      
      # Check if patches i and j are on the same side of the barrier.
      # If both are on the same side (≤ barrier_node or > barrier_node):
      if (((i <= barrier_node) && (j <= barrier_node)) || 
          ((i > barrier_node) && (j > barrier_node))) {
        
        # Movement weight without barrier penalty.
        Denominator <- Denominator + exp(-lambda * distMat[i, j])
        
      } else {
        # If crossing the barrier, apply the barrier_strength penalty.
        Denominator <- Denominator + (1 - barrier_strength) * exp(-lambda * distMat[i, j])
      }
    }
    
    # Second loop: assign normalized movement probabilities from i to j.
    for (j in 1:numPatches) {
      
      if (((i <= barrier_node) && (j <= barrier_node)) || 
          ((i > barrier_node) && (j > barrier_node))) {
        
        # Same-side movement probability: weight divided by the total.
        M[i, j] <- exp(-lambda * distMat[i, j]) / Denominator
        
      } else {
        # Cross-barrier movement probability: penalized weight divided by the total.
        M[i, j] <- (1 - barrier_strength) * exp(-lambda * distMat[i, j]) / Denominator
      }
    }
  }
  
  # Next, calculate the transition matrices for 1, ..., T_A days, which describe the 
  # probability that, given a mosquito is at node i on one day, it is at node j 1, ...,
  # T_A days later. Given a daily transition matrix, M, we have:
  
  rho <- replicate(T_A,M)
  for (i in 2:T_A) {
    # Multiply matrices recursively to get transition probabilities over multiple days
    rho[,,i] <- rho[,,(i-1)] %*% M
  }
  
  #####################################################################################
  ## FATHER-OFFSPRING (LARVA) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################
  
  # Given a larva sampled at time t2 and location x2, this is the probability that an 
  # adult male sampled at time t1 and location x1 is their father.
  
  # First, calculate the denominator:
  # This is the expected number of surviving larvae at time t2 from adult females at any
  # consistent time (assuming a constant population size, this is independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - (T_L-1)):(t2 - T_E)) {
    # Accumulate expected larvae contributions over time windows from males
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of larvae at location x2 on day t2 from an adult male
  # sampled at location x1 on day t1.
  # * By default, let t1 = 0, as the same equation will apply at all times.
  # * Earliest possible t2 is (-(T_A-1) + T_E), if the father was caught at the end of
  #   his life, but mated at the beginning of his life, & if the mother gave birth soon
  #   after mating, & the larva was caught soon after emergence.
  # * Latest possible t2 is ((T_A-1) + T_E + (T_L-1)), if the the father was caught soon
  #   after mating, the mother mated at the beginning of her life & laid eggs at the
  #   end of her life, & the larva was caught very soon before developing into a pupa.
  # * So we will explore (-(T_A-1) + T_E) <= t2 <= ((T_A-1) + T_E + (T_L-1))
  
  # Probability of adult surviving from 0 to (T_A-1) days:
  AdultSurvivalProbability <- rep(0, T_A)
  AdultAge <- rep(0, T_A)
  for (i in 1:T_A){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and (T_A-1) days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  

  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-(T_A-1) + T_E), ((T_A-1) + T_E + (T_L-1)), by=1)
  
  Numerator <- replicate(length(t2), matrix(rep(0, numSampledPatches^2), 
                                            nrow = numSampledPatches, byrow = TRUE))
  
  # Loop over all pairs of sampled patches to calculate the expected number of
  # individuals from the two locations are of father-larval offspring kinship:
  for (x1 in 1:numSampledPatches) {
    for (x2 in 1:numSampledPatches) {
      for (i in 1:length(t2)) { 
        for (tk in (t1 - (T_A-1)):t1) {
          for (y2 in tk:(tk + (T_A-1))) {
            if ((t2[i] >= (y2 + T_E)) && (t2[i] <= (y2 + T_E + (T_L-1)))) {
              # First calculate the probability that, given a father sampled at location
              # x1, the larval offspring is sampled at location x2 (rho_x1_x2):
              rho_xk_x1 <- rho[, sampledPatches[x1], (1 + t1 - tk)]
              sum_rho_xk_x1 <- sum(rho_xk_x1)
              rho_xk_x2 <- rho[, sampledPatches[x2], (1 + y2 - tk)]
              rho_x1_x2 <- as.numeric((rho_xk_x1 / sum_rho_xk_x1) %*% rho_xk_x2)
              
              # Then calculate the numerator:
              # survival * movement * fecundity * egg & larval survival
              Numerator[x1, x2, i] <- Numerator[x1, x2, i] + (AdultAgeProbability[which(AdultAge==(t1 - tk))] 
                                                              * rho_x1_x2
                                                              * ((1 - mu_A)^(y2 - tk))
                                                              * beta * ((1 - mu_E)^T_E)
                                                              * ((1 - mu_L)^(t2[i] - y2 - T_E)))
            }
          }
        }
      }
    }
  }
  
  # Father-larval offspring probability:
  # Given a larva sampled at time t2 and location x2, the probability that an adult male
  # sampled at time t1 and location x1 is their father is given by:
  
  P_FOL <- Numerator/Denominator
  
  # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
  # in general, the indices relate to times t2-t1, i.e.:
  
  t2Minust1 <- t2
  
  #####################################################################################
  ## FIND FATHER-OFFSPRING (LARVA) PAIRS IN THE DATA:                                ##
  #####################################################################################
  
  larvaID <- 0
  fatherID <- 0
  larvaSamplingTime <- 0
  fatherSamplingTime <- 0
  larvaSamplingPatch <- 0
  fatherSamplingPatch <- 0
  j <- 0
  
  # Loop over sampled larvae and match to fathers by ID
  for (i in 1:numSampledLarvae) {
    if (any(sampledAdultMales$myID == sampledLarvae[i, "dadID"])) {
      j <- j + 1
      larvaID[j] <- sampledLarvae[i, "myID"] # Larva ID
      fatherID[j] <- sampledLarvae[i, "dadID"] # Father ID
      larvaSamplingTime[j] <- sampledLarvae[i, "Time"] # Day larva sampled
      larvaSamplingPatch[j] <- sampledLarvae[i, "Patch"] # Patch larva sampled from
      # Day father sampled & patch father sampled from:
      fatherSamplingTime[j] <- sampledAdultMales$Time[sampledAdultMales$myID==fatherID[j]]
      fatherSamplingPatch[j] <- sampledAdultMales$Patch[sampledAdultMales$myID==fatherID[j]]
    }
  }
  
  # Array of father-larval offspring pairs keeping track of: i) larva ID, ii) father ID,
  # iii) day larva sampled, iv) patch larva sampled from, v) day father sampled & 
  # vi) patch father sampled from:
  FOL_Pairs_Data <- cbind(larvaID, fatherID, larvaSamplingTime, larvaSamplingPatch, 
                          fatherSamplingTime, fatherSamplingPatch)
  
  # Record number of adult females sampled on each day at each sampling patch:
  dailySampledAdultMales <- matrix(rep(0, numSampledPatches*numSamplingDays), 
                                   nrow = numSampledPatches, byrow = TRUE)
  
  for (x1 in 1:numSampledPatches) {
    for (i in 1:numSamplingDays) {
      # Count number of larvae sampled at patch x1 on day i
      dailySampledAdultMales[x1,i] <- sum(sampledAdultMales[
        which(sampledAdultMales$Patch==sampledPatches[x1]),]$Time==samplingDays[i])
    }
  }
  
  # Record number of larvae sampled on each day at each sampling patch:
  dailySampledLarvae <- matrix(rep(0, numSampledPatches*numSamplingDays), 
                               nrow = numSampledPatches, byrow = TRUE)
  
  for (x1 in 1:numSampledPatches) {
    for (i in 1:numSamplingDays) {
      dailySampledLarvae[x1,i] <- sum(sampledLarvae[
        which(sampledLarvae$Patch==sampledPatches[x1]),]$Time==samplingDays[i])
    }
  }
  
  #####################################################################################
  ## CALCULATE LIKELIHOOD OF FATHER-OFFSPRING (LARVA) PAIR DATA:                     ##
  #####################################################################################
  
  logLike <- 0
  
  for (x1 in 1:numSampledPatches) {
    for (x2 in 1:numSampledPatches) {
      
      # Calculate the loglikelihood that, for x sampled larvae at location x2 on day t2,
      # w have mothers among the y sampled adult females at location x1 on day t1:
      
      for (i in 1:length(samplingDays)) {
        
        larvaSamplingTimeI <- samplingDays[i]
        x <- dailySampledLarvae[x2,i] # Number of sampled larvae at x2 on day t2
        
        for (j in 1:length(t2Minust1)) {
          fatherSamplingTimeJ <- larvaSamplingTimeI - t2Minust1[j]
          
          if ((fatherSamplingTimeJ >= tSamplingStart) & (fatherSamplingTimeJ <= tSamplingEnd)) {
            # Number of sampled adult males at x1 on day t1:
            y <- dailySampledAdultMales[x1, which(samplingDays==fatherSamplingTimeJ)]
            
            # Number of sampled larvae at x2 on day t2 that have fathers among the sampled adult 
            # males at x1 on day t1:
            w <- sum(FOL_Pairs_Data[,"fatherSamplingTime"]
                     [((FOL_Pairs_Data[,"larvaSamplingTime"]==larvaSamplingTimeI) & 
                         ((FOL_Pairs_Data[,"larvaSamplingPatch"]==sampledPatches[x2]) &
                            (FOL_Pairs_Data[,"fatherSamplingPatch"]==sampledPatches[x1])))]==fatherSamplingTimeJ)
            
            if ((y > 0) & (P_FOL[x1,x2,j] > 0)) {
              # Probability that a given sampled larva on day t2 has a father among the y
              # sampled adult males on day t1:
              z <- 1 - ((1 - P_FOL[x1,x2,j])^y)
              
              # Now calculate the log likelihood that w sampled larvae at location x2 on 
              # day t2 have fathers among the sampled adult males at location x1 on day t1:
              logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
            }
          }
        }
      }
    }
    print(x1/numSampledPatches)
  }
  
  -logLike
  
}

#######################################################################################
## FATHER-OFFSPRING (ADULT) KINSHIP PROBABILITIES & LIKELIHOOD:                      ##
#######################################################################################

logLike_FOA <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                        beta, mu_E, mu_L, mu_P, 
                        lambda, barrier_strength, numPatches, distMat) {
  
  #####################################################################################
  ## CALCULATE TRANSITION PROBABILITIES GIVEN DISPERSAL PARAMETERS:                  ##
  #####################################################################################
  
  # First, calculate the daily transition matrix, which describes the probability that,
  # given a mosquito is at node i on one day, it is at node j on the next. Assuming a
  # mean daily dispersal distance, 1/lambda, and a barrier to movement of strength,
  # delta, we have:
  
  # Initialize an empty square matrix M (size: numPatches x numPatches) 
  # This will store movement probabilities between patches.
  M <- matrix(rep(0, numPatches^2), nrow = numPatches, byrow = TRUE)
  
  # Loop over each origin patch i.
  for (i in 1:numPatches) {
    
    # Initialize the denominator (normalization constant) for row i.
    Denominator <- 0
    
    # First loop: calculate the denominator as the sum of movement weights from patch i to all j.
    for (j in 1:numPatches) {
      
      # Check if patches i and j are on the same side of the barrier.
      # If both are on the same side (≤ barrier_node or > barrier_node):
      if (((i <= barrier_node) && (j <= barrier_node)) || 
          ((i > barrier_node) && (j > barrier_node))) {
        
        # Movement weight without barrier penalty.
        Denominator <- Denominator + exp(-lambda * distMat[i, j])
        
      } else {
        # If crossing the barrier, apply the barrier_strength penalty.
        Denominator <- Denominator + (1 - barrier_strength) * exp(-lambda * distMat[i, j])
      }
    }
    
    # Second loop: assign normalized movement probabilities from i to j.
    for (j in 1:numPatches) {
      
      if (((i <= barrier_node) && (j <= barrier_node)) || 
          ((i > barrier_node) && (j > barrier_node))) {
        
        # Same-side movement probability: weight divided by the total.
        M[i, j] <- exp(-lambda * distMat[i, j]) / Denominator
        
      } else {
        # Cross-barrier movement probability: penalized weight divided by the total.
        M[i, j] <- (1 - barrier_strength) * exp(-lambda * distMat[i, j]) / Denominator
      }
    }
  }
  
  # Next, calculate the transition matrices for 1, ..., (2*T_A) days, which describe the
  # probability that, given a mosquito is at node i on one day, it is at node j 1, ...,
  # (2*T_A) days later. Given the single-day transition matrix, M, we have:
  
  rho <- replicate((2*T_A),M)
  for (i in 2:(2*T_A)) {
    # Multiply matrices recursively to get transition probabilities over multiple days
    rho[,,i] <- rho[,,(i-1)] %*% M
  }
  
  #####################################################################################
  ## FATHER-OFFSPRING (ADULT) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################
  
  # Given an adult sampled at location x2 and time t2, this is the probability that an 
  # adult male sampled at location x1 and time t1 is their father.
  
  # First, calculate the denominator:
  # This is the expected number of surviving adults at location x2 and time t2 from 
  # adult females at any consistent time (assuming a constant population size, this is
  # independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - T_L - T_P - (T_A-1)):(t2 - T_E - T_L - T_P)) {
    # Accumulate expected adults contributions over time windows from males
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                  * ((1 - mu_P)^T_P) 
                                  * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of adult offspring at location x2 on day t2 from an 
  # adult male sampled at location x1 on day t1.
  # * By default, let t1 = 0, as the same equation will apply at all times.
  # * Earliest possible t2 is (-(T_A-1) + T_E + T_L + T_P), if the father was caught at
  #   the end of his life, but mated at the beginning of his life, & if the mother gave
  #   birth soon after mating, & the adult offspring was caught soon after emergence.
  # * Latest possible t2 is ((T_A-1) + T_E + T_L + T_P + (T_A-1)), if the the father was caught 
  #   soon after mating, the mother mated at the beginning of her life & laid eggs at the
  #   end of her life, & the adult offspring was caught at the end of its life.
  # * So we will explore (-(T_A-1) + T_E + T_L + T_P) <= t2 <= 
  #   ((T_A-1) + T_E + T_L + T_P + (T_A-1))
  
  # Probability of adult surviving from 0 to (T_A-1) days:
  AdultSurvivalProbability <- rep(0, T_A)
  AdultAge <- rep(0, T_A)
  for (i in 1:T_A){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and (T_A-1) days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-(T_A-1) + T_E + T_L + T_P), ((T_A-1) + T_E + T_L + T_P + (T_A-1)), by=1)
  
  Numerator <- replicate(length(t2), matrix(rep(0, numSampledPatches^2), 
                                            nrow = numSampledPatches, byrow = TRUE))
  
  # Loop over all pairs of sampled patches to calculate the expected number of
  # individuals from the two locations are of the father-adult offspring kinship
  for (x1 in 1:numSampledPatches) {
    for (x2 in 1:numSampledPatches) {
      for (i in 1:length(t2)) { 
        for (tk in (t1 - (T_A-1)):t1) {
          for (y2 in tk:(tk + (T_A-1))) {
            if ((t2[i] >= (y2 + T_E + T_L + T_P)) && (t2[i] <= (y2 + T_E + T_L + T_P + (T_A-1)))) {
              # First calculate the probability that, given a father sampled at location
              # x1, the adult offspring is sampled at location x2 (rho_x1_x2):
              rho_xk_x1 <- rho[, sampledPatches[x1], (1 + t1 - tk)]
              sum_rho_xk_x1 <- sum(rho_xk_x1)
              rho_xk_x2 <- rho[, sampledPatches[x2], ((1 + y2 - tk) + (1 + t2[i] - (y2 + T_E + T_L + T_P)))]
              rho_x1_x2 <- as.numeric((rho_xk_x1 / sum_rho_xk_x1) %*% rho_xk_x2)
              
              # Then calculate the numerator:
              # Add contribution: survival * movement * fecundity * egg & larval survival
              Numerator[x1, x2, i] <- Numerator[x1, x2, i] + (AdultAgeProbability[which(AdultAge==(t1 - tk))] 
                                                              * ((1 - mu_A)^(y2 - tk)) * rho_x1_x2
                                                              * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                                              * ((1 - mu_P)^T_P)
                                                              * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
            }
          }
        }
      }
    }
    print(c("CKMR probability: ", x1/numSampledPatches))
  }
  
  # Father-adult offspring probability:
  # Given an adult sampled at location x2 on day t2, the probability that an adult male
  # sampled at location x1 on day t1 is their father is given by:
  
  P_FOA <- Numerator/Denominator
  
  # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
  # in general, the indices relate to times t2-t1, i.e.:
  
  t2Minust1 <- t2
  
  #####################################################################################
  ## FIND FATHER-OFFSPRING (ADULT) PAIRS IN THE DATA:                                ##
  #####################################################################################
  
  adultID <- 0
  fatherID <- 0
  adultSamplingTime <- 0
  fatherSamplingTime <- 0
  adultSamplingPatch <- 0
  fatherSamplingPatch <- 0
  j <- 0
  
  
  # Loop over sampled adults and match to fathers by ID
  for (i in 1:numSampledAdults) {
    if (any(sampledAdultMales$myID == sampledAdults[i, "dadID"])) {
      j <- j + 1
      adultID[j] <- sampledAdults[i, "myID"] # Adult ID
      fatherID[j] <- sampledAdults[i, "dadID"] # Father ID
      adultSamplingTime[j] <- sampledAdults[i, "Time"] # Day adult sampled
      adultSamplingPatch[j] <- sampledAdults[i, "Patch"] # Patch adult sampled from
      # Day father sampled & patch father sampled from:
      fatherSamplingTime[j] <- sampledAdultMales$Time[sampledAdultMales$myID==fatherID[j]]
      fatherSamplingPatch[j] <- sampledAdultMales$Patch[sampledAdultMales$myID==fatherID[j]]
    }
  }
  
  # Array of father-adult offspring pairs keeping track of: i) adult ID, ii) father ID,
  # iii) day adult sampled, iv) patch adult sampled from, v) day father sampled & 
  # vi) patch father sampled from:
  FOA_Pairs_Data <- cbind(adultID, fatherID, adultSamplingTime, adultSamplingPatch, 
                          fatherSamplingTime, fatherSamplingPatch)
  
  # Record number of adult males sampled on each day at each sampling patch:
  dailySampledAdultMales <- matrix(rep(0, numSampledPatches*numSamplingDays), 
                                   nrow = numSampledPatches, byrow = TRUE)
  
  for (x1 in 1:numSampledPatches) {
    for (i in 1:numSamplingDays) {
      dailySampledAdultMales[x1,i] <- sum(sampledAdultMales[
        which(sampledAdultMales$Patch==sampledPatches[x1]),]$Time==samplingDays[i])
    }
  }
  
  # Record number of adults sampled on each day at each sampling patch:
  dailySampledAdults <- matrix(rep(0, numSampledPatches*numSamplingDays), 
                               nrow = numSampledPatches, byrow = TRUE)
  
  for (x1 in 1:numSampledPatches) {
    for (i in 1:numSamplingDays) {
      dailySampledAdults[x1,i] <- sum(sampledAdults[
        which(sampledAdults$Patch==sampledPatches[x1]),]$Time==samplingDays[i])
    }
  }
  
  #####################################################################################
  ## CALCULATE LIKELIHOOD OF FATHER-OFFSPRING (ADULT) PAIR DATA:                     ##
  #####################################################################################
  
  logLike <- 0
  
  for (x1 in 1:numSampledPatches) {
    for (x2 in 1:numSampledPatches) {
      
      # Calculate the loglikelihood that, for x sampled adults at location x2 on day
      # t2, w have fathers among the y sampled adult males at location x1 on day t1:
      
      for (i in 1:length(samplingDays)) {
        
        adultSamplingTimeI <- samplingDays[i]
        x <- dailySampledAdults[x2,i] # Number of sampled adults at location x2 on day t2
        
        for (j in 1:length(t2Minust1)) {
          fatherSamplingTimeJ <- adultSamplingTimeI - t2Minust1[j]
          
          if ((fatherSamplingTimeJ >= tSamplingStart) & (fatherSamplingTimeJ <= tSamplingEnd)) {
            # Number of sampled adult males at location x1 on day t1:
            y <- dailySampledAdultMales[x1, which(samplingDays==fatherSamplingTimeJ)]
            if ((j == which(t2Minust1==0)) && (x1 == x2)) { y <- (y - 1) }
            
            # Number of sampled adults at x2 on day t2 that have fathers among the sampled adult 
            # males at x1 on day t1:
            w <- sum(FOA_Pairs_Data[,"fatherSamplingTime"]
                     [((FOA_Pairs_Data[,"adultSamplingTime"]==adultSamplingTimeI) & 
                         ((FOA_Pairs_Data[,"adultSamplingPatch"]==sampledPatches[x2]) &
                            (FOA_Pairs_Data[,"fatherSamplingPatch"]==sampledPatches[x1])))]==fatherSamplingTimeJ)
            
            if ((y > 0) & (P_FOA[x1,x2,j] > 0)) {
              # Probability that a given sampled adult at location x2 on day t2 has a 
              # father among the y sampled adult males at location x1 on day t1:
              z <- 1 - ((1 - P_FOA[x1,x2,j])^y)
              
              # Now calculate the log likelihood that w sampled adults at location x2 on day
              # t2 have fathers among the sampled adult males at location x1 on day t1:
              logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
            }
          }
        }
      }
    }
    print(c("Log likelihood: ", x1/numSampledPatches))
  }
  
  -logLike
}

#######################################################################################
## FULL-SIBLING KINSHIP PROBABILITIES & LIKELIHOOD:                                  ##
#######################################################################################

logLike_Sibs <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                         beta, mu_E, mu_L, mu_P, 
                         lambda, barrier_strength, numPatches, distMat) {
  
  #####################################################################################
  ## CALCULATE TRANSITION PROBABILITIES GIVEN DISPERSAL PARAMETERS:                  ##
  #####################################################################################
  
  # First, calculate the daily transition matrix, which describes the probability that,
  # given a mosquito is at node i on one day, it is at node j on the next. Assuming a
  # mean daily dispersal distance, 1/lambda, and a barrier to movement of strength,
  # delta, we have:
  
  # Initialize an empty square matrix M (size: numPatches x numPatches) 
  # This will store movement probabilities between patches.
  M <- matrix(rep(0, numPatches^2), nrow = numPatches, byrow = TRUE)
  
  # Loop over each origin patch i.
  for (i in 1:numPatches) {
    
    # Initialize the denominator (normalization constant) for row i.
    Denominator <- 0
    
    # First loop: calculate the denominator as the sum of movement weights from patch i to all j.
    for (j in 1:numPatches) {
      
      # Check if patches i and j are on the same side of the barrier.
      # If both are on the same side (≤ barrier_node or > barrier_node):
      if (((i <= barrier_node) && (j <= barrier_node)) || 
          ((i > barrier_node) && (j > barrier_node))) {
        
        # Movement weight without barrier penalty.
        Denominator <- Denominator + exp(-lambda * distMat[i, j])
        
      } else {
        # If crossing the barrier, apply the barrier_strength penalty.
        Denominator <- Denominator + (1 - barrier_strength) * exp(-lambda * distMat[i, j])
      }
    }
    
    # Second loop: assign normalized movement probabilities from i to j.
    for (j in 1:numPatches) {
      
      if (((i <= barrier_node) && (j <= barrier_node)) || 
          ((i > barrier_node) && (j > barrier_node))) {
        
        # Same-side movement probability: weight divided by the total.
        M[i, j] <- exp(-lambda * distMat[i, j]) / Denominator
        
      } else {
        # Cross-barrier movement probability: penalized weight divided by the total.
        M[i, j] <- (1 - barrier_strength) * exp(-lambda * distMat[i, j]) / Denominator
      }
    }
  }
  
  # Next, calculate the transition matrices for 1, ..., T_A days or (2*T_A) days ((2*T_A) 
  # days is for the case where adults are sampled, to account for the possibility that the 
  # sampled adults may move in addition to their mother). This describes the probability 
  # that, given a mosquito is at node i on one day, it is at node j 1, ..., T_A days 
  # later. Given the daily transition probability matrix, M, we have:
  
  if (numSampledAdults == 0) {
    rho <- replicate(T_A,M)
    for (i in 2:T_A) {
      rho[,,i] <- rho[,,(i-1)] %*% M
    }
  }
  if (numSampledAdults > 0) {
    rho <- replicate((2*T_A),M)
    for (i in 2:(2*T_A)) {
      rho[,,i] <- rho[,,(i-1)] %*% M
    }
  }
  
  #####################################################################################
  ## FULL-SIBLING (LARVA-LARVA) KINSHIP PROBABILITIES:                               ##
  #####################################################################################
  
  # Given a larva sampled at time t1 and location x1, this is the probability that a 
  # larva sampled at time t2 and location x2 is its full sibling.
  
  # Probability of larva surviving from 0 to (T_L-1) days:
  LarvaSurvivalProbability <- rep(0, T_L)
  LarvaAge <- rep(0, T_L)
  for (i in 1:T_L) {
    LarvaSurvivalProbability[i] <- (1 - mu_L)^(i-1)
    LarvaAge[i] <- i - 1
  }
  
  # Probability of larva having an age between 0 and (T_L-1) days:
  LarvaAgeProbability <- LarvaSurvivalProbability / sum(LarvaSurvivalProbability)
  
  # Probability of adult surviving from 0 to (T_A-1) days:
  AdultSurvivalProbability <- rep(0, T_A)
  AdultAge <- rep(0, T_A)
  for (i in 1:T_A) {
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1) # adult survival by age
    AdultAge[i] <- i - 1 # age index, starts at 0
  }
  
  # Probability of adult having an age between 0 and (T_A-1) days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  # From here, the function calculates kinship probabilities for all pair types (larva-larva,
  # adult-adult, larva-adult, etc.) and handles them in separate conditional blocks.
  # Each block calculates:
  #   - Denominator: expected number of offspring or surviving individuals under neutral expectation
  #   - Numerator: expected number of kin pairs matching the observed type and timing
  #   - Ratio: kinship probability for that pair type
  
  if (numSampledLarvae > 0) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving larvae at location x2 and time t2 from adult
    # females at any consistent time (assuming a constant population size, this is 
    # independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - (T_L-1)):(t2 - T_E)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of larvae at location x2 and day t2 that are full
    # siblings of a larva sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (-T_A - T_L), if the mother laid egg 1 at the end of her
    #   life, egg 2 at the beginning of her life, larva 1 was caught at the end of its
    #   life, & larva 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_A + T_L), if the mother laid egg 1 at the beginning of her
    #   life, egg 2 at the end of her life, larva 1 was caught at the beginning of its
    #   life, & larva 2 was caught at the end of its life.
    # * So we will explore (-T_A - T_L) <= t2 <= (T_A + T_L)
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((-(T_A-1) - (T_L-1)), ((T_A-1) + (T_L-1)), by=1)
    
    Numerator <- replicate(length(t2), matrix(rep(0, numSampledPatches^2), 
                                              nrow = numSampledPatches, byrow = TRUE))
    
    print("Calculating full-sibling larva-larva probabilities:")
    
    # Loop over possible larva 2 sampling times
    for (i in 1:length(t2)) { 
      # y1: maternal birth time for larva 1 (sampled at t1)
      for (y1 in (t1 - T_E - (T_L-1)):(t1 - T_E)) {
        # y2: maternal birth time for larva 2 (sampled at t2[i])
        for (y2 in (y1 - (T_A-1)):(y1 + (T_A-1))) {
          # Check if the time offsets are compatible for two siblings.
          if ((y2 >= (t2[i] - T_E - (T_L-1))) && (y2 <= (t2[i] - T_E))) {
            # Compute the contribution to the numerator:
            # - Larva 1 age probability × half-sibling factor (1/2 for full sib) ×
            # - mother’s survival between laying events × birth rate × 
            # - embryo survival × larva 2 survival.
            NumeratorTemp <- (LarvaAgeProbability[which(LarvaAge==(t1 - y1 - T_E))]
                              * (1/2) * (1 - mu_A)^(abs(y2 - y1)) 
                              * beta * ((1 - mu_E)^T_E)
                              * ((1 - mu_L)^(t2[i] - y2 - T_E)))
            for (x1 in 1:numSampledPatches) {
              for (x2 in 1:numSampledPatches) {
                # First calculate the probability that, given a larva sampled at location
                # x1, a larva sampled at location x2 is its sibling (rho_x1_x2):
                rho_x1_x2 <- 0
                if (y2 == y1) {
                  # If birth times equal, siblings are in same patch
                  rho_x1_x2 <- (x1==x2)*1
                }
                if (y2 > y1) {
                  # If y2 later, use forward movement over y2 - y1 days
                  rho_x1_x2 <- rho[sampledPatches[x1],sampledPatches[x2], (y2 - y1)]
                }
                if (y1 > y2) {
                  # If y1 later, adjust for backward movement and normalize
                  rho_x1_x2 <- ((rho[sampledPatches[x2],sampledPatches[x1], (y1 - y2)]) / 
                                  (sum(rho[,sampledPatches[x1], (y1 - y2)])))
                }
                
                # Then calculate the numerator:
                Numerator[x1, x2, i] <- Numerator[x1, x2, i] + NumeratorTemp * rho_x1_x2
              }
            }
          }
        }
      }
      print(i/length(t2))
    }
    
    # Full sibling larva-larva probability:
    # Given a larva sampled at location x2 and time t2, the probability that a larva 
    # sampled at location x1 and time t1 is its full sibling is given by:
    
    P_FSLL <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSLL <- seq((-(T_A-1) - (T_L-1)), ((T_A-1) + (T_L-1)), by=1)
  }
  
  #####################################################################################
  ## FULL-SIBLING (ADULT-ADULT) KINSHIP PROBABILITIES:                               ##
  #####################################################################################
  
  # Given an adult sampled at location x2 and time t2, this is the probability that an
  # adult sampled at location x1 and time t1 is its full sibling.
  
  if (numSampledAdults > 0) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving adults at location x2 and time t2 from
    # adult females at any consistent time (assuming a constant population size, this
    # is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P - (T_A-1)):(t2 - T_E - T_L - T_P)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^T_P) 
                                    * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of adults at location x2 and  day t2 that are full 
    # siblings of an adult sampled at location x1 and time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (-T_A - T_A), if the mother laid egg 1 at the end of her
    #   life, egg 2 at the beginning of her life, adult 1 was caught at the end of its
    #   life, & adult 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_A + T_A), if the mother laid egg 1 at the beginning of her
    #   life, egg 2 at the end of her life, adult 1 was caught at the beginning of its
    #   life, & adult 2 was caught at the end of its life.
    # * So we will explore (-T_A - T_A) <= t2 <= (T_A + T_A)
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((-(T_A-1) - (T_A-1)), ((T_A-1) + (T_A-1)), by=1)
    
    Numerator <- replicate(length(t2), matrix(rep(0, numSampledPatches^2), 
                                              nrow = numSampledPatches, byrow = TRUE))
    
    print("Calculating full-sibling adult-adult probabilities:")
    
    i <- 1:length(t2)
    y1 <- (t1 - T_E - T_L - T_P - (T_A-1)):(t1 - T_E - T_L - T_P)# possible birth times for t1 adults
    
    # Create a dataframe pairing each i (t2 index) with y1 birth times
    data = as.data.frame(cbind(i=rep(i, each=length(y1)), y1 = rep(y1, length(i))))
    
    # For each y1, generate corresponding y2 ranges
    data$y2 <-  lapply(data$y1, function(input) {return((input - (T_A-1)):(input + (T_A-1)))})
    data <- data %>% unnest(y2)
    # Add the specific t2 value for each row
    data$t2i <- t2[data$i]
    # Filter combinations satisfying time consistency (i.e., compatible with t2 event)
    data <- data %>% filter((y2 >= (t2i - T_E - T_L - T_P - (T_A-1))) & 
                              (y2 <= (t2i - T_E - T_L - T_P)))
    
    # Find age index for adult 1 under survival distribution
    data$Numerator = mapply(data$y1, FUN=function(input) {return(which(AdultAge==(t1 - input - T_E - T_L - T_P)))})
    # Calculate numerator contribution for each pair:
    # - adult age probability × (1/2) × mother survival between births × 
    #   egg production × embryo × larva × pupa survival × adult survival to t2
    data$NumeratorTemp <- AdultAgeProbability[data$Numerator] * (1/2) * (1 - mu_A)^(abs(data$y2 - data$y1)) * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) * ((1 - mu_P)^T_P) * ((1 - mu_A)^(data$t2i - data$y2 - T_E - T_L - T_P))
    
    # Repeat patch x1 and x2 combinations for each data row
    x1 = 1:numSampledPatches
    x1 = rep(x1, nrow(data))
    x2 = 1:numSampledPatches
    x2 = rep(x2, nrow(data))
    
    data = data[rep(seq_len(nrow(data)), each = numSampledPatches), ]
    data = cbind(data, x1)
    data_x2 = cbind(data, x2)
    
    # Map sampled patch indices
    data$samp_x1 <- sampledPatches[data$x1]
    data$samp_col <- (data$y1-data$y2)*(data$y2<data$y1)+ (1 + t1 - (data$y1 + T_E + T_L + T_P))
    data_x2$samp_x2 <- sampledPatches[data_x2$x2]
    data_x2$samp_col_x2 <- (data_x2$y2-data_x2$y1)*(data_x2$y2>data_x2$y1)+ (1 + data_x2$t2i - (data_x2$y2 + T_E + T_L + T_P))
    
    # First calculate the probability that, given an adult sampled at location
    # x1, an adult sampled at location x2 is its sibling (rho_x1_x2):
    for (u in unique(data$i)) {
      # Slice data per t2 index
      slice_x1 = data %>% filter(i==u)
      slice_x2 = data_x2 %>% filter(i==u)
      
      # Build matrices holding spatial transition probabilities
      ma_x1 <- matrix(0, nrow=nrow(slice_x1), ncol=numPatches)
      for (k in 1:nrow(slice_x1)) {
        x1 = slice_x1$x1[k]
        i = slice_x1$i[k]
        temp_x1 = rho[, slice_x1$samp_x1[k], slice_x1$samp_col[k]]
        temp_x1 = temp_x1/sum(temp_x1)
        ma_x1[k, ] = temp_x1 # normalize across destinations
      }
      ma_x2 <- matrix(0, nrow=nrow(slice_x2), ncol=numPatches)
      for (j in 1:nrow(slice_x2)) {
        temp_x2 =  rho[, slice_x2$samp_x2[j], slice_x2$samp_col_x2[j]]*slice_x2$NumeratorTemp[j]
        ma_x2[j, ] = temp_x2
      }
      
      rho_x1_x2 = array(0, dim=c(numSampledPatches, ncol=numSampledPatches))
      rep_times = sum(slice_x1$x1==1) # number of replicate blocks
       
      for (rep in 1:rep_times) {
        partial_x1 = ma_x1[(1+numSampledPatches*(rep-1)):(numSampledPatches*rep), ]
        partial_x2 = ma_x2[(1+numSampledPatches*(rep-1)):(numSampledPatches*rep), ]
        
        # Sum up: probability mass across all patch combinations
        rho_x1_x2 =  rho_x1_x2 + as.numeric(partial_x1 %*% t(partial_x2))
      }
      
      Numerator[ , , u] <- rho_x1_x2
      print(u/length(unique(data$i)))
    }
    
    # Full sibling adult-adult probability:
    # Given an adult sampled at location x2 and time t2, the probability that an adult 
    # sampled at location x1 and time t1 is its full sibling is given by:
    
    P_FSAA <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSAA <- seq((-2*(T_A-1)), (2*(T_A-1)), by=1)
  }
  
  #####################################################################################
  ## FULL-SIBLING (ADULT-LARVA) KINSHIP PROBABILITIES:                               ##
  #####################################################################################
  
  # Given a larva sampled at location x2 and time t2, this is the probability that an
  # adult sampled at location x1 and time t1 is its full sibling.
  
  if ((numSampledLarvae > 0) && (numSampledAdults > 0)) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving larvae at location x2 and time t2 from
    # adult females at any consistent time (assuming a constant population size, this
    # is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    # Loop over possible birth times y2 (when mothers laid eggs producing larvae at t2)
    for (y2 in (t2 - T_E - (T_L-1)):(t2 - T_E)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of larvae at location x2 and day t2 that are full
    # siblings of an adult sampled at location x1 and day t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is -(T_L + T_P + 2*T_A), if the mother laid egg 1 at the end
    #   of her life, egg 2 at the beginning of her life, adult 1 was caught at the end of
    #   its life, & larva 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_A - T_P), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, adult 1 was caught at the 
    #   beginning of its life, & larva 2 was caught at the end of its life.
    # * So we will explore -(T_L + T_P + 2*T_A) <= t2 <= (T_A - T_P)
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- T_L - T_P - 2*T_A + 2), (T_A - T_P - 2), by=1)
    
    # Initialize a 3D array: [x1, x2, time_index]
    Numerator <- replicate(length(t2), matrix(rep(0, numSampledPatches^2), 
                                              nrow = numSampledPatches, byrow = TRUE))
    # Triple nested loop over:
    # i: t2 values
    # y1: adult birth time
    # y2: larva birth time
    print("Calculating full-sibling adult-larva probabilities:")
    col_count = (length(t2)*length((t1 - T_E - (T_L-1)):(t1 - T_E))*length((y1 - (T_A-1)):(y1 + (T_A-1))))
    vec_ma = matrix(numeric(), ncol=4, nrow=col_count)
    count = 1
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P - (T_A-1)):(t1 - T_E - T_L - T_P)) {
        for (y2 in (y1 - (T_A-1)):(y1 + (T_A-1))) {
          if ((y2 >= (t2[i] - T_E - (T_L-1))) && (y2 <= (t2[i] - T_E))) {
            NumeratorTemp <- (AdultAgeProbability[which(AdultAge==(t1 - y1 - T_E - T_L - T_P))] 
                              * (1/2) * (1 - mu_A)^(abs(y2 - y1)) 
                              * beta * ((1 - mu_E)^T_E) 
                              * ((1 - mu_L)^(t2[i] - y2 - T_E)))
            
            if (y1 >= y2) {
              third <- (y1 - y2) + (1 + t1 - (y1 + T_E + T_L + T_P))
              ma_2d <- rho[,,third]
              for (x1 in 1:numSampledPatches) {
                vec <- ma_2d[, sampledPatches[x1]]
                sum_vec <- sum(vec)
                
                for (x2 in 1:numSampledPatches) {
                  # First calculate the probability that, given a larva sampled at location
                  # x1, a larva sampled at location x2 is its sibling (rho_x1_x2):
                  rho_x1_x2 <- (vec[sampledPatches[x2]] / sum_vec)
                  
                  # Then calculate the numerator:
                  Numerator[x1, x2, i] <- Numerator[x1, x2, i] + NumeratorTemp * rho_x1_x2
                }
              }
            } else {
              vec_ma[count,] = c(i, y1, y2, NumeratorTemp)
              count <- count + 1
            }
          }
        }
      }
      print(i/length(t2))
    }
    
    data = data.frame(vec_ma[1:count,])
    colnames(data) = c('i', 'y1', 'y2', 'NumeratorTemp')
    data = data[is.na(data$i)==F,]
    data$t2i = t2[data$i]
    
    # Expand patch combinations (x1, x2) for all rows
    x1 = 1:numSampledPatches
    x1 = rep(x1, nrow(data))
    x2 = 1:numSampledPatches
    x2 = rep(x2, nrow(data))
    
    data = data[rep(seq_len(nrow(data)), each = numSampledPatches), ]
    data_x2 = cbind(data, x2)
    data = cbind(data, x1)
    
    # Add patch and time offset columns
    data$samp_x1 <- sampledPatches[data$x1]
    data$samp_col <- (1 + t1 - (data$y1 + T_E + T_L + T_P))
    data_x2$samp_x2 <- sampledPatches[data_x2$x2]
    data_x2$samp_col_x2 <- data$y2-data$y1
    
    # Calculation for each value of t2 at a time
    for (u in unique(data$i)) {
      slice_x1 = data %>% filter(i==u)
      slice_x2 = data_x2 %>% filter(i==u)
      ma_x1 <- matrix(0, nrow=nrow(slice_x1), ncol=numPatches)
      for (k in 1:nrow(slice_x1)) {
        x1 = slice_x1$x1[k]
        i = slice_x1$i[k]
        temp_x1 = rho[, slice_x1$samp_x1[k], slice_x1$samp_col[k]]
        temp_x1 = temp_x1/sum(temp_x1)
        ma_x1[k, ] = temp_x1
      }
      ma_x2 <- matrix(0, nrow=nrow(slice_x2), ncol=numPatches)
      for (j in 1:nrow(slice_x2)) {
        temp_x2 =  rho[, slice_x2$samp_x2[j], slice_x2$samp_col_x2[j]]*slice_x2$NumeratorTemp[j]
        ma_x2[j, ] = temp_x2
      }
      
      # First calculate the probability that, given a larva sampled at location
      # x1, a larva sampled at location x2 is its sibling (rho_x1_x2):
      rho_x1_x2 = array(0, dim=c(numSampledPatches, ncol=numSampledPatches))
      rep_times = sum(slice_x1$x1==1)
      
      for (rep in 1:rep_times) {
        partial_x1 = ma_x1[(1+numSampledPatches*(rep-1)):(numSampledPatches*rep), ]
        partial_x2 = ma_x2[(1+numSampledPatches*(rep-1)):(numSampledPatches*rep), ]
        rho_x1_x2 =  rho_x1_x2 + as.numeric(partial_x1 %*% t(partial_x2))
      }
      
      Numerator[ , , u] <- Numerator[ , , u] + rho_x1_x2
    }
    print('adult-larva prob')
    print(sum(Numerator))
    # Full sibling adult-larva probability:
    # Given a larva sampled at location x2 and time t2, the probability that an adult 
    # sampled at location x1 and time t1 is its full sibling is given by:
    
    P_FSAL <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSAL <- seq((- T_L - T_P - 2*T_A + 2), (T_A - T_P - 2), by=1)
  }
  
  #####################################################################################
  ## FULL-SIBLING (LARVA-ADULT) KINSHIP PROBABILITIES:                               ##
  #####################################################################################
  
  # Given an adult sampled at location x2 and time t2, this is the probability that a
  # larva sampled at location x1 and time t1 is its full sibling.
  
  if ((numSampledLarvae > 0) && (numSampledAdults > 0)) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving adults at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    # Denominator: expected number of surviving adults at time t2
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P - (T_A-1)):(t2 - T_E - T_L - T_P)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^T_P) 
                                    * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of adults at location x2 and day t2 that are full 
    # siblings of a larva sampled at location x1 and day t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (T_P - T_A), if the mother laid egg 1 at the end of her
    #   life, egg 2 at the beginning of her life, larva 1 was caught at the end of its
    #   life, & adult 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_L + T_P + 2*T_A), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, larva 1 was caught at the 
    #   beginning of its life, & adult 2 was caught at the end of its life.
    # * So we will explore (T_P - T_A) <= t2 <= (T_L + T_P + 2*T_A)
    
    # Time window to explore:
    # Earliest t2: when adult is youngest, larva is oldest → (T_P - T_A)
    # Latest t2: when adult is oldest, larva is youngest → (T_L + T_P + 2*T_A)
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((T_P - T_A + 2), (T_L + T_P + 2*T_A - 2), by=1)
    
    # Initialize numerator: stores kinship counts per patch pair and time offset
    Numerator <- replicate(length(t2), matrix(rep(0, numSampledPatches^2), 
                                              nrow = numSampledPatches, byrow = TRUE))
    
    print("Calculating full-sibling larva-adult probabilities:")
    
    # Iterate the outer loops:
    col_count = (length(t2)*length((t1 - T_E - (T_L-1)):(t1 - T_E))*length((y1 - (T_A-1)):(y1 + (T_A-1))))
    vec_ma = matrix(numeric(), ncol=4, nrow=col_count)
    count = 1
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - (T_L-1)):(t1 - T_E)){
        for (y2 in (y1 - (T_A-1)):(y1 + (T_A-1))) {
          # Check if y2 is valid relative to t2
          if ((y2 >= (t2[i] - T_E - T_L - T_P - (T_A-1))) && (y2 <= (t2[i] - T_E - T_L - T_P))) {
            # Compute survival and reproduction contribution for this (y1, y2, t2) combo
            NumeratorTemp <- (LarvaAgeProbability[which(LarvaAge==(t1 - y1 - T_E))] 
                              * (1/2) * (1 - mu_A)^(abs(y2 - y1)) 
                              * beta * ((1 - mu_E)^T_E)  * ((1 - mu_L)^T_L) * ((1 - mu_P)^T_P)
                              * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
            if (y2 >= y1) {
              # Adult born later: look up dispersal probability forward in time
              third <- (y2 - y1) + (1 + t2[i] - (y2 + T_E + T_L + T_P))
              ma_2d <- rho[,,third]
              for (x1 in 1:numSampledPatches) {
                vec <- ma_2d[sampledPatches[x1], ]
                
                
                for (x2 in 1:numSampledPatches) {
                  # First calculate the probability that, given a larva sampled at location
                  # x1, a larva sampled at location x2 is its sibling (rho_x1_x2):
                  rho_x1_x2 <- vec[sampledPatches[x2]]
                  
                  # Then calculate the numerator:
                  Numerator[x1, x2, i] <- Numerator[x1, x2, i] + NumeratorTemp * rho_x1_x2
                }
              }
            } else {
              # Larva born later: store to handle later with vectorized processing
              vec_ma[count,] = c(i, y1, y2, NumeratorTemp)
              count <- count + 1
            }
          }
        }
      }
      print(i/length(t2))
    }

    data = data.frame(vec_ma[1:count,])
    colnames(data) = c('i', 'y1', 'y2', 'NumeratorTemp')
    data = data[is.na(data$i)==F,]
    data$t2i = t2[data$i]
    x1 = 1:numSampledPatches
    x1 = rep(x1, nrow(data))
    x2 = 1:numSampledPatches
    x2 = rep(x2, nrow(data))
    
    data = data[rep(seq_len(nrow(data)), each = numSampledPatches), ]
    data_x2 = cbind(data, x2)
    data = cbind(data, x1)
    
    data$samp_x1 <- sampledPatches[data$x1]
    data$samp_col <- data$y1 - data$y2
    data_x2$samp_x2 <- sampledPatches[data_x2$x2]
    data_x2$samp_col_x2 <- (1 + data_x2$t2i - (data_x2$y2 + T_E + T_L + T_P))
    
    # For each time slice, assemble spatial matrices
    for (u in unique(data$i)) {
      slice_x1 = data %>% filter(i==u)
      slice_x2 = data_x2 %>% filter(i==u)
      ma_x1 <- matrix(0, nrow=nrow(slice_x1), ncol=numPatches)
      for (k in 1:nrow(slice_x1)) {
        x1 = slice_x1$x1[k]
        i = slice_x1$i[k]
        temp_x1 = rho[, slice_x1$samp_x1[k], slice_x1$samp_col[k]]
        temp_x1 = temp_x1/sum(temp_x1)
        ma_x1[k, ] = temp_x1
      }
      ma_x2 <- matrix(0, nrow=nrow(slice_x2), ncol=numPatches)
      for (j in 1:nrow(slice_x2)) {
        temp_x2 =  rho[, slice_x2$samp_x2[j], slice_x2$samp_col_x2[j]]*slice_x2$NumeratorTemp[j]
        ma_x2[j, ] = temp_x2
      }
      
      rho_x1_x2 = array(0, dim=c(numSampledPatches, ncol=numSampledPatches))
      # First calculate the probability that, given an adult sampled at location
      # x1, an adult sampled at location x2 is its sibling (rho_x1_x2):
      rep_times = sum(slice_x1$x1==1)
      
      for (rep in 1:rep_times) {
        partial_x1 = ma_x1[(1+numSampledPatches*(rep-1)):(numSampledPatches*rep), ]
        partial_x2 = ma_x2[(1+numSampledPatches*(rep-1)):(numSampledPatches*rep), ]
        rho_x1_x2 =  rho_x1_x2 + as.numeric(partial_x1 %*% t(partial_x2))
      }
      
      Numerator[ , , u] <- Numerator[ , , u] + rho_x1_x2
    }
    
    # Full sibling larva-adult probability:
    # Given an adult sampled at location x2 and time t2, the probability that a larva
    # sampled at location x1 and time t1 is its full sibling is given by:
    print('larva-adult prob')
    print(sum(Numerator))
    
    P_FSLA <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSLA <- seq((T_P - T_A + 2), (T_L + T_P + 2*T_A - 2), by=1)
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF FULL-SIBLING PAIR DATA:                             ##
  #####################################################################################
  
  # Record number of adults sampled on each day at each sampling patch:
  if (numSampledAdults > 0) {
    # Initialize a matrix to store the number of adults sampled per patch and day.
    # Rows = sampled patches, Columns = sampling days
    dailySampledAdults <- matrix(rep(0, numSampledPatches*numSamplingDays), 
                                 nrow = numSampledPatches, byrow = TRUE)
    
    # Loop over each patch (x1)
    for (x1 in 1:numSampledPatches) {
      # Loop over each sampling day (i)
      for (i in 1:numSamplingDays) {
        # Count how many adults were sampled at patch x1 on day i
        dailySampledAdults[x1,i] <- sum(sampledAdults[
          which(sampledAdults$Patch==sampledPatches[x1]),]$Time==samplingDays[i])
      }
    }
  }
  
  # Record number of larvae sampled on each day at each sampling patch:
  if (numSampledLarvae > 0) {
    # Initialize a matrix to store the number of larvae sampled per patch and day:
    dailySampledLarvae <- matrix(rep(0, numSampledPatches*numSamplingDays), 
                                 nrow = numSampledPatches, byrow = TRUE)
    
    # Loop over each patch
    for (x1 in 1:numSampledPatches) {
      # Loop over each sampling day
      for (i in 1:numSamplingDays) {
        # Count how many larvae were sampled at patch x1 on day i
        dailySampledLarvae[x1,i] <- sum(sampledLarvae[
          which(sampledLarvae$Patch==sampledPatches[x1]),]$Time==samplingDays[i])
      }
    }
  }
  
  # Initialize log likelihood:
  logLike <- 0
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF LARVA-ADULT FULL-SIBLING PAIR DATA:                 ##
  #####################################################################################
  
  if ((numSampledAdults > 0) && (numSampledLarvae > 0)) {
    print("Computing log likelihood of larva-adult pair data:")
    
    # Define all possible time differences between larva and adult samples
    allRelativeTimes <- seq((T_P - T_A + 2), (T_L + T_P + 2*T_A - 2), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    # Create a dataframe tracking each larva (index i)
    i = 1:numSampledLarvae
    data <- data.frame(i=i)
    # For each larva, find all adults with the same mother (full siblings)
    data$fullSibs = mapply(FUN=function(a){which(sampledAdults$momID == sampledLarvae[a, "momID"])}, data$i)
    data$larvaeTime = sampledLarvae[data$i, ]$Time
    
    # Determine the patch index for each larva
    calc_x1 <- function(input) {
      which(sampledPatches==sampledLarvae[input,]$Patch)
    }
    data$x1 = mapply(FUN=calc_x1, data$i)
    
    # Prepare combinations of larva x adult patch pairs
    x2 <- 1:numSampledPatches
    temp_x2 = data.frame(x2)
    temp_x2$patchX2 = mapply(FUN=function(a){which(sampledAdults$Patch == sampledPatches[a])}, temp_x2$x2)
    temp_x2 = temp_x2[rep(seq_len(nrow(temp_x2)), nrow(data)), ]
    data = data[rep(seq_len(nrow(data)), each = numSampledPatches), ]
    data = cbind(data, temp_x2)
    
    # For each larva-adult combination, tally full siblings by relative time
    calc_intersect <- function(a, b, c) {
      if (identical(a, integer(0))) {
        return (c())
      } else {
        union = intersect(a, b)
        fullSibRelativeTimes  = sampledAdults[union,]$Time - c
        fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
        for (k in 1:numRelativeTimes) {
          fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
        }
        return (fullSibRelativeTimesTally)
      }
    }
    
    data$fullSibRelativeTimesTally = mapply(FUN=calc_intersect, data$fullSibs, data$patchX2, data$larvaeTime)
    # Expand data across all time points
    m = 1:numRelativeTimes
    m = rep(m, nrow(data))
    
    data = data[rep(seq_len(nrow(data)), each = numRelativeTimes), ]
    data = cbind(data, m)
    data$relativeTime = allRelativeTimes[data$m]
    # Filter out cases outside sampling window
    data = data %>% filter ((relativeTime + larvaeTime >= tSamplingStart) & (relativeTime + larvaeTime <= tSamplingEnd))
    # Compute number of full sibs and probabilities
    data$numFullSibs <- mapply(FUN=calc_numFullSib, data$fullSibRelativeTimesTally, data$m)
    data$dailySampledAdultsM <- mapply(FUN=function(a,b,c) {return (dailySampledAdults[a, (b + c - tSamplingStart + 1)])}, data$x2, data$relativeTime, data$larvaeTime)
    data = data %>% filter(dailySampledAdultsM > 0)
    data$P_FSLA_value = mapply(FUN=function(a,b,c) {return(P_FSLA[a,b,c])}, data$x1, data$x2, data$m)
    data = data %>% filter(P_FSLA_value > 0)
    # Calculate log-likelihood contribution
    data$result = data$numFullSibs*log(data$P_FSLA_value) + (data$dailySampledAdultsM - data$numFullSibs)*log(1-data$P_FSLA_value)
    logLike <- logLike + sum(data$result)
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF ADULT-LARVA SIBLING PAIR DATA:                      ##
  #####################################################################################
  
  if ((numSampledAdults > 0) && (numSampledLarvae > 0)) {
    print("Computing log likelihood of adult-larva pair data:")
    
    # Define all possible time differences between adult and larva samples
    allRelativeTimes <- seq((- T_L - T_P - 2*T_A + 2), (T_A - T_P - 2), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    # Create a dataframe tracking each adult (index i)
    i = 1:numSampledAdults
    data <- data.frame(i=i)
    # For each adult, find all adults with the same mother (full siblings)
    data$fullSibs = mapply(FUN=function(a){which(sampledLarvae$momID == sampledAdults[a, "momID"])}, data$i)
    data$adultTime = sampledAdults[data$i,]$Time
    
    # Determine the patch index for each adult
    calc_x1 <- function(input) {
      which(sampledPatches==sampledAdults[input,]$Patch)
    }
    data$x1 = mapply(FUN=calc_x1, data$i)
    
    # Prepare combinations of adult x larva patch pairs
    x2 <- 1:numSampledPatches
    temp_x2 = data.frame(x2)
    temp_x2$patchX2 = mapply(FUN=function(a){which(sampledLarvae$Patch == sampledPatches[a])}, temp_x2$x2)
    temp_x2 = temp_x2[rep(seq_len(nrow(temp_x2)), nrow(data)), ]
    data = data[rep(seq_len(nrow(data)), each = numSampledPatches), ]
    data = cbind(data, temp_x2)
    
    # For each adult-larva combination, tally full siblings by relative time
    calc_intersect <- function(a, b, c) {
      if (identical(a, integer(0))) {
        return (c())
      } else {
        union = intersect(a, b)
        fullSibRelativeTimes  = sampledLarvae[union,]$Time - c
        fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
        for (k in 1:numRelativeTimes) {
          fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
        }
        return (fullSibRelativeTimesTally)
      }
    }
    
    data$fullSibRelativeTimesTally = mapply(FUN=calc_intersect, data$fullSibs, data$patchX2, data$adultTime)
    # Expand data across all time points
    m = 1:numRelativeTimes
    m = rep(m, nrow(data))
    data = data[rep(seq_len(nrow(data)), each = numRelativeTimes), ]
    data = cbind(data, m)
    data$relativeTime = allRelativeTimes[data$m]
    # Filter out cases outside sampling window
    data = data %>% filter ((relativeTime + adultTime >= tSamplingStart) & (relativeTime + adultTime <= tSamplingEnd))
    
    # Compute number of full sibs and probabilities
    data$numFullSibs <- mapply(FUN=calc_numFullSib, data$fullSibRelativeTimesTally, data$m)
    data$dailySampledLarvaeM <- mapply(FUN=function(a,b) {return (dailySampledLarvae[a + b - tSamplingStart + 1])}, data$relativeTime, data$adultTime)
    data = data %>% filter(dailySampledLarvaeM > 0)
    data$P_FSAL_value = mapply(FUN=function(a,b,c) {return(P_FSAL[a,b,c])}, data$x1, data$x2, data$m)
    data = data %>% filter(P_FSAL_value > 0)
    # Calculate log-likelihood contribution
    data$result = data$numFullSibs*log(data$P_FSAL_value) + (data$dailySampledLarvaeM - data$numFullSibs)*log(1-data$P_FSAL_value)
    logLike <- logLike + sum(data$result)
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF LARVA-LARVA SIBLING PAIR DATA:                      ##
  #####################################################################################
  
  if (numSampledLarvae > 0) {
    print("Computing log likelihood of larva-larva pair data:")
    
    # Define all possible time differences between larva and larva samples
    allRelativeTimes <- seq((-(T_A-1) - (T_L-1)), ((T_A-1) + (T_L-1)), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    # Create a dataframe tracking each larva (index i)
    i = 1:(numSampledLarvae-1)
    data <- data.frame(i=i)
    # For each larva, find all larva with the same mother (full siblings)
    data$fullSibs = mapply(FUN=function(a){return(which(sampledLarvae[(a+1):numSampledLarvae,]$momID == sampledLarvae[a, "momID"]))}, data$i)
    data$larvaeTime = sampledLarvae[data$i, ]$Time
    
    # Extract sampled patch for ith sampled larva (as referred to in sampledPatches
    # vector):
    data$x1 = mapply(FUN=function(input) {return (which(sampledPatches==sampledLarvae[input,]$Patch))}, data$i)
    
    # Prepare combinations of larva x larva patch pairs
    x2 <- 1:numSampledPatches
    x2 <- rep(x2, nrow(data))
    data = data[rep(seq_len(nrow(data)), each = numSampledPatches), ]
    data = cbind(data, x2)
    data$patchX2 = mapply(FUN=function(a, b){which(sampledLarvae[(a+1):numSampledLarvae,]$Patch == sampledPatches[b])}, data$i, data$x2)
    
    # Tally number of full-siblings a given sampled larva has at location x2 on each
    # day relative to its capture:
    calc_intersect <- function(a, b, c) {
      if (length(a)==0) {
        return (c())
      } else {
        union = intersect(a, b)
        fullSibRelativeTimes  = (sampledLarvae[(union + c),]$Time - sampledLarvae[c,]$Time)
        fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
        for (k in 1:numRelativeTimes) {
          fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
        }
        return (fullSibRelativeTimesTally)
      }
    }
    # Apply the function across all combinations
    data$fullSibRelativeTimesTally = mapply(FUN=calc_intersect, data$fullSibs, data$patchX2, data$i)
    # Calculate likelihood of observed full-sibling pair data:
    m = 1:numRelativeTimes
    m = rep(m, nrow(data))
    data = data[rep(seq_len(nrow(data)), each = numRelativeTimes), ]
    data = cbind(data, m)
    data$relativeTime = allRelativeTimes[data$m]
    # Filter out cases outside sampling window
    data = data %>% filter ((relativeTime + larvaeTime >= tSamplingStart) & (relativeTime + larvaeTime <= tSamplingEnd))
    
    # Compute number of full sibs and probabilities
    data$numFullSibs <- mapply(FUN=calc_numFullSib, data$fullSibRelativeTimesTally, data$m)
    data$dailySampledLarvaeM <- mapply(FUN=function(a,b,c) {return (dailySampledLarvae[a, (b + c - tSamplingStart + 1)])}, data$x2, data$relativeTime, data$larvaeTime)
    relativeTime0 <- which(allRelativeTimes==0)
    
    # Get total number of larvae sampled at x2, at adjusted times
    data$dailySampledLarvaeM = data$dailySampledLarvaeM - ((data$m==relativeTime0) & (data$x1==data$x2))
    
    data = data %>% filter(dailySampledLarvaeM > 0)
    data$P_FSLL_value = mapply(FUN=function(a,b,c) {return(P_FSLL[a,b,c])}, data$x1, data$x2, data$m)
    # Calculate log-likelihood contribution
    data$result = data$numFullSibs*log(data$P_FSLL_value) + (data$dailySampledLarvaeM - data$numFullSibs)*log(1-data$P_FSLL_value)
    logLike <- logLike + sum(data$result)
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF ADULT-ADULT SIBLING PAIR DATA:                      ##
  #####################################################################################
  
  if (numSampledAdults > 0) {
    print("Computing log likelihood of adult-adult pair data:")
    
    # Define all possible time differences between adult and adult samples
    allRelativeTimes <- seq((-2*(T_A-1)), (2*(T_A-1)), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    allRelativeTimes <- seq((-2*(T_A-1)), (2*(T_A-1)), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    # Create a dataframe tracking each adult (index i)
    i = 1:(numSampledAdults-1)
    data <- data.frame(i=i)
    # For each adult, find all adults with the same mother (full siblings)
    data$fullSibs = mapply(FUN=function(a){return( which(sampledAdults[(a+1):numSampledAdults,]$momID == sampledAdults[a, "momID"]))}, data$i)
    data$adultTime = sampledAdults[data$i, ]$Time
    # Extract sampled patch for ith sampled larva (as referred to in sampledPatches
    # vector):
    data$x1 = mapply(FUN=function(input) {
      which(sampledPatches==sampledAdults[input,]$Patch)
    }, data$i)
    # Tally number of full-siblings a given sampled larva has at location x2 on each
    # day relative to its capture:
    x2 <- 1:numSampledPatches
    x2 <- rep(x2, nrow(data))
    data = data[rep(seq_len(nrow(data)), each = numSampledPatches), ]
    data = cbind(data, x2)
    data$patchX2 = mapply(FUN=function(a, b){which(sampledAdults[(a+1):numSampledAdults,]$Patch == sampledPatches[b])}, data$i, data$x2)
    calc_intersect <- function(a, b, c, d) {
      if (length(a)==0) {
        return (c())
      } else {
        union = intersect(a, b)
        fullSibRelativeTimes  = (sampledAdults[(union + c),]$Time - d)
        fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
        for (k in 1:numRelativeTimes) {
          fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
        }
        return (fullSibRelativeTimesTally)
      }
    }
    data$fullSibRelativeTimesTally = mapply(FUN=calc_intersect, data$fullSibs, data$patchX2, data$i, data$adultTime)
    # Calculate likelihood of observed full-sibling pair data:
    m = 1:numRelativeTimes
    m = rep(m, nrow(data))
    data = data[rep(seq_len(nrow(data)), each = numRelativeTimes), ]
    data = cbind(data, m)
    data$relativeTime = allRelativeTimes[data$m]
    # Filter to time windows within the sampling period
    data = data %>% filter ((relativeTime + adultTime >= tSamplingStart) & (relativeTime + adultTime <= tSamplingEnd))
    # Count full siblings in each time bin
    data$numFullSibs <- mapply(FUN=calc_numFullSib, data$fullSibRelativeTimesTally, data$m)
    # Get total number of adults sampled at x2, adjusted for time
    data$dailySampledAdultsM <- mapply(FUN=function(a,b,c) {return (dailySampledAdults[a, (b + c - tSamplingStart + 1)])}, data$x2, data$relativeTime, data$adultTime)
    # Adjust for self-comparisons (exclude diagonal counts when time diff is zero)
    relativeTime0 <- which(allRelativeTimes==0)
    data$dailySampledAdultsM = data$dailySampledAdultsM - ((data$m==relativeTime0) & (data$x1==data$x2))
    data = data %>% filter(dailySampledAdultsM > 0)
    data$P_FSAA_value = mapply(FUN=function(a,b,c) {return(P_FSAA[a,b,c])}, data$x1, data$x2, data$m)
    data$result = data$numFullSibs*log(data$P_FSAA_value) + (data$dailySampledAdultsM - data$numFullSibs)*log(1-data$P_FSAA_value)
    logLike <- logLike + sum(data$result)
  }
  
  -logLike
}

#######################################################################################
# MAXIMUM-LIKELIHOOD PARAMETER INFERENCE:                                             #
#######################################################################################

# Log likelihood of all full-sibling pair data:

logLike_all <- function(DispersalPars) {
  
  # Extract parameters that we are varying from the DispersalPars vector:  
  lambda <- DispersalPars[1] * OptimScale
  barrier_strength <- DispersalPars[2]
  
  # Parameters that we are setting constant in the life history model (comment
  # out the ones you are trying to estimate):
  N_A <- 50 # Total adult mosquito population size (females & males)
  mu_A <- 0.09 # Daily mortality of adult mosquitoes
  N_F <- N_A/2 # Total adult female population size
  
  T_E <- 2 # Duration of the egg stage (days)
  T_L <- 5 # Duration of the larval stage (days)
  T_P <- 1 # Duration of the pupal stage (days)
  T_A <- 40 # Maximum adult lifespan considered (days)
  
  rM <- 1.175 # Daily mosquito population growth rate
  beta <- 20 # Daily egg production rate
  G <- T_E + T_L + T_P + (1/mu_A) # Mosquito generation time
  R_M <- rM^G # Population growth rate per generation
  
  # Daily mortality of density-independent juvenile life stages:
  mu_J <- 1 - ((2 * R_M * mu_A)/(beta * (1 - mu_A)))^(1/(T_E + T_L + T_P))
  
  mu_E <- mu_J # Daily mortality of egg stage
  mu_P <- mu_J # Daily mortality of pupal stage
  
  # Daily mortality of density-dependent larval life stage at equilibrium:
  mu_L <- 1 - (2*mu_A/(beta*(1-mu_A)*((1-mu_J)^(T_E+T_P))))^(1/T_L)
  
  # Dispersal parameters:
  # lambda <- 1/16.2 # Inverse of mean daily dispersal distance (in meters)

  # Generate the landscape for the case of a 2D grid:
  landscapeDim <- 19
  landscapeDimX <- landscapeDim
  landscapeDimY <- landscapeDim
  
  # Total number of patches:
  numPatches <- landscapeDimX * landscapeDimY
  
  # Patch coordinates (household units):
  patchCoords <- expand.grid(1:landscapeDimX, 1:landscapeDimY)
  
  # Patch coordinates (units of meters):
  patchCoords <- patchCoords * 16.6 # Taking into account spacing between houses
  
  # Matrix of distances between patches:
  distMat <- matrix(rep(0, numPatches^2), nrow = numPatches, byrow = TRUE)
  for (i in 1:numPatches) {
    for (j in 1:numPatches) {
      distMat[i, j] <- sqrt((patchCoords[i,1] - patchCoords[j,1])^2 + 
                              (patchCoords[i,2] - patchCoords[j,2])^2)
    }
  }
  
  logLike <- 0 # Initialize the log-likelihood
  
  if ((numSampledLarvae > 0) && (numSampledAdultFemales > 0)) { 
    logLike <- logLike + logLike_MOL(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                     beta, mu_E, mu_L, mu_P, 
                                     lambda, barrier_strength, numPatches, distMat) 
  }
  if ((numSampledLarvae > 0) && (numSampledAdultMales > 0)) { 
    logLike <- logLike + logLike_FOL(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                     beta, mu_E, mu_L, mu_P, 
                                     lambda, barrier_strength, numPatches, distMat) 
  }
  if (numSampledAdultFemales > 0) { 
    logLike <- logLike + logLike_MOA(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                     beta, mu_E, mu_L, mu_P, 
                                     lambda, barrier_strength, numPatches, distMat) 
  }
  if (numSampledAdultMales > 0) { 
    logLike <- logLike + logLike_FOA(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                     beta, mu_E, mu_L, mu_P, 
                                     lambda, barrier_strength, numPatches, distMat) 
  }
  logLike <- logLike + logLike_Sibs(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                    beta, mu_E, mu_L, mu_P, 
                                    lambda, barrier_strength, numPatches, distMat)
  

  logLike
}

# Choose initial values of lambda & delta for the optimization algorithm:
barrier_node <- 209
barrier_strength <- 0.15
lambda <- 1/15
OptimScale <- lambda / barrier_strength
DispersalPars <- c(lambda / OptimScale, barrier_strength)

# Choose lower & upper limits of lambda & delta for optimization algorithm:
barrier_strength_lower <- 0.01
lambda_lower <- 0.001
DispersalPars_lower <- c(lambda_lower / OptimScale, barrier_strength_lower)

barrier_strength_upper <- 0.99
lambda_upper <- 1.00
DispersalPars_upper <- c(lambda_upper / OptimScale, barrier_strength_upper)

# Find the values of lambda & delta that maximize the likelihood of the parent-offspring
# data:
DispersalPars_fit <- optimx(par=DispersalPars, fn=logLike_all, method = "nlminb", 
                            control=list(reltol = 1e-5),
                            upper = DispersalPars_upper, lower = DispersalPars_lower)
  
# The estimated values of the dispersal parameters given the kinship data are:
lambda_fit <- DispersalPars_fit$p1 * OptimScale
barrier_strength_fit <- DispersalPars_fit$p2

1 / lambda_fit
barrier_strength_fit
