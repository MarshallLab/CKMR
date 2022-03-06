#######################################################################################
## CLOSE-KIN MARK-RECAPTURE CODE FOR MOSQUITOES:                                     ##
##                                                                                   ##
## This code calculates the kinship probabilities and likelihood of observed close-  ##
## kin pairs considering mother-offspring, father-offspring, full-sibling and half-  ##
## sibling pairs where either member of a pair is larval, pupal or adult. The code   ##
## then uses an optimization algorithm to find the adult parameters - adult census   ##
## population size (N_A) and daily adult mortality rate (mu_A) that maximize the     ##
## likelihood given the data.                                                        ##
##                                                                                   ##
## Code written by John Marshall: john.marshall@berkeley.edu                         ##
## Date: February 21st, 2022                                                         ##
## Reference: Sharma Y, Bennett JB, Rasic G, Marshall JM (2022) Close-kin mark-      ##
## recapture methods to estimate demographic parameters of mosquitoes. bioRxiv doi:  ##
## https://www.biorxiv.org/content/10.1101/2022.02.19.481126v1                       ##
#######################################################################################

# Load required libraries:
library(optimx)

# Clear stored variables from the environment:
rm(list=ls())

# Directory where the data is:
setwd("...")

# Initialize the number of sampled individuals of each life stage:
numSampledLarvae <- 0
numSampledPupae <- 0
numSampledAdultMales <- 0
numSampledAdultFemales <- 0
numSampledAdults <- 0

# Load data corresponding to sampled individuals of each life stage:
if (file.exists("cut_L.csv")) { 
  sampledLarvae <- read.csv("cut_L.csv") 
  numSampledLarvae <- length(sampledLarvae[,1])
}

if (file.exists("cut_P.csv")) { 
  sampledPupae <- read.csv("cut_P.csv") 
  numSampledPupae <- length(sampledPupae[,1])
}

if (file.exists("cut_M.csv")) { 
  sampledAdultMales <- read.csv("cut_M.csv") 
  numSampledAdultMales <- length(sampledAdultMales[,1])
}

if (file.exists("cut_F.csv")) { 
  sampledAdultFemales <- read.csv("cut_F.csv") 
  numSampledAdultFemales <- length(sampledAdultFemales[,1])
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

# Record start & stop time of sampling of mosquitoes of any life stage:
if ((numSampledAdults > 0) && (numSampledLarvae == 0) && (numSampledPupae == 0)) {
  tSamplingStart <- min(sampledAdults$Time)
  tSamplingEnd <- max(sampledAdults$Time)
}
if ((numSampledAdults == 0) && (numSampledLarvae > 0) && (numSampledPupae == 0)) {
  tSamplingStart <- min(sampledLarvae$Time)
  tSamplingEnd <- max(sampledLarvae$Time)
}
if ((numSampledAdults == 0) && (numSampledLarvae == 0) && (numSampledPupae > 0)) {
  tSamplingStart <- min(sampledPupae$Time)
  tSamplingEnd <- max(sampledPupae$Time)
}
if ((numSampledAdults > 0) && (numSampledLarvae > 0) && (numSampledPupae == 0)) {
  tSamplingStart <- min(min(sampledAdults$Time), min(sampledLarvae$Time))
  tSamplingEnd <- max(max(sampledAdults$Time), max(sampledLarvae$Time))
}
if ((numSampledAdults > 0) && (numSampledLarvae == 0) && (numSampledPupae > 0)) {
  tSamplingStart <- min(min(sampledAdults$Time), min(sampledPupae$Time))
  tSamplingEnd <- max(max(sampledAdults$Time), max(sampledPupae$Time))
}
if ((numSampledAdults == 0) && (numSampledLarvae > 0) && (numSampledPupae > 0)) {
  tSamplingStart <- min(min(sampledLarvae$Time), min(sampledPupae$Time))
  tSamplingEnd <- max(max(sampledLarvae$Time), max(sampledPupae$Time))
}
if ((numSampledAdults > 0) && (numSampledLarvae > 0) && (numSampledPupae > 0)) {
  tSamplingStart <- min(min(sampledAdults$Time), min(sampledLarvae$Time), min(sampledPupae$Time))
  tSamplingEnd <- max(max(sampledAdults$Time), max(sampledLarvae$Time), max(sampledPupae$Time))
}

# Time period over which likelihood calculations will be performed:
samplingDays <- seq(tSamplingStart, tSamplingEnd, by=1)

# Total number of days over which sampling occurred:
numSamplingDays <- tSamplingEnd - tSamplingStart + 1

#######################################################################################
## MOTHER-OFFSPRING (LARVA) KINSHIP PROBABILITIES & LIKELIHOOD:                      ##
#######################################################################################

logLike_MOL <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                        beta, mu_E, mu_L, mu_P) {
  
  #####################################################################################
  ## MOTHER-OFFSPRING (LARVA) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################
  
  # Given a larva sampled at time t2, this is the probability that an adult female
  # sampled at time t1 is their mother.
  
  # First, calculate the denominator:
  # This is the expected number of surviving larvae at time t2 from adult females at any
  # consistent time (assuming a constant population size, this is independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - T_L):(t2 - T_E)) {
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of larvae at day t2 from an adult female sampled at 
  # time t1.
  # * By default, let t1 = 0, as the same equation will apply at all times.
  # * Earliest possible t2 is (-T_A + T_E), if the mother was caught at the end of her
  #   life & gave birth to the offspring soon after emergence.
  # * Latest possible t2 is (T_E + T_L), if the the mother gave birth at the time of 
  #   sampling & the larva was caught at the end of its life.
  # * So we will explore (-T_A + T_E) <= t2 <= (T_E + T_L)
  
  Numerator <- rep(0, (abs(-T_A + T_E) + (T_E + T_L) + 1))
  
  # Probability of adult surviving from 0 to T_A days:
  AdultSurvivalProbability <- rep(0, (T_A + 1))
  AdultAge <- rep(0, (T_A + 1))
  for (i in 1:(T_A + 1)){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and T_A days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-T_A + T_E), (T_E + T_L), by=1)
  
  for (i in 1:length(t2)) { 
    for (y2 in (t2[i] - T_E - T_L):(t2[i] - T_E)) {
      if ((y2 >= (t1 - T_A)) && (y2 <= t1)) {
        Numerator[i] <- Numerator[i] + ((1 - mu_A)^(t1 - y2) 
                                        * beta * ((1 - mu_E)^T_E)
                                        * ((1 - mu_L)^(t2[i] - y2 - T_E)))
      }
    }
  }
  
  # Mother-larval offspring probability:
  # Given a larva sampled at time t2, the probability that an adult female sampled at
  # time t1 is their mother is given by:
  
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
  j <- 0
  
  for (i in 1:numSampledLarvae) {
    if (any(sampledAdultFemales$myID == sampledLarvae[i, "momID"])) {
      j <- j + 1
      larvaID[j] <- sampledLarvae[i, "myID"] # Larva ID
      motherID[j] <- sampledLarvae[i, "momID"] # Mother ID
      larvaSamplingTime[j] <- sampledLarvae[i, "Time"] # Day larva sampled
      # Day mother sampled:
      motherSamplingTime[j] <- sampledAdultFemales$Time[sampledAdultFemales$myID==motherID[j]]
    }
  }
  
  # Array of mother-larval offspring pairs keeping track of: i) larva ID, ii) mother ID,
  # iii) day larva sampled, & iv) day mother sampled:
  MOL_Pairs_Data <- cbind(larvaID, motherID, larvaSamplingTime, motherSamplingTime)
  
  # Record number of adult females sampled on each day:
  dailySampledAdultFemales <- 0
  for (i in 1:numSamplingDays) {
    dailySampledAdultFemales[i] <- sum(sampledAdultFemales$Time == samplingDays[i])
  }
  
  # Record number of larvae sampled on each day:
  dailySampledLarvae <- 0
  for (i in 1:numSamplingDays) {
    dailySampledLarvae[i] <- sum(sampledLarvae$Time == samplingDays[i])
  }
  
  # Array of sampled adult females & larvae: i) sampling day, ii) number of sampled
  # adult females, & iii) number of sampled larvae:
  dailySamples <- cbind(samplingDays, dailySampledAdultFemales, dailySampledLarvae)
  
  #####################################################################################
  ## CALCULATE LIKELIHOOD OF MOTHER-OFFSPRING (LARVA) PAIR DATA:                     ##
  #####################################################################################
  
  # Format data into a 2D array containing the number of sampled larvae, w, on day t2 
  # (first index) that have mothers among the sampled adult females on day t2-t2Minust1
  # (second index). Also calculate the loglikelihood that, for x sampled larvae on day
  # t2, w have mothers among the y sampled adult females on day t1:
  
  MOL_Data_Array <- array(rep(NA, length(samplingDays) * length(t2Minust1)),
                          dim=c(length(samplingDays), length(t2Minust1)))
  logLike <- 0
  
  for (i in 1:length(samplingDays)) {
    larvaSamplingTimeI <- samplingDays[i]
    x <- dailySampledLarvae[i] # Number of sampled larvae on day t2
    
    for (j in 1:length(t2Minust1)) {
      motherSamplingTimeJ <- larvaSamplingTimeI - t2Minust1[j]
      
      if ((motherSamplingTimeJ >= tSamplingStart) & (motherSamplingTimeJ <= tSamplingEnd)) {
        # Number of sampled adult females on day t1:
        y <- dailySampledAdultFemales[which(samplingDays==motherSamplingTimeJ)]
        
        # Number of sampled larvae on day t2 that have mothers among the sampled adult 
        # females on day t1:
        w <- sum(MOL_Pairs_Data[,"motherSamplingTime"]
                 [MOL_Pairs_Data[,"larvaSamplingTime"]==
                     larvaSamplingTimeI]==motherSamplingTimeJ)
        MOL_Data_Array[i, j] <- w
        
        if (y > 0) {
          # Probability that a given sampled larva on day t2 has a mother among the y
          # sampled adult females on day t1:
          z <- 1 - ((1 - P_MOL[j])^y)
          
          # Now calculate the log likelihood that w sampled larvae on day t2 have mothers 
          # among the sampled adult females on day t1:
          logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
        }
      }
    }
  }
  
  -logLike
}

#######################################################################################
## MOTHER-OFFSPRING (ADULT) KINSHIP PROBABILITIES & LIKELIHOOD:                      ##
#######################################################################################

logLike_MOA <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                        beta, mu_E, mu_L, mu_P) {
  
  #####################################################################################
  ## MOTHER-OFFSPRING (ADULT) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################
  
  # Given an adult sampled at time t2, this is the probability that an adult female
  # sampled at time t1 is their mother.
  
  # First, calculate the denominator:
  # This is the expected number of surviving adults at time t2 from adult females at any
  # consistent time (assuming a constant population size, this is independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - T_L - T_P - T_A):(t2 - T_E - T_L - T_P)) {
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                  * ((1 - mu_P)^T_P) 
                                  * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of adults at day t2 from an adult female sampled at 
  # time t1.
  # * By default, let t1 = 0, as the same equation will apply at all times (it's the
  #   difference between t1 & t2 that matters).
  # * Earliest possible t2 is (-T_A + T_E + T_L + T_P), if the mother was caught at the 
  #   end of her life & gave birth to the offspring soon after emergence.
  # * Latest possible t2 is (T_E + T_L + + T_P + T_A), if the the mother gave birth at 
  #   the time of sampling & the adult offspring was caught at the end of its life.
  # * So we will explore (-T_A + T_E + T_L + T_P) <= t2 <= (T_E + T_L + T_P + T_A).
  
  Numerator <- rep(0, (abs(-T_A + T_E + T_L + T_P) + (T_E + T_L + T_P + T_A) + 1))
  
  # Probability of adult surviving from 0 to T_A days:
  AdultSurvivalProbability <- rep(0, (T_A + 1))
  AdultAge <- rep(0, (T_A + 1))
  for (i in 1:(T_A + 1)){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and T_A days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-T_A + T_E + T_L + T_P), (T_E + T_L + T_P + T_A), by=1)
  
  for (i in 1:length(t2)) { 
    for (y2 in (t2[i] - T_E - T_L - T_P - T_A):(t2[i] - T_E - T_L - T_P)) {
      if ((y2 >= (t1 - T_A)) && (y2 <= t1)) {
        Numerator[i] <- Numerator[i] + ((1 - mu_A)^(t1 - y2) 
                                        * beta * ((1 - mu_E)^T_E)
                                        * ((1 - mu_L)^T_L) * ((1 - mu_P)^T_P)
                                        * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
      }
    }
  }
  
  # Mother-adult offspring probability:
  # Given an adult sampled at time t2, the probability that an adult female sampled at
  # time t1 is their mother is given by:
  
  P_MOA <- Numerator/Denominator
  
  # Note that indices here relate to times in the vector t2 for which t1 = 0 and hence,
  # in general, the indices relate to times t2-t1, i.e.:
  
  t2Minust1 <- t2
  
  #####################################################################################
  ## FIND MOTHER-OFFSPRING (ADULT) PAIRS IN THE DATA:                                ##
  #####################################################################################
  
  adultID <- 0
  motherID <- 0
  adultSamplingTime <- 0
  motherSamplingTime <- 0
  j <- 0
  
  for (i in 1:numSampledAdults) {
    if (any(sampledAdultFemales$myID == sampledAdults[i, "momID"])) {
      j <- j + 1
      adultID[j] <- sampledAdults[i, "myID"] # Adult ID
      motherID[j] <- sampledAdults[i, "momID"] # Mother ID
      adultSamplingTime[j] <- sampledAdults[i, "Time"] # Day adult sampled
      # Day mother sampled:
      motherSamplingTime[j] <- sampledAdultFemales$Time[sampledAdultFemales$myID==motherID[j]]
    }
  }
  
  # Array of mother-adult offspring pairs keeping track of: i) adult ID, ii) mother ID,
  # iii) day adult sampled, & iv) day mother sampled:
  MOA_Pairs_Data <- cbind(adultID, motherID, adultSamplingTime, motherSamplingTime)
  
  # Record number of adult females sampled on each day:
  dailySampledAdultFemales <- 0
  for (i in 1:numSamplingDays) {
    dailySampledAdultFemales[i] <- sum(sampledAdultFemales$Time == samplingDays[i])
  }
  
  # Record number of adults sampled on each day:
  dailySampledAdults <- 0
  for (i in 1:numSamplingDays) {
    dailySampledAdults[i] <- sum(sampledAdults$Time == samplingDays[i])
  }
  
  # Array of sampled adult females & adults: i) sampling day, ii) number of
  # sampled adult females, & iii) number of sampled adults:
  dailySamples <- cbind(samplingDays, dailySampledAdultFemales, dailySampledAdults)
  
  #####################################################################################
  ## CALCULATE LIKELIHOOD OF MOTHER-OFFSPRING (ADULT) PAIR DATA:                     ##
  #####################################################################################
  
  # Format data into a 2D array containing the number of sampled adults, w, on day t2 
  # (first index) that have mothers among the sampled adult females on day t2-t2Minust1
  # (second index). Also calculate the loglikelihood that, for x sampled adults on day
  # t2, w have mothers among the y sampled adult females on day t1:
  
  MOA_Data_Array <- array(rep(NA, length(samplingDays) * length(t2Minust1)),
                          dim=c(length(samplingDays), length(t2Minust1)))
  logLike <- 0
  
  for (i in 1:length(samplingDays)) {
    adultSamplingTimeI <- samplingDays[i]
    x <- dailySampledAdults[i] # Number of sampled adults on day t2
    
    for (j in 1:length(t2Minust1)) {
      motherSamplingTimeJ <- adultSamplingTimeI - t2Minust1[j]
      
      if ((motherSamplingTimeJ >= tSamplingStart) & (motherSamplingTimeJ <= tSamplingEnd)) {
        # Number of sampled adult females on day t1:
        y <- dailySampledAdultFemales[which(samplingDays==motherSamplingTimeJ)]
        if (j == which(t2Minust1==0)) { y <- (y - 1) } 
        
        # Number of sampled adults on day t2 that have mothers among the sampled adult 
        # females on day t1:
        w <- sum(MOA_Pairs_Data[,"motherSamplingTime"]
                 [MOA_Pairs_Data[,"adultSamplingTime"]==
                     adultSamplingTimeI]==motherSamplingTimeJ)
        MOA_Data_Array[i, j] <- w
        
        if (y > 0) {
          # Probability that a given sampled adult on day t2 has a mother among the y
          # sampled adult females on day t1:
          z <- 1 - ((1 - P_MOA[j])^y)
          
          # Now calculate the log likelihood that w sampled adults on day t2 have mothers 
          # among the sampled adult females on day t1:
          logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
        }
      }
    }
  }
  
  -logLike
}

#######################################################################################
## FATHER-OFFSPRING (LARVA) KINSHIP PROBABILITIES & LIKELIHOOD:                      ##
#######################################################################################

logLike_FOL <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                        beta, mu_E, mu_L, mu_P) {
  
  #####################################################################################
  ## FATHER-OFFSPRING (LARVA) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################
  
  # Given a larva sampled at time t2, this is the probability that an adult male
  # sampled at time t1 is their father.
  
  # First, calculate the denominator:
  # This is the expected number of surviving larvae at time t2 from adult females at any
  # consistent time (assuming a constant population size, this is independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - T_L):(t2 - T_E)) {
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of larvae at day t2 from an adult male sampled at 
  # time t1.
  # * By default, let t1 = 0, as the same equation will apply at all times.
  # * Earliest possible t2 is (-T_A + T_E), if the father was caught at the end of his
  #   life, but mated at the beginning of his life, & if the mother gave birth soon after
  #   mating, & the larva was caught soon after emergence.
  # * Latest possible t2 is (T_A + T_E + T_L), if the the father was caught soon after 
  #   mating, the mother mated at the beginning of her life & laid eggs at the end of her
  #   life, & the larva was caught very soon before developing into a pupa.
  # * So we will explore (-T_A + T_E) <= t2 <= (T_A + T_E + T_L)
  
  Numerator <- rep(0, (abs(-T_A + T_E) + (T_A + T_E + T_L) + 1))
  
  # Probability of adult surviving from 0 to T_A days:
  AdultSurvivalProbability <- rep(0, (T_A + 1))
  AdultAge <- rep(0, (T_A + 1))
  for (i in 1:(T_A + 1)){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and T_A days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-T_A + T_E), (T_A + T_E + T_L), by=1)
  
  for (i in 1:length(t2)) { 
    for (tk in (t1 - T_A):t1) {
      for (y2 in tk:(tk + T_A)) {
        if ((t2[i] >= (y2 + T_E)) && (t2[i] <= (y2 + T_E + T_L))) {
          Numerator[i] <- Numerator[i] + (AdultAgeProbability[which(AdultAge==(t1 - tk))] 
                                          * ((1 - mu_A)^(y2 - tk))
                                          * beta * ((1 - mu_E)^T_E)
                                          * ((1 - mu_L)^(t2[i] - y2 - T_E)))
        }
      }
    }
  }
  
  # Father-larval offspring probability:
  # Given a larva sampled at time t2, the probability that an adult male sampled at time
  # t1 is their father is given by:
  
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
  j <- 0
  
  for (i in 1:numSampledLarvae) {
    if (any(sampledAdultMales$myID == sampledLarvae[i, "dadID"])) {
      j <- j + 1
      larvaID[j] <- sampledLarvae[i, "myID"] # Larva ID
      fatherID[j] <- sampledLarvae[i, "dadID"] # Father ID
      larvaSamplingTime[j] <- sampledLarvae[i, "Time"] # Day larva sampled
      # Day father sampled:
      fatherSamplingTime[j] <- sampledAdultMales$Time[sampledAdultMales$myID==fatherID[j]]
    }
  }
  
  # Array of father-larval offspring pairs keeping track of: i) larva ID, ii) father ID,
  # iii) day larva sampled, & iv) day father sampled:
  FOL_Pairs_Data <- cbind(larvaID, fatherID, larvaSamplingTime, fatherSamplingTime)
  
  # Record number of adult males sampled on each day:
  dailySampledAdultMales <- 0
  for (i in 1:numSamplingDays) {
    dailySampledAdultMales[i] <- sum(sampledAdultMales$Time == samplingDays[i])
  }
  
  # Record number of larvae sampled on each day:
  dailySampledLarvae <- 0
  for (i in 1:numSamplingDays) {
    dailySampledLarvae[i] <- sum(sampledLarvae$Time == samplingDays[i])
  }
  
  # Array of sampled adult males & larvae: i) sampling day, ii) number of sampled
  # adult males, & iii) number of sampled larvae:
  dailySamples <- cbind(samplingDays, dailySampledAdultMales, dailySampledLarvae)
  
  #####################################################################################
  ## CALCULATE LIKELIHOOD OF FATHER-OFFSPRING (LARVA) PAIR DATA:                     ##
  #####################################################################################
  
  # Format data into a 2D array containing the number of sampled larvae, w, on day t2 
  # (first index) that have fathers among the sampled adult males on day t2-t2Minust1
  # (second index). Also calculate the loglikelihood that, for x sampled larvae on day
  # t2, w have fathers among the y sampled adult males on day t1:
  
  FOL_Data_Array <- array(rep(NA, length(samplingDays) * length(t2Minust1)),
                          dim=c(length(samplingDays), length(t2Minust1)))
  logLike <- 0
  
  for (i in 1:length(samplingDays)) {
    larvaSamplingTimeI <- samplingDays[i]
    x <- dailySampledLarvae[i] # Number of sampled larvae on day t2
    
    for (j in 1:length(t2Minust1)) {
      fatherSamplingTimeJ <- larvaSamplingTimeI - t2Minust1[j]
      
      if ((fatherSamplingTimeJ >= tSamplingStart) & (fatherSamplingTimeJ <= tSamplingEnd)) {
        # Number of sampled adult males on day t1:
        y <- dailySampledAdultMales[which(samplingDays==fatherSamplingTimeJ)]
        
        # Number of sampled larvae on day t2 that have fathers among the sampled adult 
        # males on day t1:
        w <- sum(FOL_Pairs_Data[,"fatherSamplingTime"]
                 [FOL_Pairs_Data[,"larvaSamplingTime"]==
                     larvaSamplingTimeI]==fatherSamplingTimeJ)
        FOL_Data_Array[i, j] <- w
        
        if (y > 0) {
          # Probability that a given sampled larva on day t2 has a father among the y
          # sampled adult males on day t1:
          z <- 1 - ((1 - P_FOL[j])^y)
          
          # Now calculate the log likelihood that w sampled larvae on day t2 have fathers 
          # among the sampled adult males on day t1:
          logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
        }
      }
    }
  }
  
  -logLike
}

#######################################################################################
## FATHER-OFFSPRING (ADULT) KINSHIP PROBABILITIES & LIKELIHOOD:                      ##
#######################################################################################

logLike_FOA <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                        beta, mu_E, mu_L, mu_P) {
  
  #####################################################################################
  ## FATHER-OFFSPRING (ADULT) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################
  
  # Given an adult sampled at time t2, this is the probability that an adult male
  # sampled at time t1 is their father.
  
  # First, calculate the denominator:
  # This is the expected number of surviving adults at time t2 from adult females at any
  # consistent time (assuming a constant population size, this is independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - T_L - T_P - T_A):(t2 - T_E - T_L - T_P)) {
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                  * ((1 - mu_P)^T_P) 
                                  * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of adult offspring at day t2 from an adult male sampled at 
  # time t1.
  # * By default, let t1 = 0, as the same equation will apply at all times.
  # * Earliest possible t2 is (-T_A + T_E + T_L + T_P), if the father was caught at the
  #   end of his life, but mated at the beginning of his life, & if the mother gave birth
  #   soon after mating, & the adult offspring was caught soon after emergence.
  # * Latest possible t2 is (T_A + T_E + T_L + T_P + T_A), if the the father was caught 
  #   soon after mating, the mother mated at the beginning of her life & laid eggs at the
  #   end of her life, & the adult offspring was caught at the end of its life.
  # * So we will explore (-T_A + T_E + T_L + T_P) <= t2 <= (T_A + T_E + T_L + T_P + T_A)
  
  Numerator <- rep(0, (abs(-T_A + T_E + T_L + T_P) + (T_A + T_E + T_L + T_P + T_A) + 1))
  
  # Probability of adult surviving from 0 to T_A days:
  AdultSurvivalProbability <- rep(0, (T_A + 1))
  AdultAge <- rep(0, (T_A + 1))
  for (i in 1:(T_A + 1)){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and T_A days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-T_A + T_E + T_L + T_P), (T_A + T_E + T_L + T_P + T_A), by=1)
  
  for (i in 1:length(t2)) { 
    for (tk in (t1 - T_A):t1) {
      for (y2 in tk:(tk + T_A)) {
        if ((t2[i] >= (y2 + T_E + T_L + T_P)) && (t2[i] <= (y2 + T_E + T_L + T_P + T_A))) {
          Numerator[i] <- Numerator[i] + (AdultAgeProbability[which(AdultAge==(t1 - tk))] 
                                          * ((1 - mu_A)^(y2 - tk))
                                          * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                          * ((1 - mu_P)^T_P) 
                                          * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
        }
      }
    }
  }
  
  # Father-larval offspring probability:
  # Given a larva sampled at time t2, the probability that an adult male sampled at time
  # t1 is its father is given by:
  
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
  j <- 0
  
  for (i in 1:numSampledAdults) {
    if (any(sampledAdultMales$myID == sampledAdults[i, "dadID"])) {
      j <- j + 1
      adultID[j] <- sampledAdults[i, "myID"] # Adult offspring ID
      fatherID[j] <- sampledAdults[i, "dadID"] # Mother ID
      adultSamplingTime[j] <- sampledAdults[i, "Time"] # Day adult offspring sampled
      # Day father sampled:
      fatherSamplingTime[j] <- sampledAdultMales$Time[sampledAdultMales$myID==fatherID[j]]
    }
  }
  
  # Array of father-adult male offspring pairs keeping track of: i) adult ID, ii) father
  # ID, iii) day adult sampled, & iv) day father sampled:
  FOA_Pairs_Data <- cbind(adultID, fatherID, adultSamplingTime, fatherSamplingTime)
  
  # Record number of adult males sampled on each day:
  dailySampledAdultMales <- 0
  for (i in 1:numSamplingDays) {
    dailySampledAdultMales[i] <- sum(sampledAdultMales$Time == samplingDays[i])
  }
  
  # Record number of adults sampled on each day:
  dailySampledAdults <- 0
  for (i in 1:numSamplingDays) {
    dailySampledAdults[i] <- sum(sampledAdults$Time == samplingDays[i])
  }
  
  # Array of sampled adult males & adults: i) sampling day, ii) number of sampled
  # adult males, & iii) number of sampled adults:
  dailySamples <- cbind(samplingDays, dailySampledAdultMales, dailySampledAdults)
  
  #####################################################################################
  ## CALCULATE LIKELIHOOD OF FATHER-OFFSPRING (ADULT) PAIR DATA:                     ##
  #####################################################################################
  
  # Format data into a 2D array containing the number of sampled adults, w, on day t2 
  # (first index) that have fathers among the sampled adult males on day t2-t2Minust1
  # (second index). Also calculate the loglikelihood that, for x sampled adults on day
  # t2, w have fathers among the y sampled adult males on day t1:
  
  FOA_Data_Array <- array(rep(NA, length(samplingDays) * length(t2Minust1)),
                          dim=c(length(samplingDays), length(t2Minust1)))
  logLike <- 0
  
  for (i in 1:length(samplingDays)) {
    adultSamplingTimeI <- samplingDays[i]
    x <- dailySampledAdults[i] # Number of sampled adults on day t2
    
    for (j in 1:length(t2Minust1)) {
      fatherSamplingTimeJ <- adultSamplingTimeI - t2Minust1[j]
      
      if ((fatherSamplingTimeJ >= tSamplingStart) & (fatherSamplingTimeJ <= tSamplingEnd)) {
        # Number of sampled adult males on day t1:
        y <- dailySampledAdultMales[which(samplingDays==fatherSamplingTimeJ)]
        if (j == which(t2Minust1==0)) { y <- (y - 1) } 
        
        # Number of sampled adults on day t2 that have fathers among the sampled adult 
        # males on day t1:
        w <- sum(FOA_Pairs_Data[,"fatherSamplingTime"]
                 [FOA_Pairs_Data[,"adultSamplingTime"]==
                     adultSamplingTimeI]==fatherSamplingTimeJ)
        FOA_Data_Array[i, j] <- w
        
        if (y > 0) {
          # Probability that a given sampled adult on day t2 has a father among the y
          # sampled adult males on day t1:
          z <- 1 - ((1 - P_FOA[j])^y)
          
          # Now calculate the log likelihood that w sampled adults on day t2 have fathers 
          # among the sampled adult males on day t1:
          logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
        }
      }
    }
  }
  
  -logLike
}

#######################################################################################
## MOTHER-OFFSPRING (PUPA) KINSHIP PROBABILITIES & LIKELIHOOD:                       ##
#######################################################################################

logLike_MOP <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                        beta, mu_E, mu_L, mu_P) {
  
  #####################################################################################
  ## MOTHER-OFFSPRING (PUPA) KINSHIP PROBABILITIES:                                  ##
  #####################################################################################
  
  # Given a pupa sampled at time t2, this is the probability that an adult female
  # sampled at time t1 is their mother.
  
  # First, calculate the denominator:
  # This is the expected number of surviving pupae at time t2 from adult females at any
  # consistent time (assuming a constant population size, this is independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - T_L - T_P):(t2 - T_E - T_L)) {
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) 
                                  * ((1 - mu_P)^(t2 - y2 - T_E - T_L)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of pupae at day t2 from an adult female sampled at 
  # time t1.
  # * By default, let t1 = 0, as the same equation will apply at all times.
  # * Earliest possible t2 is (-T_A + T_E + T_P), if the mother was caught at the end 
  #   of her life & gave birth to the offspring soon after emergence.
  # * Latest possible t2 is (T_E + T_L + T_P), if the the mother gave birth at the time
  #   of sampling & the pupa was caught at the end of its life.
  # * So we will explore (-T_A + T_E + T_L) <= t2 <= (T_E + T_L + T+P)
  
  Numerator <- rep(0, (abs(-T_A + T_E + T_L) + (T_E + T_L + T_P) + 1))
  
  # Probability of adult surviving from 0 to T_A days:
  AdultSurvivalProbability <- rep(0, (T_A + 1))
  AdultAge <- rep(0, (T_A + 1))
  for (i in 1:(T_A + 1)){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and T_A days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-T_A + T_E + T_L), (T_E + T_L + T_P), by=1)
  
  for (i in 1:length(t2)) { 
    for (y2 in (t2[i] - T_E - T_L - T_P):(t2[i] - T_E - T_L)) {
      if ((y2 >= (t1 - T_A)) && (y2 <= t1)) {
        Numerator[i] <- Numerator[i] + ((1 - mu_A)^(t1 - y2) 
                                        * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                        * ((1 - mu_P)^(t2[i] - y2 - T_E - T_L)))
      }
    }
  }
  
  # Mother-pupal offspring probability:
  # Given a pupa sampled at time t2, the probability that an adult female sampled at
  # time t1 is its mother is given by:
  
  P_MOP <- Numerator/Denominator
  
  # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
  # in general, the indices relate to times t2-t1, i.e.:
  
  t2Minust1 <- t2
  
  #####################################################################################
  ## FIND MOTHER-OFFSPRING (PUPA) PAIRS IN THE DATA:                                 ##
  #####################################################################################
  
  pupaID <- 0
  motherID <- 0
  pupaSamplingTime <- 0
  motherSamplingTime <- 0
  j <- 0
  
  for (i in 1:numSampledPupae) {
    if (any(sampledAdultFemales$myID == sampledPupae[i, "momID"])) {
      j <- j + 1
      pupaID[j] <- sampledPupae[i, "myID"] # Pupa ID
      motherID[j] <- sampledPupae[i, "momID"] # Mother ID
      pupaSamplingTime[j] <- sampledPupae[i, "Time"] # Day pupa sampled
      # Day mother sampled:
      motherSamplingTime[j] <- sampledAdultFemales$Time[sampledAdultFemales$myID==motherID[j]]
    }
  }
  
  # Array of mother-pupal offspring pairs keeping track of: i) pupa ID, ii) mother ID,
  # iii) day pupa sampled, & iv) day mother sampled:
  MOP_Pairs_Data <- cbind(pupaID, motherID, pupaSamplingTime, motherSamplingTime)
  
  # Record number of adult females sampled on each day:
  dailySampledAdultFemales <- 0
  for (i in 1:numSamplingDays) {
    dailySampledAdultFemales[i] <- sum(sampledAdultFemales$Time == samplingDays[i])
  }
  
  # Record number of pupae sampled on each day:
  dailySampledPupae <- 0
  for (i in 1:numSamplingDays) {
    dailySampledPupae[i] <- sum(sampledPupae$Time == samplingDays[i])
  }
  
  # Array of sampled adult females & pupae: i) sampling day, ii) number of sampled
  # adult females, & iii) number of sampled pupae:
  dailySamples <- cbind(samplingDays, dailySampledAdultFemales, dailySampledPupae)
  
  #####################################################################################
  ## CALCULATE LIKELIHOOD OF MOTHER-OFFSPRING (PUPA) PAIR DATA:                      ##
  #####################################################################################
  
  # Format data into a 2D array containing the number of sampled pupae, w, on day t2 
  # (first index) that have mothers among the sampled adult females on day t2-t2Minust1
  # (second index). Also calculate the loglikelihood that, for x sampled pupae on day
  # t2, w have mothers among the y sampled adult females on day t1:
  
  MOP_Data_Array <- array(rep(NA, length(samplingDays) * length(t2Minust1)),
                          dim=c(length(samplingDays), length(t2Minust1)))
  logLike <- 0
  
  for (i in 1:length(samplingDays)) {
    pupaSamplingTimeI <- samplingDays[i]
    x <- dailySampledPupae[i] # Number of sampled pupae on day t2
    
    for (j in 1:length(t2Minust1)) {
      motherSamplingTimeJ <- pupaSamplingTimeI - t2Minust1[j]
      
      if ((motherSamplingTimeJ >= tSamplingStart) & (motherSamplingTimeJ <= tSamplingEnd)) {
        # Number of sampled adult females on day t1:
        y <- dailySampledAdultFemales[which(samplingDays==motherSamplingTimeJ)]
        
        # Number of sampled pupae on day t2 that have mothers among the sampled adult 
        # females on day t1:
        w <- sum(MOP_Pairs_Data[,"motherSamplingTime"]
                 [MOP_Pairs_Data[,"pupaSamplingTime"]==
                     pupaSamplingTimeI]==motherSamplingTimeJ)
        MOP_Data_Array[i, j] <- w
        
        if (y > 0) {
          # Probability that a given sampled pupa on day t2 has a mother among the y
          # sampled adult females on day t1:
          z <- 1 - ((1 - P_MOP[j])^y)
          
          # Now calculate the log likelihood that w sampled pupae on day t2 have mothers 
          # among the sampled adult females on day t1:
          logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
        }
      }
    }
  }
  
  -logLike
}

#######################################################################################
## FATHER-OFFSPRING (PUPA) KINSHIP PROBABILITIES & LIKELIHOOD:                       ##
#######################################################################################

logLike_FOP <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                        beta, mu_E, mu_L, mu_P) {
  
  #####################################################################################
  ## FATHER-OFFSPRING (PUPA) KINSHIP PROBABILITIES:                                  ##
  #####################################################################################
  
  # Given a pupa sampled at time t2, this is the probability that an adult male
  # sampled at time t1 is their father.
  
  # First, calculate the denominator:
  # This is the expected number of surviving pupae at time t2 from adult females at any
  # consistent time (assuming a constant population size, this is independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - T_L - T_P):(t2 - T_E - T_L)) {
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) 
                                  * ((1 - mu_P)^(t2 - y2 - T_E - T_L)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of pupae at day t2 from an adult male sampled at 
  # time t1.
  # * By default, let t1 = 0, as the same equation will apply at all times.
  # * Earliest possible t2 is (-T_A + T_E + T_L), if the father was caught at the end of 
  #   his life, but mated at the beginning of his life, & if the mother gave birth soon
  #   after mating, & the pupa was caught soon after emergence.
  # * Latest possible t2 is (T_A + T_E + T_L + T_P), if the the father was caught soon 
  #   after mating, the mother mated at the beginning of her life & laid eggs at the end 
  #   of her life, & the pupa was caught very soon before developing into an adult.
  # * So we will explore (-T_A + T_E + T_L) <= t2 <= (T_A + T_E + T_L + T_P)
  
  Numerator <- rep(0, (abs(-T_A + T_E + T_L) + (T_A + T_E + T_L + T_P) + 1))
  
  # Probability of adult surviving from 0 to T_A days:
  AdultSurvivalProbability <- rep(0, (T_A + 1))
  AdultAge <- rep(0, (T_A + 1))
  for (i in 1:(T_A + 1)){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and T_A days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-T_A + T_E + T_L), (T_A + T_E + T_L + T_P), by=1)
  
  for (i in 1:length(t2)) { 
    for (tk in (t1 - T_A):t1) {
      for (y2 in tk:(tk + T_A)) {
        if ((t2[i] >= (y2 + T_E + T_L)) && (t2[i] <= (y2 + T_E + T_L + T_P))) {
          Numerator[i] <- Numerator[i] + (AdultAgeProbability[which(AdultAge==(t1 - tk))] 
                                          * ((1 - mu_A)^(y2 - tk))
                                          * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                          * ((1 - mu_P)^(t2[i] - y2 - T_E - T_L)))
        }
      }
    }
  }
  
  # Father-pupal offspring probability:
  # Given a pupa sampled at time t2, the probability that an adult male sampled at time
  # t1 is its father is given by:
  
  P_FOP <- Numerator/Denominator
  
  # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
  # in general, the indices relate to times t2-t1, i.e.:
  
  t2Minust1 <- t2
  
  #####################################################################################
  ## FIND FATHER-OFFSPRING (PUPA) PAIRS IN THE DATA:                                 ##
  #####################################################################################
  
  pupaID <- 0
  fatherID <- 0
  pupaSamplingTime <- 0
  fatherSamplingTime <- 0
  j <- 0
  
  for (i in 1:numSampledPupae) {
    if (any(sampledAdultMales$myID == sampledPupae[i, "dadID"])) {
      j <- j + 1
      pupaID[j] <- sampledPupae[i, "myID"] # Pupa ID
      fatherID[j] <- sampledPupae[i, "dadID"] # Father ID
      pupaSamplingTime[j] <- sampledPupae[i, "Time"] # Day pupa sampled
      # Day father sampled:
      fatherSamplingTime[j] <- sampledAdultMales$Time[sampledAdultMales$myID==fatherID[j]]
    }
  }
  
  # Array of father-pupal offspring pairs keeping track of: i) pupa ID, ii) father ID,
  # iii) day pupa sampled, & iv) day father sampled:
  FOP_Pairs_Data <- cbind(pupaID, fatherID, pupaSamplingTime, fatherSamplingTime)
  
  # Record number of adult males sampled on each day:
  dailySampledAdultMales <- 0
  for (i in 1:numSamplingDays) {
    dailySampledAdultMales[i] <- sum(sampledAdultMales$Time == samplingDays[i])
  }
  
  # Record number of pupae sampled on each day:
  dailySampledPupae <- 0
  for (i in 1:numSamplingDays) {
    dailySampledPupae[i] <- sum(sampledPupae$Time == samplingDays[i])
  }
  
  # Array of sampled adult males & pupae: i) sampling day, ii) number of sampled
  # adult males, & iii) number of sampled pupae:
  dailySamples <- cbind(samplingDays, dailySampledAdultMales, dailySampledPupae)
  
  #####################################################################################
  ## CALCULATE LIKELIHOOD OF FATHER-OFFSPRING (PUPA) PAIR DATA:                      ##
  #####################################################################################
  
  # Format data into a 2D array containing the number of sampled pupae, w, on day t2 
  # (first index) that have fathers among the sampled adult males on day t2-t2Minust1
  # (second index). Also calculate the loglikelihood that, for x sampled pupae on day
  # t2, w have fathers among the y sampled adult males on day t1:
  
  FOP_Data_Array <- array(rep(NA, length(samplingDays) * length(t2Minust1)),
                          dim=c(length(samplingDays), length(t2Minust1)))
  logLike <- 0
  
  for (i in 1:length(samplingDays)) {
    pupaSamplingTimeI <- samplingDays[i]
    x <- dailySampledPupae[i] # Number of sampled pupae on day t2
    
    for (j in 1:length(t2Minust1)) {
      fatherSamplingTimeJ <- pupaSamplingTimeI - t2Minust1[j]
      
      if ((fatherSamplingTimeJ >= tSamplingStart) & (fatherSamplingTimeJ <= tSamplingEnd)) {
        # Number of sampled adult males on day t1:
        y <- dailySampledAdultMales[which(samplingDays==fatherSamplingTimeJ)]
        
        # Number of sampled pupae on day t2 that have fathers among the sampled adult 
        # males on day t1:
        w <- sum(FOP_Pairs_Data[,"fatherSamplingTime"]
                 [FOP_Pairs_Data[,"pupaSamplingTime"]==
                     pupaSamplingTimeI]==fatherSamplingTimeJ)
        FOP_Data_Array[i, j] <- w
        
        if (y > 0) {
          # Probability that a given sampled pupa on day t2 has a father among the y
          # sampled adult males on day t1:
          z <- 1 - ((1 - P_FOP[j])^y)
          
          # Now calculate the log likelihood that w sampled pupae on day t2 have fathers 
          # among the sampled adult males on day t1:
          logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
        }
      }
    }
  }
  
  -logLike
}

#######################################################################################
## FULL-SIBLING & HALF-SIBLING KINSHIP PROBABILITIES & LIKELIHOOD:                   ##
#######################################################################################

logLike_Sibs <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                         beta, mu_E, mu_L, mu_P) {
  
  #####################################################################################
  ## FULL-SIBLING (LARVA-LARVA) KINSHIP PROBABILITIES:                               ##
  #####################################################################################

  # Given a larva sampled at time t2, this is the probability that a larva sampled at
  # time t1 is its full sibling.

  # Probability of larva surviving from 0 to T_L days:
  LarvaSurvivalProbability <- rep(0, (T_L + 1))
  LarvaAge <- rep(0, (T_L + 1))
  for (i in 1:(T_L + 1)){
    LarvaSurvivalProbability[i] <- (1 - mu_L)^(i-1)
    LarvaAge[i] <- i - 1
  }
  
  # Probability of larva having an age between 0 and T_L days:
  LarvaAgeProbability <- LarvaSurvivalProbability / sum(LarvaSurvivalProbability)
  
  # Probability of adult surviving from 0 to T_A days:
  AdultSurvivalProbability <- rep(0, (T_A + 1))
  AdultAge <- rep(0, (T_A + 1))
  for (i in 1:(T_A + 1)){
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and T_A days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  if (numSampledLarvae > 0) {
  
    # First, calculate the denominator:
    # This is the expected number of surviving larvae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L):(t2 - T_E)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
    }
  
    # Next, calculate the numerator:
    # This is the expected number of larvae at day t2 that are full siblings of a larva 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (-T_A - T_L), if the mother laid egg 1 at the end of her
    #   life, egg 2 at the beginning of her life, larva 1 was caught at the end of its
    #   life, & larva 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_A + T_L), if the mother laid egg 1 at the beginning of her
    #   life, egg 2 at the end of her life, larva 1 was caught at the beginning of its
    #   life, & larva 2 was caught at the end of its life.
    # * So we will explore (-T_A - T_L) <= t2 <= (T_A + T_L)
  
    Numerator <- rep(0, (abs(-T_A - T_L) + (T_A + T_L) + 1))
  
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((-T_A - T_L), (T_A + T_L), by=1)
  
    print("Calculating full-sibling larva-larva probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L):(t1 - T_E)) {
        for (ym in (y1 - T_A):y1) {
          for (y2 in ym:(ym + T_A)) {
            if ((y2 >= (t2[i] - T_E - T_L)) && (y2 <= (t2[i] - T_E))) {
              Numerator[i] <- Numerator[i] + (LarvaAgeProbability[which(LarvaAge==(t1 - y1 - T_E))] 
                                              * AdultAgeProbability[which(AdultAge==(y1 - ym))]
                                              * (1 - mu_A)^(y2 - ym) 
                                              * beta * ((1 - mu_E)^T_E)
                                              * ((1 - mu_L)^(t2[i] - y2 - T_E)))
            } 
          }
        }
      }
      print(i/length(t2))
    }
  
    # Full sibling larva-larva probability:
    # Given a larva sampled at time t2, the probability that a larva sampled at time t1 is
    # its full sibling is given by:
    
    P_FSLL <- c(rep(0, T_A), Numerator/Denominator, rep(0, T_A))
    
    # Note that additional zeroes have been added to the beginning & end of the P_FSLL
    # vector. This is for alignment with the half-sibling relationships at the stage of
    # likelihood calculation, since half-siblings can be detected over a wider time
    # window.
    # Also note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSLL <- seq((-2*T_A - T_L), (2*T_A + T_L), by=1)
  }
  
  #####################################################################################
  ## FULL-SIBLING (ADULT-ADULT) KINSHIP PROBABILITIES:                               ##
  #####################################################################################
  
  # Given an adult sampled at time t2, this is the probability that an adult sampled at
  # time t1 is its full sibling.
  
  if (numSampledAdults > 0) {
  
    # First, calculate the denominator:
    # This is the expected number of surviving adults at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P - T_A):(t2 - T_E - T_L - T_P)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^T_P) 
                                    * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of adults at day t2 that are full siblings of an adult 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (-T_A - T_A), if the mother laid egg 1 at the end of her
    #   life, egg 2 at the beginning of her life, adult 1 was caught at the end of its
    #   life, & adult 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_A + T_A), if the mother laid egg 1 at the beginning of her
    #   life, egg 2 at the end of her life, adult 1 was caught at the beginning of its
    #   life, & adult 2 was caught at the end of its life.
    # * So we will explore (-T_A - T_A) <= t2 <= (T_A + T_A)
    
    Numerator <- rep(0, (abs(-T_A - T_A) + (T_A + T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((-T_A - T_A), (T_A + T_A), by=1)
    
    print("Calculating full-sibling adult-adult probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P - T_A):(t1 - T_E - T_L - T_P)) {
        for (ym in (y1 - T_A):y1) {
          for (y2 in ym:(ym + T_A)) {
            if ((y2 >= (t2[i] - T_E - T_L - T_P - T_A)) && (y2 <= (t2[i] - T_E - T_L - T_P))) {
              Numerator[i] <- Numerator[i] + (AdultAgeProbability[which(AdultAge==(t1 - y1 - T_E - T_L - T_P))] 
                                              * AdultAgeProbability[which(AdultAge==(y1 - ym))]
                                              * (1 - mu_A)^(y2 - ym) 
                                              * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) * ((1 - mu_P)^T_P)
                                              * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Full sibling adult-adult probability:
    # Given an adult sampled at time t2, the probability that an adult sampled at time t1 is
    # its full sibling is given by:
    
    P_FSAA <- c(rep(0, T_A), Numerator/Denominator, rep(0, T_A))
    
    # Note that additional zeroes have been added to the beginning & end of the P_FSAA
    # vector. This is for alignment with the half-sibling relationships at the stage of
    # likelihood calculation, since half-siblings can be detected over a wider time
    # window.
    # Also note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSAA <- seq((-3*T_A), (3*T_A), by=1)
  }
  
  #####################################################################################
  ## FULL-SIBLING (ADULT-LARVA) KINSHIP PROBABILITIES:                               ##
  #####################################################################################
  
  # Given a larva sampled at time t2, this is the probability that an adult sampled at
  # time t1 is its full sibling.
  
  if ((numSampledLarvae > 0) && (numSampledAdults > 0)) {
  
    # First, calculate the denominator:
    # This is the expected number of surviving larvae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L):(t2 - T_E)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of larvae at day t2 that are full siblings of an adult 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is -(T_L + T_P + 2*T_A), if the mother laid egg 1 at the end
    #   of her life, egg 2 at the beginning of her life, adult 1 was caught at the end of
    #   its life, & larva 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_A - T_P), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, adult 1 was caught at the 
    #   beginning of its life, & larva 2 was caught at the end of its life.
    # * So we will explore -(T_L + T_P + 2*T_A) <= t2 <= (T_A - T_P)
    
    Numerator <- rep(0, (abs(T_A - T_P) + (T_L + T_P + 2*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- T_L - T_P - 2*T_A), (T_A - T_P), by=1)
    
    print("Calculating full-sibling adult-larva probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P - T_A):(t1 - T_E - T_L - T_P)) {
        for (ym in (y1 - T_A):y1) {
          for (y2 in ym:(ym + T_A)) {
            if ((y2 >= (t2[i] - T_E - T_L)) && (y2 <= (t2[i] - T_E))) {
              Numerator[i] <- Numerator[i] + (AdultAgeProbability[which(AdultAge==(t1 - y1 - T_E - T_L - T_P))] 
                                              * AdultAgeProbability[which(AdultAge==(y1 - ym))]
                                              * (1 - mu_A)^(y2 - ym) 
                                              * beta * ((1 - mu_E)^T_E) 
                                              * ((1 - mu_L)^(t2[i] - y2 - T_E)))
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Full sibling adult-larva probability:
    # Given a larva sampled at time t2, the probability that an adult sampled at time t1 is
    # its full sibling is given by:
    
    P_FSAL <- c(rep(0, T_A), Numerator/Denominator, rep(0, T_A))
    
    # Note that additional zeroes have been added to the beginning & end of the P_FSAL
    # vector. This is for alignment with the half-sibling relationships at the stage of
    # likelihood calculation, since half-siblings can be detected over a wider time
    # window.
    # Also note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSAL <- seq((- T_L - T_P - 3*T_A), (2*T_A - T_P), by=1)
  }
  
  #####################################################################################
  ## FULL-SIBLING (LARVA-ADULT) KINSHIP PROBABILITIES:                               ##
  #####################################################################################
  
  # Given an adult sampled at time t2, this is the probability that a larva sampled at
  # time t1 is its full sibling.

  if ((numSampledLarvae > 0) && (numSampledAdults > 0)) {
  
    # First, calculate the denominator:
    # This is the expected number of surviving adults at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P - T_A):(t2 - T_E - T_L - T_P)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^T_P) 
                                    * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of adults at day t2 that are full siblings of a larva 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (T_P - T_A), if the mother laid egg 1 at the end of her
    #   life, egg 2 at the beginning of her life, larva 1 was caught at the end of its
    #   life, & adult 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_L + T_P + 2*T_A), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, larva 1 was caught at the 
    #   beginning of its life, & adult 2 was caught at the end of its life.
    # * So we will explore (T_P - T_A) <= t2 <= (T_L + T_P + 2*T_A)
    
    Numerator <- rep(0, (abs(T_P - T_A) + (T_L + T_P + 2*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((T_P - T_A), (T_L + T_P + 2*T_A), by=1)
    
    print("Calculating full-sibling larva-adult probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L):(t1 - T_E)) {
        for (ym in (y1 - T_A):y1) {
          for (y2 in ym:(ym + T_A)) {
            if ((y2 >= (t2[i] - T_E - T_L - T_P - T_A)) && (y2 <= (t2[i] - T_E - T_L - T_P))) {
              Numerator[i] <- Numerator[i] + (LarvaAgeProbability[which(LarvaAge==(t1 - y1 - T_E))] 
                                              * AdultAgeProbability[which(AdultAge==(y1 - ym))]
                                              * (1 - mu_A)^(y2 - ym) 
                                              * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) * ((1 - mu_P)^T_P)
                                              * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Full sibling larva-adult probability:
    # Given an adult sampled at time t2, the probability that a larva sampled at time t1 is
    # its full sibling is given by:
    
    P_FSLA <- c(rep(0, T_A), Numerator/Denominator, rep(0, T_A))
    
    # Note that additional zeroes have been added to the beginning & end of the P_FSLA
    # vector. This is for alignment with the half-sibling relationships at the stage of
    # likelihood calculation, since half-siblings can be detected over a wider time
    # window.
    # Also note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSLA <- seq((T_P - 2*T_A), (T_L + T_P + 3*T_A), by=1)
  }
  
  #####################################################################################
  ## HALF-SIBLING (LARVA-LARVA) KINSHIP PROBABILITIES:                               ##
  #####################################################################################

  # Given a larva sampled at time t2, this is the probability that a larva sampled at
  # time t1 is its half-sibling.

  if (numSampledLarvae > 0) {
  
    # First, calculate the denominator:
    # This is the expected number of surviving larvae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L):(t2 - T_E)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of larvae at day t2 that are half siblings of a larva 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (-2*T_A - T_L), if the father mated with mother 1 at the
    #   end of his life & mother 2 at the beginning of his life, mother 1 laid egg 1 
    #   after living a full T_A days while mother 2 laid egg 2 immediately, larva 1 was 
    #   caught at the end of its life, & larva 2 was caught at the beginning of its life.
    # * Latest possible t2 is (2*T_A + T_L), if the father mated with mother 1 at the
    #   beginning of his life & mother 2 at the end of his life, mother 1 laid egg 1 
    #   immediately while mother 2 laid egg 2 after living a full T_A days, larva 1 was 
    #   caught at the beginning of its life, & larva 2 was caught at the end of its life.
    # * So we will explore (-2*T_A - T_L) <= t2 <= (2*T_A + T_L)
    
    Numerator <- rep(0, (abs(-2*T_A - T_L) + (2*T_A + T_L) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((-2*T_A - T_L), (2*T_A + T_L), by=1)
    
    print("Calculating half-sibling larva-larva probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L):(t1 - T_E)) {
        for (y3 in (y1 - T_A):y1) {
          for (y5 in (y3 - T_A):y3) {
            for (y4 in y5:(y5 + T_A)) {
              for (y2 in y4:(y4 + T_A)) {
                if ((y2 >= (t2[i] - T_E - T_L)) && (y2 <= (t2[i] - T_E))) {
                  Numerator[i] <- Numerator[i] + (LarvaAgeProbability[which(LarvaAge==(t1 - y1 - T_E))] 
                                                  * AdultAgeProbability[which(AdultAge==(y1 - y3))]
                                                  * AdultAgeProbability[which(AdultAge==(y3 - y5))]
                                                  * (1 - mu_A)^(y4 - y5) * (1 - mu_A)^(y2 - y4) 
                                                  * mu_A * beta * ((1 - mu_E)^T_E)
                                                  * ((1 - mu_L)^(t2[i] - y2 - T_E)))
                }
              }
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Half sibling larva-larva probability:
    # Given a larva sampled at time t2, the probability that a larva sampled at time t1 is
    # its half-sibling is given by:
    
    P_HSLL <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
    # in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_HSLL <- t2
  }
  
  #####################################################################################
  ## HALF-SIBLING (ADULT-ADULT) KINSHIP PROBABILITIES:                               ##
  #####################################################################################
  
  # Given an adult sampled at time t2, this is the probability that an adult sampled at
  # time t1 is its half-sibling.
  
  if (numSampledAdults > 0) {
  
    # First, calculate the denominator:
    # This is the expected number of surviving adults at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P - T_A):(t2 - T_E - T_L - T_P)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^T_P) 
                                    * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of adults at day t2 that are half-siblings of an adult 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is -3*T_A, if the father mated with mother 1 at the end of
    #   his life & mother 2 at the beginning of his life, mother 1 laid egg 1 after
    #   living a full T_A days while mother 2 laid egg 2 immediately, adult 1 was caught
    #   at the end of its life, & adult 2 was caught at the beginning of its life.
    # * Latest possible t2 is 3*T_A, if the father mated with mother 1 at the beginning
    #   of his life & mother 2 at the end of his life, mother 1 laid egg 1 immediately
    #   while mother 2 laid egg 2 after living a full T_A days, adult 1 was caught at 
    #   the beginning of its life, & adult 2 was caught at the end of its life.
    # * So we will explore -3*T_A <= t2 <= 3*T_A
    
    Numerator <- rep(0, (abs(-3*T_A) + (3*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((-3*T_A), (3*T_A), by=1)
    
    print("Calculating half-sibling adult-adult probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P - T_A):(t1 - T_E - T_L - T_P)) {
        for (y3 in (y1 - T_A):y1) {
          for (y5 in (y3 - T_A):y3) {
            for (y4 in y5:(y5 + T_A)) {
              for (y2 in y4:(y4 + T_A)) {
                if ((y2 >= (t2[i] - T_E - T_L - T_P - T_A)) && (y2 <= (t2[i] - T_E - T_L - T_P))) {
                  Numerator[i] <- Numerator[i] + (AdultAgeProbability[which(AdultAge==(t1 - y1 - T_E - T_L - T_P))] 
                                                  * AdultAgeProbability[which(AdultAge==(y1 - y3))]
                                                  * AdultAgeProbability[which(AdultAge==(y3 - y5))]
                                                  * (1 - mu_A)^(y4 - y5) * (1 - mu_A)^(y2 - y4) 
                                                  * mu_A * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                                  * ((1 - mu_P)^T_P) * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
                }
              }
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Half sibling adult-adult probability:
    # Given an adult sampled at time t2, the probability that an adult sampled at time t1 is
    # its half sibling is given by:
    
    P_HSAA <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
    # in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_HSAA <- t2
  }
  
  #####################################################################################
  ## HALF-SIBLING (LARVA-ADULT) KINSHIP PROBABILITIES:                               ##
  #####################################################################################

  # Given an adult sampled at time t2, this is the probability that a larva sampled at
  # time t1 is its half-sibling.

  if ((numSampledLarvae > 0) && (numSampledAdults > 0)) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving adults at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P - T_A):(t2 - T_E - T_L - T_P)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^T_P) 
                                    * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of adults at day t2 that are half siblings of a larva 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (T_P - 2*T_A), if the father mated with mother 1 at the end of his
    #   life & mated with mother 2 at the beginning of his life, mother 1 laid egg 1 at the end of her
    #   life, mother 2 laid egg 2 at the beginning of her life, larva 1 was caught at the end of its
    #   life, & adult 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_L + T_P + 3*T_A), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, larva 1 was caught at the 
    #   beginning of its life, & adult 2 was caught at the end of its life.
    # * So we will explore (T_P - 2*T_A) <= t2 <= (T_L + T_P + 3*T_A)
    
    Numerator <- rep(0, (abs(T_P - 2*T_A) + (T_L + T_P + 3*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((T_P - 2*T_A), (T_L + T_P + 3*T_A), by=1)
    
    print("Calculating half-sibling larva-adult probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L):(t1 - T_E)) {
        for (y3 in (y1 - T_A):y1) {
          for (y5 in (y3 - T_A):y3) {
            for (y4 in y5:(y5 + T_A)) {
              for (y2 in y4:(y4 + T_A)) {
                if ((y2 >= (t2[i] - T_E - T_L - T_P - T_A)) && (y2 <= (t2[i] - T_E - T_L - T_P))) {
                  Numerator[i] <- Numerator[i] + (LarvaAgeProbability[which(LarvaAge==(t1 - y1 - T_E))] 
                                                  * AdultAgeProbability[which(AdultAge==(y1 - y3))]
                                                  * AdultAgeProbability[which(AdultAge==(y3 - y5))]
                                                  * (1 - mu_A)^(y4 - y5) * (1 - mu_A)^(y2 - y4) 
                                                  * mu_A * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                                  * ((1 - mu_P)^T_P) * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
                }
              }
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Half sibling larva-adult probability:
    # Given an adult sampled at time t2, the probability that a larva sampled at time t1 is
    # its half sibling is given by:
    
    P_HSLA <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
    # in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_HSLA <- t2
  }
  
  #####################################################################################
  ## HALF-SIBLING (ADULT-LARVA) KINSHIP PROBABILITIES:                               ##
  #####################################################################################

  # Given a larva sampled at time t2, this is the probability that an adult sampled at
  # time t1 is its half sibling.

  if ((numSampledLarvae > 0) && (numSampledAdults > 0)) {
  
    # First, calculate the denominator:
    # This is the expected number of surviving larvae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L):(t2 - T_E)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of larvae at day t2 that are half siblings of an adult
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is - (T_L + T_P + 3*T_A), if the mother laid egg 1 at the end
    #   of her life, egg 2 at the beginning of her life, adult 1 was caught at the end of
    #   its life, & larva 2 was caught at the beginning of its life.
    # * Latest possible t2 is (2*T_A - T_P), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, adult 1 was caught at the 
    #   beginning of its life, & larva 2 was caught at the end of its life.
    # * So we will explore -(T_L + T_P + 3*T_A) <= t2 <= (2*T_A - T_P)
    
    Numerator <- rep(0, (abs(2*T_A - T_P) + (T_L + T_P + 3*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- T_L - T_P - 3*T_A), (2*T_A - T_P), by=1)
    
    print("Calculating half-sibling adult-larva probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P - T_A):(t1 - T_E - T_L - T_P)) {
        for (y3 in (y1 - T_A):y1) {
          for (y5 in (y3 - T_A):y3) {
            for (y4 in y5:(y5 + T_A)) {
              for (y2 in y4:(y4 + T_A)) {
                if ((y2 >= (t2[i] - T_E - T_L)) && (y2 <= (t2[i] - T_E))) {
                  Numerator[i] <- Numerator[i] + (AdultAgeProbability[which(AdultAge==(t1 - y1 - T_E - T_L - T_P))] 
                                                  * AdultAgeProbability[which(AdultAge==(y1 - y3))]
                                                  * AdultAgeProbability[which(AdultAge==(y3 - y5))]
                                                  * (1 - mu_A)^(y4 - y5) * (1 - mu_A)^(y2 - y4) 
                                                  * mu_A * beta * ((1 - mu_E)^T_E) 
                                                  * ((1 - mu_L)^(t2[i] - y2 - T_E)))
                }
              }
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Half sibling adult-larva probability:
    # Given a larva sampled at time t2, the probability that an adult sampled at time t1 is
    # its half sibling is given by:
    
    P_HSAL <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
    # in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_HSAL <- t2
  }
  
  #####################################################################################
  ## FULL-SIBLING (PUPA-PUPA) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################

  # Given a pupa sampled at time t2, this is the probability that a pupa sampled at
  # time t1 is its full-sibling.

  if (numSampledPupae > 0) {
  
    # First, calculate the denominator:
    # This is the expected number of surviving pupae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P):(t2 - T_E - T_L)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) * ((1 - mu_P)^(t2 - y2 - T_E - T_L)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of pupae at day t2 that are full siblings of a pupa 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (-T_A - T_L), if the mother laid egg 1 at the end of her
    #   life, egg 2 at the beginning of her life, pupa 1 was caught at the end of its
    #   life, & pupa 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_A + T_L), if the mother laid egg 1 at the beginning of her
    #   life, egg 2 at the end of her life, pupa 1 was caught at the beginning of its
    #   life, & pupa 2 was caught at the end of its life.
    # * So we will explore (-T_A - T_L) <= t2 <= (T_A + T_L)
    
    Numerator <- rep(0, (abs(-T_A - T_P) + (T_A + T_P) + 1))
    
    # Probability of larva surviving from 0 to T_L days:
    LarvaSurvivalProbability <- rep(0, (T_L + 1))
    LarvaAge <- rep(0, (T_L + 1))
    for (i in 1:(T_L + 1)){
      LarvaSurvivalProbability[i] <- (1 - mu_L)^(i-1)
      LarvaAge[i] <- i - 1
    }
    
    # Probability of larva having an age between 0 and T_L days:
    LarvaAgeProbability <- LarvaSurvivalProbability / sum(LarvaSurvivalProbability)
    
    # Probability of pupa surviving from 0 to T_P days:
    PupaSurvivalProbability <- rep(0, (T_P + 1))
    PupaAge <- rep(0, (T_P + 1))
    for (i in 1:(T_P + 1)){
      PupaSurvivalProbability[i] <- (1 - mu_P)^(i-1)
      PupaAge[i] <- i - 1
    }
    
    # Probability of pupa having an age between 0 and T_P days:
    PupaAgeProbability <- PupaSurvivalProbability / sum(PupaSurvivalProbability)
    
    # Probability of adult surviving from 0 to T_A days:
    AdultSurvivalProbability <- rep(0, (T_A + 1))
    AdultAge <- rep(0, (T_A + 1))
    for (i in 1:(T_A + 1)){
      AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
      AdultAge[i] <- i - 1
    }
    
    # Probability of adult having an age between 0 and T_A days:
    AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((-T_A - T_P), (T_A + T_P), by=1)
    
    print("Calculating full-sibling pupa-pupa probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P):(t1 - T_E - T_L)) {
        for (ym in (y1 - T_A):y1) {
          for (y2 in ym:(ym + T_A)) {
            if ((y2 >= (t2[i] - T_E - T_L - T_P)) && (y2 <= (t2[i] - T_E - T_L))) {
              Numerator[i] <- Numerator[i] + (PupaAgeProbability[which(PupaAge==(t1 - y1 - T_E - T_L))] 
                                              * AdultAgeProbability[which(AdultAge==(y1 - ym))]
                                              * (1 - mu_A)^(y2 - ym) 
                                              * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                              * ((1 - mu_P)^(t2[i] - y2 - T_E - T_L)))
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Full sibling pupa-pupa probability:
    # Given a pupa sampled at time t2, the probability that a pupa sampled at time t1 is
    # its full sibling is given by:
    
    P_FSPP <- c(rep(0, T_A), Numerator/Denominator, rep(0, T_A))
    
    # Note that additional zeroes have been added to the beginning & end of the P_FSPP
    # vector. This is for alignment with the half-sibling relationships at the stage of
    # likelihood calculation, since half-siblings can be detected over a wider time
    # window.
    # Also note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSPP <- seq((-2*T_A - T_P), (2*T_A + T_P), by=1)
  }
      
  #####################################################################################
  ## HALF-SIBLING (PUPA-PUPA) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################

  # Given a pupa sampled at time t2, this is the probability that a pupa sampled at
  # time t1 is its half sibling.

  if (numSampledPupae > 0) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving pupae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P):(t2 - T_E - T_L)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) * ((1 - mu_P)^(t2 - y2 - T_E - T_L)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of pupae at day t2 that are half-siblings of a pupa 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (-2*T_A - T_L), if the father mated with mother 1 at the
    #   end of his life & mother 2 at the beginning of his life, mother 1 laid egg 1 
    #   after living a full T_A days while mother 2 laid egg 2 immediately, pupa 1 was 
    #   caught at the end of its life, & pupa 2 was caught at the beginning of its life.
    # * Latest possible t2 is (2*T_A + T_L), if the father mated with mother 1 at the
    #   beginning of his life & mother 2 at the end of his life, mother 1 laid egg 1 
    #   immediately while mother 2 laid egg 2 after living a full T_A days, pupa 1 was 
    #   caught at the beginning of its life, & pupa 2 was caught at the end of its life.
    # * So we will explore (-2*T_A - T_L) <= t2 <= (2*T_A + T_L)
    
    Numerator <- rep(0, (abs(-2*T_A - T_P) + (2*T_A + T_P) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((-2*T_A - T_P), (2*T_A + T_P), by=1)
    
    print("Calculating half-sibling pupa-pupa probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P):(t1 - T_E - T_L)) {
        for (y3 in (y1 - T_A):y1) {
          for (y5 in (y3 - T_A):y3) {
            for (y4 in y5:(y5 + T_A)) {
              for (y2 in y4:(y4 + T_A)) {
                if ((y2 >= (t2[i] - T_E - T_L - T_P)) && (y2 <= (t2[i] - T_E - T_L))) {
                  Numerator[i] <- Numerator[i] + (PupaAgeProbability[which(PupaAge==(t1 - y1 - T_E - T_L))] 
                                                  * AdultAgeProbability[which(AdultAge==(y1 - y3))]
                                                  * AdultAgeProbability[which(AdultAge==(y3 - y5))]
                                                  * (1 - mu_A)^(y4 - y5) * (1 - mu_A)^(y2 - y4) 
                                                  * mu_A * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                                  * ((1 - mu_P)^(t2[i] - y2 - T_E - T_L)))
                }
              }
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Half sibling pupa-pupa probability:
    # Given a pupa sampled at time t2, the probability that a pupa sampled at time t1 is
    # its half sibling is given by:
    
    P_HSPP <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
    # in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_HSPP <- t2
  }
  
  #####################################################################################
  ## FULL-SIBLING (PUPA-ADULT) KINSHIP PROBABILITIES:                                ##
  #####################################################################################
  
  # Given an adult sampled at time t2, this is the probability that a pupa sampled at
  # time t1 is its full-sibling.

  if ((numSampledPupae > 0) && (numSampledAdults > 0)) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving adults at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P - T_A):(t2 - T_E - T_L - T_P)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^T_P) 
                                    * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of adults at day t2 that are full siblings of a pupa 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (-T_A), if the mother laid egg 1 at the end of her
    #   life, egg 2 at the beginning of her life, pupa 1 was caught at the end of its
    #   life, & pupa 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_P + 2*T_A), if the mother laid egg 1 at the beginning of her
    #   life, egg 2 at the end of her life, pupa 1 was caught at the beginning of its
    #   life, & pupa 2 was caught at the end of its life.
    # * So we will explore (-T_A) <= t2 <= (T_P + 2*T_A)
    
    Numerator <- rep(0, (abs(- T_A) + (T_P + 2*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- T_A), (T_P + 2*T_A), by=1)
    
    print("Calculating full-sibling pupa-adult probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P):(t1 - T_E - T_L)) {
        for (ym in (y1 - T_A):y1) {
          for (y2 in ym:(ym + T_A)) {
            if ((y2 >= (t2[i] - T_E - T_L - T_P - T_A)) && (y2 <= (t2[i] - T_E - T_L - T_P))) {
              Numerator[i] <- Numerator[i] + (PupaAgeProbability[which(PupaAge==(t1 - y1 - T_E - T_L))] 
                                              * AdultAgeProbability[which(AdultAge==(y1 - ym))]
                                              * (1 - mu_A)^(y2 - ym) 
                                              * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) * ((1 - mu_P)^T_P)
                                              * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Full sibling pupa-adult probability:
    # Given an adult sampled at time t2, the probability that a pupa sampled at time t1 is
    # its full sibling is given by:
    
    P_FSPA <- c(rep(0, T_A), Numerator/Denominator, rep(0, T_A))
    
    # Note that additional zeroes have been added to the beginning & end of the P_FSPA
    # vector. This is for alignment with the half-sibling relationships at the stage of
    # likelihood calculation, since half-siblings can be detected over a wider time
    # window.
    # Also note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSPA <- seq((- 2*T_A), (T_P + 3*T_A), by=1)
  }
  
  #####################################################################################
  ## FULL-SIBLING (ADULT-PUPA) KINSHIP PROBABILITIES:                                ##
  #####################################################################################

  # Given a pupa sampled at time t2, this is the probability that an adult sampled at  #
  # time t1 is its full sibling.                                                        #

  if ((numSampledPupae > 0) && (numSampledAdults > 0)) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving pupae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P):(t2 - T_E - T_L)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^(t2 - y2 - T_E - T_L)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of pupae at day t2 that are full siblings of an adult 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (- T_P - 2*T_A), if the mother laid egg 1 at the end
    #   of her life, egg 2 at the beginning of her life, adult 1 was caught at the end of
    #   its life, & pupa 2 was caught at the beginning of its life.
    # * Latest possible t2 is T_A, if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, adult 1 was caught at the 
    #   beginning of its life, & pupa 2 was caught at the end of its life.
    # * So we will explore (- T_P - 2*T_A) <= t2 <= T_A
    
    Numerator <- rep(0, ((T_A) + (T_P + 2*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- T_P - 2*T_A), (T_A), by=1)
    
    print("Calculating full-sibling adult-pupa probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P - T_A):(t1 - T_E - T_L - T_P)) {
        for (ym in (y1 - T_A):y1) {
          for (y2 in ym:(ym + T_A)) {
            if ((y2 >= (t2[i] - T_E - T_L - T_P)) && (y2 <= (t2[i] - T_E - T_L))) {
              Numerator[i] <- Numerator[i] + (AdultAgeProbability[which(AdultAge==(t1 - y1 - T_E - T_L - T_P))] 
                                              * AdultAgeProbability[which(AdultAge==(y1 - ym))]
                                              * (1 - mu_A)^(y2 - ym) 
                                              * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                              * ((1 - mu_P)^(t2[i] - y2 - T_E - T_L)))
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Full sibling adult-pupa probability:
    # Given a pupa sampled at time t2, the probability that an adult sampled at time t1 is
    # its full sibling is given by:
    
    P_FSAP <- c(rep(0, T_A), Numerator/Denominator, rep(0, T_A))
    
    # Note that additional zeroes have been added to the beginning & end of the P_FSAP
    # vector. This is for alignment with the half-sibling relationships at the stage of
    # likelihood calculation, since half-siblings can be detected over a wider time
    # window.
    # Also note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSAP <- seq((- T_P - 3*T_A), (2*T_A), by=1)
  }
      
  #####################################################################################
  ## HALF-SIBLING (PUPA-ADULT) KINSHIP PROBABILITIES:                                ##
  #####################################################################################

  # Given an adult sampled at time t2, this is the probability that a pupa sampled at
  # time t1 is its half sibling.

  if ((numSampledPupae > 0) && (numSampledAdults > 0)) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving adults at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P - T_A):(t2 - T_E - T_L - T_P)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^T_P) 
                                    * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of adults at day t2 that are half siblings of a pupa 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (- 2*T_A), if the father mated with mother 1 at the end of his
    #   life & mated with mother 2 at the beginning of his life, mother 1 laid egg 1 at the end of her
    #   life, mother 2 laid egg 2 at the beginning of her life, pupa 1 was caught at the end of its
    #   life, & adult 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_P + 3*T_A), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, pupa 1 was caught at the 
    #   beginning of its life, & adult 2 was caught at the end of its life.
    # * So we will explore (- 2*T_A) <= t2 <= (T_P + 3*T_A)
    
    Numerator <- rep(0, (abs(- 2*T_A) + (T_P + 3*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- 2*T_A), (T_P + 3*T_A), by=1)
    
    print("Calculating half-sibling pupa-adult probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P):(t1 - T_E - T_L)) {
        for (y3 in (y1 - T_A):y1) {
          for (y5 in (y3 - T_A):y3) {
            for (y4 in y5:(y5 + T_A)) {
              for (y2 in y4:(y4 + T_A)) {
                if ((y2 >= (t2[i] - T_E - T_L - T_P - T_A)) && (y2 <= (t2[i] - T_E - T_L - T_P))) {
                  Numerator[i] <- Numerator[i] + (PupaAgeProbability[which(PupaAge==(t1 - y1 - T_E - T_L))] 
                                                  * AdultAgeProbability[which(AdultAge==(y1 - y3))]
                                                  * AdultAgeProbability[which(AdultAge==(y3 - y5))]
                                                  * (1 - mu_A)^(y4 - y5) * (1 - mu_A)^(y2 - y4) 
                                                  * mu_A * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                                  * ((1 - mu_P)^T_P) * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
                }
              }
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Half sibling pupa-adult probability:
    # Given an adult sampled at time t2, the probability that a pupa sampled at time t1 is
    # its half sibling is given by:
    
    P_HSPA <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
    # in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_HSPA <- t2
  }
  
  #####################################################################################
  ## HALF-SIBLING (ADULT-PUPA) KINSHIP PROBABILITIES:                                ##
  #####################################################################################
  
  # Given a pupa sampled at time t2, this is the probability that an adult sampled at
  # time t1 is its half sibling.

  if ((numSampledPupae > 0) && (numSampledAdults > 0)) {
  
    # First, calculate the denominator:
    # This is the expected number of surviving pupae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P):(t2 - T_E - T_L)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) * ((1 - mu_P)^(t2 - y2 - T_E - T_L)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of pupae at day t2 that are half siblings of an adult 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (- T_P - 3*T_A), if the mother laid egg 1 at the end
    #   of her life, egg 2 at the beginning of her life, adult 1 was caught at the end of
    #   its life, & pupa 2 was caught at the beginning of its life.
    # * Latest possible t2 is (2*T_A), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, adult 1 was caught at the 
    #   beginning of its life, & pupa 2 was caught at the end of its life.
    # * So we will explore (- T_P - 3*T_A) <= t2 <= (2*T_A)
    
    Numerator <- rep(0, (abs(2*T_A) + (T_P + 3*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- T_P - 3*T_A), (2*T_A), by=1)
    
    print("Calculating half-sibling adult-pupa probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P - T_A):(t1 - T_E - T_L - T_P)) {
        for (y3 in (y1 - T_A):y1) {
          for (y5 in (y3 - T_A):y3) {
            for (y4 in y5:(y5 + T_A)) {
              for (y2 in y4:(y4 + T_A)) {
                if ((y2 >= (t2[i] - T_E - T_L - T_P)) && (y2 <= (t2[i] - T_E - T_L))) {
                  Numerator[i] <- Numerator[i] + (AdultAgeProbability[which(AdultAge==(t1 - y1 - T_E - T_L - T_P))] 
                                                  * AdultAgeProbability[which(AdultAge==(y1 - y3))]
                                                  * AdultAgeProbability[which(AdultAge==(y3 - y5))]
                                                  * (1 - mu_A)^(y4 - y5) * (1 - mu_A)^(y2 - y4) 
                                                  * mu_A * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) 
                                                  * ((1 - mu_P)^(t2[i] - y2 - T_E - T_L)))
                }
              }
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Half sibling adult-pupa probability:
    # Given a pupa sampled at time t2, the probability that an adult sampled at time t1 is
    # its half sibling is given by:
    
    P_HSAP <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
    # in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_HSAP <- t2
  }
  
  #####################################################################################
  ## FULL-SIBLING (LARVA-PUPA) KINSHIP PROBABILITIES:                                ##
  #####################################################################################
  
  # Given a pupa sampled at time t2, this is the probability that a larva sampled at
  # time t1 is its full sibling.

  if ((numSampledLarvae > 0) && (numSampledPupae > 0)) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving pupae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P):(t2 - T_E - T_L)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^(t2 - y2 - T_E - T_L)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of pupae at day t2 that are full siblings of a larva 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (- T_A), if the mother laid egg 1 at the end of her
    #   life, egg 2 at the beginning of her life, larva 1 was caught at the end of its
    #   life, & pupa 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_L + T_P + T_A), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, larva 1 was caught at the 
    #   beginning of its life, & pupa 2 was caught at the end of its life.
    # * So we will explore (- T_A) <= t2 <= (T_L + T_P + T_A)
    
    Numerator <- rep(0, (abs(- T_A) + (T_L + T_P + T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- T_A), (T_L + T_P + T_A), by=1)
    
    print("Calculating full-sibling larva-pupa probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L):(t1 - T_E)) {
        for (ym in (y1 - T_A):y1) {
          for (y2 in ym:(ym + T_A)) {
            if ((y2 >= (t2[i] - T_E - T_L - T_P)) && (y2 <= (t2[i] - T_E - T_L))) {
              Numerator[i] <- Numerator[i] + (LarvaAgeProbability[which(LarvaAge==(t1 - y1 - T_E))] 
                                              * AdultAgeProbability[which(AdultAge==(y1 - ym))]
                                              * (1 - mu_A)^(y2 - ym) 
                                              * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) 
                                              * ((1 - mu_P)^(t2[i] - y2 - T_E - T_L)))
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Full sibling larva-pupa probability:
    # Given a pupa sampled at time t2, the probability that a larva sampled at time t1 is
    # its full sibling is given by:
    
    P_FSLP <- c(rep(0, T_A), Numerator/Denominator, rep(0, T_A))
    
    # Note that additional zeroes have been added to the beginning & end of the P_FSLP
    # vector. This is for alignment with the half-sibling relationships at the stage of
    # likelihood calculation, since half-siblings can be detected over a wider time
    # window.
    # Also note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSLP <- seq((- 2*T_A), (T_L + T_P + 2*T_A), by=1)
  }
  
  #####################################################################################
  ## FULL-SIBLING (PUPA-LARVA) KINSHIP PROBABILITIES:                                ##
  #####################################################################################
  
  # Given a larva sampled at time t2, this is the probability that a pupa sampled at
  # time t1 is its full sibling.
  
  if ((numSampledLarvae > 0) && (numSampledPupae > 0)) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving larvae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L):(t2 - T_E)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of larvae at day t2 that are full siblings of a pupa 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (- T_L - T_P - T_A), if the mother laid egg 1 at the end
    #   of her life, egg 2 at the beginning of her life, pupa 1 was caught at the end of
    #   its life, & larva 2 was caught at the beginning of its life.
    # * Latest possible t2 is T_A, if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, pupa 1 was caught at the 
    #   beginning of its life, & larva 2 was caught at the end of its life.
    # * So we will explore (- T_L - T_P - T_A) <= t2 <= T_A
    
    Numerator <- rep(0, ((T_A) + (T_L + T_P + T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- T_L - T_P - T_A), (T_A), by=1)
    
    print("Calculating full-sibling pupa-larva probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P):(t1 - T_E - T_L)) {
        for (ym in (y1 - T_A):y1) {
          for (y2 in ym:(ym + T_A)) {
            if ((y2 >= (t2[i] - T_E - T_L)) && (y2 <= (t2[i] - T_E))) {
              Numerator[i] <- Numerator[i] + (PupaAgeProbability[which(PupaAge==(t1 - y1 - T_E - T_L))] 
                                              * AdultAgeProbability[which(AdultAge==(y1 - ym))]
                                              * (1 - mu_A)^(y2 - ym) 
                                              * beta * ((1 - mu_E)^T_E) 
                                              * ((1 - mu_L)^(t2[i] - y2 - T_E)))
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Full sibling pupa-larva probability:
    # Given a larva sampled at time t2, the probability that a pupa sampled at time t1 is
    # its full sibling is given by:
    
    P_FSPL <- c(rep(0, T_A), Numerator/Denominator, rep(0, T_A))
    
    # Note that additional zeroes have been added to the beginning & end of the P_FSPL
    # vector. This is for alignment with the half-sibling relationships at the stage of
    # likelihood calculation, since half-siblings can be detected over a wider time
    # window.
    # Also note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSPL <- seq((- T_L - T_P - 2*T_A), (2*T_A), by=1)
  }
  
  #####################################################################################
  ## HALF-SIBLING (LARVA-PUPA) KINSHIP PROBABILITIES:                                ##
  #####################################################################################
  
  # Given a pupa sampled at time t2, this is the probability that a larva sampled at
  # time t1 is its half sibling.

  if ((numSampledLarvae > 0) && (numSampledPupae > 0)) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving pupae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P):(t2 - T_E - T_L)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^(t2 - y2 - T_E - T_L)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of pupae at day t2 that are half siblings of a larva 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (- 2*T_A), if the father mated with mother 1 at the end of his
    #   life & mated with mother 2 at the beginning of his life, mother 1 laid egg 1 at the end of her
    #   life, mother 2 laid egg 2 at the beginning of her life, larva 1 was caught at the end of its
    #   life, & pupa 2 was caught at the beginning of its life.
    # * Latest possible t2 is (T_L + T_P + 2*T_A), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, larva 1 was caught at the 
    #   beginning of its life, & pupa 2 was caught at the end of its life.
    # * So we will explore (- 2*T_A) <= t2 <= (T_L + T_P + 2*T_A)
    
    Numerator <- rep(0, (abs(- 2*T_A) + (T_L + T_P + 2*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- 2*T_A), (T_L + T_P + 2*T_A), by=1)
    
    print("Calculating half-sibling larva-pupa probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L):(t1 - T_E)) {
        for (y3 in (y1 - T_A):y1) {
          for (y5 in (y3 - T_A):y3) {
            for (y4 in y5:(y5 + T_A)) {
              for (y2 in y4:(y4 + T_A)) {
                if ((y2 >= (t2[i] - T_E - T_L - T_P)) && (y2 <= (t2[i] - T_E - T_L))) {
                  Numerator[i] <- Numerator[i] + (LarvaAgeProbability[which(LarvaAge==(t1 - y1 - T_E))] 
                                                  * AdultAgeProbability[which(AdultAge==(y1 - y3))]
                                                  * AdultAgeProbability[which(AdultAge==(y3 - y5))]
                                                  * (1 - mu_A)^(y4 - y5) * (1 - mu_A)^(y2 - y4) 
                                                  * mu_A * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                                  * ((1 - mu_P)^(t2[i] - y2 - T_E - T_L)))
                }
              }
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Half sibling larva-pupa probability:
    # Given a pupa sampled at time t2, the probability that a larva sampled at time t1 is
    # its half-sibling is given by:
    
    P_HSLP <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
    # in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_HSLP <- t2
  }
  
  #####################################################################################
  ## HALF-SIBLING (PUPA-LARVA) KINSHIP PROBABILITIES:                                ##
  #####################################################################################
  
  # Given a larva sampled at time t2, this is the probability that a pupa sampled at
  # time t1 is its half-sibling.

  if ((numSampledLarvae > 0) && (numSampledPupae > 0)) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving larvae at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L):(t2 - T_E)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of larvae at day t2 that are half siblings of a pupa
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (- T_L - T_P - 2*T_A), if the mother laid egg 1 at the end
    #   of her life, egg 2 at the beginning of her life, pupa 1 was caught at the end of
    #   its life, & larva 2 was caught at the beginning of its life.
    # * Latest possible t2 is (2*T_A), if the mother laid egg 1 at the
    #   beginning of her life, egg 2 at the end of her life, pupa 1 was caught at the 
    #   beginning of its life, & larva 2 was caught at the end of its life.
    # * So we will explore (- T_L - T_P - 2*T_A) <= t2 <= (2*T_A)
    
    Numerator <- rep(0, (abs(2*T_A) + (T_L + T_P + 2*T_A) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((- T_L - T_P - 2*T_A), (2*T_A), by=1)
    
    print("Calculating half-sibling pupa-larva probabilities:")
    for (i in 1:length(t2)) {
      for (y1 in (t1 - T_E - T_L - T_P):(t1 - T_E - T_L)) {
        for (y3 in (y1 - T_A):y1) {
          for (y5 in (y3 - T_A):y3) {
            for (y4 in y5:(y5 + T_A)) {
              for (y2 in y4:(y4 + T_A)) {
                if ((y2 >= (t2[i] - T_E - T_L)) && (y2 <= (t2[i] - T_E))) {
                  Numerator[i] <- Numerator[i] + (PupaAgeProbability[which(PupaAge==(t1 - y1 - T_E - T_L))] 
                                                  * AdultAgeProbability[which(AdultAge==(y1 - y3))]
                                                  * AdultAgeProbability[which(AdultAge==(y3 - y5))]
                                                  * (1 - mu_A)^(y4 - y5) * (1 - mu_A)^(y2 - y4) 
                                                  * mu_A * beta * ((1 - mu_E)^T_E) 
                                                  * ((1 - mu_L)^(t2[i] - y2 - T_E)))
                }
              }
            } 
          }
        }
      }
      print(i/length(t2))
    }
    
    # Half-sibling pupa-larva probability:
    # Given a larva sampled at time t2, the probability that a pupa sampled at time t1 is
    # its half-sibling is given by:
    
    P_HSPL <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
    # in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_HSPL <- t2
  }  
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF FULL & HALF-SIBLING PAIR DATA:                      ##
  #####################################################################################
  
  # Record number of adults sampled on each day:
  if (numSampledAdults > 0) {
    dailySampledAdults <- 0
    for (i in 1:numSamplingDays) {
      dailySampledAdults[i] <- sum(sampledAdults$Time == samplingDays[i])
    }
  }
  
  # Record number of larvae sampled on each day:
  if (numSampledLarvae > 0) {
    dailySampledLarvae <- 0
    for (i in 1:numSamplingDays) {
      dailySampledLarvae[i] <- sum(sampledLarvae$Time == samplingDays[i])
    }
  }
  
  # Record number of pupae sampled on each day:
  if (numSampledPupae > 0) {
    dailySampledPupae <- 0
    for (i in 1:numSamplingDays) {
      dailySampledPupae[i] <- sum(sampledPupae$Time == samplingDays[i])
    }
  }
  
  # Initialize log likelihood:
  logLike <- 0
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF LARVA-ADULT SIBLING PAIR DATA:                      ##
  #####################################################################################
  
  if ((numSampledAdults > 0) && (numSampledLarvae > 0)) {
    allRelativeTimes <- seq((T_P - 2*T_A), (T_L + T_P + 3*T_A), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of larva-adult pair data:")
    for (i in 1:numSampledLarvae) {
      # Tally number of full-siblings a given sampled larva has on each day relative to
      # its capture:
      fullSibRelativeTimes <- (sampledAdults[(which(sampledAdults$momID
                                                    == sampledLarvae[i, "momID"])),]$Time
                               - sampledLarvae[i,]$Time)
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Tally number of half-siblings a given sampled larva has on each day relative to
      # its capture:
      sameDad <- which(sampledAdults$dadID == sampledLarvae[i, "dadID"])
      sameMum <- which(sampledAdults$momID == sampledLarvae[i, "momID"])
      halfSibs <- setdiff(sameDad, sameMum)
      halfSibRelativeTimes <- (sampledAdults[halfSibs,]$Time - sampledLarvae[i,]$Time)
      halfSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        halfSibRelativeTimesTally[k] <- sum(halfSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Calculate likelihood of observed full-sibling and half-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledLarvae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledLarvae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          numHalfSibsM <- halfSibRelativeTimesTally[m]
          
          dailySampledAdultsM <- dailySampledAdults[allRelativeTimes[m] + 
                                                      sampledLarvae[i,]$Time -
                                                      tSamplingStart + 1]
          
          if (dailySampledAdultsM > 0) {
            if (P_FSLA[m]==0) {
              logLike <- (logLike + numHalfSibsM * log(P_HSLA[m]) + 
                            (dailySampledAdultsM - numHalfSibsM)
                          * log(1 - P_HSLA[m]))
            } else {
              logLike <- (logLike + numFullSibsM * log(P_FSLA[m]) + 
                            numHalfSibsM * log(P_HSLA[m]) + 
                            (dailySampledAdultsM - numFullSibsM - numHalfSibsM)
                          * log(1 - P_FSLA[m] - P_HSLA[m]))
            }
          }
        }
      }
      if (i%%100 == 0) { print(i/numSampledLarvae) }
    }
  }
      
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF ADULT-LARVA SIBLING PAIR DATA:                      ##
  #####################################################################################
  
  if ((numSampledAdults > 0) && (numSampledLarvae > 0)) {
    allRelativeTimes <- seq((- T_L - T_P - 3*T_A), (2*T_A - T_P), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of adult-larva pair data:")
    for (i in 1:numSampledAdults) {
      # Tally number of full-siblings & half-siblings a given sampled adult has on 
      # each day relative to its capture:
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      halfSibRelativeTimesTally <- rep(0,numRelativeTimes)
      
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledAdults[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledAdults[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          numHalfSibsM <- halfSibRelativeTimesTally[m]
          
          dailySampledLarvaeM <- dailySampledLarvae[allRelativeTimes[m] + 
                                                      sampledAdults[i,]$Time -
                                                      tSamplingStart + 1]

          # Calculate likelihood of observed full-sibling and half-sibling pair data:
          if (dailySampledLarvaeM > 0) {
            if (P_FSAL[m]==0) {
              logLike <- (logLike + numHalfSibsM * log(P_HSAL[m]) + 
                            (dailySampledLarvaeM - numHalfSibsM)
                          * log(1 - P_HSAL[m]))
            } else {
              logLike <- (logLike + numFullSibsM * log(P_FSAL[m]) + 
                            numHalfSibsM * log(P_HSAL[m]) + 
                            (dailySampledLarvaeM - numFullSibsM - numHalfSibsM)
                          * log(1 - P_FSAL[m] - P_HSAL[m]))
            }
          }
        }
      }
      if (i%%100 == 0) { print(i/numSampledAdults) }
    }
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF LARVA-LARVA SIBLING PAIR DATA:                      ##
  #####################################################################################
  
  if (numSampledLarvae > 0) {
    allRelativeTimes <- seq((-2*T_A - T_L), (2*T_A + T_L), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of larva-larva pair data:")
    for (i in 1:(numSampledLarvae-1)) {
      # Tally number of full-siblings a given sampled larva has on each day relative to
      # its capture:
      fullSibRelativeTimes <- (sampledLarvae[(which(sampledLarvae[(i+1):numSampledLarvae,]$momID
                                                    == sampledLarvae[i, "momID"]) + i),]$Time
                               - sampledLarvae[i,]$Time)
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Tally number of half-siblings a given sampled larva has on each day relative to
      # its capture:
      sameDad <- which(sampledLarvae[(i+1):numSampledLarvae,]$dadID == sampledLarvae[i, "dadID"])
      sameMum <- which(sampledLarvae[(i+1):numSampledLarvae,]$momID == sampledLarvae[i, "momID"])
      halfSibs <- setdiff(sameDad, sameMum)
      halfSibRelativeTimes <- (sampledLarvae[(halfSibs + i),]$Time - sampledLarvae[i,]$Time)
      halfSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        halfSibRelativeTimesTally[k] <- sum(halfSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Calculate likelihood of observed full-sibling and half-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledLarvae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledLarvae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          numHalfSibsM <- halfSibRelativeTimesTally[m]
          
          dailySampledLarvaeM <- dailySampledLarvae[allRelativeTimes[m] + 
                                                      sampledLarvae[i,]$Time -
                                                      tSamplingStart + 1]
          # Remove entry of self as a sibling:
          relativeTime0 <- which(allRelativeTimes==0)
          if (m == relativeTime0) {
            dailySampledLarvaeM <- dailySampledLarvaeM - 1
          }
          
          if (dailySampledLarvaeM > 0) {
            if (P_FSLL[m]==0) {
              logLike <- (logLike + numHalfSibsM * log(P_HSLL[m]) + 
                            (dailySampledLarvaeM - numHalfSibsM)
                          * log(1 - P_HSLL[m]))
            } else {
              logLike <- (logLike + numFullSibsM * log(P_FSLL[m]) + 
                            numHalfSibsM * log(P_HSLL[m]) + 
                            (dailySampledLarvaeM - numFullSibsM - numHalfSibsM)
                          * log(1 - P_FSLL[m] - P_HSLL[m]))
            }
          }
        }
      }
      if (i%%100 == 0) { print(i/numSampledLarvae) }
    }
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF ADULT-ADULT SIBLING PAIR DATA:                      ##
  #####################################################################################
  
  if (numSampledAdults > 0) {
    allRelativeTimes <- seq((-3*T_A), (3*T_A), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of adult-adult pair data:")
    for (i in 1:(numSampledAdults-1)) {
      # Tally number of full-siblings a given sampled adult has on each day relative to
      # its capture:
      fullSibRelativeTimes <- (sampledAdults[(which(sampledAdults[(i+1):numSampledAdults,]$momID
                                                    == sampledAdults[i, "momID"]) + i),]$Time
                               - sampledAdults[i,]$Time)
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Tally number of half-siblings a given sampled adult has on each day relative to
      # its capture:
      sameDad <- which(sampledAdults[(i+1):numSampledAdults,]$dadID == sampledAdults[i, "dadID"])
      sameMum <- which(sampledAdults[(i+1):numSampledAdults,]$momID == sampledAdults[i, "momID"])
      halfSibs <- setdiff(sameDad, sameMum)
      halfSibRelativeTimes <- (sampledAdults[(halfSibs + i),]$Time - sampledAdults[i,]$Time)
      halfSibRelativeTimesTally <- rep(0,numRelativeTimes)
      
      for (k in 1:numRelativeTimes) {
        halfSibRelativeTimesTally[k] <- sum(halfSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Calculate likelihood of observed full-sibling and half-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledAdults[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledAdults[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          numHalfSibsM <- halfSibRelativeTimesTally[m]
          
          dailySampledAdultsM <- dailySampledAdults[allRelativeTimes[m] + 
                                                      sampledAdults[i,]$Time -
                                                      tSamplingStart + 1]
          # Remove entry of self as a sibling:
          relativeTime0 <- which(allRelativeTimes==0)
          if (m == relativeTime0) {
            dailySampledAdultsM <- dailySampledAdultsM - 1
          }
          
          if (dailySampledAdultsM > 0) {
            if (P_FSAA[m]==0) {
              logLike <- (logLike + numHalfSibsM * log(P_HSAA[m]) + 
                            (dailySampledAdultsM - numHalfSibsM)
                          * log(1 - P_HSAA[m]))
            } else {
              logLike <- (logLike + numFullSibsM * log(P_FSAA[m]) + 
                            numHalfSibsM * log(P_HSAA[m]) + 
                            (dailySampledAdultsM - numFullSibsM - numHalfSibsM)
                          * log(1 - P_FSAA[m] - P_HSAA[m]))
            }
          }
        }
      }
      if (i%%100 == 0) { print(i/numSampledAdults) }
    }
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF LARVA-PUPA SIBLING PAIR DATA:                       ##
  #####################################################################################
  
  if ((numSampledLarvae > 0) && (numSampledPupae > 0)) {
    allRelativeTimes <- seq((- 2*T_A), (T_L + T_P + 2*T_A), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of larva-pupa pair data:")
    for (i in 1:numSampledLarvae) {
      # Tally number of full-siblings a given sampled larva has on each day relative to
      # its capture:
      fullSibRelativeTimes <- (sampledPupae[(which(sampledPupae$momID
                                                   == sampledLarvae[i, "momID"])),]$Time
                               - sampledLarvae[i,]$Time)
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Tally number of half-siblings a given sampled larva has on each day relative to
      # its capture:
      sameDad <- which(sampledPupae$dadID == sampledLarvae[i, "dadID"])
      sameMum <- which(sampledPupae$momID == sampledLarvae[i, "momID"])
      halfSibs <- setdiff(sameDad, sameMum)
      halfSibRelativeTimes <- (sampledPupae[halfSibs,]$Time - sampledLarvae[i,]$Time)
      halfSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        halfSibRelativeTimesTally[k] <- sum(halfSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Calculate likelihood of observed full-sibling and half-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledLarvae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledLarvae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          numHalfSibsM <- halfSibRelativeTimesTally[m]
          
          dailySampledPupaeM <- dailySampledPupae[allRelativeTimes[m] + 
                                                    sampledLarvae[i,]$Time -
                                                    tSamplingStart + 1]
          
          if (dailySampledPupaeM > 0) {
            if (P_FSLP[m]==0) {
              logLike <- (logLike + numHalfSibsM * log(P_HSLP[m]) + 
                            (dailySampledPupaeM - numHalfSibsM)
                          * log(1 - P_HSLP[m]))
            } else {
              logLike <- (logLike + numFullSibsM * log(P_FSLP[m]) + 
                            numHalfSibsM * log(P_HSLP[m]) + 
                            (dailySampledPupaeM - numFullSibsM - numHalfSibsM)
                          * log(1 - P_FSLP[m] - P_HSLP[m]))
            }
          }
        }
      }
      if (i%%100 == 0) { print(i/numSampledLarvae) }
    }
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF PUPA-LARVA SIBLING PAIR DATA:                       ##
  #####################################################################################
  
  if ((numSampledLarvae > 0) && (numSampledPupae > 0)) {
    allRelativeTimes <- seq((- T_L - T_P - 2*T_A), (2*T_A), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of pupa-larva pair data:")
    for (i in 1:numSampledPupae) {
      # Tally number of full-siblings & half-siblings a given sampled pupa has on 
      # each day relative to its capture:
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      halfSibRelativeTimesTally <- rep(0,numRelativeTimes)
      
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledPupae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledPupae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          numHalfSibsM <- halfSibRelativeTimesTally[m]
          
          dailySampledLarvaeM <- dailySampledLarvae[allRelativeTimes[m] + 
                                                      sampledPupae[i,]$Time -
                                                      tSamplingStart + 1]
          
          # Calculate likelihood of observed full-sibling and half-sibling pair data:
          if (dailySampledLarvaeM > 0) {
            if (P_FSPL[m]==0) {
              logLike <- (logLike + numHalfSibsM * log(P_HSPL[m]) + 
                            (dailySampledLarvaeM - numHalfSibsM)
                          * log(1 - P_HSPL[m]))
            } else {
              logLike <- (logLike + numFullSibsM * log(P_FSPL[m]) + 
                            numHalfSibsM * log(P_HSPL[m]) + 
                            (dailySampledLarvaeM - numFullSibsM - numHalfSibsM)
                          * log(1 - P_FSPL[m] - P_HSPL[m]))
            }
          }
        }
      }
      if (i%%100 == 0) { print(i/numSampledPupae) }
    }
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF PUPA-PUPA SIBLING PAIR DATA:                        ##
  #####################################################################################
  
  if (numSampledPupae > 0) {
    allRelativeTimes <- seq((-2*T_A - T_P), (2*T_A + T_P), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of pupa-pupa pair data:")
    for (i in 1:(numSampledPupae-1)) {
      # Tally number of full-siblings a given sampled pupa has on each day relative to
      # its capture:
      fullSibRelativeTimes <- (sampledPupae[(which(sampledPupae[(i+1):numSampledPupae,]$momID
                                                   == sampledPupae[i, "momID"]) + i),]$Time
                               - sampledPupae[i,]$Time)
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Tally number of half-siblings a given sampled pupa has on each day relative to
      # its capture:
      sameDad <- which(sampledPupae[(i+1):numSampledPupae,]$dadID == sampledPupae[i, "dadID"])
      sameMum <- which(sampledPupae[(i+1):numSampledPupae,]$momID == sampledPupae[i, "momID"])
      halfSibs <- setdiff(sameDad, sameMum)
      halfSibRelativeTimes <- (sampledPupae[(halfSibs + i),]$Time - sampledPupae[i,]$Time)
      halfSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        halfSibRelativeTimesTally[k] <- sum(halfSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Calculate likelihood of observed full-sibling and half-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledPupae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledPupae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          numHalfSibsM <- halfSibRelativeTimesTally[m]
          
          dailySampledPupaeM <- dailySampledPupae[allRelativeTimes[m] + 
                                                    sampledPupae[i,]$Time -
                                                    tSamplingStart + 1]
          # Remove entry of self as a sibling:
          relativeTime0 <- which(allRelativeTimes==0)
          if (m == relativeTime0) {
            dailySampledPupaeM <- dailySampledPupaeM - 1
          }
          
          if (dailySampledPupaeM > 0) {
            if (P_FSPP[m]==0) {
              logLike <- (logLike + numHalfSibsM * log(P_HSPP[m]) + 
                            (dailySampledPupaeM - numHalfSibsM)
                          * log(1 - P_HSPP[m]))
            } else {
              logLike <- (logLike + numFullSibsM * log(P_FSPP[m]) + 
                            numHalfSibsM * log(P_HSPP[m]) + 
                            (dailySampledPupaeM - numFullSibsM - numHalfSibsM)
                          * log(1 - P_FSPP[m] - P_HSPP[m]))
            }
          }
        }
      }
      if (i%%100 == 0) { print(i/numSampledPupae) }
    }
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF PUPA-ADULT SIBLING PAIR DATA:                       ##
  #####################################################################################
  
  if ((numSampledAdults > 0) && (numSampledPupae > 0)) {
    allRelativeTimes <- seq((- 2*T_A), (T_P + 3*T_A), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of pupa-adult pair data:")
    for (i in 1:numSampledPupae) {
      # Tally number of full-siblings a given sampled pupa has on each day relative to
      # its capture:
      fullSibRelativeTimes <- (sampledAdults[(which(sampledAdults$momID
                                                    == sampledPupae[i, "momID"])),]$Time
                               - sampledPupae[i,]$Time)
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Tally number of half-siblings a given sampled pupa has on each day relative to
      # its capture:
      sameDad <- which(sampledAdults$dadID == sampledPupae[i, "dadID"])
      sameMum <- which(sampledAdults$momID == sampledPupae[i, "momID"])
      halfSibs <- setdiff(sameDad, sameMum)
      halfSibRelativeTimes <- (sampledAdults[halfSibs,]$Time - sampledPupae[i,]$Time)
      halfSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        halfSibRelativeTimesTally[k] <- sum(halfSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Calculate likelihood of observed full-sibling and half-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledPupae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledPupae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          numHalfSibsM <- halfSibRelativeTimesTally[m]
          
          dailySampledAdultsM <- dailySampledAdults[allRelativeTimes[m] + 
                                                      sampledPupae[i,]$Time -
                                                      tSamplingStart + 1]
          
          if (dailySampledAdultsM > 0) {
            if (P_FSPA[m]==0) {
              logLike <- (logLike + numHalfSibsM * log(P_HSPA[m]) + 
                            (dailySampledAdultsM - numHalfSibsM)
                          * log(1 - P_HSPA[m]))
            } else {
              logLike <- (logLike + numFullSibsM * log(P_FSPA[m]) + 
                            numHalfSibsM * log(P_HSPA[m]) + 
                            (dailySampledAdultsM - numFullSibsM - numHalfSibsM)
                          * log(1 - P_FSPA[m] - P_HSPA[m]))
            }
          }
        }
      }
      if (i%%100 == 0) { print(i/numSampledPupae) }
    }
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF ADULT-PUPA SIBLING PAIR DATA:                       ##
  #####################################################################################
  
  if ((numSampledAdults > 0) && (numSampledPupae > 0)) {
    allRelativeTimes <- seq((- T_P - 3*T_A), (2*T_A), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of adult-pupa pair data:")
    for (i in 1:numSampledAdults) {
      # Tally number of full-siblings & half-siblings a given sampled adult has on 
      # each day relative to its capture:
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      halfSibRelativeTimesTally <- rep(0,numRelativeTimes)
      
      # Calculate likelihood of observed full-sibling and half-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledAdults[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledAdults[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          numHalfSibsM <- halfSibRelativeTimesTally[m]
          
          dailySampledPupaeM <- dailySampledPupae[allRelativeTimes[m] + 
                                                    sampledAdults[i,]$Time -
                                                    tSamplingStart + 1]
          
          if (dailySampledPupaeM > 0) {
            if (P_FSAP[m]==0) {
              logLike <- (logLike + numHalfSibsM * log(P_HSAP[m]) + 
                            (dailySampledPupaeM - numHalfSibsM)
                          * log(1 - P_HSAP[m]))
            } else {
              logLike <- (logLike + numFullSibsM * log(P_FSAP[m]) + 
                            numHalfSibsM * log(P_HSAP[m]) + 
                            (dailySampledPupaeM - numFullSibsM - numHalfSibsM)
                          * log(1 - P_FSAP[m] - P_HSAP[m]))
            }
          }
        }
      }
      if (i%%100 == 0) { print(i/numSampledAdults) }
    }
  }  
  
  -logLike
}

#######################################################################################
# MAXIMUM-LIKELIHOOD PARAMETER INFERENCE:                                             #
#######################################################################################

# Log likelihood of all parent-offspring & full-sibling pair data:

logLike_all <- function(AdultPars) {
  
  # Extract parameters that we are varying from the AdultPars vector:  
  N_A <- AdultPars[1] * OptimScale
  mu_A <- AdultPars[2]
  
  # Parameters that we are setting constant in the life history model (comment
  # out the ones you are trying to estimate):
  # N_A <- 3000 # Total adult mosquito population size (females & males)
  # mu_A <- 0.09 # Daily mortality of adult mosquitoes
  N_F <- N_A/2 # Total adult female population size
  
  T_E <- 2 # Duration of the egg stage (days)
  T_L <- 5 # Duration of the larval stage (days)
  T_P <- 1 # Duration of the pupal stage (days)
  T_A <- 30 # Maximum adult lifespan considered (days)
  
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
  
  logLike <- 0 # Initialize the log-likelihood
  
  if ((numSampledLarvae > 0) && (numSampledAdultFemales > 0)) { 
    logLike <- logLike + logLike_MOL(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                     beta, mu_E, mu_L, mu_P) 
  }
  if ((numSampledLarvae > 0) && (numSampledAdultMales > 0)) { 
    logLike <- logLike + logLike_FOL(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                     beta, mu_E, mu_L, mu_P)
  }
  if ((numSampledPupae > 0) && (numSampledAdultFemales > 0)) { 
    logLike <- logLike + logLike_MOP(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                     beta, mu_E, mu_L, mu_P)
  }
  if ((numSampledPupae > 0) && (numSampledAdultMales > 0)) { 
    logLike <- logLike + logLike_FOP(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                     beta, mu_E, mu_L, mu_P)
  }
  if (numSampledAdultFemales > 0) { 
    logLike <- logLike + logLike_MOA(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                     beta, mu_E, mu_L, mu_P) 
  }
  if (numSampledAdultMales > 0) { 
    logLike <- logLike + logLike_FOA(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                     beta, mu_E, mu_L, mu_P)
  }
  logLike <- logLike + logLike_Sibs(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                                    beta, mu_E, mu_L, mu_P)
  
  logLike
}

# Choose initial values of N_A & mu_A for the optimization algorithm:
N_A <- 2500
mu_A <- 0.1
OptimScale <- N_A/mu_A
AdultPars <- c(N_A / OptimScale, mu_A)

# Choose lower & upper limits of N_A & mu_A for optimization algorithm:
N_A_lower <- 1000
mu_A_lower <- 0.05
AdultPars_lower <- c(N_A_lower / OptimScale, mu_A_lower)

N_A_upper <- 5000
mu_A_upper <- 0.15
AdultPars_upper <- c(N_A_upper / OptimScale, mu_A_upper)

# Find the values of N_A & mu_A that maximize the likelihood of the parent-offspring
# data:
AdultPars_fit <- optimx(par=AdultPars, fn=logLike_all, method = "nlminb", 
                        upper = AdultPars_upper, lower = AdultPars_lower)

# The estimated values of the adult parameters given the kinship data are:
N_A_fit <- AdultPars_fit$p1 * OptimScale
mu_A_fit <- AdultPars_fit$p2

N_A_fit
mu_A_fit
