#######################################################################################
## CLOSE-KIN MARK-RECAPTURE CODE FOR MOSQUITOES:                                     ##
##                                                                                   ##
## This code calculates the kinship probabilities and likelihood of observed close-  ##
## kin pairs considering mother-larval offspring and mother-adult offspring          ##
## combinations when the day of capture is only known within the interval between    ##
## samples. The code then implements a grid search algorithm to estimate the daily   ##
## larval mortality rate (mu_L) and duration of the larval life stage (T_L) that     ##
## maximize the likelihood given the data.                                           ##
##                                                                                   ##
## Code written by John Marshall: john.marshall@berkeley.edu                         ##
## Date: July 15th, 2022                                                             ##
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
  samplingDaysLarvae <- sort(unique(sampledLarvae$Time))
  numSamplingDaysLarvae <- length(samplingDaysLarvae)
}

# if (file.exists("cut_P.csv")) { 
#  sampledPupae <- read.csv("cut_P.csv") 
#  numSampledPupae <- length(sampledPupae[,1])
#  samplingDaysPupae <- sort(unique(sampledPupae$Time))
#  numSamplingDaysPupae <- length(samplingDaysPupae)
# }

if (file.exists("cut_M.csv")) { 
  sampledAdultMales <- read.csv("cut_M.csv") 
  numSampledAdultMales <- length(sampledAdultMales[,1])
  samplingDaysAdultMales <- sort(unique(sampledAdultMales$Time))
  numSamplingDaysAdultMales <- length(samplingDaysAdultMales)
}

if (file.exists("cut_F.csv")) { 
  sampledAdultFemales <- read.csv("cut_F.csv") 
  numSampledAdultFemales <- length(sampledAdultFemales[,1])
  samplingDaysAdultFemales <- sort(unique(sampledAdultFemales$Time))
  numSamplingDaysAdultFemales <- length(samplingDaysAdultFemales)
}

numSampledAdults <- numSampledAdultMales + numSampledAdultFemales

# Concatenate adult female & adult male data into one adult data set:
if (numSampledAdultMales > 0) {
  sampledAdultMales$Mate <- NA
}
if ((numSampledAdultMales > 0) && (numSampledAdultFemales > 0)) {
  sampledAdults <- rbind(sampledAdultFemales,sampledAdultMales)
  samplingDaysAdults <- sort(unique(c(sampledAdultMales$Time, sampledAdultFemales$Time)))
  numSamplingDaysAdults <- length(samplingDaysAdults)
}
if ((numSampledAdultMales > 0) && (numSampledAdultFemales == 0)) {
  sampledAdults <- sampledAdultMales
  samplingDaysAdults <- samplingDaysAdultMales
  numSamplingDaysAdults <- length(samplingDaysAdults)
}
if ((numSampledAdultMales == 0) && (numSampledAdultFemales > 0)) {
  sampledAdults <- sampledAdultFemales
  samplingDaysAdults <- samplingDaysAdultFemales
  numSamplingDaysAdults <- length(samplingDaysAdults)
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

logLike_MOL <- function(mu_L) {
  
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
  for (y2 in (t2 - T_E - (T_L-1)):(t2 - T_E)) {
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^(t2 - y2 - T_E)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of larvae at day t2 from an adult female sampled at 
  # time t1.
  # * By default, let t1 = 0, as the same equation will apply at all times.
  # * Earliest possible t2 is (-(T_A-1) + T_E), if the mother was caught at the end of her
  #   life & gave birth to the offspring soon after emergence.
  # * Latest possible t2 is (T_E + (T_L-1)), if the the mother gave birth at the time of 
  #   sampling & the larva was caught at the end of its life.
  # * So we will explore (-(T_A-1) + T_E) <= t2 <= (T_E + (T_L-1))
  
  Numerator <- rep(0, (abs(-(T_A-1) + T_E) + (T_E + (T_L-1)) + 1))
  
  # Probability of adult surviving from 0 to T_A days:
  AdultSurvivalProbability <- rep(0, T_A)
  AdultAge <- rep(0, T_A)
  for (i in 1:T_A) {
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and T_A days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-(T_A-1) + T_E), (T_E + (T_L-1)), by=1)
  
  for (i in 1:length(t2)) { 
    for (y2 in (t2[i] - T_E - (T_L-1)):(t2[i] - T_E)) {
      if ((y2 >= (t1 - (T_A-1))) && (y2 <= t1)) {
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

  # Buffer days account for uncertainty over day of sampling:
  bufferDays_L <- (samplingDaysLarvae[3] - samplingDaysLarvae[2] - 1)
  bufferDays_M <- (samplingDaysAdultFemales[3] - samplingDaysAdultFemales[2] - 1)
  
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
          
          # Earliest day that the larva could have been sampled:
          if (larvaSamplingTimeI == min(samplingDaysLarvae)) {
            tEarliest_L <- larvaSamplingTimeI - bufferDays_L
          } else {
            tEarliest_L <- max(samplingDaysLarvae[samplingDaysLarvae<larvaSamplingTimeI]) + 1
          }
          
          # Latest day that the larva could have been sampled:
          tLatest_L <- larvaSamplingTimeI
          
          # Number of days over which the larva could have been sampled:
          numDaysSampled_L <- tLatest_L - tEarliest_L + 1
          
          # Earliest day that the adult female (mother) could have been sampled:
          if (motherSamplingTimeJ == min(samplingDaysAdultFemales)) {
            tEarliest_M <- motherSamplingTimeJ - bufferDays_M
          } else {
            tEarliest_M <- max(samplingDaysAdultFemales[samplingDaysAdultFemales<motherSamplingTimeJ]) + 1
          }
          
          # Latest day that the adult female (mother) could have been sampled:
          tLatest_M <- motherSamplingTimeJ
          
          # Number of days over which the adult female (mother) could have been sampled:
          numDaysSampled_M <- tLatest_M - tEarliest_M + 1
          
          P_MOL_mean <- 0
          for (tM in tEarliest_M:tLatest_M) {
            for (tL in tEarliest_L:tLatest_L) {
              if ((tL-tM) %in% t2Minust1) {
                P_MOL_mean <- P_MOL_mean + (P_MOL[which(t2Minust1==(tL-tM))]/(numDaysSampled_M * numDaysSampled_L))
              }
            }
          }
          
          # Probability that a given sampled larva on day t2 has a mother among the y
          # sampled adult females on day t1:
          z <- 1 - ((1 - P_MOL_mean)^y)
          
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

logLike_MOA <- function(T_L) {
  
  #####################################################################################
  ## MOTHER-OFFSPRING (ADULT) KINSHIP PROBABILITIES:                                 ##
  #####################################################################################
  
  # Given an adult sampled at time t2, this is the probability that an adult female
  # sampled at time t1 is their mother.
  
  # First, calculate the denominator:
  # This is the expected number of surviving adults at time t2 from adult females at any
  # consistent time (assuming a constant population size, this is independent of time).
  
  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2 (use 0 as an example)
  for (y2 in (t2 - T_E - T_L - T_P - (T_A-1)):(t2 - T_E - T_L - T_P)) {
    Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                  * ((1 - mu_P)^T_P) 
                                  * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
  }
  
  # Next, calculate the numerator:
  # This is the expected number of adults at day t2 from an adult female sampled at 
  # time t1.
  # * By default, let t1 = 0, as the same equation will apply at all times (it's the
  #   difference between t1 & t2 that matters).
  # * Earliest possible t2 is (-(T_A-1) + T_E + T_L + T_P), if the mother was caught at the 
  #   end of her life & gave birth to the offspring soon after emergence.
  # * Latest possible t2 is (T_E + T_L + + T_P + (T_A-1)), if the the mother gave birth at 
  #   the time of sampling & the adult offspring was caught at the end of its life.
  # * So we will explore (-(T_A-1) + T_E + T_L + T_P) <= t2 <= (T_E + T_L + T_P + (T_A-1)).
  
  Numerator <- rep(0, (abs(-(T_A-1) + T_E + T_L + T_P) + (T_E + T_L + T_P + (T_A-1)) + 1))
  
  # Probability of adult surviving from 0 to T_A days:
  AdultSurvivalProbability <- rep(0, T_A)
  AdultAge <- rep(0, T_A)
  for (i in 1:T_A) {
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and T_A days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-(T_A-1) + T_E + T_L + T_P), (T_E + T_L + T_P + (T_A-1)), by=1)
  
  for (i in 1:length(t2)) { 
    for (y2 in (t2[i] - T_E - T_L - T_P - (T_A-1)):(t2[i] - T_E - T_L - T_P)) {
      if ((y2 >= (t1 - (T_A-1))) && (y2 <= t1)) {
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
      adultID[j] <- sampledAdults[i, "myID"] # Adult offspring ID
      motherID[j] <- sampledAdults[i, "momID"] # Mother ID
      adultSamplingTime[j] <- sampledAdults[i, "Time"] # Day adult offspring sampled
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
  
  # Buffer days account for uncertainty over day of sampling:
  bufferDays_M <- (samplingDaysAdultFemales[3] - samplingDaysAdultFemales[2] - 1)
  bufferDays_A <- (samplingDaysAdults[3] - samplingDaysAdults[2] - 1)
  
  # Format data into a 2D array containing the number of sampled adults, w, on day t2 
  # (first index) that have mothers among the sampled adult females on day t2-t2Minust1
  # (second index). Also calculate the log likelihood that, for x sampled adults on day
  # t2, w have mothers among the y sampled adult females on day t1:
  
  MOA_Data_Array <- array(rep(NA, length(samplingDays) * length(t2Minust1)),
                          dim=c(length(samplingDays), length(t2Minust1)))
  logLike <- 0
  
  for (i in (min(which(dailySampledAdults>0))):length(samplingDays)) {
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
          
          # Earliest day that the adult offspring could have been sampled:
          if (adultSamplingTimeI == min(samplingDaysAdults)) {
            tEarliest_A <- adultSamplingTimeI - bufferDays_A
          } else {
            tEarliest_A <- max(samplingDaysAdults[samplingDaysAdults<adultSamplingTimeI]) + 1
          }
          
          # Latest day that the adult offspring could have been sampled:
          tLatest_A <- adultSamplingTimeI
          
          # Number of days over which the adult offsrping could have been sampled:
          numDaysSampled_A <- tLatest_A - tEarliest_A + 1
          
          # Earliest day that the adult female (mother) could have been sampled:
          if (motherSamplingTimeJ == min(samplingDaysAdultFemales)) {
            tEarliest_M <- motherSamplingTimeJ - bufferDays_M
          } else {
            tEarliest_M <- max(samplingDaysAdultFemales[samplingDaysAdultFemales<motherSamplingTimeJ]) + 1
          }
          
          # Latest day that the adult female (mother) could have been sampled:
          tLatest_M <- motherSamplingTimeJ
          
          # Number of days over which the adult female (mother) could have been sampled:
          numDaysSampled_M <- tLatest_M - tEarliest_M + 1
          
          P_MOA_mean <- 0
          for (tM in tEarliest_M:tLatest_M) {
            for (tA in tEarliest_A:tLatest_A) {
              if ((tA-tM) %in% t2Minust1) {
                P_MOA_mean <- P_MOA_mean + (P_MOA[which(t2Minust1==(tA-tM))]/(numDaysSampled_M * numDaysSampled_A))
              }
            }
          }
          
          # Probability that a given sampled adult on day t2 has a mother among the y
          # sampled adult females on day t1:
          z <- 1 - ((1 - P_MOA_mean)^y)
          
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
# MAXIMUM-LIKELIHOOD PARAMETER INFERENCE:                                             #
#######################################################################################

# Vary T_L from 1 to 10. For each value of T_L, find the value of mu_L that maximizes
# the likelihood of the mother-larval offspring data. Then, calculate the likelihood of 
# the mother-adult offspring data given that value of T_L & the optimal value of mu_L.

num_T_L_values <- 10
T_L_values <- seq(from=1, to=10, length.out=num_T_L_values)
mu_L_values <- rep(0, 10)
logLike_values <- rep(0, num_T_L_values)

for (iTL in 1:num_T_L_values) {
  print("Iteration:")
  print(iTL)
  
  # The value of T_L for this iteration:
  T_L <- T_L_values[iTL]
  
  # Parameters that we are setting constant in the life history model (comment
  # out the ones you are trying to estimate):
  N_A <- 3000 # Total adult mosquito population size (females & males)
  mu_A <- 0.09 # Daily mortality of adult mosquitoes
  N_F <- N_A/2 # Total adult female population size
  
  T_E <- 2 # Duration of the egg stage (days)
  # T_L <- 5 # Duration of the larval stage (days)
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
  # mu_L <- 1 - (2*mu_A/(beta*(1-mu_A)*((1-mu_J)^(T_E+T_P))))^(1/T_L)
  
  # Choose an initial value of mu_L for the optimization algorithm:
  mu_L <- 0.5
  
  # Find the value of mu_L that maximizes the likelihood of the mother-larval offspring data:
  mu_L_optim <- optimx(par=mu_L, fn=logLike_MOL, method = "nlminb", upper = 0.96, lower = 0.01)
  mu_L_values[iTL] <- mu_L_optim$p1
  
  # Update mu_L as the maximum-likelihood estimate given this value of T_L:
  mu_L <- mu_L_optim$p1
  
  # Calculate the likelihood of the mother-adult offspring data given this value 
  # of T_L & the optimal value of mu_L:
  logLike_values[iTL] <- logLike_MOA(T_L_values[iTL])
}

# The estimated values of the larval parameters given the kinship data are:
T_L_fit <- T_L_values[which.min(logLike_values)]
mu_L_fit <- mu_L_values[which.min(logLike_values)]

T_L_fit
mu_L_fit
 