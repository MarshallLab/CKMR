#######################################################################################
## CLOSE-KIN MARK-RECAPTURE CODE FOR MOSQUITOES:                                     ##
##                                                                                   ##
## This code calculates the kinship probabilities and likelihood of observed close-  ##
## kin pairs considering mother-offspring (adult) and full-sibling (adult-adult)     ##
## pairs when the day of capture is only known within the interval between samples.  ##
## The code then uses an optimization algorithm to find the adult parameters - adult ##
## census population size (N_A) and daily adult mortality rate (mu_A) that maximize  ##
## the likelihood given the data.                                                    ##
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

if (file.exists("cut_P.csv")) { 
  sampledPupae <- read.csv("cut_P.csv") 
  numSampledPupae <- length(sampledPupae[,1])
  samplingDaysPupae <- sort(unique(sampledPupae$Time))
  numSamplingDaysPupae <- length(samplingDaysPupae)
}

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
  
  # Buffer days account for uncertainty over day of sampling:
  bufferDays <- (samplingDaysAdultFemales[2] - samplingDaysAdultFemales[1] - 1)
  
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
          
          # Earliest day that the adult offspring could have been sampled:
          if (adultSamplingTimeI == min(samplingDaysAdults)) {
            tEarliest_A <- adultSamplingTimeI - bufferDays
          } else {
            tEarliest_A <- max(samplingDaysAdults[samplingDaysAdults<adultSamplingTimeI]) + 1
          }
          
          # Latest day that the adult offspring could have been sampled:
          tLatest_A <- adultSamplingTimeI
          
          # Number of days over which the adult offsrping could have been sampled:
          numDaysSampled_A <- tLatest_A - tEarliest_A + 1
          
          # Earliest day that the adult female (mother) could have been sampled:
          if (motherSamplingTimeJ == min(samplingDaysAdultFemales)) {
            tEarliest_M <- motherSamplingTimeJ - bufferDays
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
## FULL-SIBLING KINSHIP PROBABILITIES & LIKELIHOOD:                                  ##
#######################################################################################

logLike_Sibs <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
                         beta, mu_E, mu_L, mu_P) {
  
  #####################################################################################
  ## FULL-SIBLING (ADULT-ADULT) KINSHIP PROBABILITIES:                               ##
  #####################################################################################
  
  # Given an adult sampled at time t2, this is the probability that an adult sampled at
  # time t1 is its full sibling.
  
  # Probability of larva surviving from 0 to T_L days:
  LarvaSurvivalProbability <- rep(0, T_L)
  LarvaAge <- rep(0, T_L)
  for (i in 1:T_L) {
    LarvaSurvivalProbability[i] <- (1 - mu_L)^(i-1)
    LarvaAge[i] <- i - 1
  }
  
  # Probability of larva having an age between 0 and T_L days:
  LarvaAgeProbability <- LarvaSurvivalProbability / sum(LarvaSurvivalProbability)
  
  # Probability of adult surviving from 0 to T_A days:
  AdultSurvivalProbability <- rep(0, T_A)
  AdultAge <- rep(0, T_A)
  for (i in 1:T_A) {
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i-1)
    AdultAge[i] <- i - 1
  }
  
  # Probability of adult having an age between 0 and T_A days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  
  if (numSampledAdults > 0) {
    
    # First, calculate the denominator:
    # This is the expected number of surviving adults at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).
    
    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P - (T_A-1)):(t2 - T_E - T_L - T_P)) {
      Denominator <- Denominator + (N_F * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L)
                                    * ((1 - mu_P)^T_P) 
                                    * ((1 - mu_A)^(t2 - y2 - T_E - T_L - T_P)))
    }
    
    # Next, calculate the numerator:
    # This is the expected number of adults at day t2 that are full siblings of an adult 
    # sampled at time t1.
    # * By default, let t1 = 0, as the same equation will apply at all times.
    # * Earliest possible t2 is (-(T_A-1) - (T_A-1)), if the mother laid egg 1 at the end of her
    #   life, egg 2 at the beginning of her life, adult 1 was caught at the end of its
    #   life, & adult 2 was caught at the beginning of its life.
    # * Latest possible t2 is ((T_A-1) + (T_A-1)), if the mother laid egg 1 at the beginning of her
    #   life, egg 2 at the end of her life, adult 1 was caught at the beginning of its
    #   life, & adult 2 was caught at the end of its life.
    # * So we will explore (-(T_A-1) - (T_A-1)) <= t2 <= ((T_A-1) + (T_A-1))
    
    Numerator <- rep(0, (abs(-(T_A-1) - (T_A-1)) + ((T_A-1) + (T_A-1)) + 1))
    
    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((-(T_A-1) - (T_A-1)), ((T_A-1) + (T_A-1)), by=1)
    
    print("Calculating full-sibling adult-adult probabilities:")
    for (i in 1:length(t2)) { 
      for (y1 in (t1 - T_E - T_L - T_P - (T_A-1)):(t1 - T_E - T_L - T_P)) {
        for (ym in (y1 - (T_A-1)):y1) {
          for (y2 in ym:(ym + (T_A-1))) {
            if ((y2 >= (t2[i] - T_E - T_L - T_P - (T_A-1))) && (y2 <= (t2[i] - T_E - T_L - T_P))) {
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
    
    P_FSAA <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSAA <- seq((-2*(T_A-1)), (2*(T_A-1)), by=1)
  }
  
  #####################################################################################
  ## CALCULATE THE LIKELIHOOD OF FULL-SIBLING PAIR DATA:                             ##
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
  ## CALCULATE THE LIKELIHOOD OF ADULT-ADULT SIBLING PAIR DATA:                      ##
  #####################################################################################
  
  # Buffer days account for uncertainty over day of sampling:
  bufferDays <- (samplingDaysAdults[2] - samplingDaysAdults[1] - 1)
  
  if (numSampledAdults > 0) {
    allRelativeTimes <- seq((-2*(T_A-1)), (2*(T_A-1)), by=1)
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
      
      # Calculate likelihood of observed full-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledAdults[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledAdults[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          
          dailySampledAdultsM <- dailySampledAdults[allRelativeTimes[m] + 
                                                      sampledAdults[i,]$Time -
                                                      tSamplingStart + 1]
          # Remove entry of self as a sibling:
          relativeTime0 <- which(allRelativeTimes==0)
          if (m == relativeTime0) {
            dailySampledAdultsM <- dailySampledAdultsM - 1
          }
          
          if (dailySampledAdultsM > 0) {
            
            # Recorded sampling time of first adult:
            adultSamplingTime1 <- sampledAdults[i,]$Time
            
            # Earliest day that first adult offspring could have been sampled:
            if (adultSamplingTime1 == min(samplingDaysAdults)) {
              tEarliest_A1 <- adultSamplingTime1 - bufferDays
            } else {
              tEarliest_A1 <- max(samplingDaysAdults[samplingDaysAdults<adultSamplingTime1]) + 1
            }
            
            # Latest day that the adult offspring could have been sampled:
            tLatest_A1 <- adultSamplingTime1
            
            # Number of days over which first adult could have been sampled:
            numDaysSampled_A1 <- tLatest_A1 - tEarliest_A1 + 1
            
            # Recorded sampling time of second adult:
            adultSamplingTime2 <- allRelativeTimes[m] + sampledAdults[i,]$Time

            # Earliest day that second adult could have been sampled:
            if (adultSamplingTime2 == min(samplingDaysAdults)) {
              tEarliest_A2 <- adultSamplingTime2 - bufferDays
            } else {
              tEarliest_A2 <- max(samplingDaysAdults[samplingDaysAdults<adultSamplingTime2]) + 1
            }
            
            # Latest day that second adult could have been sampled:
            tLatest_A2 <- adultSamplingTime2
            
            # Number of days over which second adult could have been sampled:
            numDaysSampled_A2 <- tLatest_A2 - tEarliest_A2 + 1
            
            P_FSAA_mean <- 0
            for (tA1 in tEarliest_A1:tLatest_A1) {
              for (tA2 in tEarliest_A2:tLatest_A2) {
                if ((tA2-tA1) %in% t2Minust1_FSAA) {
                  P_FSAA_mean <- P_FSAA_mean + (P_FSAA[which(t2Minust1_FSAA==(tA2-tA1))]/(numDaysSampled_A1 * numDaysSampled_A2))
                }
              }
            }
            
            logLike <- (logLike + numFullSibsM * log(P_FSAA_mean) + 
                          (dailySampledAdultsM - numFullSibsM)
                        * log(1 - P_FSAA_mean))
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
  
  logLike <- 0 # Initialize the log-likelihood
  
  if (numSampledAdultFemales > 0) { 
    logLike <- logLike + logLike_MOA(N_A, mu_A, N_F, T_E, T_L, T_P, T_A, 
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