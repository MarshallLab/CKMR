#######################################################################################
## CLOSE-KIN MARK-RECAPTURE CODE FOR MOSQUITOES:                                     ##
##                                                                                   ##
## This code calculates the kinship probabilities and likelihood of observed close-  ##
## kin pairs considering full-sibling pairs where either member of a pair is larval, ## 
## pupal or adult. The code then uses an optimization algorithm to find the adult    ##
## parameters - adult census population size (N_A) and daily adult mortality rate    ##
## (mu_A) that maximize the likelihood given the data.                               ##
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
## FULL-SIBLING KINSHIP PROBABILITIES & LIKELIHOOD:                                  ##
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
    
    P_FSLL <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSLL <- seq((-T_A - T_L), (T_A + T_L), by=1)
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
    
    P_FSAA <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSAA <- seq((-2*T_A), (2*T_A), by=1)
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
    
    P_FSAL <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSAL <- seq((- T_L - T_P - 2*T_A), (T_A - T_P), by=1)
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
    
    P_FSLA <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSLA <- seq((T_P - T_A), (T_L + T_P + 2*T_A), by=1)
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
    
    P_FSPP <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSPP <- seq((-T_A - T_P), (T_A + T_P), by=1)
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
    
    P_FSPA <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSPA <- seq((- T_A), (T_P + 2*T_A), by=1)
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
    
    P_FSAP <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSAP <- seq((- T_P - 2*T_A), (T_A), by=1)
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
    
    P_FSLP <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSLP <- seq((- T_A), (T_L + T_P + T_A), by=1)
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
    
    P_FSPL <- Numerator/Denominator
    
    # Note that indices here relate to times in the vector, t2, where t1 = 0, and 
    # hence, in general, the indices relate to times t2-t1, i.e.:
    
    t2Minust1_FSPL <- seq((- T_L - T_P - T_A), (T_A), by=1)
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
  ## CALCULATE THE LIKELIHOOD OF LARVA-ADULT FULL-SIBLING PAIR DATA:                 ##
  #####################################################################################
  
  if ((numSampledAdults > 0) && (numSampledLarvae > 0)) {
    allRelativeTimes <- seq((T_P - T_A), (T_L + T_P + 2*T_A), by=1)
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
      
      # Calculate likelihood of observed full-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledLarvae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledLarvae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          
          dailySampledAdultsM <- dailySampledAdults[allRelativeTimes[m] + 
                                                      sampledLarvae[i,]$Time -
                                                      tSamplingStart + 1]
          
          if (dailySampledAdultsM > 0) {
            logLike <- (logLike + numFullSibsM * log(P_FSLA[m]) + 
                          (dailySampledAdultsM - numFullSibsM)
                        * log(1 - P_FSLA[m]))
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
    allRelativeTimes <- seq((- T_L - T_P - 2*T_A), (T_A - T_P), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of adult-larva pair data:")
    for (i in 1:numSampledAdults) {
      # Tally number of full-siblings a given sampled adult has on each day relative to
      # its capture:
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      
      # Calculate likelihood of observed full-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledAdults[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledAdults[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          
          dailySampledLarvaeM <- dailySampledLarvae[allRelativeTimes[m] + 
                                                      sampledAdults[i,]$Time -
                                                      tSamplingStart + 1]
          
          if (dailySampledLarvaeM > 0) {
            logLike <- (logLike + numFullSibsM * log(P_FSAL[m]) + 
                          (dailySampledLarvaeM - numFullSibsM)
                        * log(1 - P_FSAL[m]))
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
    allRelativeTimes <- seq((-T_A - T_L), (T_A + T_L), by=1)
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
      
      # Calculate likelihood of observed full-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledLarvae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledLarvae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          
          dailySampledLarvaeM <- dailySampledLarvae[allRelativeTimes[m] + 
                                                      sampledLarvae[i,]$Time -
                                                      tSamplingStart + 1]
          # Remove entry of self as a sibling:
          relativeTime0 <- which(allRelativeTimes==0)
          if (m == relativeTime0) {
            dailySampledLarvaeM <- dailySampledLarvaeM - 1
          }
          
          if (dailySampledLarvaeM > 0) {
            logLike <- (logLike + numFullSibsM * log(P_FSLL[m]) + 
                          (dailySampledLarvaeM - numFullSibsM)
                        * log(1 - P_FSLL[m]))
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
    allRelativeTimes <- seq((-2*T_A), (2*T_A), by=1)
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
            logLike <- (logLike + numFullSibsM * log(P_FSAA[m]) + 
                          (dailySampledAdultsM - numFullSibsM)
                        * log(1 - P_FSAA[m]))
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
    allRelativeTimes <- seq((- T_A), (T_L + T_P + T_A), by=1)
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
      
      # Calculate likelihood of observed full-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledLarvae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledLarvae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          
          dailySampledPupaeM <- dailySampledPupae[allRelativeTimes[m] + 
                                                    sampledLarvae[i,]$Time -
                                                    tSamplingStart + 1]
          
          if (dailySampledPupaeM > 0) {
            logLike <- (logLike + numFullSibsM * log(P_FSLP[m]) + 
                          (dailySampledPupaeM - numFullSibsM)
                        * log(1 - P_FSLP[m]))
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
    allRelativeTimes <- seq((- T_L - T_P - T_A), (T_A), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of pupa-larva pair data:")
    for (i in 1:numSampledPupae) {
      # Tally number of full-siblings a given sampled pupa has on each day relative to
      # its capture:
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      
      # Calculate likelihood of observed full-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledPupae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledPupae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          
          dailySampledLarvaeM <- dailySampledLarvae[allRelativeTimes[m] + 
                                                      sampledPupae[i,]$Time -
                                                      tSamplingStart + 1]
          
          if (dailySampledLarvaeM > 0) {
            logLike <- (logLike + numFullSibsM * log(P_FSPL[m]) + 
                          (dailySampledLarvaeM - numFullSibsM)
                        * log(1 - P_FSPL[m]))
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
    allRelativeTimes <- seq((-T_A - T_P), (T_A + T_P), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of pupa-pupa pair data:")
    for (i in 1:(numSampledPupae-1)) {
      # Tally number of full-siblings a given sampled larva has on each day relative to
      # its capture:
      fullSibRelativeTimes <- (sampledPupae[(which(sampledPupae[(i+1):numSampledPupae,]$momID
                                                   == sampledPupae[i, "momID"]) + i),]$Time
                               - sampledPupae[i,]$Time)
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      for (k in 1:numRelativeTimes) {
        fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes==allRelativeTimes[k])
      }
      
      # Calculate likelihood of observed full-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledPupae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledPupae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          
          dailySampledPupaeM <- dailySampledPupae[allRelativeTimes[m] + 
                                                    sampledPupae[i,]$Time -
                                                    tSamplingStart + 1]
          # Remove entry of self as a sibling:
          relativeTime0 <- which(allRelativeTimes==0)
          if (m == relativeTime0) {
            dailySampledPupaeM <- dailySampledPupaeM - 1
          }
          
          if (dailySampledPupaeM > 0) {
            logLike <- (logLike + numFullSibsM * log(P_FSPP[m]) + 
                          (dailySampledPupaeM - numFullSibsM)
                        * log(1 - P_FSPP[m]))
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
    allRelativeTimes <- seq((- T_A), (T_P + 2*T_A), by=1)
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
      
      # Calculate likelihood of observed full-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledPupae[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledPupae[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          
          dailySampledAdultsM <- dailySampledAdults[allRelativeTimes[m] + 
                                                      sampledPupae[i,]$Time -
                                                      tSamplingStart + 1]
          
          if (dailySampledAdultsM > 0) {
            logLike <- (logLike + numFullSibsM * log(P_FSPA[m]) + 
                          (dailySampledAdultsM - numFullSibsM)
                        * log(1 - P_FSPA[m]))
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
    allRelativeTimes <- seq((- T_P - 2*T_A), (T_A), by=1)
    numRelativeTimes <- length(allRelativeTimes)
    
    print("Computing log likelihood of adult-pupa pair data:")
    for (i in 1:numSampledAdults) {
      # Tally number of full-siblings a given sampled adult has on each day relative to
      # its capture:
      fullSibRelativeTimesTally <- rep(0,numRelativeTimes)
      
      # Calculate likelihood of observed full-sibling pair data:
      for (m in 1:numRelativeTimes) {
        if (((allRelativeTimes[m] + sampledAdults[i,]$Time) >= tSamplingStart) && 
            ((allRelativeTimes[m] + sampledAdults[i,]$Time) <= tSamplingEnd)) {
          
          numFullSibsM <- fullSibRelativeTimesTally[m]
          
          dailySampledPupaeM <- dailySampledPupae[allRelativeTimes[m] + 
                                                    sampledAdults[i,]$Time -
                                                    tSamplingStart + 1]
          
          if (dailySampledPupaeM > 0) {
            logLike <- (logLike + numFullSibsM * log(P_FSAP[m]) + 
                          (dailySampledPupaeM - numFullSibsM)
                        * log(1 - P_FSAP[m]))
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

# Log likelihood of all full-sibling pair data:

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
