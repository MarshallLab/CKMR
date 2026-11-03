# Parent-offspring and full-sibling negative log-likelihoods, followed by their sum.
# These functions read prepared sampling objects from their enclosing environment;
# run_example.R sources all stages into the same private analysis environment.
# Each likelihood constructs the daily transition probabilities it needs.

# 1. Mother-offspring (adult) negative log-likelihood -----------------------------
logLike_MOA <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A,
                        beta, mu_E, mu_L, mu_P,
                        lambda, lambda_ver, barrier_strength, numPatches, distMat, kin_df) {
  # Calculate transition probabilities given dispersal parameters ------------
  # First, calculate the daily transition matrix, which describes the probability that,
  # given a mosquito is at node i on one day, it is at node j on the next.
  # Normalize exponential weights with horizontal/vertical rates and a barrier term.

  # Grid-index barrier: cells <= 98 versus cells > 98 within each layer.
  barrier_node <- 98

  M <- matrix(rep(0, numPatches^2), nrow = numPatches, byrow = TRUE)
  for (i in 1:numPatches) {
    Denominator <- 0

    floor_i <- (i - 1) %/% patches_per_floor
    patch_i <- (i - 1) %% patches_per_floor + 1

    for (j in 1:numPatches) {
      floor_j <- (j - 1) %/% patches_per_floor
      patch_j <- (j - 1) %% patches_per_floor + 1

      if (floor_i == floor_j) {
        if (((patch_i <= barrier_node) && (patch_j <= barrier_node)) ||
          ((patch_i > barrier_node) && (patch_j > barrier_node))) {
          Denominator <- Denominator + exp(-lambda * distMat[i, j])
        } else {
          Denominator <- Denominator + (1 - barrier_strength) * exp(-lambda * distMat[i, j])
        }
      } else {
        floor_diff <- abs(floor_i - floor_j)
        total_dist <- distMat[i, j]
        horiz_dist <- distMat[patch_i, patch_j]
        ver_dist <- sqrt(total_dist^2 - horiz_dist^2)

        if (((patch_i <= barrier_node) && (patch_j <= barrier_node)) ||
          ((patch_i > barrier_node) && (patch_j > barrier_node))) {
          Denominator <- Denominator + exp(-lambda * horiz_dist - lambda_ver * ver_dist)
        } else {
          Denominator <- Denominator + (1 - barrier_strength) * exp(-lambda * horiz_dist - lambda_ver * ver_dist)
        }
      }
    }

    for (j in 1:numPatches) {
      floor_j <- (j - 1) %/% patches_per_floor
      patch_j <- (j - 1) %% patches_per_floor + 1

      floor_diff <- abs(floor_i - floor_j)
      total_dist <- distMat[i, j]
      horiz_dist <- distMat[patch_i, patch_j]
      ver_dist <- sqrt(total_dist^2 - horiz_dist^2)

      if (floor_i == floor_j) {
        if (((patch_i <= barrier_node) && (patch_j <= barrier_node)) ||
          ((patch_i > barrier_node) && (patch_j > barrier_node))) {
          M[i, j] <- exp(-lambda * distMat[i, j]) / Denominator
        } else {
          M[i, j] <- (1 - barrier_strength) * exp(-lambda * distMat[i, j]) / Denominator
        }
      } else {
        if (((patch_i <= barrier_node) && (patch_j <= barrier_node)) ||
          ((patch_i > barrier_node) && (patch_j > barrier_node))) {
          M[i, j] <- exp(-lambda * horiz_dist - lambda_ver * ver_dist) / Denominator
        } else {
          M[i, j] <- (1 - barrier_strength) * exp(-lambda * horiz_dist - lambda_ver * ver_dist) / Denominator
        }
      }
    }
  }

  # Next, calculate the transition matrices for 1, ..., T_A days, which describe the
  # probability that, given a mosquito is at node i on one day, it is at node j 1, ...,
  # T_A days later. The first slice is M (one day), followed by its powers:

  # Initialize a 3D array rho: [origin, destination, time lag]
  rho <- replicate(T_A, M)
  for (i in 2:T_A) {
    # Multiply matrices recursively to get transition probabilities over multiple days
    rho[, , i] <- rho[, , (i - 1)] %*% M
  }

  # Probability of adult surviving from 0 to (T_A-1) days:
  AdultSurvivalProbability <- rep(0, T_A)
  AdultAge <- rep(0, T_A)
  for (i in 1:T_A) {
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i - 1)
    AdultAge[i] <- i - 1
  }

  # Probability of adult having an age between 0 and (T_A-1) days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)

  # Mother-offspring (adult) kinship probabilities ------------
  # Given an adult sampled at time t2, this is the probability that an adult female
  # sampled at time t1 is their mother.

  # First, calculate the denominator:
  # This is the expected number of surviving adults at time t2 from adult females at any
  # consistent time (assuming a constant population size, this is independent of time).

  Denominator <- 0
  t2 <- 0 # The denominator should be the same for any t2
  for (y2 in (t2 - T_E - T_L - T_P - (T_A - 1)):(t2 - T_E - T_L - T_P)) {
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
  for (i in 1:T_A) {
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i - 1)
    AdultAge[i] <- i - 1
  }

  # Probability of adult having an age between 0 and (T_A-1) days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)

  t1 <- 0 # The relative difference between t1 & t2 is what matters
  t2 <- seq((-(T_A - 1) + T_E + T_L + T_P), (T_E + T_L + T_P + (T_A - 1)), by = 1)

  Numerator <- replicate(length(t2), matrix(rep(0, numSampledPatches^2),
    nrow = numSampledPatches, byrow = TRUE
  ))

  fixed_surv <- beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) * ((1 - mu_P)^T_P)
  # Loop over all pairs of sampled patches to calculate the expected number
  # of individuals from the two locations are of MOA kinship
  for (x1 in 1:numSampledPatches) {
    for (x2 in 1:numSampledPatches) {
      for (i in 1:length(t2)) {
        t2i <- t2[i]
        for (y2 in (t2i - T_E - T_L - T_P - (T_A - 1)):(t2i - T_E - T_L - T_P)) {
          if ((y2 >= (t1 - (T_A - 1))) && (y2 <= t1)) {
            # First calculate the probability that, given a mother sampled at location
            # x1, the adult offspring is sampled at location x2 (rho_x1_x2):
            rho_xi_x1 <- rho[, sampledPatches[x1], (1 + t1 - y2)]
            sum_rho_xi_x1 <- sum(rho_xi_x1)
            rho_xi_x2 <- rho[, sampledPatches[x2], (1 + t2i - (y2 + T_E + T_L + T_P))]
            rho_x1_x2 <- as.numeric((rho_xi_x1 / sum_rho_xi_x1) %*% rho_xi_x2)

            # Then calculate the numerator:
            Numerator[x1, x2, i] <- Numerator[x1, x2, i] + ((1 - mu_A)^(t1 - y2)
              * rho_x1_x2 * fixed_surv
              * ((1 - mu_A)^(t2[i] - y2 - T_E - T_L - T_P)))
          }
        }
      }
    }
  }

  # Mother-adult offspring probability:
  # Given an adult sampled at time t2 & location x2, the probability that an adult
  # female sampled at time t1 & location x1 is their mother is given by:

  P_MOA <- Numerator / Denominator

  # Note that indices here relate to times in the vector, t2, where t1 = 0, and hence,
  # in general, the indices relate to times t2-t1, i.e.:

  t2Minust1 <- t2

  # Find mother-offspring (adult) pairs in the data ------------
  is_D1 <- kin_df$PO_direction == "D1"

  kin_df <- kin_df %>%
    mutate(
      motherID = ifelse(is_D1, D1_indiv, D2_indiv),
      adultID = ifelse(is_D1, D2_indiv, D1_indiv),
      motherSamplingTime = ifelse(is_D1, D1_date, D2_date),
      adultSamplingTime = ifelse(is_D1, D2_date, D1_date),
      motherSamplingPatch = ifelse(is_D1, D1_Trap_ID, D2_Trap_ID),
      adultSamplingPatch = ifelse(is_D1, D2_Trap_ID, D1_Trap_ID)
    )

  # Array of mother-adult offspring pairs keeping track of: i) adult ID, ii) mother ID,
  # iii) day adult offspring sampled, iv) patch adult offspring sampled from, v) day
  # mother sampled & vi) patch mother sampled from:
  MOA_Pairs_Data <- kin_df[, c(
    "adultID", "motherID", "adultSamplingTime", "adultSamplingPatch",
    "motherSamplingTime", "motherSamplingPatch"
  )]

  # Record number of adult females sampled on each day at each sampling patch:
  dailySampledAdultFemales <- matrix(rep(0, numSampledPatches * numSamplingDays),
    nrow = numSampledPatches, byrow = TRUE
  )

  each <- c()
  for (x1 in 1:numSampledPatches) {
    for (i in 1:numSamplingDays) {
      dailySampledAdultFemales[x1, i] <- sum(sampledAdultFemales[
        which(sampledAdultFemales$patch_id == sampledPatches[x1]),
      ]$Date_of_collection == samplingDays[i])
    }
  }

  # Record number of adults sampled on each day at each sampling patch:
  dailySampledAdults <- matrix(rep(0, numSampledPatches * numSamplingDays),
    nrow = numSampledPatches, byrow = TRUE
  )

  # Count number of adults sampled at patch x1 on day i
  for (x1 in 1:numSampledPatches) {
    for (i in 1:numSamplingDays) {
      dailySampledAdults[x1, i] <- sum(sampledAdults[
        which(sampledAdults$patch_id == sampledPatches[x1]),
      ]$Date_of_collection == samplingDays[i])
    }
  }

  # Calculate likelihood of mother-offspring (adult) pair data ------------
  logLike <- 0

  for (x1 in 1:numSampledPatches) {
    for (x2 in 1:numSampledPatches) {
      # Calculate the loglikelihood that, for x sampled adults on day t2 at location x2,
      # w have mothers among the y sampled adult females on day t1 at location x1:

      for (i in 1:length(samplingDays)) {
        adultSamplingTimeI <- samplingDays[i]

        if (!(adultSamplingTimeI %in% sampledAdults$Date_of_collection)) {
          next
        }

        x <- dailySampledAdults[x2, i] # Number of sampled adults on day t2 at location x2

        for (j in 1:length(t2Minust1)) {
          motherSamplingTimeJ <- adultSamplingTimeI - t2Minust1[j]

          if ((motherSamplingTimeJ >= tSamplingStart) & (motherSamplingTimeJ <= tSamplingEnd)) {
            # Number of sampled adult females on day t1 at location x1:
            y <- dailySampledAdultFemales[x1, which(samplingDays == motherSamplingTimeJ)]
            if (length(y) == 0) {
              next
            }

            if ((j == which(t2Minust1 == 0)) && (x1 == x2)) {
              y <- (y - 1)
            }

            if ((y > 0) & (P_MOA[x1, x2, j] > 0)) {
              # Number of sampled adults at location x2 on day t2 that have mothers among
              # the sampled adult females at x1 on day t1:
              w <- sum(MOA_Pairs_Data[, "motherSamplingTime"]
              [((MOA_Pairs_Data[, "adultSamplingTime"] == adultSamplingTimeI) &
                  ((MOA_Pairs_Data[, "adultSamplingPatch"] == sampledPatches[x2]) &
                    (MOA_Pairs_Data[, "motherSamplingPatch"] == sampledPatches[x1])))] == motherSamplingTimeJ)
              # Probability that a given sampled adult on day t2 at location x2 has a
              # mother among the y sampled adult females on day t1 at location x1:
              z <- 1 - ((1 - P_MOA[x1, x2, j])^y)

              # Now calculate the log likelihood that w sampled adults on day t2 at
              # location x2 have mothers among the sampled adult females at location x1
              # on day t1:
              logLike <- logLike + (w * log(z)) + ((x - w) * log(1 - z))
            }
          }
        }
      }
    }
  }

  -logLike
}

# 2. Full-sibling (adult-adult) negative log-likelihood -------------------------
logLike_Sibs <- function(N_A, mu_A, N_F, T_E, T_L, T_P, T_A,
                         beta, mu_E, mu_L, mu_P,
                         lambda, lambda_ver, barrier_strength, numPatches, distMat, kin_df) {
  # Calculate transition probabilities given dispersal parameters ------------
  # First, calculate the daily transition matrix, which describes the probability that,
  # given a mosquito is at node i on one day, it is at node j on the next. Assuming a
  # normalized exponential weights, with separate horizontal and vertical rates
  # and a multiplier of (1 - barrier_strength) across the barrier:

  # Grid-index barrier: cells <= 98 versus cells > 98 within each layer.
  barrier_node <- 98

  M <- matrix(rep(0, numPatches^2), nrow = numPatches, byrow = TRUE)
  for (i in 1:numPatches) {
    Denominator <- 0

    floor_i <- (i - 1) %/% patches_per_floor
    patch_i <- (i - 1) %% patches_per_floor + 1

    for (j in 1:numPatches) {
      floor_j <- (j - 1) %/% patches_per_floor
      patch_j <- (j - 1) %% patches_per_floor + 1

      if (floor_i == floor_j) {
        if (((patch_i <= barrier_node) && (patch_j <= barrier_node)) ||
          ((patch_i > barrier_node) && (patch_j > barrier_node))) {
          Denominator <- Denominator + exp(-lambda * distMat[i, j])
        } else {
          Denominator <- Denominator + (1 - barrier_strength) * exp(-lambda * distMat[i, j])
        }
      } else {
        floor_diff <- abs(floor_i - floor_j)
        total_dist <- distMat[i, j]
        horiz_dist <- distMat[patch_i, patch_j]
        ver_dist <- sqrt(total_dist^2 - horiz_dist^2)

        if (((patch_i <= barrier_node) && (patch_j <= barrier_node)) ||
          ((patch_i > barrier_node) && (patch_j > barrier_node))) {
          Denominator <- Denominator + exp(-lambda * horiz_dist - lambda_ver * ver_dist)
        } else {
          Denominator <- Denominator + (1 - barrier_strength) * exp(-lambda * horiz_dist - lambda_ver * ver_dist)
        }
      }
    }

    for (j in 1:numPatches) {
      floor_j <- (j - 1) %/% patches_per_floor
      patch_j <- (j - 1) %% patches_per_floor + 1

      floor_diff <- abs(floor_i - floor_j)
      total_dist <- distMat[i, j]
      horiz_dist <- distMat[patch_i, patch_j]
      ver_dist <- sqrt(total_dist^2 - horiz_dist^2)

      if (floor_i == floor_j) {
        if (((patch_i <= barrier_node) && (patch_j <= barrier_node)) ||
          ((patch_i > barrier_node) && (patch_j > barrier_node))) {
          M[i, j] <- exp(-lambda * distMat[i, j]) / Denominator
        } else {
          M[i, j] <- (1 - barrier_strength) * exp(-lambda * distMat[i, j]) / Denominator
        }
      } else {
        if (((patch_i <= barrier_node) && (patch_j <= barrier_node)) ||
          ((patch_i > barrier_node) && (patch_j > barrier_node))) {
          M[i, j] <- exp(-lambda * horiz_dist - lambda_ver * ver_dist) / Denominator
        } else {
          M[i, j] <- (1 - barrier_strength) * exp(-lambda * horiz_dist - lambda_ver * ver_dist) / Denominator
        }
      }
    }
  }

  # Next, calculate the transition matrices for 1, ..., T_A days or (2*T_A) days ((2*T_A)
  # days is for the case where adults are sampled, to account for the possibility that the
  # sampled adults may move in addition to their mother). This describes the probability
  # that, given a mosquito is at node i on one day, it is at node j 1, ..., T_A days
  # later. The first slice is M (one day), followed by its powers:

  if (numSampledAdults == 0) {
    rho <- replicate(T_A, M)
    for (i in 2:T_A) {
      rho[, , i] <- rho[, , (i - 1)] %*% M
    }
  }
  if (numSampledAdults > 0) {
    rho <- replicate((2 * T_A), M)
    for (i in 2:(2 * T_A)) {
      rho[, , i] <- rho[, , (i - 1)] %*% M
    }
  }

  # Full-sibling (adult-adult) kinship probabilities ------------
  # Probability of adult surviving from 0 to (T_A-1) days:
  AdultSurvivalProbability <- rep(0, T_A)
  AdultAge <- rep(0, T_A)
  for (i in 1:T_A) {
    AdultSurvivalProbability[i] <- (1 - mu_A)^(i - 1)
    AdultAge[i] <- i - 1
  }

  # Probability of adult having an age between 0 and (T_A-1) days:
  AdultAgeProbability <- AdultSurvivalProbability / sum(AdultSurvivalProbability)
  # Given an adult sampled at time t2, this is the probability that an adult sampled at
  # time t1 is its full sibling.

  if (numSampledAdults > 0) {
    # First, calculate the denominator:
    # This is the expected number of surviving adults at time t2 from adult females at any
    # consistent time (assuming a constant population size, this is independent of time).

    Denominator <- 0
    t2 <- 0 # The denominator should be the same for any t2
    for (y2 in (t2 - T_E - T_L - T_P - (T_A - 1)):(t2 - T_E - T_L - T_P)) {
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

    t1 <- 0 # The relative difference between t1 & t2 is what matters
    t2 <- seq((-(T_A - 1) - (T_A - 1)), ((T_A - 1) + (T_A - 1)), by = 1)

    Numerator <- replicate(length(t2), matrix(rep(0, numSampledPatches^2),
      nrow = numSampledPatches, byrow = TRUE
    ))

    print("Calculating full-sibling adult-adult probabilities:")

    i <- 1:length(t2)
    y1 <- (t1 - T_E - T_L - T_P - (T_A - 1)):(t1 - T_E - T_L - T_P)
    data <- as.data.frame(cbind(i = rep(i, each = length(y1)), y1 = rep(y1, length(i))))

    data$y2 <- lapply(data$y1, function(input) {
      return((input - (T_A - 1)):(input + (T_A - 1)))
    })
    data <- data %>% unnest(y2)
    data$t2i <- t2[data$i]
    data <- data %>% filter((y2 >= (t2i - T_E - T_L - T_P - (T_A - 1))) &
      (y2 <= (t2i - T_E - T_L - T_P)))

    data$Numerator <- mapply(data$y1, FUN = function(input) {
      return(which(AdultAge == (t1 - input - T_E - T_L - T_P)))
    })
    data$NumeratorTemp <- AdultAgeProbability[data$Numerator] * (1 / 2) * (1 - mu_A)^(abs(data$y2 - data$y1)) * beta * ((1 - mu_E)^T_E) * ((1 - mu_L)^T_L) * ((1 - mu_P)^T_P) * ((1 - mu_A)^(data$t2i - data$y2 - T_E - T_L - T_P))

    x1 <- 1:numSampledPatches
    x1 <- rep(x1, nrow(data))
    x2 <- 1:numSampledPatches
    x2 <- rep(x2, nrow(data))

    data <- data[rep(seq_len(nrow(data)), each = numSampledPatches), ]
    data <- cbind(data, x1)
    data_x2 <- cbind(data, x2)

    data$samp_x1 <- sampledPatches[data$x1]
    data$samp_col <- (data$y1 - data$y2) * (data$y2 < data$y1) + (1 + t1 - (data$y1 + T_E + T_L + T_P))
    data_x2$samp_x2 <- sampledPatches[data_x2$x2]
    data_x2$samp_col_x2 <- (data_x2$y2 - data_x2$y1) * (data_x2$y2 > data_x2$y1) + (1 + data_x2$t2i - (data_x2$y2 + T_E + T_L + T_P))

    # First calculate the probability that, given an adult sampled at location
    # x1, an adult sampled at location x2 is its sibling (rho_x1_x2):
    for (u in unique(data$i)) {
      slice_x1 <- data %>% filter(i == u)
      slice_x2 <- data_x2 %>% filter(i == u)
      ma_x1 <- matrix(0, nrow = nrow(slice_x1), ncol = numPatches)
      for (k in 1:nrow(slice_x1)) {
        x1 <- slice_x1$x1[k]
        i <- slice_x1$i[k]
        temp_x1 <- rho[, slice_x1$samp_x1[k], slice_x1$samp_col[k]]
        temp_x1 <- temp_x1 / sum(temp_x1)
        ma_x1[k, ] <- temp_x1
      }
      ma_x2 <- matrix(0, nrow = nrow(slice_x2), ncol = numPatches)
      for (j in 1:nrow(slice_x2)) {
        temp_x2 <- rho[, slice_x2$samp_x2[j], slice_x2$samp_col_x2[j]] * slice_x2$NumeratorTemp[j]
        ma_x2[j, ] <- temp_x2
      }

      rho_x1_x2 <- array(0, dim = c(numSampledPatches, ncol = numSampledPatches))
      rep_times <- sum(slice_x1$x1 == 1)

      for (rep in 1:rep_times) {
        partial_x1 <- ma_x1[(1 + numSampledPatches * (rep - 1)):(numSampledPatches * rep), ]
        partial_x2 <- ma_x2[(1 + numSampledPatches * (rep - 1)):(numSampledPatches * rep), ]
        rho_x1_x2 <- rho_x1_x2 + as.numeric(partial_x1 %*% t(partial_x2))
      }

      Numerator[, , u] <- rho_x1_x2
    }

    # Full sibling adult-adult probability:
    # Given an adult sampled at time t2, the probability that an adult sampled at time t1 is
    # its full sibling is given by:

    P_FSAA <- Numerator / Denominator
    rownames(P_FSAA) <- sampledPatches
    colnames(P_FSAA) <- sampledPatches

    # Note that indices here relate to times in the vector, t2, where t1 = 0, and
    # hence, in general, the indices relate to times t2-t1, i.e.:

    t2Minust1_FSAA <- seq((-2 * (T_A - 1)), (2 * (T_A - 1)), by = 1)
  }

  # Calculate the likelihood of full-sibling pair data ------------
  # Record number of adults sampled on each day at each sampling patch:
  if (numSampledAdults > 0) {
    dailySampledAdults <- matrix(rep(0, numSampledPatches * numSamplingDays),
      nrow = numSampledPatches, byrow = TRUE
    )

    for (x1 in 1:numSampledPatches) {
      for (i in 1:numSamplingDays) {
        dailySampledAdults[x1, i] <- sum(sampledAdults[
          which(sampledAdults$patch_id == sampledPatches[x1]),
        ]$Date_of_collection == samplingDays[i])
      }
    }
    rownames(dailySampledAdults) <- sampledPatches
  }

  # Initialize log likelihood:
  logLike <- 0

  # Calculate the likelihood of adult-adult sibling pair data ------------
  if (numSampledAdults > 0) {
    print("Computing log likelihood of adult-adult pair data:")
    allRelativeTimes <- seq((-2 * (T_A - 1)), (2 * (T_A - 1)), by = 1)
    numRelativeTimes <- length(allRelativeTimes)
    allRelativeTimes <- seq((-2 * (T_A - 1)), (2 * (T_A - 1)), by = 1)
    numRelativeTimes <- length(allRelativeTimes)
    i <- 1:(numSampledAdults - 1)
    #  Make a data frame with ID and list of full siblings
    data <- data.frame(ID = sampledAdults$ID, adultTime = sampledAdults$Date_of_collection, x1 = sampledAdults$patch_id)
    # Step 1: Build the undirected graph of FS pairs
    fs_edges <- kin_df[, c("D1_indiv", "D2_indiv")]
    fs_edges <- fs_edges[complete.cases(fs_edges), ]
    g <- graph_from_data_frame(fs_edges, directed = FALSE)

    # Step 2: Get connected components (each is a full-sibling family)
    comp <- components(g)

    # Step 3: Build a lookup table: individual → family index
    family_index_lookup <- comp$membership # Named vector: names = IDs, values = family indices

    # Assign a full-sibling family to each sampled adult; unmatched IDs remain NA.
    sampledAdults$FS_family_index <- family_index_lookup[as.character(sampledAdults$ID)]

    data$fullSibs <- mapply(FUN = function(a) {
      return(which(sampledAdults[(a + 1):numSampledAdults, ]$FS_family_index == sampledAdults[a, "FS_family_index"]))
    }, 1:nrow(data))

    # Tally full siblings of a sampled adult at location x2 on each
    # day relative to its capture:
    x2 <- sort(unique(sampledAdults$patch_id))
    x2 <- rep(x2, nrow(data))
    data <- data[rep(seq_len(nrow(data)), each = numSampledPatches), ]
    data <- cbind(data, x2)

    data$patchX2 <- mapply(FUN = function(a, b) {
      which(sampledAdults[(which(sampledAdults$ID == a) + 1):numSampledAdults, ]$patch_id == b)
    }, data$ID, data$x2)
    calc_intersect <- function(a, b, c, d) {
      if (length(a) == 0) {
        return(c())
      } else {
        union <- intersect(a, b)
        fullSibRelativeTimes <- (sampledAdults[(union + which(sampledAdults$ID == c)), ]$Date_of_collection - d)
        fullSibRelativeTimesTally <- rep(0, numRelativeTimes)
        for (k in 1:numRelativeTimes) {
          fullSibRelativeTimesTally[k] <- sum(fullSibRelativeTimes == allRelativeTimes[k])
        }
        return(fullSibRelativeTimesTally)
      }
    }
    data$fullSibRelativeTimesTally <- mapply(FUN = calc_intersect, data$fullSibs, data$patchX2, data$ID, data$adultTime)
    # Calculate likelihood of observed full-sibling pair data:
    m <- 1:numRelativeTimes
    m <- rep(m, nrow(data))
    data <- data[rep(seq_len(nrow(data)), each = numRelativeTimes), ]
    data <- cbind(data, m)
    data$relativeTime <- allRelativeTimes[data$m]
    data <- data %>% filter((relativeTime + adultTime >= tSamplingStart) & (relativeTime + adultTime <= tSamplingEnd))
    data$numFullSibs <- mapply(FUN = calc_numFullSib, data$fullSibRelativeTimesTally, data$m)
    data$dailySampledAdultsM <- mapply(
      FUN = function(a, b, c) {
        return(dailySampledAdults[which(sampledPatches == a), (b + c - tSamplingStart + 1)])
      },
      data$x2, data$relativeTime, data$adultTime
    )
    relativeTime0 <- which(allRelativeTimes == 0)
    data$dailySampledAdultsM <- data$dailySampledAdultsM - ((data$m == relativeTime0) & (data$x1 == data$x2))
    data <- data %>% filter(dailySampledAdultsM > 0)
    data$P_FSAA_value <- mapply(FUN = function(a, b, c) {
      return(P_FSAA[which(sampledPatches == a), which(sampledPatches == b), c])
    }, data$x1, data$x2, data$m)
    data$result <- data$numFullSibs * log(data$P_FSAA_value) + (data$dailySampledAdultsM - data$numFullSibs) * log(1 - data$P_FSAA_value)
    logLike <- logLike + sum(data$result)
  }

  -logLike
}

# 3. Combined objective and fixed life-history parameters ----------------------
logLike_all <- function(DispersalPars) {
  # Extract parameters that we are varying from the DispersalPars vector:
  lambda <- DispersalPars[1]

  # Parameters that we are setting constant in the life history model (comment
  # out the ones you are trying to estimate):
  numPatches <- nrow(distMat)
  # The likelihood uses the fixed per-node population.
  # Total population is calculated separately for the output summary.
  N_A <- node_population_size
  lambda_ver <- DispersalPars[2] # Vertical rate applied to metric distances in distMat
  barrier_strength <- DispersalPars[3]
  mu_A <- 0.09 # Daily mortality of adult mosquitoes
  N_F <- N_A / 2 # Adult female population per node

  T_E <- 2 # Duration of the egg stage (days)
  T_L <- 5 # Duration of the larval stage (days)
  T_P <- 1 # Duration of the pupal stage (days)
  T_A <- 40 # Maximum adult lifespan considered (days)

  rM <- 1.175 # Daily mosquito population growth rate
  beta <- 20 # Daily egg production rate
  G <- T_E + T_L + T_P + (1 / mu_A) # Mosquito generation time
  R_M <- rM^G # Population growth rate per generation

  # Daily mortality of density-independent juvenile life stages:
  mu_J <- 1 - ((2 * R_M * mu_A) / (beta * (1 - mu_A)))^(1 / (T_E + T_L + T_P))

  mu_E <- mu_J # Daily mortality of egg stage
  mu_P <- mu_J # Daily mortality of pupal stage

  # Daily mortality of density-dependent larval life stage at equilibrium:
  mu_L <- 1 - (2 * mu_A / (beta * (1 - mu_A) * ((1 - mu_J)^(T_E + T_P))))^(1 / T_L)

  sampledPatchesAdultFemales <- sort(unique(sampledAdultFemales$patch_id))
  numSampledPatchesAdultFemales <- length(sampledPatchesAdultFemales)
  numSampledPatches <- numSampledPatchesAdultFemales

  logLike <- 0 # Initialize the log-likelihood

  if (numSampledAdultFemales > 0) {
    logLike <- logLike + logLike_MOA(
      N_A, mu_A, N_F, T_E, T_L, T_P, T_A,
      beta, mu_E, mu_L, mu_P,
      lambda, lambda_ver, barrier_strength, numPatches, distMat, kinship_PO_barrier
    )
  }
  logLike <- logLike + logLike_Sibs(
    N_A, mu_A, N_F, T_E, T_L, T_P, T_A,
    beta, mu_E, mu_L, mu_P,
    lambda, lambda_ver, barrier_strength, numPatches, distMat, kinship_FS_barrier
  )

  print('N_A_fixed, lambda, lambda_ver, barrier_strength, logLike:')
  print(c(N_A, lambda, lambda_ver, barrier_strength, logLike))
  logLike
}
