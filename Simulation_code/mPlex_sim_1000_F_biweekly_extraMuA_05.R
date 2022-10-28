###############################################################################
#                            ____  __          ______
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/
#                                               /_/   /_/
###############################################################################
# mPlex simulation for 1000 adult females sampled biweekly over 90 days with
# supplemental lethal sampling equivalent to a mortality rate of 0.05 per adult
# per day.
###############################################################################
# Clean environment and source files
###############################################################################
rm(list=ls());gc()
library(CKMR)

# Directory where the combineFiles.R file is:
# source("~/.../combineFiles.R")

set.seed(seed = 9)
simTime <- 595
numThreads <- 1

###############################################################################
# Setup Directories
###############################################################################

# Directory where the output data goes:
# topDirectory <- "~/.../OUTPUT/..."

if(!dir.exists(paths = topDirectory)){
  dir.create(path = topDirectory)
}

###############################################################################
# Setup Parameters for Network
###############################################################################
# single patch
migration <- as.matrix(x = 1)
numPatch <- nrow(migration)

# batch migration
# set to 0
migrationBatch <- basicBatchMigration(batchProbs = 0, numPatches = numPatch)

#setup alleles to initiate patches
reference <- list('eta'=numeric(0), 'phi'=numeric(0), 'omega'=numeric(0),
                  'xiF'=numeric(0), 'xiM'=numeric(0), 's'=numeric(0))

###############################################################################
# Release Setup
###############################################################################
# Create release list
#  there are no releases, this is just null
patchReleases = replicate(n = numPatch,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      eggReleases = NULL),
                          simplify = FALSE)

###############################################################################
# Sampling Setup
###############################################################################
# Sampling is different for every patch, and every life stage
# final object is a list with 2 matrices in it, one for when to sample, one for
# how many to sample.

# Sampling time:
sampDays <- array(data = FALSE, dim = c(simTime, 5, numPatch))

# Don't sample eggs, only larvae, pupae, adult male, adult female:
sampDays[ ,1, ] <- FALSE

# Setup sampling days:
# Generate a sequence from beginning to end of simulation, then round to whole
# values we will use these as the indices to set sampling days:
fogDayIndex<-round(x=seq.int(from=501, to =simTime, by = 1))  # daily
sampDays[fogDayIndex, 2, ] <- TRUE
sampDays[fogDayIndex, 3, ] <- TRUE
sampDays[fogDayIndex, 4, ] <- TRUE
sampDays[fogDayIndex, 5, ] <- TRUE

# Sampling coverage:
#  Additional lethality for males 5% per day
#  Additional lethality for females 5% per day
sampCov <- array(data = 0, dim = c(simTime, 5, numPatch))
sampCov[fogDayIndex ,2, ] <- 0
sampCov[fogDayIndex ,3, ] <- 0
sampCov[fogDayIndex ,4, ] <- 0.05
sampCov[fogDayIndex ,5, ] <- 0.05

# Sampling for females 9.2% per week, converted to biweekly (to obtain total of
# 1000 sampled adult females):
sampDayIndex<-ceiling(x=seq.int(from=501, to =simTime, by = 3.5)) # biweekly
sampCov[sampDayIndex ,2, ] <- sampCov[sampDayIndex ,2, ] + 0/2
sampCov[sampDayIndex ,3, ] <- sampCov[sampDayIndex ,3, ] + 0/2
sampCov[sampDayIndex ,4, ] <- sampCov[sampDayIndex ,4, ] + 0/2
sampCov[sampDayIndex ,5, ] <- sampCov[sampDayIndex ,5, ] + 0.092/2

# Sample entire adult population at end of simulation (to count population size):
sampDays[simTime,2, ] <- TRUE
sampDays[simTime,3, ] <- TRUE
sampDays[simTime,4, ] <- TRUE
sampDays[simTime,5, ] <- TRUE

sampCov[simTime, 2, ] <- 0
sampCov[simTime, 3, ] <- 0
sampCov[simTime, 4, ] <- 1
sampCov[simTime, 5, ] <- 1

# List to pass to mPlex:
samplingScheme <- list("samplingDays"=sampDays, "samplingCoverage"=sampCov)

###############################################################################
# Setup sweep
###############################################################################
# Setup desired parameter ranges here
#  This then creates a dataframe with every combination

# Sweep over cube changing parameters:
# nRep: indices of reps to do
# nPop: population sizes to test
paramCombo <- as.matrix(expand.grid('nRep' = 1:100,
                                    'nPop' = c(3000) ))
numPC <- NROW(paramCombo)

########################################
# Loop over parameters
########################################
for(i in 1:numPC){

  ####################
  # Setup folder
  ####################
  # width = 7 handles up to a million population size
  simDir <- file.path(topDirectory,
                      formatC(x = paramCombo[i,'nPop'], width = 7, format = "d", flag = "0"),
                      formatC(x = paramCombo[i,'nRep'], width = 3, format = "d", flag = "0"))
  dir.create(path = simDir, recursive = TRUE)

  ####################
  # Network parameters
  ####################
  netPar = NetworkParameters(nPatch = numPatch,
                             simTime = simTime,
                             AdPopEQ = paramCombo[i,'nPop'],
                             runID = 1L,
                             dayGrowthRate = 1.175,
                             beta = 20, tEgg = 2, tLarva = 5, tPupa = 1,
                             muAd = 0.09)

  ####################
  # Run sim!
  ####################
  # Force evaluation of the sampling
  seed <- sample(x = (-.Machine$integer.max):.Machine$integer.max, size = 1, replace = FALSE)
  runCKMR(seed = seed,
          numThreads = numThreads,
          networkParameters = netPar,
          reproductionReference = reference,
          patchReleases = patchReleases,
          migrationMale = migration,
          migrationFemale = migration,
          migrationBatch = migrationBatch,
          samplingParameters = samplingScheme,
          outputDirectory = simDir,
          verbose = FALSE)

  ####################
  # Remove empty files
  ####################
  # combineFiles.R fails when there are empty files, so remove them
  allFiles <- list.files(path = simDir, full.names = TRUE)
  emptyFiles <- which(file.size(allFiles) < 100)
  file.remove(allFiles[emptyFiles])

  ####################
  # Combine files
  ####################
  combineFiles(mainDir=simDir, workIndicator=25, fPattern=c("F","M"))

  ####################
  # Remove initial population
  ####################
  readFiles <- c('/000_F.csv','/000_M.csv')
  writeFiles <- c('/cut_F.csv','/cut_M.csv')

  for(cFile in 1:length(writeFiles)){
    
    # build command,
    cmd <- file.path("-F ',' '(NR==1) || ($1 > 500)' ", simDir, readFiles[cFile],
                     ' > ', simDir, writeFiles[cFile], fsep = '' )
    
    # Pass down to cmdline
    system2(command = 'awk', args = cmd)
    
  } # end loop
  
  ####################
  # Cleanup
  ####################
  # Only need the final cut pops
  # May as well remove the rest
  allFiles <- list.files(path = simDir, full.names = TRUE)
  trashFiles <- grep(pattern = 'cut', x = allFiles, fixed = TRUE, value = TRUE, invert = TRUE)
  file.remove(trashFiles)
  
} # end parameter sweep

detach("package:CKMR", unload=TRUE)