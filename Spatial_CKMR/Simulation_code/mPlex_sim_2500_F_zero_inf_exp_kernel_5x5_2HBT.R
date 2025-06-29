###############################################################################
#                            ____  __          ______
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/
#                                               /_/   /_/
###############################################################################
###############################################################################
#                        _____         _   _____ _ _
#                       |_   _|__  ___| |_|  ___(_) | ___
#                         | |/ _ \/ __| __| |_  | | |/ _ \
#                         | |  __/\__ \ |_|  _| | | |  __/
#                         |_|\___||___/\__|_|   |_|_|\___|
#
###############################################################################
# 20220312
#  This is a copy and update of 20220908.R
#  The sampling methodology was updated for greater flexibility.
#  The updated sampling setup comes from 20220204.R
#
#  19x19 landscape
#  16.6m distance between nodes on the grid
#  exp. rate from MGDrivE
#  Sampling Places:  61,  64,  67,  70,  73,
#                   118, 121, 124, 127, 130,
#                   175, 178, 181, 184, 187,
#                   232, 235, 238, 241, 244,
#                   289, 292, 295, 298, 301
#
###############################################################################
# Clean environment and source files
###############################################################################
rm(list=ls());gc()
library(CKMR)

# Directory where the combineFiles.R file is:
# source("~/.../combineFiles.R")

set.seed(seed = 9)
simTime <- 591
numThreads <- 1 

###############################################################################
# Setup Directories
###############################################################################

# Directory where the output data goes:
# topDirectory <- "~/.../OUTPUT/..."

# if(!dir.exists(paths = topDirectory)){
#   dir.create(path = topDirectory)
# } else {
#   unlink(x = list.files(topDirectory, full.names = TRUE), recursive = TRUE, force = TRUE)
# }

###############################################################################
# Setup Parameters for Network
###############################################################################

# Directory where migration matrix is stored:
# migration <- as.matrix(x = read.csv(file = "/.../OUTPUT/...",header = FALSE))

numPatch <- nrow(migration)

# batch migration
#  set to 0
migrationBatch <- basicBatchMigration(batchProbs = 0, numPatches = numPatch)

#setup alleles to initiate patches
reference <- list('eta'=numeric(0), 'phi'=numeric(0), 'omega'=numeric(0),
                  'xiF'=numeric(0), 'xiM'=numeric(0), 's'=numeric(0))

###############################################################################
# Release Setup
###############################################################################
# create Release List
#  there are no releases, this is just null
patchReleases = replicate(n = numPatch,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      eggReleases = NULL),
                          simplify = FALSE)

###############################################################################
# Sampling Setup
###############################################################################
# sampling is different for every patch, and every life stage
# final object is a list with 2 arrays in it, one for when to sample, one for 
#  how many to sample

# Sampling Time
#  Setup default object - no sampling of any stage on any day
sampDays <- array(data = FALSE, dim = c(simTime, 5, numPatch))

#  We will sample the adult stages (stages 4 and 5)
sampStages <- c(4, 5)

#  We only sample specific patches in the landscape ("sampPlaces" variable)
sampPlaces <- c( 61,  64,  67,  70,  73,
                 118, 121, 124, 127, 130,
                 175, 178, 181, 184, 187,
                 232, 235, 238, 241, 244,
                 289, 292, 295, 298, 301)

#  We sample the stages and places biweekly
sdIndex<-ceiling(x=seq.int(from=501, to =simTime, by = 3.5)) # biweekly
sampDays[sdIndex , sampStages, sampPlaces] <- TRUE

# Sampling Coverage
#  Setup default object - 0 sampling rate
sampCov <- array(data = 0, dim = c(simTime, 5, numPatch))

sampCov[ ,1, ] <- 0
sampCov[ , 2, sampPlaces] <- 0
sampCov[ ,3, ] <- 0
sampCov[ ,4, ] <- 0

#  We sample adult females at a rate required to achieve the desired sample size:
sampCov[ , 5, sampPlaces] <- (2.5/2) * (36/25) * 1.5 * (4/3) * 0.19 / 2

# Sample entire adult population at end of simulation:
sampDays[simTime,2, ] <- TRUE
sampDays[simTime,3, ] <- TRUE
sampDays[simTime,4, ] <- TRUE
sampDays[simTime,5, ] <- TRUE

sampCov[simTime, 2, ] <- 0
sampCov[simTime, 3, ] <- 0
sampCov[simTime, 4, ] <- 1
sampCov[simTime, 5, ] <- 1

# list to pass to mPlex
samplingScheme <- list("samplingDays"=sampDays, "samplingCoverage"=sampCov)

###############################################################################
# Setup Sweep
###############################################################################
# Setup desired parameter ranges here
#  This then creates a dataframe with every combination

# sweep over cube changing parameters
# nRep: indices of reps to do
# nPop: population sizes to test
paramCombo <- as.matrix(expand.grid('nRep' = 1:125,
                                    'nPop' = c(25) ))
numPC <- NROW(paramCombo)

########################################
# Loop over parameters
########################################
for(i in 1:numPC){
  
  ####################
  # Setup Folder
  ####################
  # width = 7 handles up to a million population size
  simDir <- file.path(topDirectory,
                      formatC(x = paramCombo[i,'nPop'], width = 7, format = "d", flag = "0"),
                      formatC(x = paramCombo[i,'nRep'], width = 3, format = "d", flag = "0"))
  dir.create(path = simDir, recursive = TRUE)
  
  ####################
  # Network Parameters
  ####################
  netPar = NetworkParameters(nPatch = numPatch,
                             simTime = simTime,
                             AdPopEQ = rep.int(x = paramCombo[i,'nPop'], times = numPatch),
                             runID = 1L,
                             dayGrowthRate = 1.175,
                             beta = 20, tEgg = 2, tLarva = 5, tPupa = 1,
                             muAd = 0.09)
  
  ####################
  # Run Sim!
  ####################
  # force evaluation of the sampling
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
  # remove empty files
  ####################
  # combineFiles.R fails when there are empty files, so remove them
  allFiles <- list.files(path = simDir, full.names = TRUE)
  emptyFiles <- which(file.size(allFiles) < 100)
  file.remove(allFiles[emptyFiles])
  
  ####################
  # Combine Files
  ####################
  # no output
  combineFiles(mainDir=simDir, workIndicator=25, fPattern=c("F","M"))
  
  ####################
  # Remove initial population
  ####################
  # Remove any days with "0" parents
  # this is to get rid of "0" parents, the ones who started the population so there was 
  #  no structure before that.
  # I will check male and female files, then remove both.
  # Or, screw it, remove first 100 days, set maxVal = 100.
  #  
  # 
  # -F: set column delimiter
  # $5/$6: column numbers
  # print if they equal 0
  # NR: number record, it is the line number, so keep header
  # then print first column from those columns.
  # then, bye eye, get max value
  # 
  # This one takes off anyone with 0 parents
  # awk -F "," '($5 == "0") && ($6 == "0")' 000_F.csv | awk -F "," '{print $1}'
  # 
  # This one takes off first 100 days
  # awk -F "," '(NR==1) || ($1 > 100)' 000_F.csv > newFemFile.csv
  
  readFiles <- c('/000_F.csv', '/000_M.csv')
  writeFiles <- c('/cut_F.csv', '/cut_M.csv')
  for(cFile in 1:length(writeFiles)){
    
    # build command, 
    cmd <- file.path("-F ',' '(NR==1) || ($1 > 100)' ", simDir, readFiles[cFile], 
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
