# Dispersal parameters:
lambda <- 1/28.9 # Inverse of dispersal distance conditional upon movement (in meters)
p0 = 0.72 # 

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

# Calculate transition probabilities given dispersal parameters:
M <- matrix(rep(0, numPatches^2), nrow = numPatches, byrow = TRUE)
for (i in 1:numPatches) {
  Denominator <- 0
  for (j in 1:numPatches) {
    if (i!=j) {
      Denominator <- Denominator + (exp(-lambda*distMat[i, j]))
    }
  }
  for (j in 1:numPatches) {
    if (i==j) {
      M[i, j] <- p0
    }
    else {
      M[i, j] <- (1-p0) * (exp(-lambda*distMat[i, j])) / Denominator
    }
  }
}

# Write output:
write.table(x = M, file = "~/Desktop/OUTPUT/mPlex/zero_inf_exp_kernel_19x19.csv",
            sep = ",", row.names = FALSE, col.names = FALSE)