# Close-kin mark-recapture-based parameter inference for mosquitoes

This repository contains kinship probabilities and pseudo-likelihood equations used to infer mosquito demographic parameters based on close-kin sample data. The code accompanies the mosquito close-kin mark-recapture (CKMR) paper available [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010755). Mosquito close-kin data was generated using the [mPlex](https://github.com/GilChrist19/mPlex) simulation framework. A spatial extension of the mosquito CKMR framework has recently been published [here](https://www.biorxiv.org/content/10.1101/2025.05.11.653364v1), with a corresponding linked repository below.

## Repository
  * [Kinship probabilities](./Kinship_probabilities) - Kinship probabilities and pseudo-likelihood equations used to infer mosquito demographic parameters
  * [Simulation code](./Simulation_code) - Code used to generate mosquito close-kin sample data using [mPlex](https://github.com/GilChrist19/mPlex)
  * [Sample data](./Sample_data) - Close-kin sample data generated using [mPlex](https://github.com/GilChrist19/mPlex)
  * [Spatial CKMR](./Spatial_CKMR) Kinship probabilities, simulation code, dispersal kernels and sample data for the spatial extension of the framework

## Authors
[John M. Marshall](https://www.marshalllab.com/) (kinship probabilities), [Jared B. Bennett](https://github.com/GilChrist19) (simulation code)

## Notes and References
 1. [Sharma Y, Bennett JB, Rašić G, Marshall JM (2022) Close-kin mark-recapture methods to estimate demographic parameters of mosquitoes. PLoS Computational Biology 18: e1010755](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010755)
 2. [Bravington MV, Skaug HJ, Anderson EC (2016) Close-kin mark-recapture. Statistical Science 31: 259-274](https://projecteuclid.org/journals/statistical-science/volume-31/issue-2/Close-Kin-Mark-Recapture/10.1214/16-STS552.pdf)
 3. [Bennett JB, Wu SL, Sánchez C. HM, Marshall JM (2019) mPlex: Agent-based modeling of mosquitoes. https://github.com/GilChrist19/mPlex](https://github.com/GilChrist19/mPlex)
 4. [Marshall JM, Yang S, Bennett JB, Filipović I, Rašić G (2025) Spatial close-kin mark-recapture methods to estimate dispersal parameters and barrier strength for mosquitoes. bioRxiv doi: https://doi.org/10.1101/2025.05.11.653364](https://doi.org/10.1101/2025.05.11.653364)