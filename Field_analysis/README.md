# CKMR estimation of dispersal and barrier effects

This example estimates horizontal dispersal rate (`lambda`), vertical dispersal rate (`lambda_ver`), and barrier strength from parent–offspring and full-sibling pairs in field-collected mosquito data. The adult population is fixed at five individuals per movement node.

The analysis uses season-2 records from blocks `943:953`, `835`, and `837:850`. Observed floors 2, 7, and 12 correspond to sampled layers 1, 3, and 5 in a five-layer movement model.

| Observed floor | Model layer | Sampled adults |
| --- | --- | ---: |
| 2 | 1 | 408 |
| 7 | 3 | 93 |
| 12 | 5 | 39 |

With the study input files, preparation yields 540 sampled adults in 42 sampled patches, with 49 parent–offspring pair records and 45 full-sibling pair records after filtering. The movement matrix has dimensions **910 × 910**, representing 182 horizontal grid cells across five layers. Layers 2 and 4 are unsampled intermediate movement layers.

## Requirements and execution

The required R packages are `ggplot2`, `sf`, `dplyr`, `igraph`, `tidyr`, `optimx`, and `stringr`. From this folder, install missing packages once:

```sh
Rscript install_dependencies.R
```

Installation of `sf` may also require system libraries appropriate to your operating system and R installation.

Input data are not distributed with the code. Before running, create a local `data/` directory and place `season2_filtered_v2.csv`, `Kinship_PO_v3.csv`, and `Kinship_FS_v3.csv` in it. File formats are described in [Input data](docs/DATA.md). The `data/` directory is excluded from Git by `.gitignore`.

To prepare and inspect the data without fitting the model:

```sh
Rscript run_example.R --prepare-only
```

To run the complete analysis:

```sh
Rscript run_example.R
```

Paths are resolved relative to the entry-point script. The script can also be called from another working directory using its full path. Fitting uses dense matrix operations and nested loops; the optimizer and subsequent Hessian calculation each require repeated likelihood evaluations.

## Code structure

| File | Purpose |
| --- | --- |
| `run_example.R` | Entry point; loads dependencies and runs the analysis stages |
| `config.R` | Paths, sample selection, floor mapping, grid settings, initial parameter values, and bounds |
| `R/helpers.R` | ID padding, sibling tally lookup, and coordinate rotation |
| `R/prepare_data.R` | Sample filtering, spatial grid, distances, kinship joins, and sampling objects |
| `R/likelihood.R` | Parent–offspring and full-sibling negative log-likelihoods, their sum, and life-history constants |
| `R/hessian_ci.R` | Bounded-logit Hessian and delta-method confidence intervals |
| `R/fit_model.R` | Bounded parameter fitting, uncertainty calculation, and result export |
| `docs/DATA.md` | Required input files and column definitions |
| `docs/MODEL_NOTES.md` | Model assumptions, spatial setup, and interpretation of diagnostics |

The calculation proceeds from sample selection and spatial preparation to kinship likelihood evaluation, parameter fitting, and uncertainty estimation. The runner loads the functions and prepared sampling objects into a shared analysis environment.

## Outputs

A complete run writes:

- `output/barrier_hessian_pop5_ver2.csv`: parameter estimates, standard errors, confidence intervals, fixed total population, convergence and Hessian diagnostics, negative log-likelihood, and elapsed time.
- `output/sessionInfo.txt`: R and package versions for the run.

Later runs overwrite these files. Preparation-only runs do not write analysis results. Inspect convergence, boundary flags, and Hessian diagnostics when interpreting fitted estimates.

For an interactive full run:

```r
analysis <- source("run_example.R")$value
analysis$fit_summary
analysis$hessian_result$summary
print(analysis$grouped_trap_map)
```

## Input integrity

The optional `checksums.sha256` file records SHA-256 fingerprints for the four input CSVs. To verify their contents:

```sh
shasum -a 256 -c checksums.sha256
```

The manifest also includes the additional season-1 file, which is not required for this example. Its check will report a missing file if it is not present locally.

## License

The analysis code is distributed under the GNU General Public License version 3. See [LICENSE](LICENSE).
