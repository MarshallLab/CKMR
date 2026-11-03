# Model specification

## Spatial representation

- The analysis selects season-2 records from blocks `943:953`, `835`, and `837:850`.
- Observed floors 2, 7, and 12 map to model layers 1, 3, and 5. Layers 2 and 4 are unsampled intermediate movement layers. Not every block has observations from every sampled floor.
- Block 844 records for units 203 and 207 are assigned to block 845 before spatial grouping.
- Block centroids are calculated from the filtered records. Each block/floor group uses its block centroid as its horizontal location.
- Centroids are rotated using `rotation_angle = -50` degrees and assigned to a grid with 45-metre square cells. The bounding box is expanded by two cells on each side.
- The grid is repeated across five layers with 7 metres between adjacent layers. A 3D patch ID is `patch_id_2d + patches_per_floor * (model_floor - 1)`.
- With these inputs, the grid contains 182 cells per layer, giving 910 movement nodes. Euclidean distances are calculated from the 3D coordinates.
- Both likelihoods use `barrier_node = 98`. Within each layer, indices `<= 98` and `> 98` define the two sides of the barrier.

The grid and centroids are calculated from the selected records. Adapting the example to another area or sample selection requires reviewing the grid and barrier threshold together.

## Movement and population

Daily movement weights are normalized by origin node. Same-layer movement uses `exp(-lambda * distance)`. Cross-layer movement uses `exp(-lambda * horizontal_distance - lambda_ver * vertical_distance)`. Movement across the barrier is multiplied by `1 - barrier_strength`.

Both distance components are measured in metres, so both dispersal rates are applied to metric distances. Successive powers of the daily transition matrix give multi-day movement probabilities.

The likelihood uses a fixed adult population per node, `N_A = 5`, with `N_F = N_A / 2`. The summary field `N_A_fixed` reports the total across nodes: `node_population_size * nrow(distMat)`.

## Life history and sampling

Life-history settings are defined in `logLike_all()`:

| Setting | Value |
| --- | --- |
| Adult daily mortality `mu_A` | 0.09 |
| Egg duration `T_E` | 2 days |
| Larval duration `T_L` | 5 days |
| Pupal duration `T_P` | 1 day |
| Maximum adult lifespan `T_A` | 40 days |
| Daily growth rate `rM` | 1.175 |
| Daily egg production `beta` | 20 |

Generation time and juvenile mortality rates are derived within the objective function. All selected sample records are treated as adult females. Parent–offspring direction is read from the pair table. Full-sibling family membership is assigned using connected components of the retained full-sibling pair graph.

## Estimation and confidence intervals

The objective is the sum of parent–offspring and full-sibling negative log-likelihoods. Parameters are ordered as `lambda`, `lambda_ver`, and `barrier_strength`, with initial values `1/15.3`, `1/5`, and `0.75`. Each parameter is bounded between `0.001` and `0.999`.

Fitting uses a single start with `optimx(method = "nlminb")` and controls `reltol = 1e-5`, `starttests = FALSE`, and `kkt = FALSE`. Some versions of `nlminb` report `reltol` as an unrecognized control element and ignore it.

Confidence intervals use `stats::optimHess()` on a bounded-logit parameter scale. The inverse Hessian is mapped back to the parameter scale by the delta method. The summary reports 95% Wald intervals and versions clipped to the parameter bounds.

`hessian_ok` checks finiteness of the covariance and standard errors. Positive definiteness and condition number are reported separately. Intervals are not automatically suppressed for an indefinite Hessian or a boundary estimate; inspect these diagnostics alongside the optimizer convergence code.

Probability-scale expressions such as `0 * log(0)` can produce non-finite objective values at some parameter settings. Sample row order also matters for the sibling-pair lookup, which searches subsequent rows for relatives.
