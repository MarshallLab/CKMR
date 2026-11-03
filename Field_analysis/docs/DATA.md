# Input data

Input data are not distributed with the code. Create a local `data/` directory containing the three files required for the default analysis: `season2_filtered_v2.csv`, `Kinship_PO_v3.csv`, and `Kinship_FS_v3.csv`.

The input files contain sampling records and inferred kinship pairs. The analysis filters and joins these tables in memory. SHA-256 fingerprints are recorded in the project-level `checksums.sha256`.

| File | Use in this example |
| --- | --- |
| `season2_filtered_v2.csv` | Sampling records used by the example; observed floors 2, 7 and 12 only |
| `Kinship_PO_v3.csv` | Parent–offspring pairs; retain pairs whose two IDs occur in the filtered sample data |
| `Kinship_FS_v3.csv` | Full-sibling pairs; retain pairs whose two IDs occur in the filtered sample data |
| `season1_filtered_v2.csv` | Additional season-1 sampling records; not used by this example |

## Columns used by the code

For season-2 sampling records:

- `ID`: individual identifier; converted to character and left-padded to at least four characters.
- `Date_of_collection`: parsed with `%Y%m%d`.
- `Block`, `Floor`: sample selection and spatial grouping.
- `Unit`: used to reassign block 844 units 203 and 207 to block 845.
- `SVY21_X`, `SVY21_Y`: projected coordinates, interpreted using EPSG:3414.

For both kinship files:

- `D1_indiv`, `D2_indiv`: individual IDs, padded in the same way as sample IDs.
- `D1_date`, `D2_date`: dates parsed with `%Y-%m-%d`.

The parent–offspring file additionally uses `PO_direction`: `D1` selects the first individual as mother; other values select the second individual as mother. Pair assignments are read from the input table without reclassification.

Patch IDs are obtained by joining filtered sample locations, rather than from kinship-file coordinates. The full-sibling likelihood uses connected components of the retained FS pair graph as family assignments.

The default floor mapping is `2 -> 1`, `7 -> 3`, `12 -> 5`. The matrix retains five movement layers. No samples from floors 6, 10 or 11 are pooled into the selected floors.
