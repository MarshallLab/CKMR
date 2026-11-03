# Run after config.R and helpers.R in the shared analysis environment.
# Prepare sample identifiers, dates, spatial locations, and kinship pair tables.

# 1. Filter samples and map observed floors to movement layers ------------------
FL_barrier <- read.csv(sample_file)
FL_barrier <- FL_barrier %>% mutate(ID = pad_id(ID))
FL_barrier$Date_of_collection <- as.Date(as.character(FL_barrier$Date_of_collection), format = "%Y%m%d")

FL_barrier <- FL_barrier %>% filter(Block %in% selected_blocks)
floor_map <- data.frame(
  observed_floor = raw_sampled_floors,
  model_floor = model_sampled_layers
)

FL_barrier <- FL_barrier %>%
  filter(Floor %in% raw_sampled_floors) %>%
  left_join(floor_map, by = c("Floor" = "observed_floor"))
# Update Block value for the specified condition
FL_barrier <- FL_barrier %>%
  mutate(Block = if_else(Block == "844" & Unit %in% c(203, 207), "845", Block))

# 2. Block centroids and rotated spatial grid ----------------------------------
FL_barrier_sf <- st_as_sf(FL_barrier, coords = c("SVY21_X", "SVY21_Y"), crs = 3414)
block_sf <- FL_barrier_sf %>%
  group_by(Block) %>%
  summarise(geometry = st_centroid(st_union(geometry))) %>%
  ungroup()

block_xy <- block_sf %>%
  mutate(
    x = st_coordinates(geometry)[, 1],
    y = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry() # keep Block, x, y

# Use block centroids as the grouped trap locations before snapping them to the standardized grid.
grouped_traps <- FL_barrier %>%
  distinct(Block, Floor, model_floor) %>%
  left_join(block_xy, by = "Block") %>%
  mutate(grouped_trap_id = paste0("B", Block, "_F", Floor)) %>%
  st_as_sf(coords = c("x", "y"), crs = 3414)

grouped_traps_rotated <- rotate_points(grouped_traps, angle_deg = rotation_angle)

bbox <- st_bbox(grouped_traps_rotated)
bbox_expanded <- bbox
bbox_expanded["xmin"] <- floor(bbox["xmin"] / grid_size) * grid_size - 2 * grid_size
bbox_expanded["ymin"] <- floor(bbox["ymin"] / grid_size) * grid_size - 2 * grid_size
bbox_expanded["xmax"] <- ceiling(bbox["xmax"] / grid_size) * grid_size + 2 * grid_size
bbox_expanded["ymax"] <- ceiling(bbox["ymax"] / grid_size) * grid_size + 2 * grid_size

grid <- st_make_grid(
  st_as_sfc(bbox_expanded),
  cellsize = grid_size,
  square = TRUE
) %>%
  st_sf()

grid <- grid %>%
  mutate(centroid = st_centroid(geometry)) %>%
  mutate(
    cx = st_coordinates(centroid)[, 1],
    cy = st_coordinates(centroid)[, 2]
  ) %>%
  arrange(desc(cy), cx) %>%
  mutate(patch_id_2d = row_number())

grouped_traps_rotated <- st_join(grouped_traps_rotated, grid, join = st_within)
grouped_traps_rotated$patch_id_2d <- st_within(grouped_traps_rotated, grid) %>% as.integer()

if (any(is.na(grouped_traps_rotated$patch_id_2d))) {
  stop("Some grouped traps were not assigned to a standardized grid cell.")
}

patches_per_floor <- nrow(grid)

comb <- expand.grid(
  Floor = floors,
  patch_id_2d = grid$patch_id_2d
) %>%
  arrange(Floor, patch_id_2d) %>%
  mutate(
    x = grid$cx[match(patch_id_2d, grid$patch_id_2d)],
    y = grid$cy[match(patch_id_2d, grid$patch_id_2d)],
    z = Floor * floor_height,
    patch_id = patch_id_2d + patches_per_floor * (Floor - 1)
  )

# 3. Stack the grid across layers and calculate 3D distances --------------------
coords_block3d <- as.matrix(comb[, c("x", "y", "z")])
rownames(coords_block3d) <- comb$patch_id

distMat <- as.matrix(dist(coords_block3d))

grouped_traps <- grouped_traps_rotated %>%
  st_drop_geometry() %>%
  mutate(
    patch_id = patch_id_2d + patches_per_floor * (model_floor - 1)
  ) %>%
  dplyr::select(Block, Floor, model_floor, grouped_trap_id, patch_id, patch_id_2d)

FL_barrier <- FL_barrier %>%
  mutate(grouped_trap_id = paste0("B", Block, "_F", Floor)) %>%
  left_join(grouped_traps %>% dplyr::select(grouped_trap_id, patch_id, patch_id_2d), by = "grouped_trap_id")

grouped_traps_plot <- FL_barrier %>%
  st_drop_geometry() %>%
  distinct(grouped_trap_id, Floor, patch_id_2d) %>%
  left_join(grid %>% st_drop_geometry() %>% dplyr::select(patch_id_2d, cx, cy), by = "patch_id_2d") %>%
  st_as_sf(coords = c("cx", "cy"), crs = 3414)

grouped_trap_map <- ggplot() +
  geom_sf(data = grid, fill = NA, color = "gray70") +
  geom_sf(data = grouped_traps_plot, aes(color = factor(Floor)), size = 2.5) +
  geom_sf_text(data = grid, aes(label = patch_id_2d), size = 2, color = "black") +
  geom_sf_text(data = grouped_traps_plot, aes(label = grouped_trap_id), size = 2.5, nudge_y = 4) +
  theme_minimal() +
  labs(
    title = "Grouped trap locations over rotated standardized grid",
    color = "Observed floor"
  )

# Inspect dimensions and the first distance entries.
print(dim(distMat)) # should be (#grid cells * #floors)^2
print(head(distMat[, 1:5])) # first few entries

# 4. Read kinship pairs and join the prepared sample locations ------------------
kinship_PO_barrier <- read.csv(po_file)
kinship_FS_barrier <- read.csv(fs_file)

kinship_PO_barrier <- kinship_PO_barrier %>% mutate(
  D1_indiv = pad_id(D1_indiv),
  D2_indiv = pad_id(D2_indiv),
  D1_date = as.Date(as.character(D1_date), format = "%Y-%m-%d"),
  D2_date = as.Date(as.character(D2_date), format = "%Y-%m-%d")
)

kinship_FS_barrier <- kinship_FS_barrier %>% mutate(
  D1_indiv = pad_id(D1_indiv),
  D2_indiv = pad_id(D2_indiv),
  D1_date = as.Date(as.character(D1_date), format = "%Y-%m-%d"),
  D2_date = as.Date(as.character(D2_date), format = "%Y-%m-%d")
)

# Prepare minimal lookup table from FL_barrier
trap_lookup <- FL_barrier %>%
  dplyr::select(ID, patch_id, grouped_trap_id) %>%
  distinct()

kinship_PO_barrier <- kinship_PO_barrier %>% filter(D1_indiv %in% FL_barrier$ID & D2_indiv %in% FL_barrier$ID)
kinship_PO_barrier <- kinship_PO_barrier %>%
  left_join(trap_lookup, by = c("D1_indiv" = "ID")) %>%
  rename(D1_Trap_ID = patch_id, D1_grouped_trap_id = grouped_trap_id) %>%
  left_join(trap_lookup, by = c("D2_indiv" = "ID")) %>%
  rename(D2_Trap_ID = patch_id, D2_grouped_trap_id = grouped_trap_id)

kinship_FS_barrier <- kinship_FS_barrier %>% filter(D1_indiv %in% FL_barrier$ID & D2_indiv %in% FL_barrier$ID)
kinship_FS_barrier <- kinship_FS_barrier %>%
  left_join(trap_lookup, by = c("D1_indiv" = "ID")) %>%
  rename(D1_Trap_ID = patch_id, D1_grouped_trap_id = grouped_trap_id) %>%
  left_join(trap_lookup, by = c("D2_indiv" = "ID")) %>%
  rename(D2_Trap_ID = patch_id, D2_grouped_trap_id = grouped_trap_id)

# 5. Sampling objects consumed by the likelihood functions ---------------------
# Preserve row order: the sibling likelihood searches later rows for relatives.
sampledAdultFemales <- FL_barrier
sampledAdultFemales <- sampledAdultFemales %>% st_drop_geometry()
sampledPatchesAdultFemales <- sort(unique(sampledAdultFemales$patch_id))
numSampledPatchesAdultFemales <- length(sampledPatchesAdultFemales)

numSampledLarvae <- 0
numSampledPupae <- 0
numSampledAdultMales <- 0
numSampledAdultFemales <- nrow(sampledAdultFemales)
numSampledAdults <- nrow(sampledAdultFemales)

numSampledAdults <- numSampledAdultMales + numSampledAdultFemales
sampledAdults <- sampledAdultFemales
sampledPatchesAdults <- sort(unique(sampledAdultFemales$patch_id))
numSampledPatchesAdults <- length(sampledPatchesAdults)
tSamplingStart <- min(sampledAdults$Date_of_collection)
tSamplingEnd <- max(sampledAdults$Date_of_collection)
sampledPatches <- sort(unique(sampledAdultFemales$patch_id))
numSampledPatches <- length(sampledPatches)
# Time period over which likelihood calculations will be performed:
samplingDays <- seq(tSamplingStart, tSamplingEnd, by = 1)
# Total number of days over which sampling occurred:
numSamplingDays <- tSamplingEnd - tSamplingStart + 1
