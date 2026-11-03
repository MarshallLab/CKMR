# Small helpers used by data preparation and the full-sibling likelihood.
# Return a full-sibling tally at a time-lag index, or zero for an empty tally.
calc_numFullSib <- function(a, b) {
  if (length(a) > 0) {
    return(a[b])
  } else {
    return(0)
  }
}

# Pad identifiers before joining samples to kinship pairs.
pad_id <- function(x) str_pad(as.character(x), width = 4, pad = "0")

# Rotate block centroids around their mean in the projected coordinate system.
rotate_points <- function(sf_obj, angle_deg, center = NULL) {
  angle_rad <- angle_deg * pi / 180
  coords <- st_coordinates(sf_obj)

  if (is.null(center)) {
    center <- colMeans(coords)
  }

  shifted <- sweep(coords, 2, center)
  rot_mat <- matrix(c(cos(angle_rad), -sin(angle_rad), sin(angle_rad), cos(angle_rad)), ncol = 2)
  rotated <- t(rot_mat %*% t(shifted))
  rotated_coords <- sweep(rotated, 2, center, "+")

  sf_obj_rotated <- st_as_sf(as.data.frame(rotated_coords), coords = 1:2, crs = st_crs(sf_obj))
  st_geometry(sf_obj_rotated) <- st_sfc(st_geometry(sf_obj_rotated), crs = st_crs(sf_obj))
  bind_cols(sf_obj_rotated, sf_obj %>% st_drop_geometry())
}
