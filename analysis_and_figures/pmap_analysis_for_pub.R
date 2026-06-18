require(terra)
require(sf)
require(sfheaders)
require(ggplot2)
require(ggpubr)

# Set islands to TRUE if dealing with lagoon islets of significant size.
# Otherwise make FALSE.
islands <- TRUE

# Load in 10-meter probability map.
pmap <- rast('~/path/to/prob_map')

# Load in the labeled training and testing data, so that we can merge
# the labeled data with the predictions. This way, the predicted probability
# map won't have gaps where the training and testing data are.
labeled <- st_read(paste0('~/shapefile_filename.shp'))
labeled <- labeled[which(!is.na(labeled$label)), ]  # Isolate labeled data.

# Fill in blank spots in probability map (from which testing and training
# data were taken) with the labels of the testing and training data.
labeled_centroids <- st_centroid(labeled)
labeled_coordinates <- st_coordinates(labeled_centroids)
labeled_cells <- cellFromXY(pmap, labeled_coordinates)
pmap[labeled_cells] <- labeled$label

# Load in the masks we made with Sentinel-2 data. The first masks both
# mangroves and water; the second masks only water.
mask <- rast(paste0('~/total_mask_10m.tif'))
water_mask <- rast(paste0('~/nir_water_mask_10m.tif'))

# Calculate total area that is not water, in m ^ 2.
water_mask_df <- as.data.frame(water_mask)
not_water_m2 <- (nrow(water_mask_df) - sum(water_mask_df)) *
  ((res(water_mask)[1] * res(water_mask)[2]))

# Load in 1-meter DTM.
dtm <- rast(paste0('~/DTM_denoised_corrected_kriging_1m.tif'))

# Use the total mask on the probability map.
pmap_mask_int <- intersect(ext(pmap), ext(mask))
pmap_cropped <- crop(pmap, pmap_mask_int)
mask_cropped <- crop(mask, pmap_mask_int)
mask_crop_resamp <- resample(mask_cropped, pmap_cropped)
pmap_masked <- mask(pmap_cropped, mask_crop_resamp, maskvalues = 1)
pmap_not_NA <- pmap_masked
pmap_not_NA[!is.na(pmap_not_NA)] <- 1

# Resample the 1-meter DTM to match the resolution and extent of the
# probability map.
dtm_10m <- resample(dtm, pmap_not_NA)

# Load in the geology shapefile and rasterize it.
geo <- st_read(paste0('~/phn_geology_detailed.shp'))
geo_vec <- vect(geo)
geo_rast <- rasterize(geo_vec, dtm_10m, field = 'Type')

# Do the same for the soils shapefile.
soil <- st_read(paste0('~/phn_soils.shp'))
soil_vec <- vect(soil)
soil_rast <- rasterize(soil_vec, dtm_10m, field = 'GENSOIL')

# Polygonize the non-NA parts of the probability map.
pmap_vects <- as.polygons(pmap_not_NA, dissolve = TRUE)
pmap_sf <- sf::st_as_sf(pmap_vects)  # Make those vectors into an sf object.
pmap_sf <- sf::st_cast(pmap_sf, 'POLYGON')  # Recast the multipolygon object as
# a polygon object.
pmap_areas <- as.numeric(sf::st_area(pmap_sf))  # Calculate the polygons' areas.
pmap_sf$area <- pmap_areas  # Add those areas to the sf object.

# Isolate polygon for Pohnpei mainland (the largest polygon).
pohnpei_mainland_sf <- pmap_sf[which(pmap_sf$area == max(pmap_sf$area)), ]

# Identify polygons for Pohnpei mainland plus all islets of significant size,
# minus airport.
significant_la_ids <- c(2, 9, 11, 14, 55, 138, 157, 254, 271, 322, 359, 370,
                        372, 659, 360, 272, 227, 85, 88, 112, 80, 660)
significant_sf <- pmap_sf[significant_la_ids, ]
surveyed_area <- as.numeric(sum(st_area(significant_sf)))

# Fill in that big data gap in the center.
pohnpei_main_filled <- sf_remove_holes(pohnpei_mainland_sf)
pohnpei_mainland_area <- as.numeric(st_area(pohnpei_main_filled))
data_gap_area <- pohnpei_mainland_area -
  as.numeric(st_area(pohnpei_mainland_sf))
mainland_surveyed <- pohnpei_mainland_area - data_gap_area
significant_sf_filled <- sf_remove_holes(significant_sf)
significant_islands_area <- as.numeric(sum(st_area(significant_sf_filled)) -
                                         st_area(pohnpei_main_filled))

# Make a polygon of the gap.
island_w_gap_geom <- st_geometry(pohnpei_mainland_sf)[[1]]
hole_indices <- 2:length(island_w_gap_geom)
hole_areas <- sapply(hole_indices, function(i) {
  st_area(st_sfc(st_polygon(list(island_w_gap_geom[[i]])), crs = st_crs(pohnpei_mainland_sf)))
})
largest_hole_ind <- hole_indices[which.max(hole_areas)]

gap_ring <- island_w_gap_geom[[largest_hole_ind]]
gap_poly <- st_polygon(list(gap_ring))
gap_poly <- st_sfc(gap_poly, crs = 'epsg:32657')
gap_poly <- st_sf(geometry = gap_poly)

# Convert our data gap polygon to a linestring, which we can use to extract
# distance-to-lagoon values for every cell that falls along the polygon
# boundary.
gap_bound <- st_boundary(gap_poly)

if (islands == TRUE) {
  hole_sf <- significant_sf
  filled_sf <- significant_sf_filled
} else {
  hole_sf <- pohnpei_mainland_sf
  filled_sf <- pohnpei_main_filled
}

# Find the island's centroid.
island_center <- st_coordinates(st_centroid(pohnpei_main_filled))

# Find the numbers of all cells for which there is a predicted probability.
pred_cells <- cells(pmap_masked)

# Get their coordinates.
pred_cell_coords <- xyFromCell(pmap_masked, pred_cells)

# Get angle between each of those coordinates and the center of the island.
angle_from_center <- function(coordinate_row) {
  lon <- pred_cell_coords[coordinate_row, 1] - island_center[1]
  lat <- pred_cell_coords[coordinate_row, 2] - island_center[2]
  theta <- 90 - (atan2(lat, lon)) * (180 / pi)
  if (theta < 0) {
    theta <- 360 + theta
  }
  return(theta)
}
thetas <- sapply(1:nrow(pred_cell_coords), angle_from_center)

# Populate the 10-meter-resolution map of Pohnpei with those angles.
angle_map <- pmap_masked
angle_map[pred_cells] <- thetas

# Mask that angle map of all but the Pohnpei mainland.
angle_masked <- mask(angle_map, vect(hole_sf))

# Create distance-to-lagoon map.
dist_map <- pmap_masked
values(dist_map) <- 1
dist_map <- mask(dist_map, vect(significant_sf_filled), inverse = TRUE)
dist_map <- distance(dist_map)
dist_map_masked <- mask(dist_map, vect(significant_sf))

# Find the maximum margin, inward from the lagoon, at which there is data
# all around the island.
dist_lagoon_to_gap <- extract(dist_map, vect(gap_bound))
uni_lag_dist <- dist_map_masked
uni_lag_dist[uni_lag_dist > min(dist_lagoon_to_gap)] <- NA

# Mask the DTM of all but Pohnpei mainland.
pohn_dtm <- mask(dtm_10m, vect(hole_sf))

# Mask the probability map of all but the Pohnpei mainland.
pmap_masked <- mask(pmap_masked, vect(hole_sf))

# Mask the geology map of all but the Pohnpei mainland.
geo_masked <- mask(geo_rast, vect(hole_sf))

# Mask the soils map of all but the Pohnpei mainland.
soil_masked <- mask(soil_rast, vect(hole_sf))

# Generate slope and aspect layers from the DTM.
dtm_prod <- terra::terrain(pohn_dtm, v = c('slope', 'aspect'))

# Generate a binary prediction map.
binary <- pmap_masked
binary[binary > 0.5] <- 1
binary[binary <= 0.5] <- 0

# Create a binary map of mainland-Pohnpei-or-not.
mainland <- pmap_masked
values(mainland) <- 1
mainland <- mask(mainland, vect(pohnpei_mainland_sf), updatevalue = 0)

# Stack the distance-to-lagoon map, probability map, and 10-meter DTM.
rstack <- c(angle_masked, dist_map_masked, geo_masked, soil_masked,
            pmap_masked, binary, pohn_dtm, dtm_prod[[1]], dtm_prod[[2]],
            mainland)
names(rstack) <- c('angle', 'dist_to_lagoon', 'geo_type', 'soil_type',
                   'prob', 'binary', 'elev', 'slope', 'aspect',
                   'mainland')

# Convert that stack to a data frame.
rstack_df <- as.data.frame(rstack)

# Cut out incomplete cases.
rstack_df_clean <- rstack_df[complete.cases(rstack_df), ]

############ OVERALL POSITIVE CLASSIFICATION ##################################
# Proportion of surveyed area classified as positive (remember, this is 
# considering only the mainland -- that is, no mangroves, lagoon, or lagoon
# islets):
pos_pred_prop <- sum(rstack_df_clean$binary == TRUE) / 
  nrow(rstack_df_clean)

# We have this many square meters/square km of positive-predicted area:
pos_pred_m2 <- sum(rstack_df_clean$binary == TRUE) * ((res(pmap)[1] * res(pmap)[2]))
pos_pred_km2 <- pos_pred_m2 / 1000000

mainland_df <- rstack_df_clean[which(rstack_df_clean$mainland == 1), ]
mainland_pos_m2 <- sum(mainland_df$binary == TRUE) * ((res(pmap)[1] * res(pmap)[2]))

key_labs <- c('mainland_area', 'mainland_area_surveyed',
              'surveyed_lagoon_islands_area', 'area_surveyed',
              'area_classified_pos', 
              'mainland_area_classified_pos',
              'prop_pos_surveyed',
              'prop_pos_mainland_plus_islets')
key_stats_m2 <- c(pohnpei_mainland_area, mainland_surveyed,
                  (surveyed_area - mainland_surveyed),
                  surveyed_area, pos_pred_m2,
                  mainland_pos_m2,
                  (pos_pred_m2 / surveyed_area),
                  (pos_pred_m2 /
                     (pohnpei_mainland_area +
                        (surveyed_area - mainland_surveyed))))
key_results <- data.frame(key_labs, key_stats_m2)

key_results_path <- paste0('~/Documents/CSRM/Further_Pohnpei_work/',
                           'RRNES_grant/Publications/lidar_methods_paper/Drafts/Figures/',
                           'key_results.csv')

if (!file.exists(key_results_path)) {
  write.csv(t(key_results), key_results_path)
}


################ SOIL #########################################################
soil_type <- levels(rstack_df_clean$soil_type)
mean_per_st <- sd_per_st <- numeric(length = length(soil_type))
for (i in 1:length(soil_type)) {
  mean_per_st[i] <- mean(rstack_df_clean[which(rstack_df_clean$soil_type == soil_type[i]), 'prob'], na.rm = TRUE)
  sd_per_st[i] <- sd(rstack_df_clean[which(rstack_df_clean$soil_type == soil_type[i]), 'prob'], na.rm = TRUE)
}
soil_df <- data.frame(soil_type, mean_per_st, sd_per_st)
soil_bars <- ggplot(soil_df, aes(soil_type, mean_per_st)) + geom_bar(stat = 'identity') + geom_errorbar(aes(ymin = mean_per_st - sd_per_st,
                                                                                                            ymax = mean_per_st + sd_per_st))
ump_mean <- rstack_df_clean[which(rstack_df_clean$soil_type == 'Umpump-Rakied'), 'prob']
not_ump_mean <- rstack_df_clean[which(rstack_df_clean$soil_type != 'Umpump-Rakied' &
                                        rstack_df_clean$soil_type != 'Naniak'), 'prob']
set.seed(2113)
not_ump_mean_sample <- sample(not_ump_mean, size = length(ump_mean))
t_test_ump <- t.test(ump_mean, not_ump_mean_sample, alternative = 'greater')

comp_df <- data.frame(ump_mean, not_ump_mean_sample)
colnames(comp_df) <- c('Umpump', 'all_else')
soil_prob_hist_overlay <- ggplot() +
  geom_histogram(aes(x = comp_df$Umpump, fill = 'Umpump'), alpha = 0.5) +
  geom_histogram(aes(x = comp_df$all_else, fill = 'all_else'), alpha = 0.5)

########## DIRECTION ##########################################################
# Define function to round to nearest 10. Apply to all angles from island center.
round_10 <- function(x) {
  if (is.na(x)) return(NA)
  if (x %% 10 != 0) {
    remainder <- x %% 10
    if (remainder >= 5) {
      x <- x + (10 - remainder)
    } else {
      x <- x - remainder
    }
  }
  return(x)
}

# Work only with data falling within the uniform margin from the lagoon -- that is,
# within the maximum distance from the lagoon for which we have data all around
# the island.
margin_stack <- c(angle_masked, dist_map_masked, binary, uni_lag_dist)
names(margin_stack) <- c('angle', 'dist_to_lagoon', 'binary', 'uni_lag_margin')
margin_df <- as.data.frame(margin_stack)
margin_df <- margin_df[complete.cases(margin_df), ]

# Cut out incomplete cases.
rstack_df_clean <- rstack_df[complete.cases(rstack_df), ]

# Round angles to nearest 10 degrees.
margin_df$theta_bin_10 <- as.factor(sapply(margin_df$angle, round_10))

# Put 0 and 360 bins in the same bin.
margin_df[which(margin_df$theta_bin_10 == 0), 'theta_bin_10'] <- 360

# Calculate proportion of positive-classified cells per 10-degree slice.
bins <- as.factor(levels(margin_df$theta_bin_10))
prop <- numeric(length = length(bins))
ind <- 1
for (i in bins) {
  cur_bin <- margin_df[which(margin_df$theta_bin_10 == i), 'binary']
  prop[ind] <- sum(cur_bin, na.rm = TRUE) / length(cur_bin)
  ind <- ind + 1
}

# Put those proportions in a data frame and clean it up.
prop_df <- data.frame(bins, prop)
prop_df <- prop_df[complete.cases(prop_df), ]

# Make a rose plot of those proportions.
rose_prob <- ggplot(prop_df, aes(bins, prop)) + geom_col() +
  coord_polar(start = 5 * (pi / 180)) +
  scale_x_discrete(limits = factor(seq(10, 360, 10))) +
  labs(y = 'Proportion of positive-classified cells',
       x = '10-degree sector',
       title = 'Proportion of positive-classified cells per\n10-degree sector')

# For plotting: Make a polygon that shows the margin for which we have 
# data all around the island.
inner_margin <- st_buffer(pohnpei_main_filled, -min(dist_lagoon_to_gap))
data_margin <- st_difference(pohnpei_main_filled, inner_margin)
st_write(data_margin, '~/Documents/CSRM/Further_Pohnpei_work/RRNES_grant/Publications/lidar_methods_paper/GIS/uniform_data_margin.shp')

# Calculate count of positive-classified cells per 10-degree slice.
bins <- as.factor(levels(margin_df$theta_bin_10))
count <- numeric(length = length(bins))
ind <- 1
for (i in bins) {
  cur_bin <- margin_df[which(margin_df$theta_bin_10 == i), 'binary']
  count[ind] <- sum(cur_bin, na.rm = TRUE)
  ind <- ind + 1
}
count_df <- data.frame(bins, count)
count_df <- count_df[complete.cases(count_df), ]

# Make rose plot similar to above, but of total positive counts, rather than of
# proportions per slice.
rose_count <- ggplot(count_df, aes(bins, count)) + geom_col() +
  coord_polar(start = 5 * (pi / 180)) +
  scale_x_discrete(limits = factor(seq(10, 360, 10))) +
  labs(y = 'Number of positive-classified cells',
       x = '10-degree slice',
       title = 'Number of positive-classified cells per 10-degree slice')

############ DISTANCE TO LAGOON ###############################################
# Bin the observations based on their distance to the lagoon.
# Define a function to cut a vector of distances into bins of a determined
# size.
des_binwidth <- 20

distance_bins <- function(input_dist, binwidth = 25) {
  (input_dist %/% binwidth) + 1
}

# Apply that function to the distances to the lagoon.
rstack_df_clean$dist_bin <- sapply(rstack_df_clean$dist_to_lagoon,
                                   distance_bins, binwidth = des_binwidth)

# Put those distances into a data frame.
bin_counts <- bin_mean_prob <- numeric(length = max(rstack_df_clean$dist_bin, na.rm = TRUE))
for (i in 1:max(rstack_df_clean$dist_bin, na.rm = TRUE)) {
  cur_sum <- sum(rstack_df_clean[which(rstack_df_clean$dist_bin == i), 'binary'], na.rm = TRUE)
  bin_counts[i] <- cur_sum
  cur_mean <- mean(rstack_df_clean[which(rstack_df_clean$dist_bin == i), 'prob'], na.rm = TRUE)
  bin_mean_prob[i] <- cur_mean
}
dist_hist_df <- data.frame(1:max(rstack_df_clean$dist_bin, na.rm = TRUE), bin_counts, bin_mean_prob)
colnames(dist_hist_df) <- c('bin', 'counts', 'mean_prob')

# Make a column plot.
dist_col <- ggplot(dist_hist_df, aes(bin * des_binwidth, mean_prob,
                                     color = mean_prob,
                                     fill = mean_prob)) +
  geom_col() + 
  scale_color_gradient(low = '#0571b0', high = '#ca0020') +
  scale_fill_gradient(low = '#0571b0', high = '#ca0020') +
  labs(y = 'Mean predicted probability', x = 'Distance from lagoon, m') + 
  xlim(c(0, 4140)) +
  guides(fill = 'none', color = 'none')

dist_col

# Define a function to rid distances of outliers (> 3 sd from the mean).
cut_outliers <- function(input_vector) {
  mu <- mean(input_vector, na.rm = TRUE)
  std <- sd(input_vector, na.rm = TRUE)
  if (is.na(std) | is.na(mu)) return(NA)
  cutoffs <- c((mu - (3 * std)), (mu + (3 * std)))
  if (cutoffs[1] < 0) cutoffs[1] <- -0.001
  return(cutoffs)
}

binned_probs <- numeric()
for (i in 1:max(rstack_df_clean$dist_bin)) {
  cur_probs <- rstack_df_clean[which(rstack_df_clean$dist_bin == i), 'prob']
  out_prob <- mean(cur_probs, na.rm = TRUE)
  binned_probs <- append(binned_probs, out_prob)
}
bin_df <- data.frame(as.numeric(c(1:max(rstack_df_clean$dist_bin))), binned_probs)
bin_df$bin <- bin_df$bin * 20  # To account for 20-meter distances
colnames(bin_df) <- c('bin', 'prob')
bin_dist_plot <- ggplot(bin_df, aes(bin)) + geom_bar(aes(weight = prob)) +
  labs(x = 'Distance from lagoon (m)',
       y = 'Mean probability')
bin_dist_point <- ggplot(bin_df, aes(bin, prob)) + geom_point(size = 0.5)

# Heatmap, mean probability against distance to lagoon.
heatmap <- dist_map_masked
heatmap <- (heatmap %/% des_binwidth)
for (i in 1:max(values(heatmap), na.rm = TRUE)) {
  heatmap[heatmap == i] <- bin_df[which(bin_df$bin == i), 'prob']
  print(paste0('Done ', i, ' of ', max(values(heatmap), na.rm = TRUE)))
}
heatmap_norm <- (heatmap - min(values(heatmap), na.rm = TRUE)) / 
  (max(values(heatmap), na.rm = TRUE) - min(values(heatmap), na.rm = TRUE))

