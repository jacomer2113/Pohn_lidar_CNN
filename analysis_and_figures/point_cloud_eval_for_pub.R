require(terra)
require(sf)
require(lidR)
require(ggplot2)

# load in water mask
water_mask <- rast('~/nir_water_mask_10m.tif')

# and DTM
dtm <- rast('DTM_denoised_corrected_kriging_1m.tif')

# Resample the DTM to match the resolution and extent of the water mask.
dtm_resamp <- resample(dtm, water_mask, method = 'bilinear')
dtm_resamp <- subst(dtm_resamp, NA, min(values(dtm_resamp), na.rm = TRUE))

# Make the resampled DTM binary.
dtm_r_bin <- dtm_resamp
dtm_r_bin[dtm_r_bin < -2164] <- 0
dtm_r_bin[dtm_r_bin != 0] <- 2

# Identify the overlap between the DTM and the water.
overlap <- dtm_r_bin + water_mask

# Isolate only the area that has both been scanned and is not water.
positive <- overlap
positive[positive != 2] <- NA

# Polygonize the surveyed land areas.
scanned_land_vect <- as.polygons(positive, dissolve = TRUE)
scanned_land <- sf::st_as_sf(scanned_land_vect)  # Make those vectors into an sf object.
scanned_land <- sf::st_cast(scanned_land, 'POLYGON')  # Recast the multipolygon object as
# a polygon object.
scanned_areas <- as.numeric(sf::st_area(scanned_land))  # Calculate the polygons' areas.
scanned_land$area <- scanned_areas  # Add those areas to the sf object.

scanned_filt <- scanned_land[which(scanned_land$area > 500), ] # Filter out
# polygons smaller than 500 m^2, many of which seem to result from small areas
# of spurious data over the lagoon.

print(paste0(sum(scanned_filt$area) / 1000000, ' square km of land surveyed.'))

# What is the area of the data gap in the center?
central_box <- ext(c(410000, 419000, 756000, 766500))
central_area <- crop(positive, central_box)
ca_area <- (nrow(central_area) * 10) * (ncol(central_area) * 10)
pos_ca_area <- expanse(central_area)
neg_ca_area <- ca_area - pos_ca_area[2]

print(paste0(neg_ca_area / 1000000, ' square km of land in island center missed.'))

# Read in our las catalog.
cat <- readLAScatalog('~/las_ellip_noise_classified/')

# Pick out the sf POLYGON for mainland Pohnpei. For some reason, the airport
# is not included here, but that's okay.
pohn_sf <- scanned_land[which(scanned_land$area == max(scanned_land$area, na.rm = TRUE)), ]
islet_sf <- scanned_land[which(scanned_land$area != max(scanned_land$area, na.rm = TRUE)), ]

# Set samp_and_plot equal to FALSE once this has been done successfully,
# because the sampling takes forever.
samp_and_plot <- TRUE
if (samp_and_plot == TRUE) {
  # Define a function to grab a random sample of 250 circular areas of 
  # 250 m^2 each from across a given input polygon.
  circ_samp <- function(input_shp) {
    samp <- st_sf(st_sample(input_shp, 1000))
    area_samp <- st_buffer(samp, dist = sqrt(250 / pi))
    success_index <- 0
    i <- 1
    while (success_index < 250) {
      cur_cat <- clip_roi(cat, area_samp[i, ])
      if (nrow(cur_cat@data) == 0) {
        i <- i + 1
        print('fail')
        next
      } else {
        cur_las <- normalize_height(cur_cat, dtm)
        cur_df <- data.frame(cur_las@data$Z, cur_las@data$ReturnNumber)
        if (success_index == 0) {
          sum_df <- cur_df
        } else {
          sum_df <- rbind(sum_df, cur_df)
        }
        i <- i + 1
        success_index <- success_index + 1
        print(paste0('Done ', success_index, ' of ', 250))
      }
    }
    colnames(sum_df) <- c('height', 'return_num')
    sum_df$return_num <- as.factor(sum_df$return_num)
    sum_df$return_num <- factor(sum_df$return_num,
                                levels = levels(sum_df$return_num)[order(
                                  levels(sum_df$return_num), decreasing = TRUE)])
    return(sum_df)
  }
  
  # Grab that sample from across the Pohnpei mainland . . . 
  set.seed(4303)
  pohn_main_samp <- circ_samp(pohn_sf)
  
  # . . . the lagoon islets only . . . 
  lagoon_islet_samp <- circ_samp(islet_sf)
  
  # . . . and from all the surveyed land.
  pohn_all_samp <- circ_samp(scanned_land)
  
  # Make a vertical histogram for mainland Pohnpei only.
  main_point_height_plot <- ggplot(subset(pohn_main_samp, return_num %in% 1:3),
                                   aes(y = height, fill = return_num)) +
    geom_histogram(binwidth = 0.5) +
    labs(x = 'Number of returns',
         y = 'Height above ground (m)',
         fill = 'Return number',
         title = 'Mainland Pohnpei') +
    ylim(y = c(-1.5, 42)) +
    xlim(y = c(0, 100000)) +
    guides(fill = 'none') + 
    theme(plot.title = element_text(size = 14),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank())
  
  main_point_height_plot
  
  # Do the same for the lagoon islets . . .
  lagoon_point_height_plot <- ggplot(subset(lagoon_islet_samp, return_num %in% 1:3),
                                     aes(y = height, fill = return_num)) +
    geom_histogram(binwidth = 0.5) +
    labs(x = 'Number of returns',
         fill = 'Return number',
         title = 'Lagoon islets') +
    ylim(y = c(-1.5, 42)) +
    xlim(y = c(0, 100000)) +
    guides(fill = 'none') + theme(axis.title.y = element_blank()) + 
    theme(plot.title = element_text(size = 14),
          axis.title = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  lagoon_point_height_plot
  
  # and for all the surveyed land.
  all_point_height_plot <- ggplot(subset(pohn_all_samp, return_num %in% 1:3),
                                  aes(y = height, fill = return_num)) +
    geom_histogram(binwidth = 0.5) +
    labs(x = 'Number of returns',
         fill = 'Return number',
         title = 'All surveyed land') +
    ylim(y = c(-1.5, 42)) +
    xlim(y = c(0, 100000)) +
#    guides(fill = 'none') +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) + 
    
#    guides(fill = guide_legend(reverse = TRUE)) + theme(axis.title.y = element_blank()) + 
    theme(plot.title = element_text(size = 14),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  all_point_height_plot
}
  
# Assess pulse and ground point densities.
pulse <- rast('~/pulse_density_1m_corrected_denoised.tif')
point <- rast('~/point_density_1m_corrected_denoised.tif')
ground <- rast('~/ground_point_density_1m_corrected_denoised.tif')

# Make a rasterized version of the sf object showing scanned land areas.
rast_scanned <- rasterize(vect(scanned_land), pulse)

# Mask the pulse, point, and ground point density layers.
pulse_masked <- mask(pulse, rast_scanned)
point_masked <- mask(point, rast_scanned)
ground_masked <- mask(ground, rast_scanned)

# Compute layers for returns/pulse, ground returns/pulse, and ground returns/return.
point_pulse <- point_masked / pulse_masked
ground_pulse <- ground_masked / pulse_masked
ground_point <- ground_masked / point_masked

# Clean up ground returns/pulse and ground returns/return, which should be <= 1.
ground_pulse[ground_pulse > 1] <- NA
ground_point[ground_point > 1] <- NA

# Define a function to calculate these statistics within an area of interest.
lidar_metrics <- function(input_shp) {
  cur_shp <- st_read(input_shp)
  cur_pulse <- crop(pulse_masked, cur_shp)
  cur_point <- crop(point_masked, cur_shp)
  cur_ground <- crop(ground_masked, cur_shp)
  cur_point_pulse <- crop(point_pulse, cur_shp)
  cur_ground_pulse <- crop(ground_pulse, cur_shp)
  cur_ground_point <- crop(ground_point, cur_shp)
  pulse_mu <- mean(values(cur_pulse), na.rm = TRUE)
  pulse_sd <- sd(values(cur_pulse), na.rm = TRUE)
  point_mu <- mean(values(cur_point), na.rm = TRUE)
  point_sd <- sd(values(cur_point), na.rm = TRUE)
  ground_mu <- mean(values(cur_ground), na.rm = TRUE)
  ground_sd <- sd(values(cur_ground), na.rm = TRUE)
  point_pulse_mu <- mean(values(cur_point_pulse), na.rm = TRUE)
  point_pulse_sd <- sd(values(cur_point_pulse), na.rm = TRUE)
  ground_pulse_mu <- mean(values(cur_ground_pulse), na.rm = TRUE)
  ground_pulse_sd <- sd(values(cur_ground_pulse), na.rm = TRUE)
  ground_point_mu <- mean(values(cur_ground_point), na.rm = TRUE)
  ground_point_sd <- sd(values(cur_ground_point), na.rm = TRUE)
  names <- c('pulse_mu', 'pulse_sd', 'point_mu', 'point_sd', 'ground_mu', 'ground_sd',
             'point_pulse_mu', 'point_pulse_sd', 'ground_pulse_mu', 'ground_pulse_sd',
             'ground_point_mu', 'ground_point_sd')
  mets <- c(pulse_mu, pulse_sd, point_mu, point_sd, ground_mu, ground_sd,
            point_pulse_mu, point_pulse_sd, ground_pulse_mu, ground_pulse_sd,
            ground_point_mu, ground_point_sd)
  output_m <- matrix(c(names, mets), nrow = 2, byrow = TRUE)
  return(output_m)
}

# Apply them to our 4 areas of particular interest.
aoi_dir <- paste0('~/report/AOIs/')
shp_names <- grep('.shp', dir(aoi_dir), value = TRUE)
for (i in 1:length(shp_names)) {
  print(paste0('Generating metrics for ', shp_names[i]))
  cur_matrix <- lidar_metrics(paste0(aoi_dir, shp_names[i]))
  if (i == 1) {
    out_matrix <- cur_matrix
  } else {
    out_matrix <- rbind(out_matrix, cur_matrix[2, ])
  }
  print(paste0('Metrics for ', shp_names[i], ' generated'))
  print('Output matrix so far:')
  print(out_matrix)
}

# Put those metrics into a data frame.
metric_df <- data.frame(out_matrix[2:nrow(out_matrix), ])
colnames(metric_df) <- c('pulse_mu', 'pulse_sd', 'point_mu', 'point_sd',
                         'ground_mu', 'ground_sd', 'point_pulse_mu',
                         'point_pulse_sd', 'ground_pulse_mu', 'ground_pulse_sd',
                         'ground_point_mu', 'ground_point_sd')

# Define a function to calculate the mean and sd of a given input layer.
mu_sd <- function(input_layer) {
  mu <- mean(values(input_layer), na.rm = TRUE)
  std <- sd(values(input_layer), na.rm = TRUE)
  return(c(mu, std))
}

# Apply that function to the entire area of surveyed land.
metric_stack <- c(pulse_masked, point_masked, ground_masked,
                  point_pulse, ground_pulse, ground_point)
total_area_metrics <- numeric()
for (i in 1:nlyr(metric_stack)) {
  print(paste0('Starting with layer ', i))
  cur_metrics <- mu_sd(metric_stack[[i]])
  total_area_metrics <- append(total_area_metrics, cur_metrics[1])
  total_area_metrics <- append(total_area_metrics, cur_metrics[2])
  print(paste0('Current vec of total area metrics: ', total_area_metrics))
  print(paste0('Done with layer ', i))
}

# Add those metrics to the metrics data frame we created above.
metric_df <- rbind(metric_df, total_area_metrics)
metric_df[] <- lapply(metric_df, as.numeric)  # Make the values numeric.
metric_df_rounded <- round(metric_df, 3)  # Round to 3 decimal places.
metric_df_rounded$area <- c('Awak', 'Lohd', 'Sapwtakai', 'Temwen', 'All')
metric_df_rounded <- metric_df_rounded[, c(13, 1:12)]

