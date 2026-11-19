require(terra)
require(sf)
require(lidR)
require(ggplot2)

# Define a function to take, from an input set of points, a circular sample 
# area of 250 m ^ 2 from the point cloud.
circ_samp <- function(samp_points, in_las_ctg, n_success = 250, area_m2 = 250) {
  area_samp <- st_buffer(samp_points, dist = sqrt(area_m2 / pi))
  results <- vector('list', n_success)
  success_index <- 0
  i <- 1
  
  while (success_index < n_success) {
    i <- i + 1
    if (i > nrow(area_samp)) {
      stop(paste0('Only ', success_index, ' of ', n_success,
                  ' samples contained points; provide more sample points.'))
    }
    
    cur_las <- clip_roi(in_las_ctg, area_samp[i, ])
    if (is.null(cur_las) || npoints(cur_las) == 0) next
    
    success_index <- success_index + 1
    results[[success_index]] <- data.frame(height = cur_las@data$Z,
                                           return_num = cur_las@data$ReturnNumber)
    message('Done ', success_index, ' of ', n_success)
  }
    sum_df <- do.call(rbind, results)
    sum_df$return_num <- factor(sum_df$return_num,
                                levels = sort(unique(sum_df$return_num),
                                              decreasing = TRUE))
    sum_df
}

# identify input dataset paths
paths <- yaml::read_yaml('../paths.yml')

# load in land groups layer
land_groups <- rast(path.expand(paths$land_groups))

# read in las catalog
ctg <- readLAScatalog(path.expand(paths$lidar_catalog))
opt_filter(ctg) <- '-drop_class 7 18'  # drop noise

# sample lidar points from one land group (keeping sample circles within
# land boundaries)
sample_group <- function(groups, values, n) {
  r <- ifel(groups %in% values, 1, NA)
  r <- focal(r, w = 3, fun = 'min', na.rm = FALSE)  # drop cells at edge
  pts <- spatSample(r, size = n, method = 'random', na.rm = TRUE,
                    as.points = TRUE)
  st_as_sf(pts)
}

set.seed(4303)
pohn_main_samp <- circ_samp(sample_group(land_groups, 1, 1000), ctg)
lagoon_islet_samp <- circ_samp(sample_group(land_groups, 2, 1000), ctg)
pohn_all_samp <- circ_samp(sample_group(land_groups, values = c(1, 2),
                                        n = 1000), ctg)

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
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12),
    axis.title.x = element_blank()
  )
main_point_height_plot

# Do the same for the lagoon islets . . .
lagoon_point_height_plot <- ggplot(subset(lagoon_islet_samp, return_num %in% 1:3),
                                   aes(y = height, fill = return_num)) +
  geom_histogram(binwidth = 0.5) +
  labs(x = 'Number of returns', fill = 'Return number', title = 'Lagoon islets') +
  ylim(y = c(-1.5, 42)) +
  xlim(y = c(0, 100000)) +
  guides(fill = 'none') + theme(axis.title.y = element_blank()) +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )
lagoon_point_height_plot

# . . . and for all the surveyed land.
all_point_height_plot <- ggplot(subset(pohn_all_samp, return_num %in% 1:3),
                                aes(y = height, fill = return_num)) +
  geom_histogram(binwidth = 0.5) +
  labs(x = 'Number of returns', fill = 'Return number', title = 'All surveyed land') +
  ylim(y = c(-1.5, 42)) +
  xlim(y = c(0, 100000)) +
  #    guides(fill = 'none') +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  
  #    guides(fill = guide_legend(reverse = TRUE)) + theme(axis.title.y = element_blank()) +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
all_point_height_plot

# Assess pulse and ground point densities.
pulse <- rast(path.expand(paths$pulse_density))
point <- rast(path.expand(paths$point_density))
ground <- rast(path.expand(paths$grnd_density))

# Resample land_groups to match the density layers.
land_1m <- resample(land_groups, pulse)

# Mask the pulse, point, and ground point density layers.
pulse_masked <- mask(pulse, land_1m)
point_masked <- mask(point, land_1m)
ground_masked <- mask(ground, land_1m)

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
aoi_dir <- path.expand(paths$aois)
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

