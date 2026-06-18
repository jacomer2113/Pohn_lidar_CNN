require(terra)

# Identify directory where S-2 download directories are stored.
sent_dir <- file.path('~/Sentinel-2/')

# Define a function to conveniently read in, from a collection of Sentinel-2
# downloads, all the layers of a particular type, mask them using the SCL 
# (Scene Classification Layer), crop them, and stack them.
grab_s2_stack <- function(sentinel_dir, identifying_string, resolution = 'R10m',
                       extent = NULL) {
  # Reads in, from a collection of Sentinel-2 downloads, all of the layers
  # of a specified type; masks them using the Scene Classification Layer (SCL),
  # keeping only those pixels classified as vegetated, non-vegetated, water, 
  # or unclassified; and stacks the resulting layers.
  #
  # Args:
  #   sentinel_dir: A directory holding all the Sentinel-2 download directories
  #      to be worked with.
  #   identifying_string: A short identifier for the particular layer type to
  #      be stacked. Obviously, this needs to constitute part of the layer
  #      type's file name -- for example, '_B08_' for band 8.
  #   resolution: The resolution of the layer type being read in. Must be one
  #      of 'R10m', 'R20m', or 'R60m.'
  #   extent: The SpatExtent to which the layers being read in should be
  #      cropped.
  #
  # Returns:
  #   A stack of Sentinel-2 layers, of the kind identified with
  #   identifying_string, masked of all SCL classifications other than 4
  #   (vegetation), 5 (not vegetated), 6 (water), and 7 (unclassified).
  #
  if ((resolution %in% c('R10m', 'R20m', 'R60m')) == FALSE) {
    stop('Resolution must be R10m, R20m, or R60m.')
  }
  index <- 1
  for (i in dir(sentinel_dir)) {
    gran_dir <- dir(paste0(sentinel_dir, i, '/GRANULE'))
    layer_file_name <- grep(identifying_string,
                            dir(paste0(sentinel_dir, i, '/GRANULE/', gran_dir,
                                       '/IMG_DATA/', resolution, '/')),
                            value = TRUE)
    layer_file_path <- paste0(sentinel_dir, i, '/GRANULE/', gran_dir,
                              '/IMG_DATA/', resolution, '/', layer_file_name)
    scl_dir <- paste0(sentinel_dir, i, '/GRANULE/', gran_dir, '/IMG_DATA/R20m/')
    scl_file_name <- grep('_SCL_', dir(scl_dir), value = TRUE)
    
    if (is.null(extent)) {
      scl <- rast(paste0(scl_dir, scl_file_name))
      stack <- rast(layer_file_path)
    } else {
      scl <- crop(rast(paste0(scl_dir, scl_file_name)), extent)
      stack <- crop(rast(layer_file_path), extent)
    }
    if (sum(res(stack) < res(scl)) > 0) {
      scl <- disagg(scl, fact = res(scl) / res(stack))
    }
    int_ext <- intersect(ext(scl), ext(stack))
    scl <- crop(scl, int_ext)
    stack <- crop(stack, int_ext)
    cloud_etc_mask <- scl == 4 | scl == 5 | scl == 6 | scl == 7
    stack_masked <- mask(stack, cloud_etc_mask, maskvalues = FALSE)
    
    plot(stack_masked)
    
    if (index == 1) {
      return_stack <- stack_masked
    } else {
      return_stack <- c(return_stack, stack_masked)
    }
    
    print(paste0('Done reading in ', index, ' of ', length(dir(sent_dir))))
    index <- index + 1
  }
  return(return_stack)
}

# Define a function to remove, from a SpatRaster, the Inf values and
# outliers (more than three standard deviations from mean).
clean_layer <- function(layer) {
  # Remove Inf values and all values more than 3 standard deviations from the
  # mean from a SpatRaster.
  #
  # Args:
  #   layer: The input SpatRaster.
  # 
  
  # Remove Inf values.
  layer_finite <- mask(layer, is.finite(layer), maskvalues = FALSE)
  
  # Remove outliers.
  mu <- mean(values(layer_finite), na.rm = TRUE)
  std <- sd(values(layer_finite), na.rm = TRUE)
  outlier_mask <- layer_finite >= (mu - (3 * std)) & layer_finite <= (mu + (3 * std))
  layer_cleaned <- mask(layer_finite, outlier_mask, maskvalues = FALSE)
  return(layer_cleaned)
}

# Identify the extent (in UTM WGS84 57N coordinates) of Pohnpei -- I'm taking
# these from the extent of the DTM derived from the lidar.
pohnpei_extent <- ext(c(401831, 429832, 749190, 775158))

# Chen et al. 2023 offer a Red-Edge Mangrove Index (REMI) calculated as
# REMI = (red edge - red) / (SWIR1 - green)
# which, for Sentinel-2, corresponds to bands
# (6 - 4) / (11 - 3).
# Baloloy et al. (2020) offer a Mangrove Vegetation Index,
# MVI = (NIR - green) / (SWIR1 - green), or
# (8 - 3) / (11 - 3).
# We'll use the function we defined above to create cloud-free composite
# images of Pohnpei for all the bands involved in these.
redge <- median(grab_s2_stack(sent_dir, '_B06_', resolution = 'R20m',
                           extent = pohnpei_extent), na.rm = TRUE)
red <- median(grab_s2_stack(sent_dir, '_B04_', extent = pohnpei_extent),
              na.rm = TRUE)
swir1 <- median(grab_s2_stack(sent_dir, '_B11_', resolution = 'R20m',
                           extent = pohnpei_extent), na.rm = TRUE)
green <- median(grab_s2_stack(sent_dir, '_B03_', extent = pohnpei_extent),
                na.rm = TRUE)
nir <- median(grab_s2_stack(sent_dir, '_B08_', extent = pohnpei_extent),
              na.rm = TRUE)

# Grab blue, too, in case we want to make an RGB image.
blue <- median(grab_s2_stack(sent_dir, '_B02_', extent = pohnpei_extent),
               na.rm = TRUE)

# Since the red edge and SWIR bands are 20-m resolution, but the red, blue, and
# green bands are 10-m, the extents of the former are different from those of
# the latter -- we have to fix that.
ext(redge) == ext(swir1)
ext(red) == ext(green)
ext(redge) == ext(red)

# Disaggregate the lower-resolution layers to match the higher-resolution
# layers.
redge_disagg <- disagg(redge, fact = res(redge) / res(red))
swir1_disagg <- disagg(swir1, fact = res(swir1) / res(red))

# Find the extent at which all 5 layers intersect.
band_int <- intersect(ext(redge), ext(red))

# Crop all 6 layers to that extent and stack them.
redge_disagg <- crop(redge_disagg, band_int)
red <- crop(red, band_int)
swir1_disagg <- crop(swir1_disagg, band_int)
green <- crop(green, band_int)
nir <- crop(nir, band_int)
blue <- crop(blue, band_int)
band_stack <- c(redge_disagg, red, swir1_disagg, green, nir, blue)
names(band_stack) <- c('redge', 'red', 'swir1', 'green', 'nir', 'blue')

# Generate REMI layer.
remi_layer <- (band_stack$redge - band_stack$red) /
  (band_stack$swir1 - band_stack$green)

# Remove Inf and outlying values.
remi_cleaned <- clean_layer(remi_layer)

# Generate MVI layer.
mvi_layer <- (band_stack$nir - band_stack$green) / 
  (band_stack$swir1 - band_stack$green)

# Remove Inf and outlying values.
mvi_cleaned <- clean_layer(mvi_layer)

# Baloloy et al. (2020) identify 4.5 as the minimum MVI threshold for 
# mangroves (see page 103). Let's just accept that threshold and apply it
# here. It's obviously not perfect, because there are areas above that 
# threshold far inland -- but if we apply this layer as a simple mask at lower
# elevations, it should be good enough for our purposes.
mangroves_only <- mask(mvi_cleaned, mvi_cleaned >= 4.5, maskvalues = FALSE)
mangrove_mask <- !is.na(mangroves_only)

# Make a water mask with the SCL layer of the Sentinel-2 data.
scl <- grab_s2_stack(sent_dir, '_SCL_', resolution = 'R20m',
                   extent = pohnpei_extent)

for (i in 1:nlyr(scl)) {
  scl[[i]] <- scl[[i]] == 6  # 6 is the value for water in the SCL layer.
}
water_cumulative <- sum(scl, na.rm = TRUE)
water_mask_scl <- water_cumulative != 0

# To give ourselves options, make another water mask, this one from
# the NIR band of Sentinel-2. I chose the reflection
# value of 3000 based only on a visual evaluation of the NIR value histogram, 
# so this mask should be considered approximate.
nir_water_mask <- band_stack$nir < 3000

# In case we want a mask for water only, set the NA values in this mask equal
# to FALSE, to avoid data gaps.
nir_water_mask[is.na(nir_water_mask)] <- FALSE

# Combine mangrove and water masks. First, we need to disaggregate and crop
# the SCL water mask so that its extent will match that of the mangrove mask.
water_mask_scl_disagg <- disagg(water_mask_scl, fact = res(water_mask_scl) /
                                res(mangrove_mask))
water_mask_scl_cropped <- crop(water_mask_scl_disagg, mangrove_mask)

# These are our two versions of the combined masks, one using NIR to mask
# water and the other using the SCL layer to mask water.
scl_combined_mask <- mangrove_mask + water_mask_scl_cropped
scl_combined_mask <- scl_combined_mask != 0

nir_combined_mask <- mangrove_mask + nir_water_mask
nir_combined_mask <- nir_combined_mask != 0

# It's not obvious that either of these masks is better than the other. I'll
# take the mask that uses NIR to identify water, smooth it with a modal
# filter, and use that one.
nir_comb_mask_filt <- focal(nir_combined_mask, w = 15, fun = 'modal')

# The last step created some gaps of NA values -- I'll set those equal to 1.
nir_comb_mask_filt <- subst(nir_comb_mask_filt, NA, 1)

# Load in our 1-meter-resolution, lidar-derived DTM to exclude all areas
# above a certain elevation from the mask, which is incorrectly applied
# to some inland areas.
dtm <- rast(paste0('~/Documents/CSRM/Further_Pohnpei_work/RRNES_grant/lidar/',
  'jake_DEMs/DTM_denoised_corrected_kriging_1m.tif'))

# Resample the DTM to match the resolution and extent of the masks we've made.
dtm_resamp <- resample(dtm, nir_comb_mask_filt, method = 'bilinear')
dtm_resamp <- subst(dtm_resamp, NA, min(values(dtm_resamp), na.rm = TRUE))

# The DTM is in ellipsoid heights, which is why our cutoff looks strangely
# high.
elev_cutoff <- mask(dtm_resamp, dtm_resamp > 75, maskvalues = FALSE)
elev_cutoff_bin <- !is.na(elev_cutoff)
elev_cutoff_filt <- focal(elev_cutoff_bin, w = 25, fun = 'modal')

# Combine the mangrove and water mask with the elevation cutoff.
masks_plus_elev <- nir_comb_mask_filt - elev_cutoff_filt
total_mask <- subst(masks_plus_elev, -1, 0)

total_mask_path <- paste0('~/Documents/CSRM/Further_Pohnpei_work/RRNES_grant/',
  'Publications/lidar_methods_paper/GIS/masks/total_mask_10m.tif')

if (!file.exists(total_mask_path)) {
  writeRaster(total_mask, total_mask_path)
}

# Do the same for the water mask.
nir_water_elev <- nir_water_mask - elev_cutoff_filt
water_mask <- subst(nir_water_elev, -1, 0)

water_mask_path <- file.path('nir_water_mask_10m.tif')

if (!file.exists(water_mask_path)) {
  writeRaster(water_mask, water_mask_path)
}


# And, while we're at it, here's our RGB image, without clouds.
rgb_s2 <- c(band_stack$red, band_stack$green, band_stack$blue)

# However, that RGB image doesn't look as good as the TCI product that 
# Sentinel-2 provides does, so we'll go with the TCI.
pohnpei_larger_extent <- ext(c(399831, 431832, 747190, 777158))
s2_image <- grab_s2_stack(sent_dir, '_TCI_', resolution = 'R10m',
                          extent = pohnpei_extent)
reds <- s2_image[[seq(1, 28, 3)]]
greens <- s2_image[[seq(2, 29, 3)]]
blues <- s2_image[[seq(3, 30, 3)]]
s2_mean_rgb <- c(mean(reds, na.rm = TRUE),
                 mean(greens, na.rm = TRUE),
                 mean(blues, na.rm = TRUE))
s2_median_rgb <- c(median(reds, na.rm = TRUE),
                 median(greens, na.rm = TRUE),
                 median(blues, na.rm = TRUE))



