require(terra)
require(lidR)
require(ranger)

# get file paths
paths <- yaml::read_yaml('../paths.yml')

# assemble Sentinel-2 data
sent_dir <- path.expand(paths$S2_directory)  # Sentinel-2 download directory

# define function to extract specific band data from S2 download
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

# identify the extent (in UTM WGS84 57N coordinates) of Pohnpei
pohnpei_extent <- ext(c(401831, 429832, 749190, 775158))

# grab necessary bands
blue <- median(grab_s2_stack(sent_dir, '_B02_', extent = pohnpei_extent),
               na.rm = TRUE)
green <- median(grab_s2_stack(sent_dir, '_B03_', extent = pohnpei_extent),
                na.rm = TRUE)
red <- median(grab_s2_stack(sent_dir, '_B04_', extent = pohnpei_extent),
              na.rm = TRUE)
redge <- median(grab_s2_stack(sent_dir, '_B06_', resolution = 'R20m',
                              extent = pohnpei_extent), na.rm = TRUE)
nir <- median(grab_s2_stack(sent_dir, '_B08_', extent = pohnpei_extent),
              na.rm = TRUE)
swir1 <- median(grab_s2_stack(sent_dir, '_B11_', resolution = 'R20m',
                              extent = pohnpei_extent), na.rm = TRUE)
swir2 <- median(grab_s2_stack(sent_dir, '_B12_', resolution = 'R20m',
                              extent = pohnpei_extent), na.rm = TRUE)

# disaggregate 20-m pixels to 10-m
to_10m <- function(in_band) {
  band_10m <- disagg(in_band, fact = 2)  # method is 'near' by default
  band_10m <- crop(band_10m, nir)  # crop to nir extent . . . 
  extend(band_10m, nir)  # . . . and extend again to ensure extents match
}

redge_10m <- to_10m(redge)
swir1_10m <- to_10m(swir1)
swir2_10m <- to_10m(swir2)

s2_stack <- c(blue, green, red, redge_10m, nir, swir1_10m, swir2_10m)
names(s2_stack) <- c('blue', 'green', 'red', 'redge', 'nir', 'swir1', 'swir2')

# calculate spectral indices
ndvi <- (nir - red) / (nir + red)
ndmi <- (nir - s2_stack$swir1) / (nir + s2_stack$swir2)
mvi <- (nir - green) / (s2_stack$swir1 - green)
mvi <- ifel(is.finite(mvi), mvi, NA)  # set Inf values in MVI to NA
indices <- c(ndvi, ndmi, mvi)
names(indices) <- c('ndvi', 'ndmi', 'mvi')

# read in las catalog
ctg <- readLAScatalog(path.expand(paths$lidar_catalog))
opt_select(ctg) <- 'xyzirc'  # keep only necessary features
opt_filter(ctg) <- '-drop_class 7 18'  # drop noise-classified points

# function to pass to pixel_metrics, gathering necessary point cloud variables
cell_metrics <- function(z, cls, rn, int) {
  first <- rn == 1L  # first return?
  ground <- cls == 2L  # ground return?
  list(
    n_pts       = as.numeric(length(z)),
    h_p95       = as.numeric(quantile(z, 0.95, names = FALSE)),
    h_sd        = as.numeric(sd(z)),
    cover_2m    = mean(z[first] > 2),  # share of first returns above 2 m, as
    # measure of canopy cover
    low_frac    = mean(z > 0.5 & z <= 1.5),  # share of all returns between
    # 0.5 and 1.5 meters, as measure of understory
    ground_frac = mean(ground),
    ground_int  = as.numeric(median(int[ground]))
  )
}

# grab abovementioned variables across grid defined by S2 10-m
lidar_metrics <- pixel_metrics(
  ctg,
  ~cell_metrics(Z, Classification, ReturnNumber, Intensity),
  res = rast(nir)   # 10 m Sentinel-2 grid
)

# ensure that the lidar metrics grid matches the S2 grid
lidar_metrics <- extend(lidar_metrics, nir)
lidar_metrics <- crop(lidar_metrics, nir)
compareGeom(lidar_metrics, nir)   # errors if the grids don't match

# read in DTM and resample to match 10-m grid
dtm <- rast(path.expand(paths$dtm))
dtm_10m <- resample(dtm, nir, method = 'average')

compareGeom(s2_stack, dtm_10m)
compareGeom(s2_stack, lidar_metrics)

# set aside data gap in center of island and all areas beyond survey region
cutoff <- 175
zone <- ifel(dtm_10m < cutoff, 1, NA)   # NA above 175 m and outside the DTM
has_dtm <- !is.na(dtm_10m)

# n_pts to 0 where lidar was flown but no points retrieved
n_pts <- ifel(is.na(lidar_metrics$n_pts) & has_dtm, 0, lidar_metrics$n_pts)

# all other flown areas where no points were recovered are set to -99
others <- setdiff(names(lidar_metrics), 'n_pts')
others_filled <- rast(lapply(others, function(n) {
  ifel(is.na(lidar_metrics[[n]]) & has_dtm, -99, lidar_metrics[[n]])
}))

lidar_filled <- c(n_pts, others_filled)
names(lidar_filled) <- c('n_pts', others)

elev <- dtm_10m
names(elev) <- 'elev'

# assemble predictor layers
predictors <- mask(c(s2_stack, indices, lidar_filled, elev), zone)

