# Load in required packages.
require(lidR)
require(terra)

# Define a function that can be used, with catalog_apply, to generate pulse
# and point density maps for the entire point cloud.
pulse_map <- function(chunk) {
  cur.las <- readLAS(chunk)
  if (lidR::is.empty(cur.las)) return(NULL)
  cur.las.pids <- retrieve_pulses(cur.las)
  cur.pulse.map <- rasterize_density(cur.las.pids, res = 1)
  return(cur.pulse.map)
}

# Indicate the catalog of las files (the directory that holds them all).
cat.dir <- file.path('drive', 'user', 'catalog_name')
cat <- readLAScatalog(cat.dir)
las_check(cat)  # Inspect the las catalog.
st_crs(cat) <- 32657  # Make sure the projection is identified (UTM, zone 57N).

# Specify a temporary directory to which output files will be written
# (in order to avoid exceeding memory limit).
opt_output_files(cat) <- paste0(tempdir(), '/{*}_noise_classified')

# Classify outlying points as noise. This is done by calculating, for each 
# point, the distance to the k nearest neighbors. Points farther away than
# the mean of these distances plus m times the standard deviation of those
# distances are considered outliers.
cat_noise_classified <- classify_noise(cat, algorithm = sor(k = 10, m = 5))
print('Noise classification done')

# Re-identify the las catalog, this time pointing to the noise-classified
# catalog.
cat <- readLAScatalog(tempdir())

# Make a directory where the las files that include the pulse IDs will
# be stored.
pulse.ID.las.dir <- '~/las_ellip_ground_density'
dir.create(pulse.ID.las.dir)
opt_output_files(cat) <- paste0(pulse.ID.las.dir, '/{*}_w_pID')
options <- list(automerge = TRUE) #  Merge all the outputs of catalog_apply
cat_w_pulse_ID <- catalog_apply(cat, pulse_map, .options = options)

# De-noise the las files in the catalog one by one. (There may be a more
# efficient way to do this than by looping through the las files,
# but if so, I don't know what it is.) Concurrently, for each of the las
# files, retrieve the pulse ID of every point. This way, we'll be able to
# make a pulse density map -- as opposed to a point density map -- later.
# First, get the names of all the noise-classified las files.
dn.las.files <- grep('.las', dir(tempdir()), value = TRUE)

# Specify a directory to hold the las files from which the noise-classified
# points have been removed. If that directory doesn't exist already, create
# it.
denoised.dir <- 'C:/Users/idumi/Desktop/denoised_las_ellip'
if (file.exists(denoised.dir) == FALSE) {
  dir.create(denoised.dir)
}

print('Begin re-saving las files with noise-classified points removed')
for(i in 1:length(dn.las.files)) {  # This loop throws an error where i = 
# 352, because that tile has only 11 points in it, and they're all noise (it's
# a small tile over the water). So:
  if (i == 352) {
    print('skipping file number 352, because it contains only noise')
  } else {
    cur.name <- dn.las.files[i]
    cur.las <- readLAS(paste0(tempdir(), '/', cur.name))
    cur.las.dn <- filter_poi(cur.las, Classification != 7 &
                               Classification != 18) # 7 and 18 are the values
    # to which noise is classified in las files (7 = noise, 18 = high noise).
    cur.las.dn.pID <- retrieve_pulses(cur.las.dn)  # For each point, retrieve
    # the pulse ID and assign it to that point.
    new.name <- gsub('_noise_classified', '_denoised', cur.name)
    writeLAS(cur.las.dn.pID, paste0(denoised.dir, '/', new.name))
    print(paste0('File ', i, ' of ', length(dn.las.files),
                 ' denoised and saved'))
  }
}

# Re-identify the las catalog, this time pointing to the catalog from which
# the noise-classified points have been removed. From here, since cat refers
# to the catalog from which the noise-classified
# points have been removed, the rasters we generate will be created from the
# denoised point cloud.
cat <- readLAScatalog(denoised.dir)

# What's the desired resolution of the output rasters?
des_res <- 1

# Generate a DTM (with kriging). For 1-m resolution, this takes a 
# while (like, a couple of days on my laptop).
dtm_kriging <- rasterize_terrain(cat, res = des_res, algorithm = kriging())
print('DTM generated')

# Generate a DSM (with points-to-raster).
dsm_p2r <- rasterize_canopy(cat, res = des_res, algorithm = p2r())
print('DSM generated')

dens <- rasterize_density(cat, res = des_res) # Generate a lidar point density
# map.
print('Point density map generated')

# Generate a ground point density map. To start, use
# opt_filter() to filter out all points except those classified as ground.
# This will make the process of generating the raster faster.
opt_filter(cat) <- '-keep_class 2'
ground_dens <- rasterize_density(cat, res = des_res)
print('Ground point density map generated')

# Resample the DSM and the point density maps to match the DTM.
dsm_p2r_resamp <- terra::resample(dsm_p2r, dtm_kriging, method = 'bilinear')
dens_resamp <- terra::resample(dens, dtm_kriging, method = 'bilinear')
grnd_dens_resamp <- terra::resample(ground_dens, dtm_kriging, method = 'bilinear')
