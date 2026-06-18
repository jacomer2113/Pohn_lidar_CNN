generate_images <- function(labeled_chunks, raster_layer, image_size,
                            image_output_directory = getwd()) {
  # Generates images of side length image_size based on the labels assigned,
  # by labeled_chunks, to sub-areas of raster_layer.
  # 
  # Args: 
  #   labeled_chunks: An sf object of type POLYGON, with fields named 'id_num,'
  #     'label,' 'training,' and 'testing,' in addition to the geometry column,
  #     that provides the labels that will be used to train the output model. It
  #     is intended that this sf object is the sf object returned by the
  #     function raster_chunking and to which labels have thereafter been
  #     applied.
  #   raster_layer: The SpatRaster to which the CNN model will be applied.
  #     labeled_chunks consists of contiguous square polygons that together 
  #     cover this SpatRaster and identify square sub-areas that are positive
  #     (spatially overlap with features of interest), negative (do not
  #     overlap with features of interest), or unlabeled.
  #   image_size: The side length, in meters of each square image that will be
  #     considered by the CNN model. The image size must be greater than the
  #     size of the polygons in labeled_chunks. From each square chunk,
  #     a corresponding, larger image, centered on the same two-dimensional
  #     point, will be generated. It is these images, rather than
  #     images corresponding only to the chunks themselves, that will
  #     be considered by the model.
  #   image_output_directory: A string identifying the directory to which the
  #     training and testing images will be written. This is the current 
  #     working directory by default. In this directory, two sub-directories
  #     will be created: training_dir and testing_dir. training_dir, in turn,
  #     has two sub-directories: class_neg and class_pos. The directory
  #     structure of training_dir makes it easy for
  #     keras::image_dataset_from_directory to read in the image data and
  #     assign the appropriate labels. As for the labels for the testing
  #     images, these will be contained in a vector that the function will
  #     return.
  #
  # Returns:
  #   A vector of the labels of the testing images (also, the training and
  #   testing images are written to disk in image_output_directory).
  
  # Make sure required packages are installed.
  terra_installed <- require(terra)
  sf_installed <- require(sf)
  png_installed <- require(png)
  stopifnot(terra_installed & sf_installed & png_installed)
  
  # Make sure the input labeled_chunks has the necessary fields.
  if ('id' %in% colnames(labeled_chunks)) {
    colnames(labeled_chunks)[colnames(labeled_chunks) == 'id'] <- 'id_num'
  }
  id_present <- 'id_num' %in% colnames(labeled_chunks)
  label_present <- 'label' %in% colnames(labeled_chunks)
  train_present <- 'training' %in% colnames(labeled_chunks)
  test_present <- 'testing' %in% colnames(labeled_chunks)
  geometry_present <- 'geometry' %in% colnames(labeled_chunks)
  if ((id_present & label_present & train_present &
       test_present & geometry_present) == FALSE) {
    stop(paste0('labeled_chunks needs \'id_num\', \'label\', \'training\',',
                '\'testing\', and \'geometry\' columns'))
  }
  
  # Create a bounding polygon for each image to be created.
  centroids <- st_centroid(labeled_chunks)  # Find centroid of each chunk.
  cent_matrix <- cbind(centroids$id, st_coordinates(centroids))  # Put centroid
  # coordinates, and an ID number for each, in matrix form.
  colnames(cent_matrix) <- c('id_num', 'X', 'Y')
  half_image_side <- image_size / 2
  
  # Define a function to create an image boundary around each centroid.
  make_image_bounds <- function(id_num, input_matrix) {
    left <- input_matrix[id_num, 2] - half_image_side
    right <- input_matrix[id_num, 2] + half_image_side
    top <- input_matrix[id_num, 3] + half_image_side
    bottom <- input_matrix[id_num, 3] - half_image_side
    xs <- c(left, right, right, left, left)
    ys <- c(top, top, bottom, bottom, top)
    cbind(id_num, xs, ys)
  }
  
  # Apply that function to all centroids.
  print('Creating boundary boxes for the images to be created.')
  image_bound_coords <- lapply(cent_matrix[, 1],
                               make_image_bounds,
                               input_matrix = cent_matrix)
  
  # Put those results into a data frame.
  bound_coord_df <- as.data.frame(do.call(rbind, image_bound_coords))
  
  # Convert the data frame to an sf data frame.
  print('Putting those boundary boxes into an sf data frame.')
  image_bound_sf <- sfheaders::sf_polygon(bound_coord_df,
                                          x = 'xs',
                                          y = 'ys',
                                          polygon_id = 'id_num')
  image_bound_sf$label <- labeled_chunks$label
  image_bound_sf$training <- labeled_chunks$training
  image_bound_sf$testing <- labeled_chunks$testing
  
  # Make sure the CRS matches that of the input shapefile.
  st_crs(image_bound_sf) <- st_crs(labeled_chunks)
  
  # Create directories to hold the image files that will be used to train and
  # test the model.
  training_dir <- paste0(image_output_directory, '/training_dir')
  testing_dir <- paste0(image_output_directory, '/testing_dir')
  train_neg_dir <- paste0(image_output_directory, '/training_dir/class_neg')
  train_pos_dir <- paste0(image_output_directory, '/training_dir/class_pos')
  dir.create(training_dir)
  dir.create(testing_dir)
  dir.create(train_neg_dir)
  dir.create(train_pos_dir)
  
  # Identify the labeled data flagged as testing data and grab the chunk numbers.
  # Testing data observations should be flagged with 1; all other observations
  # should have testing value 0. The training data should be identified in the
  # same manner, but with the training variable.
  testing_obs <- labeled_chunks[which(image_bound_sf$testing == 1), ]
  testing_indices <- testing_obs$id_num
  
  # Loop through the chunks, writing each corresponding section of the input
  # raster as a png image in the appropriate directory. NOTE: There's 
  # probably a faster way to do this than with a loop, but I haven't gotten
  # around to coding it.
  pos.train.ind <- 1
  neg.train.ind <- 1
  pos.test.ind <- 1
  neg.test.ind <- 1
  pos.ind <- 1
  neg.ind <- 1
  total.ind <- 1
  labeled_chunk_ids <- which(!is.na(image_bound_sf$label))
  total_num_labels <- length(labeled_chunk_ids)
  image_id <- sprintf(paste0('%0',
                             nchar(total_num_labels),
                             'd'),
                      1:total_num_labels)
  test_label_vec <- integer()
  for (i in labeled_chunk_ids) {
    lab_val <- image_bound_sf$label[i]
    cur_num <- image_bound_sf$id_num[i]
    cur.rast <- terra::crop(raster_layer, image_bound_sf$geometry[i])
    if (lab_val == 1) {
      if (cur_num %in% testing_indices) {
        out_dir <- testing_dir
        pos.test.ind <- pos.test.ind + 1
        pos.ind <- pos.test.ind
        test_label_vec <- append(test_label_vec, lab_val)
      } else {
        out_dir <- train_pos_dir
        pos.train.ind <- pos.train.ind + 1
        pos.ind <- pos.train.ind
      }
      png(filename = paste0(out_dir, '/', image_id[total.ind], '_chunk_', cur_num,
                            '_pos_image_', pos.ind - 1, '.png'),
          width = image_size, height = image_size)
      if (nlyr(raster_layer) == 1) {
        plot(cur.rast / 255, col = gray.colors(255), mar = c(0, 0, 0, 0),
             legend = FALSE, frame.plot = FALSE, axes = FALSE)
      } else if (nlyr(raster_layer) %in% c(3, 4)) {
        plotRGB(cur.rast)
      }
      dev.off()
    } else {
      if (cur_num %in% testing_indices) {
        out_dir <- testing_dir
        neg.test.ind <- neg.test.ind + 1
        neg.ind <- neg.test.ind
        test_label_vec <- append(test_label_vec, lab_val)
      } else {
        out_dir <- train_neg_dir
        neg.train.ind <- neg.train.ind + 1
        neg.ind <- neg.train.ind
      }
      png(filename = paste0(out_dir, '/', image_id[total.ind],
                            '_chunk_', cur_num,
                            '_neg_image_', neg.ind - 1, '.png'),
          width = image_size, height = image_size)
      if (nlyr(raster_layer) == 1) {
        plot(cur.rast / 255, col = gray.colors(255), mar = c(0, 0, 0, 0),
             legend = FALSE, frame.plot = FALSE, axes = FALSE)
      } else if (nlyr(raster_layer) == 3 | 4) {
        plotRGB(cur.rast)
      }
      dev.off()
    }
    print(paste0('Created ', total.ind, ' of ',
                 length(which(!is.na(image_bound_sf$label))), ' images.'))
    total.ind <- total.ind + 1
  }
  return(test_label_vec)
}
  
  