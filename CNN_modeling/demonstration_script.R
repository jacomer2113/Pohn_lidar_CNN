# Load required packages.
require(terra)
require(sf)
require(sfheaders)
require(png)
require(magrittr)
require(keras)
require(tfhub)
require(pins)
require(tfdatasets)
require(pROC)

# Source the necessary functions.
source('./raster_chunking.R')
source('./generate_images.R')
source('./train_test_cnn.R')
source('./run_cnn_predictions.R')

# Identify input data paths.
paths <- yaml::read_yaml('../paths.yml')

chunks_to_label <- raster_chunking(raster = rast(path.expand(paths$slrm)),
                                   output_map_res = 10)

# Write chunks_to_label as a shapefile.
st_write(chunks_to_label, path.expand(paths$unlabeled_chunks))

# The following step requires GIS software. shapefile_filename is opened in
# a GIS, like QGIS, and some of the square polygons, or chunks, are labeled. 
# Chunks that overlap with a feature of interest get 'label' 1. Chunks that do
# not overlap with a feature of interest get 'label' 0. Unlabeled chunks -- that
# is, the chunks that are going to be run through the model after it's
# trained -- retain NA as the 'label' value. Of the labeled chunks, the subset
# that will be used as training data gets 'training' value 1, NA otherwise.
# The subset that will be used as testing data gets 'testing' value 1, NA
# otherwise. Once this labeling process is complete, we load the chunks
# shapefile back in:

labeled_chunks <- st_read(path.expand(paths$labeled_chunks))

testing_data_labels <- generate_images(labeled_chunks,
                                       raster_layer = rast(
                                         path.expand(paths$slrm)),
                                       image_size = 75,
                                       image_output_directory = path.expand(
                                         paths$cnn_image_dir
                                       ))

model_performance <- train_test_cnn(testing_data_labels,
                                    image_input_directory = path.expand(
                                      paths$cnn_image_dir),
                                    architecture = 'VGG19',
                                    model_output_filename = 'model_name',
                                    validation_proportion = 0.1,
                                    n_epochs = 75, 
                                    fine_tune_epochs = NULL)

print(model_performance)

output_prob_map <- run_cnn_predictions(labeled_chunks,
                                       raster_layer =
                                         rast(path.expand(paths$slrm)),
                                       model =
                                         keras::load_model_tf('model_name'))

# Note that output_prob_map will be NA in all areas that were not run 
# through the model, which includes the areas that were labeled to serve as
# training and testing data -- but these blank spaces in output_prob_map can,
# of course, be populated with their assigned labels.

