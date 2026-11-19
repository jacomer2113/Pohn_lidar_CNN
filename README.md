# Pohn_lidar_CNN
This repository holds code relating to the manuscript 'Mapping Landscape
Modifications on Pohnpei with Airborne Lidar and Convolutional Neural Networks.'
It

1. de-noises the lidar point cloud and generates rasterized
(gridded) data products (point_cloud_to_rasters),  
2. builds, trains, and applies a predictive model using a
convolutional neural network (CNN) to the recognition of certain features
perceptible in visualizations of a high-resolution lidar-derived DTM
(CNN_modeling), and  
3. analyzes the probability map resulting from the modeling exercise in
conjunction with other data sets, and generates some of
the figures for the manuscript (analysis_and_figures).  

The first of these directories, point_cloud_to_rasters, contains only
one script, which is relatively straightforward. The second, CNN_modeling,
is more complicated. The modeling process is implemented in four functions:
raster_chunking.R, generate_images.R, train_test_cnn.R, and
run_cnn_predictions.R. For clarity's sake I include a script,
demonstration_script.R, that calls the four functions and demonstrates how
together they constitute the modeling workflow. Last, in analysis_and_figures,
the script called pmap_analysis.R is the one
that performs most of the analysis of the model output; it also generates
some of the manuscript's figures. The script point_cloud_eval.R evaluates the point
cloud data set -- that is, for gathering measurements of pulse density,
ground point density, etc. for different parts of the survey area.
The script mangrove_predictors_for_pub.R
shows how the data necessary for classifying mangrove, non-mangrove land, and
water is assembled; it generates a raster stack at 10-m resolution of
Sentinel-2 spectral reflectance data, spectral indices, and elevation and
vegetation structure metrics from the lidar point cloud. The markdown document
mangrove_modeling_markdown.md demonstrates the mangrove/non-mangrove land/water
modeling process itself.

Data inputs are specified by reference to a .yml file. For the purposes of
demonstration the contents of this .yml file are simulated examples.
