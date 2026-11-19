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
is more complicated, for which reason I include a specific script,
demonstration_script.R, to show how the four functions in that directory
lead from one to the next and, together, constitute the modeling workflow.
Last, in analysis_and_figures, the script called pmap_analysis.R is the one
that performs most of the analysis of the model output; it also generates
some of the manuscript's figures. The script S2_mask.R generates,
from a stack of Sentinel-2 imagery, masks for the mangroves and water
around Pohnpei. The script point_cloud_eval.R evaluates the point
cloud data set -- that is, for gathering measurements of pulse density,
ground point density, etc. for different parts of the survey area
(see Table 1 of the manuscript).

Data inputs are specified by reference to a .yml file. For the purposes of
demonstration this .yml file is a simulated example.
