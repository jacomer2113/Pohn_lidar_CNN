raster_chunking <- function(raster, output_map_res) {
  # Generates a grid system of contiguous, square polygons, of a size defined
  # by output_map_res, that cover the input raster, except in areas that are NA.
  #
  # Args:
  #   raster: The SpatRaster for which the grid system will be generated. The 
  #     coordinate reference system must be a Universal Transverse Mercator
  #     (UTM) system, using meters as the linear unit.
  #   output_map_res: The desired size of each square in the resulting grid
  #     system, in meters.
  #
  # Returns:
  #   An sf object of type POLYGON consisting of contiguous, square polygons
  #     with side lengths equal to output_map_res.
  terra_installed <- require(terra)
  sf_installed <- require(sf)
  sfh_installed <- require(sfheaders)
  stopifnot(terra_installed & sf_installed & sfh_installed)
  cell_size <- res(raster)[1]
  if (cell_size != res(raster)[2]) stop('Input raster must have square cells.')
  if (cell_size > output_map_res) stop('Desired chunk size must be at least that of input raster cells.')
  if (cell_size < output_map_res) {
    agg_factor <- output_map_res / cell_size
    if (output_map_res %% cell_size != 0) {
      stop('Output map cell size must be a multiple of input raster cell size.')
    }
    print('Aggregating input raster to match desired chunk size.')
    template_raster <- aggregate(raster, agg_factor)
  }
  not_na_cells <- cells(template_raster)  # For all non-NA cells in template_raster . . .
  chunk_centers <- xyFromCell(template_raster, not_na_cells)  # . . . find the center coordinates.
  half_chunk_side <- output_map_res / 2
  make_chunks <- function(row_number, input_matrix) {
    left <- input_matrix[row_number, 1] - half_chunk_side
    right <- input_matrix[row_number, 1] + half_chunk_side
    top <- input_matrix[row_number, 2] + half_chunk_side
    bottom <- input_matrix[row_number, 2] - half_chunk_side
    xs <- c(left, right, right, left, left)
    ys <- c(top, top, bottom, bottom, top)
    id <- row_number
    cbind(xs, ys, id)
  }
  coordinate_matrices <- lapply(1:nrow(chunk_centers),
                                make_chunks,
                                input_matrix = chunk_centers)
  coordinate_df <- as.data.frame(do.call(rbind, coordinate_matrices))
  chunk_sf <- sfheaders::sf_polygon(coordinate_df,
                                    x = 'xs',
                                    y = 'ys',
                                    polygon_id = 'id')
  chunk_sf$label <- NA  # Make fields for the sf data frame to hold labels,
  chunk_sf$training <- NA  # training data flags, and
  chunk_sf$testing <- NA  # testing data flags.
  st_crs(chunk_sf) <- st_crs(raster)
  return(chunk_sf)
}
