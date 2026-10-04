# Helper functions for tile-wise prediction of very large rasters.
# Splitting a raster into tiles keeps memory use flat when a fitted model has to be
# applied to tens of millions of cells.

#' Split a raster (or raster stack) into a grid of tiles
#' @param r      a terra SpatRaster
#' @param n_row  number of tile rows
#' @param n_col  number of tile columns
#' @return list of SpatRaster tiles (cropped to their extent)
make_tiles <- function(r, n_row = 4, n_col = 4) {
  e  <- terra::ext(r)
  xs <- seq(e$xmin, e$xmax, length.out = n_col + 1)
  ys <- seq(e$ymin, e$ymax, length.out = n_row + 1)
  tiles <- list()
  for (i in seq_len(n_row)) {
    for (j in seq_len(n_col)) {
      te <- terra::ext(xs[j], xs[j + 1], ys[i], ys[i + 1])
      tiles[[length(tiles) + 1]] <- terra::crop(r, te)
    }
  }
  tiles
}

#' Apply a fitted sdm ensemble tile by tile and mosaic the result
#' @param model   fitted sdm object
#' @param predictors terra SpatRaster stack of predictors
#' @param n_row,n_col tile grid
#' @param weight  ensemble statistic used for weighting (default "AUC")
predict_ensemble_tiled <- function(model, predictors, n_row = 4, n_col = 4, weight = "AUC") {
  tiles <- make_tiles(predictors, n_row, n_col)
  out <- vector("list", length(tiles))
  for (k in seq_along(tiles)) {
    out[[k]] <- tryCatch({
      sdm::ensemble(model, raster::stack(tiles[[k]]),
                    setting = list(method = "weighted", stat = weight))
    }, error = function(err) {
      message("Tile ", k, " failed: ", conditionMessage(err))
      NULL
    })
  }
  out <- Filter(Negate(is.null), out)
  do.call(terra::mosaic, lapply(out, terra::rast))
}
