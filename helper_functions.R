# Load a forestGEO dataset
load_forestGeo <- function(trees_path, #path to raw data 
                           jitter_sd = 0.01, #amount to jitter all points to avoid duplicated points
                           remove_proportion = 0., #do you want to remove a proportion of individuals?
                           grouping_threshold = 50, #how many individuals must be in a species? If a species had less individuals than this number, these are grouped into OTHERSP
                           remove_saplings = TRUE, #remove saplings?
                           species_table, #if remove_saplings = TRUE, must provide a table of species, must have a column of ReproducibleSize or ExtrapolatedReproducibleSize to be able to remove saplings
                           xrange,
                           yrange,
                           alive = TRUE) #only use alive individuals?
{
  
  # Load trees locations
  before <- c(ls(), "before")
  load(trees_path)
  raw_data <- get(setdiff(ls(), before))
  
  # Keep only living/dead trees
  if(alive) {
    symbol <- 'A'
  } else {
    symbol <- 'P'
  }
  raw_data <- raw_data[raw_data$status == symbol, ]
  
  # Convert from factor to strings
  raw_data$sp <- as.character(raw_data$sp)
  
  # Remove rows with missing coordinates or species
  raw_data <- raw_data[!is.na(raw_data$gx) & !is.na(raw_data$gy) & !is.na(raw_data$sp), ]
  
  # Remove saplings? Need reproducible size if we want to do so
  if(remove_saplings & !missing(species_table)) {
    # Best case: we did a bit of extra work and have the extrapolated version
    if("ExtrapolatedReproducibleSize" %in% colnames(species_table)) { #check this column in species_table 
      extrapolated <- sapply(unique(raw_data$sp), function(sp) {
        species_table$ExtrapolatedReproducibleSize[species_table$sp == sp][1]  #for each species in rawdata find extrapolated size from the species_table, and store in extrapolated vector
      })
      # Otherwise, we might have a generic reproducible size
    } else if("ReproducibleSize" %in% colnames(species_table)) { #check this column in species_table 
      extrapolated <- sapply(unique(raw_data$sp), function(sp) { #for each species in rawdata find extrapolated size from the species_table, and store in extrapolated vector
        species_table$ReproducibleSize[species_table$sp == sp][1]
      })
      # If we have neither of these, no way to proceed
    } else {
      extrapolated <- sapply(unique(raw_data$sp), function(sp) NA)
    }
    
    if(!all(is.na(extrapolated))) {
      # Remove individuals with NA dbh
      raw_data <- raw_data[!is.na(raw_data$dbh), ]
      # Only keep individuals where dbh is larger than extrapolated size, therefore removing saplings 
      raw_data <- raw_data[raw_data$dbh > extrapolated[raw_data$sp], ]
    }
  }
  
  # Remove a proportion of individuals
  if(remove_proportion > 0) {
    keep <- ifelse(rbinom(nrow(raw_data), size = 1, prob = 1 - remove_proportion) == 1, TRUE, FALSE)
    raw_data <- raw_data[keep, ]
  }
  
  if(jitter_sd > 0) {
    # Jitter all points to avoid numerical instability
    nr <- nrow(raw_data)
    jit <- jitter_sd / max(xrange[2] - xrange[1], yrange[2] - yrange[1])
    raw_data$gx <- stats::rnorm(nr, mean = raw_data$gx, sd = jit)
    raw_data$gy <- stats::rnorm(nr, mean = raw_data$gy, sd = jit)
    
    # Truncate everything to the observation region (it could have moved outside because of the jitter).
    raw_data$gx[raw_data$gx > xrange[2]] <- xrange[2]
    raw_data$gx[raw_data$gx < xrange[1]] <- xrange[1]
    raw_data$gy[raw_data$gy > yrange[2]] <- yrange[2]
    raw_data$gy[raw_data$gy < yrange[1]] <- yrange[1]
  }
  
  # Group less populated species together dependent on what grouping_threshold is set to
  less_populated <- as.numeric(stats::ave(raw_data$sp, raw_data$sp, FUN = length)) < grouping_threshold
  raw_data$sp[less_populated] <- "OTHSPE"
  
  print(raw_data)
}






# Convert SpatRaster  to 'im' object 
as.im.SpatRaster1 <- function(X) {
  X <- X[[1]]
  rs <- terra::res(X)
  e <- as.vector(terra::ext(X))
  out <- list(
    v = as.matrix(X, wide=TRUE)[nrow(X):1, ],
    dim = dim(X)[1:2],
    xrange = e[1:2],
    yrange = e[3:4],
    xstep = rs[1],
    ystep = rs[2],
    xcol = e[1] + (1:ncol(X)) * rs[1] + 0.5 * rs[1],
    yrow = e[4] - (nrow(X):1) * rs[2] + 0.5 * rs[2],
    type = "real",
    units  = list(singular=units(X), plural=units(X), multiplier=1)
  )
  attr(out$units, "class") <- "unitname"
  attr(out, "class") <- "im"
  out
}


# Convert a matrix to an `im` object
matrix_to_im <- function(mat, xrange, yrange, scale = TRUE) {
  # Coerce to matrix
  mat <- as.matrix(mat)
  
  # Make sure arguments are properly formatted
  if(!is.vector(xrange) | !is.vector(yrange)) {
    stop("Expecting xrange and yrange to be vectors.")
  } else if(!length(xrange) == 2 | !length(yrange) == 2) {
    stop("Expecting xrange and yrange to be of length 2.")
  }
  
  # Should we scale the matrix?
  if(scale) {
    mat <- (mat - mean(mat)) / sd(mat)
  }
  
  # Convert
  spatstat.geom::im(mat, 
                    xrange = xrange, 
                    yrange = yrange)
}


# Convert a covariate to `im` object from its path
covariate_path_to_im <- function(covariate_path, scale = TRUE) {
  extension <- tolower(tools::file_ext(covariate_path))
  if(extension == "rdata") {
    # Load whatever is in the RData
    before <- c(ls(), "before")
    load(covariate_path)
    covariate <- get(setdiff(ls(), before))
    
    # Confirm it's in CTFS format
    if(all(c("mat", "xdim", "ydim") %in% names(covariate))) {
      xrange <- c(0, covariate$xdim)
      yrange <- c(0, covariate$ydim)
      list(result = matrix_to_im(mat = covariate$mat, 
                                 xrange = xrange, 
                                 yrange = yrange,
                                 scale = scale),
           xrange = xrange,
           yrange = yrange)
    } else if(spatstat.geom::is.im(covariate)) {
      list(result = covariate,
           xrange = covariate$xrange,
           yrange = covariate$yrange)
    } else {
      stop(paste0("Detected as an RData, but could not recognize the format of the R object: ", covariate_path))
    }
  } else if(extension == "asc") {
    # TODO?:The BCI covariates have header size 6; might need to change it for other datasets
    header_length <- 6
    header <- read.table(covariate_path, sep = "", dec = ".", nrows = header_length)
    header <- setNames(lapply(seq_len(nrow(header)), function(row) header[row, 2]), nm = header[, 1])
    mat <- as.matrix(read.table(covariate_path, sep = "", dec = ".", skip = header_length))
    mat <- mat[nrow(mat):1, ]
    mat[mat == header$NODATA_VALUE] <- NA
    
    # Check that there's not obvious issue with the data
    if(ncol(mat) != header$NCOLS | nrow(mat) != header$NROWS) {
      stop(paste0("Number of rows/columns in matrix does not correspond to header in ", covariate_path))
    }
    
    # Extract xrange and yrange from header
    step <- header$CELLSIZE
    x0 <- header$XLLCORNER + step / 2
    y0 <- header$YLLCORNER + step / 2
    xmax <- x0 + (ncol(mat) - 1) * step
    ymax <- y0 + (nrow(mat) - 1) * step
    xrange <- c(x0, xmax)
    yrange <- c(y0, ymax)
    
    list(result = matrix_to_im(mat = mat, 
                               xrange = xrange, 
                               yrange = yrange,
                               scale = scale),
         xrange = xrange, 
         yrange = yrange)
  } else if(extension == "tif") {
    # Extract raster from path
    covariate <- raster::raster(covariate_path)
    
    # Scale raster?
    if(scale) {
      covariate = raster::scale(covariate)
    }
    
    # Get xrange and yrange from raster
    ext <- raster::extent(covariate)
    xrange <- c(raster::xmin(ext), raster::xmax(ext))
    yrange <- c(raster::ymin(ext), raster::ymax(ext))
    
    list(result = as.im.RasterLayer(covariate),
         xrange = xrange,
         yrange = yrange)
  } else {
    stop(paste0("Unknown extension for covariate ", covariate_path))
  }
}
