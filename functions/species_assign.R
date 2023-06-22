
# Kernel density function
st_kde <- function(points, cellsize, bandwith, extent = NULL){

  if(is.null(extent)){
    extent_vec <- sf::st_bbox(points)[c(1,3,2,4)]
  } else{
    extent_vec <- sf::st_bbox(extent)[c(1,3,2,4)]
  }

  n_y <- ceiling((extent_vec[4]-extent_vec[3])/cellsize)
  n_x <- ceiling((extent_vec[2]-extent_vec[1])/cellsize)

  extent_vec[2] <- extent_vec[1]+(n_x*cellsize)-cellsize
  extent_vec[4] <- extent_vec[3]+(n_y*cellsize)-cellsize

  coords <- sf::st_coordinates(points)
  matrix <- MASS::kde2d(coords[,1],coords[,2],h = bandwith,n = c(n_x,n_y),lims = extent_vec)
  raster::raster(matrix)
}

nrd_new <- function (x)
{
  r <- quantile(x, c(0.25, 0.75))
  h <- (r[2L] - r[1L])/1.34
  4 * 1.06 * min(sqrt(var(x)), h, na.rm = T) * length(x)^(-1/5)
}

# Assign detections to species based on:
# species seen on camera = 1
# species not seen on camera and another species seen = 0
# species not seen on camera and no others seen, but within KDE distribution = 1
# species not seen on camera, and not within KDE distribution = 0

hif_hany <- function(..., values) {
  mtx <- do.call(cbind, list(...))
  mtx <- array(mtx %in% values, dim = dim(mtx))
  rowSums(mtx) > 0
}

detection_assign <- function(sp, species, othersp) {

  sp.d <- paste(species, "distribution")

  sf_m <- dplyr::mutate(sp, across(c(Pellet, Footprint, Rubbing, Wallow),
                                   ~ dplyr::case_when(.x == 0 ~ 0,
                                                      .x > 0 & !! sym(species) == 1 ~ 1,
                                                      .x > 0 & !! sym(species) != 1 & hif_hany(!!! syms(othersp), values = 1) ~ 0,
                                                      .x > 0 & !! sym(species) != 1 & !hif_hany(!!! syms(othersp), values = 1) & !! sym(sp.d) == 1 ~ 1,
                                                      .x > 0 & !! sym(species) != 1 & !! sym(sp.d) != 1 ~ 0)))

  transects_detection <- sf_m %>%
    dplyr::select(SiteID, Transect, Distance, Pellet, Footprint, Rubbing, Wallow) %>%
    pivot_longer(cols = c("Pellet", "Footprint", "Rubbing", "Wallow"),
                 values_to = "Presence", names_to = "Survey") %>%
    mutate(Presence = coalesce(Presence, 0))

  return(transects_detection)
}
