library(terra)
library(VicmapR)
library(sf)
library(dplyr)
# get woody vegetation
l <- listLayers()
sveg <- vicmap_query("open-data-platform:sveg100") %>%
  # head(100) %>%
  collect()

vic_bound <- vicmap_query("open-data-platform:delwp_region") %>%
  # filter(state == "VIC") %>%
  collect() %>%
  st_transform(3111) %>%
  vect()

template_raster <- rast(st_bbox(vic_bound) %>% terra::ext(), res = 1000)

# landuse <- vicmap_query("open-data-platform:landuse_2014")

sveg_forest <- sveg %>%
  filter(x_vegform %in% c("OPEN FOREST") | stringr::str_detect(x_vegform, "WOODLAND"))


# woody_raster <- sveg_forest %>%
#   mutate(veg = 1) %>%
#   st_transform(3111) %>%
#   vect() %>%
#   rasterize(y = template_stack, field = "veg", background = 0)

sveg_vect <- vect(sveg_forest %>%
                    st_transform(3111) %>%
                    st_union()) |> as.lines()

woody_forest_edges <- rasterizeGeom(sveg_vect, template_raster, "length")

writeRaster(woody_forest_edges, "data/woody_forest_edges.tif", overwrite = TRUE)

# r<-woody_raster
# #3x3 sobel filters to detect edges in x and y directions
# fx=matrix(c(-1,-2,-1,0,0,0,1,2,1), nrow=3)
# fy=matrix(c(1,0,-1,2,0,-2,1,0,-1), nrow=3)
#
# wood_edge_x<-focal(r, w=fx, fun="sum", silent=FALSE)
# wood_edge_y<-focal(r, w=fx, fun="sum", silent=FALSE)
#
# plot(wood_edge_x) #horizontal edges
# plot(wood_edge_y) #vertical edges
# #combining the x and y components of edginess to get the total magnitude of edge effect
# wood_edge<-sqrt(wood_edge_x*wood_edge_x + wood_edge_y*wood_edge_y)
#
# plot(wood_edge)
