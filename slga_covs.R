#### Compiling Predictor Variables ####

#### Predictor Set A ####

#### Raster-based Data ####
library(terra)
library(SLGACloud) # download from github
library(dplyr)
library(sf)
library(VicmapR)

#### Load in victoria ####

vic_bound <- vicmap_query("open-data-platform:delwp_region") %>%
  # filter(state == "VIC") %>%
  collect()

#######################  PLEASE NOTE  ############################################
####
####  Note the use of the "/vsicurl/" preface to all URLs
####  This tells the terra package to treat these rasters as COGs.
####  ie. use the gdal virtual file system when raeding the data.
####  If you do not prepend "/vsicurl/" to your ULs terra will treat them as
####  normal GeoTIFF files and try to download the whole raster for every operation.

vars <- getProductMetaData(Detail = 'High', Product = "90m_Covariate", DataType = "DSM_Covariate") %>%
  select(Name, Description, everything())

#View(vars)

# vars_to_download <- vars$Name

# vars_to_download <- c("Clim_ADM",
#                       "Clim_PTA",
#                       "Clim_PTS1",
#                       "Clim_TNI",
#                       "Clim_TRA",
#                       "Relief_aspect_1s",
#                       "Relief_twi_1s",
#                       "Relief_slopepct1s",
#                       "Soil_SMOS_Summary_mean",
#                       "Veg_alpsbk_aust_y2009_sd5a2",
#                       "Veg_FC_Mean_PV",
#                       "Veg_NDVI_mean_Q2")

vars_to_download <- c("Relief_elev_focalrange1000m_3s",
                      "Relief_mrvbf_3s_mosaic",
                      "Relief_roughness",
                      "Relief_slope_perc",
                      "Relief_twi_3s",
                      "Veg_AVHRR_FPAR_Mean",
                      "Veg_Modis_FractCover_Mean_BS",
                      "Veg_Modis_FractCover_Mean_NPV",
                      "Veg_Modis_FractCover_Mean_PV",
                      "Veg_LandCoverTrend_evi_mean",
                      "Clim_ADM",
                      "Clim_ADRX",
                      "Clim_EAAS",
                      "Clim_fwofs",
                      "Clim_PTA",
                      "Clim_PTI",
                      "Clim_PTRI",
                      "Clim_PTX",
                      "Clim_RSM",
                      "Clim_RSRX",
                      "Clim_TNI",
                      "Clim_TNM",
                      "Clim_TNX",
                      "Clim_TRA",
                      "Clim_WDA",
                      "Clim_WDI",
                      "Clim_WDRI",
                      "Relief_mrrtf_3s",
                      "Veg_NDVI_mean_Q1",
                      "Veg_NDVI_mean_Q2",
                      "Veg_NDVI_mean_Q3",
                      "Veg_NDVI_mean_Q4")

##### Access the COG directly from the DataStore #######


rpaths <- paste0('/vsicurl/', vars %>%
                   filter(Name %in% vars_to_download) %>%
                   pull(COGsPath))

r <- rast(rpaths)



##### Clip out a section of the COG and save locally #######

e <- st_bbox(vic_bound) %>% terra::ext()
rc <- crop(r, e)
reproj_r <- project(rc, "epsg:3111", res = 1000, method = "bilinear")
# plot(rc)

writeRaster(reproj_r, "/Volumes/Cally_Camtr/StatewideRasters/slga_stack.tif", overwrite = T)
