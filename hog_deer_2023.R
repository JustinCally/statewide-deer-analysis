library(tidyverse)
library(dbplyr)
library(sf)
library(terra)
library(stars)
library(DBI)
library(data.table)
library(Distance)
library(camtrapR)
library(cmdstanr)
library(activity)
library(bayesplot)
library(cowplot)
library(nngeo)
library(posterior)
library(tidyterra)
# library(MyFunctions)

# STAN settings
#set_cmdstan_path("C:/users/dr01/cmdstan-2.30.1")

# source("r/functions.r")

cams_curated<- read_csv("data/hog_cam_operation_2023.csv")

camtrap_records_deer<- read_csv("data/hog_camtrap_records_2023.csv")

site_vars<- read_csv("data/hog_site_vars.csv")

site_vars<- site_vars |> mutate(SiteID = factor(SiteID)) |> arrange(SiteID)

### Tidy format for the records
dcount<- camtrap_records_deer %>%
  dplyr::select(Species, SiteID, Distance, size, Date, Time, Behaviour) %>%
  full_join(cams_curated %>%
              dplyr::select(SiteID, Tkt, Effort), by="SiteID") %>%
  mutate(SiteID = factor(SiteID), object=1:nrow(.)) %>%
  mutate(size = if_else(is.na(size),0L, size)) %>%
  arrange(SiteID)



### Format data in a time format for availability

dcount<- dcount |> mutate(Distance = factor(Distance,
                                    levels = c("0 - 2.5", "2.5 - 5", "5 - 7.5", "7.5 - 10", "10+")))

summarised_time <- dcount %>%
  mutate(rtime = as.integer(format(Time, "%H")),
         time_cut = cut(rtime, breaks =  seq(from = 0, to = 24, by = 2)),
         time_midpoint = 2*as.integer(time_cut)-1)

# Basic detection model
det_formula <- ~ scale(HerbaceousUnderstoryCover)
det_model_matrix <- model.matrix(det_formula, data = site_vars)


## State Space Covariates --------------------------------------------------

pred_sf<- st_as_sf(st_as_stars(deer.window)) # entire state space

bioclim<- rast("g:/processed_rasters/bioclim_stack.tif")
EVC<- rast("g:/processed_rasters/EVC.tif")
TWIND<- rast("g:/processed_rasters/twind.tif")
trees<- rast("g:/processed_rasters/tree_density.tif")
#TWIND<- terra::scale(TWIND)
#trees<- terra::scale(trees)

evc <- exactextractr::exact_extract(EVC, pred_sf, coverage_area=TRUE)
evc<- bind_rows(evc, .id = "ID") |>
  group_by(ID, value) |>
  summarize(total_area = sum(coverage_area)) |>
  group_by(ID) |>
  mutate(proportion = total_area/sum(total_area)) |> dplyr::select(-total_area)

evc<- pivot_wider(evc, id_cols = ID, names_from = value, values_from = proportion, values_fill = 0)
evc<- evc |> dplyr::select(-`NA`) |> mutate(ID = as.integer(ID)) |> arrange(ID)
evc<- evc |> rename_with(.cols=`8`:`20`, ~ paste0("V", .x, recycle0 = TRUE))

twind<- terra::extract(TWIND, pred_sf, fun = mean, na.rm=TRUE)
twind<- twind |> mutate(TWIND = coalesce(TWIND, 0))
tree_dens<- terra::extract(trees, pred_sf, fun = mean, na.rm=TRUE)
tree_dens<- tree_dens |> mutate(TreeDensity = coalesce(TreeDensity, 0))

pred_sf<- pred_sf |> mutate(ID = row_number()) |> left_join(evc, by = "ID")
pred_sf<- pred_sf |> left_join(twind, by = "ID")
pred_sf<- pred_sf |> left_join(tree_dens, by = "ID")
pred_sf<- pred_sf |> mutate_at(c("TWIND","TreeDensity"), scale2)
pred_sf$X<- as_tibble(st_coordinates(st_centroid(pred_sf))) |> pull(X)/1000
pred_sf$Y<- as_tibble(st_coordinates(st_centroid(pred_sf))) |> pull(Y)/1000

# Find cell of for each camera location
sites_sf<- site_vars %>% st_as_sf(coords=c("Longitude","Latitude"), crs=4326) |> st_transform(3111)
cell_nos<- unlist(st_nn(sites_sf, pred_sf))
sites<- slice(pred_sf, cell_nos) |> bind_cols(site_vars) |> arrange(SiteID)

n_evc_cells <- sites |> st_set_geometry(NULL) |> group_by() |>
  summarise(across(V8:V20, ~ length(.x[.x>0])))
n_evc_cells<- sort(unlist(as.vector(n_evc_cells)))

# Assemble distance related stuff -------------------------
n_site <- length(unique(dcount$SiteID))
n_distance_bins <- 5
delta <- 2.5
midpts <- c(1.25, 3.75, 6.25, 8.75, 11.25)
max_distance <- 12.5
theta_frac <- 42/360
Tkt <- dcount %>% group_by(SiteID) %>% slice(1) %>% pull(Tkt)

# Get number of observations at each group size
n_obs <- table(dcount$SiteID, dcount$size)[,-1] # remove zeros
n_gs <- ncol(n_obs)
gs <- as.integer(colnames(n_obs))

yl <- list()

for(i in 1:n_gs) {
  yl[[i]] <- dcount %>%
    mutate(SiteID = factor(SiteID)) %>%
    group_by(SiteID, Distance, .drop = FALSE) %>%
    filter(size == gs[i]) %>%
    summarise(n = length(size)) %>%
    reshape2::dcast(SiteID ~ Distance, value.var = "n") %>%
    dplyr::select(-SiteID, -`NA`)

  yl[[i]][is.na(yl[[i]])] <- 0
}

y <- str2str::ld2a(yl)

### Adjusted availability function

pa <- matrix(data = NA, nrow = max(gs), ncol = length(midpts))

for(j in 1:max(gs)) {

  pa[j, ] <- sapply(midpts, function(x) {pr_closest(x-(0.5*delta),
                                                    x+(0.5*delta),
                                                    max_distance = max_distance,
                                                    n = j)})
}


#### Determine availability prior ####
trigger.events <- camtrap_records_deer

trigger.events$rtime <- gettime(trigger.events$Time, "%Y-%m-%d %H:%M:%OS", scale="radian")

act_result <- fitact(trigger.events$rtime, sample="model", reps=50, bw = 30)

avail <- list(creation=data.frame(rate = act_result@act[1],
                                  SE   = act_result@act[2]))

# View activity and avialbility
plot(act_result)
avail # looks to be a bit higher than other deer and more variable due to lower counts


beta.pars<- get.beta.param(avail$creation$rate, 2*avail$creation$SE)

# Intercept only abundance model: change to informative predictors

ab_formula <- ~ V2 + V19 + V9 + V11
#ab_formula <- ~ 1
ab_model_matrix <- model.matrix(ab_formula, data = sites)
m_psi =  ncol(ab_model_matrix)

### Spline-based model------------------------------------------------

library(mgcv)


ss<- sites |> st_set_geometry(NULL)

oo<- dcount |> mutate(grp = if_else(size > 0, 1, 0)) |> group_by(SiteID) |> summarise(Count=sum(grp))
oo<- oo |> left_join(ss, by = "SiteID")

form<- formula(Count ~ V2 + V19 + V9 + V11 + s(X,Y, bs = 'ds',m=c(1,.5)))
mod<- jagam(form, data=oo, family = poisson, file="jagdm.jags")

ifx<- 1:5
Z.mat<- mod$jags.data$X
Z.mat<- Z.mat[,-ifx] # remove explicit intercept
n_knots<- ncol(Z.mat)
S1<- mod$jags.data$S1
ns<- dim(S1)[2]

data = list(delta = 2.5,
            n_site = n_site,
            n_gs = n_gs,
            gs = gs,
            n_distance_bins = n_distance_bins,
            midpts = midpts,
            max_distance = max_distance,
            max_int_dist = as.integer(round(max_distance)),
            theta_frac = theta_frac,
            effort = Tkt,
            n_obs = n_obs,
            y = y,
            pa = pa, # availability varies across sites based on number of individuals
            det_model_matrix = det_model_matrix,
            det_ncb = ncol(det_model_matrix),
            bshape = beta.pars[[1]],
            bscale = beta.pars[[2]],
            m_psi =  m_psi,
            X_psi = ab_model_matrix,
            Z=Z.mat, S1 = S1, nk=n_knots)


ni <- 500
nw <- 500
nt <- 1
nc <- 3

inits = lapply(1:nc, function(i) list(beta_det=runif(2),
                                      #beta_psi = c(-5, rnorm(m_psi-1)),
                                      activ = runif(1),
                                      alpha = 1, rho=2,
                                      eps_ngs = rep(1/n_gs, n_gs),
                                      eta = rnorm(n_site),
                                      b=mod$jags.ini$b[-ifx],
                                      sdb=1))


### STAN model

model <- cmdstan_model("stan/count_only_spline.stan")

# Runs in 1-2 minutes
out <- model$sample(data = data,
                    chains = nc,
                    init=inits,
                    parallel_chains = nc,
                    show_messages = FALSE,
                    save_warmup = FALSE,
                    iter_sampling = ni,
                    iter_warmup = nw, max_treedepth = 14)

out$summary(c("beta_psi","beta_det"))

out$summary(c("sdb","b"))

dens <- out$summary("N_site")
dens2 <- out$summary("Site_lambda")

### Random effects model -------------------------------------------------------

  data = list(delta = 2.5,
              n_site = n_site,
              n_gs = n_gs,
              gs = gs,
              n_distance_bins = n_distance_bins,
              midpts = midpts,
              max_distance = max_distance,
              max_int_dist = as.integer(round(max_distance)),
              theta_frac = theta_frac,
              effort = Tkt,
              n_obs = n_obs,
              y = y,
              pa = pa, # availability varies across sites based on number of individuals
              det_model_matrix = det_model_matrix,
              det_ncb = ncol(det_model_matrix),
              bshape = beta.pars[[1]],
              bscale = beta.pars[[2]],
              m_psi =  ncol(ab_model_matrix),
              X_psi = ab_model_matrix)


ni <- 500
nw <- 500
nt <- 1
nc <- 3

 inits = lapply(1:nc, function(i) list(beta_det=runif(2),
                                       activ = runif(1),
                                       site_sd = runif(1),
                                       eps_ngs = rep(1/n_gs, n_gs),
                                       site_raw = rnorm(n_site)))


### STAN model

model <- cmdstan_model("stan/count_only.stan")

# Runs in 1-2 minutes
out <- model$sample(data = data,
                           chains = nc,
                           init=inits,
                           parallel_chains = nc,
                           show_messages = FALSE,
                           save_warmup = FALSE,
                           iter_sampling = ni,
                           iter_warmup = nw)

out$summary(c("beta_psi","beta_det","site_sd"))

out$summary(c("sigma"))

dens <- out$summary("N_site")
dens2 <- out$summary("Site_lambda")

# Bind the density to the camera information
density_at_sites <- cbind(cams_curated, dens) %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4236)

## Predictions
pred_sf<- st_as_sf(st_as_stars(deer.window))
pred_grid <- st_coordinates(st_centroid(pred_sf))/1000

beta_posterior<- out$draws('beta_psi', format = "df")$`beta_psi[1]`

gp_preds <- matrix(NA, nrow = 5, ncol = NROW(pred_grid))
for(i in 1:NROW(gp_preds)){
  cat('Processing posterior draw', i, 'of 5...\n')
  gp_preds[i,] <- predict_gp(obs_distances = obs_distances,
                             pred_distances = pred_distances,
                             obs_to_pred_distances = obs_to_pred_distances,
                             obs_gp_function = as.vector(as.matrix(obs_gp_function_posterior[i,])),
                             alpha = alpha_posterior[i],
                             rho = 10,#rho_posterior[i],
                             beta = beta_posterior[i],
                             kernel = "matern32")
}

preds<- apply(gp_preds,2,median)

pred_sf<- st_as_sf(st_as_stars(deer.window))
pred_sf<- st_centroid(pred_sf)
cn<- data.frame(cells(deer.window, vect(pred_sf)))

pred_rast<- deer.window
pred_rast[cn$cell]<- preds

win.graph(15,10)
ggplot() +
  geom_spatraster(aes(fill=layer), data=pred_rast) +
  scale_fill_continuous(type = "viridis") +
  geom_sf(color="black", data = sites_sf)

### Gaussian Process version------------------------------------------------

data = list(delta = 2.5,
            n_site = n_site,
            n_gs = n_gs,
            gs = gs,
            n_distance_bins = n_distance_bins,
            midpts = midpts,
            max_distance = max_distance,
            max_int_dist = as.integer(round(max_distance)),
            theta_frac = theta_frac,
            effort = Tkt,
            n_obs = n_obs,
            y = y,
            pa = pa, # availability varies across sites based on number of individuals
            det_model_matrix = det_model_matrix,
            det_ncb = ncol(det_model_matrix),
            bshape = beta.pars[[1]],
            bscale = beta.pars[[2]],
            m_psi =  ncol(ab_model_matrix),
            X_psi = ab_model_matrix,
            coords = coords)


ni <- 500
nw <- 500
nt <- 1
nc <- 3

inits = lapply(1:nc, function(i) list(beta_det=runif(2),
                                      activ = runif(1),
                                      alpha = 1, rho=2,
                                      eps_ngs = rep(1/n_gs, n_gs),
                                      eta = rnorm(n_site)))


### STAN model

model <- cmdstan_model("stan/count_only_gp.stan")

# Runs in 1-2 minutes
out <- model$sample(data = data,
                    chains = nc,
                    init=inits,
                    parallel_chains = nc,
                    show_messages = FALSE,
                    save_warmup = FALSE,
                    iter_sampling = ni,
                    iter_warmup = nw)

out$summary(c("beta_psi","beta_det","alpha","rho"))

dens <- out$summary("N_site")
dens2 <- out$summary("Site_lambda")

##----- predictions GP process ---------------------------
pred_sf<- st_as_sf(st_as_stars(deer.window))
pred_grid <- st_coordinates(st_centroid(pred_sf))/1000

obs_distances <- calc_point_distances(coordinates = coords)
pred_distances <- calc_point_distances(coordinates = pred_grid)
obs_to_pred_distances <- calc_point_distances(coords, pred_grid)

obs_gp_function_posterior <- out$draws('gp_predict',format = "df")[,1:n_site]
rho_posterior <- out$draws('rho', format = "df")$rho
alpha_posterior <- out$draws('alpha', format = "df")$alpha
beta_posterior<- out$draws('beta_psi', format = "df")$`beta_psi[1]`

gp_preds <- matrix(NA, nrow = 5, ncol = NROW(pred_grid))
for(i in 1:NROW(gp_preds)){
  cat('Processing posterior draw', i, 'of 5...\n')
  gp_preds[i,] <- predict_gp(obs_distances = obs_distances,
                             pred_distances = pred_distances,
                             obs_to_pred_distances = obs_to_pred_distances,
                             obs_gp_function = as.vector(as.matrix(obs_gp_function_posterior[i,])),
                             alpha = alpha_posterior[i],
                             rho_posterior[i],
                             beta = beta_posterior[i],
                             kernel = "quad")
}

preds<- apply(gp_preds,2,median)

pred_sf<- st_as_sf(st_as_stars(deer.window))
pred_sf<- st_centroid(pred_sf)
cn<- data.frame(cells(deer.window, vect(pred_sf)))

pred_rast<- deer.window
pred_rast[cn$cell]<- preds

win.graph(15,10)
ggplot() +
  geom_spatraster(aes(fill=layer), data=pred_rast) +
  scale_fill_continuous(type = "viridis") +
  geom_sf(color="black", data = sites_sf)

##----- predictions Splines ---------------------------

#pred_sf<- st_as_sf(st_as_stars(deer.window))
pred_grid <- pred_sf
OS<- as.vector(st_area(pred_grid))/1e6

Xpred<- model.matrix(ab_formula, data=pred_grid)

## Jaggam
mod<- gam(form, data=oo, family = poisson)
Z.mat<- predict(mod, pred_grid, type="lpmatrix")
Z.mat<- Z.mat[,-ifx]

preds<- predict_model(samples=out, X=Xpred, Z=Z.mat, center=FALSE, Offset=OS, all.sites=TRUE)

cn<- data.frame(cells(deer.window, vect(st_centroid(pred_grid))))

pred_rast<- deer.window
pred_rast[cn$cell]<- preds$N

win.graph(15,10)
ggplot() +
  geom_spatraster(aes(fill=layer), data=pred_rast) +
  scale_fill_continuous(type = "viridis") +
  geom_sf(color="black", data = sites_sf)

