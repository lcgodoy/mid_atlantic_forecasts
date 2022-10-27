#####
## Fitting VAST model to summer flounder data for comparison with DRM 
#####


# Overview ----------------------------------------------------------------
# The objective of this script is to fit an annaul time step VAST model to summer flounder data from the Northeast Shelf Large Marine Ecosystem, focusing on validating predictive performance within latitudinal bands and comparing with an alternative size-based dynamic range model. To accomplish this objective, the script walks through X different stages -- from initial data prep, through to making predictions from a fitted VAST model. 

# IF YOU ALREADY HAVE THE FITTED MODEL -- run the preliminaries section and then skip to line 193: "Read in fitted model"


# Preliminaries -----------------------------------------------------------
library(devtools)
library(TMB)
library(VAST)
library(tidyverse)
library(lubridate)
library(sf)
library(gstat)
library(raster)
library(here) 
library(splines)
library(patchwork)
library(PresenceAbsence)
library(MLmetrics)

sf_use_s2(FALSE)
run_vast <- FALSE

# Functions -- sourced from aallyn/TargetsSDM. This includes a lot of helpful "VAST" functions.
devtools::source_url("https://raw.github.com/aallyn/TargetsSDM/main/R/vast_functions.r")
devtools::source_url("https://raw.github.com/aallyn/TargetsSDM/main/R/SDM_PredValidation_Functions.R")
devtools::source_url("https://raw.github.com/James-Thorson-NOAA/VAST/dev/R/reload_model.R")

# Initial data prep: trawl and sample data --------------------------------
flounder_train <- read_csv(here("processed-data","flounder_biomass_fall_training.csv")) %>% filter(!is.na(btemp))
flounder_test <- read_csv(here("processed-data","flounder_biomass_fall_testing.csv"))%>% filter(!is.na(btemp))


if(run_vast=TRUE){
  # VAST data --------------------------------------
  # First, a sample dataframe. With the VAST framework, we can have all the data (training and testing) and then switch Pred_TF to "1" for the "testing" data years. This means they will not be used in the likelihood, but we will predict to them.
  flounder_train$Pred_TF<- 0
  flounder_test$Pred_TF<- 1
  flounder_all<- bind_rows(flounder_train, flounder_test)
  #flounder_all<- bind_rows(flounder_train)
  vast_samp_dat<- data.frame("Year" = flounder_all$year, "Lat" = flounder_all$lat, "Lon" = flounder_all$lon, "Biomass" = flounder_all$wtcpue, "Swept" = rep(1, nrow(flounder_all)), "Pred_TF" = flounder_all$Pred_TF)
  
  # Next, a covariate dataframe and then rescaling bottom temperature
  vast_cov_dat<- data.frame("Year" = flounder_all$year, "btemp" = flounder_all$btemp, "Lat" = flounder_all$lat, "Lon" = flounder_all$lon) 
  vast_cov_dat$btemp<- as.numeric(scale(vast_cov_dat$btemp), center = TRUE, scale = TRUE)
  
  # VAST model objectives and settings --------------------------------------
  # Extrapolation grid and Strata limits. Here, we want to essentially arrive at latitude patch estimates. 
  # Grabbing the grid for northwest Atlantic in FishStatsUtils
  utils::data(northwest_atlantic_grid, package = "FishStatsUtils")
  Data_Extrap<- northwest_atlantic_grid
  str(Data_Extrap)
  
  # Looking at this quickly...
  plot(Data_Extrap$Lon, Data_Extrap$Lat)
  
  # Compare with the data...
  points(flounder_all$lon, flounder_all$lat, col = "red")
  
  # Okay, clearly want to do some work here. First is the crazy outlier in Chesapeake Bay and then the extrapolation grid also extends into area now surveyed more commonly by NOAA SEFSC. For the Chesapeake one...
  #identify(flounder_all$lon, flounder_all$lat)
  plot(Data_Extrap$Lon, Data_Extrap$Lat)
  points(flounder_all$lon, flounder_all$lat, col = "red")
  points(flounder_all$lon[922], flounder_all$lat[922], col = "green")
  
  flounder_all<- flounder_all[-922,]
  points(flounder_all$lon, flounder_all$lat, col = "blue")
  
  # Remove this obs from sample and covariate data
  vast_samp_dat<- vast_samp_dat[-922,]
  vast_cov_dat<- vast_cov_dat[-922,]
  
  # Next, work on the extrapolation grid. Go based on min latitude of observations and a little buffer (0.25 degrees)
  Data_Extrap<- Data_Extrap %>%
    filter(., Lat >= floor(min(flounder_all$lat)))
  points(Data_Extrap$Lon, Data_Extrap$Lat, col = "green")
  
  # Now, we need to have a new column that has the "STRATA" of each observation. In this case, we want those strata to be the latitude patches. From `analyze_summer_flounder` grabbing the different "lat" categories...
  vast_extrap_grid<- Data_Extrap %>%
    mutate(STRATA = factor(paste0("Lat_", floor(Lat)), levels = paste0("Lat_", seq(from = min(unique(floor(Data_Extrap$Lat))), to = max(unique(floor(Data_Extrap$Lat))))))) %>%
    rename("Area_km2" = Area_in_survey_km2) %>%
    dplyr::select(., Lon, Lat, Area_km2, STRATA)
  
  # Finally, strata.limits dataframe
  #strata_use<- strata_use<- data.frame("STRATA" = as.character(unique(vast_extrap_grid$STRATA)))
  strata_use<- data.frame( "STRATA" = paste0("Lat_", seq(from = min(unique(floor(Data_Extrap$Lat))), to = max(unique(floor(Data_Extrap$Lat))))),
                           "south_border" = seq(from = min(unique(floor(Data_Extrap$Lat))), to = max(unique(floor(Data_Extrap$Lat)))),
                           "north_border" = seq(from = min(unique(floor(Data_Extrap$Lat)))+0.99, to = max(unique(floor(Data_Extrap$Lat)))+0.99))
  
  vast_extrap_info<- make_extrapolation_info(Region = "User", strata.limits = strata_use, input_grid = vast_extrap_grid, grid_dim_km = c(25, 25))
  colSums(vast_extrap_info$a_el)
  
  # Quick check...
  ggplot() +
    geom_point(data = vast_extrap_grid, aes(x = Lon, y = Lat, color = STRATA))
  
  # Field and rho configuration settings. Field config sets up the spatial/spatio-temporal components and how many factors should be estimated. Here, we are going to turn both of those "on". Rho config sets up autoregressive structure on intercepts and spatio-temporal components. Given the interest here on forecasting, going to try AR1 process for both the intercepts and spatio-temporal variability. Will make adjustments as needed given data constraints and convergence issues. 
  field_config<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
  rho_config<- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 2, "Epsilon2" = 2)
  
  vast_settings<- make_settings(n_x = 400, Region = "User", purpose = "index2", FieldConfig = field_config, RhoConfig = rho_config, ObsModel = c(2, 1), OverdispersionConfig = c(0, 0), bias.correct = TRUE, knot_method = "samples", Options = c("Calculate_Range" = TRUE), strata.limits = strata_use)
  
  # Settings Version error ends up occuring...
  vast_settings$Version<- "VAST_v14_0_1"
  
  # Fitting VAST model ------------------------------------------------------
  nice_category_names<- "Summer flounder all"
  # Model formula
  hab_formula<- ~ bs(btemp, degree = 3, intercept = FALSE) 
  
  # Build VAST model
  vast_build<- fit_model("settings" = vast_settings, "Method" = vast_settings$Method, "input_grid" = vast_extrap_grid, "Lat_i" = vast_samp_dat[, 'Lat'], "Lon_i" = vast_samp_dat[, 'Lon'], "t_i" = vast_samp_dat[, 'Year'], "c_i" = rep(0, nrow(vast_samp_dat)), "b_i" = vast_samp_dat[, 'Biomass'], "a_i" = vast_samp_dat[, 'Swept'], "PredTF_i" = vast_samp_dat[, 'Pred_TF'], "covariate_data" = vast_cov_dat, "X1_formula" = hab_formula, "X2_formula" = hab_formula, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = TRUE)
  
  # Looks good. Fit it.
  vast_fit <- fit_model("settings" = vast_settings, "Method" = vast_settings$Method, "input_grid" = vast_extrap_grid, "Data_Extrap" = vast_build$extrapolation_list$Data_Extrap, "Lat_i" = vast_samp_dat[, "Lat"], "Lon_i" = vast_samp_dat[, "Lon"], "t_i" = vast_samp_dat[, "Year"], "c_i" = rep(0, nrow(vast_samp_dat)), "b_i" = vast_samp_dat[, "Biomass"], "a_i" = vast_samp_dat[, "Swept"], "PredTF_i" = vast_samp_dat[, "Pred_TF"], "covariate_data" = vast_cov_dat, "X1_formula" = hab_formula, "X2_formula" = hab_formula, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = TRUE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = TRUE)
  saveRDS(vast_fit, file = here::here("results/vast/summer_flounder_vast.rds"))
}

# Read in fitted model --------------------------------------
vast_fit <- readRDS(here::here("results/vast/summer_flounder_vast.rds"))
#vast_fit <- reload_model(vast_fit)
nice_category_names <- "Summer flounder all"
# Summarizing VAST model predictions --------------------------------------
## Covariate effects
vast_habcovs_effs<- get_vast_covariate_effects(vast_fit = vast_fit, params_plot = c("btemp"), params_plot_levels = 100, effects_pad_values = c(), nice_category_names = nice_category_names, out_dir = here("results/vast"))

# Warnings are...interesting. I think this is because of how the sampling is done to generate the plots -- basically some sampled values are beyond the boundary knots (e.g., bottom temperature values). 
vast_habcovs_plot<- plot_vast_covariate_effects(vast_covariate_effects = vast_habcovs_effs, vast_fit = vast_fit, nice_category_names = nice_category_names, out_dir = here("results/vast"))
vast_habcovs_plot # Seems to suggest preference for warmer waters to determine general distribution, and then density (biomass) within that area shows more of a "humped" distribution and preference for intermediate btemps

## Prediction skill...this is likely more interesting relative to another model, but just showing as an example. Closest to the gold star is better.
# Get point level predictions
vast_preds <- vast_get_point_preds(vast_fit = vast_fit, use_PredTF_only = TRUE, nice_category_names = nice_category_names, out_dir = here("results/vast"))
vast_preds$Model_Name<- "Forecast Model"

td_plot_pres <- taylor_diagram_func(dat = vast_preds, obs = "Presence", mod = "Predicted_ProbPresence", group = NULL, color.cols = "black", fill.cols = "#1b9e77", shapes = 21, alpha = 0.6, out.file = paste0(here("results/vast"), "/PresenceTD.jpg"))
td_plot_pres

td_plot_bio <- taylor_diagram_func(dat = vast_preds, obs = "Biomass", mod = "Predicted_Biomass", group = NULL, color.cols = "black", fill.cols = "#1b9e77", shapes = 21, alpha = 0.6, out.file = paste0(here("results/vast"), "/BiomassTD.jpg"))
td_plot_bio

## Stratified biomass indices -- I think this is where the comparisons would be for the patches?
vast_bio<- get_vast_index_timeseries(vast_fit = vast_fit, all_times = unique(vast_fit$covariate_data$Year), nice_category_names = nice_category_names, index_scale = c("raw"), out_dir = here("results/vast"))
summary(vast_bio)

# Plot...
color_pal<- rev(viridis::viridis(n = length(unique(vast_bio$Index_Region))))
plot_vast_index_timeseries(index_res_df = vast_bio, index_scale = "raw", nice_category_names = nice_category_names, nice_xlab = "Year", nice_ylab= "Biomass index (metric tons)", paneling = "none", color_pal = color_pal, out_dir = here("results/vast"))

# Tile plot similar to AF's for the DRM
patch_bio_tile_plot<- ggplot(data = vast_bio, aes(x = Year, y = Index_Region, fill = Index_Estimate)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Biomass index (metric tons)") +
  theme_bw() +
  scale_x_continuous(name = "Year", breaks = seq(from = min(vast_bio$Year), to = max(vast_bio$Year), by = 4), expand = c(0, 0)) +
  scale_y_discrete(name = "Latitiude patch", expand = c(0, 0))
ggsave(paste0(here("results/vast/"), nice_category_names, " patch biomass tile.jpg"), patch_bio_tile_plot)
patch_bio_tile_plot

# A land shapefile and setting lon/lat limits
land_use <- st_read(here("raw-data/ne_50m_land.shp"))
xlim_use<- c(-78, -65)
ylim_use<- c(34, 46)

range_edges <- get_range_edges(vast_fit = vast_fit, all_times = unique(vast_fit$covariate_data$Year), category_names = nice_category_names, strata_names = NULL, n_samples = 100, quantiles = c(0.01, 0.5, 0.9), calculate_relative_to_average = FALSE)

range_edge_plot <- ggplot() +
  geom_errorbar(data = range_edges, aes(x = Year, ymin = (Latitude_Mean - Latitude_SD), ymax = (Latitude_Mean + Latitude_SD), color = Quantile, group = Quantile), alpha = 0.65) +
  geom_point(data = range_edges, aes(x = Year, y = Latitude_Mean, color = Quantile)) +
  geom_smooth(data = range_edges, aes(x = Year, y = Latitude_Mean, color = Quantile, group = Quantile), method = "gam", formula = y ~ s(x)) +
  scale_color_manual(values = c("Red", "Gray", "Blue")) +
  # scale_x_date(breaks = date_breaks, date_labels = "%Y") +
  xlab("Year") +
  ylab("Mean Latitude") +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
range_edge_plot
ggsave(paste0(here("results/vast/"), nice_category_names, " range_edges.jpg"), range_edge_plot)


# Extracting and plotting predicted density at knot locations, smoothed over a regular grid. This is a bit unruly given how many time steps there are.
# Need a regional grid...
region_mask <- st_read(here("raw-data/BTS_Strata.shp")) %>%
  st_transform(., crs = 4326) %>%
  st_union()

region_utm <- region_mask %>%
  st_transform(., crs = vast_fit$extrapolation_list$projargs)

cell_size = 25
region_grid <- st_make_grid(region_utm, cellsize = cell_size, what = "centers")
region_raster <- raster(crs = st_crs(region_grid), vals = 0, resolution = cell_size, ext = extent(st_bbox(region_grid))) %>%
  rasterize(as(region_mask, "Spatial"), .)
crs(region_raster) <- st_crs(region_grid)$input
region_raster <- projectRaster(region_raster, crs = 4326)
region_raster[] <- NA

pred_dens <- vast_fit_plot_spatial(vast_fit = vast_fit, manual_pred_df = NULL, pred_grid = region_raster, spatial_var = "D_gct", nice_category_names = nice_category_names, pred_label = "Pred", mask = region_mask, all_times = unique(vast_fit$covariate_data$Year), plot_times = NULL, land_sf = land_use, xlim = c(-78.5, -56), ylim = c(35, 48), lab_lat = 36, lab_lon = -60, panel_or_gif = "panel", panel_cols = 9, panel_rows = 5, out_dir = here::here("results/vast"), land_color = "#d9d9d9")
pred_dens
