library(tidyverse)
library(tidybayes)
library(here)
library(mgcv)
library(magrittr)

# copied from prep_summer_flounder and edited down; revisit this if study domain or other design decisions change 
dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))
dat_test <- read_csv(here("processed-data","flounder_catch_at_length_fall_testing.csv"))

dat_test_summ <- dat_test %>% 
  group_by(haulid, year, lat, btemp) %>% 
  summarise(dens = sum(number_at_length)) %>% 
  group_by(year) %>% 
  summarise(quant_0.05 = Hmisc::wtd.quantile(lat, weights=dens, probs=0.05),
            quant_0.5 = Hmisc::wtd.quantile(lat, weights=dens, probs=0.5),
            quant_0.95 = Hmisc::wtd.quantile(lat, weights=dens, probs=0.95))

use_patches <- dat %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor)

patches <- sort(unique(use_patches$lat_floor))
np = length(patches) 

gamdat <- dat %>% 
  group_by(haulid, year, lat, btemp) %>% 
  summarise(dens = sum(number_at_length)) %>% 
  ungroup() 

gamdat_test <- dat_test %>% 
  group_by(haulid, year, lat, btemp) %>% 
  summarise(dens = sum(number_at_length)) %>% 
  ungroup() 

# check data coverage of testing data
gamdat_test %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(year, lat_floor) %>% 
  mutate(n_tot = n(), 
         test = is.na(btemp)) %>% 
  filter(test==FALSE) %>% 
  mutate(n_pres = n(), 
         prop_btemp = n_pres / n_tot) %>% 
  select(prop_btemp, year, lat_floor) %>% 
  distinct() %>% 
  ggplot() + 
  geom_tile(aes(x=year, y=lat_floor, fill=prop_btemp, color=prop_btemp)) 
# because so much data is missing in 2008, drop it from predictions 

gamdat_test %<>% 
  filter(!year==2008)

gam1 <- gam(dens ~ s(lat) + s(btemp), data = gamdat, family = "nb") # revisit family choice. ziP, ziplss, poisson?
plot(gam1)
gam.check(gam1) # not enough basis functions? revisit 

gamdat_pred <- gamdat_test %>% 
  mutate(dens_pred = predict.gam(gam1, newdata = gamdat_test %>% select(-dens)))
# issue here--it's predicting negative densities--which is weird for a zero-inflated Poisson 

gam_out <- gamdat_pred %>% 
  filter(!is.na(dens_pred)) %>% 
  rowwise() %>% 
  mutate(dens_pred = max(dens_pred, 0),
         lat_floor = floor(lat)) %>% 
  group_by(year) %>% 
  summarise(quant_0.05 = Hmisc::wtd.quantile(lat, weights=dens_pred, probs=0.05),
            quant_0.5 = Hmisc::wtd.quantile(lat, weights=dens_pred, probs=0.5),
            quant_0.95 = Hmisc::wtd.quantile(lat, weights=dens_pred, probs=0.95))
  
results.tmp <- file.path('~/github/mid_atlantic_forecasts/results/v0.50')

observed_abund_posterior_predictive <- read_csv(file.path(results.tmp, "density_obs_proj.csv"))
range_quantiles_proj <- read_csv(file.path(results.tmp, "range_quantiles_proj.csv"))

centroid_proj <- read_csv(file.path(results.tmp, "centroid_proj.csv")) %>% 
  mutate(year = as.numeric(year))
centroid <- dat_test_dens %>% 
  group_by(year) %>% 
  summarise(lat_obs=weighted.mean(lat_floor, w=mean_dens))

load("~/github/mid_atlantic_forecasts/processed-data/stan_data_prep.Rdata")

# centroid 
centroid_summary <- centroid_proj %>% 
  left_join(centroid) %>% 
  mutate(resid = lat - lat_obs,
         resid_sq = resid^2)

## mean bias
mean(centroid_summary$resid)

## mean RMSE
sqrt(mean(centroid_summary$resid_sq))
  
## bias over time
centroid_summary %>% 
  group_by(year) %>% 
  summarise(bias = mean(resid))

## RMSE over time 
centroid_summary %>% 
  group_by(year) %>% 
  summarise(rmse = sqrt(mean(resid_sq)))

# index of abundance 

