library(mgcv)
library(here)
library(tidyverse)
set.seed(42)

load(here("processed-data","stan_data_prep.Rdata"))
dat <- read_csv(here("processed-data","flounder_catch_for_sdm_fall_training.csv")) %>% 
  mutate(pres = ifelse(abundance==0, 0, 1),
         log_abundance = log(abundance))

dat_proj <- read_csv(here("processed-data","flounder_catch_for_sdm_fall_testing.csv"))  %>% 
  mutate(pres = ifelse(abundance==0, 0, 1),
         log_abundance = log(abundance))

# what fraction of rows have NA temperature values? 
nrow(dat %>% filter(is.na(btemp)))/nrow(dat) # 12.8% 
nrow(dat_proj %>% filter(is.na(btemp)))/nrow(dat_proj) # 12.0%

# is there a pattern to what years they are in? 
quantile((dat %>% filter(is.na(btemp)))$year) # slightly more in the earlier part of the time-series 

# GAM doesn't like NA predictors, need to remove
dat <- dat %>% filter(!is.na(btemp))
dat_proj <- dat_proj %>% filter(!is.na(btemp))

# delta model
# credit to https://github.com/stephbrodie1/Projecting_SDMs/blob/main/Estimation_Model_Functions/Fitting_GAMs.R for this code 


delta1_formula <- formula("pres ~ s(btemp)")
delta2_formula <- formula("log_abundance ~ s(btemp)")

# environmental only models ("E")
gam_E_P <- gam(delta1_formula, data=dat, family=binomial)
#plot(gam_E_P, pages=1)
gam_E_N <- gam(delta2_formula, data=dat[dat$abundance>0,], family=gaussian)
#plot(gam_E_N, pages=1)

presx <- predict(gam_E_P, dat, type="response")
abundx <- exp(predict(gam_E_N, dat, type="response"))
dat$gam_E <- presx * abundx

presx <- predict(gam_E_P, dat_proj, type="response")
abundx <- exp(predict(gam_E_N, dat_proj, type="response"))
dat_proj$gam_E <- presx * abundx

# spatial only models ("S")
gam_S_P <- gam(pres ~ s(lat,lon), data=dat, family=binomial)
#plot(gam_S_P, pages=1)
gam_S_N <- gam(log_abundance ~ s(lat,lon), data=dat[dat$abundance>0,], family=gaussian)
#plot(gam_S_N, pages=1)

presx <- predict(gam_S_P, dat, type="response")
abundx <- exp(predict(gam_S_N, dat, type="response"))
dat$gam_S <- presx * abundx

presx <- predict(gam_S_P, dat_proj, type="response")
abundx <- exp(predict(gam_S_N, dat_proj, type="response"))
dat_proj$gam_S <- presx * abundx


# environment and space models ("ES")
gam_ES_P <- gam(update(delta1_formula, ~. + s(lat,lon)), data=dat, family=binomial)
#plot(gam_ES_P)
gam_ES_N <- gam(update(delta2_formula, ~. + s(lat,lon)), data=dat[dat$abundance>0,], family=gaussian)
#plot(gam_ES_N)

presx <- predict(gam_ES_P, dat, type="response")
abundx <- exp(predict(gam_ES_N, dat, type="response"))
dat$gam_ES <- presx * abundx

presx <- predict(gam_ES_P, dat_proj, type="response")
abundx <- exp(predict(gam_ES_N, dat_proj, type="response"))
dat_proj$gam_ES <- presx * abundx

# environment, space, and time models ("EST")
gam_EST_P <- gam(update(delta1_formula, ~. + te(lat,lon,year)), data=dat, family=binomial)
#plot(gam_EST_P, pages=1)
gam_EST_N <- gam(update(delta2_formula, ~. + te(lat,lon,year)), data=dat[dat$abundance>0,], family=gaussian)
#plot(gam_EST_N, pages=1)

presx <- predict(gam_EST_P, dat, type="response")
abundx <- exp(predict(gam_EST_N, dat, type="response"))
dat$gam_EST <- presx * abundx

presx <- predict(gam_EST_P, dat_proj, type="response")
abundx <- exp(predict(gam_EST_N, dat_proj, type="response"))
dat_proj$gam_EST <- presx * abundx

dat_proj <- dat_proj %>% 
  pivot_longer(cols=c(gam_E, gam_S, gam_ES, gam_EST), names_to="model", values_to="pred") %>% 
  mutate(pred = as.numeric(pred))

dat_proj_summary <- dat_proj %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(year, lat_floor, model) %>% 
  summarise(mean_pred = mean(pred) * meanpatcharea) 

write_csv(dat_proj, file=here("results","GAM_pred_full.csv"))
write_csv(dat_proj_summary, file=here("results","GAM_pred_by_patch.csv"))

gg_gam_E_tile <- dat_proj_summary %>% 
  filter(model=="gam_E") %>% 
  ggplot(aes(x=year, y=lat_floor, fill=mean_pred)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Environment Only GAM", x="Year",y="Latitude")

gg_gam_ES_tile <- dat_proj_summary %>% 
  filter(model=="gam_ES") %>% 
  ggplot(aes(x=year, y=lat_floor, fill=mean_pred)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Environment + Space", x="Year",y="Latitude")

gg_gam_EST_tile <- dat_proj_summary %>% 
  filter(model=="gam_EST") %>% 
  ggplot(aes(x=year, y=lat_floor, fill=mean_pred)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Environment + Space + Time", x="Year",y="Latitude")

gg_gam_S_tile <- dat_proj_summary %>% 
  filter(model=="gam_S") %>% 
  ggplot(aes(x=year, y=lat_floor, fill=mean_pred)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Space Only", x="Year",y="Latitude")

ggsave(gg_gam_S_tile, filename=here("results","GAM_S_proj_tileplot.png"), scale=0.9, width=6, height=5)
ggsave(gg_gam_E_tile, filename=here("results","GAM_E_proj_tileplot.png"), scale=0.9, width=6, height=5)
ggsave(gg_gam_ES_tile, filename=here("results","GAM_ES_proj_tileplot.png"), scale=0.9, width=6, height=5)
ggsave(gg_gam_EST_tile, filename=here("results","GAM_EST_proj_tileplot.png"), scale=0.9, width=6, height=5)


# average predicted hauls within a patch (the same way we pre-processed the data for the DRM)
# multiply by patch area to be in DRM units
# from there calculate summary statistics about abundance and range position 