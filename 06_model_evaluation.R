###############
# SETUP
###############

# much of this script is copied from prep_summer_flounder.R; revisit it if that script is changed! 

# load packages 
set.seed(42)
library(tidyverse)
library(magrittr)
library(here)
#library(Hmisc) # for wtd.quantile() 
library(mgcv)
library(tidybayes)
library(ggrepel)
# read in data 
convergence_checks <- read_csv(file=here("results","convergence_checks.csv"))
ctrl_file <- read_csv(file=here("control_file.csv"))
dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))
dat_test <- read_csv(here("processed-data","flounder_catch_at_length_fall_testing.csv"))
load(here("processed-data","stan_data_prep.Rdata"))
quantiles_calc <- c(0.05, 0.5, 0.95)

# make user decisions 
divergence_cutoff <- 0.05
chains_cutoff <- 3 
generate_exploratory_plots <- FALSE

# identify models that pass convergence checks 
summarydat <- convergence_checks %>% 
  left_join(ctrl_file) %>% 
  filter(mean_divergences <= divergence_cutoff, 
         successful_chains >= chains_cutoff)

# write a function to calculate patch edge positions 
calculate_range_edge <- function(patches, weights, q){
  if(length(patches) == length(weights)){
    cutoff <- q * sum(weights)
    csum <- cumsum(weights)
    p <- min(which(csum > cutoff)) 
    if(p == 1) { # special case when range edge is in the first patch 
      out <- patches[1] + cutoff / weights[1]
    } else {
      dec <- (cutoff - csum[p-1]) / weights[p] # calculate decimal position into the edge patch -- what proportion of the patch would be occupied, assuming a constant density 
      out <- patches[p] + dec 
      return(out) }
  } else {return(print("Length of patches and weights must be equal!"))}
}

# what actually happened to edge and centroid positions in the testing data?
# need to keep aggregated at patch scale so we are comparing evenly across models 
dat_test_patch <- dat_test_dens %>% 
  group_by(year) %>% 
  summarise(
    warm_edge = calculate_range_edge(patches=lat_floor, weights=mean_dens, q=0.05),
    centroid = weighted.mean(lat_floor, w=mean_dens),
    cold_edge = calculate_range_edge(patches=lat_floor, weights=mean_dens, q=0.95),
    abund = sum(mean_dens) * meanpatcharea) %>% 
  arrange(year) %>% 
  mutate(year = year + min(years_proj) - 1, 
         abund_lr = log(abund / lag(abund))) %>% 
  pivot_longer(cols=c(warm_edge:abund_lr), names_to="feature", values_to="value") 

###############
# FIT SDMS FOR COMPARISON
###############

# prep data for fitting GAM 

# fit GAM to all hauls (not patch level)
dat_train_gam <- dat %>% 
  group_by(haulid, year, lat, btemp) %>% 
  summarise(dens = sum(number_at_length))

dat_test_gam <- dat_test %>% 
  group_by(haulid, year, lat, btemp) %>% 
  summarise(dens = sum(number_at_length))

# check data coverage of testing data
if(generate_exploratory_plots==TRUE){
  dat_test_gam %>% 
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
    geom_tile(aes(x=year, y=lat_floor, fill=prop_btemp, color=prop_btemp)) }
# because so much environmental data is missing in 2008, drop it from predictions 

dat_test_gam %<>% 
  filter(!year==2008)

# fit GAMs 
# see https://github.com/pinskylab/project_velocity/blob/master/6_model_fitting_loop.R#L182

spdata <- dat_train_gam %>%
  mutate(pres = ifelse(dens>0, 1, 0), 
         logdens = log(dens)) %>% 
  filter(!is.na(btemp))
spdata_proj <- dat_test_gam %>%
  mutate(pres = ifelse(dens>0, 1, 0), 
         logdens = log(dens)) %>% 
  filter(!is.na(btemp))

mypresmodtt<-formula(pres ~ s(btemp))
myabunmodtt<-formula(logdens ~ s(btemp))

# not sure if these are correct relative to how the source code worked, need to check 
gammaPA <- log(nrow(spdata)) / 2
gammaAbun <- log(nrow(spdata[spdata$pres==1,])) / 2

mygam1tt<-gam(mypresmodtt, family="binomial",data=spdata, select=TRUE, gamma=gammaPA) 
mygam2tt<-gam(myabunmodtt, data=spdata[spdata$pres==1,], select=TRUE, gamma=gammaAbun) 

preds1tt <- predict(mygam1tt, newdata = spdata_proj, type="response") 
preds2tt <- exp(predict(mygam2tt, newdata = spdata_proj, type='response'))

predstt <- preds1tt*preds2tt 
predstt[predstt<0] = 0
predstt[is.na(predstt)] = 0 # is this correct?

spdata_proj$predstt <- predstt

# calculate residuals by feature (centroid, edges) of forecast 
gam_summary <- spdata_proj %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor, year) %>% 
  summarise(dens_pred = mean(exp(predstt))) %>% # aggregate to patch scale for comparison to DRM 
  group_by(year) %>% # calculate summary stats 
  summarise(
    warm_edge = calculate_range_edge(patches=lat_floor, weights=dens_pred, q=0.05),
    centroid = weighted.mean(lat_floor, w=dens_pred),
    cold_edge = calculate_range_edge(patches=lat_floor, weights=dens_pred, q=0.95),
    abund = sum(dens_pred) * meanpatcharea) %>% 
  arrange(year) %>% 
  mutate(abund_lr = log(abund / lag(abund))) %>% # note that this is technically the LR of change from 2007 to 2009 because the GAM doesn't use 2008
  pivot_longer(cols=warm_edge:abund_lr, names_to="feature", values_to="value_tmp") %>% 
  left_join(dat_test_patch)%>% # compare  to true data 
  mutate(resid = value_tmp - value, 
         resid_sq = resid^2, 
         .keep = "unused", # drop all the columns used in calculations 
         id = "GAM") 

gam_time <- spdata_proj %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor, year) %>% 
  summarise(dens_pred = mean(exp(predstt))) %>% # aggregate to patch scale for comparison to DRM 
  group_by(year) %>% # calculate summary stats 
  summarise(
    warm_edge = calculate_range_edge(patches=lat_floor, weights=dens_pred, q=0.05),
    centroid = weighted.mean(lat_floor, w=dens_pred),
    cold_edge = calculate_range_edge(patches=lat_floor, weights=dens_pred, q=0.95),
    abund = sum(dens_pred) * meanpatcharea) %>% 
  arrange(year) %>% 
  mutate(abund_lr = log(abund / lag(abund))) %>% # note that this is technically the LR of change from 2007 to 2009 because the GAM doesn't use 2008
  pivot_longer(cols=warm_edge:abund_lr, names_to="feature", values_to="value_tmp") 

###############
# MAKE PERSISTENCE FORECAST
###############

persistence <- dat_train_dens %>% 
  filter(year == max(year)) %>% 
  summarise(
    warm_edge = calculate_range_edge(patches=lat_floor, weights=mean_dens, q=0.05),
    centroid = weighted.mean(lat_floor, w=mean_dens),
    cold_edge = calculate_range_edge(patches=lat_floor, weights=mean_dens, q=0.95),
    abund = sum(mean_dens) * meanpatcharea) 

persistence_dat <- data.frame(year = years_proj) %>% 
  bind_rows(persistence) %>% 
  fill(warm_edge, centroid, cold_edge, abund, .direction="up") %>% 
  filter(!is.na(year))

persistence_summary <- persistence_dat %>% 
  pivot_longer(cols=warm_edge:abund, names_to="feature", values_to="value_tmp") %>% 
  left_join(dat_test_patch)%>% # compare  to true data 
  mutate(resid = value_tmp - value, 
         resid_sq = resid^2, 
         .keep = "unused", # drop all the columns used in calculations 
         id = "Persistence") 

###############
# SUMMARIZE DRM OUTPUTS
###############

# generate summary stats for the successful DRMs

drm_out <- NULL

# pull in drm forecasts and calculate residuals by year and by feature (centroid, edges) 
# slow! try to speed this up at some point
for(i in 1:nrow(summarydat)){
  
  tmpdat <- summarydat[i,]
  
  results_path <- file.path(paste0('~/github/mid_atlantic_forecasts/results/',tmpdat$id))
  
  # get the Stan model and extract posteriors that we want for plots 
  tmp_model <-  tryCatch(read_rds(file.path(results_path, "stan_model_fit.rds")))
  observed_abund_posterior_predictive <- tidybayes::spread_draws(tmp_model, density_obs_proj[patch,year])
  centroid_proj <- tidybayes::spread_draws(tmp_model, centroid_proj[year]) 
  range_quantiles_proj <- tidybayes::spread_draws(tmp_model, range_quantiles_proj[quantile, year]) %>%
    mutate(quantile = as.factor(quantiles_calc[quantile]), .keep="unused") 
  
  # save those posteriors
  write_rds(observed_abund_posterior_predictive, file.path(results_path, "density_obs_proj.rds"))
  write_rds(range_quantiles_proj, file.path(results_path, "range_quantiles_proj.rds"))
  write_rds(centroid_proj, file.path(results_path, "centroid_proj.rds")) 
  
  # abund_tmp <- observed_abund_posterior_predictive %>%
  #   mutate(year = year + min(years_proj) - 1) %>%
  #   group_by(year, .chain, .iteration, .draw) %>%
  #   summarise(density_obs_proj_sum = sum(density_obs_proj)) %>%
  #   left_join(dat_test_patch %>% filter(feature=="centroid"))
  
  centroid_tmp <- centroid_proj %>% 
    mutate(year = year + min(years_proj) - 1) %>% 
    left_join(dat_test_patch %>% filter(feature=="centroid")) %>% 
    mutate(resid = centroid_proj - value,
           resid_sq = resid^2, 
           id = tmpdat$id) %>% 
    select(-centroid_proj, -value)
  
  warm_edge_tmp <- range_quantiles_proj %>% 
    filter(range_quantiles_proj < Inf,
           quantile == 0.05) %>%  
    mutate(
      year = year + min(years_proj) - 1, 
      range_quantiles_proj = range_quantiles_proj + min(patches)) %>% # don't need to subtract 1. we want patch 0.09 to be lat 35.09
    left_join(dat_test_patch %>% filter(feature=="warm_edge")) %>% 
    mutate(resid = range_quantiles_proj - value,
           resid_sq = resid^2, 
           id = tmpdat$id)%>% 
    select(-range_quantiles_proj, -value, -quantile)
  
  cold_edge_tmp <- range_quantiles_proj %>% 
    filter(range_quantiles_proj < Inf,
           quantile == 0.95) %>%  
    mutate(
      year = year + min(years_proj) - 1, 
      range_quantiles_proj = range_quantiles_proj + min(patches)) %>% # don't need to subtract 1. we want patch 0.09 to be lat 35.09
    left_join(dat_test_patch %>% filter(feature=="cold_edge")) %>% 
    mutate(resid = range_quantiles_proj - value,
           resid_sq = resid^2, 
           id = tmpdat$id)%>% 
    select(-range_quantiles_proj, -value, -quantile)
  
  drm_out <- bind_rows(drm_out, centroid_tmp, warm_edge_tmp, cold_edge_tmp)
}
write_csv(drm_out, file=here("processed-data","posteriors_for_model_evaluation.csv"))

drm_out <- drm_out %>% 
  filter(id %in% summarydat$id)

drm_summary <- drm_out %>% 
  filter(!is.na(resid)) %>% 
  group_by(year, feature, id) %>% 
  summarise(resid = mean(resid)) %>% 
  mutate(resid_sq = resid^2)

dat_forecasts <- drm_summary %>% 
  left_join(ctrl_file %>% select(id)) %>% 
  bind_rows(gam_summary, persistence_summary)%>% 
  filter(feature %in% c('centroid','cold_edge','warm_edge')) 

# pool across years to calculate bias and RMSE
dat_forecasts_summ <- dat_forecasts %>% 
  group_by(feature, id) %>% 
  summarise(RMSE = sqrt(mean(resid_sq)), 
            Bias = mean(resid)) %>% 
  pivot_longer(cols=c(RMSE, Bias), values_to="value", names_to="metric") 

write_csv(dat_forecasts_summ, file = here("processed-data","model_comparison_summary.csv"))

#####
# plots 
#####

gg_metrics <- dat_forecasts_summ  %>% 
  mutate(feature = case_match(feature, "centroid" ~ "Centroid", "warm_edge" ~ "Warm Edge", "cold_edge" ~ "Cold Edge", .default=feature),
         metric = case_match(metric, "rmse" ~ "RMSE", "bias" ~ "Bias", .default=metric), 
         type = ifelse(str_detect(id, "0."),"DRM", id),
         metric = factor(metric, levels = c("RMSE","Bias")))  %>%  
  ggplot(aes(x=feature, y=value)) + 
  geom_point(aes(color=type, fill=type),size=1) +
  geom_text_repel(aes(label=id), hjust = 1, max.overlaps = Inf) +
  scale_color_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
  scale_fill_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
  theme_bw() + 
  facet_wrap(~metric, scales="free_y") +
  labs(x="Feature", y="Metric") + 
  theme(legend.position = "bottom")
gg_metrics

gg_metrics2 <- dat_forecasts_summ  %>% 
  mutate(feature = case_match(feature, "centroid" ~ "Centroid", "warm_edge" ~ "Warm Edge", "cold_edge" ~ "Cold Edge", .default=feature),
         type = ifelse(str_detect(id, "0."),"DRM", id),
         metric = factor(metric, levels = c("RMSE","Bias")),
         id = gsub("v0.", "v", id))  %>%  
  pivot_wider(names_from = metric, values_from = value) %>% 
  ggplot(aes(x=RMSE, y=Bias)) + 
  geom_point(aes(color=type, fill=type),size=1) +
  geom_text_repel(aes(label=id), hjust = 1, max.overlaps = Inf) +
  scale_color_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
  scale_fill_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
  theme_bw() + 
  facet_wrap(~feature) +
  theme(legend.position = "bottom",
        legend.title=element_blank())
gg_metrics2
ggsave(gg_metrics2, filename=here("results","bias_v_rmse.png"), width=110, height=60, dpi=600, units="mm", scale=1.5)

# still too much to look at... need to trim down more 

# hacky way to get the best models 
best_drms <- dat_forecasts_summ %>% 
  filter(metric=='RMSE', !id %in% c('GAM','Persistence')) %>% 
  group_by(id) %>% 
  mutate(mean_value = mean(value)) %>% 
  filter(mean_value < 0.75)
length(unique(best_drms$id))

gg_metrics3 <- dat_forecasts_summ  %>% 
  filter(id %in% c(unique(best_drms$id), 'GAM','Persistence')) %>% 
  mutate(feature = case_match(feature, "centroid" ~ "Centroid", "warm_edge" ~ "Warm Edge", "cold_edge" ~ "Cold Edge", .default=feature),
         type = ifelse(str_detect(id, "0."),"DRM", id),
         metric = factor(metric, levels = c("RMSE","Bias")),
         id = gsub("v0.", "v", id))  %>%  
  pivot_wider(names_from = metric, values_from = value) %>% 
  ggplot(aes(x=RMSE, y=Bias)) + 
  geom_point(aes(color=type, fill=type),size=1) +
  geom_text_repel(aes(label=id), hjust = 1, max.overlaps = Inf) +
  scale_color_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
  scale_fill_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
  theme_bw() + 
  facet_wrap(~feature) +
  theme(legend.position = "bottom",
        legend.title=element_blank())
gg_metrics3
ggsave(gg_metrics3, filename=here("results","bias_v_rmse.png"), width=110, height=60, dpi=600, units="mm", scale=1.5)

ctrl_dat <- ctrl_file %>% 
  select(-id, -do_dirichlet, -exp_yn, -known_historic_f, -description) 

ctrl_dat <- cbind(id = ctrl_file$id, data.frame(ifelse(ctrl_dat == 0, "No", "Yes")))  # hacky way to convert 0s and 1s into yeses and nos 

drms_without_t <- ctrl_file %>% 
  mutate(tmp = T_dep_movement + T_dep_recruitment + T_dep_mortality) %>% 
  filter(tmp == 0) %>% 
  pull(id)

tbl_out <- ctrl_dat %>% 
  filter(id %in% best_drms$id) %>% 
  pivot_longer(cols = c('T_dep_movement','T_dep_mortality','T_dep_recruitment'), values_to='value',names_to='formulation') %>% 
  group_by(id) %>% 
  mutate(T_effect = case_when(
    value=='Yes' ~ formulation,
    id %in% drms_without_t ~ 'None'
    )) %>% 
  select(-value, -formulation) %>% 
  group_by(id) %>% 
  fill(T_effect, .direction="updown") %>% 
  distinct() %>% 
  mutate(T_effect = case_match(T_effect, 
    'T_dep_recruitment' ~ 'Recruitment',
    'T_dep_mortality' ~ 'Mortality',
    'T_dep_movement' ~ 'Movement',
    'None' ~ 'None'
  ))
  
colnames(tbl_out) <- c("ID","Fit to Length", "Spawner-Recruit Relationship","Process Error", "Known F","Temperature Effect")

write_csv(tbl_out, file = here("results","best_drm_table.csv"))

# want to plot dat_test_patch, gam_time, and persistence_dat against the posteriors 


points_for_plot <- dat_test_patch %>% 
  mutate(id = 'Observed')%>% rename(value_tmp = value) %>% 
  bind_rows(gam_time %>% mutate(id = 'GAM') )  %>% 
  bind_rows(persistence_dat %>% pivot_longer(cols=c('warm_edge','cold_edge','abund','centroid'), values_to='value_tmp', names_to='feature') %>% mutate(id='Persistence')) %>% 
  filter(feature %in% c('warm_edge','cold_edge','centroid')) %>% 
  mutate(feature = case_match(feature, "centroid" ~ "Centroid", "warm_edge" ~ "Warm Edge", "cold_edge" ~ "Cold Edge", .default=feature))

for(k in best_drms){
  results_path <- file.path(paste0('~/github/mid_atlantic_forecasts/results/',k))
  tmp_edges <- read_rds(file.path(results_path, "range_quantiles_proj.rds")) %>% 
    mutate(year = year + min(years_proj) - 1,
           range_quantiles_proj = range_quantiles_proj + min(patches) - 1) %>% 
    filter(range_quantiles_proj < Inf,
           year < 2017)
  tmp_centroids <- read_rds(file.path(results_path, "centroid_proj.rds"))%>% 
    mutate(year = year + min(years_proj) - 1) %>% 
    filter(year < 2017)
  
  gg_range_time_drm <- ggplot() +
    stat_lineribbon(data = tmp_edges[tmp_edges$quantile==0.05,], aes(x=year, y=range_quantiles_proj))+ 
    stat_lineribbon(data = tmp_edges[tmp_edges$quantile==0.95,], aes(x=year, y=range_quantiles_proj))+ 
    stat_lineribbon(data = tmp_centroids, aes(x=year, y=centroid_proj))+ 
    geom_point(data=points_for_plot %>% filter(id == "Observed") %>% rename("Feature" = feature), aes(x=year, y=value_tmp, shape=Feature), color="red", size=2) + 
    geom_line(data=points_for_plot %>% filter(id == "Observed")%>% rename("Feature" = feature), aes(x=year, y=value_tmp, group=Feature), color="red") +
    scale_fill_brewer(name = "Credible Interval") +
    scale_x_continuous(breaks =seq(2007, 2016, 1)) + 
    scale_y_continuous(breaks = seq(35, 44, 1), labels =  seq(35, 44, 1), limits=c(34, 44.5)) +
    labs(x="Year", y="Latitude", title = paste0("DRM ",k)) + 
    theme_bw()
  ggsave(gg_range_time_drm, filename=paste0(here("results"),"/range_time_",k,".png"), dpi=600, units="mm", width=130, height=75)
}

gg_range_time <-  ggplot() +
  geom_point(data=points_for_plot %>% rename("Feature" = feature), aes(x=year, y=value_tmp, color = id, shape=Feature), size=2) + 
  geom_line(data=points_for_plot %>% rename("Feature" = feature), aes(x=year, y=value_tmp, color=id, group = interaction(id, Feature))) +
  scale_x_continuous(breaks =seq(2007, 2016, 1)) + 
  scale_y_continuous(breaks = seq(35, 44, 1), labels =  seq(35, 44, 1), limits=c(34, 44.5)) +
  scale_color_manual(values = c('blue','red','green')) +
  labs(x="Year", y="Latitude", title = "Model Comparisons") + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  guides(color=guide_legend("Feature"), shape = "none") 
ggsave(gg_range_time, filename=here("results","range_time.png"), dpi=600, units="mm", width=115, height=75)

#     gg_length <- n_at_length_hat %>% 
#       ggplot(aes(length, plength)) +
#       stat_lineribbon(size = .1) +
#       geom_point(data = n_p_l_y, aes(length, plength), color = "red", alpha = 0.25) +
#       facet_grid(patch ~ year, scales = "free_y") + 
#       scale_x_continuous(limits = c(0, 50)) # generates warnings because of the close-to-zero probabilities at larger lengths

abund_p_y_proj <-  dat_test_dens %>%
  mutate(abundance = mean_dens * meanpatcharea)

gg_observed_abundance_tile <- abund_p_y_proj %>%
  mutate(Year = (year + min(years_proj) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>%
  ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Observed")

gg_real <- dat_test_patch %>% 
  filter(feature %in% c('centroid','cold_edge','warm_edge')) %>% 
  ggplot(aes(x=year, y=value )) + 
  geom_line() +
  geom_point() +
  theme_bw() + 
  labs(x="Year", y="Latitude") + 
  facet_wrap(~feature, ncol=1)
gg_real

gg_real_plus <- 
  
  gg_bias <- dat_forecasts %>% 
  ggplot(aes(x=year, y=resid,color=model_name, fill=model_name, shape = model_name )) + 
  geom_line() +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "bottom") +
  labs(x="Year", y="Residuals (° lat)") + 
  facet_wrap(~feature)
gg_bias
