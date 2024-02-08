
set.seed(42)
library(tidyverse)
library(ggridges)
library(tidybayes)
#library(Cairo)
library(here)
library(magrittr)
library(rstan)
library(Matrix)
library(rstanarm)
library(cmdstanr)
library(data.table)

rstan_options(javascript = FALSE, auto_write = TRUE)

# load fit_drm function 
# fit_drm() fits the model and writes out the model object and a plot to the results directory
funs <- list.files("functions")
sapply(funs, function(x)
  source(file.path("functions", x)))

# which range edges should be calculated?
quantiles_calc <- c(0.05, 0.5, 0.95)
ctrl_file <- read_csv("control_file.csv") 
# ctrl_file <- read_csv("control_file.csv") %>%
#   filter(eval_l_comps==0,
#          spawner_recruit_relationship==1,
#          process_error_toggle==1,
#          known_f==1
#          )
fit_drms <- TRUE
make_plots <- FALSE
write_summary <- TRUE
iters <- 5000
warmups <- 2000
chains <- 4
cores <- 4

for(k in 1:nrow(ctrl_file)){
  i = ctrl_file$id[k]  
  
  results_path <- file.path("results",i)
  
  # turn off if you just want to load already-fitted models and analyze them
  
  if (fit_drms==TRUE){
    drm_fits <-  ctrl_file %>%
      filter(id == i)
    
    drm_fits$fits <- list(tryCatch(fit_drm(
      amarel = FALSE,
      create_dir = TRUE,
      run_name = drm_fits$id,
      do_dirichlet = drm_fits$do_dirichlet,
      eval_l_comps = drm_fits$eval_l_comps,
      T_dep_movement = drm_fits$T_dep_movement,
      T_dep_mortality = drm_fits$T_dep_mortality,
      T_dep_recruitment = drm_fits$T_dep_recruitment,
      spawner_recruit_relationship = drm_fits$spawner_recruit_relationship,
      process_error_toggle = drm_fits$process_error_toggle,
      exp_yn = drm_fits$exp_yn,
      known_f = drm_fits$known_f,
      known_historic_f = drm_fits$known_historic_f,
      warmup = warmups,
      iter = iters,
      chains = chains,
      cores = cores,
      adapt_delta = 0.99, 
      run_forecast = 1,
      quantiles_calc = quantiles_calc, 
    )
    ) 
    )# as currently written this just adds a column to drm_fits that says "all done". the column `fits` used to contain the model object itself 
    
    
  } # close fit_drms 
  
  # would be more efficient to do this above without having to load in the model 
  if(write_summary==TRUE){
    
    tmp_model <-  tryCatch(read_rds(file.path(results_path, "stan_model_fit.rds")))
    # may want to add in rhat metrics later, but those are per parameter
    #   https://mc-stan.org/cmdstanr/articles/cmdstanr.html 
    
    diagnostic_ls <- tryCatch(c(list(num_chains = chains, num_cores = cores, num_iters = iters, num_warmups = warmups), tmp_model$diagnostic_summary()))
    tryCatch(saveRDS(diagnostic_ls, file = file.path(results_path,"diagnostics.rds"))    )
    
  } # close summary 
} # close loop over models
# process_stanfit <- function(run_name){
#   
#   results_path <- here::here("results", run_name)
#   
#   tmp_model <-  read_rds(file.path(results_path, "stan_model_fit.rds"))
#   
#   
#   
#   #### insert things to do  ###
#   # tidybayes is the easiest way to get things out of a cmdstan object
#   out <- list()
#   out$diagnostics <- tmp_model$diagnostic_summary()
#   
#   out$n_at_age_hat <- tidybayes::gather_draws(tmp_model, n_at_age_hat[year,patch, age],n = 100)
#   
#   rm(tmp_model)
#   
#   gc()
#   
#   return(out)
# }
# 
# drm_fits <- ctrl_file %>%
#   mutate(fits = map(id, process_stanfit, .progress = TRUE))
# 
# drm_fits$fits[[1]]$diagnostics


# load results and do something

# diagnostic_fit <- drm_fits$fits[[which(drm_fits$id == i)]]
# rstan::check_hmc_diagnostics(tmp_model)
# 
# num_divergent <- get_num_divergent(tmp_model)
# num_iters <- length(get_divergent_iterations(diagnostic_fit))
# num_max_treedepth <- get_num_max_treedepth(diagnostic_fit)
# bfmi <- get_bfmi(diagnostic_fit)
# rhat <- rstan::Rhat(diagnostic_fit)
#  bulk_ess <- rstan::ess_bulk(diagnostic_fit)
# diagnostic_ls <- list(num_divergent=num_divergent, num_iters=num_iters, num_max_treedepth=num_max_treedepth, bfmi=bfmi)
# capture.output(diagnostic_ls, file = file.path(results_path,"diagnostics.txt"))

# load results and make plots

#   if(make_plots==TRUE){
#     # process results ---------------------------------------------------------
#     
#     load(here("processed-data","stan_data_prep.Rdata"))
#     
#     # visualize abundance over time 
#     abund_p_y <-  dat_train_dens %>%
#       mutate(abundance = mean_dens * meanpatcharea) 
#     
#     abund_p_y_hat <- tidybayes::spread_draws(diagnostic_fit, density_hat[patch,year])
#     
#     abundance_v_time <- abund_p_y_hat %>% 
#       ggplot(aes(year, density_hat)) + 
#       stat_lineribbon() +
#       geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
#       facet_wrap(~patch, scales = "free_y") +
#       labs(x="Year",y="Abundance", fill="Probability") + 
#       scale_fill_brewer()
#     
#     abund_p_y_proj <-  dat_test_dens %>%
#       mutate(abundance = mean_dens * meanpatcharea) 
#     
#     observed_abund_posterior_predictive <- tidybayes::spread_draws(diagnostic_fit, density_obs_proj[patch,year])
#     
#  #    abund_posterior_predictive <- tidybayes::spread_draws(diagnostic_fit, density_proj[patch,year])
#     
#     observed_abundance_forecast <- observed_abund_posterior_predictive %>% 
#       ggplot(aes(year, density_obs_proj)) + 
#       stat_lineribbon() +
#       geom_point(data = abund_p_y_proj, aes(year, abundance), color = "red") +
#       scale_x_continuous(breaks=seq(0, 10, 2)) +
#       facet_wrap(~patch, scales = "free_y", ncol=5) +
#       labs(x="Year",y="Abundance", fill="Probability") + 
#       scale_fill_brewer()
#     
#     # abundance_forecast <- abund_posterior_predictive %>%
#     #   ggplot(aes(year, density_proj)) +
#     #   stat_lineribbon() +
#     #   geom_point(data = abund_p_y_proj, aes(year, abundance), color = "red") +
#     #   facet_wrap(~patch, scales = "free_y", ncol=5) +
#     #   labs(x="Year",y="Abundance") +
#     #   scale_fill_brewer()
#     
#     n_p_l_y <- as.data.table(len) %>% 
#       rename(year = V3, patch = V1, length = V2) %>% 
#       group_by(year, patch) %>% 
#       mutate(plength = value / sum(value, na.rm = TRUE)) %>% 
#       ungroup() %>% 
#       filter(year == min(year) | year == max(year))
#     
#     observed_abundance_tile <- abund_p_y %>% 
#       mutate(Year = (year + min(years) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
#       ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
#       geom_tile() +
#       theme_bw() +
#       scale_x_continuous(breaks=seq(min(years), max(years), 4)) +
#       scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
#       scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
#       labs(title="Observed")
#     
#     proj_observed_abundance_tile <- abund_p_y_proj %>% 
#       mutate(Year = (year + min(years_proj) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
#       ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
#       geom_tile() +
#       theme_bw() +
#       scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
#       scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
#       scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
#       labs(title="Observed")
#     
#     proj_est_patch_abund <- observed_abund_posterior_predictive %>% 
#       group_by(year, patch) %>% 
#       summarise(abundance = mean(density_obs_proj))
#     
#     estimated_abundance_tile <- abund_p_y_hat%>% 
#       group_by(year, patch) %>% 
#       summarise(abundance = mean(density_hat)) %>% 
#       mutate(Year = (year + min(years) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
#       ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
#       geom_tile() +
#       theme_bw() +
#       scale_x_continuous(breaks=seq(min(years), max(years), 4)) +
#       scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
#       scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
#       labs(title="Estimated")
#     
#     proj_estimated_abundance_tile <- proj_est_patch_abund %>% 
#       mutate(Year = (year + min(years_proj) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
#       ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
#       geom_tile() +
#       theme_bw() +
#       scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
#       scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
#       scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
#       labs(title="Estimated")
#     
#     # other parameters 
#     raw <-  rstan::extract(diagnostic_fit,"raw")[[1]]
# 
#     gg_raw <- raw %>% 
#       as.data.frame() %>% 
#       pivot_longer(cols = everything(), names_to="year", values_to="value") %>% 
#       mutate(year = gsub("V", "", year)) %>% 
#       ggplot(aes(x=value, y=as_factor(year))) +
#       geom_density_ridges(scale=4) +
#       theme_ridges() +
#       labs(x="Estimate", y="Year", title="Raw") +
#       coord_cartesian(clip = "off") 
#     
#    gg_params_small <- plot(diagnostic_fit, pars=c('sigma_obs','d','width','alpha','beta_obs','theta_d','beta_obs_int','sigma_r_raw','beta_t','beta_rec','p_length_50_sel'))
#    
#    gg_params_large <- plot(diagnostic_fit, pars=c('sel_delta','log_mean_recruits','log_r0','Topt'))
#     
#     # length frequency 
#     n_at_length_hat <- rstan::extract(diagnostic_fit,"n_at_length_hat")[[1]]
#     
#     
#     
#     n_at_length_hat <- n_at_length_hat[sample(1:dim(n_at_length_hat)[1], 100, replace = TRUE),c(1,35),,]
#     
#     n_at_length_hat <- data.table::as.data.table(n_at_length_hat) %>% 
#       rename(year = V3, patch = V2, length = V1) %>% 
#       group_by(year, patch,iterations) %>% 
#       mutate(plength = value / sum(value)) %>% 
#       ungroup()
#     
#     n_at_length_hat$year[n_at_length_hat$year ==  max(n_at_length_hat$year)] <-  max(n_p_l_y$year)
#     
#     gg_length <- n_at_length_hat %>% 
#       ggplot(aes(length, plength)) +
#       stat_lineribbon(size = .1) +
#       geom_point(data = n_p_l_y, aes(length, plength), color = "red", alpha = 0.25) +
#       facet_grid(patch ~ year, scales = "free_y") + 
#       scale_x_continuous(limits = c(0, 50)) # generates warnings because of the close-to-zero probabilities at larger lengths
#     
#     # # range edges and centroids 
#     centroid_proj <- tidybayes::spread_draws(diagnostic_fit, centroid_proj[year]) 
#     
#     range_quantiles_proj <- tidybayes::spread_draws(diagnostic_fit, range_quantiles_proj[quantile, year]) %>%
#       mutate(quantile = as.factor(quantiles_calc[quantile]), .keep="unused") 
#     
#     # boop <- dat_test_dens %>% 
#     #   group_by(year) %>%
#     #   mutate(sumdens = sum(mean_dens)) %>%
#     #   arrange(patch) %>%
#     #   mutate(csum_dens = cumsum(mean_dens) / sumdens, 
#     #          diff_p = csum_dens - 0.05) %>%
#     #   group_by(year) %>% 
#     #   filter(diff_p == min(abs(diff_p)))
#     
#     
#     
#     # centroid position by year
#     dat_centroid_proj <- abund_p_y_proj %>%
#       group_by(year) %>%
#       summarise(centroid_lat = weighted.mean(x=patch, w=abundance) + (min(patches)-1))
#     
#     # save plots of model fit 
#     ggsave(observed_abundance_tile, filename=file.path(results_path, "observed_abundance_v_time_tileplot.png"), scale=0.9, width=6, height=5)
#     ggsave(estimated_abundance_tile, filename=file.path(results_path,"estimated_abundance_v_time_tileplot.png"), scale=0.9, width=6, height=5)
#     ggsave(abundance_v_time, filename=file.path(results_path,"abundance_est_v_time_by_patch.png"), dpi=300, width=10, height=5)
#     
#     # save plots of model forecast 
#     ggsave(proj_estimated_abundance_tile, filename=file.path(results_path,"proj_estimated_abundance_v_time_tileplot.png"), scale=0.9, width=6, height=5)
#     ggsave(proj_observed_abundance_tile, filename=file.path(results_path,"proj_observed_abundance_v_time_tileplot.png"), scale=0.9, width=6, height=5)
#     ggsave(observed_abundance_forecast, filename=file.path(results_path,"abundance_est_v_time_by_patch_proj.png"), dpi=300, width=10, height=5)
#     
#     # save plots of parameter posteriors 
#     ggsave(gg_length, filename=file.path(results_path,"length_dist_first_last_year.png"), dpi=300, width=5, height=10)
#     ggsave(gg_raw, filename=file.path(results_path,"raw_posterior.png"), dpi=300, width=5, height=10)
#     ggsave(gg_params_large, filename=file.path(results_path,"param_posteriors_large.png"), dpi=300, width=6, height=5)
#     ggsave(gg_params_small, filename=file.path(results_path,"param_posteriors_small.png"), dpi=300, width=6, height=5)
#     
#     # write out data for summary stats
#     write_csv(range_quantiles_proj, file=file.path(results_path,"range_quantiles_proj.csv"))
#     write_csv(centroid_proj, file=file.path(results_path,"centroid_proj.csv"))
#     write_csv(observed_abund_posterior_predictive, file=file.path(results_path,"density_obs_proj.csv"))
#     
#     
#   } # close make plots 
#   
# } # close for loop

