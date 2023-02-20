args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)>1) {
  stop("Only one argument must be supplied", call.=FALSE)
}

# designed to be sent as parallel jobs via a bash script on a HPC like Amarel at Rutgers
# I have commented out or deleted some of the desireable features for a local machine, like analyzing an already-fitted model

dyn.load('/projects/community/gcc/9.2.0/gc563/lib64/libgfortran.so.5')

set.seed(42)
library(tidyverse)
library(tidybayes)
library(magrittr)
#library(here)
library(rstan)
library(Matrix)
library(rstanarm)
library(data.table)

# results.dir <- "/projects/f_mlp195/fredston/"
home_drm_path <- "/home/fredston/mid_atlantic_forecasts/"

rstan_options(javascript = FALSE, auto_write = TRUE)

# load fit_drm function 
# fit_drm() fits the model and writes out the model object and a plot to the results directory
funs <- list.files(paste0(home_drm_path, "functions"))
sapply(funs, function(x)
  source(file.path(paste0(home_drm_path, "functions"), x)))

# which range edges should be calculated?
quantiles_calc <- c(0.05, 0.5, 0.95)

ctrl_file <- read_csv(paste0(home_drm_path, "control_file.csv") )

fit_drms <- TRUE
make_plots <- TRUE

i = ctrl_file$id[args[1]]  

# results_path = paste0(results.dir,i,"/")
results_path = ""

# turn off if you just want to load already-fitted models and analyze them

# if (fit_drms){
# couldn't get this to work by creating a column with pmap(), got the uninformative warning "Problem with `mutate()` column" and it didn't work  
drm_fits <-  ctrl_file %>%
  filter(id == i) 

drm_fits$fits <- list(fit_drm(
  amarel = TRUE,
  run_name = drm_fits$id,
  results_path = results_path, 
  create_dir = FALSE,
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
  warmup = 10,
  iter = 50,
  chains = 4,
  cores = 4,
  run_forecast = 1,
  quantiles_calc = quantiles_calc, 
)
)


# } else {

#   drm_fits <- ctrl_file %>%
#     filter(id == i) %>% 
#     mutate(fits = pmap(list(run_name = id), ~ purrr::safely(readr::read_rds)(
#       here("results", .x, "stan_model_fit.rds")
#     )))
#   
#   #   drm_fits %<>% as.list()
#   
#   drm_worked <- map_lgl(map(drm_fits$fits,"error"), is.null)
#   
#   drm_fits <- drm_fits %>% 
#     filter(drm_worked) %>% 
#     mutate(fits = map(fits, "result"))
#   
# }

if(make_plots==TRUE){
  # process results ---------------------------------------------------------
  
  diagnostic_fit <- drm_fits$fits[[which(drm_fits$id == i)]]
  
  rstan::check_hmc_diagnostics(diagnostic_fit)
  
  num_divergent <- get_num_divergent(diagnostic_fit)
  num_iters <- length(get_divergent_iterations(diagnostic_fit))
  num_max_treedepth <- get_num_max_treedepth(diagnostic_fit)
  bfmi <- get_bfmi(diagnostic_fit)
  
  diagnostic_ls <- list(num_divergent=num_divergent, num_iters=num_iters, num_max_treedepth=num_max_treedepth, bfmi=bfmi)
  
  capture.output(diagnostic_ls, file = paste0(results_path,"diagnostics.txt"))
  
  load("/home/fredston/mid_atlantic_forecasts/processed-data/stan_data_prep.Rdata")
  
  # visualize abundance over time 
  abund_p_y <-  dat_train_dens %>%
    mutate(abundance = mean_dens * meanpatcharea) 
  
  abund_p_y_hat <- tidybayes::spread_draws(diagnostic_fit, density_hat[patch,year])
  
  abundance_v_time <- abund_p_y_hat %>% 
    ggplot(aes(year, density_hat)) + 
    stat_lineribbon() +
    geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
    facet_wrap(~patch, scales = "free_y") +
    labs(x="Year",y="Abundance", fill="Probability") + 
    scale_fill_brewer()
  
  abund_p_y_proj <-  dat_test_dens %>%
    mutate(abundance = mean_dens * meanpatcharea) 
  
  observed_abund_posterior_predictive <- tidybayes::spread_draws(diagnostic_fit, density_obs_proj[patch,year])
  
  abund_posterior_predictive <- tidybayes::spread_draws(diagnostic_fit, density_proj[patch,year])
  
  observed_abundance_forecast <- observed_abund_posterior_predictive %>% 
    ggplot(aes(year, density_obs_proj)) + 
    stat_lineribbon() +
    geom_point(data = abund_p_y_proj, aes(year, abundance), color = "red") +
    scale_x_continuous(breaks=seq(0, 10, 2)) +
    facet_wrap(~patch, scales = "free_y", ncol=5) +
    labs(x="Year",y="Abundance", fill="Probability") + 
    scale_fill_brewer()
  
  # abundance_forecast <- abund_posterior_predictive %>%
  #   ggplot(aes(year, density_proj)) +
  #   stat_lineribbon() +
  #   geom_point(data = abund_p_y_proj, aes(year, abundance), color = "red") +
  #   facet_wrap(~patch, scales = "free_y", ncol=5) +
  #   labs(x="Year",y="Abundance") +
  #   scale_fill_brewer()
  
  n_p_l_y <- as.data.table(len) %>% 
    rename(year = V3, patch = V1, length = V2) %>% 
    group_by(year, patch) %>% 
    mutate(plength = value / sum(value, na.rm = TRUE)) %>% 
    ungroup() %>% 
    filter(year == min(year) | year == max(year))
  
  observed_abundance_tile <- abund_p_y %>% 
    mutate(Year = (year + min(years) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
    ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
    geom_tile() +
    theme_bw() +
    scale_x_continuous(breaks=seq(min(years), max(years), 4)) +
    scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
    scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
    labs(title="Observed")
  
  proj_observed_abundance_tile <- abund_p_y_proj %>% 
    mutate(Year = (year + min(years_proj) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
    ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
    geom_tile() +
    theme_bw() +
    scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
    scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
    scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
    labs(title="Observed")
  
  proj_est_patch_abund <- observed_abund_posterior_predictive %>% 
    group_by(year, patch) %>% 
    summarise(abundance = mean(density_obs_proj))
  
  estimated_abundance_tile <- abund_p_y_hat%>% 
    group_by(year, patch) %>% 
    summarise(abundance = mean(density_hat)) %>% 
    mutate(Year = (year + min(years) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
    ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
    geom_tile() +
    theme_bw() +
    scale_x_continuous(breaks=seq(min(years), max(years), 4)) +
    scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
    scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
    labs(title="Estimated")
  
  proj_estimated_abundance_tile <- proj_est_patch_abund %>% 
    mutate(Year = (year + min(years_proj) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
    ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
    geom_tile() +
    theme_bw() +
    scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
    scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
    scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
    labs(title="Estimated")
  
  # other parameters 
  sigma_obs_hat <-  rstan::extract(diagnostic_fit,"sigma_obs")[[1]]
  
  theta_d <-  rstan::extract(diagnostic_fit,"sigma_obs")[[1]]
  
  # length frequency 
  n_at_length_hat <- rstan::extract(diagnostic_fit,"n_at_length_hat")[[1]]
  
  n_at_length_hat <- n_at_length_hat[sample(1:1000, 100, replace = FALSE),c(1,35),,]
  
  n_at_length_hat <- data.table::as.data.table(n_at_length_hat) %>% 
    rename(year = V3, patch = V2, length = V1) %>% 
    group_by(year, patch,iterations) %>% 
    mutate(plength = value / sum(value)) %>% 
    ungroup()
  
  n_at_length_hat$year[n_at_length_hat$year ==  max(n_at_length_hat$year)] <-  max(n_p_l_y$year)
  
  gg_length <- n_at_length_hat %>% 
    ggplot(aes(length, plength)) +
    stat_lineribbon(size = .1) +
    geom_point(data = n_p_l_y, aes(length, plength), color = "red", alpha = 0.25) +
    facet_grid(patch ~ year, scales = "free_y") + 
    scale_x_continuous(limits = c(0, 50)) # generates warnings because of the close-to-zero probabilities at larger lengths
  
  # # range edges and centroids 
  centroid_proj <- rstan::extract(diagnostic_fit,"centroid_proj")[[1]] %>%
    as.data.table() %>%
    pivot_longer(cols=everything(), names_to="year", values_to="lat") %>%
    mutate(year = as.numeric(str_replace(year, "V", "")))
  
  range_quantiles_proj <- rstan::extract(diagnostic_fit,"range_quantiles_proj")[[1]] %>%
    as.data.table() %>%
    rename(year=V1, patch=value, quantile_category = V2) %>%
    mutate(quantile = as.factor(quantiles_calc[quantile_category]), .keep="unused") %>%
    group_by(quantile, year) %>%
    summarize(lat = mean(patch))
  
  # centroid position by year
  dat_centroid_proj <- abund_p_y_proj %>%
    group_by(year) %>%
    summarise(centroid_lat = weighted.mean(x=patch, w=abundance) + (min(patches)-1))
  
  # centroid_plot <- ggplot() +
  #   geom_line(data=range_quantiles_proj, aes(x=year, y=lat, color=quantile)) + 
  #   geom_point(data=dat_centroid_proj, aes(x=year, y=centroid_lat), color="red")
  # 
  ggsave(observed_abundance_tile, filename=paste0(results_path,"observed_abundance_v_time_tileplot.pdf"), scale=0.9, width=6, height=5)
  ggsave(estimated_abundance_tile, filename=paste0(results_path,"estimated_abundance_v_time_tileplot.pdf"), scale=0.9, width=6, height=5)
  
  ggsave(proj_estimated_abundance_tile, filename=paste0(results_path,"proj_estimated_abundance_v_time_tileplot.pdf"), scale=0.9, width=6, height=5)
  ggsave(proj_observed_abundance_tile, filename=paste0(results_path,"proj_observed_abundance_v_time_tileplot.pdf"), scale=0.9, width=6, height=5)
  
  
  ggsave(observed_abundance_forecast, filename=paste0(results_path, "abundance_est_v_time_by_patch_proj.pdf"), dpi=300, width=10, height=5)
  ggsave(abundance_v_time, filename=paste0(results_path, "abundance_est_v_time_by_patch.pdf"), dpi=300, width=10, height=5)
  
  #  ggsave(abundance_forecast, filename=here(paste0("results/",i), "abundance_est_v_time_by_patch_proj.pdf"), dpi=300, width=10, height=5) 
  
  ggsave(gg_length, filename=paste0(results_path,"length_dist_first_last_year.pdf"), dpi=300, width=5, height=10)
  
  # write out data for summary stats
  
  write_csv(range_quantiles_proj, file=paste0(results_path,"range_quantiles_proj.csv"))
  write_csv(centroid_proj, file=paste0(results_path,"centroid_proj.csv"))
  write_csv(observed_abund_posterior_predictive, file=paste0(results_path,"density_obs_proj.csv"))
  
  
} # close make plots 
