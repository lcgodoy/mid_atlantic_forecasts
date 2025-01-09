
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

