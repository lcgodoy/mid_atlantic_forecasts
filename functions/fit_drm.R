fit_drm <- function(amarel = FALSE,
                    run_name = "test",
                    results_path = file.path("results", run_name),
                    create_dir = TRUE,
                    do_dirichlet = 1,
                    eval_l_comps = 0,
                    T_dep_mortality = 0,
                    T_dep_recruitment = 0,
                    T_dep_movement = 0,
                    spawner_recruit_relationship = 1,
                    run_forecast = 0,
                    process_error_toggle = 1,
                    exp_yn = 0,
                    warmup = 1000,
                    iter = 2000,
                    max_treedepth =  10,
                    chains =  1,
                    refresh = 10,
                    cores = 1,
                    adapt_delta = 0.95, 
                    drm_name = "process_sdm",
                    number_quantiles = 3,
                    quantiles_calc = c(0.05, 0.5, 0.95),
                    known_f = 0,
                    known_historic_f = 1,
                    sigma_obs_cv = 0.1,
                    h = 0.8,
                    dcap = 1 / 3) {
  if (amarel == TRUE) {
    dyn.load('/projects/community/gcc/9.2.0/gc563/lib64/libgfortran.so.5')
    
  }
  
  if (create_dir == TRUE) {
    if (!dir.exists(results_path)) {
      dir.create(results_path, recursive = TRUE)
    }
  }
  
  if (amarel == FALSE) {
    load(here::here("processed-data", "stan_data_prep.Rdata"))
  }
  
  if (amarel == TRUE) {
    load(
      '/home/fredston/mid_atlantic_forecasts/processed-data/stan_data_prep.Rdata'
    )
  }
  
  if (known_historic_f == 0) {
    f <- matrix(0.2, nrow = nrow(f), ncol = ncol(f))
    
  }
  
  if (known_f == 0) {
    f_proj <-
      matrix(rep(f[, ncol(f)], ncol(f_proj)),
             ncol = ncol(f_proj),
             nrow = nrow(f_proj))
    
  }
  
  if (T_dep_movement == 1) {
    dcap = 1
  }
  
  stan_data <- list(
    np = np,
    patches = patches,
    n_ages = n_ages,
    ny_train = ny,
    ny_proj = ny_proj,
    n_lbins = n_lbins,
    n_at_length = len,
    dens = dens,
    sbt = sbt,
    sbt_proj = sbt_proj,
    m = m,
    f = f,
    f_proj = f_proj,
    k = k,
    loo = loo,
    t0 = t0,
    cv = cv,
    length_50_sel_guess = length_50_sel_guess,
    age_sel = age_sel,
    bin_mids = bin_mids,
    sel_100 = sel_100,
    age_at_maturity = age_at_maturity,
    l_at_a_key = l_at_a_mat,
    wt_at_age = wt_at_age,
    do_dirichlet = do_dirichlet,
    eval_l_comps = eval_l_comps,
    # evaluate length composition data? 0=no, 1=yes
    T_dep_mortality = T_dep_mortality,
    # CURRENTLY NOT REALLY WORKING
    T_dep_recruitment = T_dep_recruitment,
    # think carefully before making more than one of the temperature dependencies true,
    T_dep_movement = T_dep_movement,
    spawner_recruit_relationship = spawner_recruit_relationship,
    run_forecast = run_forecast,
    exp_yn = exp_yn,
    process_error_toggle = process_error_toggle,
    number_quantiles = number_quantiles,
    quantiles_calc = quantiles_calc,
    sigma_obs_cv = sigma_obs_cv,
    h = h,
    dcap = dcap
  )
  nums <- 100 * exp(-.2 * (0:(n_ages - 1)))
  
  if (amarel == FALSE) {
    stan_file = here::here("src", paste0(drm_name, ".stan"))
  }
  
  if (amarel == TRUE) {
    stan_file = paste0('/home/fredston/mid_atlantic_forecasts/src/',
                       drm_name,
                       '.stan')
  }
  drm_model <- cmdstan_model( here::here("src","process_sdm.stan"))
  
  stan_model_fit = drm_model$sample(
    data = stan_data,
    chains = chains,
    iter_warmup = warmup,
    iter_sampling = iter - warmup,
    parallel_chains = cores,
    refresh = refresh,
    max_treedepth = max_treedepth,
    adapt_delta = 0.85,
    init = lapply(1:chains, function(x)
      list(
        Topt = jitter(12, 4),
        log_r0 = jitter(10, 5),
        beta_obs = jitter(1e-6, 4),
        beta_obs_int = jitter(-10, 2)
      ))
  )


stan_model_fit$save_object(file = file.path(results_path,
                                     "stan_model_fit.rds"))

# fit <- readRDS(file.path(results_path,
#                          "stan_model_fit.rds"))
# browser()
#   
#   stan_model_fit <- stan(
#     sample_file = file.path(results_path,paste0(run_name,".csv")),
#     diagnostic_file = file.path(results_path,paste0(run_name,"_diagnostics.csv")),
#     file = stan_file,
#     data = stan_data,
#     chains = chains,
#     warmup = warmup,
#     iter = iter,
#     cores = cores,
#     refresh = refresh,
#     control = list(max_treedepth = 8,
#                    adapt_delta = 0.95),
#     init = lapply(1:chains, function(x)
#       list(
#         Topt = jitter(12, 4),
#         log_r0 = jitter(10, 5),
#         beta_obs = jitter(1e-6, 4),
#         beta_obs_int = jitter(-10, 2)
#       ))
#   )
  rm(stan_model_fit)
  # readr::write_rds(stan_model_fit, file = file.path(results_path,
  #                                                 "stan_model_fit.rds"))
  # didn't work on HPC
  
  # readr::write_rds(stan_model_fit, file = file.path(results_path,
  #                                                   "stan_model_fit.rds"))
  # 
  # abund_p_y <- dat_train_dens %>%
  #   mutate(abundance = mean_dens * meanpatcharea)
  #
  # abund_p_y_hat <- tidybayes::spread_draws(stan_model_fit, density_hat[patch,year])
  #
  #
  # abundance_v_time <- abund_p_y_hat %>%
  #   ggplot(aes(year, density_hat)) +
  #   stat_lineribbon() +
  #   geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  #   facet_wrap(~patch, scales = "free_y") +
  #   labs(x="Year",y="Abundance") +
  #   scale_fill_brewer()
  # # abundance_v_time
  # ggsave(abundance_v_time, filename=file.path(results_path,"abundance_fits.pdf"), width=7, height=4)
  #
  gc()
  return("all done")
  
} # close fit_drm