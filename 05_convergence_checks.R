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

ctrl_file <- read_csv("control_file.csv") %>% 
  filter(!id %in% c('v0.28','v0.43'))

datalist = list()

for(k in 1:nrow(ctrl_file)){
  
  i = ctrl_file$id[k]  
  
  results_path <- file.path("results",i)
  
  diagnostic_ls <- readRDS(file.path(results_path,"diagnostics.rds"))
  
  dat <- data.frame(id = i,
    successful_chains = length(diagnostic_ls$num_divergent),
                    mean_divergences = mean(diagnostic_ls$num_divergent/(diagnostic_ls$num_iters-diagnostic_ls$num_warmups)),
                    low_ebfmi = ifelse(mean(diagnostic_ls$ebfmi)<0.3, TRUE, FALSE))
  
  datalist[[k]] <- dat
}

summarydat = do.call(rbind, datalist)

write_csv(summarydat, file=here("results","convergence_checks.csv"))
