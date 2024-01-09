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

ctrl_file <- read_csv("control_file.csv") 

for(k in 1:nrow(ctrl_file)){
  
  i = ctrl_file$id[k]  
  
  results_path <- file.path("results",i)
  
  diagnostic_ls <- readRDS(file.path(results_path,"diagnostics.rds"))
  
  
}

