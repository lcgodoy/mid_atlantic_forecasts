###################################
# READ ME:

# this script uses the 2020 release of OceanAdapt
# to re-run it, you will need to download it
# it can be downloaded from: https://zenodo.org/record/3885625
# DOI 10.5281/zenodo.3885625
# it should download a in folder called pinskylab-OceanAdapt-c913be5
# you can place that directly in this repository, mid_atlantic_forecasts, for the relative file paths below to work
# or to avoid cluttering up this repo, create a new folder inside it called pinskylab-OceanAdapt-c913be5, and move in only the needed files
# the four files you will need are:
# data_clean/dat_exploded.rds (this is a big file which is why we can't just put these all on github)
# data_raw/neus_Survdat.RData
# data_raw/neus_SVSPP.RData
# data_raw/neus_fall_svsta.csv

###################################

# load packages
library(tidyverse)
library(here)
library(stringr)
library(lubridate)
here <- here::here

OApath <- "pinskylab-OceanAdapt-c913be5"

dat_exploded <- readRDS(here(OApath,"data_clean","dat_exploded.rds")) # get zero-inflated survey data

# get haul-level variables
load(here(OApath,"data_raw","neus_Survdat.RData")) 

# note that in the NEFSC survey you may need to aggregate across sex classes ('CATCHSEX' column) for correct abundance/length values. for summer flounder that column is always 0 so it's not an issue

# get taxonomy
load(here(OApath, "data_raw","neus_SVSPP.RData"))

# get tow duration data 
towdur_raw <- read_csv(here(OApath, "data_raw","neus_fall_svsta.csv"))

spp_of_interest <- c("Paralichthys dentatus")
reg_of_interest <- c("Northeast US Fall")

dat_exploded_neus <- dat_exploded %>% 
  filter(region==reg_of_interest) 

max_yr <- 2016 # 2017 is missing, for some reason 
min_yr <- 1972 # early years have some issues; the newly filtered OA data starts here anyway
# to explore more the changes in samplng over time, make annual maps of haul locations, or a tile plot of year vs. lat band 
forecast_yr_1 <- max_yr - 9

# trim to only hauls that we know were close to standard length (30min), which we use for the CPUE conversion 
# sensu https://github.com/afredston/marine_heatwaves_trawl/blob/main/prep_trawl_data.R#L55
towdur <- towdur_raw %>% 
  filter(!STATION %in% c('8093','7046')) %>% 
  mutate(
  CRUISE6 = as.character(CRUISE6), 
  STATION = str_sub(STATION, start=2, end=4),
  STRATUM = str_sub(STRATUM, start=2, end=5), # QA/QC up to this point is just to make the haulid consistent with the other NEUS dataframe haulids 
  haulid = paste(formatC(CRUISE6, width=6, flag=0), formatC(STATION, width=3, flag=0), formatC(STRATUM, width=4, flag=0), sep='-')) %>% 
  select(haulid, TOWDUR, EST_YEAR) %>% 
  mutate(haul_dur = TOWDUR / 60) %>% 
  rename(
    year = EST_YEAR) %>% 
  mutate(flag = case_when(
    year < 2009 & haul_dur < 0.42 ~ "drop",
    year < 2009 & haul_dur > 0.58 ~ "drop",
    year >= 2009 & haul_dur < 0.25 ~ "drop",
    year >= 2009 & haul_dur > 0.42 ~ "drop",
    haul_dur == NA ~ "drop",
    .default = "keep"
  )) 

towdur %>% group_by(flag) %>% summarise(n=n()) # drops 174 of 19102 hauls, ie <1%  

towdur <- towdur %>% 
  filter(flag == 'keep') 

dat_exploded_neus <- dat_exploded_neus %>% 
  filter(haulid %in% unique(towdur$haulid)) # drop hauls that are an irregular duration 

# creating a separate haul info dataframe to get the date and btemp
hauldat <- survdat %>% 
  # create a haulid for joining
  mutate(haulid = paste(formatC(CRUISE6, width=6, flag=0), formatC(STATION, width=3, flag=0), formatC(STRATUM, width=4, flag=0), sep='-')) %>% 
  select(haulid, BOTTEMP, YEAR, EST_TOWDATE, LAT, LON) %>% 
  distinct() %>%
  rename("btemp"=BOTTEMP,
         "year"=YEAR,
         "date"=EST_TOWDATE,
         "lat"=LAT,
         "lon"=LON)
        
# tidy length data
len_flounder_prep <- survdat %>% 
  # create a haulid for joining with dat.exploded
  mutate(haulid = paste(formatC(CRUISE6, width=6, flag=0), formatC(STATION, width=3, flag=0), formatC(STRATUM, width=4, flag=0), sep='-')) %>% 
  select(haulid, SVSPP, LENGTH, NUMLEN) %>% 
  left_join(spp, by="SVSPP") %>% # get species names from species codes
  mutate(spp = str_to_sentence(SCINAME)) %>%
  select(spp, haulid, LENGTH, NUMLEN) %>%
  filter(!is.na(LENGTH),
         spp == spp_of_interest,
         LENGTH>4 # get rid of a single fish at 4cm -- the next smallest is 9cm
  ) %>% 
  rename("length"=LENGTH,
         "number_at_length"=NUMLEN)

# need to create length df that still includes zeroes
# important to use haulids from dat_exploded_neus because it has been cleaned
len_flounder <- expand.grid(haulid=unique(dat_exploded_neus$haulid), length=seq(min(len_flounder_prep$length), max(len_flounder_prep$length), 1)) %>% # get full factorial of every haul * length bin
  mutate(spp = spp_of_interest) %>% 
  # left_joining to use only the hauls in dat_exploded_neus
  left_join(len_flounder_prep, by=c('length','haulid','spp')) %>% 
  mutate(number_at_length = replace_na(number_at_length, 0)) %>%  # fill in absences with true zeroes
  left_join(hauldat)

# count up all abundances for SDMs, which don't use length data 
abund_flounder <- len_flounder %>% 
  group_by_at(vars(c(-number_at_length, -length))) %>% 
  summarise(abundance = sum(number_at_length)) %>%
  left_join(dat_exploded_neus %>% select(haulid, depth) %>% distinct(), by="haulid") # get depth in case wanted for SDM

# use dat_exploded_neus for biomass SDMs
bio_flounder <- dat_exploded_neus %>% 
  filter(spp == spp_of_interest) %>% 
  select(-common, -stratum, -stratumarea) %>%
  left_join(hauldat %>% select(haulid, btemp), by="haulid")

# split all into training/testing and save 
abund_len_flounder_train <- len_flounder %>% 
  filter(year >= min_yr,
         year < forecast_yr_1)

abund_flounder_train <- abund_flounder %>% 
  filter(year >= min_yr,
         year < forecast_yr_1)

bio_flounder_train <- bio_flounder %>% 
  filter(year >= min_yr,
         year < forecast_yr_1)

abund_len_flounder_test <- len_flounder %>% 
  filter(year >= forecast_yr_1,
         year <= max_yr)

abund_flounder_test <- abund_flounder %>% 
  filter(year >= forecast_yr_1,
         year <= max_yr)

bio_flounder_test <- bio_flounder %>% 
  filter(year >= forecast_yr_1,
         year <= max_yr)

write_csv(abund_len_flounder_train, here("processed-data","flounder_catch_at_length_fall_training.csv"))
write_csv(abund_len_flounder_test, here("processed-data","flounder_catch_at_length_fall_testing.csv"))

write_csv(abund_flounder_train, here("processed-data","flounder_catch_fall_training.csv"))
write_csv(abund_flounder_test, here("processed-data","flounder_catch_fall_testing.csv"))


write_csv(bio_flounder_train, here("processed-data","flounder_biomass_fall_training.csv"))
write_csv(bio_flounder_test, here("processed-data","flounder_biomass_fall_testing.csv"))
