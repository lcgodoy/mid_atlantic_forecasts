# load packages and data 
library(tidyverse)
library(here)

# packages for nice plots
library(tidytext)
library(ggrepel)
# data
dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))
dat_test <- read_csv(here("processed-data","flounder_catch_at_length_fall_testing.csv"))
dat_catchonly <- read_csv(here("processed-data","flounder_catch_fall_training.csv"))
dat_test_catchonly <- read_csv(here("processed-data","flounder_catch_fall_testing.csv"))
convergence_checks <- read_csv(file=here("results","convergence_checks.csv"))
dat_forecasts_summ <- read_csv(file = here("processed-data","model_comparison_summary.csv"))
ctrl_file <- read_csv(file=here("control_file.csv"))
dat_test_patch <- read_csv(file=here("processed-data","dat_test_patch.csv"))
points_for_plot <- read_csv(file=here("processed-data","points_for_plot.csv"))

drm_outputs_available_locally <- FALSE
############
# calculate all the statistics reported in-text
############

# how many hauls? 
length(unique(dat$haulid)) + length(unique(dat_test$haulid)) #12318
nrow(dat_catchonly) + nrow(dat_test_catchonly) #12318 — should be the same as above

# how many positive encounters? 
nrow(dat_catchonly[dat_catchonly$abundance>0,]) + nrow(dat_test_catchonly[dat_test_catchonly$abundance>0,]) # 2370 

# how frequently were more than 10 individuals caught? 
nrow(bind_rows(dat_catchonly, dat_test_catchonly) %>% 
       filter(abundance >= 10)) / nrow(bind_rows(dat_catchonly, dat_test_catchonly))

# how many individuals caught? 
sum(dat_catchonly$abundance) # 8145
sum(dat_test_catchonly$abundance) # 6368
# 14513 total 

# how many individuals caught per year?
bind_rows(dat_catchonly, dat_test_catchonly) %>% 
  group_by(year) %>% 
  summarise(n=sum(abundance))

# how many models passed convergence checks?
divergence_cutoff <- 0.05
chains_cutoff <- 3 
nrow(convergence_checks %>% 
       filter(mean_divergences <= divergence_cutoff, 
              successful_chains >= chains_cutoff))

############
# make figures
############
library(ggplot2)
# https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")

# area map figure

states <- ne_states(country = c("United States of America","Canada"), returnclass = "sf")
select <- dplyr::select

load(here("processed-data","stan_data_prep.Rdata"))
recdat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv")) %>% 
  mutate(lat_floor = floor(lat)) %>% 
  select(lat_floor, lat, lon) %>% 
  group_by(lat_floor) %>% 
  summarise(minlon = min(lon),
            maxlon = max(lon))

recs <- data_frame(ymin = patches, 
                   ymax = patches + 0.99, 
                   xmin = recdat$minlon + 0.1, 
                   xmax = xmin + 2,
                   x = xmin + (xmax-xmin)/2, 
                   y = ymin + 0.5
) # be sure to note that rectangle width is approximate

recs2 <- recs %>% 
  select(y, xmax, xmin)  %>%
  pivot_longer(cols = c(xmax, xmin), names_to="identity", values_to="x") %>% 
  group_by(y) %>% 
  mutate(meanx = mean(x),
         newcol = x - min(x)) # plus or minus? 

state_points <- cbind(st_centroid(states), st_coordinates(st_centroid(states$geometry))) %>% 
  select(X, Y, name, postal) %>% 
  mutate(
    X = ifelse(postal=='NJ', X+0.25, X)) %>% 
  filter(!postal == "DC") # nudge for plotting 

ggmap <- ggplot(data = states) +
  geom_rect(data=recs, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="black",  alpha=0.5) +
  geom_sf(fill = "grey92") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location = "bl", width_hint = 0.5, pad_x = unit(1, "in")) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(2.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_text(data= state_points,aes(x=X, y=Y, label=postal),
            color = "black", fontface = "bold", check_overlap = FALSE, position="dodge") +
  # geom_rect(data=recs, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="black",  alpha=0.5) +
  #  geom_tile(data = recs2, aes(x=meanx, y=y, width=2, height=1, fill=newcol), color="black") +
  scale_fill_gradient(low="darkblue", high="white") +
  coord_sf(xlim = c(-80, -65), ylim = c(34, 45.5), expand = FALSE) +
  annotate(geom = "text", x = -69.5, y = 38, label = "Atlantic Ocean", 
           fontface = "italic", color = "grey22", size = 6)  +
  NULL

ggsave(ggmap, filename=here("results","area_map.png"), width=110, height=110, dpi=600, units="mm")

# hovmoller plot 
dat_btemp <- bind_rows(dat_catchonly, dat_test_catchonly) %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(year, lat_floor) %>% 
  summarise(btemp = mean(btemp, na.rm=TRUE),
            n = length(haulid[!is.na(btemp)]))

gg_btemp <- ggplot() + 
  geom_raster(data=dat_btemp, aes(x=year, y=lat_floor, fill=btemp)) + 
  scale_fill_gradientn(colors=c("blue3","darkturquoise", "gold", "orangered", "red3"), limits = c(6,24), breaks = c(seq(6, 24, 2)),
                       guide = guide_colourbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE)) + 
  scale_x_continuous(limits = c(1971, 2017), breaks=seq(1962, 2016, 4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(34, 45), breaks=seq(35, 44, 1), expand = c(0, 0)) +
  labs(x="Year", y="Latitude") +
  theme_bw() +
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=12),
        axis.text=element_text(family="sans",size=8,color="black"), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title=element_text(family="sans",size=12,color="black"),
        plot.margin=margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")) + 
  guides(fill = guide_colourbar(barwidth = 0.75, barheight = 25, title="SBT", reverse=TRUE)) +
  NULL
ggsave(gg_btemp, filename=here("results","btemp_lat_time.png"), width=110, height=70, scale = 2.1, dpi=600, units="mm")

# model comparison plots 

# commented out this plot since it has too many points 
# gg_metrics <- dat_forecasts_summ  %>% 
#   mutate(feature = case_match(feature, "centroid" ~ "Centroid", "warm_edge" ~ "Warm Edge", "cold_edge" ~ "Cold Edge", .default=feature),
#          metric = case_match(metric, "rmse" ~ "RMSE", "bias" ~ "Bias", .default=metric), 
#          type = ifelse(str_detect(id, "0."),"DRM", id),
#          metric = factor(metric, levels = c("RMSE","Bias")))  %>%  
#   ggplot(aes(x=feature, y=value)) + 
#   geom_point(aes(color=type, fill=type),size=1) +
#   geom_text_repel(aes(label=id), hjust = 1, max.overlaps = Inf) +
#   scale_color_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
#   scale_fill_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
#   theme_bw() + 
#   facet_wrap(~metric, scales="free_y") +
#   labs(x="Feature", y="Metric") + 
#   theme(legend.position = "bottom")
# gg_metrics

# commented out this plot since it has too many points 
# gg_metrics2 <- dat_forecasts_summ  %>% 
#   mutate(feature = case_match(feature, "centroid" ~ "Centroid", "warm_edge" ~ "Warm Edge", "cold_edge" ~ "Cold Edge", .default=feature),
#          type = ifelse(str_detect(id, "0."),"DRM", id),
#          metric = factor(metric, levels = c("RMSE","Bias")),
#          id = gsub("v0.", "v", id))  %>%  
#   pivot_wider(names_from = metric, values_from = value) %>% 
#   ggplot(aes(x=RMSE, y=Bias)) + 
#   geom_point(aes(color=type, fill=type),size=1) +
#   geom_text_repel(aes(label=id), hjust = 1, max.overlaps = Inf) +
#   scale_color_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
#   scale_fill_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
#   theme_bw() + 
#   facet_wrap(~feature) +
#   theme(legend.position = "bottom",
#         legend.title=element_blank())
# gg_metrics2
# ggsave(gg_metrics2, filename=here("results","bias_v_rmse.png"), width=110, height=60, dpi=600, units="mm", scale=1.5)


# still too much to look at... need to trim down more 

# hacky way to get the best models 
best_drms <- dat_forecasts_summ %>% 
  filter(metric=='RMSE', !id %in% c('GAM','Persistence')) %>% 
  group_by(id) %>% 
  mutate(mean_value = mean(value)) %>% 
  filter(mean_value < 0.75)
length(unique(best_drms$id))

# gg_metrics3 <- dat_forecasts_summ  %>% 
#   filter(id %in% c(unique(best_drms$id), 'GAM','Persistence')) %>% 
#   mutate(feature = case_match(feature, "centroid" ~ "Centroid", "warm_edge" ~ "Warm Edge", "cold_edge" ~ "Cold Edge", .default=feature),
#          type = ifelse(str_detect(id, "0."),"DRM", id),
#          metric = factor(metric, levels = c("RMSE","Bias")),
#          id = gsub("v0.", "v", id))  %>%  
#   pivot_wider(names_from = metric, values_from = value) %>% 
#   ggplot(aes(x=RMSE, y=Bias)) + 
#   geom_point(aes(color=type, fill=type),size=1) +
#   geom_text_repel(aes(label=id), hjust = 1, max.overlaps = Inf) +
#   scale_color_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
#   scale_fill_manual(values= wesanderson::wes_palette("Darjeeling1", n = 3)) +
#   theme_bw() + 
#   facet_wrap(~feature) +
#   theme(legend.position = "bottom",
#         legend.title=element_blank())
# gg_metrics3
# ggsave(gg_metrics3, filename=here("results","bias_v_rmse.png"), width=110, height=60, dpi=600, units="mm", scale=1.5)

dat_evil_plot <- dat_forecasts_summ  %>% 
  filter(id %in% c(unique(best_drms$id), 'GAM','Persistence')) %>% 
  mutate(feature = case_match(feature, "centroid" ~ "Centroid", "warm_edge" ~ "Warm Edge", "cold_edge" ~ "Cold Edge", .default=feature),
         type = ifelse(str_detect(id, "0."),"DRM", id),
         id = gsub("v0.", "v", id))

gg_metrics4 <- dat %>%  
  mutate(
    metric = factor(metric, levels = c("RMSE","Bias")),
    id = reorder_within(id, value, list(metric, feature))) %>%  
  ggplot( )+
  geom_point(aes(id, y=value)) + 
  geom_hline(yintercept=0) + 
  coord_flip() + 
  scale_x_reordered() +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + 
  facet_grid(feature ~ metric, scales="free", axes="all_y", axis.labels = "all_y") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  NULL
gg_metrics4

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

if(drm_outputs_available_locally == TRUE) {
  for(k in unique(best_drms$id)){
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

# gg_bias <- dat_forecasts %>% 
#   ggplot(aes(x=year, y=resid,color=model_name, fill=model_name, shape = model_name )) + 
#   geom_line() +
#   geom_point() +
#   theme_bw() + 
#   theme(legend.position = "bottom") +
#   labs(x="Year", y="Residuals (° lat)") + 
#   facet_wrap(~feature)
# gg_bias