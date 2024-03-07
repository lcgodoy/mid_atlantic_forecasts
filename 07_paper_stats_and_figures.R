# load packages and data 

library(tidyverse)
library(here)

# packages for map
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")

# packages for nice plots
library(ggplot2)
library(tidytext)
library(ggrepel)
library(ggh4x)
library(ggdist)
theme_set(theme_bw())

# data
dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))
dat_test <- read_csv(here("processed-data","flounder_catch_at_length_fall_testing.csv"))
dat_catchonly <- read_csv(here("processed-data","flounder_catch_fall_training.csv"))
dat_test_catchonly <- read_csv(here("processed-data","flounder_catch_fall_testing.csv"))
convergence_checks <- read_csv(file=here("results","convergence_checks.csv"))
dat_forecasts_summ <- read_csv(file = here("processed-data","model_comparison_summary.csv"))
ctrl_file <- read_csv(file=here("control_file.csv"))
dat_test_patch <- read_csv(file=here("processed-data","dat_test_patch.csv"))
gam_out <- read_csv(file = here("processed-data","gam_abundance_time.csv"))
points_for_plot <- read_csv(file=here("processed-data","points_for_plot.csv"))
load(here("processed-data","stan_data_prep.Rdata"))

drm_outputs_available_locally <- TRUE

if(drm_outputs_available_locally){
  drm_out <- read_csv(here("processed-data","posteriors_for_model_evaluation.csv")) }
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

# did summer flounder shift north? 
lm_dat <- points_for_plot %>% 
  filter(id == 'Observed')
summary(lm(value_tmp ~ year, data = lm_dat %>% 
             filter(feature == 'Centroid')))
summary(lm(value_tmp ~ year, data = lm_dat %>% 
             filter(feature == 'Warm Edge')))
summary(lm(value_tmp ~ year, data = lm_dat %>% 
             filter(feature == 'Cold Edge')))

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

gg_mod_compare <- dat_evil_plot %>%  
  mutate(
    metric = factor(metric, levels = c("RMSE","Bias")),
    id = reorder_within(id, value, list(metric, feature))) %>%  
  ggplot( )+
  geom_point(aes(id, y=value, shape = type)) + 
  geom_hline(yintercept=0, linetype="dashed") + 
  coord_flip() + 
  scale_x_reordered() +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + 
  facet_grid(feature ~ metric, scales="free", axes="all_y", axis.labels = "all_y") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none") +
  ggh4x::facet_nested_wrap(~metric+feature, scales = "free") 
ggsave(gg_mod_compare, filename=here("results","bias_v_rmse.png"), width=110, height=80, dpi=600, units="mm", scale=1.5)

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
    
    tmp_dens <- read_rds(file.path(results_path, "density_obs_proj.rds")) %>% 
      mutate(year = year + min(years_proj) - 1,
             patch = patch + min(patches) - 1) %>% 
      filter(year < 2017)
    
    tmp_edges <- read_rds(file.path(results_path, "range_quantiles_proj.rds")) %>% 
      mutate(year = year + min(years_proj) - 1,
             range_quantiles_proj = range_quantiles_proj + min(patches) - 1) %>% 
      filter(range_quantiles_proj < Inf,
             year < 2017)
    
    tmp_centroids <- read_rds(file.path(results_path, "centroid_proj.rds"))%>% 
      mutate(year = year + min(years_proj) - 1) %>% 
      filter(year < 2017)
    
    # gg_range_time_drm <- ggplot() +
    #   stat_lineribbon(data = tmp_edges[tmp_edges$quantile==0.05,], aes(x=year, y=range_quantiles_proj))+ 
    #   stat_lineribbon(data = tmp_edges[tmp_edges$quantile==0.95,], aes(x=year, y=range_quantiles_proj))+ 
    #   stat_lineribbon(data = tmp_centroids, aes(x=year, y=centroid_proj))+ 
    #   geom_point(data=points_for_plot %>% filter(id == "Observed") %>% rename("Feature" = feature), aes(x=year, y=value_tmp, shape=Feature), color="#E20134", size=2) + 
    #   geom_line(data=points_for_plot %>% filter(id == "Observed")%>% rename("Feature" = feature), aes(x=year, y=value_tmp, group=Feature), color="#E20134") +
    #   scale_fill_brewer(name = "Credible Interval") +
    #   scale_x_continuous(breaks =seq(2007, 2016, 1)) + 
    #   scale_y_continuous(breaks = seq(35, 44, 1), labels =  seq(35, 44, 1), limits=c(34, 44.5)) +
    #   labs(x="Year", y="Latitude", title = paste0("DRM ",k)) + 
    #   theme_bw()
    # ggsave(gg_range_time_drm, filename=paste0(here("results"),"/range_time_",k,".png"), dpi=600, units="mm", width=130, height=75)
    
    gg_range_time_drm_cent <- ggplot() +
      stat_lineribbon(data = tmp_centroids, aes(x=year, y=centroid_proj))+ 
      geom_point(data=points_for_plot %>% filter(id == "Observed", feature == 'Centroid') %>% rename("Feature" = feature), aes(x=year, y=value_tmp, shape=Feature), color="#E20134", size=2) + 
      geom_line(data=points_for_plot %>% filter(id == "Observed", feature == 'Centroid')%>% rename("Feature" = feature), aes(x=year, y=value_tmp, group=Feature), color="#E20134") +
      scale_fill_brewer(name = "Credible Interval") +
      scale_x_continuous(breaks =seq(2007, 2016, 1)) + 
      scale_y_continuous(breaks = seq(37, 40, 1), labels =  seq(37, 40, 1), limits=c(36, 41)) +
      labs(x="Year", y="Latitude", title = paste0("Centroid, DRM ",gsub("0.", "", k))) + 
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, vjust=0.8))
    ggsave(gg_range_time_drm_cent, filename=paste0(here("results"),"/centroid_time_",k,".png"), dpi=600, units="mm", width=75, height=55)
    
    gg_range_time_drm_cold <- ggplot() +
      stat_lineribbon(data = tmp_edges[tmp_edges$quantile==0.95,], aes(x=year, y=range_quantiles_proj))+ 
      geom_point(data=points_for_plot %>% filter(id == "Observed", feature == 'Cold Edge') %>% rename("Feature" = feature), aes(x=year, y=value_tmp, shape=Feature), color="#E20134", size=2) + 
      geom_line(data=points_for_plot %>% filter(id == "Observed", feature == 'Cold Edge')%>% rename("Feature" = feature), aes(x=year, y=value_tmp, group=Feature), color="#E20134") +
      scale_fill_brewer(name = "Credible Interval") +
      scale_x_continuous(breaks =seq(2007, 2016, 1)) + 
      scale_y_continuous(breaks = seq(38, 44, 1), labels =  seq(38, 44, 1), limits=c(38, 44.1)) +
      labs(x="Year", y="Latitude", title = paste0("Cold edge, DRM ",gsub("0.", "", k))) + 
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, vjust=0.8))
    ggsave(gg_range_time_drm_cold, filename=paste0(here("results"),"/cold_edge_time_",k,".png"), dpi=600, units="mm", width=75, height=55)
    
    gg_range_time_drm_warm <- ggplot() +
      stat_lineribbon(data = tmp_edges[tmp_edges$quantile==0.05,], aes(x=year, y=range_quantiles_proj))+ 
      geom_point(data=points_for_plot %>% filter(id == "Observed", feature == 'Warm Edge') %>% rename("Feature" = feature), aes(x=year, y=value_tmp, shape=Feature), color="#E20134", size=2) + 
      geom_line(data=points_for_plot %>% filter(id == "Observed", feature == 'Warm Edge')%>% rename("Feature" = feature), aes(x=year, y=value_tmp, group=Feature), color="#E20134") +
      scale_fill_brewer(name = "Credible Interval") +
      scale_x_continuous(breaks =seq(2007, 2016, 1)) + 
      scale_y_continuous(breaks = seq(35, 38, 1), labels =  seq(35, 38, 1), limits=c(34, 39)) +
      labs(x="Year", y="Latitude", title = paste0("Warm edge, DRM ",gsub("0.", "", k))) + 
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, vjust=0.8))
    ggsave(gg_range_time_drm_warm, filename=paste0(here("results"),"/warm_edge_time_",k,".png"), dpi=600, units="mm", width=75, height=55)
    
    gg_est_tile_drm <- tmp_dens %>%
      group_by(patch, year) %>% 
      summarise(Abundance = mean(density_obs_proj)) %>% 
      rename("Latitude" = patch, "Year" = year) %>% 
      ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
      geom_tile() +
      scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
      scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
      #     scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
      labs(x="Year", y="Latitude", title = paste0("DRM ",gsub("0.", "", k))) + 
      theme(legend.position="none",
            axis.text.x = element_text(angle = 45, vjust=0.8)) +
      NULL
    ggsave(gg_est_tile_drm, filename=paste0(here("results"),"/tileplot_time_",k,".png"), dpi=600, units="mm", width=75, height=75)
    
  } 
}

#     gg_length <- n_at_length_hat %>% 
#       ggplot(aes(length, plength)) +
#       stat_lineribbon(size = .1) +
#       geom_point(data = n_p_l_y, aes(length, plength), color = "red", alpha = 0.25) +
#       facet_grid(patch ~ year, scales = "free_y") + 
#       scale_x_continuous(limits = c(0, 50)) # generates warnings because of the close-to-zero probabilities at larger lengths

gg_range_time_warm <- ggplot(data = points_for_plot %>% filter(feature=='Warm Edge') %>% 
                             rename("Latitude" = value_tmp, "Year" = year) %>% 
                               mutate(id = factor(id, levels=c('Observed','GAM','Persistence')))) +
  geom_point(aes(x=Year, y=Latitude, color=id, fill=id, shape = id), size=2) + 
  geom_line(aes(x=Year, y=Latitude, color=id)) + 
  scale_x_continuous(breaks =seq(2007, 2016, 1)) + 
  scale_y_continuous(breaks = seq(35, 38, 1), labels =  seq(35, 38, 1), limits=c(34, 39)) +
  scale_color_manual(values=c("#E20134","#8400CD","#009F81")) + 
  scale_fill_manual(values=c("#E20134","#8400CD","#009F81")) + 
  labs(x="Year", y="Latitude", title = "Warm edge") + 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust=0.8))
ggsave(gg_range_time_warm, filename=paste0(here("results"),"/warm_edge_time.png"), dpi=600, units="mm", width=75, height=55)

gg_range_time_cold <- ggplot(data = points_for_plot %>% filter(feature=='Cold Edge') %>% 
                               rename("Latitude" = value_tmp, "Year" = year) %>% 
                               mutate(id = factor(id, levels=c('Observed','GAM','Persistence')))) +
  geom_point(aes(x=Year, y=Latitude, color=id, fill=id, shape = id), size=2) + 
  geom_line(aes(x=Year, y=Latitude, color=id)) + 
  scale_x_continuous(breaks =seq(2007, 2016, 1)) + 
  scale_y_continuous(breaks = seq(38, 44, 1), labels =  seq(38, 44, 1), limits=c(38, 44.1)) +
  scale_color_manual(values=c("#E20134","#8400CD","#009F81")) + 
  scale_fill_manual(values=c("#E20134","#8400CD","#009F81")) + 
  labs(x="Year", y="Latitude", title = "Cold edge") + 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust=0.8))
ggsave(gg_range_time_cold, filename=paste0(here("results"),"/cold_edge_time.png"), dpi=600, units="mm", width=75, height=55)

gg_range_time_centroid <- ggplot(data = points_for_plot %>% filter(feature=='Centroid') %>% 
                               rename("Latitude" = value_tmp, "Year" = year) %>% 
                               mutate(id = factor(id, levels=c('Observed','GAM','Persistence')))) +
  geom_point(aes(x=Year, y=Latitude, color=id, fill=id, shape = id), size=2) + 
  geom_line(aes(x=Year, y=Latitude, color=id)) + 
  scale_x_continuous(breaks =seq(2007, 2016, 1)) + 
  scale_y_continuous(breaks = seq(37, 40, 1), labels =  seq(37, 40, 1), limits=c(36, 41)) +
  scale_color_manual(values=c("#E20134","#8400CD","#009F81")) + 
  scale_fill_manual(values=c("#E20134","#8400CD","#009F81")) + 
  labs(x="Year", y="Latitude", title = "Centroid") + 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust=0.8))
ggsave(gg_range_time_centroid, filename=paste0(here("results"),"/centroid_time.png"), dpi=600, units="mm", width=75, height=55)

abund_p_y_proj <-  dat_test_dens %>%
  mutate(abundance = mean_dens * (1/0.0384) * meanpatcharea)

gg_observed_abundance_tile <- abund_p_y_proj %>%
  mutate(Year = (year + min(years_proj) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>%
  ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
  geom_tile() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  #  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Relative abundance (data)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust=0.8))
ggsave(gg_observed_abundance_tile, filename=paste0(here("results"),"/tileplot_time_observed.png"), dpi=600, units="mm", width=75, height=75)

gg_persistence_tile <- expand_grid(patches, years_proj) %>% 
  rename("Year" = years_proj, "Latitude" = patches) %>% 
  left_join(dat_train_dens %>% 
              filter(year == max(year)) %>% 
              select(-year, -patch) %>% 
              rename("Latitude" = lat_floor, "Abundance"=mean_dens)) %>% 
            ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
              geom_tile() +
              scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
              scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
              #  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
              labs(title="Persistence") +
              theme(legend.position = "none",
                    axis.text.x = element_text(angle = 45, vjust=0.8))
ggsave(gg_persistence_tile, filename=paste0(here("results"),"/tileplot_time_persistence.png"), dpi=600, units="mm", width=75, height=75)

gg_gam_tile <- gam_out %>% 
  rename("Latitude" = lat_floor, "Year" = year, "Abundance" = dens_pred) %>% 
  ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
  geom_tile() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  #  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="GAM") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust=0.8))
ggsave(gg_gam_tile, filename=paste0(here("results"),"/tileplot_time_gam.png"), dpi=600, units="mm", width=75, height=75)
