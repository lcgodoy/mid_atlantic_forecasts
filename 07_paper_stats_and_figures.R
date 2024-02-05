# load packages and data 
library(tidyverse)
library(here)
dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))
dat_test <- read_csv(here("processed-data","flounder_catch_at_length_fall_testing.csv"))

dat_catchonly <- read_csv(here("processed-data","flounder_catch_fall_training.csv"))
dat_test_catchonly <- read_csv(here("processed-data","flounder_catch_fall_testing.csv"))

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
