# https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library(tidyverse)
library(here)

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
