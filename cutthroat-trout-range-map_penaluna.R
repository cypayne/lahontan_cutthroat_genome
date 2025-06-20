## Cutthroat Trout subspecies distribution map, historical and contemporary ranges
## Map is adapted from Penaluna et al 2016
## Code creates distribution map figure (Figure 1B) that appears in 
## Payne & Escalona et al 2025. JHeredity.


# Read this shape file with the rgdal library.
library(sf)
library(ggplot2)
library(ggmap) # ggplot functionality for maps
library(data.table)
library(dplyr)
library(ggspatial) # for N arrow

## prepare historic trout ranges and site coordinates

# look at the available layers
#st_layers("fish_ranges.gpkg")
st_layers("penaluna_cutthroat_map_layers/greenback_h_uc_1.shp")
range_bonneville_h <- st_read("penaluna_cutthroat_map_layers/bonneville_h_1.shp")
range_bonneville_c <- st_read("penaluna_cutthroat_map_layers/bonneville_c_1.shp")
range_coastal_c <- st_read("penaluna_cutthroat_map_layers/coastal_c_h_1.shp")
range_colorado_h <- st_read("penaluna_cutthroat_map_layers/colorardo_h_1.shp")
range_colorado_c <- st_read("penaluna_cutthroat_map_layers/colorardo_c_1.shp")
range_greenback_h <- st_read("penaluna_cutthroat_map_layers/greenback_h_1.shp")
range_greenback_h_uc <- st_read("penaluna_cutthroat_map_layers/greenback_h_uc_1.shp")
range_greenback_c <- st_read("penaluna_cutthroat_map_layers/greenback_c_1.shp")
range_lahontan_h <- st_read("penaluna_cutthroat_map_layers/lahontan_h_1.shp")
range_lahontan_c <- st_read("penaluna_cutthroat_map_layers/lahontan_c_1.shp")
range_paiute_c <- st_read("penaluna_cutthroat_map_layers/paiute_c_h_1.shp")
range_riogrande_h <- st_read("penaluna_cutthroat_map_layers/riogrande_h_1.shp")
range_riogrande_h_uc1 <- st_read("penaluna_cutthroat_map_layers/riogrande_h_uc_1.shp")
range_riogrande_h_uc2 <- st_read("penaluna_cutthroat_map_layers/riogrande_h_uc_1.shp")
range_riogrande_c <- st_read("penaluna_cutthroat_map_layers/riogrande_c_1.shp")
range_westslope_h <- st_read("penaluna_cutthroat_map_layers/westslope_h_1.shp")
range_westslope_c1 <- st_read("penaluna_cutthroat_map_layers/westslope_c_1.shp")
range_westslope_c2 <- st_read("penaluna_cutthroat_map_layers/westslope_c_2.shp")
range_westslope_c3 <- st_read("penaluna_cutthroat_map_layers/westslope_c_3.shp")
range_yellowstone_h <- st_read("penaluna_cutthroat_map_layers/yellowstone_h_1.shp")
range_yellowstone_c <- st_read("penaluna_cutthroat_map_layers/yellowstone_c_1.shp")

# transform coordinates to EPSG:4326
range_bonneville_h <- range_bonneville_h %>% st_zm %>% st_transform(., 4326) %>% mutate(group="bonneville_h")
range_bonneville_c <- range_bonneville_c %>% st_zm %>% st_transform(., 4326) %>% mutate(group="bonneville_c") 
range_coastal_c <- range_coastal_c %>% st_zm %>% st_transform(., 4326) %>% mutate(group="coastal_c") 
range_colorado_h <- range_colorado_h %>% st_zm %>% st_transform(., 4326) %>% mutate(group="colorado_h") 
range_colorado_c <- range_colorado_c %>% st_zm %>% st_transform(., 4326) %>% mutate(group="colorado_c") 
range_greenback_h <- range_greenback_h %>% st_zm %>% st_transform(., 4326) %>% mutate(group="greenback_h") 
range_greenback_h_uc <- range_greenback_h_uc %>% st_zm %>% st_transform(., 4326) %>% mutate(group="greenback_uc") 
range_greenback_c <- range_greenback_c %>% st_zm %>% st_transform(., 4326) %>% mutate(group="greenback_c") 
range_lahontan_h <- range_lahontan_h %>% st_zm %>% st_transform(., 4326) %>% mutate(group="lahontan_h")
range_lahontan_c <- range_lahontan_c %>% st_zm %>% st_transform(., 4326) %>% mutate(group="lahontan_c")
range_paiute_c <- range_paiute_c %>% st_zm %>% st_transform(., 4326) %>% mutate(group="paiute_c")
range_riogrande_h <- range_riogrande_h %>% st_zm %>% st_transform(., 4326) %>% mutate(group="riogrande_h") 
range_riogrande_c <- range_riogrande_c %>% st_zm %>% st_transform(., 4326) %>% mutate(group="riogrande_c") 
range_riogrande_h_uc1 <- range_riogrande_h_uc1 %>% st_zm %>% st_transform(., 4326) %>% mutate(group="riogrande_uc")
range_riogrande_h_uc2 <- range_riogrande_h_uc2 %>% st_zm %>% st_transform(., 4326) %>% mutate(group="riogrande_uc")
range_westslope_h <- range_westslope_h %>% st_zm %>% st_transform(., 4326) %>% mutate(group="westslope_h")
range_westslope_c1 <- range_westslope_c1 %>% st_zm %>% st_transform(., 4326) %>% mutate(group="westslope_c")
range_westslope_c2 <- range_westslope_c2 %>% st_zm %>% st_transform(., 4326) %>% mutate(group="westslope_c")
range_westslope_c3 <- range_westslope_c3 %>% st_zm %>% st_transform(., 4326) %>% mutate(group="westslope_c")
range_yellowstone_h <- range_yellowstone_h %>% st_zm %>% st_transform(., 4326) %>% mutate(group="yellowstone_h")
range_yellowstone_c <- range_yellowstone_c %>% st_zm %>% st_transform(., 4326) %>% mutate(group="yellowstone_c")

# cutthroat palette
#cutthroat_palette <- c("bonneville_h"="#116178","bonneville_c"="#003c4d",
#                       "colorado_h"="#dc6601","colorado_c"="#d24e01",
#                       "coastal_c"="#BFBDD1","paiute_c"="#cc877d",
#                       "lahontan_h"="#4c5e46","lahontan_c"="#334a36",
#                       "greenback_h"="#5c6b28","greenback_c"="#424a26","greenback_uc"="grey",
#                       "riogrande_h"="#7b6b45","riogrande_c"="#4e452e","riogrande_uc"="grey",
#                       "westslope_h"="#a03b27","westslope_c"="#5e271b",
#                       "yellowstone_h"="#d7b013","yellowstone_c"="#a36f04")

cutthroat_palette <- c("bonneville_h"="#116178","bonneville_c"="#003c4d",
                       "colorado_h"="#dc6601","colorado_c"="#d24e01",
                       "coastal_c"="#BFBDD1","paiute_c"="black",
                       "lahontan_h"="#f76","lahontan_c"="#ff5742",
                       "greenback_h"="#5c6b28","greenback_c"="#424a26","greenback_uc"="grey",
                       "riogrande_h"="#a03b27","riogrande_c"="#5e271b","riogrande_uc"="grey",
                       "westslope_h"="#4c5e46","westslope_c"="#334a36",
                       "yellowstone_h"="#d7b013","yellowstone_c"="#a36f04")


# Pyramid Lake coordinates
map_points <- data.frame(location = c("Pyramid Lake"),
                     x = c(-119.52301025),
                     y = c(40.03798276))

# add major water bodies (only US though)
#water_bodies <- st_read("USA_Detailed_Water_Bodies/USA_Detailed_Water_Bodies.shp")
#lakes <- st_read("North_America_Lakes/North_America_Lakes_and_Rivers.shp") # seems to only be lakes
#lakes <- lakes %>% st_zm %>% st_transform(., 4326) %>% mutate(group="water") %>% filter(Type=="Lake / Lago / Lac")
#rivers <- st_read("North_America_Rivers/North_America_Lakes_and_Rivers.shp") # seems to only be rivers
#rivers <- rivers %>% st_zm %>% st_transform(., 4326) %>% mutate(group="water")


## cutthroat subspecies map
us_states <- map_data("state")
alaska <- map_data("world","USA:Alaska")
canada <- map_data("world","Canada")
mexico <- map_data("world","Mexico")
cutt_base <- ggplot() +
  geom_polygon(data = us_states, aes(x = long, y = lat, group = group), color = "white", fill = "#ECF1ED") +
  geom_polygon(data = canada, aes(x = long, y = lat, group = group), color = "white", fill = "#ECF1ED") +
  geom_polygon(data = alaska, aes(x = long, y = lat, group = group), color = "white", fill = "#ECF1ED") +
  geom_polygon(data = mexico, aes(x = long, y = lat, group = group), color = "white", fill = "#ECF1ED") +
  #geom_sf(data = lakes, aes(fill=group), lwd=0.2, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  #geom_sf(data = rivers, aes(fill=group), lwd=0.2, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  theme_bw() +
  theme(legend.position="none") +
  coord_sf(xlim=c(-150,-100), ylim=c(62,28)) +
  xlab("longitude") +
  ylab("latitude") 
cutt_base
  
cutt_map <- cutt_base +
  geom_sf(data = range_bonneville_h, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_bonneville_c, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_coastal_c, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_colorado_h, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_colorado_c, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_greenback_h, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_greenback_c, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_greenback_h_uc, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.4,show.legend="polygon") +
  geom_sf(data = range_lahontan_h, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_lahontan_c, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_paiute_c, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_riogrande_h, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_riogrande_c, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_riogrande_h_uc1, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.4,show.legend="polygon") +
  geom_sf(data = range_riogrande_h_uc2, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.4,show.legend="polygon") +
  geom_sf(data = range_westslope_h, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_westslope_c1, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_westslope_c2, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_westslope_c3, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_yellowstone_h, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  geom_sf(data = range_yellowstone_c, aes(fill=group), color=NA, lwd=0, inherit.aes = FALSE, alpha=0.7,show.legend="polygon") +
  scale_fill_manual(values=cutthroat_palette) + 
  geom_point(data = map_points, aes(x = x, y = y), size=0.5, shape=16, color="white") +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  #theme(legend.position = c(.85, .75), text = element_text(size = 12)) +
  theme(legend.position="none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  coord_sf(xlim=c(-150,-100), ylim=c(62,28)) +
  xlab("Longitude") +
  ylab("Latitude") 
cutt_map
ggsave(cutt_map, filename='historical-cutthroat-trout-range-map_penaluna-adapted.pdf', width=8, height=8, bg="transparent")



