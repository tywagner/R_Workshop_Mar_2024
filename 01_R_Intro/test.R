library(dplyr) # data management
library(tidyverse)
library(ggplot2) # plot
library(lubridate) # dates
library(bslib) # a modern UI toolkit for Shiny and R Markdown
library(stringr) # manipulate character stings
library(sf) # map creation
library(spData) # spatial data
library(car)  # logit function
library(ggmap)

library(nhdplusTools)
# library(rgdal)
library(ggplot2)
library(googleway)
library(osmdata) # get waterbody data for mapping

fhc_dat <- read_csv('./02_data/FHC_length_weigtht.csv')
str(fhc_dat)

fhc_dat1 <-  fhc_dat %>%
  rename("fish_id" = `Fish ID`,
         "length_mm" = `Length (mm)`,
         "weight_g" = `Weight (kg)`,
         "fin_clip" = `Fin Clip ( Y or N )`,
         "liver_muscle" = `Liver/Muscle (Y or N)`)
head(fhc_dat1, 3)


fhc_dat2 <-  fhc_dat %>%
  rename("fish_id" = `Fish ID`,
         "length_mm" = `Length (mm)`,
         "weight_g" = `Weight (kg)`,
         "fin_clip" = `Fin Clip ( Y or N )`,
         "liver_muscle" = `Liver/Muscle (Y or N)`) %>%
  mutate(Date = mdy(Date), # mdy is from the lubridate package
         fish_id = factor(fish_id),
         Site = factor(Site),
         log_length = log(length_mm),
         log_weight = log(weight_g))
head(fhc_dat2, 3)

site_dat <- read_csv('./02_data/FHC_site_info.csv')
str(site_dat) # display structure of R object


fhc_dat_merged <- fhc_dat2 %>%
  rename("site_name" = Site) %>% # rename so they match in each data set
  select(-fin_clip, -liver_muscle) %>% # remove unwanted columns
  mutate(site_name = recode_factor(site_name, "East Donegal"="East Donegal Township")) %>% # change factor level name
  left_join(site_dat, by=c("site_name")) # merge
head(fhc_dat_merged)
dim(fhc_dat_merged)

# Grab state boundaries from spData and
# transform the coordinate references system (crs)
# crs = 4326 = WGS84; WGS84 CRS is often used for lat and long positions
us_states2 <- st_transform(us_states, crs = 4326)
# Rename column
colnames(us_states2)[2] <- "State"
# Select state(s) of interest
selectStates <- c("Pennsylvania")
# Subset data for plotting
us_state_select <- us_states2[us_states2$State %in% selectStates, ]

map.dat <- st_as_sf(site_dat, coords = c("long", "lat"), crs = 4326)

ggplot() +
  geom_sf(data = us_state_select, color = "gray30", lwd=1, fill="grey80") +
  geom_sf(data=filter(map.dat, site_name=='Bainbridge' | site_name=='Union'), shape=16,    size = 3, colour="red") +
  geom_sf(data=filter(map.dat, site_name=='Clemson Island' | site_name=='East Donegal Township' |site_name=='Shady Nook'), shape='x',    size = 5, colour="green") +
  labs(title="", y="Latitude", x="Longitude") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))

bb1 <- st_bbox(map.dat)
# plot_nhdplus(bbox = bb1, streamorder = 3,
#              stoponlargerequest = FALSE)

# NHD Flowlines


# USGS 01578310 SUSQUEHANNA RIVER AT CONOWINGO, MD
# plot_nhdplus("01578310")

# first identify the gage of interest
nldi_nwis <- list(featureSource = "nwissite", featureID = "USGS-01578310")

# next download the basin:
basin <- get_nldi_basin(nldi_feature = nldi_nwis)

# find out comid of the gage or point:
(comid_pt <- discover_nhdplus_id(nldi_feature = nldi_nwis))

huc_test <- get_huc(basin, type = "huc12")
str(huc_test)
# ----- Get Streamline Data
# We also can specify whether we want upstream or downstream mainstem or tributaries.
# UT= “Upstream Tributaries”, UM=“Upstream Mainstem”, DM=“Downstream Main”,
# or DD=“Downstream Diversions”. We can also specify a distance we want to
# travel from that initial starting point. The default is 10km.


# pull mainstems from a USGS gage
us_main <- navigate_nldi(nldi_feature = nldi_nwis,
                         mode="UM",
                         distance_km =  800)$UM %>%
  st_geometry()

# -- (not work) Subset UT comids based on stream order
# us_main2 <- navigate_nldi(nldi_feature = nldi_nwis,
#                          mode="UM",
#                          distance_km =  800)
#
# sub_main <- subset_nhdplus(us_main2$UM_flowlines$nhdplus_comid, streamorder = 3)
#


# Pull upstream tributaries
us_tribs <- navigate_nldi(nldi_feature = nldi_nwis,
                          mode="UT", distance_km =  800)$UT %>%
  st_geometry()

# ----------- Grab West Branch
# 01553500
nldi_nwis_WB <- list(featureSource = "nwissite", featureID = "USGS-01553500")
# find out comid of the gage or point:
(comid_pt_WB <- discover_nhdplus_id(nldi_feature = nldi_nwis_WB))
# Upstream
us_main_WB <- navigate_nldi(nldi_feature = nldi_nwis_WB,
                         mode="UM",
                         distance_km =  1500)$UM %>%
  st_geometry()
# Downstream
ds_main_WB <- navigate_nldi(nldi_feature = nldi_nwis_WB,
                            mode="DM",
                            distance_km =  80)$DM %>%
  st_geometry()

ggplot() +
  geom_sf(data = us_state_select, color = "gray30", lwd=1, fill="grey80") +
  geom_sf(data=basin) +
  geom_sf(data=us_tribs, col="lightblue", alpha=0.2) +
  geom_sf(data=us_main, col = 'black') +
  geom_sf(data=us_main_WB, col = 'black') +
  geom_sf(data=ds_main_WB, col = 'black') +
  geom_sf(data=filter(map.dat, site_name=='Bainbridge' | site_name=='Union'),
          shape=16,    size = 3, colour="red") +
  geom_sf(data=filter(map.dat, site_name=='Clemson Island' | site_name=='East Donegal Township' |site_name=='Shady Nook'),
          shape='x',    size = 5, colour="green") +
  labs(title="", y="Latitude", x="Longitude") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))

## ------------ Final map
selectStates2 <- c("Pennsylvania", "New York")
# Subset data for plotting
us_state_select2 <- us_states2[us_states2$State %in% selectStates2, ]

cities <- data.frame(city = c("Pittsburgh", "Harrisburg", "Philadephia"),
                     lat = c(40.440624, 40.263680, 39.952583),
                     long = c(-79.995888, -76.890739, -75.165222),
                     plat = c(40.440624, 40.263680, 39.952583),
                     plong = c(-79.995888, -76.890739, -75.165222))
citites.dat <- st_as_sf(cities, coords = c("long", "lat"), crs = 4326)

p1 <- ggplot() +
  geom_sf(data = us_state_select2, color = "gray30", lwd=1, fill="grey80") +
  geom_sf(data=basin) +
  geom_sf(data=us_tribs, col="lightblue", alpha=0.2) +
  geom_sf(data=us_main, col = 'darkblue') +
  geom_sf(data=us_main_WB, col = 'darkblue') +
  geom_sf(data=ds_main_WB, col = 'darkblue') +
  geom_sf(data=citites.dat, shape = 17, col = "olivedrab4") +
  geom_text(data = citites.dat, aes(x = plong, y = plat, label = city),
            size = 3.9, col = "black", fontface = "bold", nudge_x = 0.8) +
  geom_sf(data=filter(map.dat, site_name=='Bainbridge' | site_name=='Union'),
          shape=16,    size = 3, colour="red") +
  geom_sf(data=filter(map.dat, site_name=='Clemson Island' | site_name=='East Donegal Township' |site_name=='Shady Nook'),
          shape='x',    size = 5, colour="green") +
  labs(title="", y="Latitude", x="Longitude") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))



# Inset map
# Box for study area
bb1 <- st_as_sfc(st_bbox(us_state_select2))

ggm1 = ggplot() +
  geom_sf(data = us_states2, fill = "white") +
  geom_sf(data = bb1, fill = NA, color = "red", size = 1.2) +
  theme_void()

# gg_inset_map1 =

# Create final inset map
library(cowplot)
ggdraw() +
  draw_plot(p1) +
  draw_plot(ggm1, x = 0.04, y = 0.65, width = 0.3, height = 0.3)


## --------------- Google API maps

get_map(location = c(lon = -77.85807816278827, lat = 40.80762950782668),
               maptype = "satellite", zoom=14) %>%
  ggmap()


# Sample locations
# Grab bounding box coordinates for use in get_map
# Get bounding box for Susq. Riv Basin
bb_basin <- st_as_sfc(st_bbox(basin))


xvalues <- st_coordinates(bb_basin)[c(1,2),1] # min & max of x values
yvalues <- st_coordinates(bb_basin)[c(2,3),2] # min & max of y values

boundingbox <- c(left = xvalues[1], bottom = yvalues[1], right = xvalues[2], top = yvalues[2]-1)

m1 <- get_map(boundingbox,
        maptype = "terrain", zoom=8, source = "google") %>%
  ggmap() +
  geom_sf(data = map.dat,size = 3,colour="red",
          show.legend = "point",
          inherit.aes = FALSE) +
  geom_sf(data=us_main, col = 'darkblue',
          inherit.aes = FALSE) +
  geom_sf(data=us_main_WB, col = 'darkblue',
          inherit.aes = FALSE) +
  geom_sf(data=ds_main_WB, col = 'darkblue',
          inherit.aes = FALSE) +
  labs(title="", y="Latitude", x="Longitude") +
  theme_bw() +
  #set the CRS for the entire map (sometimes needed)
  coord_sf(crs = st_crs(4326))

m1


# # Grab openstreet waterbody informatio
# # https://www.r-bloggers.com/2018/08/how-to-quickly-enrich-a-map-with-natural-and-anthropic-details/
#
# osm_lakes.sf <-
#   opq(bbox = st_bbox(basin),  timeout = 250) %>%
#   add_osm_feature(key = 'water', value = 'lake') %>%
#   osmdata_sf()
# osm_lakes.sf <- osm_lakes.sf$osm_multipolygons
#
# osm_rivers.sf <-
#   opq(bbox = st_bbox(basin), timeout = 250) %>%
#   add_osm_feature(key = 'waterway', value = 'river') %>%
#   osmdata_sf()
# osm_rivers.sf <- osm_rivers.sf$osm_lines
#
#
# m2 <- get_map(boundingbox,
#               maptype = "toner-lite", zoom=8, source = "stamen") %>%
#   ggmap() +
#   geom_sf(data = map.dat,size = 3,colour="red",
#           show.legend = "point",
#           inherit.aes = FALSE) +
#   geom_sf(data = osm_lakes.sf, fill = '#9ecae1', colour = NA,  inherit.aes = FALSE) +
#   geom_sf(data = osm_rivers.sf, colour = '#9ecae1', size = 0.05, inherit.aes = FALSE) +
#   labs(title="", y="Latitude", x="Longitude") +
#   theme_bw() +
#   #set the CRS for the entire map (sometimes needed)
#   coord_sf(crs = st_crs(4326))
# m2



# Elevation
# set_key("AIzaSyBm043WmwUI0J4JdUo7mWUi45IvMXd14Cc")
# df <- data.frame(lat = c(40.80762950782668,40.831941019195774),
#                  lon = c(-77.85807816278827, -77.97305536595117) )
#
#
# df <- google_elevation(df_locations = df,
#                        location_type = "path",
#                        samples = 20,
#                        simplify = TRUE)
#
# df_plot <- data.frame(elevation = df$results$elevation,
#                       location = as.integer(rownames(df$results)))
#
# ggplot(data = df_plot, aes(x = location, y = elevation)) +
#   geom_line()

