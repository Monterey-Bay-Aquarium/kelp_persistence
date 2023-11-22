##Joshua G. Smith
#May 1, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan, reshape2, raster)


################################################################################
#set directories and load data
localdir <- "/Volumes/seaotterdb$/kelp_recovery/"
basedir <- here::here("output")
figdir <- here::here("figures")

# Read MURSST data
# Source
# https://coastwatch.pfeg.noaa.gov/erddap/griddap/index.html?page=1&itemsPerPage=1000
# Multi-scale Ultra-high Resolution (MUR) SST Analysis fv04.1, Global, 0.01°, 2002-present, Monthly
# data are too large for GitHub -- download and read locally
sst <- raster::brick(file.path(localdir, "data/environmental_data/MURSST/raw/2002_2022_mursst_monthly_raster.grd"))

#load monitoring site table
site_table <- read.csv(file.path(basedir, "monitoring_data/raw/MLPA_kelpforest_site_table.4.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(site, latitude, longitude, ca_mpa_name_short, mpa_class=site_designation, mpa_designation=site_status, 
                baseline_region)%>%
  distinct() #remove duplicates


################################################################################
#intersect monitoring sites with SST

# Convert sites to spatial points
sites <- site_table %>% 
  # Convert to SF
  sf::st_as_sf(coords=c("longitude", "latitude"), crs=crs(sst))


#Extract SST
sst_ts_orig <- extract(x=sst, y=sites, method = "bilinear") #note, bilinear is four-point interpolation

# Format SST extraction
sst_ts <- sst_ts_orig %>% 
  # Convert to df
  as.data.frame() %>% 
  # Add site
  mutate(site=sites$site) %>% 
  dplyr::select(site, everything()) %>% 
  # Gather
  gather(key="date", value="sst_c", 2:ncol(.)) %>% 
  # Format date
  mutate(date=gsub("X", "", date) %>% lubridate::ymd(.)) %>% 
  # Add metadata
  left_join(sites) %>% 
  # Arrange
  dplyr::select(baseline_region, site, everything()) 

#examine output for central coast

sst_NA <- sst_ts %>% filter(baseline_region == "CENTRAL") %>%
            #separate year
          mutate(year = year(as.Date(date))) %>%
          filter(is.na(sst_c))

# Export SST 
saveRDS(sst_ts, file.path(basedir, 
                          "environmental_data/processed/2002_2022_mursst_monthly_by_PISCO_site.Rds"))





