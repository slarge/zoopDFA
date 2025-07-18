library(dplyr)
library(tidyr)
library(marmap)
# library(googlesheets4)
# remotes::install_github("NOAA-EDAB/ecodata")

source("R/helper_functions.R")

## Right now the most up to date survdat is on ECSA. Soon this will migrate to ecodata
load(url("https://github.com/NOAA-EDAB/ECSA/blob/master/data/Survdat.RData?raw=true"))
usethis::use_data(survdat, overwrite = TRUE)

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
## Download and reformat EcoMon data ##
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# download.file(url = "https://www.nodc.noaa.gov/archive/arc0143/0187513/3.3/data/0-data/EcoMon_Plankton_Data_v3_8.csv",
#               destfile = "data-raw/EcoMon_Plankton_data_v3_8.csv", mode = "wb")


ecomon_long <- readr::read_csv("data-raw/EcoMon_Plankton_data_v3_10.csv") %>%
  rename_with(stringr::str_to_lower) %>%
  mutate(volume = as.double(volume_100m3),
         date = as.Date(date, format="%d-%b-%y")) %>%
  select(-ends_with("_10m2"),
         -ends_with("_10m2x"),
         -starts_with("volume_")) %>%
  pivot_longer(cols = c(ends_with("_100m3")),
               names_to = "spp",
               values_to = "abundance") %>%
  rename(station = staton)

usethis::use_data(ecomon_long, overwrite = TRUE)

# These oblique tows include a measure of total volume swept, and we divide the total number of
# zoop by volume swept and then multiply by the seafloor depth at the beginning of the tow to
# obtain vertically integrated numbers-density.
# Using vertically integrated numbers-density as response variable then allows us to predict
# vertically integrated densities across a standard survey area, where the sum across this survey
# area represents a prediction of vertical and spatially integrated abundance in numbers.


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
## Prepare EcoMon data for VAST ##
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# download.file(url = "https://www.nodc.noaa.gov/archive/arc0143/0187513/3.3/data/0-data/EcoMon_Plankton_Data_v3_8.csv",
#               destfile = "data-raw/EcoMon_Plankton_data_v3_8.csv", mode = "wb")


data("ecomon_long")
#

# prepare bottom depth
xlims <- c(-77, -65)
ylims <- c(35, 45)
res <- 1
bath_filename <- sprintf("marmap_coord_%s;%s;%s;%s_res_%s.csv",
                         xlims[1], ylims[1], xlims[2], ylims[2], res)

if(!bath_filename %in% list.files("data-raw")){
  nesbath <- marmap::getNOAA.bathy(lon1 = xlims[1], lon2 = xlims[2],
                                   lat1 = ylims[1], lat2 = ylims[2],
                                   resolution = res,
                                   keep = TRUE) %>%
    marmap::as.raster()

  file.copy(bath_filename, "data-raw")
  file.remove(bath_filename)
} else {
  nesbath <- marmap::read.bathy(sprintf("data-raw/%s", bath_filename), header = T) %>%
    marmap::as.raster()
}

ecomon_format <- ecomon_long %>%
  mutate(marmap_depth = raster::extract(nesbath, y = cbind(.$lon, .$lat)) * -1,
         depth = ifelse(depth == 9999,
                        marmap_depth,
                        depth),
         spp = gsub("_100m3", "", spp),
         date = as.Date(date, format = "%d-%b-%y"),
         id = as.factor(paste0(cruise_name, "_", station)),
         vessel = as.factor(substr(cruise_name, start = 1, stop = 2)),
         areaswept_km2 = .001,
         day = as.numeric(strftime(date, format = "%j")),
         year = as.numeric(strftime(date, format = "%Y")),
         season = case_when(day %in% 1:90 ~ "winter",
                            day %in% 91:181 ~ "spring",
                            day %in% 182:273 ~ "summer",
                            day %in% 274:365 ~ "fall",
                            TRUE ~ NA_character_)) %>%
  dplyr::filter(grepl("6B3", zoo_gear),
                !grepl("6B5", ich_gear))

crs_epu <- 4269 # NAD83 https://epsg.org/crs_4269/NAD83.html
# crs_epu <- 9311 # NAD27 https://epsg.org/crs_9311/NAD27-US-National-Atlas-Equal-Area.html

epu <- sf::read_sf(here::here("data-raw/EPU_NOESTUARIES.shp")) %>%
# epu <- ecodata::epu_sf %>%
  sf::st_transform(crs = sf::st_crs(crs_epu))

sf::sf_use_s2(FALSE)

## Post stratify data according to EPUs
ecomon_epu <- ecomon_format %>%
  sf::st_as_sf(coords = c("lon","lat"), crs = crs_epu) %>%
  sf::st_join(epu) %>%
  # sf::st_join(ecodata::epu_sf) %>%
  sfc_as_cols(names = c("lon", "lat")) %>%
  sf::st_drop_geometry() %>%
  select(id,
         spp,
         EPU,
         day,
         season,
         bottom_depth = depth,
         vessel,
         abundance,
         areaswept_km2,
         lat,
         lon,
         year,
         zoo_gear,
         ich_gear) %>%
  data.frame()

### Bring in the forage taxa info
forage_taxa <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1z0TkgN3XYQ9NCF2r0CLTYnJRFog74-fSjXHpAfbt-UI/edit?usp=sharing") %>%
  mutate(spp = gsub("_10m2|_abnd", "", `COLUMN NAME`),
         sciname = stringr::str_to_sentence(`TAXA NAME`),
         sciname = gsub(" - append", "", sciname),
         forage_group = ifelse(is.na(IchGroup), ZooGroup, IchGroup + 100)) %>%
  select(sciname, spp, forage_name = `Forage Name`, forage_group)

ecomon_epu <- ecomon_epu %>%
  left_join(forage_taxa)

usethis::use_data(ecomon_epu, overwrite = TRUE)


## Post stratify data according to BTS strata
crs_strata <- 4269 # NAD83 https://epsg.org/crs_4269/NAD83.html
# crs_epu <- 9311 # NAD27 https://epsg.org/crs_9311/NAD27-US-National-Atlas-Equal-Area.html

strata <- sf::st_read(dsn = here::here("data-raw/NES_BOTTOM_TRAWL_STRATA.shp")) %>%
  sf::st_transform(crs = sf::st_crs(crs_strata))

sf::sf_use_s2(FALSE)

## Post stratify data according to EPUs
ecomon_strata <- ecomon_format %>%
  sf::st_as_sf(coords = c("lon","lat"), crs = crs_strata) %>%
  sf::st_join(strata) %>%
  sfc_as_cols(names = c("lon", "lat")) %>%
  sf::st_drop_geometry() %>%
  select(id,
         spp,
         strata = STRATA,
         day,
         season,
         bottom_depth = depth,
         vessel,
         abundance,
         areaswept_km2,
         lat,
         lon,
         year,
         zoo_gear,
         ich_gear) %>%
  data.frame()

### Bring in the forage taxa info
forage_taxa <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1z0TkgN3XYQ9NCF2r0CLTYnJRFog74-fSjXHpAfbt-UI/edit?usp=sharing") %>%
  mutate(spp = gsub("_10m2|_abnd", "", `COLUMN NAME`),
         sciname = stringr::str_to_sentence(`TAXA NAME`),
         sciname = gsub(" - append", "", sciname),
         forage_group = ifelse(is.na(IchGroup), ZooGroup, IchGroup + 100)) %>%
  select(sciname, spp, forage_name = `Forage Name`, forage_group)

ecomon_strata <- ecomon_strata %>%
  left_join(forage_taxa)

usethis::use_data(ecomon_strata, overwrite = TRUE)


### PANGAEA Copepod data

tmp <- here::here("data-raw", "Brun-etal_2016_Copepode_trait.xlsx")
download.file(url = "https://store.pangaea.de/Publications/BrunP-etal_2016/Brun-etal_2016_Copepode_trait.xlsx", destfile = tmp)

copepod_traits <- readxl::read_xlsx(tmp, sheet = "Body size")
usethis::use_data(copepod_traits, overwrite = TRUE)
