# VAST EOF model with seasonal effects
# modified from https://github.com/James-Thorson-NOAA/VAST/wiki/

# install.packages('TMB', type = 'source')
# remotes::install_github("james-thorson/VAST")
library(dplyr)
library(tidyr)
library(ggplot2)
library(VAST)

#
# ## For some reason I need to make sure Rtools has path properly set, else TMB won't compile
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# Load Data -----
data("ecomon_epu")

# # # # Number of stations per year
# total_stations <- ecomon_epu %>%
#   group_by(year) %>%
#   summarize(total = n_distinct(id))
#
#
# ### Proportion of stations with positive tows per year. Cutoff used to "keep"
# spps <- ecomon_epu %>%
#   mutate(present = ifelse(abundance > 0,
#                           1, 0)) %>%
#   filter(year >= 1994,
#          spp != "euph1") %>%
#   left_join(total_stations) %>%
#   group_by(year, spp) %>%
#   summarize(positive_stations = sum(present, na.rm = TRUE)/total,
#             keep = ifelse(positive_stations > .1,
#                           1, 0),
#             forage_group = forage_group) %>%
#   distinct()
#
# # #
# # # ## List of spp where the total number of years by season is greater than 5
# # # ## and each spp group has more than cutoff of positive tows
# spp_list <- spps %>%
#   dplyr::filter(#season %in% season_list,
#     forage_group >= 100) %>%
#   group_by(spp) %>%
#   summarize(count = sum(keep)) %>%
#   filter(count > 1) %>%
#   select(spp) %>%
#   distinct() %>%
#   pull(spp)

spp_list <- c("calfin", "chaeto", "cham", "clauso", "ctyp",
              "euph", "gas", "hyper", "larvaceans",
              "mlucens", "oithspp", "para", "pseudo", "tlong")
n_x = 50

## Seasonal model -----
working_dir <- here::here("analysis/vast_seasonal_EOF3/")

if(!dir.exists(working_dir)) {
  dir.create(working_dir, recursive  = TRUE)
}
#
# ## Attempt to create a log file
# my_log <- file(sprintf("%s/%s_log-%s.txt", working_dir, sp, n_x)) # File name of output log
#
# sink(my_log, append = TRUE, type = "output") # Writing console output to log file
# on.exit(sink(file = NULL), add = TRUE, after = TRUE)
#
# sink(my_log, append = TRUE, type = "message")
# on.exit(sink(file = NULL), add = TRUE, after = TRUE)
#
year_set = 1994:2017
season_set = c("winter", "spring", "summer", "fall")
#
# zoop_dat <- ecomon_epu %>%
#   dplyr::filter(spp %in% spp_list,
#                 EPU %in% c("GB", "GOM", "MAB"),
#                 as.numeric(year) >= 1994) %>%
#   # right_join(expand_grid(year = year_set,
#   #                       season = season_set,
#   #                       spp = spp_list)) %>%
#   dplyr::mutate(species_number = as.numeric(factor(spp)),
#                 areaswept_km2 = 1,
#                 year_season = factor( paste(year, season, sep = "_"))) %>% #,
#                 # lat = ifelse(is.na(lat),
#                 #              mean(lat, na.rm = TRUE),
#                 #              lat),
#                 # lon = ifelse(is.na(lon),
#                 #              mean(lon, na.rm = TRUE),
#                 #              lon)) %>%
#   select(species_number,
#          spp,
#          year,
#          year_season,
#          season,
#          abundance,
#          areaswept_km2,
#          lat,
#          lon) %>%
#   data.frame()

zoop_dat <- ecomon_epu %>%
  dplyr::filter(spp %in% spp_list,
                EPU %in% c("GB", "GOM", "MAB"),
                as.numeric(year) >= 1994) %>%
  dplyr::mutate(areaswept_km2 = 1,
                year_season = factor( paste(year, season, sep = "_")),
                species_number = as.numeric(factor(spp)) - 1) %>%
  select(species_number,
         spp,
         year,
         year_season,
         abundance,
         areaswept_km2,
         lat,
         lon) %>%
  data.frame()



# total_stations <- ecomon_epu %>%
#   dplyr::filter(spp %in% spp_list,
#                 EPU %in% c("GB", "GOM", "MAB"),
#                 as.numeric(year) >= 2010) %>%
#   right_join(expand_grid(year = year_set,
#                          season = season_set,
#                          spp = spp_list)) %>%
#   mutate(year_season = factor(paste(year, season, sep = "_")),
#          present = ifelse(abundance > 0, 1, 0)) %>%
#   # select(year, season, spp, present) %>%
#   group_by(year_season, spp) %>%
#   summarize(total = n_distinct(id),
#             prop = sum(present)/total)

# Load data and quick exploration of structure
# Set of years and seasons. Order accordingly.
#
# # Create a grid with all unique combinations of seasons and years and then combine these into one "year_season" variable
# yearseason_grid = expand.grid("season" = season_set, "year" = year_set)
# yearseason_levels = apply(yearseason_grid[,2:1], MARGIN = 1, FUN = paste, collapse = "_")
# yearseason_labels = round(yearseason_grid[,'year'] + (as.numeric(factor(yearseason_grid[,'season'],
#                                                                         levels = season_set))-1)/length(season_set), digits=1)
#
# # Similar process, but for the observations
# yearseason_i = apply(zoop_dat[,c("year","season")], MARGIN = 1, FUN = paste, collapse = "_")
# yearseason_i = factor(yearseason_i, levels = yearseason_levels)
#
# # Add the year_season factor column to our sampling_data data set
# zoop_dat$year_season = yearseason_i
# zoop_dat$season = factor(zoop_dat$season, levels = season_set)


# ggplot(zoop_dat, aes(x = lon, y = lat)) +
#   geom_point(data = zoop_dat %>% filter(abundance == 0), color = "black", fill = "black", shape = 21) +
#   geom_point(data = zoop_dat %>% filter(abundance > 0), aes(color = spp, size = abundance), alpha = 0.5) +
#   facet_wrap(~ year_season, ncol = 4) +
#   # labs(title = i) +
#   NULL


# Some last processing steps
# zoop_dat = zoop_dat[, c("species_number", "spp", "year", "season", "year_season", "lat", "lon", "areaswept_km2", "abundance")]

# # Make dummy observation for each season-year combination
# dummy_data = data.frame(
#   # species_number = zoop_dat[, "species_number"],
#   # spp = zoop_dat[, "spp"],
#   year = yearseason_grid[,'year'],
#   season = yearseason_grid[,'season'],
#   year_season = yearseason_levels,
#   lat = mean(zoop_dat[,'lat']),
#   lon = mean(zoop_dat[,'lon']),
#   areaswept_km2 = mean(zoop_dat[,'areaswept_km2']),
#   abundance = 0,
#   dummy = TRUE) %>%
#   expand_grid(spp = unique(zoop_dat[, "spp"])) %>%
#   mutate(species_number =  as.numeric(factor(spp)) - 1)
#
# # Combine with sampling data
# full_data = rbind(cbind(zoop_dat, dummy = FALSE), dummy_data)
#
# # Create sample data
# samp_dat = data.frame(
#   "spp" =  full_data$spp,
#   "species_number" = full_data$species_number,
#   "year_season" = as.numeric(full_data$year_season)-1,
#   "Lat" = full_data$lat,
#   "Lon" = full_data$lon,
#   "abundance" = full_data$abundance,
#   "areaswept_km2" = full_data$areaswept_km2,
#   "Dummy" = full_data$dummy )
#
# # Covariate data. Note here, case sensitive!
# cov_dat = data.frame(
#   "Year" = as.numeric(full_data$year_season)-1,
#   "Year_Cov" = factor(full_data$year, levels = year_set),
#   "Season" = full_data$season,
#   "Lat" = full_data$lat,
#   "Lon" = full_data$lon)

# Inspect
# table("year_season"=cov_dat$Year, "Actual_year"=cov_dat$Year_Cov)
# table("year_season"=cov_dat$Year, "Actual_season"=cov_dat$Season)

#####
## Model settings
#####


##  Random Fields -----

## Control the random fields part of the model.
## Omega = X is the number of random spatial fields to apply
## and Epsilon = X is the number of random spatio-temporal
## fields to apply. Omega1 is for the probability of
## occurrence, and Omega2 is for the density given occurrence,
## similarly for Epsilon.

## 0 = off
## "AR1" = AR1 process
## >1 = number of elements in a factor-analysis covariance
## "IID" = random effect following an IID distribution

FieldConfig <- c("Omega1" = "IID", "Epsilon1" = "IID",
                 "Omega2" = "IID", "Epsilon2" = "IID")


## Autoregressive structure -----

## Control autoregressive structure for parameters over time
## Changing the settings here creates different
## autoregressive models for the intercept (Beta) and
## spatio-temporal process (Epsilon).

## 0 = each year is a fixed effect
## 1 = random effect
## 2 = random walk
## 3 = fixed effect that is constant over time
## 4 = AR1 process

RhoConfig <- c("Beta1"    = c(0, 1, 2, 3, 4)[1],
               "Beta2"    = c(0, 1, 2, 3, 4)[1],
               "Epsilon1" = c(0, 1, 2, 3, 4)[1],
               "Epsilon2" = c(0, 1, 2, 3, 4)[1])

## Correlated overdispersion -----

## Control correlated overdispersion among categories
## for each level of v_i, where eta1 is for encounter
## probability, and eta2 is for positive catch rates
# eta1 = vessel effects on prey encounter rate
# eta2 = vessel effects on prey weight

## 0 = off,
## "AR1" = AR1 process,
## >0 = number of elements in a factor-analysis covariance
OverdispersionConfig <- c("eta1" = 0,
                          "eta2" = 0)


## Observation model -----
# Control observation model structure. The first
# component sets the distribution of the positive
# distribution component. ?VAST::make_data()

ObsModel <- c("PosDist" = 1, # Delta-Gamma; Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number
              "Link"    = 3)

# Make settings
settings = make_settings(n_x = n_x,
                         Region = "northwest_atlantic",
                         strata.limits = "EPU",
                         purpose = "EOF3",
                         n_categories = 2,
                         # FieldConfig = FieldConfig,
                         RhoConfig = RhoConfig,
                         bias.correct = TRUE,
                         Options = c('treat_nonencounter_as_zero' = TRUE),
                         ObsModel = ObsModel)#,

settings$epu_to_use <- c("Georges_Bank", "Gulf_of_Maine", "Mid_Atlantic_Bight")

#####
## Model fit -- make sure to use new functions
#####

fit = fit_model(settings = settings,
                Lat_i = zoop_dat$lat,
                Lon_i = zoop_dat$lon,
                # t_i = zoop_dat$year,
                t_i = as.numeric(zoop_dat$year_season) - 1,
                c_i = zoop_dat$species_number,
                b_i = as_units(zoop_dat$abundance, "count"),
                a_i = as_units(zoop_dat$areaswept_km2, "km^2"),
                epu_to_use = settings$epu_to_use,
                working_dir = working_dir,
                Use_REML = FALSE,
                getsd = FALSE,
                newtonsteps = 1,
                Options = c('treat_nonencounter_as_zero' = TRUE),
                optimize_args = list("lower" = -Inf,
                                     "upper" = Inf))
# Creating model formula
# X1_formula = ~ Season + Year_Cov
# X2_formula = ~ Season + Year_Cov

# X1_formula = ~ Season #+ Year_Cov
# X2_formula = ~ Season #+ Year_Cov
#

# Implement corner constraint for linear effect but not spatially varying effect:
# * one level for each term is 2 (just spatially varying) spatially varying, zero-centered linear effect on 1st linear predictor for category c
# * all other levels for each term is 3 (spatially varying plus linear effect), spatially varying linear effect on 1st linear predictor for category c
# X1config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat$Season)-1), 2, rep(3,nlevels(cov_dat$Year_Cov)-1) ), nrow=1 )
# X2config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat$Season)-1), 2, rep(3,nlevels(cov_dat$Year_Cov)-1) ), nrow=1 )

# X1config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat$Season)-1)), nrow=1 )
# X2config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat$Season)-1)), nrow=1 )

# X1config_cp_use = matrix( rep(c(2, rep(3, nlevels(cov_dat$Season)-1)), length(spp_list)), nrow = length(spp_list), byrow = TRUE)
# X2config_cp_use = matrix( rep(c(2, rep(3, nlevels(cov_dat$Season)-1)), length(spp_list)), nrow = length(spp_list), byrow = TRUE)
#


#####
## Model fit -- make sure to use new functions
# #####
#
# fit = fit_model(settings = settings,
#                 Lat_i = zoop_dat$lat,
#                 Lon_i = zoop_dat$lon,
#                 # t_i = as.numeric(zoop_dat$year_season) - 1,
#                 t_i = zoop_dat$year,
#                 c_i = zoop_dat$species_number,
#                 b_i = as_units(zoop_dat$abundance, "count"),
#                 a_i = as_units(zoop_dat$areaswept_km2, "km^2"),
#                 epu_to_use = settings$epu_to_use,
#                 working_dir = working_dir,
#                 # X1config_cp = X1config_cp_use,
#                 # X2config_cp = X2config_cp_use,
#                 # covariate_data = cov_dat,
#                 # X1_formula = X1_formula,
#                 # X2_formula = X2_formula,
#                 # X_contrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE)#,
#                 # Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)
#                 # ),
#                 Use_REML = FALSE,
#                 # PredTF_i = samp_dat$Dummy,
#                 getsd = FALSE,
#                 newton_steps = 0,
#                 Options = c('treat_nonencounter_as_zero' = TRUE),
#                 optimize_args = list("lower" = -Inf,
#                                      "upper" = Inf))



# Adjust mapping for log_sigmaXi and fitting final model -- pool variance for all seasons and then set year's to NA
Map_adjust = fit_orig$tmb_list$Map

# Pool variances for each term to a single value
# Map_adjust$log_sigmaXi1_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), nlevels(cov_dat$Season)),
#                                       rep(as.numeric(Map_adjust$log_sigmaXi1_cp[nlevels(cov_dat$Season)+1]), nlevels(cov_dat$Year_Cov))))
# Map_adjust$log_sigmaXi2_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), nlevels(cov_dat$Season)),
# rep(as.numeric(Map_adjust$log_sigmaXi2_cp[nlevels(cov_dat$Season)+1]), nlevels(cov_dat$Year_Cov))))
#
#

# Map_adjust$log_sigmaXi1_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), nlevels(cov_dat$Season))))#,
# Map_adjust$log_sigmaXi2_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), nlevels(cov_dat$Season))))#,

Map_adjust$log_sigmaXi1_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), 20)))#,
Map_adjust$log_sigmaXi2_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), 20)))#,


# Fit final model with new mapping
fit = fit_model(settings = settings,
                Lat_i = samp_dat$Lat,
                Lon_i = samp_dat$Lon,
                t_i = samp_dat$year_season,
                c_i = samp_dat$species_number,
                b_i =  as_units(samp_dat$abundance, "count"),
                a_i = as_units(samp_dat$areaswept_km2, "km^2"),
                epu_to_use = settings$epu_to_use,
                working_dir = working_dir,
                X1config_cp = X1config_cp_use,
                X2config_cp = X2config_cp_use,
                covariate_data = cov_dat,
                X1_formula = X1_formula,
                X2_formula = X2_formula,
                X_contrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE)#,
                                   # Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)
                ),
                Use_REML = TRUE,
                PredTF_i = samp_dat$Dummy,
                newtonsteps = 1,
                run_model = TRUE,
                Map = Map_adjust,
                getsd = FALSE,
                Options = c('treat_nonencounter_as_zero' = TRUE),
                optimize_args = list("lower" = -Inf,
                                     "upper" = Inf))

# Calculate predictive distribution for a Index_ctl
index_ctl_sd_array <- sample_variable(Sdreport = fit$par$SD,
                                      Obj = fit$tmb_list$Obj,
                                      variable_name = "Index_ctl",
                                      n_samples = 100,
                                      sample_fixed = FALSE) ## Only sample random effects

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
saveRDS(index_ctl_sd_array, file = paste0(working_dir, "/index_ctl_sd_array.rds"))



possibly_vast_wrapper <- purrr::quietly(vast_wrapper)

# plan(multisession, workers = 4)
vast_runs <- purrr::map(.x = spp_list, .f = possibly_vast_wrapper)
# plan(sequential)