# VAST univariate model with seasonal effects
# modified from https://github.com/James-Thorson-NOAA/VAST/wiki/Index-standardization

# install.packages('TMB', type = 'source')
# remotes::install_github("james-thorson/VAST")

library(dplyr)
library(tidyr)
library(ggplot2)
library(VAST)
# library(furrr)


#
# ## For some reason I need to make sure Rtools has path properly set, else TMB won't compile
# Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# Load Data -----
data("ecomon_epu")

# # # Number of stations per year
total_stations <- ecomon_epu %>%
  group_by(year) %>%
  summarize(total = n_distinct(id))


### Proportion of stations with positive tows per year. Cutoff used to "keep"
spps <- ecomon_epu %>%
  mutate(present = ifelse(abundance > 0,
                          1, 0)) %>%
  filter(year >= 1994,
         spp != "euph1") %>%
  left_join(total_stations) %>%
  group_by(year, spp) %>%
  summarize(positive_stations = sum(present, na.rm = TRUE)/total,
            keep = ifelse(positive_stations > .1,
                          1, 0),
            forage_group = forage_group) %>%
  distinct()

# #
# # ## List of spp where the total number of years by season is greater than 5
# # ## and each spp group has more than cutoff of positive tows
spp_list <- spps %>%
  dplyr::filter(#season %in% season_list,
    forage_group >= 100) %>%
  group_by(spp) %>%
  summarize(count = sum(keep)) %>%
  filter(count > 1) %>%
  select(spp) %>%
  distinct() %>%
  pull(spp)

# n_x = 1000
# td <- vast_wrapper(n_x = 500)

vast_wrapper <- function(n_x = 50) {

  ## Seasonal model -----
  working_dir <- here::here("analysis/vast_EOF3/zoop_ichthyo")

  if(!dir.exists(working_dir)) {
    dir.create(working_dir, recursive  = TRUE)
  }
  #
  ## Attempt to create a log file
  my_log <- file(sprintf("%s/log-%s.txt", working_dir, n_x)) # File name of output log

  sink(my_log, append = TRUE, type = "output") # Writing console output to log file
  on.exit(sink(file = NULL), add = TRUE, after = TRUE)

  sink(my_log, append = TRUE, type = "message")
  on.exit(sink(file = NULL), add = TRUE, after = TRUE)


  zoop_dat <- ecomon_epu %>%
    dplyr::filter(spp %in% spp_list,
                  EPU %in% c("GB", "GOM", "MAB"),
                  # season == !!season,
                  as.numeric(year) >= 1999) %>%
    dplyr::mutate(areaswept_km2 = 1,
                  species_number = as.numeric(factor(spp)) - 1) %>%
    select(species_number,
           spp,
           year,
           abundance,
           areaswept_km2,
           lat,
           lon) %>%
    data.frame()

#
#   ggplot(zoop_dat, aes(x = lon, y = lat)) +
#     geom_point(data = zoop_dat %>% filter(abundance == 0), color = "black", fill = "black", shape = 21) +
#     geom_point(data = zoop_dat %>% filter(abundance > 0), aes(color = season, size = abundance), alpha = 0.5) +
#     facet_wrap(~year) +
#     # labs(title = i) +
#     NULL


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

  # FieldConfig_eof3 <- matrix(c("IID", "Identity", "IID", 2, 0, 0, "IID", "Identity"), ncol = 2, nrow = 4,
  #                            dimnames = list(c("Omega", "Epsilon", "Beta", "Epsilon_year"),
  #                                            c("Component_1", "Component_2")))


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
                "Link"    = 1)

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
                  t_i = zoop_dat$year,
                  c_i = zoop_dat$species_number,
                  b_i = as_units(zoop_dat$abundance, "count"),
                  a_i = as_units(zoop_dat$areaswept_km2, "km^2"),
                  epu_to_use = settings$epu_to_use,
                  working_dir = working_dir,
                  Use_REML = TRUE,
                  # run_model = FALSE,
                  # build_model = TRUE,
                  # test_fit = FALSE,
                  getsd = TRUE,
                  newtonsteps = 1,
                  Options = c('treat_nonencounter_as_zero' = TRUE),
                  optimize_args = list("lower" = -Inf,
                                       "upper" = Inf))


#
#   results = plot( fit,
#                   check_residuals=FALSE,
#                   plot_set=c(3,16) )
#
#   results = plot( fit,
#                   check_residuals=FALSE,
#                   plot_set=c(2,3,6:9,11:16))
#
#   index1_tf = results$Factors$Rotated_loadings$EpsilonTime1
#   cpe_vf = ecodist::vf(index1_tf, data.frame("CPE"=CPE), nperm=1000)
#   png( "year_ordination.png", width=6, height=6, res=200, units="in")
#   plot( index1_tf, type="n" )
#   text( index1_tf, labels=rownames(index1_tf) )
#   plot( cpe_vf )
#   dev.off()

  saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
  # saveRDS(index_ctl_sd_array, file = paste0(working_dir, "/index_ctl_sd_array.rds"))
  return(fit)
}

fit <- vast_wrapper(n_x = 500)


fit$data_list

# Plot results, including spatial term Omega1
results = plot( fit,
                check_residuals=FALSE,
                plot_set=c(3,16),
                category_names = spp_list)

