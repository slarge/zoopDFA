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
# Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# setwd(R'(C:\Users\James.Thorson\Desktop\Git\zoopDFA)')

# Load Data -----
data("ecomon_epu")


#
# # # # # Number of stations per year
# total_stations <- ecomon_epu %>%
#   group_by(year) %>%
#   summarize(total = n_distinct(id))
# #
# #
# # ### Proportion of stations with positive tows per year. Cutoff used to "keep"
# spps <- ecomon_epu %>%
#   mutate(present = ifelse(abundance > 0,
#                           1, 0)) %>%
#   filter(year >= 1994,
#          spp != "euph1") %>%
#   left_join(total_stations) %>%
#   group_by(year, season, spp) %>%
#   summarize(positive_stations = sum(present, na.rm = TRUE)/total,
#             keep = ifelse(positive_stations > .1,
#                           1, 0),
#             forage_group = forage_group) %>%
#   distinct()

# #
# # ## List of spp where the total number of years by season is greater than 5
# # ## and each spp group has more than cutoff of positive tows
# spp_list <- spps %>%
#   dplyr::filter(#season %in% season_list,
#     forage_group < 100) %>%
#   group_by(spp) %>%
#   summarize(count = sum(keep)) %>%
#   filter(count > 26) %>%
#   select(spp) %>%
#   distinct() %>%
#   pull(spp)

# spp_list <- c("acarspp", "calfin", "chaeto", "cham", "cirr", "clauso", "ctyp",
#               "echino", "euph", "evadnespp", "gas", "hyper", "larvaceans",
#               "mlucens", "oithspp", "para", "penilia", "pseudo", "salps", "tlong")
#

spp_list <- c("calfin", "chaeto", "cham", "clauso", "ctyp",
              "euph", "gas", "hyper", "larvaceans",
              "mlucens", "oithspp", "para",  "pseudo", "tlong")


# spp_list <- c("calfin", "cham", "ctyp", "tlong")
# n_x = 500

vast_wrapper <- function(n_x = 50){
  ## Seasonal model -----
  working_dir <- here::here("analysis/vast_seasonal_EOF3/")

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
#                  EPU %in% c("GB", "GOM", "MAB"),
#                  as.numeric(year) >= 2010,
#                  as.numeric(year) < 2015) %>%
                  as.numeric(year) >= 1992,
                  as.numeric(year) <= 2021,
                  EPU %in% c("GB", "GOM", "MAB")) %>%
    dplyr::mutate(areaswept_km2 = 1,
                  year_level = factor(year),
                  season = factor(season, levels = c("winter", "spring", "summer", "fall")),
                  year_labels = factor(paste(year, season, sep = "_")),
                  year_season = factor(year_labels, levels = paste(rep(levels(year_level),
                                                                       each = nlevels(season)),
                                                                   levels(season),
                                                                   sep="_")),
                  species_number = as.numeric(factor(spp)) - 1) %>%
    select(year,
           year_season,
           season,
           lat,
           lon,
           areaswept_km2,
           species_number,
           spp,
           abundance) %>%
    droplevels() %>%
    # mutate(t_i = as.numeric(year_season)- 1) %>%
    arrange(year_season) %>%
    data.frame()


  # table( zoop_dat$year, zoop_dat$season )

  # ggplot(zoop_dat, aes(x = lon, y = lat)) +
  #   geom_point(data = zoop_dat %>% filter(abundance == 0), color = "black", fill = "black", shape = 21) +
  #   geom_point(data = zoop_dat %>% filter(abundance > 0), aes(color = spp, size = abundance), alpha = 0.5) +
  #   facet_wrap(~ year_season, ncol = 4) +
  #   # labs(title = i) +
  #   NULL


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

  # FieldConfig <- c("Omega1" = "IID", "Epsilon1" = "IID",
  #                  "Omega2" = "IID", "Epsilon2" = "IID")

  FieldConfig = matrix(c("IID", "Identity", "IID", 2, 0, 0, "IID", "Identity"), ncol = 2, nrow = 4,
                       dimnames = list(c("Omega", "Epsilon", "Beta", "Epsilon_year"),
                                       c("Component_1", "Component_2")))
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
#
  ObsModel <- c("PosDist" = 1, # Delta-Gamma; Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number
                "Link"    = 4)

  # ObsModel <- c("PosDist" = 1, # Delta-Gamma; Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number
  #               "Link"    = 3)

  # Make settings
  settings = make_settings(n_x = n_x,
                           Region = "northwest_atlantic",
                           strata.limits = "EPU",
                           Version = "VAST_v14_0_1",
                           purpose = "EOF3",
                           n_categories = 2,
                           ObsModel = ObsModel,
                           RhoConfig = RhoConfig#,
                           # use_anisotropy = FALSE,
                           # FieldConfig = FieldConfig,
                           # bias.correct = FALSE,
                           # Options = c('treat_nonencounter_as_zero' = FALSE)
                           )

  settings$epu_to_use <- c("Georges_Bank", "Gulf_of_Maine", "Mid_Atlantic_Bight")

  #####
  ## Model fit -- make sure to use new functions
  #####

  # To run the model with missing year_seasons, need to get the indexing right in the spatially varying season indicator
  # it wouldn't run for season-year combos that have no data, but you can also drop those with essentially no impact on results
  # as long as you then do the plotting right, and get the season-indicator Q term inputted correctly
  # it's just a matter of getting the covariate_data right for earlier years with missing season-year combos ...
  # I think it's getting the Year column (which actually represents Season-Year) aligned with the right Season

  cov_dat <- zoop_dat %>%
    mutate(Year = as.numeric(year_season) - 1) %>%
    select(Lat = lat,
           Lon = lon,
           season,
           Year) %>%
    data.frame()

  # fit_og <- readRDS(paste0(working_dir, "/fit.rds"))
    # Don't use X_contrasts, so that fixed season-year slope isn't confounded with beta term
  fit = fit_model(settings = settings,
                  Lat_i = zoop_dat$lat,
                  Lon_i = zoop_dat$lon,
                  # t_i = zoop_dat$t_i,
                  t_i = as.numeric(zoop_dat$year_season) - 1,
                  c_i = zoop_dat$species_number,
                  b_i = as_units(zoop_dat$abundance, "count"),
                  a_i = as_units(zoop_dat$areaswept_km2, "km^2"),
                  epu_to_use = settings$epu_to_use,
                  newtonsteps = 0,
                  covariate_data = cov_dat,
                  X1_formula = ~ season,
                  X1config_cp = matrix(2, nrow = length(spp_list), ncol = 4),
                  X_contrasts = list(season = contrasts(zoop_dat$season, contrasts = FALSE)),
                  # Parameters = fit_og$ParHat,
                  #hessian.fixed = H,
                  getsd = TRUE,
                  Use_REML = TRUE,
                  run_model = TRUE,
                  working_dir = paste0(working_dir, "/"),
                  optimize_args = list("lower" = -Inf,
                                       "upper" = Inf))
                  # category_names = spp_list,
                  # year_labels = levels(zoop_dat$year_season))

  results = plot_results( fit,
                          check_residuals = FALSE,
                          plot_set= c(3,14,16,18), #c(3,16,18),
                          working_dir = "analysis/vast_seasonal_EOF3/figures/",
                          category_names = spp_list,
                          Version = "VAST_v14_0_1",
                          # year_labels = levels(zoop_dat$year_season),
                          strata_names =  c("Georges Bank", "Gulf of Maine", "Mid-Atlantic Bight"))


  saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
  saveRDS(results, file = paste0(working_dir, "/results.rds"))

  return(fit)
}

## Dual ordination  ##
# settings$FieldConfig["Epsilon",1] = 4

working_dir <- here::here("analysis/vast_seasonal_EOF3/")
fit <- readRDS(paste0(working_dir, "/fit.rds"))

vast_runs <- vast_wrapper(n_x = 500)



MapDetails_List = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list)

MapDetails<-MapDetails_List[["MappingDetails"]]

fit$year_labels <- levels(zoop_dat$year_season)

expand_grid(year = 1992:2021, season = c("winter", "spring", "summer", "fall")) %>%
  mutate(year_season = paste(year, season, sep = "_")) %>%
  pull(year_season)



ffit <- plot_factors(fit = fit, Report = fit$Report, ParHat = fit$ParHat, Data = fit$data_list,
             SD = fit$parameter_estimates$SD,
             # category_names = spp_list,
             # mapdetails_list = MapDetails_List,
             # n_cells = 2000,
             plotdir = paste0(working_dir, "/"))

# ffit <- plot_factors(fit, working_dir = paste0(working_dir, "/")) #)c(3,14,16,18) )
#
td <- data.frame(factor(levels(zoop_dat$year_season), levels = levels(zoop_dat$year_season)), ffit$Rotated_loadings$EpsilonTime1)
colnames(td) <- c("year_season", "Factor_1", "Factor_2")


## Reimagine the example data as seasonal sampling instead of annual sampling
season_dat <- expand.grid(year = min(zoop_dat$year):max(zoop_dat$year),
                          season = factor(c("winter", "spring", "summer", "fall"), levels = c("winter", "spring", "summer", "fall"))) %>%
  mutate(year_level = factor(year),
         season = factor(season, levels = c("winter", "spring", "summer", "fall")),
         year_labels = factor(paste(year, season, sep = "_")),
         year_season = factor(year_labels, levels = paste(rep(levels(year_level),
                                                              each = nlevels(season)),
                                                          levels(season),
                                                          sep="_")))


season_td <- td %>%
  right_join(season_dat) %>%
  pivot_longer(cols = starts_with("Factor"), names_to = "var", values_to = "val") %>%
  arrange(year_season) %>%
  mutate(year_season = factor(year_season, levels = levels(year_labels)),
         val = val * -1)

# levels(season_td$tt)
#
# f1_ts <- decompose(ts(data = td$Factor_1, start = c(1994, 1), end = c(2015, 4), frequency = 4))$trend
# f2_ts <- decompose(ts(data = td$Factor_2, start = c(1994, 1), end = c(2015, 4), frequency = 4))$trend %>%
#   as.data.frame() %>%
#   mutate(x = ifelse(is.na(x), 0, x))

#
ggplot(season_td, aes(x = year_season, xend = year_season, y = val, color = season)) +
  geom_point() +
  geom_segment( aes(y=0, yend = val)) +
  facet_wrap(~var, nrow = 2) +
  scale_color_brewer(palette = "Set1") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_minimal() +
  NULL


# possibly_vast_wrapper <- purrr::quietly(vast_wrapper)

# plan(multisession, workers = 4)
# vast_runs <- possibly_vast_wrapper(n_x = 750)

# fit <- vast_runs$result$tmb_list

# fit <- readRDS(here::here("analysis/vast_seasonal_EOF3/fit.rds"))

results = plot_results( fit,
                check_residuals = FALSE,
                plot_set= c(3,14,16,18), #c(3,16,18),
                working_dir = "analysis/vast_seasonal_EOF3/",
                category_names = spp_list,
                Version = "VAST_v14_0_1",
                year_labels = levels(zoop_dat$year_season),
                strata_names =  c("Georges Bank", "Gulf of Maine", "Mid-Atlantic Bight"))




# use `ecodist` to display ordination
index1_tf = results$Factors$Rotated_loadings$EpsilonTime1
ft_ts_vf = ecodist::vf(index1_tf, data.frame("F2_trend" = data.frame(f2_ts)), nperm=1000)
png( "year_ordination.png", width=6, height=6, res=200, units="in")
plot( index1_tf, type="n" )
text( index1_tf, labels=rownames(index1_tf) )
plot( ft_ts_vf )
dev.off()



# Plot against cold-pool extent index
index2 = results$Factors$Rotated_loadings$EpsilonTime1[,2]
index2 = sign(cor(index2,as.numeric(f2_ts$x))) * index2
png( "EOF_index.png", width=6, height=6, res=200, units="in")
matplot( x=fit$year_labels, y=scale(cbind(f2_ts$x,index2)),
         type="l", lty="solid", col=c("blue","black"), lwd=2, ylab="Index", xlab="Year" )
legend( "bottom", ncol=2, fill=c("blue","black"), legend=c("Trend","factor-2"), bty="n")
dev.off()


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

#
#
# # Adjust mapping for log_sigmaXi and fitting final model -- pool variance for all seasons and then set year's to NA
# Map_adjust = fit_orig$tmb_list$Map
#
# # Pool variances for each term to a single value
# # Map_adjust$log_sigmaXi1_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), nlevels(cov_dat$Season)),
# #                                       rep(as.numeric(Map_adjust$log_sigmaXi1_cp[nlevels(cov_dat$Season)+1]), nlevels(cov_dat$Year_Cov))))
# # Map_adjust$log_sigmaXi2_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), nlevels(cov_dat$Season)),
# # rep(as.numeric(Map_adjust$log_sigmaXi2_cp[nlevels(cov_dat$Season)+1]), nlevels(cov_dat$Year_Cov))))
# #
# #
#
# # Map_adjust$log_sigmaXi1_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), nlevels(cov_dat$Season))))#,
# # Map_adjust$log_sigmaXi2_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), nlevels(cov_dat$Season))))#,
#
# Map_adjust$log_sigmaXi1_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), 20)))#,
# Map_adjust$log_sigmaXi2_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), 20)))#,
#
#
# # Fit final model with new mapping
# fit = fit_model(settings = settings,
#                 Lat_i = samp_dat$Lat,
#                 Lon_i = samp_dat$Lon,
#                 t_i = samp_dat$year_season,
#                 c_i = samp_dat$species_number,
#                 b_i =  as_units(samp_dat$abundance, "count"),
#                 a_i = as_units(samp_dat$areaswept_km2, "km^2"),
#                 epu_to_use = settings$epu_to_use,
#                 working_dir = working_dir,
#                 X1config_cp = X1config_cp_use,
#                 X2config_cp = X2config_cp_use,
#                 covariate_data = cov_dat,
#                 X1_formula = X1_formula,
#                 X2_formula = X2_formula,
#                 X_contrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE)#,
#                                    # Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)
#                 ),
#                 Use_REML = TRUE,
#                 PredTF_i = samp_dat$Dummy,
#                 newtonsteps = 1,
#                 run_model = TRUE,
#                 Map = Map_adjust,
#                 getsd = FALSE,
#                 Options = c('treat_nonencounter_as_zero' = TRUE),
#                 optimize_args = list("lower" = -Inf,
#                                      "upper" = Inf))
#
# # Calculate predictive distribution for a Index_ctl
# index_ctl_sd_array <- sample_variable(Sdreport = fit$par$SD,
#                                       Obj = fit$tmb_list$Obj,
#                                       variable_name = "Index_ctl",
#                                       n_samples = 100,
#                                       sample_fixed = FALSE) ## Only sample random effects
#
# saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
# saveRDS(index_ctl_sd_array, file = paste0(working_dir, "/index_ctl_sd_array.rds"))



# plan(sequential)