### Load packages
# remotes::install_github("james-thorson-NOAA/FishStatsUtils@dev")
library(VAST)

### load data set
example = load_example( data_set="five_species_ordination" )

working_dir <- here::here("analysis/eof3_test/")
### Make settings:
### including modifications from default settings to match
### analysis in original paper
settings = make_settings( n_x=25,
                          Region=example$Region,
                          purpose="EOF3",
                          n_categories=2,
                          mesh_package = "fmesher",
                          ObsModel=c(1,1),
                          RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0) )

### Run model (including settings to speed up run)
fit = fit_model( settings=settings,
                 Lat_i=example$sampling_data[,'Lat'],
                 Lon_i=example$sampling_data[,'Lon'],
                 t_i=example$sampling_data[,'Year'],
                 c_i=example$sampling_data[,'species_number']-1,
                 b_i=example$sampling_data[,'Catch_KG'],
                 a_i=example$sampling_data[,'AreaSwept_km2'],
                 newtonsteps=0,
                 getsd = TRUE,
                 Use_REML = TRUE,
                 run_model = TRUE,
                 working_dir = paste0(working_dir, "/"),
                 optimize_args = list("lower" = -Inf,
                                      "upper" = Inf))
saveRDS(fit, file = paste0(working_dir, "/fit.rds"))

### Plot results, including spatial term Omega1

fit <- readRDS(here::here(working_dir, "fit.rds"))
fit <- reload_model(fit)

results = plot( fit,
                check_residuals=FALSE,
                plot_set=c(3,16),
                category_names = c("pollock", "cod", "arrowtooth", "snow_crab", "yellowfin") )

str(results$Factors)

sessionInfo()

