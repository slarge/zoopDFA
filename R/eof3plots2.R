# VAST EOF model with seasonal effects

# Install packages ----
# modified from https://github.com/James-Thorson-NOAA/VAST/wiki/

# install.packages('TMB', type = 'source')
# remotes::install_github("james-thorson/VAST")
library(dplyr)
library(tidyr)
library(ggplot2)
library(VAST)
library(dendextend)
library(patchwork)


# Load data ----

# working_dir <- here::here("analysis/20250410_vast_seasonal_EOF3_310/")
working_dir <- here::here("analysis/20250417_vast_seasonal_EOF3_dual_310/")

fit <- readRDS(here::here(working_dir, "/fit.rds"))
results <- readRDS(here::here(working_dir, "/results.rds"))

load(here::here("data/ecomon_epu.rda"))
spp_list <- c("calfin", "chaeto", "cham", "clauso", "ctyp",
              "euph", "gas", "hyper", "larvaceans",
              "mlucens", "oithspp", "para", "pseudo", "tlong")
start_year <- 1992
end_year <- 2023
zoop_dat <- ecomon_epu %>%
  dplyr::filter(spp %in% spp_list,
                as.numeric(year) >= start_year,
                as.numeric(year) <= end_year,
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

# results = plot_results( fit,
#                         check_residuals = FALSE,
#                         plot_set= c(3,16), #c(3,14,16,18), #c(3,16,18),
#                         working_dir = working_dir,
#                         # Version = "VAST_v14_0_1",
#                         category_names = spp_list,
#                         # year_labels = levels(zoop_dat$year_season),
#                         year_labels = as.numeric(1:length(fit$year_labels)),
#                         strata_names =  c("Georges Bank", "Gulf of Maine", "Mid-Atlantic Bight"))
# saveRDS(results, file = paste0(working_dir, "/results.rds"))
# fit$settings

MapDetails_List <- make_map_info( "Region"= fit$settings$Region,
                                  "spatial_list" = fit$spatial_list,
                                  "Extrapolation_List" = fit$extrapolation_list)


# factor_fit <- plot_factors(fit = fit, Report = fit$Report, ParHat = fit$ParHat, Data = fit$data_list,
#                            SD = fit$parameter_estimates$SD,
#                            category_names = spp_list,
#                            year_labels = as.numeric(1:length(fit$year_labels)),
#                            mapdetails_list = MapDetails_List,
#                            # n_cells = 2000,
#                            plotdir = paste0(working_dir, "/"))

factor_fit <- results$Factors
# Plots ----
## Cluster spp ---- (Only for dual ordination)
hc_c = hclust(dist(fit$Report[["L_epsilon1_cf"]], diag=TRUE, upper=TRUE))
hc_c$labels <- spp_list
hc_c_df <- data.frame(group = cutree(hc_c, k = 2))
hc_c_df$spp <- spp_list

## Hierarchical clustering of Epsilon Time
hc_t <- hclust(dist(fit$Report[["Ltime_epsilon1_tf"]], diag=TRUE, upper=TRUE))
hc_t_df <- data.frame(group = cutree(hc_t, k = 2))
hc_t_df$year_season <- gsub(" ", "_", rownames(hc_t_df))
hc_t_df$year_season <- gsub(" ", "_", levels(zoop_dat$year_season))
hc_t$labels <- gsub(pattern = "_", " ", levels(zoop_dat$year_season))


## Index plot ----
index_years <- zoop_dat %>%
  select(year_season) %>%
  unique() %>%
  mutate(Time = 1:n() -1)


yearseason_set <- expand.grid(year_season = paste0(rep(start_year:end_year, each = 4), "_",
                                                   c("winter", "spring", "summer", "fall")),
                              epu = c("Georges Bank", "Gulf of Maine", "Mid-Atlantic Bight"),
                              spp = spp_list) %>%
  arrange(year_season) %>%
  left_join(index_years)


index <- read.csv(file = here::here(working_dir, "/index.csv")) %>%
  rename(#year_season = Time,
    epu = Stratum,
    value = Estimate,
    value_sd = Std..Error.for.Estimate,
    value_lnsd = Std..Error.for.ln.Estimate.,
    spp = Category) %>%
  left_join(yearseason_set) %>%
  mutate( year_season = factor(year_season, levels = levels(zoop_dat$year_season)),
          spp = as.factor(spp))
# str(index)
# index$year_season <- as.factor(levels(zoop_dat$year_season))
# index$spp <- as.factor(spp_list)
# index$epu <- c("Georges Bank", "Gulf of Maine", "Mid-Atlantic Bight")

index_full <- index %>%
  left_join(yearseason_set) %>%
  mutate(year = as.numeric(gsub("_.*", "", year_season)),
         season = factor(gsub(".*_", "", year_season),
                         levels = c("winter", "spring", "summer", "fall")),
         time = case_when(season == "winter" ~ 0,
                          season == "spring" ~ 0.25,
                          season == "summer" ~ 0.5,
                          season == "fall" ~ 0.75,
                          TRUE ~ NA_integer_),
         time = year + time,
         value_up = value + value_sd,
         value_lo = value - value_sd,
         value_lo = ifelse(value_lo < 0, 0, value_lo)) #%>%
  # left_join(hc_c_df)

spp_names <- c("calfin" = "italic(Calanus~~finmarchicus)",
               "cham" = "italic(Centropages~~hamatus)",
               "chaeto" = "Chaetognatha",
               "clauso" = "italic(Clausocalanus~~arcuicornis)",
               "ctyp" = "italic(Centropages~~typicus)",
               "euph" = "Euphausiacea",
               "gas"= "Gastropoda",
               "hyper" = "Hyperiidea",
               "larvaceans" = "larvaceans",
               "mlucens" = "italic(Metridia~~lucens)",
               "oithspp" = "italic(Oithona)~~spp.",
               "para" = "italic(Paracalanus~~parvus)",
               "pseudo" = "italic(Pseudocalanus)~~spp.",
               "tlong" = "italic(Temora~~longicornis)")
season_names <- c("winter" = "winter",
                  "spring" = "spring",
                  "summer" = "summer",
                  "fall" = "fall")



spp_labels <- as_labeller(spp_names, label_parsed)

epu_colors <- c("Mid-Atlantic Bight" = "#5145b9",
                "Gulf of Maine" =  "#FF8400",
                "Georges Bank" = "#007582")



index_plot <- ggplot(index_full, aes(x = time, y = value/1000000000, ymin = value_lo/1000000000, ymax = value_up/1000000000,
                                color = epu, fill = epu)) +
  geom_point(na.rm = FALSE, alpha = 0.25, size = 0.5, shape = 19) +
  geom_errorbar(na.rm = FALSE, alpha = 0.5, linewidth = 0.5) +
  facet_wrap(.~spp, scales = "free_y", labeller = spp_labels, ncol = 3) +
  labs(x = "", y = expression("billion/km"^2), color = "", fill = "") +
  theme_minimal() +
  scale_color_manual(values = epu_colors) +
  theme(legend.position = "bottom",
        text = element_text(size = 8)) +
  NULL


ggsave(filename = here::here(working_dir, "index_plot.jpeg"), index_plot, bg = "white",
       width = 170, units = "mm", device = jpeg, dpi = 500)

### Cluster Loadings plot ----


## CCA
ndim_cca <- 2 # number of dimension to filter in the EOFs
ndim_cat <- 4
# names(fit$Report)
E = list(
  "v" = fit$Report$Ltime_epsilon1_tf, # temporal loadings (p * p)
  "u" = fit$Report$Epsiloninput1_gff[, 1:ndim_cat,] # spatial loadings (n * p)
)

## Construct the seasonal variable
t_series <- expand_grid(year = start_year:end_year,
                        season = c("winter", "spring", "summer", "fall")) %>%
  mutate(season = factor(season,  levels = c("winter", "spring", "summer", "fall")),
         year_season = factor(paste0(year, "_", season))) %>%
  left_join(data.frame(year_season = levels(zoop_dat$year_season),
                       Time = seq(1, nrow(E$v))))

p <- nrow(t_series)
ts <- seq(1, p, 1) # sequence of time steps
component.strength <- 2 # amplitude of the signal
component.freqs <- 1/4 # frequency of signal (quarter^-1)
f.0 <- 2 * pi # fundamental frequency (month^-1)
component.delay <- 0 # delay of signal components (radians)

# Create the seasonal signal
signal_vec <- - component.strength *
  sin(component.freqs * f.0 * ts + component.delay)

t_sig <- t_series %>%
  mutate(signal_vec = ifelse(is.na(Time), NA, signal_vec),
         time = case_when(season == "winter" ~ 0,
                          season == "spring" ~ 0.25,
                          season == "summer" ~ 0.5,
                          season == "fall" ~ 0.75,
                          TRUE ~ NA_integer_),
         time = year + time) %>%
  drop_na()


## Filter EOFs dimensions and apply CCA
EOFset1 <- E$v[, 1:ndim_cca] * sqrt(p - 1)# filter the loadings
cca_res <- stats::cancor(EOFset1, t_sig$signal_vec) # perform CCA

## Rotate the temporal loadings with CCA results
ccavar1_t <- data.frame(Time = 1:length(fit$year_labels),
                        ccavar1 = scale((as.matrix(EOFset1) %*% cca_res$xcoef)[,1], scale = TRUE, center = TRUE)[,1],
                        eofvar1 =  EOFset1[,1],
                        eofvar2 =  EOFset1[,2])

cca_ts <- t_series %>%
  mutate(Ancillary = signal_vec) %>%
  left_join(ccavar1_t) %>%
  mutate(dummy_var = scale(signal_vec, scale = TRUE, center = TRUE)[,1],
         Ancillary = scale(signal_vec, scale = TRUE, center = TRUE)[,1],
         eofvar1 = scale(eofvar1, scale = TRUE, center = TRUE)[,1],
         eofvar2 = scale(eofvar2, scale = TRUE, center = TRUE)[,1]) %>%
  tidyr::pivot_longer(Ancillary:dummy_var) %>%
  mutate(type = ifelse(name %in% c("ccavar1", "dummy_var"), "CCA variable", "EOF variable"),
         name = case_when(name == "ccavar1" ~ "CCA",
                          name %in% c("Ancillary", "dummy_var") ~ "Ancillary",
                          name == "eofvar1" ~ "EOF Dimension 1",
                          name == "eofvar2" ~ "EOF Dimension 2",
                          TRUE ~ NA_character_),
         type = factor(type, levels = c("EOF variable","CCA variable")),
         name = factor(name, levels = c("EOF Dimension 1","EOF Dimension 2","Ancillary","CCA")),
         time = case_when(season == "winter" ~ 0,
                          season == "spring" ~ 0.25,
                          season == "summer" ~ 0.5,
                          season == "fall" ~ 0.75,
                          TRUE ~ NA_integer_),
         time = year + time)

cca_ts_plot <- ggplot(data = cca_ts, aes(x = time, y = value, color = name)) +
  geom_vline(xintercept = cca_ts$time[which(stringr::str_detect(cca_ts$season,"winter"))],
             # linetype = "dashed",
             color = "grey90", alpha = .75) +
  geom_point(size = 0.5, shape = 19) +
  geom_path() +
  facet_wrap(.~type, nrow = 2) +
  labs(y = "", x = "") +
  scale_color_manual(values = c("EOF Dimension 1" = "#3e86ffff",
                                "EOF Dimension 2" = "#66C2A5",
                                "Ancillary" = "#F8766DFF",
                                "CCA" = "#00ba38ff")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  NULL
# cca_ts_plot

cor_mat <- cca_ts %>%
  filter(!(name == "Ancillary" & type == "EOF variable")) %>%
  select(time, name, value) %>%
  mutate(name = factor(name, levels = c("Ancillary", "CCA", "EOF Dimension 1", "EOF Dimension 2"))) %>%
  arrange(time, name) %>%
  pivot_wider(id_cols = time, names_from = name) %>%
  select(-time) %>%
  cor(use = "complete.obs")

cor_mat_plot <- ggcorrplot::ggcorrplot((cor_mat),
                                       lab = TRUE,
                                       colors = c("blue", "white", "red"))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_blank())
# cor_mat_plot

cor_p <- (plot_spacer() + cor_mat_plot) / cca_ts_plot + plot_layout(nrow = 2, widths = c(1, 1, 2), heights = c(2, 2, 4))

ggsave(filename = here::here(working_dir, "correlation_plot.jpeg"), cor_p, bg = "white",
       width = 170, height = 191, units = "mm", device = jpeg, dpi = 500)


### PCA

pca_ts <- rotate_factors(L_pj = fit$Report$Ltime_epsilon1_tf, RotationMethod = "PCA", testcutoff = 1e-10)
# compare_pca <- data.frame(ffit$Rotated_loadings$EpsilonTime1,
#                           tt$L_pj_rot)




factor_t <- t_series %>%
  left_join(data.frame(pca_ts$L_pj_rot, #results$Factors$Rotated_loadings$EpsilonTime1,
                       year_season = levels(zoop_dat$year_season)), by = join_by(year_season)) %>%
  left_join(cca_ts %>% filter(type == "CCA variable", name == "CCA") %>% select(year_season, cca = value), by = join_by(year_season)) %>%
  left_join(cca_ts %>% filter(type == "CCA variable", name == "Ancillary") %>% select(year_season, ancillary = value), by = join_by(year_season)) %>%
  left_join(hc_t_df, by = join_by(year_season))


# factor_t <- data.frame(pca_ts$L_pj_rot, #results$Factors$Rotated_loadings$EpsilonTime1,
#                        year_season = levels(zoop_dat$year_season)) %>%
#   right_join(yearseason_set, by = join_by(year_season)) %>%
#
#

factor_time_all <- factor_t %>%
  mutate(group = 0)

### The problem is ancillary needs to be a full timeseries with no NAs and join to the clusters that might have NAs. Ancillary should have no missing values

factor_time <- bind_rows(factor_t,
                         factor_time_all) %>%
  mutate(year = as.numeric(gsub("_.*", "", year_season)),
         season = factor(gsub(".*_", "", year_season),
                         levels = c("winter", "spring", "summer", "fall")),
         time = case_when(season == "winter" ~ 0,
                          season == "spring" ~ 0.25,
                          season == "summer" ~ 0.5,
                          season == "fall" ~ 0.75,
                          TRUE ~ NA_integer_),
         time = year + time,
         group = case_when(group == 0 ~ "All data",
                           group > 0 ~ paste0("Cluster ", group),
                           TRUE ~ NA_character_)) %>%
  rename(`Factor 1` = X1,
         `Factor 2` = X2) %>%
  pivot_longer(cols = c(`Factor 1`, `Factor 2`,
                        "cca", "ancillary"))

labels_pca = factor_time %>%
  drop_na(group) %>%
  filter(name %in% c("Factor 1", "Factor 2")) %>%
  select(name, group) %>%
  distinct %>%
  group_by(group, name) %>%
  mutate(label = letters[cur_group_id()])

labels_cca = factor_time %>%
  drop_na(group) %>%
  filter(name == "cca") %>%
  select(name, group) %>%
  distinct %>%
  group_by(group, name) %>%
  mutate(label = letters[cur_group_id()])


factor_labels <- as_labeller(c(`factor_1` = "Factor 1",
                               `factor_2` = "Factor 2",
                               `hclust_1` = "Cluster 1",
                               `hclust_2` = "Cluster 2"))

factor_cluster_pca_plot <- ggplot(factor_time %>% filter(name %in% c("Factor 1", "Factor 2")) %>% drop_na(group), aes(x = time, y = value, color = season)) +
  geom_segment(aes(xend = time, yend = 0), na.rm = FALSE, linewidth = 0.5) +
  geom_point(na.rm = FALSE, size = 0.5) +
  geom_text(data = labels_pca, aes(label = label), x = start_year + 3, y = 2.5, hjust = 0, vjust = 0, inherit.aes = FALSE) +
  geom_hline(yintercept = 0) +
  # geom_line() +
  labs(x = "", y = "Factor Loadings", color = "Season") +
  scale_color_brewer(palette = "Set1") +
  facet_grid(group~name, drop = TRUE) +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 8)) +
  NULL
factor_cluster_pca_plot

ancillary_time <- factor_time %>%
  filter(name == "ancillary") %>%
  select(time, value, season) %>%
  expand_grid(group = c("All data", "Cluster 1", "Cluster 2"))

cca_time <- factor_time %>%
  filter(name == "cca") %>%
  drop_na(group)

factor_cluster_cca_plot <- ggplot() +
  # geom_point(data = ancillary_time, aes(x = time, y = value), color = "#F8766DFF", alpha = 0.5) +
  geom_line(data = ancillary_time, aes(x = time, y = value), color = "#F8766DFF", alpha = 0.5) +
  geom_segment(data = cca_time, aes(x = time, y = value, xend = time, yend = 0, color = season), na.rm = FALSE, linewidth = 1) +
  geom_point(data = cca_time, aes(x = time, y = value, color = season), na.rm = FALSE, size = 1) +
  geom_text(data = labels_cca, aes(label = label), x = start_year, y = 1.25, hjust = 0, vjust = 0, inherit.aes = FALSE) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Factor Loadings", color = "Season") +
  scale_color_brewer(palette = "Set1") +
  facet_grid(group~name) +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 8)) +
  NULL
factor_cluster_cca_plot


dend_year <- as.dendrogram(hc_t, type = "rectangle") %>%
  color_labels(k = 2, col = c("black", "grey50")) %>%
  color_branches(k = 2, col = c("black", "grey50")) %>%
  set("labels_cex", .4) %>%
  ladderize()

tree_plot_year <- function(main = "g") {
  par(mar = c(0,1,0,0))
  tree_plot <- wrap_elements(full = ~plot_horiz.dendrogram(dend_year, side = TRUE,
                                                           center = FALSE, main = main,
                                                           cex.main = .75, font.main = 1, yaxt = "n", text_pos = 4, yaxs = "i"))
  return(tree_plot)
}

cluster_loadings_pca <- factor_cluster_pca_plot + tree_plot_year(main = "g")
cluster_loadings_cca <- factor_cluster_cca_plot + tree_plot_year(main = "d")


ggsave(filename = here::here(working_dir, "/cluster_loadings_pca.jpeg"), cluster_loadings_pca, bg = "white",
       width = 170, height = 191, units = "mm", device = jpeg, dpi = 500)
ggsave(filename = here::here(working_dir, "/cluster_loadings_cca.jpeg"), cluster_loadings_cca, bg = "white",
       width = 170, height = 191, units = "mm", device = jpeg, dpi = 500)



### ## Cluster spp ---- (Only for dual ordination)
dend_spp <- as.dendrogram(hc_c, type = "rectangle") %>%
  color_labels(k = 2, col = c("black", "grey50")) %>%
  color_branches(k = 2, col = c("black", "grey50")) %>%
  set("labels_cex", .8) %>%
  ladderize()

par(mar = c(0,1,0,0))
tree_plot_spp <- wrap_elements(full = ~plot_horiz.dendrogram(dend_spp, side = TRUE,
                                                             center = FALSE,
                                                             cex.main = .85, font.main = 1, yaxt = "n", text_pos = 4, yaxs = "i"))
# cluster_loadings <- factor_cluster_plot + tree_plot_year

ggsave(filename = here::here(working_dir, "/cluster_loadings_spp.jpeg"), tree_plot_spp, bg = "white",
       width = 170, units = "mm", device = jpeg, dpi = 500)


## (Only for dual ordination)
epsilon1_loadings_se <- data.frame(factor_fit$Rotated_loadings_SE$Epsilon1) %>% #results$Factors$Rotated_loadings_SE$Epsilon1) %>%
  tibble::rownames_to_column(var = "spp") %>%
  pivot_longer(!spp, names_to = "factor", values_to = "se")
#
#
epsilon1_loadings <- data.frame(results$Factors$Rotated_loadings$Epsilon1) %>%
  tibble::rownames_to_column(var = "spp") %>%
  pivot_longer(!spp, names_to = "factor", values_to = "loadings") %>%
  left_join(epsilon1_loadings_se) %>%
  arrange(desc(spp)) %>%
  mutate(spp = factor(spp, levels = unique(spp)),
         upper = loadings + se,
         lower = loadings - se,
         factor = gsub("X", "Factor ", factor),
         season = case_when(grepl(1, factor) ~ "winter",
                            grepl(2, factor) ~ "spring",
                            grepl(3, factor) ~ "summer",
                            grepl(4, factor) ~ "fall",
                            TRUE ~ NA_character_),
         season = factor(season, levels = c("winter", "spring", "summer", "fall")))
#
factor_taxa <- ggplot(data = epsilon1_loadings, aes(x = spp, y = loadings, ymin = lower, ymax = upper, color = season)) +
  geom_hline(yintercept = 0) +
  geom_linerange(show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  scale_x_discrete(labels = \(x) scales::label_parse()(spp_names[x])) +
  facet_wrap(. ~ season, nrow = 1) +
  labs(x = "", y = "factor loadings") +
  theme_minimal()
#
ggsave(filename = here::here(working_dir, "/factor_loadings_spp.jpeg"), factor_taxa, bg = "white",
       width = 170, units = "mm", device = jpeg, dpi = 500)

### Temporal variation (β) for each variable, where these intercepts are specified to follow a random-walk.
### Spatial variation (ω) in the expected value for each variable, representing long-term spatial patterns.
### Spatio-temporal variation (ε), estimated as one or more dominant modes of ecosystem variability as well as
### a map representing the spatial response for each variable to these estimated modes of variation.


## EpsilonTime 1

## Download data layers
## 1) North America layer
xlims = c(-77, -65)
ylims = c(35, 45)
# crs <- sf::st_crs("+proj=utm +zone=19 +datum=WGS84 +units=km")
# crs <-  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs <- 4269
trans_crs <- 32619
ne_countries <- rnaturalearth::ne_countries(scale = 10,
                                            continent = "North America",
                                            returnclass = "sf") %>%
  sf::st_transform(crs = trans_crs)

## 2) State layer
ne_states <- rnaturalearth::ne_states(country = "united states of america",
                                      returnclass = "sf") %>%
  sf::st_transform(crs = trans_crs)


# cca

cca_epsilontime_1 <- data.frame(season = 1:4) %>%
  mutate(epsilontime = purrr::map(.x = season, ~ (as.matrix(E$u[,.x,]) %*% cca_res$xcoef[,1])[,1]),
         season = case_when(season == 1 ~ "winter",
                            season == 2 ~ "spring",
                            season == 3 ~ "summer",
                            season == 4 ~ "fall",
                            TRUE ~ NA_character_)) %>%
  tidyr::unnest(cols = epsilontime) %>%
  mutate(Lat = rep.int(fit$spatial_list$latlon_g[,"Lat"], times = ndim_cat),
         Lon = rep.int(fit$spatial_list$latlon_g[, "Lon"], times = ndim_cat),
         season = factor(season, levels = c("winter", "spring", "summer", "fall"))) %>%
  sf::st_as_sf(coords = c("Lon","Lat")) %>%
  sf::st_set_crs(crs) %>%
  sf::st_transform(crs = trans_crs) ## project to UTM Zone 19N

# pca
dimnames(factor_fit$Rotated_projected_factors$EpsilonTime1) <- list(1:2000, season_names, c("Factor_1", "Factor_2"))
pca_epsilontime_1 <- data.frame(factor_fit$Rotated_projected_factors$EpsilonTime1,# Factors$Rotated_projected_factors$EpsilonTime1,
                                fit$spatial_list$latlon_g) %>%
  tidyr::pivot_longer(cols =  !c("Lat", "Lon"),
                      names_to = "season_factor",
                      values_to = "epsilontime") %>%
  tidyr::separate_wider_delim(delim = ".", names = c("season", "factor"), cols = "season_factor") %>%
  mutate(epsilontime = as.numeric(epsilontime),
         season = factor(season, levels = c("winter", "spring", "summer", "fall")),
         factor = as.factor(gsub("_", " ", factor))) %>%
  sf::st_as_sf(coords = c("Lon","Lat")) %>%
  sf::st_set_crs(crs) %>%
  sf::st_transform(crs = trans_crs) ## project to UTM Zone 19N


region <- NEFSCspatial::epu_sf %>%
  sf::st_transform(trans_crs) %>%
  sf::st_make_valid() %>%
  filter(EPU != "SS")

# Get the coordinates from the sf object
coords <- sf::st_coordinates(cca_epsilontime_1)

# Compute the tessellation using the deldir package
voronoi <- deldir::deldir(coords[,1], coords[,2])

# Convert the Voronoi tiles to sf polygons
voronoi_polygons_cca <- sf::st_collection_extract(
  sf::st_voronoi(sf::st_union(cca_epsilontime_1))
)
voronoi_polygons_pca <- sf::st_collection_extract(
  sf::st_voronoi(sf::st_union(pca_epsilontime_1))
)

# Create a boundary (the convex hull of the points) to clip the polygons
boundary <- sf::st_union(region)

# Create an sf data frame from the polygons
voronoi_sf_cca <- sf::st_sf(geometry = voronoi_polygons_cca)
voronoi_sf_pca <- sf::st_sf(geometry = voronoi_polygons_pca)

# Clip the Voronoi polygons to the boundary
voronoi_clipped_cca <- sf::st_intersection(voronoi_sf_cca, boundary)
voronoi_clipped_pca <- sf::st_intersection(voronoi_sf_pca, boundary)

# Join the original point data to the new polygons
final_polygons_cca <- sf::st_join(voronoi_clipped_cca, cca_epsilontime_1)
final_polygons_pca <- sf::st_join(voronoi_clipped_pca, pca_epsilontime_1)

# pca
epsilontime_plot_pca <- ggplot() +
  geom_sf(data = final_polygons_pca %>% drop_na(season), aes(color = epsilontime, fill = epsilontime)) +
  geom_sf(data = ne_countries, color = "grey40", size = 0.05) +
  geom_sf(data = ne_states, color = "grey40", size = 0.05) +
  facet_grid(season~factor) +
  coord_sf(crs = crs, xlim = xlims, ylim = ylims) +
  scale_color_viridis_c(option = "inferno", na.value = "transparent") +
  scale_fill_viridis_c(option = "inferno", na.value = "transparent") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "bottom",
        # legend.position.inside = c(1, 0),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0)) +
  NULL
# epsilontime_plot

ggsave(filename = here::here(working_dir, "/epsilontime_season_pca.jpeg"), epsilontime_plot, bg = "white",
       width = 85, units = "mm", device = jpeg, dpi = 500)

epsilontime_cca_plot <- ggplot() +
  geom_sf(data = final_polygons_cca %>% drop_na(season), aes(fill = epsilontime, color = epsilontime)) +#, size = 0.5, stroke=0, shape=16) +
  facet_wrap(~season, drop = TRUE, ncol = 1) +
  geom_sf(data = ne_countries, color = "grey40", size = 0.05) +
  geom_sf(data = ne_states, color = "grey40", size = 0.05) +
  coord_sf(crs = crs, xlim = xlims, ylim = ylims) +
  scale_color_viridis_c(option = "rocket", na.value = "transparent") +
  scale_fill_viridis_c(option = "rocket", na.value = "transparent") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "bottom",
        # legend.position.inside = c(1, 0),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0)) +
  NULL
# epsilontime_cca_plot

ggsave(filename = here::here(working_dir, "/epsilontime_season_cca.jpeg"), epsilontime_cca_plot, bg = "white",
       width = 85, units = "mm", device = jpeg, dpi = 500)


cor_p <- (epsilontime_cca_plot + cor_mat_plot) / cca_ts_plot + plot_layout(nrow = 2, widths = c(1, 1, 2), heights = c(2, 2, 4))

ggsave(filename = here::here(working_dir, "correlation_plot.jpeg"), cor_p, bg = "white",
       width = 170, height = 191, units = "mm", device = jpeg, dpi = 500)



### CCA Epsiloninput1_gff

c_index = 1
E = list(
  "v" = fit$Report$Ltime_epsilon1_tf,
  "u" = fit$Report$Epsiloninput1_gff[,c_index,]
)









str(fit$Report$Epsiloninput1_gff)
dimnames(fit$Report$Epsiloninput1_gff) <- list(1:2000, season_names, c("Factor_1", "Factor_2"))

epsiloninput_1 <- data.frame(season = 1:4) %>%
  mutate(dat = purrr::map(.x = season, ~ (fit$Report$Epsiloninput1_gff[, .x, ] %*% cca_res$xcoef[, 1])[,1])) %>%
  tidyr::unnest(cols = c(dat)) %>%
  mutate(Lat = rep(fit$spatial_list$latlon_g[,1], 4),
         Lon = rep(fit$spatial_list$latlon_g[,2], 4),
         season = case_when(season == 1 ~ "winter",
                            season == 2 ~ "spring",
                            season == 3 ~ "summer",
                            season == 4 ~ "fall",
                            TRUE ~ NA_character_),
         season =  factor(season, levels = c("winter", "spring", "summer", "fall"))) %>%
  rename("epsilontime" = 2) %>%
  sf::st_as_sf(coords = c("Lon","Lat")) %>%
  sf::st_set_crs(crs)



# epsiloninput_1 <- data.frame(fit$Report$Epsiloninput1_gff,
#                             fit$spatial_list$latlon_g) %>%
#   tidyr::pivot_longer(cols =  !c("Lat", "Lon"),
#                       names_to = "season_factor",
#                       values_to = "epsilontime") %>%
#   tidyr::separate_wider_delim(delim = ".", names = c("season", "factor"), cols = "season_factor") %>%
#   mutate(epsilontime = as.numeric(epsilontime),
#          season = factor(season, levels = c("winter", "spring", "summer", "fall")),
#          factor = as.factor(gsub("_", " ", factor))) %>%
#   sf::st_as_sf(coords = c("Lon","Lat")) %>%
#   sf::st_set_crs(crs)


epsilontime_plot <- ggplot() +
  geom_sf(data = epsiloninput_1, aes(color = epsilontime), alpha = 0.5) +
  geom_sf(data = ne_countries, color = "grey40", size = 0.05) +
  geom_sf(data = ne_states, color = "grey40", size = 0.05) +
  facet_wrap(season~.) +
  coord_sf(crs = crs, xlim = xlims, ylim = ylims) +
  scale_color_viridis_c(option = "inferno", na.value = "transparent") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text = element_text(size = rel(.75))) +
  # theme(legend.position = "bottom") +
  NULL
epsilontime_plot

ggsave(filename = here::here("analysis/vast_seasonal_EOF3_dual/epsilontime_season_cca.jpeg"), epsilontime_plot, bg = "white",
       width = 85, units = "mm", device = jpeg, dpi = 500)





## Start the GAMS ----
# yearseason_set <- expand.grid(year_season = paste0(rep(1992:2021, each = 4), "_",
#                                                    c("winter", "spring", "summer", "fall")))
#
# # results$Factors$Rotated_loadings$EpsilonTime1
# factor_time <- data.frame(results$Factors$Rotated_loadings$EpsilonTime1,
#                           results$Factors$Rotated_loadings_SE$EpsilonTime1,#fit$Report$Ltime_epsilon1_tf,
#                           year_season = levels(zoop_dat$year_season)) %>%
#   right_join(yearseason_set) %>%
#   mutate(year = as.numeric(gsub("_.*", "", year_season)),
#          season = factor(gsub(".*_", "", year_season),
#                          levels = c("winter", "spring", "summer", "fall")),
#          nseason = case_when(season == "winter" ~ 0,
#                              season == "spring" ~ 0.25,
#                              season == "summer" ~ 0.50,
#                              season == "fall" ~ 0.75),
#          time = year + nseason) %>%
#   rename(`Factor_1` = X1,
#          `Factor_2` = X2,
#          `Factor_1SE` = X1.1,
#          `Factor_2SE` = X2.1)

# ## Timeseries analysis
# factor_dat <- factor_time %>%
#   select(year, season, value, name) %>%
#   distinct() %>%
#   group_by(name) %>%
#   complete(year, season = season_names, fill = list(value = 0),
#            explicit = TRUE) %>%
#   mutate(date = case_when(season == "winter" ~ paste0(year, "-01-01"),
#                           season == "spring" ~ paste0(year, "-04-01"),
#                           season == "summer" ~ paste0(year, "-07-01"),
#                           season == "fall" ~ paste0(year, "-10-01"),
#                           TRUE ~ NA_character_),
#          date = as.Date(date, format = "%Y-%m-%d"),
#          yq = tsibble::yearquarter(date))
#
# tt <- as_tsibble(factor_dat, index = yq, key = name) %>%
#   select(yq, value, name) %>%
#   model(STL(value ~ season(period = "year"), robust = TRUE)) %>%
#   components()
#
# ggplot() +
#   geom_path(data = tt, aes(x = yq, y = trend, color = name))
#
#   autoplot() +
#   theme_minimal()
#


index1_tf <- results$Factors$Rotated_loadings$EpsilonTime1
colSums(index1_tf^2)/sum(index1_tf^2)

rowSums(results$Factors$Rotated_loadings$Epsilon1^2)/sum(results$Factors$Rotated_loadings$Epsilon1^2)
colSums(results$Factors$Rotated_loadings$Epsilon1^2)/sum(results$Factors$Rotated_loadings$Epsilon1^2)


factor_season_plot <- ggplot(factor_time %>% filter(group == "All data"), aes(x = year, y = value, color = season)) +
  geom_segment(aes(xend = year, yend = 0), na.rm = FALSE, linewidth = 0.5) +
  geom_point(na.rm = FALSE, size = 0.5) +
  # geom_text(data = labels, aes(label = label), x = 1992, y = 1.8, hjust = 0, vjust = 0, inherit.aes = FALSE) +
  geom_hline(yintercept = 0) +
  # geom_line() +
  labs(x = "", y = "Factor Loadings", color = "Season") +
  scale_color_brewer(palette = "Set1") +
  facet_grid(.~name) +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 8)) +
  NULL


## ----
## MGCV weight argument with 1/SE^2


# factor_t <- data.frame(fit$Report$Ltime_epsilon1_tf,
#                        year_season = levels(zoop_dat$year_season)) %>%
#   right_join(yearseason_set)

factor_t <- data.frame(results$Factors$Rotated_loadings$EpsilonTime1,
                       # fit$Report$Ltime_epsilon1_tf,
                       year_season = levels(zoop_dat$year_season)) %>%
  right_join(yearseason_set) %>%
  mutate(year = as.numeric(gsub("_.*", "", year_season)),
         season = factor(gsub(".*_", "", year_season),
                         levels = c("winter", "spring", "summer", "fall")),
         time = case_when(season == "winter" ~ 0,
                          season == "spring" ~ 0.25,
                          season == "summer" ~ 0.5,
                          season == "fall" ~ 0.75,
                          TRUE ~ NA_integer_),
         time = year + time) %>%
  rename(`Factor 1` = X1,
         `Factor 2` = X2) %>%
  pivot_longer(cols = c(`Factor 1`, `Factor 2`))

labels = factor_time %>%
  select(name, group) %>%
  distinct %>%
  group_by(group, name) %>%
  mutate(label = letters[cur_group_id()])

plot(fit$Report$Ltime_epsilon1_tf[,1], type = "l")
plot(fit$Report$Ltime_epsilon1_tf[,1]~ results$Factors$Rotated_loadings$EpsilonTime1[,1], type = "p")


factor_labels <- as_labeller(c(`factor_1` = "Factor 1",
                               `factor_2` = "Factor 2"))

factor_loadings <- ggplot(factor_t, aes(x = time, y = value, color = season)) +
  geom_hline(yintercept = 0) +
  geom_segment(aes(xend = time, yend = 0), na.rm = FALSE, size = 0.5) +
  geom_point(na.rm = FALSE, size = 0.5) +
  # geom_text(data = labels, aes(label = label), x = 1992, y = 2, hjust = 0, vjust = 0, inherit.aes = FALSE) +
  # geom_line() +
  labs(x = "", y = "Factor Loadings", color = "Season") +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~name, ncol = 1) +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 10))


ggsave(filename = here::here("analysis/vast_seasonal_EOF3_dual/factor_loadings.jpeg"), factor_loadings, bg = "white",
       width = 170, units = "mm", device = jpeg, dpi = 500)

##

zoop_dat <- ecomon_epu %>%
  dplyr::filter(spp %in% spp_list,
                #                  EPU %in% c("GB", "GOM", "MAB"),
                #                  as.numeric(year) >= 2010,
                #                  as.numeric(year) < 2015) %>%
                as.numeric(year) >= end_year,
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



ffit <- plot_factors(fit = fit, Report = fit$Report, ParHat = fit$ParHat, Data = fit$data_list,
                     SD = fit$parameter_estimates$SD,
                     category_names = spp_list)

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

levels(season_dat$year_season)

season_td <- td %>%
  right_join(season_dat) %>%
  pivot_longer(cols = starts_with("Factor"), names_to = "var", values_to = "val") %>%
  arrange(year_season) %>%
  mutate(year_season = factor(year_season, levels = levels(year_labels)))

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

phase_td <- td %>%
  right_join(season_dat) %>%
  # pivot_longer(cols = starts_with("Factor"), names_to = "var", values_to = "val") %>%
  arrange(year_season) %>%
  mutate(year_season = factor(year_season, levels = levels(year_labels)))

library(ggnewscale)
library(gghighlight)

ggplot(phase_td, aes(x = Factor_1, y = Factor_2)) +
  geom_path(aes(color = as.numeric(year))) +
  scale_color_viridis_c(option = "magma") +
  new_scale_color() +
  geom_point(aes(color = as.factor(season))) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  facet_wrap(~season, scales = "fixed")


# geom_line() +
geom_label(aes(label = year_season, color = season))


# geom_segment( aes(y=0, yend = val)) +
facet_wrap(~var, nrow = 2) +
  scale_color_brewer(palette = "Set1") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_minimal() +
  NULL

tt <- zoop_dat %>%
  right_join(season_dat) %>%
  arrange(year_season) %>%
  mutate(year_season = factor(year_season, levels = levels(season_dat$year_season)))

levels(tt$year_season)

raw_plot <- ggplot(tt, aes(x = lon, y = lat)) +
  geom_point(data = tt %>% filter(abundance == 0), color = "black", fill = "black", shape = 21) +
  geom_point(data = tt %>% filter(abundance > 0), aes(color = spp, size = abundance), alpha = 0.5) +
  facet_wrap(~ year_season, ncol = 4) +
  theme_minimal() +
  # labs(title = i) +
  NULL

ggsave(filename = here::here(working_dir, "figures/raw_plot.png"), plot = raw_plot, bg = "white", width = 1200, height = 2000, units = "px")


results = plot_results( fit,
                        check_residuals = FALSE,
                        plot_set= c(3,14,16,18), #c(3,16,18),
                        working_dir = "analysis/vast_seasonal_EOF3/",
                        category_names = spp_list,
                        Version = "VAST_v14_0_1",
                        year_labels = levels(zoop_dat$year_season),
                        strata_names =  c("Georges Bank", "Gulf of Maine", "Mid-Atlantic Bight"))


results <- readRDS(here::here(working_dir, "results.rds"))


f2_ts = fit$covariate_data[match(fit$year_labels, fit$covariate_data$Year),'AREA_SUM_KM2_LTE2']

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

