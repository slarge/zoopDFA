########
#zoop plots
library(TMB)
library(VAST)
library(FishStatsUtils)
library(dplyr)
library(tidyverse)
library(ggrepel)

dir = "~/Desktop/indicators/indicators_working_output/zoop/GOM-zoop-sal-sfc-common species-spat_var"
dir2<-"~/Desktop/indicators/indicators_working_output/zoop_output"

spp_list <- c("calfin", "ctyp", "cham", "tlong", "oithspp")


#plot of GOM time series, spatial predictions, slopewater comparison

load(file.path(dir,"/output.RData"))

results = plot( fit, category_names=spp_list,
                check_residuals=FALSE,
                plot_set=3,
                working_dir=moddir)

#run extract_factor_things.R

loadings<-read.table(file.path(dir2, "/eps_t_loadings.txt"), header=T)%>%
  dplyr::rename("Factor 1" = X1, "Factor 2" = X2)%>%
  dplyr::mutate(Year = 2000:2017)%>%
  tidyr::pivot_longer(cols = c("Factor 1","Factor 2"),
                      names_to = "Factor",
                      values_to = "Loading")
spat<-readRDS(file.path(dir2, "GOM-zoop-sal-sfc-common species-spat_var/spat_eps_t.txt"))
ts_length <- length(fit$year_labels)


SEs<-read.table(file.path(dir2, "/eps_t_SE.txt"), header=T)%>%
  dplyr::rename("Factor 1" = X1, "Factor 2" = X2)%>%
  dplyr::mutate(Year = 2000:2017)%>%
  tidyr::pivot_longer(cols = c("Factor 1","Factor 2"),
                      names_to = "Factor",
                      values_to = "SE")


#factor loadings plot
graphlabs<-data.frame(Factor = c("Factor 1", "Factor 2"),
                 PropExp= c("Proportion of Variance: 70.6%",
                                      "Proportion of Variance: 29.4%"))

loadings$SE<-SEs$SE

r<-ggplot(loadings,aes(x = Year, y = Loading ))+
  geom_bar(stat="identity", width = 0.5)+
  geom_errorbar(aes(ymin=Loading-SE, ymax=Loading+SE), width=.2,
                position=position_dodge(.9))+
  facet_wrap(~Factor)+
  theme_bw()+
  theme(strip.background=element_blank(),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
r

ggsave(filename="zoop_fact_load.tiff", plot=r, device='tiff', dpi=700)




crs <-  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ne_countries <- rnaturalearth::ne_countries(scale = 10,
                                            continent = "North America",
                                            returnclass = "sf") %>%
  sf::st_transform(crs = crs)

ne_states <- rnaturalearth::ne_states(country = "united states of america",
                                      returnclass = "sf")%>%
  sf::st_transform(crs = crs)

ne_strata <- sf::read_sf("~/Desktop/indicators/shapefiles/BTS_Strata.shp") %>%
  sf::st_transform(crs = crs)

xmin<- min(fit$spatial_list$latlon_g[,2])
xmax<- max(fit$spatial_list$latlon_g[,2])
ymin<- 41
ymax<- max(fit$spatial_list$latlon_g[,1])

xlims<-c(xmin,xmax)
ylims<-c(ymin,ymax)



epsilon1_gct_1<- data.frame(spat$Epsilon1_sct[,,1],
                                     fit$spatial_list$latlon_g)%>%
  dplyr::rename_with(~ spp_list, .cols = !c("Lat", "Lon")) %>%
  tidyr::pivot_longer(cols =  !c("Lat", "Lon"),
                      names_to = "spp",
                      values_to = "val")%>%
  sf::st_as_sf(coords = c("Lon","Lat")) %>%
  sf::st_set_crs(crs)

epsilon1_gct_2<- data.frame(spat$Epsilon1_sct[,,2],
                            fit$spatial_list$latlon_g)%>%
  dplyr::rename_with(~ spp_list, .cols = !c("Lat", "Lon")) %>%
  tidyr::pivot_longer(cols =  !c("Lat", "Lon"),
                      names_to = "spp",
                      values_to = "val")%>%
  sf::st_as_sf(coords = c("Lon","Lat")) %>%
  sf::st_set_crs(crs)

spp_names<- as_labeller(c(
  "calfin"="cfin",
  "cham"="cham",
  "ctyp"="ctyp",
  "oithspp"="oith",
  "tlong"="tlong"))



epsilon1 <-  ggplot() +
  geom_sf(data = epsilon1_gct_1 , aes(color = val)) +
  geom_sf(data = ne_countries, color = "grey40", size = 0.05) +
  geom_sf(data = ne_states, color = "grey40", size = 0.05) +
  facet_wrap(~spp, labeller=spp_names) +
  theme_bw() +
  coord_sf(crs = crs, xlim = xlims, ylim = ylims) +
  theme(legend.position = c(0.8, 0.05),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.title = element_text(size=7), legend.text=element_text(size=7),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle=90, hjust=1),
        strip.background = element_blank())+
  scale_color_viridis_c(option = "cividis") +
  labs( title = "Spatial map of Factor 1",
        x = "",
        y = "",
        color = "log(density)")

epsilon2 <-  ggplot() +
  geom_sf(data = epsilon1_gct_2 , aes(color = val)) +
  geom_sf(data = ne_countries, color = "grey40", size = 0.05) +
  geom_sf(data = ne_states, color = "grey40", size = 0.05) +
  facet_wrap(~spp, labeller=spp_names) +
  theme_bw() +
  coord_sf(crs = crs, xlim = xlims, ylim = ylims) +
  theme(legend.position = c(0.8, 0.05),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.title = element_text(size=7), legend.text=element_text(size=7),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle=90, hjust=1),
        strip.background = element_blank())+
  scale_color_viridis_c(option = "cividis") +
  labs( title = "Spatial map of Factor 2",
        x = "",
       y = "",
       color = "log(density)")

epsilon1
epsilon2

#quartz()
sp<-gridExtra::grid.arrange(epsilon1, epsilon2, nrow=1)

ggsave(filename="zoop_spat_map.tiff", plot=sp, height=9, width=11.5, device='tiff', dpi=700)






