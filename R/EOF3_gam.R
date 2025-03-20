# remotes::install_github("NOAA-EDAB/ecodata")
# library(ecodata)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(mgcv)
library(gratia)
##

# Load data ----

working_dir <- here::here("analysis/vast_seasonal_EOF3_dual/")
results <- readRDS(here::here("analysis/vast_seasonal_EOF3_dual/results.rds"))
fit <- readRDS(here::here("analysis/vast_seasonal_EOF3_dual/fit.rds"))

data("ecomon_epu")
spp_list <- c("calfin", "chaeto", "cham", "clauso", "ctyp",
              "euph", "gas", "hyper", "larvaceans",
              "mlucens", "oithspp", "para", "pseudo", "tlong")

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


## Start the GLMS ----
library(ggfortify)
m <- glm(value ~ season + year, data = f1)
m1 <- glm(value ~ year * season, data = f1)
m2 <- glm(value ~ year * season, data = f1, weights = weight)


summary(m2)
acf(resid(m1))
autoplot(m2, label.size = 3)
AIC(m, m1, m2)


## Start the GAMS ----
yearseason_set <- expand.grid(year_season = paste0(rep(1992:2021, each = 4), "_",
                                                   c("winter", "spring", "summer", "fall")))

# results$Factors$Rotated_loadings$EpsilonTime1
factor_time <- data.frame(results$Factors$Rotated_loadings$EpsilonTime1,
                          results$Factors$Rotated_loadings_SE$EpsilonTime1,#fit$Report$Ltime_epsilon1_tf,
                          year_season = levels(zoop_dat$year_season)) %>%
  right_join(yearseason_set) %>%
  mutate(year = as.numeric(gsub("_.*", "", year_season)),
         season = factor(gsub(".*_", "", year_season),
                         levels = c("winter", "spring", "summer", "fall")),
         nseason = case_when(season == "winter" ~ 0,
                             season == "spring" ~ 0.25,
                             season == "summer" ~ 0.50,
                             season == "fall" ~ 0.75),
         time = year + nseason) %>%
  rename(`Factor_1` = X1,
         `Factor_2` = X2,
         `Factor_1SE` = X1.1,
         `Factor_2SE` = X2.1)

# fit <- NULL
# gc()
# factor_time <- factor_t


# # Location	ecoregions
# # Y	yr
# # M	mo
# # D	day
# # DayStart	j date
# # Date	doy
# # Percent	% possible pixels
# # Summary.Type	mean
# # Mission	GB is ocean color CHL, AV is SST
# # Mean	mean value
# # SE	SE of mean
# # n.pixels	number of pixels sampled
# remote_data <- read.csv(here::here("data/remote_data.csv"),
#                         stringsAsFactors = FALSE) %>%
#   dplyr::rename(EPU = Location,
#                 Year = Y,
#                 Month = M,
#                 Value = Mean,
#                 Var = Mission) %>%
#   dplyr::mutate(Var = dplyr::case_when(Var == "GB" ~ "CHL",
#                                        Var == "AV" ~ "SST",
#                                        TRUE ~ NA_character_),
#                 EPU = dplyr::case_when(EPU == "SCS" ~ "SS",
#                                        EPU == "NES" ~ "All",
#                                        EPU == "GBK" ~ "GB",
#                                        EPU == "MAB" ~ "MAB",
#                                        EPU == "GOM" ~ "GOM",
#                                        TRUE ~ NA_character_),
#                 EPU = as.character(EPU),
#                 nMonth = as.integer(Month)) %>%
#   dplyr::select(Year, nMonth, Var, Value, EPU) %>%
#   tidyr::pivot_wider(names_from = Var,
#                      values_from = Value)

#
# join_dat <- ecodata::chl_pp %>%
#   filter(#EPU %in% c("GOM"),
#          str_detect(Var, "MONTHLY_PPD_MEDIAN")) %>%
#   tidyr::separate(.,Time, into = c("Year","Month"), sep = 4, remove = FALSE) %>%
#   mutate(nMonth = as.integer(Month),
#          EPU = as.character(EPU),
#          Year = as.integer(Year)) %>%
#   dplyr::left_join(remote_data, by = c("Year", "EPU", "nMonth"))
#
# out_pp <- join_dat %>%
#   dplyr::mutate(log_value = log(Value + 1),
#          fMonth = plyr::mapvalues(Month, from = c("01","02","03","04","05","06",
#                                                  "07","08","09","10","11","12"),
#                                  to = c(month.abb)),
#          Date = as.Date(paste(Year, Month, "15", sep = "-")),# format = "%Y-%b-%d"),
#          Time = as.numeric(Date)/1000,
#          sYear = scale(Year, center = TRUE, scale = TRUE))


# out_pp$fMonth <- factor(out_pp$fMonth, levels = month.abb)

knots <- list(nSeason = seq(from = 0.5, to = 4.5, length.out = 4))


f1 <- factor_time %>%
  mutate(nseason = as.numeric(factor_time$season),
         time = time/1000,
         weight = 1/(Factor_1SE^2)) %>%
  select(value = Factor_1,
         weight,
         time,
         year,
         season,
         nseason)
#
ggplot(f1, aes(x = time, y = value, ymin = value - weight, ymax = value + weight, color = season)) +
  geom_point() +
  geom_errorbar()

library(mgcv)
library(gratia)
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

m1 <- gamm(value ~ s(year, bs = "cr", k = 20) + s(nseason, bs = "cc", k = 4) +
             te(year, nseason, bs = c("cr","cc"), k = c(20, 4)),
           data = f1, method = "ML", control = ctrl, knots = knots,
           correlation = corARMA(fixed = TRUE, form = ~ 1 | year, p = 1))
gratia::draw(m1)
gam.check(m1$gam)
appraise(m1$gam)
acf(resid(m1$lme))

pdat2 <- with(f1,
              data.frame(year = rep(1991:2021, each = 4),
                         nseason = rep(1:4, times = 31),
                         season = as.factor(rep(c("winter", "spring", "summer", "fall"), times = 31))))

pred2 <- predict(m1$gam, newdata = pdat2, se.fit = TRUE)
crit <- qt(0.975, df = df.residual(m1$gam)) # ~95% interval critical t

## add predictions & SEs to the new data ready for plotting
pdat2 <- transform(pdat2,
                   fitted = pred2$fit,  # predicted values
                   se = pred2$se.fit)#,   # standard errors
                   # season = factor(month.abb[nMonth], # month as a factor
                                   # levels = month.abb))
pdat2 <- transform(pdat2,
                   upper = fitted + (crit * se), # upper and...
                   lower = fitted - (crit * se)) # lower confidence bounds


p2 <- ggplot(pdat2, aes(x = year, y = fitted, group = season)) +
  geom_line(aes(colour = season)) +   # draw trend lines
  geom_ribbon(aes(ymin=lower, ymax = upper, fill = season), alpha=0.3) +
  geom_point(data = f1, aes(x = year, y = value, color = season)) +
  theme_bw() +                        # minimal theme
  theme(legend.position = "none") +   # no legend
  labs(caption = "sensu https://www.fromthebottomoftheheap.net/2015/11/23/are-some-seasons-warming-more-than-others/",
       subtitle = "GAMM tensor product interaction smooth",
       # y = expression("PP (gC m"^-2*" d"^-1*")"),
       x = "") +
  facet_wrap(~ season, ncol = 4) +    # facet on month
  # scale_y_continuous(breaks = seq(4, 17, by = 1),
  #                    minor_breaks = NULL) # nicer ticks
  NULL
p2



## basic gam ----
basic_model <- gam(value ~ season + s(year), data = f1, method = "REML")
plot(basic_model, all.terms = TRUE, pages = 1)
# summary(basic_model)$s.table
# acf(resid(basic_model))

factor_interact <- gam(value ~ season + s(year, bs = "ts", k = 25, by = season), data = f1, method = "REML")

plot(factor_interact, all.terms = TRUE, pages = 1)
summary(basic_model)$s.table
summary(factor_interact)$s.table
pacf(resid(factor_interact))

# acf(resid(factor_interact))
appraise(factor_interact)
gam.check(factor_interact)

## uncorrelated errors ----
# seasonal continuous
knots <- list(nSeason = seq(from = 0.5, to = 4.5, length.out = 4))
m <- gamm(value ~ s(year, bs = "ts") + s(nseason, k = 4, bs = "cc"),
          weights = weight,
          control = ctrl, knots = knots,
          method = "ML", data = f1)

# Shrink wiggliness to zero
m1 <- gamm(value ~ season + s(year, bs = "tp"),
           weights = weight,
           control = ctrl,
           method = "ML", data = f1)

# add factor interaction
m2 <- gamm(value ~ season + s(year, bs = "tp", by = season),
           weights = weight,
           control = ctrl,
           method = "ML", data = f1)

# add factor interaction and increase knots
m3 <- gamm(value ~ season + s(year, bs = "ts", k = 25, by = season),
           weights = weight,
           control = ctrl,
           method = "ML", data = f1)

# add random factor interaction and increase knots
m4 <- gamm(value ~ season + s(year, bs = "fs", k = 25, by = season, m = 4),
           weights = weight,
           method = "ML", data = f1)

AIC(m$lme, m1$lme, m2$lme, m3$lme, m4$lme)
## Based on the AIC
gam.check(m4$gam)
draw(m4)
summary(m4$gam)
appraise(m4$gam)
acf(resid(m4$gam))
acf(resid(m4$lme,type = "normalized"))

# add factor interaction and increase knots
m3ar1 <- gamm(value ~ season + s(time, bs = "ts", k = 25, by = season),
           weights = weight,
           control = ctrl,
           correlation = corARMA(form = ~ 1|year, p = 1),
           method = "ML", data = f1)

m3ar2 <- gamm(value ~ season + s(time, bs = "ts", k = 25, by = season),
              weights = weight,
              control = ctrl,
              correlation = corARMA(form = ~ 1|year, p = 2),
              method = "ML", data = f1)
m3ar3 <- gamm(value ~ season + s(time, bs = "ts", k = 25, by = season),
              weights = weight,
              control = ctrl,
              correlation = corARMA(form = ~ 1|year, p = 3),
              method = "ML", data = f1)
m3ar4 <- gamm(value ~ season + s(time, bs = "ts", k = 25, by = season),
              weights = weight,
              control = ctrl,
              correlation = corARMA(form = ~ 1|year, p = 4),
              method = "ML", data = f1)

AIC(m3ar1$lme, m3ar2$lme, m3ar3$lme, m3ar4$lme)
summary(m3ar1$gam)
appraise(m3ar1$gam)
draw(m3ar1$gam)
acf(residuals(m3ar1$lme, type = "normalized"))
acf(residuals(m2ar4$gam))
# add factor interaction
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

m2ar1 <- gamm(value ~ season + s(time, bs = "tp", by = season), method = "ML", weights = weight,
           data = f1,
           correlation = corARMA(form = ~ 1|year, p = 1),
           control = ctrl)

m2ar2 <- gamm(value ~ season + s(time, bs = "tp", by = season), method = "ML", weights = weight,
             data = f1,
             correlation = corARMA(form = ~ 1|year, p = 2),
             control = ctrl)

m2ar3 <- gamm(value ~ season + s(time, bs = "tp", by = season), method = "ML", weights = weight,
              data = f1,
              correlation = corARMA(form = ~ 1|year, p = 3),
              control = ctrl)

m2ar4 <- gamm(value ~ season + s(time, bs = "tp", by = season), method = "ML", weights = weight,
              data = f1,
              correlation = corARMA(form = ~ 1|year, p = 4),
              control = ctrl)


AIC(m2ar1$lme, m2ar2$lme, m2ar3$lme, m2ar4$lme)

summary(m2ar4$gam)
appraise(m2ar4$gam)
draw(m2ar4$gam)
acf(residuals(m2ar1$lme, type = "normalized"))
acf(residuals(m2ar4$gam))



ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

m2ar1 <- gamm(value ~ season + s(time, bs = "tp", by = season), method = "ML", weights = weight,
              data = f1,
              correlation = corARMA(form = ~ 1|season, p = 1),
              control = ctrl)

m2ar2 <- gamm(value ~ season + s(time, bs = "tp", by = season), method = "ML", weights = weight,
              data = f1,
              correlation = corARMA(form = ~ 1|season, p = 2),
              control = ctrl)

m2ar3 <- gamm(value ~ season + s(time, bs = "tp", by = season), method = "ML", weights = weight,
              data = f1,
              correlation = corARMA(form = ~ 1|season, p = 3),
              control = ctrl)

m2ar4 <- gamm(value ~ season + s(time, bs = "tp", by = season), method = "ML", weights = weight,
              data = f1,
              correlation = corARMA(form = ~ 1|season, p = 4),
              control = ctrl)


summary(m2ar4$gam)
appraise(m2ar4$gam)
draw(m2ar4$gam)
acf(residuals(m2ar1$lme, type = "normalized"))
acf(residuals(m2ar4$gam))


AIC(m2ar1$lme, m2ar2$lme, m2ar3$lme, m2ar4$lme)


ma1 <- gamm(value ~ s(time, bs = "ts") + s(nseason, k = 4, bs = "cc"), weights = weight,
          method = "ML", data = f1,
          correlation = corARMA(form = ~ 1|year, p = 1),
          control = ctrl)

ma2 <- gamm(value ~ s(time, bs = "ts") + s(nseason, k = 4, bs = "cc"), weights = weight,
            method = "ML", data = f1,
            correlation = corARMA(form = ~ 1|year, p = 2),
            control = ctrl)


AIC(ma1$lme, ma2$lme)
acf(residuals(ma1$lme, type = "normalized"))

#AR(1)
m1 <- gamm(value ~ s(nseason, k = 4, bs = "cc") + s(time, k = 20),
           data = f1, correlation = corARMA(form = ~ 1|year, p = 1),
           control = ctrl)

m3w <- gamm(value ~ season + s(year, bs = "tp", k = 25, by = season), method = "ML", weights = weight,
           data = f1)
AIC(m$lme, m3$lme, m3w$lme)
summary(m3w$lme)
plot(m3w$gam, scale = 0, all.terms = TRUE, pages = 1)
acf(resid(m3w$lme))
gam.check(m2w$gam)

summary(m2$gam)
plot(m4$gam, scale = 0, all.terms = TRUE, pages = 1)
acf(resid(m2$gam))
AIC(m$lme, m1$lme, m2$lme, m3$lme, m4$lme)

gam.check(m3$gam)
gam.check(m4$gam)
acf(resid(m2$lme, type = "normalized"))

pdat2 <- with(f1,
              data.frame(year = rep(1991:2021, each = 4),
                         nseason = rep(1:4, times = 31),
                         season = factor(rep(c("winter", "spring", "summer", "fall"), times = 31),
                                         levels = c("winter", "spring", "summer", "fall"))))

pred2 <- predict(factor_interact, newdata = pdat2, se.fit = TRUE)
crit <- qt(0.975, df = df.residual(factor_interact)) # ~95% interval critical t

## add predictions & SEs to the new data ready for plotting
pdat2 <- transform(pdat2,
                   fitted = pred2$fit,  # predicted values
                   se = pred2$se.fit)#,   # standard errors
# season = factor(month.abb[nMonth], # month as a factor
# levels = month.abb))
pdat2 <- transform(pdat2,
                   upper = fitted + (crit * se), # upper and...
                   lower = fitted - (crit * se)) # lower confidence bounds

levels(pdat2$season)





p2 <- ggplot(pdat2, aes(x = year, y = fitted, group = season)) +
  geom_line(aes(colour = season)) +   # draw trend lines
  geom_ribbon(aes(ymin=lower, ymax = upper, fill = season), alpha=0.3) +
  geom_point(data = f1, aes(x = year, y = value, color = season)) +
  theme_bw() +                        # minimal theme
  theme(legend.position = "none") +   # no legend
  labs(caption = "sensu https://www.fromthebottomoftheheap.net/2015/11/23/are-some-seasons-warming-more-than-others/",
       subtitle = "GAMM tensor product interaction smooth",
       # y = expression("PP (gC m"^-2*" d"^-1*")"),
       x = "") +
  facet_wrap(~ season, ncol = 1) +    # facet on month
  # scale_y_continuous(breaks = seq(4, 17, by = 1),
  #                    minor_breaks = NULL) # nicer ticks
  NULL
p2


## correlated errors ----

ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

#AR(1)
m1 <- gamm(value ~ s(nseason, k = 4, bs = "cc") + s(time, k = 20),
           data = f1, correlation = corARMA(form = ~ 1|year, p = 1),
           control = ctrl)


#AR(2)
m2 <- gamm(value ~ s(nseason, k = 4, bs = "cc") + s(time, k = 20),
           data = f1, correlation = corARMA(form = ~ 1|year, p = 2),
           control = ctrl)

#AR(3)
m3 <- gamm(value ~ s(nseason, k = 4, bs = "cc") + s(time, k = 20),
           data = f1, correlation = corARMA(form = ~ 1|year, p = 3),
           control = ctrl)


anova(m$lme, m1$lme, m2$lme, m3$lme)

layout(matrix(1:2, ncol = 2))
plot(m2$gam, scale = 0)
layout(1)

layout(matrix(1:2, ncol = 2))
acf(resid(m0$lme, type = "normalized"))
layout(1)

## Spline interactions ----

knots <- list(nSeason = seq(from = 0.5, to = 4.5, length.out = 4))

m0 <- gamm(value ~ te(year, nseason, bs = c("cr","cc"), k = c(10, 4)),
           data = f1, method = "REML", knots = knots)

ctrl <- list(niterEM = 0, optimMethod="L-BFGS-B", maxIter = 100, msMaxIter = 100)
m <- gamm(value ~ te(year, nseason, bs = c("cr","cc"), k = c(20, 4)),
          data = f1, method = "REML", control = ctrl, knots = knots,
          correlation = corARMA(form = ~ 1 | year, p = 4))

phi <- unname(intervals(m$lme, which = "var-cov")$corStruct[, 2])

m1 <- gamm(value ~ s(year, bs = "cr", k = 20) + s(nseason, bs = "cc", k = 4) +
            ti(year, nseason, bs = c("cr","cc"), k = c(20, 4)),
          data = f1, method = "ML", control = ctrl, knots = knots,
          correlation = corARMA(form = ~ 1 | year, p = 1))


m0 <- gamm(value ~ s(year, bs = "cr", k = 20) + s(nseason, bs = "cc", k = 4),
           data = f1, method = "ML", control = ctrl, knots = knots,
           correlation = corARMA(form = ~ 1 | year, p = 1))
anova(m0$lme, m1$lme)


m0 <- gam(value ~ ti(year, k = 10, bs = "cr") + ti(nseason, k = 4, bs = "cc") +
           ti(year, nseason, k = c(10, 4), bs = c("cr", "cc")),
         data = f1,
         knots = knots)

ctrl <- list(niterEM = 0, optimMethod="L-BFGS-B", maxIter = 100, msMaxIter = 100)
for (i in 1:8) {
  m <- gamm(value ~ te(year, nseason, bs = c("cr","cc"), k = c(10, 4)),
            data = f1, method = "REML", control = ctrl, knots = knots,
            correlation = corARMA(form = ~ 1 | year, p = i))
  assign(paste0("m", i), m)
}



m1 <- gamm(value ~ ti(year, k = 10, bs = "cr") + ti(nseason, k = 4, bs = "cc") +
            ti(year, nseason, k = c(10, 4), bs = c("cr", "cc")),
           correlation = corARMA(form = ~ 1|year, p = 1),
          # family = Gamma(link = "log"),
          data = f1,
          control = ctrl,
          knots = knots)

layout(matrix(1:2, ncol = 2))
acf(resid(m$lme), lag.max = 36, main = "ACF")
pacf(resid(m$lme), lag.max = 36, main = "pACF")
layout(1)

m2 <- gamm(value ~ ti(year, k = 10, bs = "cr") + ti(nseason, k = 4, bs = "cc") +
             ti(year, nseason, k = c(10, 4), bs = c("cr", "cc")),
           correlation = corARMA(form = ~ 1|year, p = 2),
           # family = Gamma(link = "log"),
           control = ctrl,
           data = f1,
           knots = knots)

layout(matrix(1:2, ncol = 2))
acf(resid(m2$lme), lag.max = 36, main = "ACF")
pacf(resid(m2$lme), lag.max = 36, main = "pACF")
layout(1)


layout(matrix(1:3, nrow = 1))
plot(m1$gam, scheme = 1, theta = 40)
plot(m1$gam, scheme = 1, theta = 80)
plot(m1$gam, scheme = 1, theta = 120)


m3 <- gamm(log_value ~ s(SST, k = 10, bs = "ts") + ti(Year, k = 10, bs = "cr") + ti(nMonth, k = 12, bs = "cc") +
             ti(Year, nMonth, k = c(10, 12), bs = c("cr", "cc")),
           correlation = corARMA(form = ~ 1|Year, p = 2),
           # family = Gamma(link = "log"),
           control = ctrl,
           data = out_pp,
           knots = knots)

summary(m3$gam)

gam.check(m3$gam)
layout(matrix(1:2, ncol = 2))
acf(resid(m3$lme), lag.max = 36, main = "ACF")
pacf(resid(m3$lme), lag.max = 36, main = "pACF")
layout(1)

anova(m2$lme, m3$lme)

plot(out_pp$SST, out_pp$nMonth, type = "p")
concurvity(m3$gam)
library(gratia)
gratia::draw(m3$gam)

gratia::derivatives(m2$gam)

#
# m <- gamm(Value ~ s(nMonth, bs = "cc", k = 12) +  s(Time),
#           data = out_pp)
#
# layout(matrix(1:2, ncol = 2))
# acf(resid(m$lme), lag.max = 36, main = "ACF")
# pacf(resid(m$lme), lag.max = 36, main = "pACF")
# layout(1)
#
# ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
#
#
# ## AR(1)
# m1 <- gamm(Value ~ s(nMonth, bs = "cc", k = 12)  + s(Time, k = 20),
#            data = out_pp, correlation = corARMA(form = ~ 1|Year, p = 1),
#            control = ctrl)
#
# ## AR(2)
# m2 <- gamm(Value ~ s(nMonth, bs = "cc", k = 12) + s(Time, k = 20),
#            data = out_pp, correlation = corARMA(form = ~ 1|Year, p = 2),
#            control = ctrl)
#
# ## AR(3)
# m3 <- gamm(Value ~ s(nMonth, bs = "cc", k = 12) + s(Time, k = 20),
#            data = out_pp, correlation = corARMA(form = ~ 1|Year, p = 3),
#            control = ctrl)
#
# anova(m$lme, m1$lme, m2$lme, m3$lme)
#
#
# arma_res <- forecast::auto.arima(resid(m$lme, type = "normalized"),
#                        stationary = TRUE, seasonal = FALSE)


# ## AR(1)
# m1ma <- gamm(Value ~ s(nMonth, bs = "cc", k = 12)  + s(Time),
#            data = out_pp, correlation = corARMA(form = ~ 1|Year, p = 1),
#            control = ctrl)
#
# gam.check(m1ma$gam)
#
# anova(m$lme, m1$lme, m1ma$lme)


## Climate change and spline interactions
## https://www.fromthebottomoftheheap.net/2015/11/21/climate-change-and-spline-interactions/


## naive model that assumes independence of observations
# m0 <- gamm(Value ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10, 12)),
#            data = out_pp, method = "REML", knots = knots)
# plot(acf(resid(m0$lme, type = "normalized")))
#
# layout(matrix(1:2, ncol = 2))
# plot(m0$gam)
# layout(1)
#

# ctrl <- list(niterEM = 0, optimMethod="L-BFGS-B", maxIter = 100, msMaxIter = 100)
#
# m1 <- gamm(Value ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10,12)),
#            data = out_pp, method = "REML", control = ctrl, knots = knots,
#            correlation = corARMA(form = ~ 1 | Year, p = 1))
#
# m2 <- gamm(Value ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10,12)),
#            data = out_pp, method = "REML", control = ctrl, knots = knots,
#            correlation = corARMA(form = ~ 1 | Year, p = 2))
#
# m3 <- gamm(Value ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10,12)),
#            data = out_pp, method = "REML", control = ctrl, knots = knots,
#            correlation = corARMA(form = ~ 1 | Year, p = 3))
#
# anova(m0$lme, m1$lme, m2$lme, m3$lme)
#

# m <- m1
# plot(m$gam, pers = TRUE)
#
# range(out_pp$Year)
#
#
#
# p3 <- ggplot(pdat2, aes(x = Year, y = fitted, group = fMonth)) +
#   geom_line(aes(colour = fMonth)) +   # draw trend lines
#   theme_bw() +                        # minimal theme
#   theme(legend.position = "top") +    # legend on top
#   scale_fill_discrete(name = "Month") + # nicer legend title
#   scale_colour_discrete(name = "Month") +
#   labs(y = expression(Temperature ~ (degree*C)), x = NULL) +
#   facet_grid(Quarter ~ ., scales = "free_y") # facet by Quarter
# p3

m0 <- gamm(log(Value) ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10,12)),
           data = out_pp, method = "REML", control = ctrl, knots = knots)

m1 <- gamm(log(Value) ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10,12)),
           data = out_pp, method = "REML", control = ctrl, knots = knots,
           correlation = corARMA(form = ~ 1 | Year, p = 1))

m2 <- gamm(log(Value) ~ te(Year, nMonth, bs = c("cr","cc"), k = c(10,12)),
           data = out_pp, method = "REML", control = ctrl, knots = knots,
           correlation = corARMA(form = ~ 1 | Year, p = 1, q = 1))

anova(m0$lme, m1$lme, m2$lme)

arma_res <- forecast::auto.arima(resid(m1$lme, type = "normalized"),
                                 stepwise = FALSE,
                                 approximation = FALSE,
                       stationary = TRUE, seasonal = FALSE)
vis.gam(m1$gam, theta=35)
arma_res$coef

layout(matrix(1:2, ncol = 2))
acf(resid(m1$lme), lag.max = 36, main = "ACF")
pacf(resid(m1$lme), lag.max = 36, main = "pACF")
layout(1)



knots <- list(nMonth = c(0.5, seq(1, 12, length = 10), 12.5))

m0 <- gam(log(Value) ~ s(Year, bs = "cr", k = 10) + s(nMonth, bs = "cc", k = 12),# +
             # ti(Year, nMonth, bs = c("cr","cc"), k = c(10, 12)),
           data = out_pp, method = "ML",
          # control = ctrl,
          knots = knots#,
           # correlation = corARMA(value = phi, fixed = TRUE, form = ~ 1 | Year, p = 1)
           )
summary(m0)
gam.check(m0)
gratia::appraise(m0)
forecast::ggPacf(m0$residuals)
forecast::ggAcf(m0$residuals)
forecast::auto.arima(m0$residuals)

m1 <- gamm(log(Value) ~ s(Year, bs = "cr", k = 10) + s(nMonth, bs = "cc", k = 12) +
             ti(Year, nMonth, bs = c("cr","cc"), k = c(10, 12)),
           data = out_pp, method = "ML", control = ctrl, knots = knots,
           correlation = corARMA(value = phi, fixed = TRUE, form = ~ 1 | Year, p = 1))

# phi <- unname(intervals(m1$lme, which = "var-cov")$corStruct[, 2])
## above pulls an error, so just use phi from summary(m$lme)
phi <- 0.5131951

m1 <- gamm(log(Value) ~ s(Year, bs = "cr", k = 10) + s(nMonth, bs = "cc", k = 12) +
             ti(Year, nMonth, bs = c("cr","cc"), k = c(10, 12)),
           data = out_pp, method = "ML", control = ctrl, knots = knots,
           correlation = corARMA(value = phi, fixed = TRUE, form = ~ 1 | Year, p = 1))
m0 <- gamm(log(Value) ~ s(Year, bs = "cr", k = 20) + s(nMonth, bs = "cc", k = 12),
           data = out_pp, method = "ML", control = ctrl, knots = knots,
           correlation = corARMA(value = phi, fixed = TRUE, form = ~ 1 | Year, p = 1))


arma_res <- forecast::auto.arima(resid(m1$lme, type = "normalized"),
                                 stationary = FALSE, seasonal = TRUE)


arma_res$coef

anova(m0$lme, m1$lme)


plot(m1$gam, scheme = TRUE)
m1$gam$model

phi <- 0.5131951

m1 <- gamm(log_value ~ s(Year, bs = "cr", k = 10) + s(nMonth, bs = "cc", k = 12) +
             ti(Year, nMonth, bs = c("cr","cc"), k = c(10, 12)),
           data = out_pp, method = "ML", control = ctrl, knots = knots,
           correlation = corARMA(value = phi, fixed = TRUE, form = ~ 1 | Year, p = 1))


pdat2 <- with(out_pp,
              data.frame(Year = rep(1998:2020, each = 12),
                         nMonth = rep(1:12, times = 23)))


pred2 <- predict(m1$gam, newdata = pdat2, se.fit = TRUE)
crit <- qt(0.975, df = df.residual(m1$gam)) # ~95% interval critical t

## add predictions & SEs to the new data ready for plotting
pdat2 <- transform(pdat2,
                   fitted = pred2$fit,  # predicted values
                   se = pred2$se.fit,   # standard errors
                   fMonth = factor(month.abb[nMonth], # month as a factor
                                   levels = month.abb))
pdat2 <- transform(pdat2,
                   upper = fitted + (crit * se), # upper and...
                   lower = fitted - (crit * se)) # lower confidence bounds


p2 <- ggplot(pdat2, aes(x = Year, y = fitted, group = fMonth)) +
  geom_line(aes(colour = fMonth)) +   # draw trend lines
  geom_ribbon(aes(ymin=lower, ymax = upper, fill = fMonth), alpha=0.3) +
  geom_point(data = out_pp, aes(x = Year, y = log_value, color = fMonth)) +
  theme_bw() +                        # minimal theme
  theme(legend.position = "none") +   # no legend
  labs(title, "GOM monthly median PPD",
       caption = "sensu https://www.fromthebottomoftheheap.net/2015/11/23/are-some-seasons-warming-more-than-others/",
       subtitle = "GAMM tensor product interaction smooth",
       y = expression("PP (gC m"^-2*" d"^-1*")"),
       x = "") +
  facet_wrap(~ fMonth, ncol = 6) +    # facet on month
  # scale_y_continuous(breaks = seq(4, 17, by = 1),
  #                    minor_breaks = NULL) # nicer ticks
  NULL
p2
 ggsave(p2, filename = "~/pp_GOM.pdf", width = 11, height = 8)


pdat2$Quarter <- NA
pdat2$Quarter[pdat2$nMonth %in% c(12,1,2)] <- "Winter"
pdat2$Quarter[pdat2$nMonth %in% 3:5] <- "Spring"
pdat2$Quarter[pdat2$nMonth %in% 6:8] <- "Summer"
pdat2$Quarter[pdat2$nMonth %in% 9:11] <- "Autumn"
pdat2 <- transform(pdat2,
                   Quarter = factor(Quarter,
                                    levels = c("Spring","Summer","Autumn","Winter")))




p2 <- ggplot(pdat2, aes(x = Year, y = fitted, group = fMonth)) +
  geom_line(aes(colour = fMonth)) +   # draw trend lines
  geom_ribbon(aes(ymin=lower, ymax = upper, fill = fMonth), alpha=0.3) +
  geom_point(data = out_pp, aes(x = Year, y = Value, color = fMonth)) +
  theme_bw() +                        # minimal theme
  theme(legend.position = "none") +   # no legend
  labs(y = expression("PP (gC m"^-2*" d"^-1*")")) +
  facet_wrap(~ fMonth, ncol = 6) +    # facet on month
  # scale_y_continuous(breaks = seq(4, 17, by = 1),
  #                    minor_breaks = NULL) # nicer ticks
  NULL
p2



layout(matrix(1:3, ncol = 3))
op <- par(mar = rep(4, 4) + 0.1)
plot(m1$gam, scheme = TRUE, scale = 0)
par(op)
layout(1)
