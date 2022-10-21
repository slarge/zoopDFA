library(MARSS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

source(here::here("R/helper_functions.r"))


dfa_out <- readRDS(here::here("analysis/zoop_seasonal_dfa_epu.rds"))


aicc_mod <- dfa_out %>%
  select(-mod,
         -data) %>%
  tidyr::unnest(AICc)

aicc_plot <- ggplot(aicc_mod, aes(x = m, y = AICc, color = covariate, shape = R)) +
  geom_point() +
  geom_line() +
  facet_wrap(~epu, scales = "free_y", ncol = 1) +
  theme_minimal()

# all_mod <- dfa_out %>%
#   select(-mod,
#          -data) %>%
#   tidyr::unnest(AICc) %>%
#   group_by(EPU) %>%
#   arrange(AICc, .by_group = TRUE) %>%
#   slice_min(order_by = AICc)

table_mod <- dfa_out %>%
  select(-mod,
         -data) %>%
  tidyr::unnest(AICc) %>%
  group_by(epu) %>%
  arrange(AICc, .by_group = TRUE) %>%
  slice_min(order_by = AICc)

best_mod <- dfa_out %>%
  tidyr::unnest(AICc) %>%
  group_by(epu) %>%
  arrange(AICc, .by_group = TRUE) %>%
  slice_min(order_by = AICc)

plot_tibble <- best_mod %>%
  mutate(fitted.ytT = purrr::map(.x = mod, function(x) autoplot.marssMLE(x,
                                                                         plot.type = "fitted.ytT",
                                                                         form = "dfa",
                                                                         rotate = TRUE,
                                                                         silent = TRUE)),
         dfa_plot = purrr::pmap(.l = list(object = mod, EPU = epu),
                                .f = function(object, EPU) dfa_plot(object = object, EPU = epu))) %>%
  select(-mod,
         -data)

ggsave(filename = "analysis/figures/gom_dfa.png", plot_tibble$dfa_plot[[1]])
ggsave(filename = "analysis/figures/gb_dfa.png", plot_tibble$dfa_plot[[2]])
ggsave(filename = "analysis/figures/mab_dfa.png", plot_tibble$dfa_plot[[3]])
