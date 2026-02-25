
# Load functions, packages, & data -----

# Define the packages
packages <- c("dplyr", "lubridate", "magrittr", "tidyr", "ggplot2", "ggpubr")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
# if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages])
# }

## Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Suppress dplyr summarise warning
options(dplyr.summarise.inform = FALSE)

# Make sure using dplyr select
select <- dplyr::select



# Load  the environmental data ----

data_ker <- read.csv("Data_inputs/bbal_ker_weekly_sst.csv")
#range(year(data_ker$date)) # 1981-2025

data_bi <- read.csv('Data_inputs/bbal_bi_weekly_sst.csv')
#range(year(data_bi$date)) # 1981-2025


# Compute variation ---- 

# For each spatial cell, compute mean, SD and CV of SST through time

## Create spatial cell ID
data_bi %<>% 
  mutate(cell = paste(x, y))

data_ker %<>% 
  mutate(cell = paste(x, y))


## Compute variability
cell_var_bi <- data_bi %>%
  group_by(cell) %>%
  summarise(
    mean_sst = mean(sst, na.rm = TRUE),
    sd_sst   = sd(sst, na.rm = TRUE),
    cv_sst   = ifelse(abs(mean_sst) < 1, NA, sd_sst / mean_sst), 
    n_obs    = n() ) %>%
  ungroup()

cell_var_ker <- data_ker %>%
  group_by(cell) %>%
  summarise(
    mean_sst = mean(sst, na.rm = TRUE),
    sd_sst   = sd(sst, na.rm = TRUE),
    cv_sst   = ifelse(abs(mean_sst) < 1, NA, sd_sst / mean_sst), 
    n_obs    = n() ) %>%
  ungroup()


## Summarise variability across space
site_summary_bi <- cell_var_bi %>%
  summarise(
    mean_sd_sst = mean(sd_sst, na.rm = TRUE), # mean temporal SD across cells 
    sd_sd_sst   = sd(sd_sst, na.rm = TRUE), # spatial hetereogeneity in temporal SD
    mean_cv_sst = mean(cv_sst, na.rm = TRUE), # mean coefficient of variation across spatial cells
    n_cells = n() ) %>%
  mutate(colony = "Bird Island")

site_summary_ker <- cell_var_ker %>%
  summarise(
    mean_sd_sst = mean(sd_sst, na.rm = TRUE), # mean temporal SD across cells 
    sd_sd_sst   = sd(sd_sst, na.rm = TRUE), # spatial hetereogeneity in temporal SD
    mean_cv_sst = mean(cv_sst, na.rm = TRUE), # mean coefficient of variation across spatial cells
    n_cells = n() ) %>%
  mutate(colony = "Kerguelen")

site_summary <- rbind(site_summary_bi, site_summary_ker)

save(site_summary, file = "Data_outputs/sst_variability_comparison")

# mean_sd_sst -> bird island has higher absolute temporal variability: higher absolute fluctuations in SST
# sd_sd_sst -> bird island has more heterogeneity across space; kerguelen is more spatially uniform
# mean_cv_sst -> bird island is much more variable relative to local average temperature; more variable thermal conditions



# Model variability across sites ----

cell_var_bi %<>% mutate(colony = "birdisland")
cell_var_ker %<>% mutate(colony = "kerguelen")

cell_var_all <- rbind(cell_var_bi, cell_var_ker)


lm_mean <- lm(mean_sst ~ colony, data = cell_var_all)
summary(lm_mean)

lm_var <- lm(sd_sst ~ colony, data = cell_var_all)
summary(lm_var)

lm_cv <- lm(cv_sst ~ colony, data = cell_var_all)
summary(lm_cv)


# Plot variability ----

## Split cells into coords
cell_var_xy.bi <- cell_var_bi %>%
  separate(cell, into = c("x", "y"), sep = " ") %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y) )

cell_var_xy.ker <- cell_var %>%
  separate(cell, into = c("x", "y"), sep = " ") %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y) )


## Plot gridmap

# bi range way too big - fix
plot_grid_bi <- ggplot(cell_var_xy.bi, aes(x = x, y = y, fill = cv_sst)) +
  geom_tile() +
  ylim(c(-60, -45)) +
  coord_equal() +
  scale_fill_viridis_c(
    name = "SST variability (CV)",
    option = "C") +
  labs(x = "Longitude", y = "Latitude", tag = "(A)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        legend.position = "none")

plot_grid_ker <- ggplot(cell_var_xy.ker, aes(x = x, y = y, fill = cv_sst)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis_c(
    name = "SST variability (CV)",
    option = "C") +
  labs(x = "Longitude", y = "Latitude", tag = "(B)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        plot.margin = margin(0, 0, 0, 0))

png(file = "Figures/FIGURES3.png", width = 12, height = 4, units = "in", res = 600)
ggarrange(plot_grid_bi, plot_grid_ker,
          nrow = 1, 
          ncol = 2,
          align = "hv",
          widths = c(1,1))
dev.off()
