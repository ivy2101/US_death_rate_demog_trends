library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggforce)
library(patchwork)
library(colorspace)
library(sf)
library(viridis)
library(usmap) 
library(shadowtext)

# Female 
f_mod <- readRDS("/Users/irisyang/Bayes_env/full_inla.rds") 
subsetf <- readRDS("/Users/irisyang/Bayes_env/full_data_female.rds")

# Male 
m_mod <- readRDS("/Users/irisyang/Bayes_env/full_inla_male.rds") 
subsetm <- readRDS("/Users/irisyang/Bayes_env/full_data_male.rds")

age_levels <- c("<1", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
                "35-39", "40-44", "45-49", "50-54", "55-59", "60-64",
                "65-69", "70-74", "75-79", "80-84", "85+")
race_levels <- c("White", "Black or African American", "Asian or Pacific Islander",
                 "American Indian or Alaska Native", "All other races")

allraces_fitted_oneage <- function(target_year, target_age, subset, sex, mod, age_levels, race_levels, 
                                   return_data = FALSE, global_min = NULL, global_max = NULL) {
  
  state_levels <- unique(subset$state)
  n_age <- length(age_levels)
  n_state <- length(state_levels)
  n_race <- length(race_levels)
  
  intercept <- mod$summary.fixed["(Intercept)", "mean"]
  global_slope <- mod$summary.fixed["x", "mean"]
  
  # Random effects
  state_eff <- mod$summary.random$state %>% mutate(state = state_levels) %>% rename(state_intercept = mean)
  state_slope <- mod$summary.random$state2 %>% mutate(state = state_levels) %>% rename(state_slope = mean)
  age_eff <- mod$summary.random$age %>% mutate(age = 1:n_age) %>% rename(baseline = mean)
  age_slope <- mod$summary.random$age2 %>% mutate(age = 1:n_age) %>% rename(age_slope = mean)
  race_eff <- mod$summary.random$race %>% mutate(race = race_levels) %>% rename(race_intercept = mean)
  race_slope <- mod$summary.random$race2 %>% mutate(race = race_levels) %>% rename(race_slope = mean)
  age_state_eff <- mod$summary.random$age5 %>% mutate(state = rep(state_levels, each = n_age), age = rep(1:n_age, times = n_state)) %>% rename(age_state_intercept = mean)
  age_state_slope <- mod$summary.random$age6 %>% mutate(state = rep(state_levels, each = n_age), age = rep(1:n_age, times = n_state)) %>% rename(age_state_slope = mean)
  age_race_eff <- mod$summary.random$age3 %>% mutate(race = rep(race_levels, each = n_age), age = rep(1:n_age, times = n_race)) %>% rename(age_race_intercept = mean)
  age_race_slope <- mod$summary.random$age4 %>% mutate(race = rep(race_levels, each = n_age), age = rep(1:n_age, times = n_race)) %>% rename(age_race_slope = mean)
  race_state_eff <- mod$summary.random$race5 %>% mutate(state = rep(state_levels, each = n_race), race = rep(race_levels, times = n_state)) %>% rename(race_state_intercept = mean)
  race_state_slope <- mod$summary.random$race6 %>% mutate(state = rep(state_levels, each = n_race), race = rep(race_levels, times = n_state)) %>% rename(race_state_slope = mean)
  
  id_3way_eff <- mod$summary.random$id_3way %>% rename(threeway_intercept = mean)
  id_3way_slope <- mod$summary.random$id_3way2 %>% rename(threeway_slope = mean)
  
  x_val <- target_year - 1972
  age_index <- which(age_levels == target_age)
  
  # Population weights
  weights <- subset %>%
    filter(year == target_year, age == target_age) %>%
    group_by(state, race) %>%
    summarise(population = sum(population, na.rm = TRUE), .groups = "drop") %>%
    group_by(state) %>%
    mutate(prop = population / sum(population, na.rm = TRUE)) %>%
    ungroup()
  
  if (nrow(weights) == 0) stop("No population data for selected age group.")
  
  # Expand grid
  fitted_all <- expand.grid(race = race_levels, state = state_levels) %>%
    mutate(age = age_index) %>%
    mutate(id = paste(state, age, race, sep = "_")) %>%
    left_join(age_eff, by = "age") %>%
    left_join(age_slope, by = "age") %>%
    left_join(state_eff, by = "state") %>%
    left_join(state_slope, by = "state") %>%
    left_join(race_eff, by = "race") %>%
    left_join(race_slope, by = "race") %>%
    left_join(age_state_eff, by = c("age", "state")) %>%
    left_join(age_state_slope, by = c("age", "state")) %>%
    left_join(age_race_eff, by = c("age", "race")) %>%
    left_join(age_race_slope, by = c("age", "race")) %>%
    left_join(race_state_eff, by = c("race", "state")) %>%
    left_join(race_state_slope, by = c("race", "state")) %>%
    left_join(
      tibble(
        id = paste(
          rep(state_levels, times = n_race),
          age_index,
          rep(race_levels, each = n_state),
          sep = "_"
        ),
        row = 1:(n_race * n_state)
      ) %>%
        left_join(id_3way_eff %>% mutate(row = row_number()), by = "row") %>%
        left_join(id_3way_slope %>% mutate(row = row_number()), by = "row") %>%
        select(id, threeway_intercept, threeway_slope),
      by = "id"
    ) %>%
    mutate(
      threeway_intercept = replace_na(threeway_intercept, 0),
      threeway_slope = replace_na(threeway_slope, 0)
    ) %>%
    left_join(weights, by = c("state", "race")) %>%
    mutate(
      prop = if_else(is.na(prop), 0, prop),
      linear_pred = intercept + baseline + state_intercept + race_intercept +
        age_race_intercept + race_state_intercept + age_state_intercept + threeway_intercept +
        x_val * (global_slope + state_slope + age_slope + race_slope +
                   age_state_slope + age_race_slope + race_state_slope + threeway_slope),
      fitted_rate = exp(linear_pred) * 100000
    )
  
  fitted_combined <- fitted_all %>%
    group_by(state) %>%
    summarise(fitted_rate = sum(fitted_rate * prop, na.rm = TRUE), .groups = "drop")
  
  if (return_data) {
    return(fitted_combined)
  }
  
  us_map_data <- us_map("states")
  
  joined_fitted <- us_map_data %>%
    left_join(fitted_combined, by = c("abbr" = "state"))
  
  # Use global min/max if provided
  min_val <- if (!is.null(global_min)) global_min else min(fitted_combined$fitted_rate, na.rm = TRUE)
  max_val <- if (!is.null(global_max)) global_max else max(fitted_combined$fitted_rate, na.rm = TRUE)
  
  p_fitted <- plot_usmap(data = joined_fitted, values = "fitted_rate", labels = FALSE) +
    scale_fill_continuous_sequential(
      palette = "Mako", p1 = 1, p2 = 1.1,
      limits = c(min_val, max_val),
      name = "Death Rate", na.value = "grey70"
    ) +
    geom_sf_text(aes(label = abbr), size = 3, color = "black", check_overlap = TRUE) +
    labs(
      title = paste0(target_age, " ", sex, " in ", target_year)
    ) +
    theme_void() +
    theme(plot.title = element_text(size = 18, hjust = 0.5))
  
  return(p_fitted)
}

# Get the fitted data separately
fit_1973 <- allraces_fitted_oneage(target_year = 1973, target_age = "65-69", subset = subsetf,
                                   sex = "Females", mod = f_mod, age_levels = age_levels, race_levels = race_levels,
                                   return_data = TRUE)

fit_2022 <- allraces_fitted_oneage(target_year = 2022, target_age = "65-69", subset = subsetf,
                                   sex = "Females", mod = f_mod, age_levels = age_levels, race_levels = race_levels,
                                   return_data = TRUE)

# Find shared global scale
global_min <- min(fit_1973$fitted_rate, fit_2022$fitted_rate, na.rm = TRUE)
global_max <- max(fit_1973$fitted_rate, fit_2022$fitted_rate, na.rm = TRUE)

# Now plot both using SAME scale
p1973 <- allraces_fitted_oneage(target_year = 1973, target_age = "65-69", subset = subsetf,
                                sex = "Females", mod = f_mod, age_levels = age_levels, race_levels = race_levels,
                                global_min = global_min, global_max = global_max)

p2022 <- allraces_fitted_oneage(target_year = 2022, target_age = "65-69", subset = subsetf,
                                sex = "Females", mod = f_mod, age_levels = age_levels, race_levels = race_levels,
                                global_min = global_min, global_max = global_max)

# View side by side
p1973 + p2022


# fit_1973$year <- 1973
# fit_2022$year <- 2022
# combined_fit <- bind_rows(fit_1973, fit_2022)
# 
# # Find the row with the maximum fitted rate
# max_row <- combined_fit %>%
#   filter(fitted_rate == max(fitted_rate, na.rm = TRUE))
# 
# # View the result
# print(max_row)

# Plots by Race 
allages_fitted_raceonly <- function(target_year, subset, mod, sex, age_levels, race_levels, 
                                    target_race, return_data = FALSE, global_min = NULL, global_max = NULL) {
  
  x_val <- target_year - 1972
  
  age_levels <- levels(subset$age)
  race_levels <- levels(subset$race)
  state_levels <- unique(subset$state)
  
  n_age <- length(age_levels)
  n_race <- length(race_levels)
  n_state <- length(state_levels)
  
  # --- MODEL COMPONENTS ---
  intercept <- mod$summary.fixed["(Intercept)", "mean"]
  global_slope <- mod$summary.fixed["x", "mean"]
  
  state_eff <- mod$summary.random$state %>% mutate(state = state_levels) %>% rename(state_intercept = mean)
  state_slope <- mod$summary.random$state2 %>% mutate(state = state_levels) %>% rename(state_slope = mean)
  age_eff <- mod$summary.random$age %>% mutate(age = 1:n_age) %>% rename(baseline = mean)
  age_slope <- mod$summary.random$age2 %>% mutate(age = 1:n_age) %>% rename(age_slope = mean)
  race_eff <- mod$summary.random$race %>% mutate(race = race_levels) %>% rename(race_intercept = mean)
  race_slope <- mod$summary.random$race2 %>% mutate(race = race_levels) %>% rename(race_slope = mean)
  age_state_eff <- mod$summary.random$age5 %>% mutate(state = rep(state_levels, each = n_age), age = rep(1:n_age, times = n_state)) %>% rename(age_state_intercept = mean)
  age_state_slope <- mod$summary.random$age6 %>% mutate(state = rep(state_levels, each = n_age), age = rep(1:n_age, times = n_state)) %>% rename(age_state_slope = mean)
  age_race_eff <- mod$summary.random$age3 %>% mutate(race = rep(race_levels, each = n_age), age = rep(1:n_age, times = n_race)) %>% rename(age_race_intercept = mean)
  age_race_slope <- mod$summary.random$age4 %>% mutate(race = rep(race_levels, each = n_age), age = rep(1:n_age, times = n_race)) %>% rename(age_race_slope = mean)
  race_state_eff <- mod$summary.random$race5 %>% mutate(state = rep(state_levels, each = n_race), race = rep(race_levels, times = n_state)) %>% rename(race_state_intercept = mean)
  race_state_slope <- mod$summary.random$race6 %>% mutate(state = rep(state_levels, each = n_race), race = rep(race_levels, times = n_state)) %>% rename(race_state_slope = mean)
  
  id_3way_eff <- mod$summary.random$id_3way %>% rename(threeway_intercept = mean)
  id_3way_slope <- mod$summary.random$id_3way2 %>% rename(threeway_slope = mean)
  
  us_map_data <- us_map("states")
  
  # --- Find index of selected race ---
  race_index <- which(race_levels == target_race)
  
  # --- POPULATION WEIGHTS ---
  weights <- subset %>%
    filter(year == target_year) %>%
    group_by(state, age) %>%
    summarise(population = sum(population, na.rm = TRUE), .groups = "drop") %>%
    group_by(state) %>%
    mutate(prop = population / sum(population, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(age = as.integer(factor(age, levels = age_levels)))
  
  # --- COMPUTE FITTED RATES ---
  fitted_all <- expand.grid(
    age = age_levels,
    state = state_levels
  ) %>%
    mutate(
      age = as.integer(factor(age, levels = age_levels)),
      race = race_levels[race_index]
    ) %>%
    mutate(id = paste(state, age, race, sep = "_")) %>%
    left_join(age_eff, by = "age") %>%
    left_join(age_slope, by = "age") %>%
    left_join(state_eff, by = "state") %>%
    left_join(state_slope, by = "state") %>%
    left_join(race_eff, by = "race") %>%
    left_join(race_slope, by = "race") %>%
    left_join(age_state_eff, by = c("age", "state")) %>%
    left_join(age_state_slope, by = c("age", "state")) %>%
    left_join(age_race_eff, by = c("age", "race")) %>%
    left_join(age_race_slope, by = c("age", "race")) %>%
    left_join(race_state_eff, by = c("race", "state")) %>%
    left_join(race_state_slope, by = c("race", "state")) %>%
    left_join(
      tibble(
        id = paste(
          rep(state_levels, times = n_age),
          rep(1:n_age, each = n_state),
          race_levels[race_index],
          sep = "_"
        ),
        row = 1:(n_age * n_state)
      ) %>%
        left_join(id_3way_eff %>% mutate(row = row_number()), by = "row") %>%
        left_join(id_3way_slope %>% mutate(row = row_number()), by = "row") %>%
        select(id, threeway_intercept, threeway_slope),
      by = "id"
    ) %>%
    mutate(
      threeway_intercept = replace_na(threeway_intercept, 0),
      threeway_slope = replace_na(threeway_slope, 0)
    ) %>%
    left_join(weights, by = c("state", "age")) %>%
    mutate(
      prop = if_else(is.na(prop), 0, prop),
      linear_pred = intercept + baseline + state_intercept + age_state_intercept +
        race_intercept + age_race_intercept + race_state_intercept + threeway_intercept +
        x_val * (global_slope + state_slope + age_slope + age_state_slope +
                   race_slope + age_race_slope + race_state_slope + threeway_slope),
      fitted_rate = exp(linear_pred) * 100000
    )
  
  fitted_combined <- fitted_all %>%
    group_by(state) %>%
    summarise(fitted_rate = sum(fitted_rate * prop, na.rm = TRUE), .groups = "drop")
  
  if (return_data) {
    return(fitted_combined)
  }
  
  # Otherwise, make a plot
  fitted_plot <- us_map_data %>%
    left_join(fitted_combined, by = c("abbr" = "state")) %>%
    ggplot(aes(fill = fitted_rate)) +
    geom_sf(color = "black") +
    scale_fill_continuous_sequential(
      palette = "Mako", p1 = 1, p2 = 1.1,
      limits = if (!is.null(global_min) && !is.null(global_max)) c(global_min, global_max) else NULL,
      name = "Death Rate", na.value = "grey70"
    ) +
    geom_sf_text(aes(label = abbr), size = 3, color = "black", check_overlap = TRUE) +
    labs(
      title = paste0(target_race, " ", sex, " in ", target_year)
    ) +
    theme_void() +
    theme(plot.title = element_text(size = 16, hjust = 0.5))
  
  return(fitted_plot)
}

# Get fitted data for Black and White
black_fitted <- allages_fitted_raceonly(target_year = 2022, subset = subsetm, mod = m_mod, 
                                        sex = "Males", age_levels = age_levels, 
                                        race_levels = race_levels,
                                        target_race = "Black or African American",
                                        return_data = TRUE)

white_fitted <- allages_fitted_raceonly(target_year = 2022, subset = subsetm, mod = m_mod, 
                                        sex = "Males", age_levels = age_levels, 
                                        race_levels = race_levels,
                                        target_race = "White",
                                        return_data = TRUE)

# Find global min and max
global_min <- min(black_fitted$fitted_rate, white_fitted$fitted_rate, na.rm = TRUE)
global_max <- max(black_fitted$fitted_rate, white_fitted$fitted_rate, na.rm = TRUE)

# Now plot with SAME scale
black_plot <- allages_fitted_raceonly(target_year = 2022, subset = subsetm, mod = m_mod, 
                                      sex = "Males", age_levels = age_levels, 
                                      race_levels = race_levels,
                                      target_race = "Black or African American",
                                      global_min = global_min, global_max = global_max)

white_plot <- allages_fitted_raceonly(target_year = 2022, subset = subsetm, mod = m_mod, 
                                      sex = "Males", age_levels = age_levels, 
                                      race_levels = race_levels,
                                      target_race = "White",
                                      global_min = global_min, global_max = global_max)

# Display side by side

black_plot + white_plot
