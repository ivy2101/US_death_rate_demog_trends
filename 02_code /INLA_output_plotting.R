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

# Collapsed and Standardized over Age and Race 
national_plot_weighted <- function(target_year, subset, mod, sex, age_levels, race_levels) {
  
  x_val <- target_year - 1972
  
  state_levels <- unique(subset$state)
  n_age <- length(age_levels)
  n_state <- length(state_levels)
  n_race <- length(race_levels)
  
  # Fixed Effects
  intercept <- mod$summary.fixed["(Intercept)", "mean"]
  global_slope <- mod$summary.fixed["x", "mean"]
  
  # Random Effects
  state_eff <- mod$summary.random$state %>%
    mutate(state = state_levels) %>%
    rename(state_intercept = mean)
  
  state_slope <- mod$summary.random$state2 %>%
    mutate(state = state_levels) %>%
    rename(state_slope = mean)
  
  age_eff <- mod$summary.random$age %>%
    mutate(age = 1:n_age) %>%
    rename(baseline = mean)
  
  age_slope <- mod$summary.random$age2 %>%
    mutate(age = 1:n_age) %>%
    rename(age_slope = mean)
  
  race_eff <- mod$summary.random$race %>%
    mutate(race = race_levels) %>%
    rename(race_intercept = mean)
  
  race_slope <- mod$summary.random$race2 %>%
    mutate(race = race_levels) %>%
    rename(race_slope = mean)
  
  age_state_eff <- mod$summary.random$age5 %>%
    mutate(state = rep(state_levels, each = n_age),
           age = rep(1:n_age, times = n_state)) %>%
    rename(age_state_intercept = mean)
  
  age_state_slope <- mod$summary.random$age6 %>%
    mutate(state = rep(state_levels, each = n_age),
           age = rep(1:n_age, times = n_state)) %>%
    rename(age_state_slope = mean)
  
  age_race_eff <- mod$summary.random$age3 %>%
    mutate(race = rep(race_levels, each = n_age),
           age = rep(1:n_age, times = n_race)) %>%
    rename(age_race_intercept = mean)
  
  age_race_slope <- mod$summary.random$age4 %>%
    mutate(race = rep(race_levels, each = n_age),
           age = rep(1:n_age, times = n_race)) %>%
    rename(age_race_slope = mean)
  
  race_state_eff <- mod$summary.random$race5 %>%
    mutate(state = rep(state_levels, each = n_race),
           race = rep(race_levels, times = n_state)) %>%
    rename(race_state_intercept = mean)
  
  race_state_slope <- mod$summary.random$race6 %>%
    mutate(state = rep(state_levels, each = n_race),
           race = rep(race_levels, times = n_state)) %>%
    rename(race_state_slope = mean)
  
  id_3way_eff <- mod$summary.random$id_3way %>%
    rename(threeway_intercept = mean)
  
  id_3way_slope <- mod$summary.random$id_3way2 %>%
    rename(threeway_slope = mean)
  
  # --- POPULATION WEIGHTS ---
  weights <- subset %>%
    filter(year == target_year) %>%
    group_by(state, age, race) %>%
    summarise(population = sum(population, na.rm = TRUE), .groups = "drop") %>%
    group_by(state) %>%
    mutate(prop = population / sum(population, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      age = as.integer(factor(age, levels = age_levels)),
      race = factor(race, levels = race_levels)
    ) %>%
    select(state, age, race, prop)
  
  # --- COMPUTE FITTED ---
  fitted_all <- expand.grid(
    race = race_levels,
    age = age_levels,
    state = state_levels
  ) %>%
    mutate(
      age = as.integer(factor(age, levels = age_levels)),
      race = factor(race, levels = race_levels),
      id = paste(state, age, race, sep = "_")
    ) %>%
    left_join(age_eff %>% select(age, baseline), by = "age") %>%
    left_join(age_slope %>% select(age, age_slope), by = "age") %>%
    left_join(state_eff %>% select(state, state_intercept), by = "state") %>%
    left_join(state_slope %>% select(state, state_slope), by = "state") %>%
    left_join(race_eff %>% select(race, race_intercept), by = "race") %>%
    left_join(race_slope %>% select(race, race_slope), by = "race") %>%
    left_join(age_state_eff %>% select(age, state, age_state_intercept), by = c("age", "state")) %>%
    left_join(age_state_slope %>% select(age, state, age_state_slope), by = c("age", "state")) %>%
    left_join(age_race_eff %>% select(age, race, age_race_intercept), by = c("age", "race")) %>%
    left_join(age_race_slope %>% select(age, race, age_race_slope), by = c("age", "race")) %>%
    left_join(race_state_eff %>% select(race, state, race_state_intercept), by = c("race", "state")) %>%
    left_join(race_state_slope %>% select(race, state, race_state_slope), by = c("race", "state")) %>%
    left_join(
      tibble(
        id = paste(
          rep(state_levels, each = n_age * n_race),
          rep(rep(1:n_age, each = n_race), times = n_state),
          rep(race_levels, times = n_age * n_state),
          sep = "_"
        ),
        row = 1:(n_state * n_age * n_race)
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
    left_join(weights, by = c("state", "age", "race")) %>%
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
  
  # --- RAW DATA ---
  raw_data <- subset %>%
    filter(year == target_year) %>%
    mutate(
      age = as.integer(factor(age, levels = age_levels)),
      race = factor(race, levels = race_levels)
    ) %>%
    left_join(weights, by = c("state", "age", "race")) %>%
    mutate(
      deaths = if_else(is.na(deaths), 0, deaths),
      population = if_else(is.na(population), 0, population),
      prop = if_else(is.na(prop), 0, prop),
      raw_rate = if_else(population > 0, (deaths / population) * 100000, 0)
    ) %>%
    group_by(state) %>%
    summarise(
      obs_rate = sum(raw_rate * prop, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- PLOT ---
  min_val <- min(c(fitted_combined$fitted_rate, raw_data$obs_rate), na.rm = TRUE)
  max_val <- max(c(fitted_combined$fitted_rate, raw_data$obs_rate), na.rm = TRUE)
  
  us_map_data <- us_map("states")
  
  joined_fitted <- us_map_data %>%
    left_join(fitted_combined, by = c("abbr" = "state"))
  
  joined_raw <- us_map_data %>%
    left_join(raw_data, by = c("abbr" = "state"))
  
  p_fitted <- plot_usmap(data = joined_fitted, values = "fitted_rate", labels = FALSE) +
    scale_fill_continuous_sequential(
      palette = "Mako", p1 = 1, p2 = 1.1,
      limits = c(min_val, max_val),
      name = "Death Rate",
      na.value = "grey70"
    ) +
    geom_sf_text(aes(label = abbr), size = 1, color = "black", check_overlap = TRUE) +
    labs(title = paste0("Fitted")) +
    theme_void()
  
  p_raw <- plot_usmap(data = joined_raw, values = "obs_rate", labels = FALSE) +
    scale_fill_continuous_sequential(
      palette = "Mako", p1 = 1, p2 = 1.1,
      limits = c(min_val, max_val),
      name = "Death Rate",
      na.value = "grey70"
    ) +
    geom_sf_text(aes(label = abbr), size = 1.6, color = "black", check_overlap = TRUE) +
    labs(title = paste0("Raw")) +
    theme_void()
  
  (p_raw | p_fitted) +
    plot_annotation(
      title = paste0(" ", sex, " In ", target_year, " — Age-Standardized"),
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
    )
}

# Collapsed and Standardized by Race 
allraces_per_age<- function(target_year, subset, sex, mod, age_levels, race_levels) {
  
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
  
  for (i in seq_along(age_levels)) {
    target_age <- age_levels[i]
    age_index <- i
    
    # Population weights
    weights <- subset %>%
      filter(year == target_year, age == target_age) %>%
      group_by(state, race) %>%
      summarise(population = sum(population, na.rm = TRUE), .groups = "drop") %>%
      group_by(state) %>%
      mutate(prop = population / sum(population, na.rm = TRUE)) %>%
      ungroup()
    
    if (nrow(weights) == 0) next  # skip if no data
    
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
        prop = if_else(is.na(prop), 0, prop),  # Add this to protect against missing weights
        linear_pred = intercept + baseline + state_intercept + race_intercept + 
          age_race_intercept + race_state_intercept + age_state_intercept + threeway_intercept +
          x_val * (global_slope + state_slope + age_slope + race_slope +
                     age_state_slope + age_race_slope + race_state_slope + threeway_slope),
        fitted_rate = exp(linear_pred) * 100000
      )
    
    fitted_combined <- fitted_all %>%
      group_by(state) %>%
      summarise(fitted_rate = sum(fitted_rate * prop, na.rm = TRUE), .groups = "drop")
    
    raw_data <- subset %>%
      filter(year == target_year, age == target_age) %>%
      group_by(state, race) %>%
      summarise(
        deaths = sum(deaths, na.rm = TRUE),
        population = sum(population, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        raw_rate = if_else(population > 0, (deaths / population) * 100000, 0)
      ) %>%
      left_join(weights, by = c("state", "race")) %>%  
      mutate(
        prop = if_else(is.na(prop), 0, prop),           # in case any NA
        weighted_rate = raw_rate * prop
      ) %>%
      group_by(state) %>%
      summarise(
        obs_rate = sum(weighted_rate, na.rm = TRUE),
        .groups = "drop"
      )
    
    
    min_val <- min(c(fitted_combined$fitted_rate, raw_data$obs_rate), na.rm = TRUE)
    max_val <- max(c(fitted_combined$fitted_rate, raw_data$obs_rate), na.rm = TRUE)
    
    us_map_data <- us_map("states")
    
    joined_fitted <- us_map_data %>% left_join(fitted_combined, by = c("abbr" = "state"))
    joined_raw <- us_map_data %>% left_join(raw_data, by = c("abbr" = "state"))
    
    p_fitted <- plot_usmap(data = joined_fitted, values = "fitted_rate", labels = FALSE) +
      scale_fill_continuous_sequential(
        palette = "Mako", p1 = 1, p2 = 1.1,
        limits = c(min_val, max_val),
        name = "Death Rate",
        na.value = "grey70"
      ) +
      geom_sf_text(aes(label = abbr), size = 1, color = "black", check_overlap = TRUE) +
      labs(title = paste0("Fitted")) +
      theme_void()
    
    p_raw <- plot_usmap(data = joined_raw, values = "obs_rate", labels = FALSE) +
      scale_fill_continuous_sequential(
        palette = "Mako", p1 = 1, p2 = 1.1,
        limits = c(min_val, max_val),
        name = "Death Rate",
        na.value = "grey70"
      ) +
      geom_sf_text(aes(label = abbr), size = 1.6, color = "black", check_overlap = TRUE) +
      labs(title = paste0("Observed")) +
      theme_void()
    
    p_final <- (p_raw | p_fitted) +
      plot_annotation(
        title = paste0(target_age, " " ,sex, " In ", target_year, " — Weighted"),
        theme = theme(plot.title = element_text(size = 18, hjust = 0.5))
      )
    
    print(p_final)
    
  }
}

# Collapsed and Standardized by Age 
allages_perrace <- function(target_year, subset, mod, sex, age_levels, race_levels) {
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
  age_state_eff <- mod$summary.random$age5 %>% mutate(state = rep(state_levels, each = n_age),
                                                      age = rep(1:n_age, times = n_state)) %>% rename(age_state_intercept = mean)
  age_state_slope <- mod$summary.random$age6 %>% mutate(state = rep(state_levels, each = n_age),
                                                        age = rep(1:n_age, times = n_state)) %>% rename(age_state_slope = mean)
  age_race_eff <- mod$summary.random$age3 %>% mutate(race = rep(race_levels, each = n_age),
                                                     age = rep(1:n_age, times = n_race)) %>% rename(age_race_intercept = mean)
  age_race_slope <- mod$summary.random$age4 %>% mutate(race = rep(race_levels, each = n_age),
                                                       age = rep(1:n_age, times = n_race)) %>% rename(age_race_slope = mean)
  race_state_eff <- mod$summary.random$race5 %>% mutate(state = rep(state_levels, each = n_race),
                                                        race = rep(race_levels, times = n_state)) %>% rename(race_state_intercept = mean)
  race_state_slope <- mod$summary.random$race6 %>% mutate(state = rep(state_levels, each = n_race),
                                                          race = rep(race_levels, times = n_state)) %>% rename(race_state_slope = mean)
  id_3way_eff <- mod$summary.random$id_3way %>% rename(threeway_intercept = mean)
  id_3way_slope <- mod$summary.random$id_3way2 %>% rename(threeway_slope = mean)
  
  us_map_data <- us_map("states")
  
  plots <- list() # Store plots here
  
  # --- LOOP THROUGH RACES ---
  for (i in seq_along(race_levels)) {
    target_race <- race_levels[i]
    race_index <- i
    
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
            race_levels[race_index],  # single fixed race
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
        prop = if_else(is.na(prop), 0, prop),  # <<< ADDED this line for safety
        linear_pred = intercept + baseline + state_intercept + age_state_intercept +
          race_intercept + age_race_intercept + race_state_intercept + threeway_intercept +
          x_val * (global_slope + state_slope + age_slope + age_state_slope +
                     race_slope + age_race_slope + race_state_slope + threeway_slope),
        fitted_rate = exp(linear_pred) * 100000
      )
    
    
    fitted_combined <- fitted_all %>%
      group_by(state) %>%
      summarise(fitted_rate = sum(fitted_rate * prop, na.rm = TRUE), .groups = "drop")
    
    raw_data <- subset %>%
      filter(year == target_year, race == target_race) %>%
      mutate(age = as.integer(factor(age, levels = age_levels))) %>%
      group_by(state, age) %>%
      summarise(
        deaths = sum(deaths, na.rm = TRUE),
        population = sum(population, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        raw_rate = if_else(population > 0, (deaths / population) * 100000, 0)
      ) %>%
      left_join(weights, by = c("state", "age")) %>%
      mutate(
        prop = if_else(is.na(prop), 0, prop),   # safety in case of missing weights
        weighted_rate = raw_rate * prop
      ) %>%
      group_by(state) %>%
      summarise(
        obs_rate = sum(weighted_rate, na.rm = TRUE),
        .groups = "drop"
      )
    
    
    min_val <- min(c(fitted_combined$fitted_rate, raw_data$obs_rate), na.rm = TRUE)
    max_val <- max(c(fitted_combined$fitted_rate, raw_data$obs_rate), na.rm = TRUE)
    
    joined_fitted <- us_map_data %>%
      left_join(fitted_combined, by = c("abbr" = "state"))
    
    joined_raw <- us_map_data %>%
      left_join(raw_data, by = c("abbr" = "state"))
    
    p_fitted <- plot_usmap(data = joined_fitted, values = "fitted_rate", labels = FALSE) +
      scale_fill_continuous_sequential(
        palette = "Mako", p1 = 1, p2 = 1.1,
        limits = c(min_val, max_val),
        name = "Death Rate",
        na.value = "grey70"
      ) +
      geom_sf_text(aes(label = abbr), size = 1, color = "black", check_overlap = TRUE) +
      labs(title = paste0("Fitted")) +
      theme_void()
    
    p_raw <- plot_usmap(data = joined_raw, values = "obs_rate", labels = FALSE) +
      scale_fill_continuous_sequential(
        palette = "Mako", p1 = 1, p2 = 1.1,
        limits = c(min_val, max_val),
        name = "Death Rate",
        na.value = "grey70"
      ) +
      geom_sf_text(aes(label = abbr), size = 1.6, color = "black", check_overlap = TRUE) +
      labs(title = paste0("Observed")) +
      theme_void()
    
    plots[[target_race]] <- (p_raw | p_fitted) +
      plot_annotation(
        title = paste0(target_race, " ", sex, " in ", target_year, " — Weighted"),
        theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
      )
  }
  
  return(plots)
}

# Aggregated Plots  
individual_4pplots <- function(data, mod, sex, age = age_levels, race = race_levels) {
  # Add fitted and observed rates
  data$fitted <- mod$summary.fitted.values$mean * 100000
  data$obs_rate <- (data$deaths / data$population) * 100000
  
  us_map_data <- us_map("states")
  
  for (r in race) {
    for (a in age) {
      
      combined_data <- data %>%
        filter(year %in% c(1973, 2022), race == r, age == a) %>%
        group_by(year, state) %>%
        summarise(
          mean_fitted = mean(fitted, na.rm = TRUE),
          obs_fitted = mean(obs_rate, na.rm = TRUE),
          .groups = "drop"
        )
      
      if (nrow(combined_data) == 0) next  # Skip if no data
      
      min_val <- min(combined_data$obs_fitted, combined_data$mean_fitted, na.rm = TRUE)
      max_val <- max(combined_data$obs_fitted, combined_data$mean_fitted, na.rm = TRUE)
      
      make_plot <- function(data, val_col, title_suffix) {
        joined <- us_map_data %>%
          left_join(data, by = c("abbr" = "state"))
        
        plot_usmap(data = joined, values = val_col) +
          scale_fill_gradientn(
            name = "Death Rate",
            colours = c("blue", "yellow"),
            limits = c(min_val, max_val),
            na.value = "grey70"
          ) +
          geom_sf_text(aes(label = abbr), size = 1.6, color = "black", check_overlap = TRUE) +
          labs(title = title_suffix) +
          theme_void()
      }
      
      data_1973 <- combined_data %>% filter(year == 1973)
      data_2022 <- combined_data %>% filter(year == 2022)
      
      p1 <- make_plot(data_1973, "obs_fitted", "1973 Raw")
      p2 <- make_plot(data_1973, "mean_fitted", "1973 Fitted")
      p3 <- make_plot(data_2022, "obs_fitted", "2022 Raw")
      p4 <- make_plot(data_2022, "mean_fitted", "2022 Fitted")
      
      top_row <- (p1 | p2) + plot_layout(guides = "collect")
      bottom_row <- (p3 | p4) + plot_layout(guides = "collect")
      
      final_plot <- (top_row / bottom_row) +
        plot_annotation(
          title = paste(sex, r, "- Age", a),
          caption = "Note: Rates are per 100,000",
          theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
        )
      
      print(final_plot)
    }
  }
}

