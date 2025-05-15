# All-cause mortality in the United States from 1973-2022

## Project description

This project utilizes Bayesian hierarchical modeling to the NCHS Multiple Cause of Death dataset to better understand changes in death rates across the United States from 1973 to 2022, with a focus on specific subgroups such as race and age.

## 1. Data

## 2.Code 

## 3. Figures

## Directory structure

```md
├── 01_data
│   └── summary_1970_2022.rds
├── 02_code
│   ├── 2a_data_exploration
│   │   └── EDA.Rmd
│   ├── 2b_model
│   │   └── INLA_model_script.R
│   ├── 2c_model_plotting
│   │   ├── raw_fitted_plotting
│   │   │   ├── national.Rmd
│   │   │   ├── age_specific.Rmd
│   │   │   └── race_specific.Rmd
│   │   └── smoothed_same_scale
│   │   │   └── example.R
├── 03_figures
└── README.md
```

## Data Availability 
Mortality data used in this analysis is found here: https://www.nber.org/research/data/mortality-data-vital-statistics-nchs-multiple-cause-death-data
