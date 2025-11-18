---
title: "m6APrediction: An R Package for Predicting m6A RNA Modification Sites"
author: "Xiaoyu Zhang"
date: "2025-11-18"
output:
  html_document:
    toc: true        
    toc_depth: true     
    number_sections: true
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Package Overview

`m6APrediction` is an R package designed to predict N6-methyladenosine (m6A) RNA modification sites using a Random Forest model. It supports both batch prediction and single-site prediction based on extracted sequence and structural features.

This tool was developed as part of a practical coursework project in BIO215 and can be integrated into Shiny-based bioinformatics web applications.

# Installation

To install the package:

```{r}

# Install devtools if not already installed
install.packages("devtools")

# Install m6APrediction from GitHub
devtools::install_github("zxy2360145-commits/m6APrediction")

# Load the package
library(m6APrediction)

```

# Usage

## Use Batch Prediction

The `prediction_multiple()` function takes a data frame of input features and predicts the m6A methylation status for each site.

```{r}
model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
input_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
prediction_multiple(model, input_df, positive_threshold = 0.6)
```

## Use Single Prediction

The prediction_single() function allows users to input site-specific features to predict whether a single RNA site is methylated.

```{r}
model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
prediction_single(
   ml_fit = model,
   gc_content = 0.6,
   RNA_type = "mRNA",
   RNA_region = "CDS",
   exon_length = 120,
   distance_to_junction = 5,
   evolutionary_conservation = 0.8,
   DNA_5mer = "ATCGAT",
   positive_threshold = 0.5
 )
```

# Model Performance

The ROC and PRC curves below show the model’s performance based on validation data (from Practical 4):

```{r}

knitr::include_graphics("m6APrediction/man/figure/ROC.png")
knitr::include_graphics("m6APrediction/man/figure/PRC.png")

```

# Project Information
	1. GitHub Repository: https://github.com/zxy2360145-commits/m6APrediction
	2. Bug Reports: https://github.com/zxy2360145-commits/m6APrediction/issues
	3. License: MIT

```{r}

sessionInfo()

```
