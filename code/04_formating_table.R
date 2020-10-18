# Formatting AIC table
library(dplyr)
library(knitr)
library(kableExtra)
library(stringr)

# Reading output tables
aic.files <- list.files(path = "results",
                    pattern = "aic_table", full.names = TRUE)
aic.list <- lapply(aic.files, read.csv, row.names = 1)

# Fixing model names
model.names <- c("Trait-mediated Selection & Drift",
                 "Selection & Drift",
                 "Drift",
                 "Selection",
                 "Trait-mediated Selection",
                 "Idiosyncratic")
names(model.names) <- rownames(aic.list[[1]])

aic.list <- lapply(aic.list, function(x) {
  x$models = str_replace_all(rownames(x), model.names)
  return(x)
})

names(aic.list) <- c("Core", "Occasional")

# Binding and selecting colums to add in the manuscript
aic.df <- bind_rows(aic.list, .id = "Group") %>%
  select(Group, models, dAICc,  AICc, df, weight)

kable(aic.df, "latex", digits = 2) %>%
  collapse_rows(columns = 1)

