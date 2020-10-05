# Formatting AIC table
library(dplyr)
library(knitr)
library(kableExtra)

aic.files <- list.files(path = "results",
                    pattern = "aic_table", full.names = TRUE)
aic.list <- lapply(aic.files, read.csv)
names(aic.list) <- c("Core", "Occasional")

aic.df <- bind_rows(aic.list, .id = "Group") %>%
  select(-logLik, -dLogLik)

table <- kable(aic.df, "latex", digits = 2)

