library(ggplot2)
library(dplyr)
#library(colortools)

# Code to read csv simulation files
path.res <- "results/simulations"
files.prop <- list.files(path.res, "*._prop", full.names = TRUE)
files.effects <- list.files(path.res, "*._effects", full.names = TRUE)

df.prop <- lapply(files.prop, read.csv) %>%
  bind_rows(.id = "type")

df.effects <- lapply(files.effects, read.csv) %>%
  bind_rows(.id = "type")

# First set of plots: one round of simulations #################################
t1 <- read.csv("results/test1_tutorial.csv")
t2 <- read.csv("results/test2_tutorial.csv")

all.t1 <- tapply(t1$Value,
                 list(t1$Community),
                 sum, na.rm = TRUE)

all.t2 <- tapply(t2$Value,
                 list(t2$Community),
                 sum, na.rm = TRUE)

table(t1$Community)
table(t2$Community)

t1$Total <- rep(all.t1, each = 6)
t2$Total <- rep(all.t2, each = 6)

t1$Proportion <- t1$Value/t1$Total
t2$Proportion <- t2$Value/t2$Total

neutral <- "#091346"
niche <- "#801411"
nineu <- "#D6A838" #"#e1b727"
idiosyncratic <- "#B4B1AC"

paleta <- c(niche, neutral, nineu, grey)

# other stuff w/ colors
# library(colortools)
# niche.cols <- analogous(niche)
# neutral.cols <- analogous(neutral)
# wheel(niche, num = 10)
# wheel(neutral, num = 10)
# wheel(nineu, num = 10)

p1 <- ggplot(t1, aes(x = Community,
                     y = Proportion,
                     fill = Effect)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic()


p1
p2 <- ggplot(t2, aes(x = Community,
                     y = Proportion,
                     fill = Effect)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic()

p2

grid.arrange(p1, p2, nrow = 1)

png("figures/test1.png")
p1
dev.off()

png("figures/test2.png")
p2
dev.off()

p1
p2
