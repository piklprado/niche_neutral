library(ggplot2)
library(dplyr)

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

cores <- c(ni = rgb(202, 0, 32, maxColorValue = 255),
           neu = rgb(4, 4, 160, maxColorValue = 255),
           id = rgb(180, 177, 172,  maxColorValue = 255),
           nineu = rgb(245, 219 ,58, maxColorValue = 255))

cores



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

png("figures/test1.png")
p1
dev.off()

png("figures/test2.png")
p2
dev.off()

p1
p2
cores
