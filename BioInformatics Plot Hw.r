######################################
# 10/18/2016
# Multiple Linear Regression (MLR) for 1-Mer vs 1-Mer vs Shape
# BISC 481
######################################

## Install and initialize packages
install.packages("ggplot2")
install.packages("grid")
library(ggplot2)
library(grid)

## Theme
my.theme <- theme(
  plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm"),
  axis.text = element_text(colour="black", size=12),
  axis.title.x = element_text(colour="black", size=12),
  axis.title.y = element_text(colour="black", size=12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(colour ="black"),
  axis.ticks = element_line(colour = "black")
)
## Data
## For 1-Mer
# R^2 for Mad= .7752117
# R^2 for Max= .7856751
# R^2 for MyC= .7780480
##1-Mer+ 1-shape
## R^2 for Mad= .863343
## R^2 for Max= .8644034
## R^2 for MyC= .854405

## Data preparation
"1-Mer"=data1 <- c(0.78, 0.79, 0.78)
"1-Mer+Shape"=data2 <- c(0.87, 0.86, 0.85)

## Ploting
ggplot() +
  geom_point(aes(x = data1, y = data2), color = "red", size=1) +
  geom_abline(slope=1) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  my.theme  

