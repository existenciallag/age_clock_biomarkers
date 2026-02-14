# 1. Instalar BioAge y cargar librerías
install.packages("devtools")
library(devtools)



install.packages("devtools", repos = "http://cran.us.r-project.org")
devtools::install_github("dayoonkwon/BioAge")
library(BioAge)
library(dplyr)
library(ggplot2)
