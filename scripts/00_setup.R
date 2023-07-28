
# PROJECT-WIDE SET-UP

library(tidyverse)
library(lubridate)
library(tictoc)
library(furrr)
library(stars)
library(units)

sf_use_s2(F)

options(future.fork.enable = T)