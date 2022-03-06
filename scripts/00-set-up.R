# PROJECT-WIDE SET-UP

# *************************************************************************************************

library(tidyverse)
library(lubridate)
library(tictoc)
library(furrr)
library(stars)
library(units)

sf_use_s2(F)

options(future.fork.enable = T)
options(future.globals.maxSize= 10000*1024^2)

# system("sudo sysctl -w vm.dirty_bytes=50331648")
# system("sudo sysctl -w vm.dirty_background_bytes=16777216")

# *************************************************************************************************

# MOUNT BUCKETS AND DRIVE

# risk team bucket
# dir_bucket_risk <- "/home/cdobler/bucket_risk/"

# my bucket
# dir_bucket_mine <- "/home/cdobler/bucket_mine/"

# persistent disk
# dir_disk <- "/home/cdobler/pers_disk"

# *************************************************************************************************

# OTHER DIRS
# dir_data <- "/home/cdobler/bucket_mine/misc_data/"
# dir_output <- "/home/cdobler/bucket_mine/results/global_scpdsi_ww/"
