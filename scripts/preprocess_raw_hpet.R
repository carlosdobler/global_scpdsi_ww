
# SCRIPT TO TRANSFORM RAW hPET FILES (0.1 RES, DAILY, -180:180) TO
# MONTHLY ERA5 FORMAT (0.25 RES, MONTHLY SUMS, 0:360) 


library(tidyverse)



# specify grid
# based on <cdo griddes raw_file>

c("gridtype = lonlat",
  "xsize = 3600",
  "ysize = 1801",
  "xname = longitude",
  "yname = latitude",
  "xfirst = -180",
  "xinc = 0.1",
  "yfirst = 90",
  "yinc = -0.1") %>%
  
  write_lines("grid.txt")


# loop through years

walk(1981:2021, function(yr){
  
  print(yr)
  
  infile <- "/mnt/bucket_mine/era/daily/hPET_raw/{yr}_daily_pet.nc" %>% 
    str_glue()
  
  outfile <- "/mnt/bucket_mine/era/monthly/total-pet-v-hPET/era5_mon_total-pet_{yr}.nc" %>% 
    str_glue()
  
  template <- "/mnt/bucket_mine/era/monthly/mean-daily-precip/era5_mon_mean-daily-precip_1959.nc"
  
  # setgrid > aggregate to monthly > aggregate to 0.25 deg (nearest neighbor)
  
  "cdo remapnn,{template} -monsum -setgrid,grid.txt {infile} {outfile}" %>% 
    str_glue() %>% 
    system()
  
})



