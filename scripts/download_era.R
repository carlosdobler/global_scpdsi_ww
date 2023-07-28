library(reticulate)
use_condaenv("my_first_env")

source_python("scripts/era_download.py")



# cdsapi <- import("cdsapi")
# cds <- cdsapi$Client()
# 
# seq(1989, 1997) -> x
# split(x, ceiling(seq_along(x)/(length(x)/2))) -> x_yrs
# unname(x_yrs) -> x_yrs
# x_yrs[[3]] <- c(2019, 2020)
# 
# 
# for(yrs in x_yrs){
#   
#   print(str_glue("Downloading {first(yrs)} - {last(yrs)}"))
#   
#   cds$retrieve('reanalysis-era5-single-levels',
#                
#                list(product_type = "reanalysis",
#                     variable = "total_precipitation",
# 
#                     year = as.character(yrs),
# 
#                     month = seq_len(12) %>% str_pad(2, "left", "0"),
#                     day = seq_len(31) %>% str_pad(2, "left", "0"),
#                     time = seq(0, 23) %>% str_pad(2, "left", pad = "0") %>% str_c(":00"),
# 
#                     format = "netcdf"),
#                
#                str_glue("/home/cdobler/bucket_mine/era/hourly/total_precipitation/era5_totprecip_hourly_{yrs}.nc"))
#   
# }
# 
# cds$retrieve('reanalysis-era5-single-levels-monthly-means',
#              
#              list(product_type = "monthly_averaged_reanalysis",
#                   variable = "total_precipitation",
#                   
#                   year = seq(1979,2020) %>% as.character(),
#                   
#                   month = seq_len(12) %>% str_pad(2, "left", "0"),
#                   time = '00:00',
#                   
#                   format = "netcdf"),
#              
#              str_glue("/home/cdobler/bucket_mine/era/monthly/era5_totprecip_monthlymeans_rawfile.nc"))
# 
# 
# cds$retrieve('reanalysis-era5-single-levels-monthly-means',
#              
#              list(product_type = "monthly_averaged_reanalysis",
#                   variable = "potential_evaporation",
#                   
#                   year = seq(1979,2020) %>% as.character(),
#                   
#                   month = seq_len(12) %>% str_pad(2, "left", "0"),
#                   time = '00:00',
#                   
#                   format = "netcdf"),
#              
#              str_glue("/home/cdobler/bucket_mine/era/monthly/era5_potevap_dailymeans_monthlyres_rawfile.nc"))
# 
# 
# cds$retrieve('reanalysis-era5-single-levels-monthly-means',
#              
#              list(product_type = "monthly_averaged_reanalysis_by_hour_of_day",
#                   variable = "2m_temperature",
#                   
#                   year = 1980, #seq(1979,2020) %>% as.character(),
#                   
#                   month = seq_len(12) %>% str_pad(2, "left", "0"),
#                   time = seq(0, 23) %>% str_pad(2, "left", pad = "0") %>% str_c(":00"),
#                   
#                   format = "netcdf"),
#              
#              "/home/cdobler/bucket_mine/era/monthly_by_hr_day/test.nc")
# 
# # *********
# 
# cds$retrieve('reanalysis-era5-single-levels-monthly-means',
#              
#              list(product_type = "monthly_averaged_reanalysis",
#                   variable = "2m_temperature",
#                   
#                   year = seq(1979,2020) %>% as.character(),
#                   
#                   month = seq_len(12) %>% str_pad(2, "left", "0"),
#                   time = '00:00',
#                   
#                   format = "netcdf"),
#              
#              "/home/cdobler/bucket_mine/era/monthly/era5_2mtemp_dailymeans_monthlyres_rawfile.nc")
# 
# # ************
# 
# cds$retrieve('reanalysis-era5-single-levels-monthly-means',
#              
#              list(product_type = 'monthly_averaged_reanalysis',
#                   variable = c('surface_solar_radiation_downwards', 'toa_incident_solar_radiation'),
#                   
#                   year = seq(1979,2020) %>% as.character(),
#                   
#                   month = seq_len(12) %>% str_pad(2, "left", "0"),
#                   time = '00:00',
#                   
#                   format = "netcdf"),
#              
#              "/home/cdobler/bucket_mine/era/monthly/era5_radiation_dailymeans_monthlyres_rawfile.nc")
# 
# # ***********
# 
# cds$retrieve('reanalysis-era5-single-levels-monthly-means',
#              
#              list(product_type = "monthly_averaged_reanalysis",
#                   variable = "10m_wind_speed",
#                   
#                   year = seq(1979,2020) %>% as.character(),
#                   
#                   month = seq_len(12) %>% str_pad(2, "left", "0"),
#                   time = '00:00',
#                   
#                   format = "netcdf"),
#              
#              "/home/cdobler/bucket_mine/era/monthly/era5_monthly_mean_10mwindspeed.nc"
#              )
# 
# # ***********
# 
# cds$retrieve('reanalysis-era5-single-levels-monthly-means',
#              
#              list(product_type = "monthly_averaged_reanalysis",
#                   variable = "2m_dewpoint_temperature",
#                   
#                   year = seq(1979,2020) %>% as.character(),
#                   
#                   month = seq_len(12) %>% str_pad(2, "left", "0"),
#                   time = '00:00',
#                   
#                   format = "netcdf"),
#              
#              "/home/cdobler/bucket_mine/era/monthly/era5_monthly_mean_dewpointtemp.nc"
# )
# 
# # *************
# 
# cds$retrieve('reanalysis-era5-single-levels-monthly-means',
#              
#              list(product_type = "monthly_averaged_reanalysis",
#                   variable = "surface_pressure",
#                   
#                   year = seq(1979,2020) %>% as.character(),
#                   
#                   month = seq_len(12) %>% str_pad(2, "left", "0"),
#                   time = '00:00',
#                   
#                   format = "netcdf"),
#              
#              "/home/cdobler/bucket_mine/era/monthly/era5_monthly_mean_surfpressure.nc"
# )
# 
# 
# 
# 
# 
# 
# # *************
# 
# nc_dir <- str_glue("{dir_bucket_mine}/era/hourly/total_precipitation")
# read_ncdf(str_c(nc_dir, "/test.nc")) -> t
# t %>% st_get_dimension_values("time")
# 
# nc_meta(str_c())
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
