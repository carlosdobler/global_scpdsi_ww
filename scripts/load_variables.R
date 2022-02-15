
# SETUP --------------------------------------------------------------------------------------------

seq(as_date("1979-01-01"), as_date("2020-12-01"), by = "months") -> date_vector
days_in_month(date_vector) %>% unname() -> days_mth

# mex
# c(st_point(c(-120, 33)),
#   st_point(c(-80, 13))) %>%
#   st_bbox() %>%
#   st_set_crs(4326) -> lim

# central europe
c(st_point(c(0.5, 58)),
  st_point(c(20, 43))) %>% 
  st_bbox() %>% 
  st_set_crs(4326) -> lim

# "~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip.nc" %>% 
#   read_stars(proxy = T) %>%
#   slice(time, 1) -> s_proxy

"~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip.nc" %>% 
  read_ncdf(ncsub = cbind(start = c(1,1,1), count = c(NA,NA,1))) %>%
  adrop() -> s_proxy

# lon
s_proxy %>% 
  st_get_dimension_values(1, where = "start") -> s_proxy_lon 

which.min(abs(s_proxy_lon - lim[1])) -> lon_start
which.min(abs(s_proxy_lon - lim[3])) - lon_start -> lon_count

# lat
s_proxy %>% 
  st_get_dimension_values(2, where = "start") -> s_proxy_lat

which.min(abs(s_proxy_lat - lim[4])) -> lat_start
which.min(abs(s_proxy_lat - lim[2])) - lat_start -> lat_count

sub <- cbind(start = c(lon_start, lat_start, 1),
             count = c(lon_count, lat_count, NA))

# PRECIPITATION ------------------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip.nc" %>%
  read_ncdf(ncsub = sub) -> s

st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]

s %>%
  mutate(tp = tp %>% units::set_units("mm"),
         tp = tp %>% units::drop_units()) %>% 
  setNames("X") -> var_pr

var_pr %>% 
  st_apply(c(1,2), function(x){
    x * days_mth
  },
  .fname = "time") %>% 
  
  st_set_dimensions("time", values = date_vector) %>% 
  aperm(c(2,3,1)) -> var_pr


# RADIATION ----------------------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_monthly_mean_daily_radiation.nc" %>% 
  read_ncdf(ncsub = sub) -> s

st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]

# external / extraterrestrial
s %>%   
  select(tisr) %>% # toa incident solar radiation
  mutate(tisr = units::set_units(tisr, MJ/m^2),
         tisr = units::drop_units(tisr)) %>% 
  st_set_dimensions("time", values = date_vector) %>% 
  setNames("X") -> var_ra

# internal / incoming
s %>% 
  select(ssrd) %>% # surface solar radiation downwards
  mutate(ssrd = units::set_units(ssrd, MJ/m^2),
         ssrd = units::drop_units(ssrd)) %>% 
  st_set_dimensions("time", values = date_vector) %>% 
  setNames("X") -> var_rs


# POTENTIAL EVAPORATION ---------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_monthly_mean_daily_potevap.nc" %>% 
  read_ncdf(ncsub = sub) -> s

st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]  

s %>%   
  mutate(pev = pev %>% units::set_units("mm"),
         pev = pev %>% units::drop_units()) %>% 
  setNames("X") -> var_potevap

var_potevap %>% 
  st_apply(c("longitude", "latitude"), function(x){
    x * days_mth
  },
  .fname = "time") %>% 
  
  st_set_dimensions("time", values = date_vector) %>% 
  aperm(c(2,3,1)) -> var_potevap


# MIN MAX TEMPERATURE ------------------------------------------------------------------------------
# (based on monthly by hr averages)

# file.copy("~/bucket_mine/era/monthly_by_hr_day/2mtemp.nc", "~/pers_disk/")
# 
# seq(1, 12384-(24*12), by = 24) -> day_breaks
# 
# imap(day_breaks, function(db, i){
# 
#   tic(str_glue("Done {i} / 504"))
# 
#   # "~/bucket_mine/era/monthly_by_hr_day/2mtemp.nc" %>%
#   "~/pers_disk/2mtemp.nc" %>%
#     read_ncdf(ncsub = cbind(start = c(1,1,db),
#                             count = c(NA, NA, 24))) %>%
#     suppressMessages() -> s_mth
# 
#   s_mth %>%
#     setNames("X") %>%
#     mutate(X = X %>% units::set_units("degC"),
#            X = X %>% units::drop_units()) -> s_mth
# 
#   s_mth %>%
#     st_apply(c(1,2), function(x){
#       c(max = max(x),
#         min = min(x))
#     },
#     .fname = "func") %>%
#     split("func") -> s_mth_range
# 
#   toc()
# 
#   return(s_mth_range)
# 
# }) -> s_mth_range
# 
# s_mth_range %>% 
#   {do.call(c, c(., along = "time"))} -> s_mth_range
# 
# s_mth_range %>% 
#   st_set_dimensions("time", values = date_vector) -> s_mth_range
# 
# saveRDS(s_mth_range, "~/bucket_mine/era/monthly/era5_monthly_mean_minmaxtemp.rds")

"~/bucket_mine/era/monthly/era5_monthly_mean_minmaxtemp.rds" %>% 
  readRDS() -> s
  
s %>% 
  st_crop(var_pr) -> var_t

# max
var_t %>% 
  select(max) %>% 
  setNames("X") -> var_tmax

# min
var_t %>% 
  select(min) %>% 
  setNames("X") -> var_tmin


# MIN MAX TEMPERATURE 2 ----------------------------------------------------------------------------
# (real one)

# # max
# map(seq_len(42), function(y){
#   
#   list.files("~/bucket_mine/era/monthly/2m_temperature_max/", full.names = T) %>% 
#     .[y] %>% 
#     readRDS() %>% 
#     st_crop(lim)
#   
# }) %>% 
#   do.call(c, .) -> var_tmax
# 
# # min
# map(seq_len(42), function(y){
#   
#   list.files("~/bucket_mine/era/monthly/2m_temperature_min/", full.names = T) %>% 
#     .[y] %>% 
#     readRDS() %>% 
#     st_crop(lim)
#   
# }) %>% 
#   do.call(c, .) -> var_tmin


# AVERAGE TEMPERATURE ------------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_monthly_mean_meantemp.nc" %>% 
  read_ncdf(ncsub = sub) -> s

st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]

s %>%   
  mutate(t2m = units::set_units(t2m, degC),
         t2m = units::drop_units(t2m)) %>%
  st_set_dimensions("time", value = date_vector) %>% 
  setNames("X") -> var_tave


# LATITUDE -----------------------------------------------------------------------------------------

var_tave %>%
  slice(time, 1) %>%
  st_dim_to_attr(2) -> var_lat


# WIND SPEED ---------------------------------------------------------------------------------------
# (need to convert to 2m)

"~/bucket_mine/era/monthly/era5_monthly_mean_10mwindspeed.nc" %>% 
  read_ncdf(ncsub = sub) -> s 
  
st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]  

s %>%
  mutate(si10 = si10 %>% units::drop_units()) %>% 
  st_set_dimensions("time", value = date_vector) %>% 
  setNames("X") -> var_wind


# DEWPOINT TEMP ------------------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_monthly_mean_dewpointtemp.nc" %>% 
  read_ncdf(ncsub = sub) -> s

st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]
  
s %>% 
  mutate(d2m = d2m %>% units::set_units(degC),
         d2m = d2m %>% units::drop_units()) %>% 
  st_set_dimensions("time", value = date_vector) %>% 
  setNames("X") -> var_dewpoint
  

# SURFACE PRESSURE ---------------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_monthly_mean_surfpressure.nc" %>% 
  read_ncdf(ncsub = sub) -> s

st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]  

s %>%
  mutate(sp = sp %>% units::set_units(kPa),
         sp = sp %>% units::drop_units()) %>% 
  st_set_dimensions("time", value = date_vector) %>% 
  setNames("X") -> var_pressure


# ELEVATION ----------------------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_geopotential.nc" %>% 
  read_ncdf(ncsub = sub) -> s
  
st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]
  
s %>% 
  slice(time, 1) %>%
  
  mutate(z = z %>% units::drop_units(),
         z = z/9.80665,
         z = ifelse(z < 0, 0, z)) %>% 
  setNames("X") -> var_z


# VAN DER SCHRIER ----------------------------------------------------------------------------------

here::here("data", "scPDSI.cru_ts4.05early1.1901.2020.cal_1901_20.bams.2021.GLOBAL.1901.2020.nc") %>% 
  read_ncdf() %>%
  st_warp(var_pr) %>% 
  filter(year(time) >= 1979) %>% 
  st_set_dimensions("time", values = date_vector) -> vds

vds %>% 
  setNames("m") %>% 
  mutate(m = ifelse(is.na(m), NA, 1)) -> mask
  
  



