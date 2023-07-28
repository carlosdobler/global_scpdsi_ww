
 #  _    _         _  
 # |_)  |_  |\/|  / \ 
 # | \  |_  |  |  \_/ 
 #                    

# DOWNLOAD FILES

# chmod 744 <wget_script.sh>
# bash <wget_script.sh> -H


# **********************************************************

# SET LIMITS

"~/bucket_mine/remo/monthly/hurs_EUR-22_MOHC-HadGEM2-ES_historical_r1i1p1_GERICS-REMO2015_v1_mon_200101-200512_regrid.nc" %>% 
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

which.min(abs(s_proxy_lat - lim[2])) -> lat_start
which.min(abs(s_proxy_lat - lim[4])) - lat_start -> lat_count

sub <- cbind(start = c(lon_start, lat_start, 1),
             count = c(lon_count, lat_count, NA))


# ************************************************************

# IMPORT

variable <- c("tasmax", "tasmin", "hurs", "rsds", "sfcWind", "pr")

map(c("Had", "MPI", "Nor") %>% set_names(), function(model){
  
  case_when(model == "Had" ~ 2014,
            model == "MPI" ~ 2006,
            model == "Nor" ~ 2019) -> mid_yr
  
  date_vector_1deg <- seq(as_date(str_c(mid_yr-10,"-01-01")), 
                          as_date(str_c(mid_yr+10,"-12-01")),
                          by = "months")
  
  map(variable, function(var){
    
    print(str_glue("Processing model {model} - variable {var}"))
    
    list.files("~/bucket_mine/remo/monthly/", full.names = T) %>% 
      .[str_detect(., model)] %>% 
      .[str_detect(., "regrid")] %>% 
      .[str_detect(., var)] %>%
      
      map(function(f){
        
        f %>% 
          read_ncdf(ncsub = sub) %>% 
          suppressMessages() -> s
        
        st_dimensions(s)$lon$offset <- s_proxy_lon[lon_start]
        st_dimensions(s)$lat$offset <- s_proxy_lat[lat_start]
        
        s %>%
          filter(year(time) >= mid_yr-10,
                 year(time) <= mid_yr+10)
        
      }) -> s
    
    s %>% 
      do.call(c, .)
    
  }) -> s_vars
  
  s_vars %>% 
    do.call(c, .) %>%
    st_set_dimensions("time", values = date_vector_1deg)
  
}) -> r_models_vars

rm(s_proxy, sub, lat_count, lat_start, lon_count, 
   lon_start, s_proxy_lat, s_proxy_lon, variable)


# *********************************************************

# LAND MASK

# st_read("~/bucket_mine/misc_data/ne_50m_land/ne_50m_land.shp") %>%
#   mutate(a = 1) %>% 
#   select(a) %>% 
#   # st_crop(st_bbox(s_models_vars[[1]])) %>% 
#   st_rasterize(st_as_stars(st_bbox(), dx = 0.22, dy = 0.22, values = NA)) -> r_land_rast
# 
# st_warp(r_land_rast %>% 
#           st_set_dimensions(which = c(1,2),
#                             names = c("lon", "lat")),
#         r_models_vars[[1]] %>% 
#           slice(time, 1)) -> r_land_rast

# ********************************************************

# MODIFY

func_dewpoint <- function(hurs, tasmean){
  
  hurs <- set_units(hurs, NULL)
  tasmean <- set_units(tasmean, NULL)
  
  243.04*(log(hurs/100)+((17.625*tasmean)/(243.04+tasmean)))/
    (17.625-log(hurs/100)-((17.625*tasmean)/(243.04+tasmean)))
  
}

r_models_vars %>% 
  imap(function(s, i){
    
    print(str_glue("Processing {i}"))
    
    s %>%
      select(-pr) %>% 
      mutate(tasmax = tasmax %>% set_units(degC) %>% set_units(NULL),
             tasmin = tasmin %>% set_units(degC) %>% set_units(NULL),
             rsds = rsds %>% set_units(MJ/d/m^2) %>% set_units(NULL),
             sfcWind = sfcWind %>% set_units(NULL),
             tasmean = (tasmax + tasmin)/2,
             dewpoint = func_dewpoint(hurs, tasmean)) %>% 
      select(-hurs, -tasmean) -> s_no_pr
    
    if(i == "Had"){
      
      s %>%
        select(pr) %>% 
        mutate(pr = pr %>% 
                 set_units(kg/m^2/d) %>% 
                 set_units(NULL) %>% 
                 {.*30}) -> s_pr
      
    } else {
      
      s %>% 
        st_get_dimension_values("time") -> dates
      
      dates %>%   
        days_in_month() %>% 
        unname() -> days_mth
      
      s %>% 
        select(pr) %>% 
        mutate(pr = pr %>% 
                 set_units(kg/m^2/d) %>% 
                 set_units(NULL)) %>%
        st_apply(c(1,2),
                 function(x) x * days_mth,
                 .fname = "time") %>% 
        st_set_dimensions("time", values = dates, point = FALSE) %>% 
        aperm(c(2,3,1)) -> s_pr
      
    }
      
    c(s_no_pr, s_pr) %>% 
      st_set_dimensions("time", point = FALSE)
      
  }) -> r_models_vars

# *****************************************

# Z

"~/bucket_mine/era/monthly/era5_geopotential_shifted.nc" %>% 
  read_ncdf() %>% 
  slice(time, 1) -> r_var_z

r_var_z %>%
  mutate(z = z %>% set_units(NULL) %>% {./9.80665},
         z = ifelse(z < 0, 0, z)) %>% 
  setNames("z") -> r_var_z

r_var_z %>% 
  st_warp(r_models_vars[[1]] %>% slice(time, 1)) -> r_var_z

# list(var_z_remo) %>% rep(252) %>% {do.call(c, c(., along = "time"))} -> var_z_remo_252

# *****************************************

# LAT

r_var_z %>% 
  st_dim_to_attr(2) -> r_var_lat

# *****************************************

# seq(as_date(15), as_date(364-15), by = "1 month") %>% 
#   as.integer() %>% 
#   map(function(J){
#     
#     matrix(J, nrow = 507, ncol = 227) -> m
#     names(dim(m)) <- c("x", "y")
#     st_as_stars(m)
#     
#   }) -> m_one_year
# 
# m_one_year %>% rep(21) %>% {do.call(c, c(., along = "time"))} %>% 
#   setNames("J") -> m_all_yrs
# 
# st_dimensions(m_all_yrs) <- st_dimensions(var_lat_remo_252)


# *****************************************

# PET

source(here::here("scripts/penman_mine.R"))

r_models_vars %>%
  imap(function(s, i){
    
    print(str_glue("Processing {i}"))
    
    abind(
      s[[1]], # tmax
      s[[2]], # tmin
      s[[3]], # rsds
      s[[4]], # sfcWind
      s[[5]], # dewpoint
      r_var_z[[1]], # z
      r_var_lat[[1]], # lat
      along = 3) -> s_array
    
    names(dim(s_array)) <- c("lon", "lat", "time")
      
    s_array %>% 
      st_as_stars() -> s_array
    
    c(dim(s[[1]])[3],
      dim(s[[2]])[3],
      dim(s[[3]])[3],
      dim(s[[4]])[3],
      dim(s[[5]])[3],
      1,
      1) %>% 
      unname() %>% 
      cumsum() -> index
    
    tic("applying")
    s_array %>% 
      st_apply(c(1,2), function(x){
        
        penman_mine(Tmax = x[1:index[1]],
                    Tmin = x[(index[1]+1):index[2]],
                    Rs = x[(index[2]+1):index[3]],
                    uz = x[(index[3]+1):index[4]],
                    Tdew = x[(index[4]+1):index[5]],
                    z = x[index[6]],
                    lat = x[index[7]],
                    model_i = i)
        
      }, 
      .fname = "time") -> s_pet
    toc()
    
    s_pet %>% 
      aperm(c(2,3,1)) %>% 
      setNames("pet") -> s_pet
    
    st_dimensions(s_pet) <- st_dimensions(s) 
      
    return(s_pet)           
    
  }) -> r_models_pet

# *************************************

# REMO PET MSC ----

r_models_pet %>% 
  map(function(s){
    
    map(seq_len(12), function(mth){
      
      s %>% 
        filter(month(time) == mth) %>% 
        st_apply(c(1,2),
                 mean,
                 na.rm = T,
                 .fname = month.abb[mth])
      
    }) %>% 
      do.call(c, .) %>% 
      merge(name = "month") %>% 
      setNames("pet")
    
  }) -> r_models_pet_msc

r_models_pet_msc %>% 
  {do.call(c, c(., along = "model"))} %>% 
  st_apply(c(1,2,3), mean, na.rm = T) -> r_models_pet_msc[["ENS"]]
  
    
# *******************************************************

# REMO PRECIP MSC ----

r_models_vars %>% 
  map(select, "pr") %>% 
  map(function(s){
    
    map(seq_len(12), function(mth){
      
      s %>% 
        filter(month(time) == mth) %>% 
        st_apply(c(1,2),
                 mean,
                 na.rm = T,
                 .fname = month.abb[mth])
      
    }) %>% 
      do.call(c, .) %>% 
      merge(name = "month") %>% 
      setNames("pr")
    
  }) -> r_models_pr_msc

r_models_pr_msc %>% 
  {do.call(c, c(., along = "model"))} %>% 
  st_apply(c(1,2,3), mean, na.rm = T) -> r_models_pr_msc[["ENS"]]


# ******************************************************

# REMO PET ANN ----

r_models_pet_msc %>% 
  {do.call(c, c(., along = "model"))} %>% 
  st_apply(c(1,2,4), mean, na.rm = T) -> r_models_pet_ann

# REMO PRECIP ANN ----

r_models_pr_msc %>% 
  {do.call(c, c(., along = "model"))} %>% 
  st_apply(c(1,2,4), mean, na.rm = T) -> r_models_pr_ann


# ******************************************************

# REMO SPEI ----

# plan(multicore, workers = 6)
# 
# library(SPEI)
# 
# c("Had", "MPI", "Nor") %>% set_names() %>%
#   map(function(model){
# 
#     print(str_glue("Processing {model}"))
# 
#     r_models_pet %>%
#       pluck(model) -> s_pet
# 
#     r_models_vars %>%
#       pluck(model) %>%
#       select(pr) -> s_pr
# 
#     c(s_pet, s_pr) %>%
#       mutate(wb = pr - pet) %>%
#       select(wb) -> s_wb
# 
#     s_wb %>%
#       st_apply(c(1,2),
#                function(x){
# 
#                  spei(x, scale = 3, na.rm = T)$fitted %>%
#                    as.vector() %>%
#                    ifelse(is.infinite(.), NA, .)
# 
#                },
#                FUTURE = T,
#                .fname = "time") -> s_spei
# 
#     s_spei %>%
#       st_set_dimensions("time",
#                         values = st_get_dimension_values(s_pet, "time")) %>%
#       aperm(c(2,3,1)) %>%
#       setNames("spei")
# 
#   }) -> r_models_spei
# 
# saveRDS(r_models_spei, "output/r_models_spei3.rds")

readRDS(here::here("output/r_models_spei3.rds")) -> r_models_spei

# **************************************

# REMO SPEI MSC ----

r_models_spei %>% 
  map(function(s){
    
    map(seq_len(12), function(mth){
      
      s %>% 
        filter(month(time) == mth) %>% 
        st_apply(c(1,2),
                 function(x){
                   sum(x < -1.5, na.rm = T)
                 },
                 .fname = month.abb[mth])
      
    }) %>% 
      do.call(c, .) %>% 
      merge(name = "month") %>% 
      setNames("spei")
    
  }) -> r_models_spei_msc

r_models_spei_msc %>% 
  {do.call(c, c(., along = "model"))} %>% 
  st_apply(c(1,2,3), mean, na.rm = T, .fname = "spei") -> r_models_spei_msc[["ENS"]]

# **************************************

# REMO SPEI ANN ----

r_models_spei %>% 
  map(function(s){
    
    map2(seq(1, 252, 12),
         seq(12, 253, 12),
         function(st, en){
           
           s %>% 
             slice(time, st:en) %>% 
             st_apply(c(1,2),
                      function(x){
                        sum(x < -1.5, na.rm = T)
                      },
                      .fname = "ann_count")
           
         }) %>% 
      do.call(c, .) %>% 
      merge() %>% 
      st_apply(c(1,2), mean, na.rm = T, .fname = "spei")
    
    }) -> r_models_spei_ann

r_models_spei_ann %>% 
  {do.call(c, c(., along = "model"))} %>% 
  st_apply(c(1,2), mean, na.rm = T, .fname = "spei") -> r_models_spei_ann[["ENS"]]
      
r_models_spei_ann %>% 
  {do.call(c, c(., along = "model"))} -> r_models_spei_ann


# ******************************************************
# ******************************************************

#  ___  __       
# |__  |__)  /\  
# |___ |  \ /~~\ 
#                

# SET LIMITS

"~/bucket_mine/era/monthly/era5_geopotential_shifted.nc" %>% 
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

sub <- cbind(start = c(lon_start, lat_start, 253),
             count = c(lon_count, lat_count, NA))

# ****************************************************

# ERA IMPORT -----

# precip
"~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip_shifted.nc" %>%
  read_ncdf(ncsub = sub) -> s

st_dimensions(s)$longitude$offset <- s_proxy_lon[lon_start]
st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]

s %>%
  mutate(tp = tp %>% set_units("mm") %>% set_units(NULL)) -> e_var_pr

e_var_pr %>% 
  st_get_dimension_values(3) %>% 
  as_date() -> dates_1deg

dates_1deg %>%   
  days_in_month() %>% 
  unname() -> days_mth_1deg

e_var_pr %>% 
  st_apply(c(1,2), function(x){
    x * days_mth_1deg
  },
  .fname = "time") %>% 
  
  st_set_dimensions("time", values = dates_1deg) %>% 
  aperm(c(2,3,1)) -> e_var_pr


# radiation
"~/bucket_mine/era/monthly/era5_monthly_mean_daily_radiation_shifted.nc" %>%
  read_ncdf(ncsub = sub) -> s

st_dimensions(s)$longitude$offset <- s_proxy_lon[lon_start]
st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]

s %>% 
  select(ssrd) %>% # surface solar radiation downwards
  mutate(ssrd = set_units(ssrd, MJ/m^2) %>% set_units(NULL)) %>% 
  st_set_dimensions("time", values = dates_1deg) -> e_var_rs


# wind speed
"~/bucket_mine/era/monthly/era5_monthly_mean_10mwindspeed_shifted.nc" %>%
  read_ncdf(ncsub = sub) -> s

st_dimensions(s)$longitude$offset <- s_proxy_lon[lon_start]
st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]

s %>%
  mutate(si10 = si10 %>% set_units(NULL)) %>% 
  st_set_dimensions("time", value = dates_1deg) -> e_var_wind


# dewpoint
"~/bucket_mine/era/monthly/era5_monthly_mean_dewpointtemp_shifted.nc" %>%
  read_ncdf(ncsub = sub) -> s

st_dimensions(s)$longitude$offset <- s_proxy_lon[lon_start]
st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]

s %>% 
  mutate(d2m = d2m %>% set_units(degC) %>% set_units(NULL)) %>% 
  st_set_dimensions("time", value = dates_1deg) -> e_var_dewpoint


# elevation
"~/bucket_mine/era/monthly/era5_geopotential_shifted.nc" %>% 
  read_ncdf(ncsub = sub[1:2,] %>% rbind(1)) %>% adrop() -> s

st_dimensions(s)$longitude$offset <- s_proxy_lon[lon_start]
st_dimensions(s)$latitude$offset <- s_proxy_lat[lat_start]

s %>%
  mutate(z = z %>% set_units(NULL),
         z = z/9.80665,
         z = ifelse(z < 0, 0, z)) -> e_var_z


# tasmax
list.files("~/bucket_mine/era/monthly/2m_temperature_max/", full.names = T) %>% 
  .[str_detect(., seq(2000,2020) %>% as.character())] %>% 
  imap(function(y, i){
    
    print(str_glue("Processing {i} / 21"))
    
    y %>% 
      readRDS() %>% 
      st_warp(e_var_z)
    
  }) %>% 
  do.call(c, .) %>% 
  setNames("tasmax") -> e_var_tmax


# tasmin
list.files("~/bucket_mine/era/monthly/2m_temperature_min/", full.names = T) %>% 
  .[str_detect(., seq(2000,2020) %>% as.character())] %>% 
  imap(function(y, i){
    
    print(str_glue("Processing {i} / 21"))
    
    y %>% 
      readRDS() %>% 
      st_warp(e_var_z)
    
  }) %>% 
  do.call(c, .) %>% 
  setNames("tasmin") -> e_var_tmin


# latitude
e_var_z %>%
  st_dim_to_attr(2) -> e_var_lat

rm(s, s_proxy, sub, lat_count, lat_start, lon_count, 
   lon_start, s_proxy_lat, s_proxy_lon)

# **********************************************

# PET

abind(
  e_var_tmax[[1]], # tmax
  e_var_tmin[[1]], # tmin
  e_var_rs[[1]], # rsds
  e_var_wind[[1]], # sfcWind
  e_var_dewpoint[[1]], # dewpoint
  e_var_z[[1]], # z
  e_var_lat[[1]], # lat
  along = 3) -> s_array

names(dim(s_array)) <- c("longitude", "latitude", "time")

s_array %>% 
  st_as_stars() -> s_array

c(dim(e_var_tmax[[1]])[3],
  dim(e_var_tmin[[1]])[3],
  dim(e_var_rs[[1]])[3],
  dim(e_var_wind[[1]])[3],
  dim(e_var_dewpoint[[1]])[3],
  1,
  1) %>% 
  unname() %>% 
  cumsum() -> index

tic("applying")
s_array %>% 
  st_apply(c(1,2), function(x){
    
    penman_mine(Tmax = x[1:index[1]],
                Tmin = x[(index[1]+1):index[2]],
                Rs = x[(index[2]+1):index[3]],
                uz = x[(index[3]+1):index[4]],
                uz_2 = F,
                Tdew = x[(index[4]+1):index[5]],
                z = x[index[6]],
                lat = x[index[7]])
    
  }, 
  .fname = "time") -> e_var_pet
toc()

e_var_pet %>% 
  aperm(c(2,3,1)) %>% 
  setNames("pet") -> e_var_pet

st_dimensions(e_var_pet) <- st_dimensions(e_var_tmax) 

# *********************************

# ERA PET MSC ----

map(seq_len(12), function(mth){
  
  e_var_pet %>% 
    filter(month(time) == mth) %>% 
    st_apply(c(1,2),
             mean,
             na.rm = T,
             .fname = month.abb[mth])
  
}) %>% 
  do.call(c, .) %>% 
  merge(name = "month") %>% 
  setNames("pet") -> e_var_pet_msc

# *******************************************************

# PRECIP MSC

map(seq_len(12), function(mth){
  
  e_var_pr %>% 
    filter(month(time) == mth) %>% 
    st_apply(c(1,2),
             mean,
             na.rm = T,
             .fname = month.abb[mth])
  
}) %>% 
  do.call(c, .) %>% 
  merge(name = "month") %>% 
  setNames("pr") -> e_var_pr_msc

# ******************************************************

# PET ANN

e_var_pet_msc %>% 
  st_apply(c(1,2), mean, na.rm = T, .fname = "ERA") -> e_var_pet_ann

# PRECIP ANN

e_var_pr_msc %>% 
  st_apply(c(1,2), mean, na.rm = T, .fname = "ERA") -> e_var_pr_ann

# ****************************************************

# SPEI MSC ----

# c(e_var_pr, e_var_pet) %>%
#   mutate(wb = tp - pet) %>%
#   select(wb) %>%
# 
#   st_apply(c(1,2),
#            function(x){
# 
#              spei(x, scale = 3, na.rm = T)$fitted %>%
#                as.vector() %>%
#                ifelse(is.infinite(.), NA, .)
# 
#            },
#            FUTURE = T,
#            .fname = "time") %>%
# 
#         st_set_dimensions("time",
#                           values = st_get_dimension_values(e_var_pr, "time")) %>%
#   aperm(c(2,3,1)) -> e_var_spei
# 
# saveRDS(e_var_spei, "output/e_var_spei3.rds")

readRDS(here::here("output/e_var_spei3.rds")) -> e_var_spei

# *************************************************

# ERA SPEI MSC ----

map(seq_len(12), function(mth){
  
  e_var_spei %>%
    filter(month(time) == mth) %>% 
    st_apply(c(1,2),
             function(x){
               sum(x < -1.5, na.rm = T)
             },
             .fname = month.abb[mth])
      
    }) %>% 
      do.call(c, .) %>% 
      merge(name = "month") %>% 
      setNames("spei") -> e_var_spei_msc

# **************************************

# ERA SPEI ANN ----

map2(seq(1, 252, 12),
     seq(12, 253, 12),
     function(st, en){
       
       e_var_spei %>% 
         slice(time, st:en) %>% 
         st_apply(c(1,2),
                  function(x){
                    sum(x < -1.5, na.rm = T)
                  })
       
     }) %>% 
  do.call(c, .) %>% 
  merge() %>% 
  st_apply(c(1,2), mean, na.rm = T, .fname = "ERA") -> e_var_spei_ann


# *************************************************

rm(s_array, dates_1deg, days_mth, days_mth_1deg, index, J)

# *************************************************

st_read("~/bucket_mine/misc_data/ne_50m_land/ne_50m_land.shp") %>%
  mutate(a = 1) %>%
  select(a) %>%
  st_rasterize(st_as_stars(st_bbox(), dx = 0.25, dy = 0.25, values = NA)) -> r_land_rast

st_warp(r_land_rast %>%
          st_set_dimensions(which = c(1,2),
                            names = c("longitude", "latitude")),
        e_var_pet_ann) -> r_land_rast




         
