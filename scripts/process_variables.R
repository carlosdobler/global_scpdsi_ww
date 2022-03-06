
source("scripts/00-set-up.R")
plan(multicore, workers = 6)

# datasrs <- c("ERA", "REMO", "vdS")[1]

# mex
# c(st_point(c(-120, 33)),
#   st_point(c(-80, 13))) %>%
#   st_bbox() %>%
#   st_set_crs(4326) -> lim

# central europe
c(st_point(c(-11, 65)),
  st_point(c(39, 32))) %>% 
  st_bbox() %>% 
  st_set_crs(4326) -> lim

source("scripts/functions_to_process_vars.R")
source("scripts/penman_mine.R")
library(scPDSI)



# **************************************************************************************************

#  ___  __       
# |__  |__)  /\  
# |___ |  \ /~~\ 
#                

# **************************************************************************************************

# ERA SET LIMITS ----

"~/bucket_mine/era/monthly/era5_geopotential_shifted.nc" %>%
  read_ncdf(ncsub = cbind(start = c(1,1,1), count = c(NA,NA,1))) %>%
  adrop() -> s_proxy

lims <-
  func_ncsub(s_proxy, lim)


 
# ERA IMPORT VARS ----

# all except tasmax and min:
variable <- c("meantemp", "dewpointtemp", "radiation", "10mwindspeed", "precip")

map(variable, function(var){

  list.files("~/bucket_mine/era/monthly/", full.names = T) %>%
    .[str_detect(., var)] %>%
    .[str_detect(., "shifted")] %>%

    read_ncdf(ncsub = lims$sub) -> e_vars

  st_dimensions(e_vars)$longitude$offset <- lims$lon_ini
  st_dimensions(e_vars)$latitude$offset <- lims$lat_ini

  return(e_vars)

}) -> e_vars

e_vars %>%
  do.call(c, .) %>%
  select(-tisr) %>%
  st_set_dimensions("time",
                    values = st_get_dimension_values(., "time") %>%
                      as_date()) -> e_vars


# tasmax and min:
# (no longer needed; using mean temp for PM PET)
# map(c("max", "min"), function(var){
# 
#   print(str_glue("Processing {var}"))
# 
#   list.files(str_glue("~/bucket_mine/era/monthly/2m_temperature_{var}/"), full.names = T) %>%
#     map(readRDS) %>%
#     do.call(c, .) %>%
#     setNames(str_glue("tas{var}")) %>%
#     st_warp(e_vars)
# 
# }) %>%
#   do.call(c, .) %>%
#   c(e_vars) -> e_vars



# ERA MODIFY VARS ----

e_vars %>%
  mutate(t2m = t2m %>% set_units(degC) %>% set_units(NULL),
         d2m = d2m %>% set_units(degC) %>% set_units(NULL),
         ssrd = ssrd %>% set_units(MJ/m^2) %>% set_units(NULL),
         si2m = si10 %>% set_units(NULL) %>% {. * 4.87 / log(67.8 * 10 - 5.42)},
         tp = tp %>% set_units(mm) %>% set_units(NULL)) %>% 
  select(-si10) -> e_vars

e_vars %>%
  st_get_dimension_values("time") %>%
  days_in_month() %>%
  unname() -> days_mth

e_vars %>%
  select(tp) %>%
  st_apply(c(1,2), function(x){
    x * days_mth
  },
  .fname = "time") %>%
  aperm(c(2,3,1)) %>%
  st_set_dimensions("time", values = st_get_dimension_values(e_vars, "time")) %>%
  {c(e_vars %>% select(-tp), .)} -> e_vars

rm(days_mth)

saveRDS(e_vars, here::here("output/e_vars_1979_2020.rds"))



# ERA Z AND LAT ----

"~/bucket_mine/era/monthly/era5_geopotential_shifted.nc" %>% 
  read_ncdf(ncsub = lims$sub[1:2,] %>% rbind(1)) %>% adrop() -> e_var_z

st_dimensions(e_var_z)$longitude$offset <- lims$lon_ini
st_dimensions(e_var_z)$latitude$offset <- lims$lat_ini

e_var_z %>%
  mutate(z = z %>% set_units(NULL) %>% {./9.80665},
         z = ifelse(z < 0, 0, z)) -> e_var_z

saveRDS(e_var_z, here::here("output/e_var_z.rds"))

e_var_z %>%
  st_dim_to_attr(2) -> e_var_lat

saveRDS(e_var_lat, here::here("output/e_var_lat.rds"))



# ERA OTHER VARS ----

# potential evaporation

"~/bucket_mine/era/monthly/era5_monthly_mean_daily_potevap_shifted.nc" %>%
  read_ncdf(ncsub = lims$sub) -> s

st_dimensions(s)$longitude$offset <- lims$lon_ini
st_dimensions(s)$latitude$offset <- lims$lat_ini

s %>%   
  mutate(pev = pev %>% set_units("mm") %>% set_units(NULL)) -> e_var_pev

e_vars %>%
  st_get_dimension_values("time") %>%
  days_in_month() %>%
  unname() -> days_mth

e_var_pev %>% 
  st_apply(c("longitude", "latitude"), function(x){
    x * days_mth
  },
  .fname = "time") %>% 
  
  st_set_dimensions("time", values = st_get_dimension_values(e_vars, "time")) %>%
  aperm(c(2,3,1)) -> e_var_pev

saveRDS(e_var_pev, here::here("output/e_var_pev_1979_2020.rds"))


# surface pressure

# "~/bucket_mine/era/monthly/era5_monthly_mean_surfpressure_shifted.nc" %>%
#   read_ncdf(ncsub = lims$sub) -> s
# 
# st_dimensions(s)$longitude$offset <- lims$lon_ini
# st_dimensions(s)$latitude$offset <- lims$lat_ini
# 
# s %>%   
#   mutate(sp = sp %>% set_units(kPa) %>% set_units(NULL)) %>% 
#   st_set_dimensions("time", values = st_get_dimension_values(e_vars, "time")) -> e_var_press



# ERA PET PM ----

abind(
  e_vars[[1]], # tmean
  e_vars[[3]], # ssrd (radiation)
  e_vars[[4]], # si10 (wind speed)
  e_vars[[2]], # dewpoint
  e_var_z[[1]], # z
  e_var_lat[[1]], # lat
  along = 3) -> s_array

names(dim(s_array)) <- c("longitude", "latitude", "time")

s_array %>% 
  st_as_stars() -> s_array

c(1,
  dim(e_vars[[1]])[3],
  dim(e_vars[[3]])[3],
  dim(e_vars[[4]])[3],
  dim(e_vars[[2]])[3],
  1,
  1) %>% 
  unname() %>% 
  cumsum() -> index

cbind(index[-length(index)],
      index[-1]-1) -> index

s_array %>% 
  st_apply(c(1,2), function(x){
    
    penman_mine(Tmean = x[index[1,1]:index[1,2]],
                Rs = x[index[2,1]:index[2,2]],
                u2 = x[index[3,1]:index[3,2]],
                Tdew = x[index[4,1]:index[4,2]],
                z = x[index[5,1]:index[5,2]],
                lat = x[index[6,1]:index[6,2]])
    
  }, 
  .fname = "time") %>%
  aperm(c(2,3,1)) %>% 
  setNames("pet") -> e_pet_pm

st_dimensions(e_pet_pm) <- st_dimensions(e_vars)

rm(s_array, index)

saveRDS(e_pet_pm, here::here("output/e_pet_pm_1979_2020.rds"))



# ERA PDSI PM ----

c(e_vars %>% select(tp), e_pet_pm, along = 3) %>% 
  st_apply(c(1,2), function(x){
    
    pdsi(P = x[1:504], 
         PE = x[505:1008], 
         sc = T)$X %>% as.vector()
    
  },
  FUTURE = T,
  future.seed = NULL,
  .fname = "time") %>% 
  
  st_set_dimensions("time", values = st_get_dimension_values(e_vars, "time")) %>% 
  aperm(c(2,3,1)) %>% 
  setNames("pdsi") -> e_pdsi_pm

saveRDS(e_pdsi_pm, here::here("output/e_pdsi_pm_1979_2020.rds"))



# **************************************************************************************************

#                    __   ___  __      __   __        __     ___  __  
# \  /  /\  |\ |    |  \ |__  |__)    /__` /  ` |__| |__) | |__  |__) 
#  \/  /~~\ | \|    |__/ |___ |  \    .__/ \__, |  | |  \ | |___ |  \ 
#                                                                     

# **************************************************************************************************

# VDS PDSI ----

here::here(
  "data", "scPDSI.cru_ts4.05early1.1901.2020.cal_1901_20.bams.2021.GLOBAL.1901.2020.nc"
  ) %>% 
  #read_stars() %>% st_get_dimension_values(3) %>% {which(. == "1979-01-01")}
  read_ncdf(ncsub = cbind(start = c(1,1,937),
                          count = c(NA, NA, NA))) -> v_pdsi

v_pdsi %>% 
  st_set_dimensions("time", values = st_get_dimension_values(e_vars, "time")) %>% 
  st_warp(e_var_z) -> v_pdsi

saveRDS(v_pdsi, here::here("output/v_pdsi_1979_2020.rds"))



# **************************************************************************************************

#  _    _         _
# |_)  |_  |\/|  / \
# | \  |_  |  |  \_/
# 


# REMO SET LIMITS ----------------------------------------------------------------------------------

list.files("~/bucket_mine/remo/monthly", full.names = T) %>% 
  .[str_detect(., "regrid")] %>% 
  first() %>%
  read_ncdf(ncsub = cbind(start = c(1,1,1), count = c(NA,NA,1))) %>%
  adrop() -> s_proxy

lims <- 
  func_ncsub(s_proxy, lim)




# REMO IMPORT VARS ---------------------------------------------------------------------------------

variable <- c("tasmax", "tasmin", "hurs", "rsds", "sfcWind", "pr")

map(c("Had", "MPI", "Nor") %>% set_names(), function(model){

  map(variable, function(var){

    print(str_glue("Processing model {model} - variable {var}"))

    list.files("~/bucket_mine/remo/monthly/", full.names = T) %>%
      .[str_detect(., model)] %>%
      .[str_detect(., "regrid")] %>%
      .[str_detect(., var)] %>%

      map(function(f){

        f %>%
          read_ncdf(ncsub = lims$sub) %>%
          suppressMessages() -> s

        st_dimensions(s)$lon$offset <- lims$lon_ini
        st_dimensions(s)$lat$offset <- lims$lat_ini

        return(s)

      }) %>%
      do.call(c, .) -> s_var

    s_var %>%
      st_get_dimension_values("time") %>%
      as.character() %>%
      as_date() %>%
      {st_set_dimensions(s_var, "time", values = .)} # %>%
      # filter(year(time) >= 1979,
      #        year(time) <= 2020)

  }) %>% 
    do.call(c, .)

}) -> r_models_vars

rm(variable)




# REMO MODIFY --------------------------------------------------------------------------------------

func_dewpoint <- function(hurs, tasmean){

  hurs <- set_units(hurs, NULL)
  tasmean <- set_units(tasmean, NULL)

  243.04*(log(hurs/100)+((17.625*tasmean)/(243.04+tasmean)))/
    (17.625-log(hurs/100)-((17.625*tasmean)/(243.04+tasmean)))

}

r_models_vars %>%
  imap(function(s_vars, i){

    print(str_glue("Processing {i}"))

    s_vars %>%
      select(-pr) %>%
      mutate(tasmax = tasmax %>% set_units(degC) %>% set_units(NULL),
             tasmin = tasmin %>% set_units(degC) %>% set_units(NULL),
             rsds = rsds %>% set_units(MJ/d/m^2) %>% set_units(NULL),
             sfcWind = sfcWind %>% set_units(NULL),
             tasmean = (tasmax + tasmin)/2,
             dewpoint = func_dewpoint(hurs, tasmean)) %>%
      select(-hurs, -tasmax, -tasmin) -> s_no_pr

    if(i == "Had"){

      s_vars %>%
        select(pr) %>%
        mutate(pr = pr %>%
                 set_units(kg/m^2/d) %>%
                 set_units(NULL) %>%
                 {.*30}) -> s_pr

    } else {

      s_vars %>%
        st_get_dimension_values("time") -> dates

      dates %>%
        days_in_month() %>%
        unname() -> days_mth

      s_vars %>%
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

saveRDS(r_models_vars, "output/r_models_vars_1971_20XX.rds")
# readRDS("output/r_models_vars_1971_20XX.rds") -> r_models_vars



# REMO Z AND LAT -----------------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_geopotential_shifted.nc" %>% 
  read_ncdf() %>% 
  slice(time, 1) -> r_var_z

r_var_z %>%
  mutate(z = z %>% set_units(NULL) %>% {./9.80665},
         z = ifelse(z < 0, 0, z)) %>% 
  setNames("z") -> r_var_z

r_var_z %>% 
  st_warp(r_models_vars[[1]] %>% slice(time, 1)) -> r_var_z

saveRDS(r_var_z, "output/r_var_z.rds")
# readRDS("output/r_var_z.rds") -> r_var_z

r_var_z %>% 
  st_dim_to_attr(2) -> r_var_lat

saveRDS(r_var_lat, "output/r_var_lat.rds")
# readRDS("output/r_var_lat.rds") -> r_var_lat




# REMO PET 1979-2020 -------------------------------------------------------------------------------

r_models_vars %>%
  imap(function(s_vars, i){
    
    s_vars %>% 
      filter(year(time) >= 1979,
             year(time) <= 2020) -> s_vars
    
    print(str_glue("Processing {i}"))
    
    abind(
      s_vars[[3]], # tasmean
      s_vars[[1]], # rsds
      s_vars[[2]], # sfcWind
      s_vars[[4]], # dewpoint
      r_var_z[[1]], # z
      r_var_lat[[1]], # lat
      along = 3) -> s_array
    
    names(dim(s_array)) <- c("lon", "lat", "time")
    
    s_array %>% 
      st_as_stars() -> s_array
    
    c(1,
      dim(s_vars[[3]])[3],
      dim(s_vars[[1]])[3],
      dim(s_vars[[2]])[3],
      dim(s_vars[[4]])[3],
      1,
      1) %>% 
      unname() %>% 
      cumsum() -> index
    
    cbind(index[-length(index)],
          index[-1]-1) -> index
    
    tic("applying")
    s_array %>% 
      st_apply(c(1,2), function(x){
        
        penman_mine(Tmean = x[index[1,1]:index[1,2]],
                    Rs = x[index[2,1]:index[2,2]],
                    u2 = x[index[3,1]:index[3,2]],
                    Tdew = x[index[4,1]:index[4,2]],
                    z = x[index[5,1]:index[5,2]],
                    lat = x[index[6,1]:index[6,2]],
                    model_i = i)
        
      }, 
      .fname = "time") -> s_pet
    toc()
    
    s_pet %>% 
      aperm(c(2,3,1)) %>% 
      setNames("pet") -> s_pet
    
    st_dimensions(s_pet) <- st_dimensions(s_vars) 
    
    return(s_pet)           
    
  }) -> r_models_pet

saveRDS(r_models_pet, "output/r_models_pet_1979_2020.rds")




# REMO PDSI 1979-2020 ------------------------------------------------------------------------------

r_models_vars %>% 
  map(select, "pr") %>% 
  map(filter, year(time) >= 1979, year(time) <= 2020) -> r_models_pr

map2(r_models_pr, r_models_pet, function(s_pr, s_pet){
  
  c(s_pr, s_pet, along = 3) %>% 
    st_apply(c(1,2),
             function(x){
               
               c(1,
                 dim(s_pr)[3],
                 dim(s_pet)[3]) %>% 
                 unname() %>% 
                 cumsum() -> index
               
               cbind(index[-length(index)],
                     index[-1]-1) -> index
               
               pdsi(P = x[index[1,1]:index[1,2]], 
                    PE = x[index[2,1]:index[2,2]], 
                    sc = T)$X %>% as.vector()
               
             },
             FUTURE = T,
             future.seed = NULL,
             .fname = "time") %>% 
    
    st_set_dimensions("time", values = st_get_dimension_values(s_pet, "time")) %>% 
    aperm(c(2,3,1)) %>% 
    setNames("pdsi")
  
}) -> r_models_pdsi

saveRDS(r_models_pdsi, "output/r_models_pdsi_1979_2020.rds")



# REMO PET 1 DEG -----------------------------------------------------------------------------------


r_models_vars %>%
  imap(function(s_vars, i){
    
    print(str_glue("Processing {i}"))
    
    case_when(i == "Had" ~ 2014,
              i == "MPI" ~ 2006,
              i == "Nor" ~ 2019) -> mid_yr
    
    # date_vector_1deg <- seq(as_date(str_c(mid_yr-10,"-01-01")), 
    #                         as_date(str_c(mid_yr+10,"-12-01")),
    #                         by = "months")
    
    s_vars %>% 
      filter(year(time) >= mid_yr-10,
             year(time) <= mid_yr+10) -> s_vars
    
    abind(
      s_vars[[3]], # tasmean
      s_vars[[1]], # rsds
      s_vars[[2]], # sfcWind
      s_vars[[4]], # dewpoint
      r_var_z[[1]], # z
      r_var_lat[[1]], # lat
      along = 3) -> s_array
    
    names(dim(s_array)) <- c("lon", "lat", "time")
    
    s_array %>% 
      st_as_stars() -> s_array
    
    c(1,
      dim(s_vars[[3]])[3],
      dim(s_vars[[1]])[3],
      dim(s_vars[[2]])[3],
      dim(s_vars[[4]])[3],
      1,
      1) %>% 
      unname() %>% 
      cumsum() -> index
    
    cbind(index[-length(index)],
          index[-1]-1) -> index
    
    tic("applying")
    s_array %>% 
      st_apply(c(1,2), function(x){
        
        penman_mine(Tmean = x[index[1,1]:index[1,2]],
                    Rs = x[index[2,1]:index[2,2]],
                    u2 = x[index[3,1]:index[3,2]],
                    Tdew = x[index[4,1]:index[4,2]],
                    z = x[index[5,1]:index[5,2]],
                    lat = x[index[6,1]:index[6,2]],
                    model_i = i)
        
      }, 
      .fname = "time") -> s_pet
    toc()
    
    s_pet %>% 
      aperm(c(2,3,1)) %>% 
      setNames("pet") -> s_pet
    
    st_dimensions(s_pet) <- st_dimensions(s_vars) 
    
    return(s_pet)           
    
  }) -> r_models_pet

saveRDS(r_models_pet, "output/r_models_pet_1deg.rds")




# REMO PDSI 1 DEG ----------------------------------------------------------------------------------

c("Had", "MPI", "Nor") %>% set_names() %>% 
  map(function(i){
    
    print(str_glue("Processing {i}"))
    
    case_when(i == "Had" ~ 2014,
              i == "MPI" ~ 2006,
              i == "Nor" ~ 2019) -> mid_yr
    
    r_models_vars %>% 
      pluck(i) %>%
      select(pr) %>% 
      filter(year(time) >= mid_yr-10,
             year(time) <= mid_yr+10) -> s_pr
    
    r_models_pet %>% 
      pluck(i) %>% 
      filter(year(time) >= mid_yr-10,
             year(time) <= mid_yr+10) -> s_pet
    
    c(s_pr, s_pet, along = 3) %>% 
      st_apply(c(1,2),
               function(x){
                 
                 c(1,
                   dim(s_pr)[3],
                   dim(s_pet)[3]) %>% 
                   unname() %>% 
                   cumsum() -> index
                 
                 cbind(index[-length(index)],
                       index[-1]-1) -> index
                 
                 pdsi(P = x[index[1,1]:index[1,2]], 
                      PE = x[index[2,1]:index[2,2]], 
                      sc = T)$X %>% as.vector()
                 
               },
               FUTURE = T,
               future.seed = NULL,
               .fname = "time") %>% 
      
      st_set_dimensions("time", values = st_get_dimension_values(s_pet, "time")) %>% 
      aperm(c(2,3,1)) %>% 
      setNames("pdsi")
    
  }) -> r_models_pdsi

saveRDS(r_models_pdsi, "output/r_models_pdsi_1deg.rds")



# REMO PDSI RISK BUCKET 1979-2030 ------------------------------------------------------------------

"~/bucket_risk/RCM_regridded_data/REMO2015/EUR/monthly/palmer_drought_severity_index/pdsi_EUR_HadGEM2-ES_1970_2099.nc" %>%
  read_ncdf(ncsub = cbind(start = c(1,1,1),
                          count = c(NA,NA,1)),
            make_time = F) %>% 
  adrop() -> s_proxy

lims <- 
  func_ncsub(s_proxy, lim)


map(c("Had", "MPI", "Nor") %>% set_names(), function(model){
  
  "~/bucket_risk/RCM_regridded_data/REMO2015/EUR/monthly/palmer_drought_severity_index/" %>% 
    list.files(full.names = T) %>% 
    .[str_detect(., model)] %>% 
    
    read_ncdf(ncsub = rbind(lims$sub[1:2,], c(1,1560)), make_time = F) %>% 
    suppressMessages() -> s
  
  s %>% 
    st_set_dimensions("lon", 
                      values = s %>% st_get_dimension_values("lon") %>% round(1)) -> s
  
  s %>% 
    st_set_dimensions("lat", 
                      values = s %>% st_get_dimension_values("lat") %>% round(1)) -> s
  
  st_set_crs(s, 4326) -> s
  
  s %>% 
    st_set_dimensions("time",
                      values = seq(as_date("1970-01-01"), as_date("2099-12-01"), by = "month")) %>% 
    filter(year(time) >= 1979,
           year(time) <= 2030) -> s
  
  return(s)
  
}) -> r_models_pdsi_ir

saveRDS(r_models_pdsi_ir, "output/r_models_pdsi_ir_1979_2030.rds")
    
 


# REMO PDSI RISK BUCKET 1 DEG ----------------------------------------------------------------------
r_models_pdsi_ir %>% 
  imap(function(s, i){
    
    case_when(i == "Had" ~ 2014,
              i == "MPI" ~ 2006,
              i == "Nor" ~ 2019) -> mid_yr
    
    s %>% 
      filter(year(time) >= mid_yr-10,
             year(time) <= mid_yr+10)
    
  }) -> r_models_pdsi_ir
  
saveRDS(r_models_pdsi_ir, "output/r_models_pdsi_ir_1deg.rds")  




# **************************************************************************************************

#                 __                 __       
# |     /\  |\ | |  \     |\/|  /\  /__` |__/ 
# |___ /~~\ | \| |__/     |  | /~~\ .__/ |  \ 
#                                             

# LAND MASK ----
# warped at ERA resolution

st_read("~/bucket_mine/misc_data/ne_50m_land/ne_50m_land.shp") %>%                                   # ch: land mask
  mutate(a = 1) %>%
  select(a) %>%
  st_rasterize(st_as_stars(st_bbox(), dx = 0.25, dy = 0.25, values = NA)) -> land_rast

st_warp(land_rast %>%
          st_set_dimensions(which = c(1,2),
                            names = c("longitude", "latitude")),
        e_var_z) -> land_rast                                                                        # ch: land mask

  



