


# pr %>% 
#   slice(time, 1:12) -> pr_1979
# 
# pev %>% 
#   slice(time, 1:12) -> pev_1979

c(pr, pev) %>% 
  setNames(c("pr", "pev")) -> s

c(pr, pev, along = 3) -> s2




library(scPDSI)

plan(multicore, workers = 6)

s %>% select(pr) %>% slice(longitude,1) %>% slice(latitude,1) %>% pull(1) %>% as.vector() -> pr_1
s %>% select(pev) %>% slice(longitude,1) %>% slice(latitude,1) %>% pull(1) %>% as.vector() -> pev_1
pdsi(pr_1, pev_1, sc = T) -> result_1
result_1$X %>% plot()

# s2 %>% 
#   st_apply(c(1,2), function(x){
#     
#     pdsi(x[1:504], x[505:1008], sc = T)$X %>% as.vector()
#     
#   },
#   FUTURE = T,
#   .fname = "time") -> s2_app
# 
# seq(as_date("1979-01-01"), as_date("2020-12-01"), by = "months") -> date_vector
# 
# s2_app %>% 
#   st_set_dimensions("time", values = date_vector) %>% 
#   aperm(c(2,3,1)) -> s2_app

s2_app %>% slice(longitude,1) %>% slice(latitude,1) %>% pull(1) %>% as.vector() %>% plot(type = "l")


s %>% 
  as_tibble() -> s_tbl

s_tbl %>% 
  group_by(longitude, latitude) %>% 
  nest() %>% 
  mutate(pdsi = map(data, function(df){
    
    pdsi(df$pr, df$pev)$X
    
  })) -> s_tbl


s_1979 %>% 
  st_apply(1:3, pdsi,
           FUTURE = T) -> s


1:81 %>% array(dim = c(3,3,3,3)) %>% st_as_stars() %>% .[[1]]

pdsi(s[1,1,1,], s[2,1,1,])

array(1:27, dim = c(3,3,3)) %>% st_as_stars() -> a
array(1:27, dim = c(3,3,3)) %>% st_as_stars() -> b
c(a,b) -> c

c %>% select(A1) %>% slice(X1, 1) %>% slice(X2, 1) %>% pull(1)

c %>% 
  st_apply(1:2, function(x){
    print(x)
  })


vds <- read_ncdf("data/scPDSI.cru_ts4.05early1.1901.2020.cal_1901_20.bams.2021.GLOBAL.1901.2020.nc")
vds %>% st_crop(lim) -> vds

vds %>% 
  filter(year(time) >= 1979) -> vds

vds %>% 
  st_warp(s2_app) %>% 
  {c(s2_app, ., along = 3)} %>% 
  
  st_apply(c(1,2), function(x){
    
    cor(x[1:504], x[505:1008])
    
  }, 
  FUTURE = T) -> s_cor

mapview::mapview(s_cor, at = c(-Inf, seq(-0.5, 0.5, 0.1), Inf), col.regions = colorRampPalette(diverging_hcl(palette = "Blue-Red", 11)))


# *************

library(SPEI)
library(scPDSI)
library(colorspace)

c(st_point(c(-120, 33)),
  st_point(c(-80, 13))) %>% 
  st_bbox() %>% 
  st_set_crs(4326) -> lim

# read_ncdf(str_glue("{dir_bucket_mine}/era/monthly/era5_2mtemp_dailymeans_monthlyres_rawfile.nc")) %>% 
#   st_crop(lim) -> tave
# 
# tave %>% 
#   mutate(t2m = units::set_units(t2m, degC),
#          t2m = units::drop_units(t2m)) -> tave
# 
# tave %>% 
#   # slice(time,1) %>% 
#   st_dim_to_attr(2) -> s_lat

# tave %>% 
#   split("time") %>% 
#   c(s_lat) %>%
#   merge() -> s
# 
# s %>% 
#   st_apply(c(1,2), function(x){
#     
#     x[length(x)] -> lat
#     x[-length(x)] -> x_ts
#     
#     thornthwaite(x_ts, lat)
#     
#   }, 
#   FUTURE = T,
#   .fname = "time") -> pet_thorn

# seq(as_date("1979-01-01"), as_date("2020-12-01"), by = "months") -> date_vector

# pet_thorn %>% 
#   st_set_dimensions("time", values = date_vector) %>% 
#   aperm(c(2,3,1)) -> pet_thorn

readRDS("output/starr_precip_monthlysums_monthlyres.rds") %>% 
  st_crop(lim) -> pr

# c(pr, pet_thorn, along = 3) -> s2
# 
# s2 %>% 
#   st_apply(c(1,2), function(x){
#     
#     pdsi(x[1:504], x[505:1008], sc = T)$X %>% as.vector()
#     
#   },
#   FUTURE = T,
#   .fname = "time") -> pdsi_thorn
# 
# pdsi_thorn %>% 
#   aperm(c(2,3,1)) -> pdsi_thorn

# vds <- read_ncdf("data/scPDSI.cru_ts4.05early1.1901.2020.cal_1901_20.bams.2021.GLOBAL.1901.2020.nc")
# vds %>% 
#   st_warp(pdsi_thorn) %>% 
#   st_crop(lim) %>% 
#   filter(year(time) >= 1979) -> vds
# 
# 
# vds %>% 
#   # st_warp(pdsi_thorn) %>% 
#   {c(pdsi_thorn, ., along = 3)} %>% 
#   
#   st_apply(c(1,2), function(x){
#     
#     cor(x[1:504], x[505:1008])
#     
#   }, 
#   FUTURE = T) -> s_cor
# 
# mapview(s_cor, at = c(-Inf, seq(-0.8, 0.8, 0.2), Inf), 
#         col.regions = colorRampPalette(diverging_hcl(palette = "Blue-Red", 11)))
# 
# s_cor %>% 
#   as_tibble() %>% 
#   filter(near(X, 0.5, 0.01)) %>% 
#   slice_sample(n = 1) %>% 
#   {c(.$longitude, .$latitude)} -> coords
# 
# vds %>% 
#   # st_warp(pdsi_thorn) %>% 
#   as_tibble() %>% 
#   filter(longitude == coords[1],
#          latitude == coords[2]) %>% 
#   pull(scpdsi) -> vect_vds
# 
# pdsi_thorn %>% 
#   as_tibble() %>% 
#   filter(longitude == coords[1],
#          latitude == coords[2]) %>% 
#   pull(X) -> vect_pdsi_thorn
# 
# tibble(vds = vect_vds,
#        thorn = vect_pdsi_thorn,
#        time = date_vector) %>% 
#   pivot_longer(-time, names_to = "src", values_to = "pdsi") %>% 
#   
#   ggplot(aes(x = time, y = pdsi, color = src)) +
#   geom_line()
# 
# vds %>% 
#   # st_warp(pdsi_thorn) %>% 
#   {c(pdsi_thorn, ., along = list(foo = c("bar1", "bar2")))} %>% 
#   # st_set_dimensions(3, values = date_vector) %>% 
#   split("foo") %>% 
#   as_tibble() %>% 
#   group_by(time) %>% 
#   summarize(cor = cor(bar1, bar2, use = "complete.obs")) %>% 
#   mutate(time = date_vector) -> t_cor 
#   
# t_cor %>%   
#   ggplot(aes(x = time, y = cor)) +
#   geom_line() +
#   geom_point()
# 
# t_cor %>% 
#   as_tibble() %>% 
#   filter(near(cor, 0.6, 0.01)) %>% 
#   slice_sample(n = 1) %>% 
#   pull(time) -> t_time
# 
# pdsi_thorn[is.na(adrop(vds[,,,1]))] <- NA
# 
# pdsi_thorn %>% st_set_dimensions(3, values = date_vector) %>%
#   filter(time == t_time) %>% 
#   adrop() %>% 
#   as_tibble() -> m_1
# 
# vds %>%
#   filter(time == t_time) %>% 
#   adrop() %>% 
#   as_tibble() -> m_2
#   
# left_join(m_1, m_2, by = c("longitude", "latitude")) %>% 
#   pivot_longer(c(X, scpdsi), names_to = "src", values_to = "pdsi") %>% 
#   mutate(pdsi = case_when(pdsi < -5 ~ -5,
#                           pdsi > 5 ~ 5,
#                           TRUE ~ pdsi)) %>% 
#   ggplot(aes(x = longitude, y = latitude, fill = pdsi)) +
#   geom_raster() +
#   scale_fill_continuous_diverging(palette = "Blue-Red", rev = T) +
#   facet_wrap(~src, ncol = 2) +
#   coord_cartesian()

# **********

# HARGRAVES

# library(SPEI)
# library(scPDSI)
# library(colorspace)
# 
# plan(multicore, workers = 6)
# 
# c(st_point(c(-120, 33)),
#   st_point(c(-80, 13))) %>% 
#   st_bbox() %>% 
#   st_set_crs(4326) -> lim
# 
# seq(as_date("1979-01-01"), as_date("2020-12-01"), by = "months") -> date_vector

# ***

# readRDS("output/starr_precip_monthlysums_monthlyres.rds") %>% 
#   st_crop(lim) -> pr

# ***

# EXTERNAL RADIATION
# read_ncdf(str_glue("{dir_bucket_mine}/era/monthly/era5_monthly_mean_daily_radiation.nc")) %>% 
#   st_crop(lim) -> Ra
# 
# Ra %>%
#   select(tisr) -> Ra # top (external?)
# 
# Ra %>% 
#   mutate(tisr = units::set_units(tisr, MJ/m^2),
#          tisr = units::drop_units(tisr)) -> Ra

# ***

# TMAX
# map(seq_len(42), function(y){
#   
#   list.files("~/bucket_mine/era/monthly/2m_temperature_max/", full.names = T) %>% 
#     .[y] %>% 
#     readRDS() %>% 
#     st_crop(lim)
#   
# }) -> tmax
# 
# tmax %>% 
#   do.call(c, .) -> tmax
# 
# # ***
# 
# # TMIN
# map(seq_len(42), function(y){
#   
#   print(y)
#   
#   list.files("~/bucket_mine/era/monthly/2m_temperature_min/", full.names = T) %>% 
#     .[y] %>% 
#     readRDS() %>% 
#     st_crop(lim)
#   
# }) -> tmin
# 
# tmin %>% 
#   do.call(c, .) -> tmin

# ***

# PET
# c(Ra, tmax, tmin) %>%
#   st_apply(c(1,2), function(x){
#     
#     hargreaves(Tmin = x[1009:1512],
#                Tmax = x[505:1008],
#                Ra = x[1:504],
#                na.rm = T)
#     
#   },
#   FUTURE = T,
#   .fname = "time") -> pet_harg
# 
# pet_harg %>% 
#   st_set_dimensions("time", values = date_vector) %>% 
#   aperm(c(2,3,1)) -> pet_harg

# ***

# SCPDSI
# c(pr, pet_harg, along = 3) %>%
#   st_apply(c(1,2), function(x){
#     
#     pdsi(x[1:504], x[505:1008], sc = T)$X %>% as.vector()
#     
#   },
#   FUTURE = T,
#   .fname = "time") %>% 
#   st_set_dimensions("time", values = date_vector) %>% 
#   aperm(c(2,3,1)) -> pdsi_harg

# ***

# VAN DER SCHRIER

# vds <- read_ncdf("data/scPDSI.cru_ts4.05early1.1901.2020.cal_1901_20.bams.2021.GLOBAL.1901.2020.nc")
# vds %>% 
#   st_warp(pdsi_harg) %>% 
#   st_crop(lim) %>% 
#   filter(year(time) >= 1979) -> vds
# 
# vds %>% 
#   st_set_dimensions("time", values = date_vector) -> vds

# ***

# CORRELATION MAP

# c(pdsi_harg, vds, along = 3) %>% 
#   st_apply(c(1,2), function(x){
#     
#     cor(x[1:504], x[505:1008])
#     
#   }, 
#   FUTURE = T) -> s_cor

# s_cor %>% 
#   mutate(X = case_when(X < -0.9 ~ -0.9,
#                        X > 0.9 ~ 0.9,
#                        TRUE ~ X)) %>% 
#   {
#     ggplot() +
#       geom_stars(data = .) +
#       scale_fill_continuous_diverging(palette = "Blue-Red", rev = T, na.value = "grey80") +
#       coord_fixed()
#   }

# ts of HIGH cor
# s_cor %>% 
#   as_tibble() %>% 
#   filter(near(X, 0.85, 0.01)) %>% 
#   slice_sample(n = 1) %>% 
#   {c(.$longitude, .$latitude)} -> coords
# 
# vds %>%
#   as_tibble() %>% 
#   filter(longitude == coords[1],
#          latitude == coords[2]) %>% 
#   pull(scpdsi) -> vect_vds
# 
# pdsi_harg %>% 
#   as_tibble() %>% 
#   filter(longitude == coords[1],
#          latitude == coords[2]) %>% 
#   pull(X) -> vect_pdsi
# 
# tibble(vds = vect_vds,
#        era5 = vect_pdsi,
#        time = date_vector) %>% 
#   pivot_longer(-time, names_to = "src", values_to = "pdsi") %>% 
#   
#   ggplot(aes(x = time, y = pdsi, color = src)) +
#   geom_line()

# ts of LOW cor
# s_cor %>% 
#   as_tibble() %>% 
#   filter(near(X, 0.2, 0.01)) %>% 
#   slice_sample(n = 1) %>% 
#   {c(.$longitude, .$latitude)} -> coords
# 
# vds %>%
#   as_tibble() %>% 
#   filter(longitude == coords[1],
#          latitude == coords[2]) %>% 
#   pull(scpdsi) -> vect_vds
# 
# pdsi_harg %>% 
#   as_tibble() %>% 
#   filter(longitude == coords[1],
#          latitude == coords[2]) %>% 
#   pull(X) -> vect_pdsi
# 
# tibble(vds = vect_vds,
#        era5 = vect_pdsi,
#        time = date_vector) %>% 
#   pivot_longer(-time, names_to = "src", values_to = "pdsi") %>% 
#   
#   ggplot(aes(x = time, y = pdsi, color = src)) +
#   geom_line()

# ***

# TS OF COR

# vds %>% 
#   {c(pdsi_harg, ., along = list(foo = c("era5", "vds")))} %>% 
#   split("foo") %>% 
#   as_tibble() %>% 
#   group_by(time) %>% 
#   summarize(cor = cor(era5, vds, use = "complete.obs")) -> t_cor 
# 
# t_cor %>%   
#   ggplot(aes(x = time, y = cor)) +
#   geom_point() +
#   geom_smooth()

# comparison of HIGH cor maps
# pdsi_harg[is.na(adrop(vds[,,,1]))] <- NA

# t_cor %>% 
#   as_tibble() %>% 
#   filter(near(cor, 0.6, 0.01)) %>% 
#   slice_sample(n = 1) %>% 
#   pull(time) -> t_time
# 
# pdsi_harg %>%
#   filter(time == t_time) %>% 
#   adrop() %>% 
#   as_tibble() -> m_1
# 
# vds %>%
#   filter(time == t_time) %>% 
#   adrop() %>% 
#   as_tibble() -> m_2
# 
# left_join(m_1, m_2, by = c("longitude", "latitude")) %>%
#   rename(era5 = 3,
#          vds = 4) %>% 
#   pivot_longer(3:4, names_to = "src", values_to = "pdsi") %>% 
#   mutate(pdsi = case_when(pdsi < -5 ~ -5,
#                           pdsi > 5 ~ 5,
#                           TRUE ~ pdsi)) %>% 
#   ggplot(aes(x = longitude, y = latitude, fill = pdsi)) +
#   geom_raster() +
#   scale_fill_continuous_diverging(palette = "Blue-Red", rev = T, na.value = "grey80") +
#   facet_wrap(~src, ncol = 2) +
#   coord_equal()

# # comparison of LOW cor maps
# t_cor %>% 
#   as_tibble() %>% 
#   filter(near(cor, 0.2, 0.01)) %>% 
#   slice_sample(n = 1) %>% 
#   pull(time) -> t_time
# 
# pdsi_harg %>%
#   filter(time == t_time) %>% 
#   adrop() %>% 
#   as_tibble() -> m_1
# 
# vds %>%
#   filter(time == t_time) %>% 
#   adrop() %>% 
#   as_tibble() -> m_2
# 
# left_join(m_1, m_2, by = c("longitude", "latitude")) %>%
#   rename(era5 = 3,
#          vds = 4) %>% 
#   pivot_longer(3:4, names_to = "src", values_to = "pdsi") %>% 
#   mutate(pdsi = case_when(pdsi < -5 ~ -5,
#                           pdsi > 5 ~ 5,
#                           TRUE ~ pdsi)) %>% 
#   ggplot(aes(x = longitude, y = latitude, fill = pdsi)) +
#   geom_raster() +
#   scale_fill_continuous_diverging(palette = "Blue-Red", rev = T, na.value = "grey80") +
#   facet_wrap(~src, ncol = 2) +
#   coord_equal()





 