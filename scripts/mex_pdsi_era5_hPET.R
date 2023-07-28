
source("scripts/00_setup.R")
library(scPDSI)

plan(multicore)



# foo <- array(rep(1:3, each = 9), dim = c(3,3,3)) %>% 
#   st_as_stars() %>% 
#   setNames("foo")
# 
# bar <- array(rep(4:6, each = 9), dim = c(3,3,3)) %>% 
#   st_as_stars() %>% 
#   setNames("bar")
# 
# c(foo, bar) %>% 
#   st_apply(c(1,2), sum)







dir_cat <- "/mnt/pers_disk/cat"
# dir.create(dir_cat)

# concatenate PET

# ff <- 
#   "/mnt/bucket_mine/era/monthly/total-pet-v-hPET/" %>% 
#   list.files(full.names = T) %>% 
#   str_flatten(" ")
# 
# "cdo cat {ff} {dir_cat}/pet.nc" %>% 
#   str_glue() %>% 
#   system()
# 
# # concatenate PRECIP
# 
# ff <- 
#   "/mnt/bucket_mine/era/monthly/mean-daily-precip/" %>% 
#   list.files(full.names = T) %>%
#   str_subset(str_flatten(1981:2021, "|")) %>% 
#   str_flatten(" ")
# 
# "cdo cat {ff} {dir_cat}/precip.nc" %>% 
#   str_glue() %>% 
#   system()




# 

fun_apply_pdsi <- function(x){
  
  if(any(is.na(x))){
    
    rep(NA, length(x)/2)
  
  } else {
    
    x_P <- x[1:length(time_vector)] %>% 
      {. * 1000} %>% 
      {. * days_mth}
    
    x_PE <- x[(length(time_vector)+1):(length(time_vector)*2)]
    
    
    res <- pdsi(P = x_P,
                PE = x_PE,
                sc = T,
                cal_start = 1,
                cal_end = 20
    ) %>% 
      .$X
      
    case_when(res > 5 ~ 5,
              res < -5 ~ -5,
              TRUE ~ as.numeric(res))
    
  }
}




#

pet_proxy <- 
  "/mnt/pers_disk/cat/pet.nc" %>% 
  read_ncdf(proxy = T)

precip_proxy <- 
  "/mnt/pers_disk/cat/precip.nc" %>% 
  read_ncdf(proxy = T)

time_vector <-
  st_get_dimension_values(precip_proxy, "time") %>% 
  as_date()

time_vector %>% 
  days_in_month() %>% 
  unname() -> days_mth
  

s_precip <- 
  precip_proxy[, 950:1100, 220:310, ] %>%
  st_as_stars() %>%
  drop_units()

s_pet <-
  pet_proxy[, 950:1100, 220:310, ] %>%
  st_as_stars()

tic()
s_pdsi <- 
  c(s_precip,
  s_pet,
  along = 3) %>%
  
  st_apply(c(1,2),
           fun_apply_pdsi,
           FUTURE = T,
           .fname = "time") %>%
  aperm(c(2,3,1)) %>%
  st_set_dimensions("time", values = time_vector)
toc()




prob_extr_drought <- 
  s_pdsi %>% 
  st_apply(c(1,2), function(x){
    
    sum(x <= -4)/length(x)*100
    
  },
  FUTURE = T,
  .fname = "prob")


write_stars(prob_extr_drought,
            "/mnt/bucket_mine/results/mexico_drought_ww/observational_prob_extr_drought.tif")



historic <-
  "/mnt/bucket_cmip5/Wellington/Drought/pdsi_cmip6_historical_6month_neg4_monthly_percentage_1971_2000_basd.tif" %>%
  read_stars() %>%
  drop_units() %>%
  .[,121:192,110:158] %>%
  setNames("hist") %>%
  mutate(hist = ifelse(hist == 0, 0.0001, hist))

ssp585 <- 
  "/mnt/bucket_cmip5/Wellington/Drought/pdsi_cmip6_ssp585_6month_neg4_monthly_percentage_2021_2050_basd.tif" %>% 
  read_stars() %>% 
  drop_units() %>% 
  .[,121:192,110:158] %>% 
  setNames("future") %>% 
  mutate(future = ifelse(future == 0, 0.0001, future))

s_diff <- ssp585-historic

s_diff %>% plot()
  setNames("change") %>% 
  mutate(change = case_when(change > 20 ~ 20,
                            change < -20 ~ -20,
                            TRUE ~ change))


write_stars(s_diff,
            "/mnt/bucket_mine/results/mexico_drought_ww/cmip6_change_prob_extr_drought.tif")







prob_extr_drought %>% 
  as_tibble() %>% 
  mutate(prob = case_when(prob > 40 ~ 40,
                          TRUE ~ prob)) %>% 
  ggplot(aes(longitude, latitude, fill = prob)) +
  geom_raster() +
  colorspace::scale_fill_binned_sequential("viridis",
                                           na.value = "transparent",
                                           n.breaks = 7) +
  coord_equal()




land <- "/mnt/bucket_mine/misc_data/ne_110m_land/ne_110m_land.shp" %>% 
  st_read() %>% 
  select(1)

ggplot() +
  geom_raster(data = s_diff %>% #plot(axes = T)
                as_tibble() %>% 
                mutate(prob = case_when(future > 40 ~ 40,
                                        TRUE ~ future)),
              aes(x, y, fill = prob)) +
  colorspace::scale_fill_binned_diverging(na.value = "transparent",
                                           n.breaks = 7,
                                          limits = c(-10, 40)) +
  geom_sf(data = land, fill = NA) +
  coord_sf(xlim = c(-120, -85), ylim = c(12, 35))



tif = system.file("tif/L7_ETMs.tif", package = "stars")
r = read_stars(tif, proxy = TRUE)
tiles = st_tile(nrow(r), ncol(r), 150, 150)

