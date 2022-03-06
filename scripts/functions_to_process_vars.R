
# to shift central meridian ----

func_shift_mer <- function(variable){

  variable %>%
    str_split("/") %>%
    unlist() -> dir_var

  dir_var %>%
    .[-length(.)] %>%
    str_flatten("/") -> dir

  dir_var %>%
    .[length(.)] %>%
    str_remove(".nc") -> var

  system(str_glue("cdo sellonlatbox,-180,180,-90.25,90 {dir}/{var}.nc {dir}/{var}_shifted.nc"))
}




# to set limits ----

func_ncsub <- function(s_proxy, lim){
  
  # lon
  s_proxy %>% 
    st_get_dimension_values(1, where = "start") -> s_proxy_lon
  
  which.min(abs(s_proxy_lon - lim[1])) -> lon_start
  which.min(abs(s_proxy_lon - lim[3])) - lon_start -> lon_count
  
  # lat
  s_proxy %>% 
    st_get_dimension_values(2, where = "start") -> s_proxy_lat
  
  if(s_proxy_lat[1] > s_proxy_lat[2]){
    
    which.min(abs(s_proxy_lat - lim[4])) -> lat_start
    which.min(abs(s_proxy_lat - lim[2])) - lat_start -> lat_count
    
  } else {
    
    which.min(abs(s_proxy_lat - lim[2])) -> lat_start
    which.min(abs(s_proxy_lat - lim[4])) - lat_start -> lat_count
    
  }
  
  sub <- cbind(start = c(lon_start, lat_start, 1),
               count = c(lon_count, lat_count, NA))
  
  return(list(sub = sub,
              lon_ini = s_proxy_lon[lon_start],
              lat_ini = s_proxy_lat[lat_start]))
  
}





