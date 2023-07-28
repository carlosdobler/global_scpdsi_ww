
# PRECIPITATION ----

# CONVERT DAILY MEANS INTO MONTHLY TOTALS

dir <- str_glue("{dir_bucket_mine}/era/monthly/")

# raw_file <- str_glue("{dir}/era5_totprecip_monthlymeans_rawfile.nc")
raw_file <- str_glue("{dir}/era5_potevap_dailymeans_monthlyres_rawfile.nc")

read_ncdf(raw_file) -> s

s %>% 
  st_get_dimension_values(3) %>% 
  as_date() -> date_vector

tibble(time = date_vector,
       yr = year(time),
       mo = month(time)) %>%
  rowid_to_column("pos") %>% 
         
  group_by(yr) %>% 
  summarize(pos_i = min(pos),
            pos_f = max(pos)) -> tb_time_pos

tb_time_pos %>%
  pmap(function(yr, pos_i, pos_f){
    
    print(str_glue("Processing: {yr}"))
    
    s %>% 
      slice(time, pos_i:pos_f) -> s_year
    
    map(seq_len(12), function(mth){
      
      str_glue("{yr}-{mth}-01") %>% 
        days_in_month() %>% 
        unname() -> days_mth
      
      s_year %>% 
        slice(time, mth) %>% 
        rename(var = 1) %>% 
        mutate(var = var %>% units::set_units("mm"),
               var = var * days_mth) -> s_monthlytotals
      
    }) %>% 
      do.call(c, .) %>% 
      merge(name = "time") %>% 
      st_set_dimensions("time", 
                        values = date_vector[str_detect(date_vector, as.character(yr))])
    
  }) %>% 
  do.call(c, .) -> foo

rm(s)

# saveRDS(foo, "output/starr_precip_monthlysums_monthlyres.rds")
saveRDS(foo, "output/starr_potevap_monthlysums_monthlyres.rds")




# ************

# T MAX AND MIN

# yr = 1979

walk(1979:2020, function(yr){
  
  print(yr)
  
  file.copy(str_glue("{dir_bucket_risk}/Reanalysis_data/ERA5_raw_data/daily/minimum_temperature/daily_min_temperature_{yr}.nc"),
            dir_disk)
  
  read_ncdf(str_glue("{dir_disk}/daily_min_temperature_{yr}.nc")) -> s_daily
  
  seq(as_date(str_glue("{yr}-01-01")),
      as_date(str_glue("{yr}-12-31")),
      by = "1 day") -> date_vector
  
  s_daily %>% 
    st_set_dimensions("time", values = date_vector) %>% 
    st_set_dimensions(1, names = "lon") %>% 
    st_set_dimensions(2, names = "lat") -> s_daily
  
  s_daily %>% 
    mutate(t2m = units::set_units(t2m, degC)) %>% 
    aggregate(by = "months", mean, na.rm = T) %>% 
    aperm(c(2,3,1)) -> s_monthly
  
  saveRDS(s_monthly,
          str_glue("{dir_bucket_mine}/era/monthly/2m_temperature_min/era5_monthly_mean_daily_mintemp_{yr}.rds"))
  
  unlink(str_glue("{dir_disk}/daily_min_temperature_{yr}.nc"))
  
  
})


