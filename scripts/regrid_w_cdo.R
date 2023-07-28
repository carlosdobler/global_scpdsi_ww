
# REGRID WITH CDO

library(RNetCDF)

# reference grid:

# list.files("~/bucket_mine/remo/monthly/", full.names = T) %>%
#   .[str_detect(., "regrid", negate = T)] %>%
#   .[1] %>%
#   open.nc -> nc
ff <- "tasmax_SAM-22_NCC-NorESM1-M_historical_r1i1p1_GERICS-REMO2015_v1_day_19700101-19701231.nc"
nc <- open.nc(ff)

print.nc(nc)

x <- var.get.nc(nc, "lon") %>% ifelse(. > 180, .-360, .)
y <- var.get.nc(nc, "lat")

c("gridtype = latlon",
  str_c("xfirst = ", round(min(x), 1)),
  "xinc = 0.2",
  str_c("xsize = ", round((round(max(x), 1) - round(min(x), 1))/0.2)),
  str_c("yfirst = ", round(min(y), 1)),
  "yinc = 0.2",
  str_c("ysize = ", round((round(max(y), 1) - round(min(y), 1))/0.2))
  ) %>%
  write_lines("grid.txt")

# regrid with cdo:

list.files("~/bucket_mine/remo/monthly/new", full.names = T) %>%
  .[str_detect(., "wget", negate = T)] %>%

  walk(function(f){

    f %>%
      str_split("/") %>%
      unlist() %>%
      last() %>%
      print()

    f %>%
      str_remove(".nc") %>%
      str_c("_regrid.nc") -> f_out

    system(str_glue("cdo remapbil,grid.txt {f} {f_out}"))

  })

# ******

"~/bucket_risk/RCM_raw_data/REMO2015/lm_files/lm_SAM.nc" -> dom_grid
system(str_glue("cdo remapbil,{dom_grid} {ff} {f_out}"))
