
J <- seq(as_date(15), as_date(364-15), by = "1 month") %>% 
  as.integer()

days_mth <- seq(as_date(0), as_date(364), by = "1 month") %>% 
  days_in_month() %>% 
  unname()

# https://edis.ifas.ufl.edu/pdf%5CAE%5CAE459%5CAE459-4802269.pdf

penman_mine <- function(Tmean, Rs, u2, Tdew, z, lat, co2adj = F, model_i = NA){
  
  # STEP 1: Mean daily temp
  
  # Tmean 
  
  # ***********************************************
  
  # STEP 2: Mean daily solar radiation
  
  # Rs
  
  # ***********************************************
  
  # STEP 3: Wind speed
  
  # u2
  
  # ***********************************************
  
  # ***********************************************
  
  # STEP 10: Mean saturation vapor pressure
  
  # original approach:
  # eTmax <-
  #   0.6108 * exp((17.27 * Tmax) / (Tmax + 237.3))
  # eTmin <-   
  #   0.6108 * exp((17.27 * Tmin) / (Tmin + 237.3))
  # 
  # es <- 
  #   (eTmax + eTmin) / 2
  
  # Bliss' approach:
  es <- 0.618 * exp((17.27 * Tmean) / (Tmean + 237.3))
  
  # ***********************************************
  
  # STEP 4: Slope of saturation vapor pressure curve
  
  Delta <- 
    4098 * es /
    (Tmean + 237.3)^2
  
  # ***********************************************
  
  # STEP 5: Atmospheric pressure
  # (available in ERA5)
  
  P <-
    101.3 * ((293 - 0.0065 * z) / 293)^5.26
  
  # ***********************************************
  
  # STEP 6: Psychometric constant
  
  gamma <- 
    1.013e-3 * P / (0.622 * 2.45)
  
  # ***********************************************
  
  # STEP 7: Delta term
  
  # DT <-
  #   Delta / 
  #   (Delta + gamma * (1 + 0.34 * u2))
  
  # ***********************************************
  
  # STEP 8: Psi term
  
  # PT <- 
  #   gamma / 
  #   (Delta + gamma * (1 + 0.34 * u2))
  
  # ***********************************************
  
  # STEP 9: Temperature term
  
  # TT <- 
  #   900 / (Tmean + 273) * u2
  
  # ***********************************************
  
  # STEP 11: Actual vapor pressure
  
  ea <- 
    0.6108 * exp((17.27 * Tdew) / (Tdew + 273.3))
  
  # ***********************************************
  
  # STEP 12: 
  # inverse relative distance Earth-Sun
  dr <-
    1 + 0.033 * cos(2 * pi / 365 * J)
  
  # solar declination
  delta <-
    0.409 * sin(2 * pi / 365 * J - 1.39)
  
  # ***********************************************
  
  # STEP 13: latitude in radians
  
  phi <- 
    pi / 180 * lat
  
  # ***********************************************
  
  # STEP 14: sunset hour angle
  
  omega_s <- 
    acos(-tan(phi) * tan(delta))
  
  # ***********************************************
  
  # STEP 15: Extraterrestrial radiation
  # (available in ERA5)
  
  Ra <- 
    (24*60/pi) * 0.0820 * dr * ((omega_s * sin(phi) * sin(delta)) + (cos(phi) * cos(delta) * sin(omega_s)))
  
  # ***********************************************
  
  # STEP 16: Clear sky solar radiation
  
  Rso <-
    (0.75 + 2e-10 * z) * Ra
  
  # ***********************************************
  
  # Step 17: Net solar radiation
  
  Rns <- 
    (1 - 0.23) * Rs
  
  # ***********************************************
  
  # STEP 18: Net outgoing long wave solar radiation
  
  Rnl <- 
    4.903e-9 * 
    # (((Tmax + 273.16)^4 + (Tmin + 273.16)^4) / 2) *
    (Tmean + 273.16)^4 *
    (0.34 - 0.14 * sqrt(ea)) *
    ((1.35 * Rs / Rso) - 0.35)
  
  # ***********************************************
  
  # STEP 19: Net radiation
  
  Rn <-
    Rns - Rnl

  # Rng <- 
  #   0.408 * Rn
  
  # ***********************************************
  
  # G (Bliss')
  G <- 
    0.3 * Rn
  
  # ***********************************************
  
  # FINAL STEP
  
  # ETrad <- 
  #   DT * Rng
  # 
  # ETwind <- 
  #   PT * TT * (es - ea)
  # 
  # ET0 <- ETwind + ETrad
  
  # *********************************************

  if(isFALSE(co2adj)){
    
    ET0 <- 
      ((0.408 * Delta * (Rn - G)) + (gamma * (900 / (Tmean + 273)) * u2 * (es - ea))) /
      (Delta + gamma * (1 + 0.34 * u2))
    
  } else {
    
    # with CO2 adjustment
    ET0 <- 
      (0.408 * Delta * (Rn - G)) + (gamma * (900 / (Tmean + 273)) * u2 * (es - ea)) /
      (Delta + gamma * (1 + u2 * (0.34 + 2.4e-4 * (co2_ppm - 300))))
    
  }
  
  if(is.na(model_i) | model_i != "Had"){
    
    ET0 <- ET0 * days_mth
    
  } else {
    
    ET0 <- ET0 * 30
    
  }
  
  ET0 <- ifelse(ET0 < 0, 0, ET0)
  
  return(ET0)
  
}


