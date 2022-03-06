
# to calculate pdsi ----
func_pdsi <- function(s_pr, s_pet){
  
  c(s_pr, s_pet, along = 3) %>% 
    st_apply(c(1,2), function(x){
      
      scPDSI::pdsi(P = x[1:504], 
                   PE = x[505:1008], 
                   sc = T)$X %>% as.vector()
      
    },
    FUTURE = T,
    future.seed = NULL,
    .fname = "time") %>% 
    
    st_set_dimensions("time", values = st_get_dimension_values(s_pr, "time")) %>% 
    aperm(c(2,3,1)) %>% 
    setNames("pdsi") -> s_pdsi
  
  return(s_pdsi)
  
}



# to map pdsi on a given date ----
func_map_date <- function(map_1, date){
  
  ggplot() +
    geom_stars(data = map_1 %>% 
                 filter(time == date) %>% 
                 adrop() %>% 
                 mutate(pdsi = case_when(pdsi < -5 ~ -5,
                                         pdsi > 5 ~ 5,
                                         TRUE ~ pdsi))) +
    geom_polygon(data = land %>% st_crop(land_rast) %>% as("Spatial") %>% fortify(), 
                 fill = NA,
                 color = "black",
                 aes(x = long, y = lat, group = group)) +
    scale_fill_binned_divergingx(palette = "RdYlBu", 
                                 rev = F, 
                                 na.value = "transparent",
                                 name = "PDSI",
                                 # breaks = c(-5, -2.5, 0, 2.5, 5),
                                 # labels = c("< -5", "-2.5", "0", "2.5", "> 5"
                                 limits = c(-5,5),
                                 n.breaks = 9
                                     
    ) +
    coord_fixed() +
    labs(subtitle = date,
         x = NULL,
         y = NULL) +
    guides(fill = guide_colorsteps(barwidth = 15, barheight = 0.5)) +
    theme(legend.position = "bottom")
  
}



# to map correlation between at a pixel level (ts vs ts) ----
func_t_cor_map <- function(map_1, map_2){
  
  c(map_1, map_2, along = 3) %>% 
    st_apply(c(1,2), function(x){
      
      cor(x[1:504], x[505:1008])
      
    }, 
    FUTURE = T)}

func_t_cor_map_plot <- function(mapp){
  
  ggplot() +
    geom_stars(data = mapp) +
    geom_polygon(data = land %>% st_crop(land_rast) %>% as("Spatial") %>% fortify(), 
                 fill = NA,
                 color = "black",
                 aes(x = long, y = lat, group = group)) +
    scale_fill_binned_diverging(palette = "Vik", 
                                rev = T,
                                na.value = "transparent",
                                limits = c(-1,1),
                                # breaks = c(-0.9, 0, 0.9),
                                # labels = c("-0.9", "0", "0.9"),
                                name = "r",
                                n.breaks = 11) +
    coord_fixed() +
    
    guides(fill = guide_colorsteps(barwidth = 15, barheight = 0.5)) +
    theme(legend.position = "bottom") +
    labs(x = NULL,
         y = NULL)
  
}



# to compare 2 per-pixel time series with a correlation threshold ----
func_ts_comparison_thres <- function(cor_map, map_1, src_name, thres) {
  
  cor_map %>% 
    as_tibble() %>%
    rename("X" = 3) %>% 
    filter(near(X, thres, 0.01)) %>% 
    slice_sample(n = 1) %>% 
    {c(.$longitude, .$latitude)} -> coords
  
  map_1 %>%
    as_tibble() %>% 
    filter(longitude == coords[1],
           latitude == coords[2]) %>% 
    pull(4) -> vect_1
  
  v_pdsi %>% 
    as_tibble() %>% 
    filter(longitude == coords[1],
           latitude == coords[2]) %>% 
    pull(4) -> vect_vds
  
  tibble(vdS = vect_vds,
         src = vect_1,
         time = st_get_dimension_values(v_pdsi, "time")) %>% 
    rename({{src_name}} := src) %>% 
    pivot_longer(-time, names_to = "src", values_to = "pdsi") %>% 
    
    ggplot(aes(x = time, y = pdsi, color = src)) +
    geom_hline(yintercept = 0, linetype = "2222") +
    geom_line(size = 0.6) +
    scale_color_discrete_qualitative("dark3", name = NULL) +
    theme(legend.position = "bottom") +
    labs(x = NULL)
  
}



# to compute correlations spatially ----
func_sp_correlation <- function(map_1){
  
  c(map_1, v_pdsi) %>% 
    setNames(c("map_1", "vdS")) %>% 
    as_tibble() %>% 
    group_by(time) %>% 
    summarize(cor = cor(map_1, vdS, use = "complete.obs")) -> t_cor 
  
  t_cor %>%   
    ggplot(aes(x = time, y = cor)) +
    geom_point(size = 0.7) +
    geom_smooth() +
    labs(y = "r",
         x = NULL)
  
}



# to plot two maps given a sp correlation thres ----
func_sp_comparison_maps <- function(map_1, src_name, thres){
  
  c(map_1, v_pdsi) %>% 
    setNames(c("map_1", "vdS")) %>% 
    as_tibble() %>% 
    group_by(time) %>% 
    summarize(cor = cor(map_1, vdS, use = "complete.obs")) -> t_cor 
  
  t_cor %>% 
    as_tibble() %>% 
    filter(near(cor, thres, 0.01)) %>% 
    slice_sample(n = 1) %>% 
    pull(time) -> t_time
  
  map_1 %>%
    filter(time == t_time) %>% 
    adrop() %>% 
    as_tibble() -> m_1
  
  v_pdsi %>%
    filter(time == t_time) %>% 
    adrop() %>% 
    as_tibble() -> m_2
  
  left_join(m_1, m_2, by = c("longitude", "latitude")) %>%
    rename({{src_name}} := 3,
           vdS = 4) %>% 
    pivot_longer(3:4, names_to = "src", values_to = "pdsi") %>% 
    mutate(pdsi = case_when(pdsi < -5 ~ -5,
                            pdsi > 5 ~ 5,
                            TRUE ~ pdsi)) %>% 
    
    {
      ggplot() +
        geom_raster(data = ., aes(x = longitude, y = latitude, fill = pdsi)) +
        geom_polygon(data = land %>% st_crop(land_rast) %>% as("Spatial") %>% fortify(), 
                     fill = NA,
                     color = "black",
                     aes(x = long, y = lat, group = group)) +
        scale_fill_binned_divergingx(palette = "RdYlBu", 
                                     rev = F, 
                                     na.value = "transparent",
                                     name = "PDSI",
                                     # breaks = c(-5, -2.5, 0, 2.5, 5),
                                     # labels = c("< -5", "-2.5", "0", "2.5", "> 5"),
                                     limits = c(-5,5),
                                     n.breaks = 9) +
        facet_wrap(~src, ncol = 2) +
        coord_equal() +
        labs(subtitle = t_time,
             x = NULL,
             y = NULL) +
        guides(fill = guide_colorsteps(barwidth = 15, barheight = 0.5)) +
        theme(legend.position = "bottom")
    }
}

ggplot() +
  geom_polygon(data = land %>% st_crop(land_rast) %>% as("Spatial") %>% fortify(), 
               fill = NA,
               color = "black",
               aes(x = long, y = lat, group = group))

# to compare +2 per-pixel randomly chosen time series ----
func_random_ts_comparison <- function(vars){

  all_pdsi %>%
    select(vars) %>%
    as_tibble() %>%
    filter(!is.na(!!sym(vars[1]))) %>%
    group_by(longitude, latitude) %>%
    mutate(gp = cur_group_id()) %>%
    ungroup() %>%
    filter(gp == sample(max(gp), 1)) -> tb

  tibble(ds1 = seq_along(vars)) %>%
    expand(ds1, ds2 = ds1) %>% 
    filter(ds1 < ds2) -> ds
  
  pmap(ds, function(ds1, ds2){
    
    v = str_glue("{vars[ds1]} vs {vars[ds2]}")
    cor = cor(tb[[vars[ds1]]], tb[[vars[ds2]]]) %>% 
      round(3)
    
    str_glue("{v} = {cor}")

  }) -> cor_coefficients
  
  tb %>%
    pivot_longer(vars, names_to = "model", values_to = "pdsi") %>%

    ggplot(aes(x = time, y = pdsi, color = model)) +
    geom_hline(yintercept = 0, linetype = "2222") +
    geom_line(size = 0.8, show.legend = F) +
    scale_color_discrete_qualitative("Dark3") +#, name = NULL) +
    # theme(legend.position = "bottom")
    labs(caption = str_glue("Correlation coefficients:\n
                            {str_c(cor_coefficients) %>% str_flatten(collapse = '\n')}"),
         x = NULL) +
    facet_wrap(~model, ncol = 1)

}



# to compare spatially +2 datasets ----
func_sp_correlation_over2 <- function(vars){
  
  tibble(ds1 = seq_along(vars)) %>%
    expand(ds1, ds2 = ds1) %>% 
    filter(ds1 < ds2) -> ds
  
  pmap_dfr(ds, function(ds1, ds2){
    
    all_pdsi %>%
      select(vars) %>%
      as_tibble() %>%
      group_by(time) %>% 
      summarize(r = cor(!!sym(vars[ds1]), !!sym(vars[ds2]), use = "complete.obs")) %>% 
      mutate(comp = str_glue("{vars[ds1]} vs {vars[ds2]}"))
    
  }) -> t_cor
  
  t_cor %>%   
    ggplot(aes(x = time, y = r, color = comp)) +
    geom_point(show.legend = F, size = 0.6) +
    geom_smooth(color = "black", linetype = "2222", size = 0.6, se = F) +
    labs(y = "r",
         x = NULL) +
    facet_wrap(~comp, ncol = 1)
  
}



# ----
func_random_sp_comparison <- function(vars, vars_comp){
  
  sample(st_get_dimension_values(all_pdsi, "time"), 1) -> d

  all_pdsi %>%
    filter(time == d) %>% 
    adrop() -> ss
  
  tibble(ds1 = seq_along(vars_comp)) %>%
    expand(ds1, ds2 = ds1) %>%
    filter(ds1 < ds2) -> ds

  ss %>%
    select(vars_comp) %>%
    as_tibble() -> tb

  pmap(ds, function(ds1, ds2){

    v = str_glue("{vars_comp[ds1]} vs {vars_comp[ds2]}")
    cor = cor(tb[[vars_comp[ds1]]], tb[[vars_comp[ds2]]], use = "complete.obs") %>%
      round(3)

    str_glue("{v} = {cor}")

  }) -> cor_coefficients
  
  ss %>% 
    select(vars) %>% 
    merge(name = "model") %>%
    mutate(X = case_when(X < -5 ~ -5,
                         X > 5 ~ 5,
                         TRUE ~ X)) %>% 
    
    {
      
      ggplot() +
        geom_stars(data = .) +
        geom_polygon(data = land %>% st_crop(land_rast) %>% as("Spatial") %>% fortify(), 
                     fill = NA,
                     color = "black",
                     aes(x = long, y = lat, group = group)) +
        scale_fill_binned_divergingx(palette = "RdYlBu", 
                                     rev = F, 
                                     na.value = "transparent",
                                     name = "PDSI",
                                     limits = c(-5, 5),
                                     n.breaks = 7
                                         
        ) +
        coord_fixed() +
        labs(subtitle = d,
             
             caption = str_glue("Correlation coefficients:\n
                            {str_c(cor_coefficients) %>% str_flatten(collapse = '\n')}"),
             
             x = NULL,
             y = NULL) +
        
        facet_wrap(~model, ncol = 3) +
        theme(legend.position = "bottom") +
        guides(fill = guide_colorsteps(barwidth = 15, barheight = 0.5))
      
    }
}


