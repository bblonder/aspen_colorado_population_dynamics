library(terra)

fn_rasters_swe = dir('data/rasters/',pattern='*_SWE_*',full.names = TRUE)
rasters_swe = rast(fn_rasters_swe)
names(rasters_swe) = paste("SWE",2016:2023,sep=".")

fn_rasters_stb = dir('data/rasters/',pattern='*_short_term_blend_*',full.names = TRUE)
rasters_stb = rast(fn_rasters_stb)
names(rasters_stb) = paste("STB",2016:2023,sep=".")

fn_rasters_tmax = dir('data/rasters/',pattern='*_tmax_*',full.names = TRUE) # mean tmax
rasters_tmax = rast(fn_rasters_tmax)
names(rasters_tmax) = paste("Tmax",2016:2023,sep=".")

get_climate_ts_at_location <- function(lat, lon)
{
  stopifnot(length(lat)==1)
  p = data.frame(lon=lon,lat=lat) %>% st_as_sf(coords=c('lon','lat'),crs=4326)

  ts_swe = terra::extract(rasters_swe, p) %>% 
    select(-ID) %>%
    t %>%
    as.data.frame %>%
    rename(SWE=1)
  
  ts_stb = terra::extract(rasters_stb, p) %>% 
    select(-ID) %>%
    t %>%
    as.data.frame %>%
    rename(STB=1)
  
  ts_tmax = terra::extract(rasters_tmax, p) %>% 
    select(-ID) %>%
    t %>%
    as.data.frame %>%
    rename(Tmax=1)
  
  return(cbind(ts_swe, ts_stb, ts_tmax))
}

make_climate_ts_at_location <- function(lat, lon, num_time_points=100, weight_1=3)
{
  climate_this = get_climate_ts_at_location(lat, lon)
  
  tmpFiles(remove=TRUE)
  
  # # of years to skip along series, relative weights
  sequence = 1 + (cumsum(sample(c(1,2),size=num_time_points,replace=TRUE,prob=c(weight_1,1))) %% nrow(climate_this)) 
  
  return(climate_this[sequence,])
}


# df_climate_summary_SWE= df_climate %>% 
#                         group_by(site_code, metric) %>%
#                         filter(metric=='SWE') %>%
#   # only a few values are lost
#                         filter(value>0) %>%
#                         reframe(t(fitdistr(sqrt(value),'normal')$estimate) %>% as.data.frame) %>%
#   select(-metric) %>%
#   rename(SWE.meansqrt=mean,SWE.sdsqrt=sd)
# 
# df_climate_summary_STB = df_climate %>% 
#   group_by(site_code, metric) %>%
#   filter(metric=='STB') %>%
#   reframe(t(fitdistr(value,'normal')$estimate) %>% as.data.frame) %>%
#   select(-metric) %>%
#   rename(STB.mean=mean,STB.sd=sd)


