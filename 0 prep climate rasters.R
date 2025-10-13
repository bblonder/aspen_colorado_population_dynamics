library(dplyr)
library(ggplot2)
library(terra)
library(tidyr)

# 
# fn_rasters_spi = dir('data/rasters/',pattern='*spi90d*',full.names = TRUE)
# rasters_spi = rast(fn_rasters_spi)
# names(rasters_spi) = paste("SPI",2016:2023,sep=".")

fn_rasters_swe = dir('data/rasters/',pattern='*_SWE_*',full.names = TRUE)
rasters_swe = rast(fn_rasters_swe)
names(rasters_swe) = paste("SWE",2016:2023,sep=".")

fn_rasters_stb = dir('data/rasters/',pattern='*_short_term_blend_*',full.names = TRUE)
rasters_stb = rast(fn_rasters_stb)
names(rasters_stb) = paste("STB",2016:2023,sep=".")

data_sites = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/manuscript ploidy 2019/new phyt/SI new phyt/File S1 - aspen data site-level processed 30 Mar 2020.csv') %>%
  select(site_code=Site_Code, Latitude, Longitude)

xy = data_sites %>% 
  select(lon=Longitude, lat=Latitude) %>%
  vect(crs="+proj=longlat +datum=WGS84")

swe_extracted = terra::extract(rasters_swe, xy)

stb_extracted = terra::extract(rasters_stb, xy)

df_climate = cbind(swe_extracted, stb_extracted) %>%
  cbind(data_sites %>% select(site_code)) %>%
  select(-ID) %>%
  pivot_longer(!site_code) %>%
  mutate(year=sapply(name, function(x) { as.numeric(strsplit(x,split="\\.")[[1]][2]) })) %>%
  mutate(metric=sapply(name, function(x) { strsplit(x,split="\\.")[[1]][1] }))

g_climate = ggplot(df_climate, aes(x=year,y=value,group=site_code)) +
  geom_line() +
  facet_wrap(~metric,scales='free') +
  theme_bw()
ggsave(g_climate, file='output_figures/g_climate.pdf')


g_climate_by_year = df_climate %>% filter(metric!='SPI') %>% select(-name) %>% pivot_wider(names_from=metric) %>%
  ggplot(aes(x=SWE,y=STB)) + geom_point() +
  facet_wrap(~year) +
  theme_bw()
ggsave(g_climate_by_year, file='output_figures/g_climate_by_year.pdf')



get_lag <- function(df, year_this, delta_years, metric_this)
{
  df_result = df %>%
    filter(metric==metric_this) %>%
    filter(year==(year_this - delta_years)) %>%
    mutate(year=year_this) %>%
    select(site_code, metric, year, value)
  
  names(df_result)[length(names(df_result))] = paste(metric_this,delta_years,sep=".")
  return(df_result)
}

get_lag_for_year <- function(df, year, metric)
{
  get_lag(df, year, 0, metric) %>%
    left_join(get_lag(df, year, 1, metric), by=c('site_code','metric','year')) %>%
    left_join(get_lag(df, year, 2, metric), by=c('site_code','metric','year')) %>%
    select(-metric)
}

df_lagged_climate_stb = do.call("rbind",lapply(c(2018, 2020, 2023), get_lag_for_year, df=df_climate,metric='STB'))
df_lagged_climate_swe = do.call("rbind",lapply(c(2018, 2020, 2023), get_lag_for_year, df=df_climate,metric='SWE'))

df_lagged_climate = df_lagged_climate_stb %>%
  left_join(df_lagged_climate_swe, by=c('site_code','year'))


# write out extracts
write.csv(df_lagged_climate,file='output_data/df_lagged_climate.csv', row.names=FALSE)
write.csv(df_climate, file='output_data/df_climate.csv', row.names=FALSE)
