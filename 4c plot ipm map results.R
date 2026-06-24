library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(parallel)
library(tidyr)
library(ggpubr)
library(progress)
library(Rage)
library(pbapply)
library(stars)

fn_scenarios = dir(path='output_data',pattern='ipm_grid*',full.names = TRUE)
fn_scenarios = fn_scenarios[grepl('csv', fn_scenarios)]

# load in scenarios
lapply(fn_scenarios, function(x) {
  name_final = gsub('.csv','',basename(x),fixed=TRUE)
  z = read.csv(x)
  assign(name_final, z, envir=globalenv())
  return(NULL)
})

# load in raster
r_elev_proj = rast('output_data/r_elevation_aggregated.tif')





# define plotting code
map_scenario <- function(s_this, s_base, variable, name_this)
{
  r_this = make_raster(s_this, variable)
  r_base = make_raster(s_base, variable)
  
  val_max = 100
  
  g = ggplot() +
    theme_void() +
    geom_stars(data=st_as_stars(100*(r_this - r_base)/r_base)) +
    scale_fill_gradient2(limits=c(-val_max, val_max),name='Percent change from baseline (%)') +
    coord_sf() +
    theme(axis.title = theme_bw()$axis.title) +
    xlab(variable) + ylab('')
  
  values = (100*(r_this - r_base)/r_base)[]
  
  print(sprintf('%s - %s: %.0f +/- %.0f percent', variable, name_this, median(values, na.rm=TRUE), abs(diff(quantile(values, c(0.25,0.75), na.rm=TRUE)))))
  
  return(g)
}

map_scenario_set <- function(s_this, s_base, name_this)
{
  m1 = map_scenario(s_this=s_this, 
                    s_base=s_base, 
                    variable='n_A_final', 
                    name_this=name_this)
  
  m2 = map_scenario(s_this=s_this, 
                    s_base=s_base, 
                    variable='basal_area_density_final', 
                    name_this=name_this)
  
  m3 = map_scenario(s_this=s_this, 
                    s_base=s_base, 
                    variable='n_S_final', 
                    name_this=name_this)
  
  m4 = map_scenario(s_this=s_base, 
                    s_base=s_this, 
                    variable='longevity_90_final', 
                    name_this=name_this)  
  
  m_all = list(m1, m2, m3, m4)
  
  return(m_all)
}


make_raster <- function(result, variable)
{
  pixel_means = result %>%
    group_by(pixel) %>%
    summarize(mean=mean(.data[[variable]]))
  
  r_this = r_elev_proj
  r_this[] = NA
  r_this[pixel_means$pixel] = pixel_means$mean
  names(r_this) = variable
  
  return(r_this)
}

s1 = map_scenario_set(s_this=ipm_gridded_n_S_high, 
                      s_base=ipm_gridded_base, 
                      name_this='n_S_high')

# s2 = map_scenario_set(s_this=result_scenario_tmax_high, 
#              s_base=result_scenario_base, 
#              name_this='tmax_high')
# 
# s3 = map_scenario_set(s_this=result_scenario_SWE_low, 
#                  s_base=result_scenario_base, 
#                  name_this='SWE_low')

s4 = map_scenario_set(s_this=ipm_gridded_climate_stress, 
                      s_base=ipm_gridded_base, 
                      name_this='climate_stress')

g_scenarios = ggarrange(plotlist=c(s1, 
                                   #s2, 
                                   #s3, 
                                   s4), nrow=2,ncol=4,
                        align='hv', 
                        common.legend = TRUE,
                        legend='bottom',
                        hjust = c(0,0,0),
                        labels=c('(A) Herbivore management - 200% n_S*',rep('',3), '(B) Climate change - 75% SWE, 110% Tmax',rep('',3)))

ggsave(g_scenarios, file='output_figures/g_scenarios.pdf',width=11,height=7)
ggsave(g_scenarios, file='output_figures/g_scenarios.png',width=11,height=7)



# pull in rehfeldt model
r_rehfeldt = rast('data/rasters/POTRmodel/POTR5_CGH_GMUG_change.tif')
r_rehfeldt_projected = project(r_rehfeldt, r_elev_proj, method='mode')


r_pct_change_climate_stress = (make_raster(ipm_gridded_climate_stress,'basal_area_density_final') -
  make_raster(ipm_gridded_base,'basal_area_density_final')) / make_raster(ipm_gridded_base,'basal_area_density_final')

r_classsifed_climate_stress = r_pct_change_climate_stress
r_classsifed_climate_stress[abs(as.numeric(r_classsifed_climate_stress[]))<0.1] = 0
r_classsifed_climate_stress[] = sign(r_classsifed_climate_stress[])
r_classsifed_climate_stress[is.na(r_classsifed_climate_stress)] = -2
r_classsifed_climate_stress = r_classsifed_climate_stress + 2
levels(r_classsifed_climate_stress) = c('not present','decreasing','stable','increasing')

table_xtabs = table(levels(r_classsifed_climate_stress)[[1]]$category[1+as.numeric(r_classsifed_climate_stress[])], 
      levels(r_rehfeldt_projected)[[1]]$category[1+as.numeric(r_rehfeldt_projected[])])
table_xtabs

g_ipm = ggplot() + 
  theme_void() +
  geom_stars(data=st_as_stars(r_classsifed_climate_stress)) +
  coord_sf() +
  scale_fill_manual(values=c('gray','red','black','blue'),name='Class')

g_rehfeldt = ggplot() + 
  theme_void() +
  geom_stars(data=st_as_stars(r_rehfeldt_projected)) +
  coord_sf() +
  scale_fill_manual(values=c('gray','red','blue','lightgray','darkgray'),name='Class')

ggarrange(g_rehfeldt, g_ipm, labels='auto', align='hv')
