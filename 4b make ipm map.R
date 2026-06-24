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


load('output_data/workspace for ipm.Rdata')
source('get_climate.R')

NUM_CORES = 8 # as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE')) # how many cores to use
NUM_CLIMATE_REPLICATES = 1 # number of pulls of climate
NUM_RESAMPLES = 3

r_sex_proj = rast('output_data/r_sex_projected_masked.tif') 
r_ploidy_level_proj = rast('output_data/r_ploidy_level_projected_masked.tif') 
r_cos_aspect_proj = rast('output_data/r_cos_aspect_aggregated.tif') 
r_elev_proj = rast('output_data/r_elevation_aggregated.tif')

r_predictors = c(r_sex_proj, r_ploidy_level_proj, r_cos_aspect_proj, r_elev_proj)
names(r_predictors) = c('geneticSexID','Ploidy_level','Cos.aspect','Elevation')


df_predictors = r_predictors[] %>%
  as.data.frame %>%
  mutate(skip = apply(.,1,function(x) {any(is.na(x))}))

coords_utm = crds(r_predictors, na.rm=FALSE, df=TRUE)
names(coords_utm) = c('X.UTM','Y.UTM')
coords_utm_sf = st_as_sf(coords_utm, coords=c('X.UTM','Y.UTM'), crs=32613)
coords_lat_lon_sf = coords_utm_sf %>% 
  st_transform(crs=4326) %>%
  st_coordinates %>%
  as.data.frame
names(coords_lat_lon_sf) = c('Longitude','Latitude')

df_predictors = cbind(df_predictors, coords_lat_lon_sf) %>%
  as.data.frame %>%
  mutate(pixel=row_number())

# check combinations
table(df_predictors$Ploidy_level, df_predictors$geneticSexID) 
table(df_predictors$skip) 
 
indices_predictors_non_na = which(df_predictors$skip==FALSE)







# run across all sites

df_sites_for_ipm_base = df_predictors %>%
  filter(skip==FALSE) %>%
  # rename a few things to standardize with the IPM code
  mutate(n_iterations=300,
         # assume a standard population at each location since we don't have data
         population_density_initial = 0.1, 
         size_mean = 20, 
         size_sd = 5) %>%
  mutate(cores=1, # do the looping over pixels not cores here  
         progress=FALSE) %>%
  # rename columns for the script
  mutate(geneticSexIDmale=as.character(as.numeric((geneticSexID=='male'))),
         Ploidy_leveltriploid=as.character(as.numeric((Ploidy_level=='triploid'))),
         geneticSexIDunknown=as.character(as.numeric((geneticSexID=='unknown'))),
         Ploidy_levelunknown=as.character(as.numeric((Ploidy_level=='unknown')))) %>%
  select(-geneticSexID, -Ploidy_level, -skip) %>%
  group_by(pixel) %>%
  slice(rep(1:n(), each=NUM_CLIMATE_REPLICATES)) %>% 
  mutate(prefix_this=paste('map', pixel, row_number(), sep='_'))

#n_iterations=10

run_scenario_spatial <- function(df, scenario_name, n_S_scale_factor=1, Tmax_scale_factor=1, SWE_scale_factor=1)
{
  set.seed(1)

  climate_all_this = lapply(1:nrow(df), function(i)
  {
    cat(sprintf('climate %d %.3f\n',i, i/nrow(df)))
    lat_this = df$Latitude[i]
    lon_this = df$Longitude[i]
    climate_this = make_climate_ts_at_location(lon=lon_this, lat=lat_this, 
                                               num_time_points = df$n_iterations[i] + 1)
    return(climate_this)
  })
  
  ipm_all_this = pblapply(1:nrow(df), function(i)
  {
    cat(sprintf('run %d %.3f\n',i, i/nrow(df)))
    params_this = as.list(df[i,])
    # make environment sequence
    climate_this = climate_all_this[[i]]
    params_this$SWE.sequence=climate_this$SWE * SWE_scale_factor
    params_this$Tmax.sequence=climate_this$Tmax * Tmax_scale_factor
    params_this$n_S_scale_factor = n_S_scale_factor
    params_this$Latitude = NULL
    params_this$Longitude = NULL
    params_this$pixel = NULL
    
    result = do.call("run_model_resampled",params_this)
    
    return(result)
  }, cl=NUM_CORES) # multiprocess over pixels
  
  ipm_all_this_df = do.call('rbind',lapply(1:length(ipm_all_this), function(r) { do.call('rbind',lapply(ipm_all_this[[r]], function(x) { as.data.frame(x$results) %>% mutate(index=r) })) })) %>%
    mutate(pixel=as.numeric(sapply(strsplit(prefix_this,split='_'),function(x) {x[2]}))) %>%
    mutate(replicate_climate_resample=as.numeric(sapply(strsplit(prefix_this,split='_'),function(x) { x[3] }))) %>%
    mutate(replicate_data_resample=sapply(strsplit(prefix_this,split='_'),function(x) {x[4]})) %>%
    mutate(Ploidy_level=ifelse(Ploidy_leveltriploid=='1','triploid',ifelse(Ploidy_levelunknown=='1','unknown','diploid'))) %>%
    mutate(geneticSexID=ifelse(geneticSexIDmale=='1','male',ifelse(geneticSexIDunknown=='1','unknown','female')))
  
  saveRDS(ipm_all_this, file=sprintf('output_data/ipm_gridded_%s.Rdata', scenario_name))
  write.csv(ipm_all_this_df, file=sprintf('output_data/ipm_gridded_%s.csv', scenario_name))
  
  return(list(ipm_all_this_df=ipm_all_this_df,ipm_all_this=ipm_all_this))  
}



result_scenario_base = run_scenario_spatial(df_sites_for_ipm_base,scenario_name = 'base')
result_scenario_n_S_low = run_scenario_spatial(df_sites_for_ipm_base,scenario_name = 'n_S_low', n_S_scale_factor = 0.5)
result_scenario_n_S_high = run_scenario_spatial(df_sites_for_ipm_base,scenario_name = 'n_S_high', n_S_scale_factor = 2)
result_scenario_tmax_high = run_scenario_spatial(df_sites_for_ipm_base,scenario_name = 'tmax_high', Tmax_scale_factor = 1.1) # see 2°C over 20°C
result_scenario_SWE_low = run_scenario_spatial(df_sites_for_ipm_base,scenario_name = 'SWE_low', SWE_scale_factor = 0.75) # see https://www.nrcs.usda.gov/sites/default/files/2023-03/2021GL094985Future%20Changes%20in%20Snowpack%20Snowmelt%20and%20Runoff.pdf
result_scenario_climate_stress = run_scenario_spatial(df_sites_for_ipm_base,scenario_name = 'climate_stress', SWE_scale_factor = 0.75, Tmax_scale_factor = 1.1)



#save(list=ls(pattern='result*'),file='~/Downloads/scenarios.Rdata')






# 
# run_scenario_recruitment <- function(df, scenario_name, n_iterations_burn_in, n_iterations_perturbation)
# {
#   set.seed(1)
#   
#   n_S_scale_factor.sequence_this = c(rep(1, n_iterations_burn_in+1),rep(0.5, n_iterations_perturbation+1))
#   
#   climate_all_this = lapply(1:nrow(df), function(i)
#   {
#     cat(sprintf('%d %.3f\n',i, i/nrow(df)))
#     lat_this = df$Latitude[i]
#     lon_this = df$Longitude[i]
#     climate_this = make_climate_ts_at_location(lon=lon_this, lat=lat_this, 
#                                                num_time_points = n_iterations_burn_in + n_iterations_perturbation + 1)
#     return(climate_this)
#   })
#   
#   
#   ipm_all_this = lapply(1:nrow(df), function(i)
#   {
#     cat(sprintf('%d %.3f\n',i, i/nrow(df)))
#     params_this = as.list(df[i,])
#     # make environment sequence
#     climate_this = climate_all_this[[i]]
#     params_this$SWE.sequence=climate_this$SWE
#     params_this$STB.sequence=climate_this$STB
#     params_this$n_S_scale_factor.sequence = n_S_scale_factor.sequence_this
#     params_this$Latitude = NULL
#     params_this$Longitude = NULL
#     params_this$pixel = NULL
#     
#     result = do.call("run_model_resampled",params_this)
#     
#     return(result)
#   })
#   
#   ipm_all_this_df = do.call('rbind',lapply(1:length(ipm_all_this), function(r) { do.call('rbind',lapply(ipm_all_this[[r]], function(x) { as.data.frame(x$results) %>% mutate(index=r) })) })) %>%
#     mutate(pixel=as.numeric(sapply(strsplit(prefix_this,split='_'),function(x) {x[2]}))) %>%
#     mutate(replicate_climate_resample=as.numeric(sapply(strsplit(prefix_this,split='_'),function(x) { x[3] }))) %>%
#     mutate(replicate_data_resample=sapply(strsplit(prefix_this,split='_'),function(x) {x[4]})) %>%
#     mutate(Ploidy_level=ifelse(Ploidy_leveltriploid=='1','triploid',ifelse(Ploidy_levelunknown=='1','unknown','diploid'))) %>%
#     mutate(geneticSexID=ifelse(geneticSexIDmale=='1','male',ifelse(geneticSexIDunknown=='1','unknown','female')))
#   
#   saveRDS(ipm_all_this_df, file=sprintf('output_data/ipm_%s.Rdata', scenario_name))
#   
#   return(list(ipm_all_this_df=ipm_all_this_df,ipm_all_this=ipm_all_this))
# }
# 
# z = run_scenario_recruitment(df_sites_for_ipm_base[1:1,], 'n_S_low', 
#                              n_iterations_burn_in = 200, 
#                              n_iterations_perturbation = 30)
# 
# w = lapply(z$ipm_all_this, function(z_this) { do.call('cbind',lapply(z_this, function(x) { 
#   df = data.frame(t=1:20, y=tail(x$bad,20))  
#   slope = coef(lm(y~t,data=df))[2]
#   return(slope)
# }))    })
# z$ipm_all_this[[1]][[2]]$bad %>% tail(20) %>% plot
# 
# make_raster <- function(result_scenario)
# {
#   
# }






#df_predictors$lambda[indices_predictors_non_na_ss] = unlist(lambdas)
# 
# make_ipm_spatial_scenario <- function(r_medium_trees_scenario)
# {
#   lambdas = mclapply(indices_predictors_non_na_ss, function(i) {
#     print(i)
#     ipm_this = make_ipm_for_site(
#       m_survival = m_survival_all,
#       m_growth = m_growth_all,
#       m_recruit = m_recruit_all,
#       other_vars=list(geneticSexIDM=paste(df_predictors$geneticSexID[i]),
#                       Ploidy_levelTriploid=paste(df_predictors$Ploidy_level[i]),
#                       n_medium_trees=df_predictors$n_medium_trees[i],
#                       Cos.aspect=df_predictors$Cos.aspect[i],
#                       Elevation=df_predictors$Elevation[i])
#     )
#     
#     if (!is.null(ipm_this))
#     {
#       lambda = mean(tail(as.numeric(lambda(ipm_this)),5)) # take the final 5 lambda values and average
#     }
#     else
#     {
#       lambda = NA
#     }
#     return(lambda)
#   }, mc.cores = NUM_CORES)
#   
# 
#   
#   r_lambda = r_predictors[[1]]
#   r_lambda[] = df_predictors$lambda
#   
#   return(r_lambda)
# }
# 
# r_lambda_n_medium_constant_5 = make_ipm_spatial_scenario(r_medium_trees_scenario=r_medium_trees_constant_proj)  
# 
# 
# writeRaster(r_lambda_n_medium_constant_5, 
#             file='output_data/r_lambda_n_medium_constant_5.tif',
#             overwrite=TRUE)
# 
# # currently this has no bootstrapping to get the uncertainty in dataset and/or sex propagated. 
# # i think it is too computationally intensive to do this.