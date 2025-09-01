library(terra)
library(dplyr)
library(ggplot2)
# library(RStoolbox)
# library(GGally)
library(ipmr)
library(pbapply)

NUM_PIXELS = NULL # debug amount, set to NULL to do all pixels
NUM_CORES = 1 # how many cores to use

load('output_data/workspace for ipm.Rdata')
source('ipm parameters.R')

r_sex_proj = rast('output_data/r_sex_proj.tif') %>% toMemory
r_ploidy_level_proj = rast('output_data/r_ploidy_level_proj.tif') %>% toMemory
r_cos_aspect_proj = rast('output_data/r_cos_aspect_proj.tif') %>% toMemory

r_medium_trees_constant_proj = rast('output_data/r_medium_trees_constant_proj.tif') %>% toMemory
r_medium_trees_midelevation_proj = rast('output_data/r_medium_trees_midelevation_proj.tif') %>% toMemory


make_ipm_spatial_scenario <- function(r_medium_trees_scenario)
{
  r_predictors = c(r_sex_proj, r_ploidy_level_proj, r_cos_aspect_proj, r_medium_trees_scenario)
  names(r_predictors) = c('geneticSexID','Ploidy_level','Cos.aspect','n_medium_trees')
  
  
  df_predictors = r_predictors[] %>%
    as.data.frame %>%
    mutate(skip = apply(.,1,function(x) {any(is.na(x))})) %>%
    mutate(lambda=NA)
  
  indices_predictors_non_na = which(df_predictors$skip==FALSE)
  
  # potentially subset to some pixels
  if (!is.null(NUM_PIXELS))
  {
    indices_predictors_non_na_ss = sample(indices_predictors_non_na, NUM_PIXELS)
  } else
  {
    indices_predictors_non_na_ss = indices_predictors_non_na
  }
  
  print(sprintf('%d total pixels', length(indices_predictors_non_na_ss)))
  
  lambdas = pbsapply(indices_predictors_non_na_ss, function(i) {
    ipm_this = make_ipm_for_site(
      m_survival = m_survival_all,
      m_growth = m_growth_all,
      m_recruit = m_recruit_all,
      other_vars=list(geneticSexIDM=paste(df_predictors$geneticSexID[i]),
                      Ploidy_levelTriploid=paste(df_predictors$Ploidy_level[i]),
                      n_medium_trees=df_predictors$n_medium_trees[i],
                      Cos.aspect=df_predictors$Cos.aspect[i])
    )
    
    if (!is.null(ipm_this))
    {
      lambda = mean(tail(as.numeric(lambda(ipm_this)),5)) # take the final 5 lambda values and average
    }
    else
    {
      lambda = NA
    }
    return(lambda)
  }, cl = NUM_CORES)
  
  df_predictors$lambda[indices_predictors_non_na_ss] = lambdas
  
  r_lambda = r_predictors[[1]]
  r_lambda[] = df_predictors$lambda
  
  return(r_lambda)
}

r_lambda_n_medium_constant_5 = make_ipm_spatial_scenario(r_medium_trees_scenario=r_medium_trees_constant_proj)  
r_lambda_n_medium_midelevation = make_ipm_spatial_scenario(r_medium_trees_scenario=r_medium_trees_midelevation_proj)  


writeRaster(r_lambda_n_medium_constant_5, 
            file='output_data/r_lambda_n_medium_constant_5.tif',
            overwrite=TRUE)
writeRaster(r_lambda_n_medium_midelevation, 
            file='output_data/r_lambda_n_medium_midelevation.tif',
            overwrite=TRUE)


# currently this has no bootstrapping to get the uncertainty in dataset and/or sex propagated. 
# i think it is too computationally intensive to do this.






# 
# 
# plot_map <- function(r_all, include_sex, df_sites_ss)
# {  
#   names(r_all) = c('geneticSexID', 'Ploidy_level', 'Cos.aspect', 'Slope', 'Elevation', 'n_medium_trees')
#   predictors_all = as.data.frame(r_all[])
#   
#   # pick some pseudoabsence coordinates that are not in aspen but are in terrain
#   # candidate_indices = which(is.na(r_ploidy_level_proj[]))
#   # predictors_ss_random = predictors_all[sample(candidate_indices,200),]
#   # predictors_ss_random$lambda = 0.9 # no aspen here, but the size of this value may reall
#   
#   #
#   # df_for_training = rbind(predictors_ss_random,
#   #                         df_sites_for_ipm_joined %>% select(all_of(names(predictors_ss_random)))) %>%
#   #   na.omit
#   df_for_training = df_sites_ss %>%
#     mutate(geneticSexID = as.numeric(geneticSexID=='M')) %>%# set 0=female,1=male
#     filter(!is.na(lambda)) # deal with ensemble missing values
#    
#   # make model to interpolate the computationally slow IPM and for unknown values of the other variables
#   # might be better to actually run the IPM drawing random values for the other variables, then kriging
#   # to explicitly account for assumed uncertainty
#   # really would need to see how the predictors like num_medium covary with other predictors, but this in turn
#   # requires knowing about other vegetation dynamics and how seedling recruitment and size dependent competition play out- 
#   # too many uncertainties to probably be realistically useful
#   # at least we know about the sex effect
#   if (include_sex==TRUE)
#   {
#     rf_remotely_sensed = randomForest(lambda ~ Cos.aspect + Slope + Elevation +
#                                         n_medium_trees + 
#                                         Ploidy_level + geneticSexID, 
#                                       # train on real data + pseudo absence
#                                       data=df_for_training)
#   }
#   else
#   {
#     rf_remotely_sensed = randomForest(lambda ~ Cos.aspect + Slope + Elevation +
#                                         n_medium_trees + 
#                                         Ploidy_level, 
#                                       # train on real data + pseudo absence
#                                       data=df_for_training)    
#   }
#   randomForest::importance(rf_remotely_sensed)
#   rf_remotely_sensed
#   
#   # concerned the ploidy_level response will be incorrect if we extrapolate with a nonlinear model
#   
#   lambda_predicted = predict(rf_remotely_sensed, newdata=predictors_all)
#   r_lambda = r_cos_aspect_proj
#   r_lambda[] = lambda_predicted
#   
#   candidate_indices = which(is.na(r_ploidy_level_proj[]))
#   r_lambda[candidate_indices] = NA
#   
#   # pdf(file='figures/g_ipm_map_projected.pdf',width=10,height=10)
#   # plot(r_lambda, maxcell=1e7)
#   # dev.off()
#   
#   return(r_lambda)
# }
# 
# df_sites_split = df_sites_for_ipm_joined %>%
#   group_by(replicate) %>%
#   group_split
# 
# r_scenario_medium_trees_constant_bootstrap = lapply(1:length(df_sites_split), function(i) {
#   cat('.')
#   r_scenario = plot_map(r_all = c(r_sex_proj, 
#                      r_ploidy_level_proj, 
#                      r_cos_aspect_proj, 
#                      r_slope_proj, 
#                      r_elevation, 
#                      r_medium_trees_scenario_constant), 
#            include_sex = TRUE,
#            df_sites_ss = df_sites_split[[i]])
#   return(r_scenario)
#   })
# r_scenario_medium_trees_constant_bootstrap_stack = do.call("c", r_scenario_medium_trees_constant_bootstrap)
# r_scenario_medium_trees_constant_bootstrap_stack_mean = mean(r_scenario_medium_trees_constant_bootstrap_stack)
# r_scenario_medium_trees_constant_bootstrap_stack_sd = stdev(r_scenario_medium_trees_constant_bootstrap_stack, na.rm=TRUE)
# 
# r_scenario_medium_trees_midelevation_bootstrap = lapply(1:length(df_sites_split), function(i) {
#   cat('.')
#   r_scenario = plot_map(r_all = c(r_sex_proj, 
#                                   r_ploidy_level_proj, 
#                                   r_cos_aspect_proj, 
#                                   r_slope_proj, 
#                                   r_elevation, 
#                                   r_medium_trees_scenario_midelevation), 
#                         include_sex = TRUE,
#                         df_sites_ss = df_sites_split[[i]])
#   return(r_scenario)
# })
# r_scenario_medium_trees_midelevation_bootstrap_stack = do.call("c", r_scenario_medium_trees_midelevation_bootstrap)
# r_scenario_medium_trees_midelevation_bootstrap_stack_mean = mean(r_scenario_medium_trees_midelevation_bootstrap_stack)
# r_scenario_medium_trees_midelevation_bootstrap_stack_sd = stdev(r_scenario_medium_trees_midelevation_bootstrap_stack, na.rm=TRUE)
# 
# 
# 
# 
# g_lambda_projected_scenario_medium_trees_constant_bootstrap_stack_mean = ggR(r_scenario_medium_trees_constant_bootstrap_stack_mean, geom_raster = TRUE, maxpixels=1e7) + 
#   scale_fill_gradient2(low='red',high='blue',mid='gray',midpoint=1,name='lambda') +
#   theme_bw()
# g_lambda_projected_scenario_medium_trees_constant_bootstrap_stack_sd = ggR(r_scenario_medium_trees_constant_bootstrap_stack_sd, geom_raster = TRUE, maxpixels=1e7) + 
#   scale_fill_viridis_c(name='lambda') +
#   theme_bw()
# ggsave(g_lambda_projected_scenario_medium_trees_constant_bootstrap_stack_mean, file=sprintf('output_figures/g_lambda_projected_scenario_medium_trees_constant_bootstrap_stack_mean.png'),width=10,height=10)
# ggsave(g_lambda_projected_scenario_medium_trees_constant_bootstrap_stack_sd, file=sprintf('output_figures/g_lambda_projected_scenario_medium_trees_constant_bootstrap_stack_sd.png'),width=10,height=10)
# 
# g_lambda_projected_scenario_medium_trees_midelevation_bootstrap_stack_mean = ggR(r_scenario_medium_trees_midelevation_bootstrap_stack_mean, geom_raster = TRUE, maxpixels=1e7) + 
#   scale_fill_gradient2(low='red',high='blue',mid='gray',midpoint=1,name='lambda') +
#   theme_bw()
# g_lambda_projected_scenario_medium_trees_midelevation_bootstrap_stack_sd = ggR(r_scenario_medium_trees_midelevation_bootstrap_stack_sd, geom_raster = TRUE, maxpixels=1e7) + 
#   scale_fill_viridis_c(name='lambda') +
#   theme_bw()
# ggsave(g_lambda_projected_scenario_medium_trees_midelevation_bootstrap_stack_mean, file=sprintf('output_figures/g_lambda_projected_scenario_medium_trees_midelevation_bootstrap_stack_mean.png'),width=10,height=10)
# ggsave(g_lambda_projected_scenario_medium_trees_midelevation_bootstrap_stack_sd, file=sprintf('output_figures/g_lambda_projected_scenario_medium_trees_midelevation_bootstrap_stack_sd.png'),width=10,height=10)
# 
# 
# 