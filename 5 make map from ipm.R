library(randomForest)
library(terra)
library(dplyr)
library(mgcv)
library(ggplot2)
library(RStoolbox)
library(GGally)

df_sites_for_ipm_joined = read.csv('output_data/sites_data_frame.csv') %>%
  mutate(Ploidy_level = as.numeric(factor(Ploidy_level))-1) # convert to numeric, triploid = 1



# need to add other sites out of the current range that are true absences, zero lambda

r_sex = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspenSex/aspen_sex_prediction_20250110.tif')
r_ploidy_level = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/min_phase_cytotype_medfilt-seived.tif')
r_ploidy_level[r_ploidy_level==-1] = NA # all 0 and 1
r_cos_aspect = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/r_cos_aspect.tif')
r_slope = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/r_slope.tif')
r_elevation = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/altitudes.tif')

# project to elevation (5 m resolution) for now
r_sex_proj = project(r_sex-1, r_elevation) # set male=1, female=0
r_ploidy_level_proj = project(r_ploidy_level, r_elevation)
r_cos_aspect_proj = project(r_cos_aspect, r_elevation)
r_slope_proj = project(r_slope, r_elevation)


r_medium_trees_scenario_midelevation = r_cos_aspect_proj
m_nmt = nls(n_medium_trees ~ k*exp(-1/2*(Elevation-mu)^2/sigma^2),start=c(mu=3000,sigma=500,k=3) , 
            data = df_sites_for_ipm_joined)
r_medium_trees_scenario_midelevation[] = predict(m_nmt, newdata=data.frame(Elevation=as.numeric(r_elevation[])))
# 
# #r_medium_trees_scenario[!is.na(r_medium_trees_scenario)] = 0 # cosine aspect locations to zero
# r_medium_trees_scenario[!is.na(r_cos_aspect_proj)] = sample(df_sites_for_ipm_joined$n_medium_trees, 
#                                          length(which(!is.na(r_cos_aspect_proj[]))), replace=TRUE)
# 
# r_ploidy_level_proj_scenario = r_ploidy_level_proj
# m_pl = glm(Ploidy_level~Elevation,data=df_sites_for_ipm_joined, family=binomial)
# r_ploidy_level_proj_scenario[] = predict(m_pl, newdata=data.frame(Elevation=as.numeric(r_elevation[])))

r_medium_trees_scenario_constant = r_cos_aspect_proj
r_medium_trees_scenario_constant[] = median(df_sites_for_ipm_joined$n_medium_trees)

plot_map <- function(r_all, name_this, include_sex)
{  
  names(r_all) = c('geneticSexID', 'Ploidy_level', 'Cos.aspect', 'Slope', 'Elevation', 'n_medium_trees')
  predictors_all = as.data.frame(r_all[])
  
  # pick some pseudoabsence coordinates that are not in aspen but are in terrain
  # candidate_indices = which(is.na(r_ploidy_level_proj[]))
  # predictors_ss_random = predictors_all[sample(candidate_indices,200),]
  # predictors_ss_random$lambda = 0.9 # no aspen here, but the size of this value may reall
  
  #
  # df_for_training = rbind(predictors_ss_random,
  #                         df_sites_for_ipm_joined %>% select(all_of(names(predictors_ss_random)))) %>%
  #   na.omit
  df_for_training = df_sites_for_ipm_joined %>%
    mutate(geneticSexID = as.numeric(geneticSexID=='M')) # set 0=female,1=male
  
  # make model to interpolate the computationally slow IPM and for unknown values of the other variables
  # might be better to actually run the IPM drawing random values for the other variables, then kriging
  # to explicitly account for assumed uncertainty
  # really would need to see how the predictors like num_medium covary with other predictors, but this in turn
  # requires knowing about other vegetation dynamics and how seedling recruitment and size dependent competition play out- 
  # too many uncertainties to probably be realistically useful
  # at least we know about the sex effect
  if (include_sex==TRUE)
  {
    rf_remotely_sensed = randomForest(lambda ~ Cos.aspect + Slope + Elevation +
                                        n_medium_trees + 
                                        Ploidy_level + geneticSexID, 
                                      # train on real data + pseudo absence
                                      data=df_for_training)
  }
  else
  {
    rf_remotely_sensed = randomForest(lambda ~ Cos.aspect + Slope + Elevation +
                                        n_medium_trees + 
                                        Ploidy_level, 
                                      # train on real data + pseudo absence
                                      data=df_for_training)    
  }
  randomForest::importance(rf_remotely_sensed)
  rf_remotely_sensed
  
  # concerned the ploidy_level response will be incorrect if we extrapolate with a nonlinear model
  
  lambda_predicted = predict(rf_remotely_sensed, newdata=predictors_all)
  r_lambda = r_cos_aspect_proj
  r_lambda[] = lambda_predicted
  
  candidate_indices = which(is.na(r_ploidy_level_proj[]))
  r_lambda[candidate_indices] = NA
  
  # pdf(file='figures/g_ipm_map_projected.pdf',width=10,height=10)
  # plot(r_lambda, maxcell=1e7)
  # dev.off()
  
  g_lambda_projected = ggR(r_lambda, geom_raster = TRUE, maxpixels=1e7) + 
    scale_fill_gradient2(low='red',high='blue',mid='gray',midpoint=1,name='lambda') +
    theme_bw()
  ggsave(g_lambda_projected, file=sprintf('output_figures/g_ipm_lambda_projected_%s.png',name_this),width=10,height=10)
  
  # r_lambda_clipped_classified = classify(r_lambda, c(0,0.995,1.005,Inf))
  # plot(r_lambda_clipped_classified)
  
  g_lambda_projected_hist = ggplot(data.frame(lambda=as.numeric(lambda_predicted)),aes(x=lambda)) +
    geom_histogram(binwidth=0.001) +
    geom_vline(xintercept = 1,color='purple') +
    theme_bw()
  ggsave(g_lambda_projected_hist,file=sprintf('output_figures/g_ipm_lambda_projected_hist_%s.pdf',name_this))
  
  table(lambda_predicted>1)
  
  return(r_lambda)
}


r_scenario_medium_trees_constant_include_sex = plot_map(c(r_sex_proj, r_ploidy_level_proj, r_cos_aspect_proj, r_slope_proj, r_elevation, r_medium_trees_scenario_constant), 
         'medium_trees_constant_include_sex', include_sex = TRUE)
r_scenario_medium_trees_constant_exclude_sex = plot_map(c(r_sex_proj, r_ploidy_level_proj, r_cos_aspect_proj, r_slope_proj, r_elevation, r_medium_trees_scenario_constant), 
         'medium_trees_constant_exclude_sex', include_sex = FALSE)

r_scenario_medium_trees_midelevation_include_sex = plot_map(c(r_sex_proj, r_ploidy_level_proj, r_cos_aspect_proj, r_slope_proj, r_elevation, r_medium_trees_scenario_midelevation), 
         'medium_trees_midelevation_include_sex', include_sex = TRUE)
r_scenario_medium_trees_midelevation_exclude_sex = plot_map(c(r_sex_proj, r_ploidy_level_proj, r_cos_aspect_proj, r_slope_proj, r_elevation, r_medium_trees_scenario_midelevation), 
         'medium_trees_midelevation_exclude_sex', include_sex = FALSE)

r_scenarios = c(r_scenario_medium_trees_constant_include_sex, r_scenario_medium_trees_constant_exclude_sex, r_scenario_medium_trees_midelevation_include_sex, r_scenario_medium_trees_midelevation_exclude_sex)
r_scenarios_values = r_scenarios[]

r_scenarios_values_sampled = r_scenarios_values %>% 
  na.omit
r_scenarios_values_sampled = r_scenarios_values_sampled[sample(1:nrow(r_scenarios_values_sampled),1000),] %>%
  as.data.frame
names(r_scenarios_values_sampled) = c('medium_trees_constant_include_sex','medium_trees_constant_exclude_sex','medium_trees_midelevation_include_sex','medium_trees_midelevation_exclude_sex')

g_pairs = ggpairs(r_scenarios_values_sampled) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0,color='red') +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1)
ggsave(g_pairs, file='output_figures/g_ipm_lambda_pairs.pdf',width=12,height=12)
