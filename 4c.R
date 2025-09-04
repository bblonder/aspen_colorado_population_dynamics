library(terra)
library(dplyr)
library(ggplot2)
library(ipmr)
library(parallel)

NUM_PIXELS = NULL # debug amount, set to NULL to do all pixels
NUM_CORES = as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE')) # how many cores to use

load('output_data/workspace for ipm.Rdata')
source('ipm parameters.R')

r_sex_proj = rast('output_data/r_sex_proj.tif') 
r_ploidy_level_proj = rast('output_data/r_ploidy_level_proj.tif') 
r_cos_aspect_proj = rast('output_data/r_cos_aspect_proj.tif') 

r_medium_trees_constant_proj = rast('output_data/r_medium_trees_constant_proj.tif') 
r_medium_trees_midelevation_proj = rast('output_data/r_medium_trees_midelevation_proj.tif') 


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
  
  lambdas = mclapply(indices_predictors_non_na_ss, function(i) {
    print(i)
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
  }, mc.cores = NUM_CORES)
  
  df_predictors$lambda[indices_predictors_non_na_ss] = unlist(lambdas)
  
  r_lambda = r_predictors[[1]]
  r_lambda[] = df_predictors$lambda
  
  return(r_lambda)
}

r_lambda_n_medium_midelevation = make_ipm_spatial_scenario(r_medium_trees_scenario=r_medium_trees_midelevation_proj)  

writeRaster(r_lambda_n_medium_midelevation, 
            file='output_data/r_lambda_n_medium_midelevation.tif',
            overwrite=TRUE)

# currently this has no bootstrapping to get the uncertainty in dataset and/or sex propagated. 
# i think it is too computationally intensive to do this.