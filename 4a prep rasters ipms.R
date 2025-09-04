library(terra)
library(dplyr)
library(ggplot2)
library(RStoolbox)
library(GGally)
library(ggplot2)
library(pbapply)
library(ipmr)

load('output_data/workspace for ipm.Rdata')

AGGREGATION_LEVEL = 10 # starting at 5 m, so this is multiplicative of that amount

r_sex = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_prediction_20250110.tif')
r_ploidy_level = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/min_phase_cytotype_medfilt-seived.tif')
r_ploidy_level[r_ploidy_level==-1] = NA # all 0 and 1
r_cos_aspect = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/r_cos_aspect.tif')
# r_slope = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/r_slope.tif')
r_elevation = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/altitudes.tif')

warning('need to set correct aggregation amount')
r_elevation_coarse = aggregate(r_elevation, AGGREGATION_LEVEL, fun='mean')

# project to elevation (coarser resolution) for now
warning('need to make sure these rasters have levels coded correctly')
r_sex_proj = project(!(r_sex-1), r_elevation_coarse, method='mode', threads=TRUE) # set male=1, female=0
r_ploidy_level_proj = project(!r_ploidy_level, r_elevation_coarse, method='mode', threads=TRUE) # set triploid=1, diploid=0
r_cos_aspect_proj = project(r_cos_aspect, r_elevation_coarse)
#r_slope_proj = project(r_slope, r_elevation)

# 

# r_ploidy_level_proj_scenario = r_ploidy_level_proj
# m_pl = glm(Ploidy_level~Elevation,data=df_sites_for_ipm_joined, family=binomial)
# r_ploidy_level_proj_scenario[] = predict(m_pl, newdata=data.frame(Elevation=as.numeric(r_elevation[])))

r_medium_trees_constant_proj = r_cos_aspect_proj
r_medium_trees_constant_proj[] = 5 #median(df_sites_for_ipm_joined$n_medium_trees)

r_medium_trees_midelevation_proj = r_cos_aspect_proj
m_nmt = nls(n_medium_trees ~ k*exp(-1/2*(Elevation-mu)^2/sigma^2),start=c(mu=3000,sigma=500,k=3) ,
            data = transitions_all_filtered_joined_no_na)
# do an optimistic scenario
r_medium_trees_midelevation_proj[] = 4*predict(m_nmt, newdata=data.frame(Elevation=as.numeric(r_elevation_coarse[])))


writeRaster(r_sex_proj, file='output_data/r_sex_proj.tif', overwrite=TRUE)
writeRaster(r_ploidy_level_proj, file='output_data/r_ploidy_level_proj.tif', overwrite=TRUE)
writeRaster(r_cos_aspect_proj, file='output_data/r_cos_aspect_proj.tif', overwrite=TRUE)

writeRaster(r_medium_trees_constant_proj, file='output_data/r_medium_trees_constant_proj.tif', overwrite=TRUE)
writeRaster(r_medium_trees_midelevation_proj, file='output_data/r_medium_trees_midelevation_proj.tif', overwrite=TRUE)










