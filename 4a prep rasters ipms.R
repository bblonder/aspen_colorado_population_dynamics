library(terra)
library(dplyr)

AGGREGATION_LEVEL = 100 # starting at 5 m, so this is multiplicative of that amount

r_sex = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_prediction_20250110.tif')
r_ploidy_level = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/min_phase_cytotype_medfilt-seived.tif')
r_ploidy_level[r_ploidy_level==-1] = NA # all 0 and 1
r_cos_aspect = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/r_cos_aspect.tif')
# r_slope = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/r_slope.tif')
r_elevation = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/altitudes.tif')

r_elevation_coarse = aggregate(r_elevation, AGGREGATION_LEVEL, fun='mean')

# project to elevation (coarser resolution) for now
r_sex_proj = project(!(r_sex-1), r_elevation_coarse, method='mode', threads=TRUE) # set male=1, female=0
r_ploidy_level_proj = project(!r_ploidy_level, r_elevation_coarse, method='mode', threads=TRUE) # set triploid=1, diploid=0
r_cos_aspect_proj = project(r_cos_aspect, r_elevation_coarse)


writeRaster(r_elevation_coarse, file='output_data/r_elev_proj.tif', overwrite=TRUE)
writeRaster(r_cos_aspect_proj, file='output_data/r_cos_aspect_proj.tif', overwrite=TRUE)

writeRaster(r_sex_proj, file='output_data/r_sex_proj.tif', overwrite=TRUE)
writeRaster(r_ploidy_level_proj, file='output_data/r_ploidy_level_proj.tif', overwrite=TRUE)





