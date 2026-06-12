library(terra)
library(dplyr)
library(sf)


AGGREGATION_LEVEL = 40 # starting at 5 m, so this is multiplicative of that amount


# get topography at 5m resolution
r_elevation = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/altitudes.tif')
r_cos_aspect = cos(terrain(r_elevation, v='aspect',unit='radians'))

r_cos_aspect_aggregated = aggregate(r_cos_aspect, AGGREGATION_LEVEL)
r_elevation_aggregated = aggregate(r_elevation, AGGREGATION_LEVEL)



# get genetic predictors
r_sex = rast('data/rasters/CRBU_2018_mosaic_sex_masked.tif')
r_sex_projected = project(r_sex, r_elevation_aggregated, method='mode', threads=TRUE)
levels(r_sex_projected) = data.frame(id=c(0,1),sex=c('male','female')) # 1=female


r_ploidy_level = rast('data/rasters/CRBU_2018_mosaic_ploidy_masked.tif')
r_ploidy_level_projected = project(r_ploidy_level, r_elevation_aggregated, method='mode', threads=TRUE)
levels(r_ploidy_level_projected) = data.frame(id=c(0,1),cytotype=c('triploid','diploid'))

# make mask
r_mask = project(!is.na(r_sex), r_elevation_aggregated, method='average') >= 0.5


# mask genetic predictors
r_sex_projected_masked = mask(r_sex_projected, r_mask, maskvalue=0, updatevalue=NA)
r_ploidy_level_projected_masked = mask(r_ploidy_level_projected, r_mask, maskvalue=0, updatevalue=NA)


# count nnz pixels
na.omit(r_sex_projected_masked[]) %>% length

# make sure everything lines up
c(r_sex_projected_masked, r_ploidy_level_projected_masked, r_cos_aspect_aggregated, r_elevation_aggregated)


# write out rasters
writeRaster(r_sex_projected_masked, file='output_data/r_sex_projected_masked.tif', overwrite=TRUE)
writeRaster(r_ploidy_level_projected_masked, file='output_data/r_ploidy_level_projected_masked.tif', overwrite=TRUE)

writeRaster(r_cos_aspect_aggregated, file='output_data/r_cos_aspect_aggregated.tif', overwrite=TRUE)
writeRaster(r_elevation_aggregated, file='output_data/r_elevation_aggregated.tif', overwrite=TRUE)





