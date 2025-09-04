library(terra)
library(dplyr)
library(ggplot2)
library(RStoolbox)
library(GGally)
library(ggplot2)
library(stars)
library(ggspatial)
library(ggpubr)

r_scenario_constant = rast('output_data/r_lambda_n_medium_constant_5.tif')
r_scenario_midelevation = rast('output_data/r_lambda_n_medium_midelevation.tif')

make_map <- function(raster, title)
{
  g = ggplot() + 
    theme_bw() + 
    geom_stars(data=st_as_stars(raster)) +
    scale_fill_gradient2(midpoint=1,name=expression(paste(lambda)),
                         low = 'darkorange',high='darkorchid1',mid = 'gray',
                         limits=c(0.9,1.1)) +
    annotation_scale() +
    annotation_north_arrow(location='br') +
    coord_equal() +
    xlab('Easting (m)') + ylab('Northing (m)') +
    ggtitle(title)
  
  return(g)
}

g1 = make_map(r_scenario_constant, 'Recruitment spatially constant')
g2 = make_map(r_scenario_midelevation, 'Recruitment peak at mid-elevation')

ggsave(ggarrange(g1, g2,labels='AUTO', common.legend = TRUE, legend='bottom',nrow=1,ncol=2), 
                   file='output_figures/g_lambda_map.png',width=8,height=4)



