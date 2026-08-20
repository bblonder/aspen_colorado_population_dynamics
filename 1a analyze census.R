library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(car)
library(DHARMa)
library(ggbiplot)
library(glmmTMB)
library(ggeffects)
library(MuMIn)
library(performance)
library(sjPlot)
library(elevatr)
library(sf)
library(terra)
library(stars)
library(ggnewscale)
library(ggspatial)
library(rnaturalearth)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(RColorBrewer)

if (!file.exists('output_figures'))
{
  dir.create('output_figures')
}
if (!file.exists('output_data'))
{
  dir.create('output_data')
}

if(!("PREFIX_TYPE" %in% ls()))
{
  PREFIX_TYPE = 'v2closest11'
}

data_climate = read.csv('output_data/df_lagged_climate.csv')

data_site = read.csv(sprintf('data/aspen_data_site-level_2018-2023_%s_2024-11-27.csv',PREFIX_TYPE))
data_sex = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_aug_11_2021.csv') %>%
  mutate(site_code = Site_Code) %>%
  dplyr::select(-Site_Code,-X.UTM,-Y.UTM) %>%
  mutate(geneticSexID=ifelse(is.na(geneticSexID),'unknown',geneticSexID))
data_site_info = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv') %>%
  dplyr::select(site_code=Site_Code,Cos.aspect,Elevation,Slope,X.UTM,Y.UTM,Watershed)
data_ploidy_new = read.csv('/Users/benjaminblonder/Documents/berkeley/roxy cruz postdoc/aspen greenhouse/data analysis - best version with no exclusion bugs/00 update ploidy rmbl/data_rmbl_ploidy_updated_2025-12-04.csv') %>%
  dplyr::rename(site_code=Site_Code)

data_site = data_site %>%
  left_join(data_sex, by='site_code') %>%
  left_join(data_ploidy_new, by='site_code') %>%
  left_join(data_site_info, by='site_code') %>%
  mutate(Point_Type = ifelse(nchar(site_code)==4,'Grid','Random')) %>%
  mutate(Point_Type_detailed = ifelse(nchar(site_code)>5,'Random-opportunistic','Other')) %>%
  # give more illustrative names
  mutate(geneticSexID=as.character(factor(geneticSexID, levels=c('M','F'), labels = c('male','female')))) %>%
  mutate(geneticSexID=ifelse(is.na(geneticSexID),'unknown',geneticSexID)) %>%
  left_join(data_climate, by=c('site_code','year'))
  

# distribution of ploidy level and sex
data_site %>% select(site_code, Ploidy_level, Point_Type) %>% unique %>%
  group_by(Point_Type, Ploidy_level) %>%
  tally

data_site %>% select(site_code, geneticSexID, Point_Type) %>% unique %>%
  group_by(Point_Type, geneticSexID) %>%
  tally

# get recruitment quantiles
data_site %>% group_by(Point_Type) %>% summarize(quantile(n_medium_trees,0.95,na.rm=TRUE))

# how many plots total
data_site %>% filter(Point_Type=='Random') %>% pull(site_code) %>% unique %>% length

# how many plots per year
data_site %>% select(site_code, year, Point_Type_detailed) %>% unique %>%
  group_by(Point_Type_detailed, year) %>%
  tally

# any biases in un-called
data_site %>% select(site_code, geneticSexID, Ploidy_level) %>% unique %>%
  group_by(geneticSexID, Ploidy_level) %>%
  tally


# make map
# map out where site are missed
g_counts = ggplot(data_site %>% 
                    select(site_code, X.UTM, Y.UTM, Point_Type, year) %>% 
                    filter(Point_Type=='Random') %>% 
                    group_by(site_code, X.UTM, Y.UTM) %>% 
                    tally(), aes(x=X.UTM,y=Y.UTM,color=factor(n))) +
  geom_point() +
  theme_bw() +
  scale_color_brewer(palette='Set1')



coords = data_site %>% st_as_sf(coords=c("X.UTM","Y.UTM"), crs = 32613) %>%
  filter(Point_Type=='Random')
elev = get_elev_raster(st_bbox(st_buffer(coords, dist=1000)), z=12) %>%
  trim(-50)
sl <- terrain(rast(elev), "slope", unit = "radians")
asp <- terrain(rast(elev), "aspect", unit = "radians")
hill_single <- shade(sl, asp,
                     angle = 45,
                     direction = 300,
                     normalize = TRUE
)


g_map_inset = ggplot() + 
  theme_bw() + 
  xlab('Longitude') + ylab('Latitude') +
  geom_stars(data=st_as_stars(hill_single),show.legend=FALSE,alpha=0.5) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_stars(data=st_as_stars(elev),alpha=0.4) +
  scale_fill_viridis_c(option = "cividis",name='Elevation (m)',direction = -1) +
  annotation_scale() +
  annotation_north_arrow(location='br') +
  new_scale_fill() +
  new_scale_color() +
  geom_sf(data=coords, mapping=aes(color=Ploidy_level,shape=geneticSexID),size=1.5) + 
  scale_color_manual(values=c(alpha("blue",0.2),alpha("red",0.2),alpha("gray",0.2)),na.value = 'gray',name='Ploidy level') + 
  #scale_fill_manual(values=c(alpha("blue",0.5),alpha("red",0.5),alpha("gray",0.5)),na.value = 'gray',name='Ploidy level') + 
  scale_shape_discrete(name='Sex',na.value=3)



world <- ne_countries(scale = "medium", returnclass = "sf")
world_coordinates <- map_data("world")

g_map_larger = ggplot(data = world) +
  geom_sf(fill='#EEEEEE') +
  theme_bw() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  geom_sf(data=coords %>% st_bbox %>% st_as_sfc,color='purple',linewidth=1,size=4)  +
  coord_sf(xlim = c(-125,-70),ylim=c(20,52))



g_map_final = ggdraw(g_map_inset) +
  draw_plot(g_map_larger, x=0.0,y=0.75,width=0.25,height=0.25)
ggsave(g_map_final, file='output_figures/FIG1_map.png',width=6,height=6)
ggsave(g_map_final, file='output_figures/FIG1_map.pdf',width=6,height=6)


hill_single_trimmed = crop(hill_single,st_bbox(st_buffer(coords, dist=500))) 

# make SI maps of remotely sensed predictor layers
r_sex = rast('data/rasters/CRBU_2018_mosaic_sex_masked.tif')
r_sex_downsampled = crop(aggregate(r_sex,5,fun='median',na.rm=TRUE), st_bbox(st_buffer(coords, dist=500)))
levels(r_sex_downsampled) = data.frame(id=c(0,1),sex=c('male','female')) # 1=female

r_ploidy_level = rast('data/rasters/CRBU_2018_mosaic_ploidy_masked.tif')
r_ploidy_level_downsampled = crop(aggregate(r_ploidy_level,5,fun='median',na.rm=TRUE), st_bbox(st_buffer(coords, dist=500)))
levels(r_ploidy_level_downsampled) = data.frame(id=c(0,1),cytotype=c('triploid','diploid'))

g_sex = ggplot() + 
  theme_bw() + 
  geom_stars(data=st_as_stars(hill_single_trimmed),show.legend=FALSE,alpha=0.5) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_stars(data=st_as_stars(r_sex_downsampled),alpha=0.75) +
  scale_fill_manual(values=c('purple3','green4'),name='Sex',na.value=NA) +
  annotation_scale() +
  annotation_north_arrow(location='br') +
  coord_sf() + 
  xlab("Easting (m)") + ylab("Northing (m)")

g_cytotype = ggplot() + 
  theme_bw() + 
  geom_stars(data=st_as_stars(hill_single_trimmed),show.legend=FALSE,alpha=0.5) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_stars(data=st_as_stars(r_ploidy_level_downsampled),alpha=0.75) +
  scale_fill_manual(values=c('red','blue'),name='Cytotype',na.value=NA) +
  annotation_scale() +
  annotation_north_arrow(location='br') +
  coord_sf() + 
  xlab("Easting (m)") + ylab("Northing (m)")

ggsave(ggarrange(g_cytotype, g_sex, labels='AUTO',nrow=2,ncol=1),width=5,height=8, file='output_figures/g_sex_cytotype.png')



# basic survey counts
data_site %>% group_by(year, Point_Type) %>% tally

# count up the number of sites that we visited
data_site_by_year = data_site %>%
  group_by(site_code, year) %>%
  tally() %>%
  pivot_wider(id_cols=site_code, names_from=year,values_from=n, values_fill=0) %>%
  pivot_longer(!site_code) %>%
  left_join(data_site %>% select(site_code, Point_Type) %>% unique,by='site_code')

g_which_sampled = ggplot(data_site_by_year, aes(x=name,y=site_code,fill=factor(value))) + geom_tile() +
    facet_wrap(~Point_Type,scales='free') +
    scale_fill_manual(values=c('lightgray','blue'),name='Visited') +
  xlab('Year') +
  ylab('Plot name')
ggsave(g_which_sampled, file=sprintf('output_figures/g_which_sampled_%s.pdf',PREFIX_TYPE),width=5,height=27)
ggsave(g_which_sampled, file=sprintf('output_figures/g_which_sampled_%s.png',PREFIX_TYPE),width=5,height=27, dpi=600)


data_counts_by_year = data_site %>%
  group_by(year, Point_Type) %>%
  tally

g_total_counts = ggplot(data_counts_by_year, aes(x=factor(year),y=n)) +
  geom_bar(stat='identity') +
  facet_wrap(~Point_Type) +
  theme_bw()
ggsave(g_total_counts, file=sprintf('output_figures/g_total_counts_%s.pdf',PREFIX_TYPE),width=8,height=6)



t_aspen_counts = table(data_site$n_aspen_trees_surveyed, data_site$year)
t_aspen_counts

write.csv(t_aspen_counts, file=sprintf('output_data/t_aspen_counts_%s.csv', PREFIX_TYPE), row.names=TRUE)
# check on details - this all looks OK to me
write.csv(data_site %>% 
            filter(year==2023 & n_aspen_trees_surveyed > 1 & n_aspen_trees_surveyed < 10), 
          file = sprintf('output_data/t_aspen_counts_details_%s.csv', PREFIX_TYPE), row.names = FALSE)



fit_lm <- function(df, yvar)
{
  m_lm  = NULL
  try(m_lm <- lm(formula(sprintf("%s~year",yvar)),data=df))
          
  if (!is.null(m_lm))
  {
    coefs = data.frame(t(coef(m_lm)))
    names(coefs)=c('int','slope')
  }
  else
  {
    coefs = data.frame(int=NA,slope=NA)
  }
  return(coefs)
}

# do time series analyses
slopes_mortality = data_site %>% 
  filter(Point_Type=='Random') %>%
  group_by(site_code) %>%
  do(fit_lm(.,"frac_adult_dead_w_background_mortality_estimate"))

slopes_mortality %>%
  ungroup %>%
  summarize(slope_mean = mean(slope, na.rm=TRUE), slope_sd = sd(slope, na.rm=TRUE)) * 100

slopes_recruitment = data_site %>% 
  filter(Point_Type=='Random') %>%
  group_by(site_code) %>%
  do(fit_lm(.,"n_medium_trees"))

slopes_recruitment %>%
  ungroup %>%
  summarize(slope_mean = mean(slope, na.rm=TRUE), slope_sd = sd(slope, na.rm=TRUE))

slopes_size = data_site %>% 
  mutate(dbh_mean_live = ifelse(is.na(dbh_mean_live), dbh_center_live, dbh_mean_live)) %>%
  filter(Point_Type=='Random') %>%
  group_by(site_code) %>%
  do(fit_lm(.,"dbh_mean_live"))

slopes_size %>%
  ungroup %>%
  summarize(slope_mean = mean(slope, na.rm=TRUE), slope_sd = sd(slope, na.rm=TRUE))



g_mortality_by_site = ggplot(data_site %>% filter(Point_Type=='Random'), aes(x=year,y=frac_adult_dead_w_background_mortality_estimate,group=site_code)) +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.25,size=1) +
  theme_bw() +
  facet_grid(geneticSexID~Ploidy_level) +
  ylim(0,1) + 
  xlab('Year') + ylab("Fraction adult stems dead (corrected)") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
ggsave(g_mortality_by_site, file=sprintf('output_figures/g_site_mortality_%s.pdf', PREFIX_TYPE), width=6,height=6)
ggsave(g_mortality_by_site, file=sprintf('output_figures/g_site_mortality_%s.png', PREFIX_TYPE), width=6,height=6)


g_damage_by_site = ggplot(data_site %>% filter(Point_Type=='Random'), aes(x=year,y=frac_adult_damaged,group=site_code)) +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  theme_bw() +
  facet_grid(geneticSexID~Ploidy_level) +
  ylim(0,1) + 
  xlab('Year') + ylab("Fraction adult stems damaged") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
ggsave(g_damage_by_site, file=sprintf('output_figures/g_site_damage_%s.pdf', PREFIX_TYPE), width=6,height=6)
ggsave(g_damage_by_site, file=sprintf('output_figures/g_site_damage_%s.png', PREFIX_TYPE), width=6,height=6)


g_n_medium_by_site = ggplot(data_site %>% filter(Point_Type=='Random'), aes(x=year,y=n_medium_trees/(pi*3^2),group=site_code)) +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  theme_bw() +
  facet_grid(geneticSexID~Ploidy_level) +
  #scale_y_sqrt() +
  xlab('Year') + ylab(expression(paste('n'[S], ' (m'^{-2},')'))) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
ggsave(g_n_medium_by_site, file=sprintf('output_figures/g_site_medium_%s.pdf', PREFIX_TYPE), width=6,height=6)
ggsave(g_n_medium_by_site, file=sprintf('output_figures/g_site_medium_%s.png', PREFIX_TYPE), width=6,height=6)


g_n_small_by_site = ggplot(data_site %>% filter(Point_Type=='Random'), aes(x=year,y=n_small_trees/(pi*3^2),group=site_code)) +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  theme_bw() +
  facet_grid(geneticSexID~Ploidy_level) +
  #scale_y_sqrt() +
  xlab('Year') + ylab(expression(paste('Stem density with dgh smaller than pencil (m'^{-2},')'))) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
ggsave(g_n_small_by_site, file=sprintf('output_figures/g_site_small_%s.pdf', PREFIX_TYPE), width=6,height=6)
ggsave(g_n_small_by_site, file=sprintf('output_figures/g_site_small_%s.png', PREFIX_TYPE), width=6,height=6)



#


# deal with growth rates

data_rgr = data_site %>% 
  select(site_code, year, dbh_center_live, Ploidy_level) %>% 
  pivot_wider(id_cols=site_code, names_from=year, values_from = dbh_center_live) %>% 
  mutate(rgr_2020_2018 = (`2020`-`2018`)/`2018`/2) %>% 
  mutate(rgr_2023_2020 = (`2023`-`2020`)/`2020`/3) %>% 
  mutate(gr_2020_2018 = (`2020`-`2018`)/2) %>% 
  mutate(gr_2023_2020 = (`2023`-`2020`)/3) %>%
  mutate(gr_2023_2018 = (`2023`-`2018`)/5) %>%
  mutate(outlier_growth = abs(gr_2020_2018) > 1 | abs(gr_2023_2020) > 1 | abs(gr_2023_2018) > 1)

# ggplot(data_rgr, aes(x=gr_2020_2018,y=gr_2023_2020)) +
#   geom_point()

data_site_with_growth = data_site %>%
  left_join(data_rgr,by='site_code') %>%
  mutate(gr_average = (gr_2023_2018 + gr_2023_2020)/2)

# look at stem growth
g_growth_rate_by_site = ggplot(data_site_with_growth, 
                               aes(x=year, y=dbh_center_live, group=site_code, color=outlier_growth)) + 

  geom_line(alpha=0.5) + 
  geom_point(alpha=0.5) +
  facet_grid(geneticSexID~Ploidy_level) +
  xlab('Year') + ylab("x (focal stem, cm)") + 
  theme_bw() +
  scale_color_manual(values=c('black','red')) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  theme(legend.position='none')
ggsave(g_growth_rate_by_site, file=sprintf('output_figures/g_growth_rate_by_site_%s.pdf', PREFIX_TYPE),width=6,height=6)
ggsave(g_growth_rate_by_site, file=sprintf('output_figures/g_growth_rate_by_site_%s.png', PREFIX_TYPE),width=6,height=6)


# find outliers
outliers_growth = data_rgr %>%
  filter(outlier_growth==TRUE) %>%
  select(site_code, `2018`,`2020`,`2023`, starts_with('rgr'), starts_with('gr')) 
outliers_growth %>%
  write.csv(sprintf('output_data/t_growth_rate_high_%s.csv',PREFIX_TYPE))



g_basal_area = ggplot(data_site %>% 
                        filter(Point_Type=='Random') %>%
                        mutate(outlier_growth = site_code %in% outliers_growth$site_code), 
                               aes(x=year,y=basal_area_density_live,group=site_code,color=outlier_growth)) +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  theme_bw() +
  facet_grid(geneticSexID~Ploidy_level) +
  scale_color_manual(values=c('black','red')) + 
  #scale_y_sqrt() +
  xlab('Year') + ylab(expression(paste('BAD', ' (m'^2, ' m'^{-2},')'))) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  theme(legend.position='none')
ggsave(g_basal_area, file=sprintf('output_figures/g_basal_area_by_site_%s.pdf', PREFIX_TYPE),width=6,height=6)
ggsave(g_basal_area, file=sprintf('output_figures/g_basal_area_by_site_%s.png', PREFIX_TYPE),width=6,height=6)



### show where sites are

# map out where site are missed
g_counts = ggplot() +
  geom_stars(data=st_as_stars(hill_single_trimmed),show.legend=FALSE,alpha=0.5) +
  scale_fill_distiller(palette = "Greys") +
  geom_point(data=data_site %>% 
               select(site_code, X.UTM, Y.UTM, Point_Type, year) %>% 
               group_by(site_code, X.UTM, Y.UTM, Point_Type) %>% 
               tally(), aes(x=X.UTM,y=Y.UTM,color=factor(n)),size=0.5) +
  theme_bw() +
  facet_wrap(~Point_Type) +
  scale_color_viridis_d(option='viridis',name='# re-visits') +
  coord_equal() +
  xlab('Easting (m)') + ylab('Northing (m)') +
  theme(legend.position='bottom')
ggsave(g_counts, file=sprintf('output_figures/g_counts_%s.pdf',PREFIX_TYPE),width=8,height=5)
ggsave(g_counts, file=sprintf('output_figures/g_counts_%s.png',PREFIX_TYPE),width=8,height=5)


g_counts_grid = ggplot() +
  #geom_stars(data=st_as_stars(hill_single_trimmed),show.legend=FALSE,alpha=0.5) +
  #scale_fill_distiller(palette = "Greys") +
  geom_point(data=data_site %>% 
               select(site_code, X.UTM, Y.UTM, Point_Type, Watershed, year) %>% 
               filter(Point_Type=='Grid') %>% 
               group_by(site_code, Watershed, X.UTM, Y.UTM) %>% 
               tally(), mapping=aes(x=X.UTM,y=Y.UTM,color=factor(n))) +
  theme_bw() +
  scale_color_viridis_d(option='viridis',name='# re-visits') +
  xlab('Easting (m)') + ylab('Northing (m)') +
  theme(legend.position='bottom') +
  facet_wrap(~Watershed,scales='free')
ggsave(g_counts_grid, file=sprintf('output_figures/g_counts_grid_%s.pdf',PREFIX_TYPE),width=8,height=8)
ggsave(g_counts_grid, file=sprintf('output_figures/g_counts_grid_%s.png',PREFIX_TYPE),width=8,height=8)



# make counts of plots by type
data_counts_sex_ploidy = data_site %>%
  select(site_code, Ploidy_level, geneticSexID, Point_Type) %>%
  ungroup %>%
  group_by(Ploidy_level, geneticSexID, Point_Type) %>%
  tally

g_counts_sex_ploidy = ggplot(data_counts_sex_ploidy, aes(fill=Ploidy_level,x=geneticSexID,y=n)) +
  geom_bar(stat='identity',position='dodge') +
  theme_bw() +
  ylab("Number of sites") + 
  scale_fill_manual(values=c('blue','red','gray'),na.value = 'gray',name='Cytotype') +
  facet_wrap(~Point_Type) +
  xlab("Sex")
ggsave(g_counts_sex_ploidy, file=sprintf('output_figures/g_counts_sex_ploidy_%s.pdf', PREFIX_TYPE),width=8,height=5)
ggsave(g_counts_sex_ploidy, file=sprintf('output_figures/g_counts_sex_ploidy_%s.png', PREFIX_TYPE),width=8,height=5)




data_site_with_growth_for_pca = data_site_with_growth %>% 
  mutate(dbh_mean_live = ifelse(is.na(dbh_mean_live), dbh_center_live, dbh_mean_live)) %>%
  # set growth rates for all 3 years
  mutate(growth_rate=ifelse(year==2018, gr_2020_2018,
                            ifelse(year==2020, 0.5*(gr_2023_2020+gr_2020_2018),
                                   ifelse(year==2023, gr_2023_2020, NA)))) %>%
  # use averages if needed
  mutate(growth_rate = ifelse(is.na(growth_rate), gr_2023_2018, growth_rate)) %>%
  # remove outliers
  mutate(growth_rate = ifelse(abs(growth_rate)>1, NA, growth_rate))







plot_pca <- function(pca_this, data_this, name)
{
  variance_fraction = 100*(pca_this$sdev^2)/sum(pca_this$sdev^2)
  
  ch_indices = chull(pca_this$x[,1:2])
  
  g_pca_scores = ggplot(pca_this$x %>% 
           cbind(data_this %>% select(site_code, year, Ploidy_level, geneticSexID)),
         aes(x=PC1,y=PC2,color=factor(year))) +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    geom_point() +
    facet_grid(geneticSexID~Ploidy_level) +
    theme_bw() +
    coord_equal() +
    scale_color_viridis_d(name='Year') +
    geom_polygon(data=pca_this$x[ch_indices,], inherit.aes = FALSE, mapping=aes(x=PC1,y=PC2),fill=NA,color='gray20') +
    theme(legend.position='bottom') +
    xlab(sprintf("PC1 (%.1f%%)",variance_fraction[1])) + ylab(sprintf("PC2 (%.1f%%)",variance_fraction[2]))
  
    # scale_color_manual(values=c('blue','red'),na.value = 'gray',name='Cytotype') + 
    # scale_shape_discrete(name='Sex',na.value=3)
  
  g_pca_loadings = ggplot(pca_this$rotation %>% 
                            as.data.frame %>% 
                            mutate(var=sapply(sapply(gsub("_"," ", row.names(.)),strwrap,width=20),paste,collapse='\n'))) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),arrow = arrow(length = unit(0.1, "inches")),color='gray') +
    geom_label_repel(aes(x=PC1*1.1, y=PC2*1.1, label=var),size=3,alpha=0.5) +
    theme_bw() +
    coord_equal() +
    xlab(sprintf("PC1 (%.1f%%)",variance_fraction[1])) + ylab(sprintf("PC2 (%.1f%%)",variance_fraction[2])) +
    xlim(-1,1)+ylim(-1,1)
  
  ggsave(ggarrange(g_pca_loadings, g_pca_scores, nrow=2,ncol=1,labels='AUTO'),
         file=sprintf('output_figures/g_pca_plot_%s_%s.pdf', name, PREFIX_TYPE),
         width=6,height=9)
  ggsave(ggarrange(g_pca_loadings, g_pca_scores, nrow=2,ncol=1,labels='AUTO'),
         file=sprintf('output_figures/g_pca_plot_%s_%s.png', name, PREFIX_TYPE),
         width=6,height=9)
}


# PCA of demographic variables by year (all variables)
data_pca_all = data_site_with_growth_for_pca %>% 
  filter(Point_Type=='Random') %>%
  mutate(n_S = n_medium_trees / (pi*3^2)) %>%
  select(site_code, year, Ploidy_level, geneticSexID, 
         #Cos.aspect, Elevation, Slope,
         n_pencil = n_small_trees, 
         nS=n_S,
         dbh_mean_live,
         bad_live = basal_area_density_live,
         frac_adult_dead_w_background_mortality_estimate,
         frac_adult_damaged,
         growth_rate#,
         # n_small_dead,
         # n_medium_dead,
         # n_dead_down
  )

rows_na_pca_full = data_pca_all %>%
  select(n_pencil:growth_rate) %>%
  rowSums %>%
  is.na

data_pca_all = data_pca_all[which(!rows_na_pca_full),]

pca_all = prcomp(data_pca_all %>%
                   select(n_pencil:growth_rate),
                 center=TRUE, scale=TRUE)

plot_pca(pca_all, data_pca_all, 'all_vars')







with(pca_all, 100*(sdev^2)/sum(sdev^2))





















######################## 
# demography
df_for_plot_level_model = data_site_with_growth %>%
  filter(Point_Type=='Random') %>%
  mutate(Ploidy_level=ifelse(Ploidy_level=="unknown",NA,Ploidy_level),Ploidy_level) %>%
  mutate(geneticSexID=ifelse(geneticSexID=="unknown",NA,geneticSexID),geneticSexID) %>%
  mutate(Ploidy_level = factor(Ploidy_level)) %>%
  mutate(geneticSexID = factor(geneticSexID)) %>%
  # take center tree DBH as mean if only one tree censused
  mutate(dbh_mean_live = ifelse(is.na(dbh_mean_live), dbh_center_live, dbh_mean_live)) %>%
  mutate(years_after_2018 = year - 2018)
  # filter(!is.na(Ploidy_level)) %>%
  # filter(!is.na(geneticSexID))
  # mutate(Ploidy_level = factor(Ploidy_level)) %>%
  # mutate(geneticSexID = factor(geneticSexID))
  



fit_plot_level_model <- function(response_var, family, data, zi=TRUE, title=response_var)
{
  formula_this = sprintf("%s ~ years_after_2018 + Ploidy_level * geneticSexID + (1|site_code)", response_var)
  
  m_this = glmmTMB(formula=formula(formula_this), 
                   data=data,
                   ziformula = formula(ifelse(zi==TRUE, "~ Ploidy_level * geneticSexID + years_after_2018", "~0")),
                   family=family,
                   control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
  
  dharma_this = simulateResiduals(m_this)

  plot_resid_this = ggplot(data.frame(resid=resid(m_this)),aes(x=resid)) +
    geom_density() +
    geom_vline(xintercept = 0,color='red') +
    theme_bw()
  
  r2_this = c(r2(m_this), r2_zeroinflated(m_this)) 
  
  return(list(model = m_this,
              summary = summary(m_this),
              dharma = dharma_this,
              #preds_this = preds_this,
              #plot_effect = plot_effect_this,
              plot_resid = plot_resid_this,
              r2 = r2_this))
}

results_frac_dead = fit_plot_level_model(response_var="frac_adult_dead_w_background_mortality_estimate", 
                              data=df_for_plot_level_model, family=lognormal, zi=TRUE, 
                              title = 'Fraction adult dead')

results_size = fit_plot_level_model(response_var="dbh_mean_live", 
                                         data=df_for_plot_level_model, family=lognormal, zi=FALSE, 
                                         title = 'Focal tree DBH (cm)')

results_n_medium = fit_plot_level_model(response_var="I(n_medium_trees/(pi*3^2))", 
                                    data=df_for_plot_level_model, family=lognormal, zi=TRUE, 
                                    title = 'Number medium saplings')

# 
# results_size$summary
# results_n_medium$summary
# results_frac_dead$summary





tab_model(results_size$model, file='output_figures/table_model_plot_level_size.html',encoding='UTF-16')

tab_model(results_frac_dead$model, file='output_figures/table_model_plot_level_mortality.html',encoding='UTF-16')

tab_model(results_n_medium$model, file='output_figures/table_model_plot_level_recruitment.html',encoding='UTF-16')



do_ggpredict <- function(model, terms, zi)
{
  preds_this = ggpredict(model=model,
                         terms=terms,
                         type=ifelse(zi==TRUE,'zero_inflated','fixed'),
                         condition=c(Elevation=3000,
                                     Cos.aspect=-1))
  
  plot_effect_this = plot(preds_this,
                          show_data=TRUE,
                          jitter=c(0.1,0))
  
    
  
  return(plot_effect_this)
}

show_ggpredict <- function(model, zi)
{
  g1 = do_ggpredict(model = model, 
               zi=zi, 
               terms=c('years_after_2018','Ploidy_level','geneticSexID')) +
    theme_bw() +
    scale_color_manual(values=c('blue','red','gray')) +
    scale_fill_manual(values=c('blue','red','gray')) +
    labs(color='Cytotype') + 
    labs(title = NULL) +
    xlab('Years after 2018')
  
  g1
  
  # g2 = do_ggpredict(model = model, 
  #              zi=zi, 
  #              terms=c('SWE.1','STB.1')) +
  #   theme_bw() +
  #   scale_color_gradient() +
  #   scale_fill_gradient() + 
  #   labs(title = NULL)
  # 
  # ggarrange(g1, g2)
}

gp1 = show_ggpredict(results_size$model, zi=FALSE) + ylab('DBH mean (cm)')
gp2 = show_ggpredict(results_frac_dead$model, zi=FALSE) + ylab('Fraction adult dead')
gp3 = show_ggpredict(results_n_medium$model, zi=FALSE) + ylab('Number medium stems') + scale_y_sqrt()


g_plot_level = ggarrange(gp1, gp2, gp3,
                          nrow=2,ncol=2,labels='AUTO',
                          common.legend = TRUE,
                          legend = 'bottom')

ggsave(g_plot_level, file=sprintf('output_figures/g_plot_level_temporal_trends_%s.pdf', PREFIX_TYPE),width=7,height=7)
ggsave(g_plot_level, file=sprintf('output_figures/g_plot_level_temporal_trends_%s.png', PREFIX_TYPE),width=7,height=7)







# make tables of distributions

g_dist1 = ggplot(data_site, aes(x=n_medium_trees/(pi*3^2),color=factor(year),fill=factor(year))) + 
  geom_density(alpha=0.2) +
  theme_bw() +
  facet_wrap(~Point_Type) +
  xlab(expression(paste('n'[S], ' (m'^{-2},')'))) +
  scale_color_discrete(name='Year') +
    scale_fill_discrete(name='Year') +
  ylab('Probability density')

g_dist2 = ggplot(data_site, aes(x=n_small_trees/(pi*3^2),color=factor(year),fill=factor(year))) + 
  geom_density(alpha=0.2) +
  theme_bw() +
  facet_wrap(~Point_Type) +
  xlab(expression(paste('n'['tiny']))) +
  scale_color_discrete(name='Year') +
      scale_fill_discrete(name='Year') +
  ylab('Probability density')

g_dist3 = ggplot(data_site, aes(x=dbh_center_live,color=factor(year),fill=factor(year))) + 
  geom_density(alpha=0.2) +
  theme_bw() +
  facet_wrap(~Point_Type) +
  xlab('x (focal stem, cm)') +
  scale_color_discrete(name='Year') +
      scale_fill_discrete(name='Year') +
  ylab('Probability density')

g_dist = ggarrange(g_dist3, g_dist1,
          nrow=1,ncol=2,
          common.legend = TRUE,
          legend='bottom',
          labels='AUTO',
          align='hv')
ggsave(g_dist, file=sprintf('output_figures/g_distribution_%s.pdf', PREFIX_TYPE),width=7,height=4)
ggsave(g_dist, file=sprintf('output_figures/g_distribution_%s.png', PREFIX_TYPE),width=7,height=4)



# make plot of climate data
r_sex = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_prediction_20250110.tif')
fn_stb = dir(path='data/rasters/',pattern='GRIDMET_DROUGHT_short_term_blend',full.names = TRUE)
r_stb = rast(fn_stb)
r_stb_projected = project(r_stb, aggregate(r_sex, 20))
names(r_stb_projected) = paste(2016:2023)

fn_swe = dir(path='data/rasters/',pattern='SNODAS_SWE',full.names = TRUE)
r_swe = rast(fn_swe)
r_swe_projected = project(r_swe, aggregate(r_sex, 20))
names(r_swe_projected) = paste(2016:2023)

fn_tmax = dir(path='data/rasters/',pattern='_tmax_',full.names = TRUE)
r_tmax = rast(fn_tmax)
r_tmax_projected = project(r_tmax, aggregate(r_sex, 20))
names(r_tmax_projected) = paste(2016:2023)


plots_stb = lapply(r_stb_projected - mean(r_stb_projected), function(r) {
  ggplot() + 
    geom_stars(data=st_as_stars(r)) +
    scale_fill_gradient2(name='Detrended STB (August)',limits=c(-2,2)) +
    annotation_scale() +
    theme_void() + 
    ggtitle(names(r)) +
    coord_equal()
})

plots_swe = lapply(r_swe_projected - mean(r_swe_projected), function(r) {
  ggplot() + 
    geom_stars(data=st_as_stars(r)) +
    scale_fill_continuous(palette=brewer.pal(9, 'PiYG'), name='Detrended SWE (May-June) (mm)',limits=c(-600,600)) +
    annotation_scale() +
    theme_void() + 
    ggtitle(names(r)) +
    coord_equal()
})

plots_tmax = lapply(r_tmax_projected - mean(r_tmax_projected), function(r) {
  ggplot() + 
    geom_stars(data=st_as_stars(r)) +
    scale_fill_continuous(palette=brewer.pal(9, 'RdBu'), name='Detrended mean Tmax (June - August) (°C)',limits=c(-2.5,2.5)) +
    annotation_scale() +
    theme_void() + 
    ggtitle(names(r)) +
    coord_equal()
})

#g_plot_stb = ggarrange(plotlist=plots_stb,nrow=1,common.legend = TRUE,legend='bottom')
g_plot_swe = ggarrange(plotlist=plots_swe,nrow=1,common.legend = TRUE,legend='bottom')
g_plot_tmax = ggarrange(plotlist=plots_tmax,nrow=1,common.legend = TRUE,legend='bottom')

#g_plot_stb, 
ggsave(ggarrange(g_plot_swe, g_plot_tmax, nrow=2,labels='AUTO'), file=sprintf('output_figures/g_climate_rasters_%s.pdf', PREFIX_TYPE),width=7,height=4.5)
ggsave(ggarrange(g_plot_swe, g_plot_tmax, nrow=2,labels='AUTO'), file=sprintf('output_figures/g_climate_rasters_%s.png', PREFIX_TYPE),width=7,height=4.5)
