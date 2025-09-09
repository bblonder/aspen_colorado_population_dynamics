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
  PREFIX_TYPE = 'v1allAspen'
}

data_site = read.csv(sprintf('data/aspen_data_site-level_2018-2023_%s_2024-11-27.csv',PREFIX_TYPE))
data_sex = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_aug_11_2021.csv') %>%
  mutate(site_code = Site_Code) %>%
  dplyr::select(-Site_Code,-X.UTM,-Y.UTM)
data_ploidy = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv') %>%
  dplyr::select(site_code=Site_Code,Ploidy_level,Cos.aspect,Elevation,Slope,X.UTM,Y.UTM,Watershed)
data_site = data_site %>%
  left_join(data_sex, by='site_code') %>%
  left_join(data_ploidy, by='site_code') %>%
  mutate(Point_Type = ifelse(nchar(site_code)==4,'Grid','Random'))

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

library(elevatr)
library(sf)
library(terra)
library(stars)
library(ggnewscale)
library(ggspatial)

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
  geom_sf(data=coords, mapping=aes(color=Ploidy_level,shape=geneticSexID),size=0.5) + 
  scale_color_manual(values=c('blue','red'),na.value = 'gray',name='Cytotype') + 
  scale_shape_discrete(name='Sex',na.value=3)

library("rnaturalearth")
library("ggpubr")

world <- ne_countries(scale = "medium", returnclass = "sf")
world_coordinates <- map_data("world")

g_map_larger = ggplot(data = world) +
  geom_sf(fill='#EEEEEE') +
  theme_bw() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  geom_sf(data=coords %>% st_bbox %>% st_as_sfc,color='black',linewidth=1)  +
  coord_sf(xlim = c(-125,-70),ylim=c(20,52))


library(cowplot)

g_map_final = ggdraw(g_map_inset) +
  draw_plot(g_map_larger, x=0.0,y=0.75,width=0.25,height=0.25)
ggsave(g_map_final, file='output_figures/FIG1_map.png',width=6,height=6)


hill_single_trimmed = crop(hill_single,st_bbox(st_buffer(coords, dist=500))) 

# make SI maps of remotely sensed predictor layers
r_sex = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_prediction_20250110.tif')
r_sex_downsampled = crop(aggregate(r_sex,5,fun='median',na.rm=TRUE), st_bbox(st_buffer(coords, dist=500))) - 1 # shift from 1 2 to 0 1 classification
levels(r_sex_downsampled) = c('male','female')

r_ploidy_level = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/min_phase_cytotype_medfilt-seived.tif')
crs(r_ploidy_level) = crs(hill_single)
r_ploidy_level[r_ploidy_level==-1] = NA # all 0 and 1
r_ploidy_level_downsampled = crop(aggregate(r_ploidy_level,5,fun='median',na.rm=TRUE), st_bbox(st_buffer(coords, dist=500)))
levels(r_ploidy_level_downsampled) = c('diploid','triploid')

g_sex = ggplot() + 
  theme_bw() + 
  geom_stars(data=st_as_stars(hill_single_trimmed),show.legend=FALSE,alpha=0.5) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_stars(data=st_as_stars(r_sex_downsampled)) +
  scale_fill_manual(values=c('green','purple'),name='Sex',na.value=NA) +
  annotation_scale() +
  annotation_north_arrow(location='br') +
  xlab("Easting (m)") + ylab("Northing (m)")

g_cytotype = ggplot() + 
  theme_bw() + 
  geom_stars(data=st_as_stars(hill_single_trimmed),show.legend=FALSE,alpha=0.5) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_stars(data=st_as_stars(r_ploidy_level_downsampled)) +
  scale_fill_manual(values=c('red','blue'),name='Cytotype',na.value=NA) +
  annotation_scale() +
  annotation_north_arrow(location='br') +
  xlab("Easting (m)") + ylab("Northing (m)")

ggsave(ggarrange(g_sex, g_cytotype, labels='AUTO',nrow=2,ncol=1),width=5,height=8, file='output_figures/g_sex_cytotype.png')



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


g_n_medium_by_site = ggplot(data_site %>% filter(Point_Type=='Random'), aes(x=year,y=n_medium_trees,group=site_code)) +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  theme_bw() +
  facet_grid(geneticSexID~Ploidy_level) +
  #scale_y_sqrt() +
  xlab('Year') + ylab("# medium saplings per subplot") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
ggsave(g_n_medium_by_site, file=sprintf('output_figures/g_site_medium_%s.pdf', PREFIX_TYPE), width=6,height=6)
ggsave(g_n_medium_by_site, file=sprintf('output_figures/g_site_medium_%s.png', PREFIX_TYPE), width=6,height=6)


g_n_small_by_site = ggplot(data_site %>% filter(Point_Type=='Random'), aes(x=year,y=n_small_trees,group=site_code)) +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  theme_bw() +
  facet_grid(geneticSexID~Ploidy_level) +
  #scale_y_sqrt() +
  xlab('Year') + ylab("# small saplings per subplot") + 
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
  xlab('Year') + ylab("Plot focal tree diameter at breast height (cm)") + 
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
  xlab('Year') + ylab("Basal area density") + 
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
  scale_color_brewer(palette='Set2',name='# re-visits') +
  coord_equal() +
  xlab('Easting (m)') + ylab('Northing (m)') +
  theme(legend.position='bottom')
ggsave(g_counts, file=sprintf('output_figures/g_counts_%s.pdf',PREFIX_TYPE),width=8,height=5)
ggsave(g_counts, file=sprintf('output_figures/g_counts_%s.png',PREFIX_TYPE),width=8,height=5)


g_counts_grid = ggplot(data_site %>% 
                         select(site_code, X.UTM, Y.UTM, Point_Type, Watershed, year) %>% 
                         filter(Point_Type=='Grid') %>% 
                         group_by(site_code, Watershed, X.UTM, Y.UTM) %>% 
                         tally(), aes(x=X.UTM,y=Y.UTM,color=factor(n))) +
  geom_stars(data=st_as_stars(hill_single_trimmed),show.legend=FALSE,alpha=0.5) +
  scale_fill_distiller(palette = "Greys") +
  geom_point() +
  theme_bw() +
  scale_color_viridis_d(option='viridis',name='# re-visits') +
  xlab('Easting (m)') + ylab('Northing (m)') +
  theme(legend.position='bottom') +
  facet_wrap(~Watershed,scales='free')
ggsave(g_counts_grid, file=sprintf('output_figures/g_counts_grid_%s.pdf',PREFIX_TYPE),width=8,height=8)
ggsave(g_counts_grid, file=sprintf('output_figures/g_counts_grid_%s.png',PREFIX_TYPE),width=8,height=8)


#

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
  scale_fill_manual(values=c('blue','red'),na.value = 'gray',name='Cytotype') +
  facet_wrap(~Point_Type) +
  xlab("Sex")
ggsave(g_counts_sex_ploidy, file=sprintf('output_figures/g_counts_sex_ploidy_%s.pdf', PREFIX_TYPE),width=8,height=5)
ggsave(g_counts_sex_ploidy, file=sprintf('output_figures/g_counts_sex_ploidy_%s.png', PREFIX_TYPE),width=8,height=5)



































 
# # track mortality trends
# g_mortality_boxplot = ggplot(data_site, aes(x=factor(year),y=I(frac_adult_dead_w_background_mortality_estimate),color=Ploidy_level)) +
#   geom_boxplot() +
#   facet_wrap(~Point_Type) +
#   theme_bw()
# ggsave(g_mortality_boxplot, file=sprintf('output_figures/g_mortality_boxplot_%s.pdf', PREFIX_TYPE))
# 
# 
# # track damage trends
# g_damage_boxplot = ggplot(data_site, aes(x=factor(year),y=I(frac_adult_damaged),color=Ploidy_level)) +
#   geom_boxplot() +
#   facet_wrap(~Point_Type) +
#   theme_bw()
# ggsave(g_damage_boxplot, file=sprintf('output_figures/g_damage_boxplot_%s.pdf', PREFIX_TYPE))
# 
# 
# 
# 
# # track regen trends- varies by aspect x cytotype, no temporal trend
# g_regen_boxplot = ggplot(data_site, aes(x=factor(year),y=n_small_trees,color=Ploidy_level)) +
#   geom_boxplot() +
#   facet_wrap(~cut(Cos.aspect,breaks=seq(-1,1,by=0.5)) + Point_Type) +
#   theme_bw()
# ggsave(g_regen_boxplot, file=sprintf('output_figures/g_regen_boxplot_%s.pdf', PREFIX_TYPE))
# 
# 
# # track sapling growth: no real trends except north slopes diploids do much better
# g_sapling_boxplot = ggplot(data_site, aes(x=factor(year),y=n_medium_trees,color=Ploidy_level)) +
#   geom_boxplot() +
#   facet_wrap(~cut(Cos.aspect,breaks=seq(-1,1,by=0.5)) + Point_Type) +
#   theme_bw()
# ggsave(g_sapling_boxplot, file=sprintf('output_figures/g_sapling_boxplot_%s.pdf', PREFIX_TYPE))
# 
# 
# 
# 
# 
# ########################################
# # track as individual sites
# g_mortality_by_site = ggplot(data_site, aes(x=factor(year),y=frac_adult_dead_w_background_mortality_estimate,group=site_code)) +
#   geom_line(alpha=0.5,color='red') +
#   #geom_line(stat="smooth",method = "lm",alpha=0.5,color='black') +
#   theme_bw() +
#   facet_wrap(~Point_Type+Ploidy_level) +
#   ylim(0,1)
# ggsave(g_mortality_by_site, file=sprintf('output_figures/g_mortality_by_site_%s.pdf', PREFIX_TYPE),width=12,height=8)
# # track as individual sites
# g_damage_by_site = ggplot(data_site, aes(x=factor(year),y=frac_adult_damaged,group=site_code)) +
#   geom_line(alpha=0.5,color='orange') +
#   geom_line(stat="smooth",method = "lm",alpha=0.5,color='black') +
#   theme_bw() +
#   facet_wrap(~Point_Type+Ploidy_level) +
#   ylim(0,1)
# ggsave(g_damage_by_site, file=sprintf('output_figures/g_damage_by_site_%s.pdf', PREFIX_TYPE),width=12,height=8)
# 
# g_small_trees_by_site = ggplot(data_site, aes(x=factor(year),y=n_small_trees,group=site_code)) +
#   geom_line(alpha=0.5,color='purple') +
#   geom_line(stat="smooth",method = "lm",alpha=0.5,color='black') +
#   theme_bw() +
#   facet_wrap(~Point_Type+Ploidy_level) +
#   scale_y_sqrt(limits=c(0,max(data_site$n_small_trees,na.rm=TRUE)))
# ggsave(g_small_trees_by_site, file=sprintf('output_figures/g_small_trees_by_site_%s.pdf', PREFIX_TYPE),width=12,height=8)
# 
# 
# g_medium_trees_by_site = ggplot(data_site, aes(x=factor(year),y=n_medium_trees,group=site_code)) +
#   geom_line(alpha=0.5,color='blue') +
#   geom_line(stat="smooth",method = "lm",alpha=0.5,color='black') +
#   theme_bw() +
#   facet_wrap(~Point_Type+Ploidy_level) +
#   scale_y_sqrt(limits=c(0,max(data_site$n_medium_trees,na.rm=TRUE)))
# ggsave(g_small_trees_by_site, file=sprintf('output_figures/g_medium_trees_by_site_%s.pdf', PREFIX_TYPE),width=12,height=8)
# 
# 
# 
# # map out where mortality is occurring
# g_map_mortality = ggplot(data_site %>% filter(Point_Type=='Random'), aes(x=X.UTM,y=Y.UTM,color=frac_adult_dead_w_background_mortality_estimate)) +
#   geom_point() +
#   facet_wrap(~year) +
#   theme_bw() +
#   coord_equal() + 
#   scale_color_viridis_c() # why are some coal creek adult deads NAs?
# ggsave(g_map_mortality, file=sprintf('output_figures/g_map_mortality_%s.pdf',PREFIX_TYPE),width=12,height=7)
# 




###############################



# 
# # growth rates seems to be about the same among ploidy levels, but more anticorrelated in diploids
# 
# g_growth_rate_distribution = ggplot(data_site_with_growth %>% filter(Point_Type=='Random' & outlier_growth==FALSE), aes(x=gr_2023_2020,y=gr_2020_2018,color=Ploidy_level)) +
#   facet_wrap(~Ploidy_level) +
#   coord_equal() + 
#   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
#   geom_point() +
#   theme_bw()
# ggsave(g_growth_rate_distribution, file=sprintf('output_figures/g_growth_rate_distribution_%s.pdf',PREFIX_TYPE),width=10,height=7)



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
    geom_polygon(data=pca_this$x[ch_indices,], inherit.aes = FALSE, mapping=aes(x=PC1,y=PC2),fill=NA,color='#AAAAAAAA') +
    theme(legend.position='bottom') +
    xlab(sprintf("PC1 (%.1f%%)",variance_fraction[1])) + ylab(sprintf("PC2 (%.1f%%)",variance_fraction[2]))
  
    # scale_color_manual(values=c('blue','red'),na.value = 'gray',name='Cytotype') + 
    # scale_shape_discrete(name='Sex',na.value=3)
  
  g_pca_loadings = ggplot(pca_this$rotation %>% 
                            as.data.frame %>% 
                            mutate(var=sapply(sapply(gsub("_"," ", row.names(.)),strwrap,width=20),paste,collapse='\n'))) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),arrow = arrow(length = unit(0.1, "inches")),color='gray') +
    geom_label(aes(x=PC1*1.1, y=PC2*1.1, label=var),size=3,alpha=0.5) +
    theme_bw() +
    coord_equal() +
    xlab(sprintf("PC1 (%.2f%%)",variance_fraction[1])) + ylab(sprintf("PC2 (%.2f%%)",variance_fraction[2])) +
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
  select(site_code, year, Ploidy_level, geneticSexID, 
         #Cos.aspect, Elevation, Slope,
         n_small_trees, n_medium_trees, 
         dbh_mean_live,
         basal_area_density_live,
         frac_adult_dead_w_background_mortality_estimate,
         frac_adult_damaged,
         growth_rate#,
         # n_small_dead,
         # n_medium_dead,
         # n_dead_down
  )

rows_na_pca_full = data_pca_all %>%
  select(n_small_trees:growth_rate) %>%
  rowSums %>%
  is.na

data_pca_all = data_pca_all[which(!rows_na_pca_full),]

pca_all = prcomp(data_pca_all %>%
                   select(n_small_trees:growth_rate),
                 center=TRUE, scale=TRUE)

plot_pca(pca_all, data_pca_all, 'all_vars')





data_pca_reduced = data_site_with_growth_for_pca %>% 
  filter(Point_Type=='Random') %>%
  filter(!is.na(Ploidy_level) & !is.na(geneticSexID)) %>%
  select(site_code, year, Ploidy_level, geneticSexID, 
         #Cos.aspect, Elevation, Slope,
         #n_small_trees, 
         n_medium_trees, 
         dbh_mean_live,
         #basal_area_density_live,
         frac_adult_dead_w_background_mortality_estimate,
         #frac_adult_damaged,
         growth_rate#,
         # n_small_dead,
         # n_medium_dead,
         # n_dead_down
  )

rows_na_pca_reduced = data_pca_reduced %>%
  select(n_medium_trees:growth_rate) %>%
  rowSums %>%
  is.na

data_pca_reduced = data_pca_reduced[which(!rows_na_pca_reduced),]

# this is the model we use in the manuscript
pca_reduced = prcomp(data_pca_reduced %>%
                   select(n_medium_trees:growth_rate),
                 center=TRUE, scale=TRUE)

plot_pca(pca_reduced, data_pca_reduced, 'reduced_vars')


with(pca_reduced, 100*(sdev^2)/sum(sdev^2))


######################## 
# demography
df_for_plot_level_model = data_site_with_growth %>%
  filter(Point_Type=='Random') %>%
  # mutate(Ploidy_level=ifelse(is.na(Ploidy_level), "unknown",Ploidy_level),Ploidy_level) %>%
  # mutate(geneticSexID=ifelse(is.na(geneticSexID), "unknown",geneticSexID),geneticSexID) %>%
  mutate(Ploidy_level = factor(Ploidy_level)) %>%
  mutate(geneticSexID = factor(geneticSexID)) %>%
  mutate(dbh_mean_live = ifelse(is.na(dbh_mean_live), dbh_center_live, dbh_mean_live))
  # # filter(!is.na(Ploidy_level)) %>%
  # # filter(!is.na(geneticSexID)) %>%
  # mutate(Ploidy_level = factor(Ploidy_level)) %>%
  # mutate(geneticSexID = factor(geneticSexID))
  



fit_plot_level_model <- function(response_var, family, data, zi=TRUE, title=response_var)
{
  formula_this = sprintf("%s ~ Ploidy_level * geneticSexID + year + Cos.aspect + Elevation + Slope + (1|site_code)", response_var)
  
  m_this = glmmTMB(formula=formula(formula_this), 
                   data=data,
                   ziformula = formula(ifelse(zi==TRUE, "~ Ploidy_level * geneticSexID + I(factor(year))", "~0")),
                   family=family,
                   control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
  
  dharma_this = simulateResiduals(m_this)
  
  preds_this = ggpredict(m_this, 
                         terms=c('year','Ploidy_level','geneticSexID'),
                         type=ifelse(zi==TRUE,'zero_inflated','fixed'))
                         # bias_correction = TRUE)
  plot_effect_this = plot(preds_this,
                          show_data=TRUE,
                          jitter=c(0.1,0),
                          colors=c('blue','red')) +
    labs(color='Cytotype') +
    #scale_color_discrete(name='Cytotype') +
    xlab('Year') + 
    ylab('Predicted value') +
    ggtitle(title)
  
  plot_resid_this = ggplot(data.frame(resid=resid(m_this)),aes(x=resid)) + 
    geom_density() +
    geom_vline(xintercept = 0,color='red') +
    theme_bw()
  
  r2_this = c(r2(m_this), r2_zeroinflated(m_this)) 
  
  return(list(model = m_this,
              summary = summary(m_this),
              dharma = dharma_this,
              preds_this = preds_this,
              plot_effect = plot_effect_this,
              plot_resid = plot_resid_this,
              r2 = r2_this))
}

results_frac_dead = fit_plot_level_model(response_var="frac_adult_dead_w_background_mortality_estimate", 
                              data=df_for_plot_level_model, family=lognormal, zi=TRUE, 
                              title = 'Fraction adult dead')

results_size = fit_plot_level_model(response_var="dbh_mean_live", 
                                         data=df_for_plot_level_model, family=lognormal, zi=FALSE, 
                                         title = 'Focal tree DBH (cm)')

results_n_medium = fit_plot_level_model(response_var="sqrt_n_medium_trees", 
                                    data=df_for_plot_level_model %>% mutate(sqrt_n_medium_trees = sqrt(n_medium_trees)), family=lognormal, zi=TRUE, 
                                    title = 'Number medium saplings (sqrt)')


results_size$summary
results_n_medium$summary
results_frac_dead$summary



g_plot_level = ggarrange( results_size$plot_effect, 
                          results_frac_dead$plot_effect, 
                          results_n_medium$plot_effect,
          nrow=3,labels='AUTO',
          common.legend = TRUE,
          legend = 'bottom')

ggsave(g_plot_level, file=sprintf('output_figures/g_plot_level_temporal_trends_%s.pdf', PREFIX_TYPE),width=6,height=7)
ggsave(g_plot_level, file=sprintf('output_figures/g_plot_level_temporal_trends_%s.png', PREFIX_TYPE),width=6,height=7)



tab_model(results_size$model, file='output_figures/table_model_plot_level_size.html',encoding='UTF-16')

tab_model(results_frac_dead$model, file='output_figures/table_model_plot_level_mortality.html',encoding='UTF-16')

tab_model(results_n_medium$model, file='output_figures/table_model_plot_level_recruitment.html',encoding='UTF-16')

















### OLDER CODE BELOW


# 
# df_demography = data_site_with_growth %>% 
#   filter(Point_Type=='Random') %>%
#   #mutate(growth_rate_mean = rowMeans(cbind(gr_2020_2018, gr_2023_2020, gr_2023_2018,na.rm=TRUE)))
#   select(site_code, year, Ploidy_level, n_medium_trees, n_small_trees, frac_adult_dead_w_background_mortality_estimate, frac_adult_damaged, ) %>% 
#   mutate(n_small_trees = sqrt(n_small_trees), n_medium_trees = sqrt(n_medium_trees)) %>%
#   na.omit
# 
# pca_demography = prcomp(df_demography %>% select(-site_code, -year, -Ploidy_level),center=TRUE,scale=TRUE)
# 
# g_demography_pca_1_2 = pca_demography %>% 
#   ggbiplot(choices=c(1,2), groups=factor(df_demography$Ploidy_level)) +
#   theme_bw() +
#   ggtitle('Random plots no growth outliers')
# g_demography_pca_1_3 = pca_demography %>% 
#   ggbiplot(choices=c(1,3), groups=factor(df_demography$Ploidy_level)) +
#   theme_bw() +
#   ggtitle('Random plots no growth outliers')
# ggsave(ggarrange(g_demography_pca_1_2, g_demography_pca_1_3), file=sprintf('output_figures/g_demography_pca_%s.pdf', PREFIX_TYPE),width=14,height=8)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ######################## 
# # demography
# df_for_regression = data_site_with_growth %>%
#   filter(Point_Type=='Random') %>%
#   mutate(Ploidy_level=factor(Ploidy_level)) %>%
#   mutate(geneticSexID=factor(geneticSexID))
# 
# 
# fit_model <- function(response_var, family=lognormal, data, zi=TRUE)
# {
#   formula_this = sprintf("%s ~ year*(Ploidy_level*geneticSexID) + cow + Cos.aspect + Elevation + (1|site_code)", response_var)
# 
#   m_this = glmmTMB(formula=formula(formula_this), 
#           data=data,
#           ziformula = formula(ifelse(zi==TRUE, "~1", "~0")),
#           family=family,
#           control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
# 
#   dharma_this = simulateResiduals(m_this)
# 
#   preds_this = ggpredict(m_this, terms=c('year','Ploidy_level','geneticSexID'), type='zero_inflated')
#   plot_effect_this = plot(preds_this) +
#     scale_color_manual(values=c('blue','red'))
# 
#   plot_resid_this = ggplot(data.frame(resid=resid(m_this)),aes(x=resid)) + 
#     geom_density() +
#     geom_vline(xintercept = 0,color='red') +
#     theme_bw()
#   
#   r2_this = c(r2(m_this), r2_zeroinflated(m_this)) 
#   
#   return(list(model = m_this,
#               summary = summary(m_this),
#               dharma = dharma_this,
#               plot_effect = plot_effect_this,
#               plot_resid = plot_resid_this,
#               r2 = r2_this))
# }
#   
# results_frac_dead = fit_model(response_var="frac_adult_dead_w_background_mortality_estimate", 
#           data=df_for_regression, family=lognormal, zi=TRUE)
# ggsave(results_frac_dead$plot_effect, file=sprintf('output_figures/g_demography_frac_dead_%s.pdf', PREFIX_TYPE),width=10,height=7)
# 
# results_n_dead = fit_model(response_var="n_adult_dead_w_background_mortality_estimate_per_m2", 
#                               data=df_for_regression, family=lognormal, zi=TRUE)
# ggsave(results_n_dead$plot_effect, file=sprintf('output_figures/g_demography_n_dead_%s.pdf', PREFIX_TYPE),width=10,height=7)
# 
# 
# 
# results_frac_damaged = fit_model(response_var="frac_adult_damaged", 
#                               data=df_for_regression, family=lognormal, zi=TRUE) # sqrt needed for convergence
# ggsave(results_frac_damaged$plot_effect, file=sprintf('output_figures/g_demography_frac_damaged_%s.pdf', PREFIX_TYPE),width=10,height=7)
# 
# results_n_small = fit_model(response_var="n_small_trees", 
#                                  data=df_for_regression, family=lognormal, zi=TRUE)
# ggsave(results_n_small$plot_effect, file=sprintf('output_figures/g_demography_n_small_%s.pdf', PREFIX_TYPE),width=10,height=7)
# 
# results_n_medium = fit_model(response_var="sqrt(n_medium_trees)", 
#                             data=df_for_regression, family=lognormal, zi=TRUE) # sqrt needed for convergence
# ggsave(results_n_medium$plot_effect, file=sprintf('output_figures/g_demography_sqrt_n_medium_%s.pdf', PREFIX_TYPE),width=10,height=7)
# 
# results_growth = fit_model(response_var="dbh_center_live", 
#                              data=df_for_regression %>% filter(outlier_growth==FALSE), family=gaussian, zi=FALSE) # sqrt needed for convergence
# ggsave(results_growth$plot_effect, file=sprintf('output_figures/g_demography_growth_%s.pdf', PREFIX_TYPE),width=10,height=7)
# 













# 
# # get a handle on how the mortality metrics are related
# ggplot(df_for_regression, aes(x=frac_adult_dead_w_background_mortality_estimate, 
#                               y=n_adult_dead_w_background_mortality_estimate_per_m2,color=factor(year))) +
#   geom_point()
# 
# ggplot(df_for_regression, aes(x=frac_adult_dead_w_background_mortality_estimate, 
#                               y=frac_standing_dead,color=factor(year))) +
#   geom_point()
# 
# ggplot(df_for_regression, aes(x=n_adult_dead_w_background_mortality_estimate_per_m2, 
#                               y=frac_standing_dead,color=factor(year))) +
#   geom_point()
# 
# g_mortality_1 = ggplot(df_for_regression, aes(x=year,y=frac_adult_dead_w_background_mortality_estimate, group=site_code)) +
#   geom_line(alpha=0.25) +
#   theme_bw()
# g_mortality_2 = ggplot(df_for_regression, aes(x=year,y=n_adult_dead_w_background_mortality_estimate_per_m2, group=site_code)) +
#   geom_line(alpha=0.25) +
#   theme_bw() +
#   scale_y_sqrt()
# g_mortality_3 = ggplot(df_for_regression, aes(x=year,y=frac_standing_dead, group=site_code)) +
#   geom_line(alpha=0.25) +
#   theme_bw()
# g_mortality_all = ggarrange(g_mortality_1, g_mortality_2, g_mortality_3,nrow=3)
# ggsave(g_mortality_all, file=sprintf('output_figures/g_mortality_all_%s.pdf', PREFIX_TYPE),width=10,height=12)
# 
# 
# 
# # who has the ery high value of 
# df_for_regression %>% 
#   filter(n_adult_dead_w_background_mortality_estimate_per_m2 > 0.4)
# # FOWCL in 2023 - had a lot of dead trees even wihout correct - i guess no action; this is just a biological outlier



# need to predict prob of death (maybe size dependent)
# need to predict 



# inferences
# sex maybe matters
# increased mortality, not much change in recruitment
# growth higher for triploids, and females


  # try to make an IPM
  # assume we have small, medium, largeish, large trees
  # can do growth metrics
  # survival based on counts
  # reproduction based non-size model, with env covariates?
  # do environmental IPM?
  
  
  
  # first build some regression models
  
  
  







# look at overall death/reproduction metric (ratio of recruitment to adult tree density?)

# fit a site-level metric of decline over time,
# or look for a year x climate effect, in a negbin glmmtmb?
# 
# d_random = data_site %>% filter(Point_Type=='Random')
# m_dead = lmer(sqrt(tree_dead_fraction) ~ (Ploidy_level+Slope+Cos.aspect+Elevation+year)^2 + (1|Genotype), 
#      data=d_random)
# 
# m_dead %>% Anova
# visreg(m_dead, xvar='Elevation',by='Ploidy_level',gg=TRUE, overlay=TRUE)
# simulateResiduals(m_dead) %>% plot
# 
# # for the sites with tree level data, try to fit an IPM # we have growth, mortality, recruitment # need to make estimate of whole-plot recruitment and small stem mortality
