library(dplyr)
library(tidyr)
library(ggplot2)
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
  select(-Site_Code,-X.UTM,-Y.UTM)
data_ploidy = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv') %>%
  select(site_code=Site_Code,Ploidy_level,Cos.aspect,Elevation,Slope,X.UTM,Y.UTM,Watershed)
data_site = data_site %>%
  left_join(data_sex, by='site_code') %>%
  left_join(data_ploidy, by='site_code') %>%
  mutate(Point_Type = ifelse(nchar(site_code)==4,'Grid','Random'))


# count up the number of sites that we visited
data_site_by_year = data_site %>%
  group_by(site_code, year) %>%
  tally() %>%
  pivot_wider(id_cols=site_code, names_from=year,values_from=n, values_fill=0) %>%
  pivot_longer(!site_code) %>%
  left_join(data_site %>% select(site_code, Point_Type) %>% unique,by='site_code')

g_which_sampled = ggplot(data_site_by_year, aes(x=name,y=site_code,fill=factor(value))) + geom_tile() +
    facet_wrap(~Point_Type,scales='free') +
    scale_fill_manual(values=c('red','blue'),name='Visited')
ggsave(g_which_sampled, file=sprintf('output_figures/g_which_sampled_%s.pdf',PREFIX_TYPE),width=12,height=15)

data_counts_by_year = data_site %>%
  group_by(year, Point_Type) %>%
  tally

g_total_counts = ggplot(data_counts_by_year, aes(x=factor(year),y=n)) +
  geom_bar(stat='identity') +
  facet_wrap(~Point_Type) +
  theme_bw()
ggsave(g_total_counts, file=sprintf('output_figures/g_total_counts_%s.pdf',PREFIX_TYPE),width=8,height=6)



# why are some sites not marked as having 11 trees surveyed? why 1?
# how is fraction_dead calculated?
t_aspen_counts = table(data_site$n_aspen_trees_surveyed, data_site$year)

t_aspen_counts

write.csv(t_aspen_counts, file=sprintf('output_data/t_aspen_counts_%s.csv', PREFIX_TYPE), row.names=TRUE)
# check on details - this all looks OK to me
write.csv(data_site %>% 
            filter(year==2023 & n_aspen_trees_surveyed > 1 & n_aspen_trees_surveyed < 10), 
          file = sprintf('output_data/t_aspen_counts_details_%s.csv', PREFIX_TYPE), row.names = FALSE)


# track mortality trends
g_mortality_boxplot = ggplot(data_site, aes(x=factor(year),y=I(frac_adult_dead_w_background_mortality_estimate),color=Ploidy_level)) +
  geom_boxplot() +
  facet_wrap(~Point_Type) +
  theme_bw()
ggsave(g_mortality_boxplot, file=sprintf('output_figures/g_mortality_boxplot_%s.pdf', PREFIX_TYPE))


# track damage trends
g_damage_boxplot = ggplot(data_site, aes(x=factor(year),y=I(frac_adult_damaged),color=Ploidy_level)) +
  geom_boxplot() +
  facet_wrap(~Point_Type) +
  theme_bw()
ggsave(g_damage_boxplot, file=sprintf('output_figures/g_damage_boxplot_%s.pdf', PREFIX_TYPE))




# trrack regen trends- varies by aspect x cytotype, no temporal trend
g_regen_boxplot = ggplot(data_site, aes(x=factor(year),y=n_small_trees,color=Ploidy_level)) +
  geom_boxplot() +
  facet_wrap(~cut(Cos.aspect,breaks=seq(-1,1,by=0.5)) + Point_Type) +
  theme_bw()
ggsave(g_regen_boxplot, file=sprintf('output_figures/g_regen_boxplot_%s.pdf', PREFIX_TYPE))


# track sapling growth: no real trends except north slopes diploids do much better
g_sapling_boxplot = ggplot(data_site, aes(x=factor(year),y=n_medium_trees,color=Ploidy_level)) +
  geom_boxplot() +
  facet_wrap(~cut(Cos.aspect,breaks=seq(-1,1,by=0.5)) + Point_Type) +
  theme_bw()
ggsave(g_sapling_boxplot, file=sprintf('output_figures/g_sapling_boxplot_%s.pdf', PREFIX_TYPE))





########################################
# track as individual sites
g_mortality_by_site = ggplot(data_site, aes(x=factor(year),y=frac_adult_dead_w_background_mortality_estimate,group=site_code)) +
  geom_line(alpha=0.5,color='red') +
  geom_line(stat="smooth",method = "lm",alpha=0.5,color='black') +
  theme_bw() +
  facet_wrap(~Point_Type+Ploidy_level) +
  ylim(0,1)
ggsave(g_mortality_by_site, file=sprintf('output_figures/g_mortality_by_site_%s.pdf', PREFIX_TYPE),width=12,height=8)
# track as individual sites
g_damage_by_site = ggplot(data_site, aes(x=factor(year),y=frac_adult_damaged,group=site_code)) +
  geom_line(alpha=0.5,color='orange') +
  geom_line(stat="smooth",method = "lm",alpha=0.5,color='black') +
  theme_bw() +
  facet_wrap(~Point_Type+Ploidy_level) +
  ylim(0,1)
ggsave(g_damage_by_site, file=sprintf('output_figures/g_damage_by_site_%s.pdf', PREFIX_TYPE),width=12,height=8)

g_small_trees_by_site = ggplot(data_site, aes(x=factor(year),y=n_small_trees,group=site_code)) +
  geom_line(alpha=0.5,color='purple') +
  geom_line(stat="smooth",method = "lm",alpha=0.5,color='black') +
  theme_bw() +
  facet_wrap(~Point_Type+Ploidy_level) +
  scale_y_sqrt(limits=c(0,max(data_site$n_small_trees,na.rm=TRUE)))
ggsave(g_small_trees_by_site, file=sprintf('output_figures/g_small_trees_by_site_%s.pdf', PREFIX_TYPE),width=12,height=8)


g_medium_trees_by_site = ggplot(data_site, aes(x=factor(year),y=n_medium_trees,group=site_code)) +
  geom_line(alpha=0.5,color='blue') +
  geom_line(stat="smooth",method = "lm",alpha=0.5,color='black') +
  theme_bw() +
  facet_wrap(~Point_Type+Ploidy_level) +
  scale_y_sqrt(limits=c(0,max(data_site$n_medium_trees,na.rm=TRUE)))
ggsave(g_small_trees_by_site, file=sprintf('output_figures/g_medium_trees_by_site_%s.pdf', PREFIX_TYPE),width=12,height=8)



# map out where mortality is occurring
g_map_mortality = ggplot(data_site %>% filter(Point_Type=='Random'), aes(x=X.UTM,y=Y.UTM,color=frac_adult_dead_w_background_mortality_estimate)) +
  geom_point() +
  facet_wrap(~year) +
  theme_bw() +
  coord_equal() + 
  scale_color_viridis_c() # why are some coal creek adult deads NAs?
ggsave(g_map_mortality, file=sprintf('output_figures/g_map_mortality_%s.pdf',PREFIX_TYPE),width=12,height=7)


# map out where site are missed
g_counts = ggplot(data_site %>% 
         select(site_code, X.UTM, Y.UTM, Point_Type, year) %>% 
         filter(Point_Type=='Random') %>% 
         group_by(site_code, X.UTM, Y.UTM) %>% 
         tally(), aes(x=X.UTM,y=Y.UTM,color=factor(n))) +
  geom_point() +
  theme_bw() +
  scale_color_brewer(palette='Set1')
ggsave(g_counts, file=sprintf('output_figures/g_counts_%s.pdf',PREFIX_TYPE),width=8,height=8)

g_counts_grid = ggplot(data_site %>% 
         select(site_code, X.UTM, Y.UTM, Point_Type, Watershed, year) %>% 
         filter(Point_Type=='Grid') %>% 
         group_by(site_code, Watershed, X.UTM, Y.UTM) %>% 
         tally(), aes(x=X.UTM,y=Y.UTM,color=factor(n))) +
  geom_point() +
  theme_bw() +
  scale_color_brewer(palette='Set1') +
  facet_wrap(~Watershed,scales='free')
ggsave(g_counts_grid, file=sprintf('output_figures/g_counts_grid_%s.pdf',PREFIX_TYPE),width=8,height=8)






###############################
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

ggplot(data_rgr, aes(x=gr_2020_2018,y=gr_2023_2020)) +
  geom_point()

data_site_with_growth = data_site %>%
  left_join(data_rgr,by='site_code') %>%
  mutate(gr_average = (gr_2023_2018 + gr_2023_2020)/2)

# look at stem growth
g_growth_rate_by_site = ggplot(data_site_with_growth, 
                                        aes(x=year, y=dbh_center_live, group=site_code, color=outlier_growth)) + 
  geom_point() +
  geom_line(alpha=0.5) + 
  facet_wrap(~Ploidy_level,scales='free') +
  theme_bw() +
  scale_color_manual(values=c('black','red'))
ggsave(g_growth_rate_by_site, file=sprintf('output_figures/g_growth_rate_by_site_%s.pdf', PREFIX_TYPE),width=10,height=8)

 # find outliers
data_rgr %>%
  filter(outlier_growth==TRUE) %>%
  select(site_code, `2018`,`2020`,`2023`, starts_with('rgr'), starts_with('gr')) %>%
  write.csv(sprintf('output_data/t_growth_rate_high_%s.csv',PREFIX_TYPE))




# growth rates seems to be about the same among ploidy levels, but more anticorrelated in diploids

g_growth_rate_distribution = ggplot(data_site_with_growth %>% filter(Point_Type=='Random' & outlier_growth==FALSE), aes(x=gr_2023_2020,y=gr_2020_2018,color=Ploidy_level)) +
  facet_wrap(~Ploidy_level) +
  coord_equal() + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point() +
  theme_bw()
ggsave(g_growth_rate_distribution, file=sprintf('output_figures/g_growth_rate_distribution_%s.pdf',PREFIX_TYPE),width=10,height=7)










# PCA of demographic variables by year
df_demography = data_site_with_growth %>% 
  filter(outlier_growth==FALSE & Point_Type=='Random') %>%
  select(site_code, year, Ploidy_level, n_medium_trees, n_small_trees, frac_adult_dead_w_background_mortality_estimate, frac_adult_damaged, rgr_2023_2020, rgr_2020_2018) %>% 
  mutate(n_small_trees = sqrt(n_small_trees), n_medium_trees = sqrt(n_medium_trees)) %>%
  na.omit

pca_demography = prcomp(df_demography %>% select(-site_code, -year, -Ploidy_level),center=TRUE,scale=TRUE)

g_demography_pca_1_2 = pca_demography %>% 
  ggbiplot(choices=c(1,2), groups=factor(df_demography$Ploidy_level)) +
  theme_bw() +
  ggtitle('Random plots no growth outliers')
g_demography_pca_1_3 = pca_demography %>% 
  ggbiplot(choices=c(1,3), groups=factor(df_demography$Ploidy_level)) +
  theme_bw() +
  ggtitle('Random plots no growth outliers')
ggsave(ggarrange(g_demography_pca_1_2, g_demography_pca_1_3), file=sprintf('output_figures/g_demography_pca_%s.pdf', PREFIX_TYPE),width=14,height=8)











######################## 
# demography
df_for_regression = data_site_with_growth %>%
  filter(Point_Type=='Random') %>%
  mutate(Ploidy_level=factor(Ploidy_level)) %>%
  mutate(geneticSexID=factor(geneticSexID))


fit_model <- function(response_var, family=lognormal, data, zi=TRUE)
{
  formula_this = sprintf("%s ~ year*(Ploidy_level*geneticSexID) + cow + Cos.aspect + Elevation + (1|site_code)", response_var)

  m_this = glmmTMB(formula=formula(formula_this), 
          data=data,
          ziformula = formula(ifelse(zi==TRUE, "~1", "~0")),
          family=family,
          control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))

  dharma_this = simulateResiduals(m_this)

  preds_this = ggpredict(m_this, terms=c('year','Ploidy_level','geneticSexID'), type='zero_inflated')
  plot_effect_this = plot(preds_this) +
    scale_color_manual(values=c('blue','red'))

  plot_resid_this = ggplot(data.frame(resid=resid(m_this)),aes(x=resid)) + 
    geom_density() +
    geom_vline(xintercept = 0,color='red') +
    theme_bw()
  
  r2_this = c(r2(m_this), r2_zeroinflated(m_this)) 
  
  return(list(model = m_this,
              summary = summary(m_this),
              dharma = dharma_this,
              plot_effect = plot_effect_this,
              plot_resid = plot_resid_this,
              r2 = r2_this))
}
  
results_frac_dead = fit_model(response_var="frac_adult_dead_w_background_mortality_estimate", 
          data=df_for_regression, family=lognormal, zi=TRUE)
ggsave(results_frac_dead$plot_effect, file=sprintf('output_figures/g_demography_frac_dead_%s.pdf', PREFIX_TYPE),width=10,height=7)

results_n_dead = fit_model(response_var="n_adult_dead_w_background_mortality_estimate_per_m2", 
                              data=df_for_regression, family=lognormal, zi=TRUE)
ggsave(results_n_dead$plot_effect, file=sprintf('output_figures/g_demography_n_dead_%s.pdf', PREFIX_TYPE),width=10,height=7)



results_frac_damaged = fit_model(response_var="frac_adult_damaged", 
                              data=df_for_regression, family=lognormal, zi=TRUE) # sqrt needed for convergence
ggsave(results_frac_damaged$plot_effect, file=sprintf('output_figures/g_demography_frac_damaged_%s.pdf', PREFIX_TYPE),width=10,height=7)

results_n_small = fit_model(response_var="n_small_trees", 
                                 data=df_for_regression, family=lognormal, zi=TRUE)
ggsave(results_n_small$plot_effect, file=sprintf('output_figures/g_demography_n_small_%s.pdf', PREFIX_TYPE),width=10,height=7)

results_n_medium = fit_model(response_var="sqrt(n_medium_trees)", 
                            data=df_for_regression, family=lognormal, zi=TRUE) # sqrt needed for convergence
ggsave(results_n_medium$plot_effect, file=sprintf('output_figures/g_demography_sqrt_n_medium_%s.pdf', PREFIX_TYPE),width=10,height=7)

results_growth = fit_model(response_var="dbh_center_live", 
                             data=df_for_regression %>% filter(outlier_growth==FALSE), family=gaussian, zi=FALSE) # sqrt needed for convergence
ggsave(results_growth$plot_effect, file=sprintf('output_figures/g_demography_growth_%s.pdf', PREFIX_TYPE),width=10,height=7)




data_counts_sex_ploidy = df_for_regression %>%
  ungroup %>%
  group_by(Ploidy_level, geneticSexID) %>%
  tally

g_counts_sex_ploidy = ggplot(data_counts_sex_ploidy, aes(x=Ploidy_level,fill=geneticSexID,y=n)) +
  geom_bar(stat='identity',position='dodge') +
  theme_bw() +
  ylab("Number of random sites")
ggsave(g_counts_sex_ploidy, file=sprintf('output_figures/g_counts_sex_ploidy_%s.pdf', PREFIX_TYPE),width=10,height=7)













# get a handle on how the mortality metrics are related
ggplot(df_for_regression, aes(x=frac_adult_dead_w_background_mortality_estimate, 
                              y=n_adult_dead_w_background_mortality_estimate_per_m2,color=factor(year))) +
  geom_point()

ggplot(df_for_regression, aes(x=frac_adult_dead_w_background_mortality_estimate, 
                              y=frac_standing_dead,color=factor(year))) +
  geom_point()

ggplot(df_for_regression, aes(x=n_adult_dead_w_background_mortality_estimate_per_m2, 
                              y=frac_standing_dead,color=factor(year))) +
  geom_point()

g_mortality_1 = ggplot(df_for_regression, aes(x=year,y=frac_adult_dead_w_background_mortality_estimate, group=site_code)) +
  geom_line(alpha=0.25) +
  theme_bw()
g_mortality_2 = ggplot(df_for_regression, aes(x=year,y=n_adult_dead_w_background_mortality_estimate_per_m2, group=site_code)) +
  geom_line(alpha=0.25) +
  theme_bw() +
  scale_y_sqrt()
g_mortality_3 = ggplot(df_for_regression, aes(x=year,y=frac_standing_dead, group=site_code)) +
  geom_line(alpha=0.25) +
  theme_bw()
g_mortality_all = ggarrange(g_mortality_1, g_mortality_2, g_mortality_3,nrow=3)
ggsave(g_mortality_all, file=sprintf('output_figures/g_mortality_all_%s.pdf', PREFIX_TYPE),width=10,height=12)



# who has the ery high value of 
df_for_regression %>% 
  filter(n_adult_dead_w_background_mortality_estimate_per_m2 > 0.4)
# FOWCL in 2023 - had a lot of dead trees even wihout correct - i guess no action; this is just a biological outlier



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
