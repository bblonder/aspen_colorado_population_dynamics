library(ggbiplot)
library(ggpubr)
library(dplyr)
library(tidyr)

### move to another script
source('ipm parameters.R')

params = read.csv('output_data/ipm_outcomes_lambda.csv')

g_lambda = ggplot(params, aes(color=factor(Ploidy_levelTriploid),x=geneticSexIDM,y=lambda)) + 
  geom_hline(yintercept = 1) +
  geom_point(alpha=0.5) +
  facet_grid(Cos.aspect~n_medium_trees) +
  theme_bw() +
  scale_y_log10() +
  scale_color_manual(values=c('blue','red'))

g_lambda
ggsave(g_lambda, file='output_figures/g_ipm_lambda.pdf')



# so the females and the triploids are doing worse, except in cases where triploids have higher recruitment or other factorsc cause more medium trees



# check age structure
age_structure = read.csv('output_data/sites_age_structure.csv')

age_structure_normalized = apply(age_structure, 1, function(x) {x/sum(x)}) %>% t

g_age_structure_normalized = age_structure_normalized %>% 
  as.data.frame %>%
  mutate(row=row_number()) %>% 
  pivot_longer(!row) %>%
  mutate(dbh_index=as.numeric(gsub("V","",name))) %>%
  mutate(dbh = dbh_range[dbh_index]) %>%  
  ggplot(aes(x=dbh,y=value, group=row)) +
  geom_line(alpha=0.5) +
  theme_bw() +
  scale_y_sqrt() +
  xlab("Stable age structure DBH") +
  ylab("Probability")
ggsave(g_age_structure_normalized, file='output_figures/g_ipm_age_structure_normalized.pdf')




df_sites_for_ipm_joined = read.csv('output_data/sites_data_frame.csv')

g_ipm_size_distribution = ggplot(df_sites_for_ipm_joined, aes(y=lambda, x=dbh_modal, color=geneticSexID)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 1, color='purple') +
  geom_smooth(method='lm') + 
  xlab("Stable age distribution, modal DBH (cm)")
ggsave(g_ipm_size_distribution, file='output_figures/g_ipm_size_distribution.pdf')



# 
# ggplot(df_sites_for_ipm_joined, aes(x=geneticSexID,y=lambda)) + 
#   geom_boxplot()

g_map = ggplot(df_sites_for_ipm_joined, aes(x=X.UTM,y=Y.UTM, color=lambda, shape=geneticSexID)) + 
  geom_point() +
  scale_color_gradient2(low='red',high='blue',mid='gray',midpoint=1) +
  theme_bw()
ggsave(g_map, file='output_figures/g_ipm_map.pdf',width=12,height=12)

g_map_binned = ggplot(df_sites_for_ipm_joined, aes(x=X.UTM,y=Y.UTM, color=lambda_binned, shape=geneticSexID)) + 
  geom_point() +
  scale_color_manual(values=c('red','black','blue')) +
  theme_bw()
ggsave(g_map_binned, file='output_figures/g_ipm_map_binned.pdf',width=12,height=12)

g_lambda_hist = ggplot(df_sites_for_ipm_joined, aes(x=lambda)) +
  geom_histogram(binwidth = 0.002) +
  geom_vline(xintercept = 1,color='purple') +
  theme_bw()
ggsave(g_lambda_hist, file='output_figures/g_lambda_hist.pdf')

g_lambda_binned_hist = ggplot(df_sites_for_ipm_joined, aes(x=lambda_binned,fill=lambda_binned)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(values=c('red','black','blue')) +
  theme(legend.position='none')
ggsave(g_lambda_binned_hist, file='output_figures/g_ipm_lambda_binned_hist.pdf')

# ggplot(df_sites_for_ipm_joined, aes(x=Elevation,y=lambda)) + 
#   geom_point() +
#   #scale_color_gradient2(low='red',high='blue',mid='gray',midpoint=1) +
#   theme_bw() + 
#   facet_wrap(~geneticSexID+Ploidy_level)
# 




pca_all = df_sites_for_ipm_joined %>% 
  ungroup %>% 
  select(Cos.aspect, Slope, Elevation, 
         n_medium_trees, population_density_m2, size_mean, 
         geneticSexID, Ploidy_level,
         lambda) %>% 
  mutate(SexMale = as.numeric(factor(geneticSexID)), PloidyTriploid = as.numeric(factor(Ploidy_level))) %>%
  select(-geneticSexID, -Ploidy_level) %>%
  na.omit %>% 
  prcomp(center=TRUE,scale=TRUE)


g_pc_12 = ggbiplot(pca_all, alpha=0.5, choices=c(1,2)) +
  theme_bw()
g_pc_13 = ggbiplot(pca_all, alpha=0.5, choices=c(1,3)) +
  theme_bw()
ggsave(ggarrange(g_pc_12, g_pc_13,nrow=1,ncol=2), file='output_figures/g_ipm_pca.pdf', width=15,height=7)








# 
# # can we post hoc explain with remotely sensible variables
# library(randomForest)
# randomForest(lambda ~ Cos.aspect + Slope + Elevation + Summer.Insolation, data=df_sites_for_ipm_joined)
# randomForest(lambda ~ geneticSexID + Ploidy_level, data=df_sites_for_ipm_joined)
# randomForest(lambda ~ Cos.aspect + Slope + Elevation + Summer.Insolation + 
#                geneticSexID + Ploidy_level, data=df_sites_for_ipm_joined)
# 
# # only about 45% with remotely sensible variables
# rf_remotely_sensed = randomForest(lambda ~ Cos.aspect + Slope + Elevation +
#                                     Ploidy_level, 
#                                   data=df_sites_for_ipm_joined)
# 
# # how about with everything
# randomForest(lambda ~ Cos.aspect + Slope + Elevation + Summer.Insolation + 
#                geneticSexID + Ploidy_level + size_mean + population_density_m2 + n_medium_trees, data=df_sites_for_ipm_joined)
# 
# 


# could do a separate IPM for each of the grids?


# incorporate stochastic mortality?

# we are using n_medium because the observed rate of recruitment at each site is too low to direclty observe over 2 panels
# and it turns out n_small is not a good predictor
# this in turn will allow us to make per-site models or explore scenarios of variation in herbivory

# feedbacks
# between variables, e.g. if triploids are outcompeted by diploids
# what if there are episodic mortality events
# none of these are currently included

# to consider
# it appears that large individuals have high survival, high enough that eventually all individuals progress into this size class
# need to read about truncation or think about ways to get more mortality (episodic stochastic events?)

# Does recruitment do what I think in the kernel or is the rate too high

# how to account for DBH density at a site - put it in direclty or use a covariate like cos aspect?