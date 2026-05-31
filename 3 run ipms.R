library(MASS)
library(dplyr)
library(visreg)
library(DHARMa)
library(MuMIn)
library(ggplot2)
library(car)
library(sjPlot)
library(tidyr)
library(ggeffects)
library(ggpubr)
library(pbapply)
library(insight)
library(progress)
library(Rage)
library(parallel)
library(lme4)
library(sf)
library(elevatr)
library(stars)
library(ggnewscale)
library(GGally)

#source('ipm parameters.R')
NUM_RESAMPLES = 20

set.seed(1) # reproducibility of resampling

if (!file.exists('output_figures'))
{
  dir.create('output_figures')
}

# load in data
transitions_all_filtered_joined = read.csv('output_data/transitions_all_filtered_joined_with_climate_year_final.csv')
#message('could fit only the grid plots too')

# get no NA version for model selection
transitions_all_filtered_joined_no_na = transitions_all_filtered_joined %>%
  mutate(site_type = ifelse(nchar(site_code)==4,'grid','random')) %>%
  mutate(population_density_m2 = 10 / plot_area_m2) %>%
  group_by(site_code, year) %>%
  mutate(size_mean = mean(size, na.rm=TRUE), size_sd = sd(size, na.rm=TRUE)) %>%
  dplyr::select(surv, sizeNext, recruit, size, 
                site_type, site_code, 
                delta_years,
                Ploidy_level, geneticSexID, 
                Cos.aspect, Elevation, 
                STB.0:SWE.2,
                n_medium_trees, 
                population_density_m2, size_mean, size_sd) %>%
  # estimate a size_sd for the cases we don't have one
  ungroup %>%
  mutate(size_sd = ifelse(is.na(size_sd), mean(size_sd, na.rm=TRUE), size_sd)) %>%
  mutate(Ploidy_level=factor(Ploidy_level), geneticSexID=factor(geneticSexID))
#%>%
#  na.omit %>%
#  filter(site_type == 'random')


# transitions_all_filtered_joined_no_na = transitions_all_filtered_joined_no_na %>%
#   left_join(counts_by_plot_year %>% dplyr::select(year, site_code),
#             by=c('year','site_code'))




# get resampled dataset for bootstrapping
transitions_resampled = lapply(1:NUM_RESAMPLES, function(x) {
  transitions_grouped = transitions_all_filtered_joined_no_na %>%
    #filter(site_type!='grid') %>% 
    group_by(Ploidy_level, geneticSexID)
  
  counts_this = transitions_grouped %>% 
    tally
  print(counts_this)
  
  transitions_this_resampled = transitions_grouped %>%
    sample_n(200,replace = TRUE) %>%
    ungroup
  
  return(transitions_this_resampled)
  
  #[sample(1:nrow(transitions_all_filtered_joined_no_na), size=nrow(transitions_all_filtered_joined_no_na), replace=TRUE),]
  })

# utility function
# get_coefficient_distribution <- function(model_list)
# {
#   do.call('rbind',lapply(model_list, function(x) { as.matrix(coef(x)) %>% as.data.frame %>% t} ))
# }

get_coefficient_vector <- function(model)
{
  coefs = coef(model)
  result = coefs; 
  #names(result) = names_this
  return(result)
}

get_coefficient_table <- function(models)
{
  sapply(models, get_coefficient_vector) %>% 
    t %>% 
    as.data.frame %>% 
    mutate(rep=row_number()) %>% 
    pivot_longer(!rep)
}

plot_coefficient_table <- function(models)
{
  get_coefficient_table(models) %>% 
    ggplot(aes(x=value,fill=name)) + 
    geom_density(alpha=0.5) + 
    facet_wrap(~name,scales='free') +
    theme_bw() +
    theme(legend.position='none') +
    geom_vline(xintercept = 0)
}



# fit vital rate models
formula_survival = formula(factor(surv) ~ size + (Ploidy_level + geneticSexID + Cos.aspect + Elevation + I(Elevation^2) + STB.0 + SWE.0 + STB.1 + SWE.1 + population_density_m2))

m_survival_all = glm(formula=formula_survival,
                 data = transitions_all_filtered_joined_no_na,
                 family = 'binomial')
#m_survival_all_reduced = stepAIC(m_survival_all)
#summary(m_survival_all_reduced) 
# 
pdf(file='output_figures/g_ipm_survival.pdf',width=10,height=10)
plot_model(m_survival_all, sort.est=TRUE) +
  theme_bw() +
  geom_hline(yintercept = 1) +
  ggtitle('survival, standardized effect')
simulateResiduals(m_survival_all) %>% plot
dev.off()

tab_model(m_survival_all, file='output_figures/table_m_survival_all.html')
# 
# Anova(m_survival_all_reduced,type=3)

# get resampled version using whatever the final formula was
# m_survival_all_resampled = lapply(transitions_resampled, function(df) {
#   glm(formula=m_survival_all_reduced$formula,
#       data = df,
#       family = binomial())
#   })
m_survival_resampled = lapply(transitions_resampled, function(df) {
  glm(formula=formula_survival,
      data = df,
      family = 'binomial')
})

g_survival_resampled = plot_coefficient_table(m_survival_resampled)
ggsave(g_survival_resampled, file='output_figures/g_survival_resampled.pdf')
ggsave(g_survival_resampled, file='output_figures/g_survival_resampled.png')





#formula_growth = formula(sizeNext ~ size * (Ploidy_level + geneticSexID))

formula_growth = formula(sizeNext ~ size + (Ploidy_level + geneticSexID + Cos.aspect + Elevation + I(Elevation^2) + STB.0 + SWE.0 + STB.1 + SWE.1 + population_density_m2))
m_growth_all = glm(formula=formula_growth, 
                         family='gaussian',
               data = transitions_all_filtered_joined_no_na)
# m_growth_all_reduced = stepAIC(m_growth_all)
# summary(m_growth_all_reduced)
# 
# 
pdf(file='output_figures/g_ipm_growth.pdf',width=10,height=10)
plot_model(m_growth_all, type='est', sort.est=TRUE, rm.terms='size') +
  theme_bw() +
  geom_hline(yintercept = 0) +
  ggtitle('growth, standardized effect')
simulateResiduals(m_growth_all) %>% plot
dev.off()

tab_model(m_growth_all, file='output_figures/table_m_growth_all.html')

# Anova(m_growth_all_reduced,type=3)

# m_growth_all_resampled = lapply(transitions_resampled, function(df) {
#   glm(formula=m_growth_all_reduced$formula,
#       data = df)
# })
m_growth_resampled = lapply(transitions_resampled, function(df) {
  glm(formula=formula_growth,
         family='gaussian',
      data = df)
})

g_growth_resampled = plot_coefficient_table(m_growth_resampled)
ggsave(g_growth_resampled,file='output_figures/g_growth_resampled.pdf')
ggsave(g_growth_resampled,file='output_figures/g_growth_resampled.png')






counts_recruit = transitions_all_filtered_joined_no_na %>%
  group_by(site_code, site_type, delta_years, year, n_medium_trees, Ploidy_level, geneticSexID, Cos.aspect, Elevation, STB.0, SWE.0, STB.1, SWE.1, population_density_m2) %>%
  tally(recruit) %>%
  #mutate(log_one_plus_n_med_per_pop_dens = log(1+n_medium_trees/population_density_m2)) %>%
  #mutate(has_nonzero_n_med = factor(n_medium_trees > 0)) %>%
  as.data.frame %>%
  mutate(has_medium = factor(n_medium_trees > 0))

counts_recruit_resampled = lapply(1:NUM_RESAMPLES, function(x) {
  counts_recruit_grouped = counts_recruit %>%
    #filter(site_type!='grid') %>% 
    group_by(Ploidy_level, geneticSexID)
  
  counts_recruit_this = counts_recruit_grouped %>% 
    tally
  print(counts_recruit_this)

  counts_recruit_this_resampled = counts_recruit_grouped %>%
    sample_n(50,replace = TRUE) %>%
    ungroup

  return(counts_recruit_this_resampled)
  #counts_recruit[sample(1:nrow(counts_recruit), size=nrow(counts_recruit), replace=TRUE),]
})


formula_recruit_count = formula(n ~ n_medium_trees + (Ploidy_level + geneticSexID + Cos.aspect + Elevation + I(Elevation^2) + STB.0 + SWE.0 + STB.1 + SWE.1 + population_density_m2))

m_recruit_count_all = glm(formula = formula_recruit_count,
                          data=counts_recruit,
                          family = 'poisson')
# m_recruit_count_all_reduced = stepAIC(m_recruit_count_all)
# summary(m_recruit_count_all_reduced)
# 
pdf(file='output_figures/g_ipm_recruit_count.pdf',width=10,height=10)
plot_model(m_recruit_count_all, sort.est=TRUE) +
  theme_bw() +
  geom_hline(yintercept = 1) +
  ggtitle('recruit count, standardized effect')
simulateResiduals(m_recruit_count_all) %>% plot
dev.off()

tab_model(m_recruit_count_all, file='output_figures/table_m_recruit_count_all.html')
# 
# Anova(m_recruit_count_all_reduced,type=3)


# m_recruit_count_all_resampled = lapply(counts_recruit_all_resampled, function(df) {
#   glm(formula=m_recruit_count_all_reduced$formula,
#       family=poisson,
#       data = df)
# })
m_recruit_count_resampled = lapply(counts_recruit_resampled, function(df) {
  glm(formula=formula_recruit_count,
      family='poisson',
      data = df)
})

g_recruit_count_resampled = plot_coefficient_table(m_recruit_count_resampled)
ggsave(g_recruit_count_resampled, file='output_figures/g_recruit_count_resampled.pdf')
ggsave(g_recruit_count_resampled, file='output_figures/g_recruit_count_resampled.png')




# add quadratic elevation term, lose the big outliers
# m_n_medium_all = glm.nb(I(n_medium_trees / population_density_m2) ~ (Ploidy_level + geneticSexID + Cos.aspect + Elevation + I(Elevation^2) + STB.0 + SWE.0 + STB.1 + SWE.1 + population_density_m2),
#                         data=counts_recruit_all)


formula_prob_medium = formula(has_medium ~ (Ploidy_level + geneticSexID + Cos.aspect + Elevation + I(Elevation^2) + STB.0 + SWE.0 + STB.1 + SWE.1 + population_density_m2))

m_prob_medium_all = glm(formula_prob_medium, data=counts_recruit, family='binomial')


formula_n_medium_count = formula(I(n_medium_trees / population_density_m2) ~ (Ploidy_level + geneticSexID + Cos.aspect + Elevation + I(Elevation^2) + STB.0 + SWE.0 + STB.1 + SWE.1 + population_density_m2))

# get a theta estimate
# formula_n_medium = formula(I(n_medium_trees / population_density_m2) ~ (Ploidy_level + geneticSexID + Cos.aspect + Elevation + I(Elevation^2) + STB.0 + SWE.0 + STB.1 + SWE.1 + population_density_m2))
# m_n_medium_all = glm.nb(formula_n_medium,
#                         data=counts_recruit)

m_n_medium_count_all = glm.nb(formula_n_medium_count,
                        data=counts_recruit %>% filter(n_medium_trees>0))







# m_n_medium_all = cv.glmnet(formula = formula_n_medium,
#           data=counts_recruit,
#           family = negative.binomial(theta=m_n_medium_all_base$theta))
# 
# # this is effectively a no-intercept model for n_medium
# counts_recruit_all_nonzero = counts_recruit_all %>%
#   filter(has_nonzero_n_med==TRUE)
# 
# formula_n_medium = formula(log_one_plus_n_med_per_pop_dens ~ (Ploidy_level + geneticSexID + Cos.aspect + Elevation + STB.0 + SWE.0 + STB.1 + SWE.1 + population_density_m2))
# m_n_medium_all = glm(formula_n_medium,
#                  family=Gamma,
#                  data=counts_recruit_all_nonzero)
# 
# m_n_medium_all_reduced = stepAIC(m_n_medium_all)
# summary(m_n_medium_all_reduced)
# 
# pdf(file='output_figures/g_ipm_n_medium.pdf',width=10,height=10)
# plot_model(m_n_medium_all, sort.est=TRUE) +
#   theme_bw() +
#   geom_hline(yintercept = 1) +
#   ggtitle('n medium, standardized effect')
# simulateResiduals(m_n_medium_all) %>% plot
# dev.off()
# 
# tab_model(m_n_medium_all, file='output_figures/table_m_n_medium_all.html')

pdf(file='output_figures/g_ipm_prob_medium.pdf',width=10,height=10)
plot_model(m_prob_medium_all, sort.est=TRUE) +
  theme_bw() +
  geom_hline(yintercept = 1) +
  ggtitle('prob medium, standardized effect')
simulateResiduals(m_prob_medium_all) %>% plot
dev.off()

tab_model(m_prob_medium_all, file='output_figures/table_m_prob_medium_all.html')

pdf(file='output_figures/g_ipm_n_medium_count.pdf',width=10,height=10)
plot_model(m_n_medium_count_all, sort.est=TRUE) +
  theme_bw() +
  geom_hline(yintercept = 1) +
  ggtitle('prob medium, standardized effect')
simulateResiduals(m_n_medium_count_all) %>% plot
dev.off()

tab_model(m_n_medium_count_all, file='output_figures/table_m_n_medium_count_all.html')



# Anova(m_n_medium_all_reduced,type=3)
# 
# # m_n_medium_all_resampled = lapply(counts_recruit_all_resampled, function(df) {
# #   glm(formula=m_n_medium_all_reduced$formula,
# #       family=poisson,
# #       data = df %>%
# #         filter(has_nonzero_n_med==TRUE))
# # })


# m_n_medium_all_resampled = lapply(counts_recruit_all_resampled, function(df) {
#   glm.nb(formula=formula(m_n_medium_all_reduced),
#       data = df)
# })

# GOOD CODE BELOW
# m_n_medium_resampled = lapply(counts_recruit_resampled, function(df) {
#   m_n_medium_all_base = glm.nb(formula_n_medium,
#                                data=df)
# })
# 
# g_n_medium_resampled = plot_coefficient_table(m_n_medium_resampled)
# ggsave(g_n_medium_resampled, file='output_figures/g_n_medium_resampled.pdf')
# ggsave(g_n_medium_resampled, file='output_figures/g_n_medium_resampled.png')

m_prob_medium_resampled = lapply(counts_recruit_resampled, function(df) {
  m = glm(formula_prob_medium, data=df, family='binomial')
})

g_prob_medium_resampled = plot_coefficient_table(m_prob_medium_resampled)
ggsave(g_prob_medium_resampled, file='output_figures/g_prob_medium_resampled.pdf')
ggsave(g_prob_medium_resampled, file='output_figures/g_prob_medium_resampled.png')

m_n_medium_count_resampled = lapply(counts_recruit_resampled, function(df) {
  m = glm.nb(formula_n_medium_count,
                               data=df %>% filter(n_medium_trees>0))
})

g_n_medium_count_resampled = plot_coefficient_table(m_n_medium_count_resampled)
ggsave(g_n_medium_count_resampled, file='output_figures/g_n_medium_count_resampled.pdf')
ggsave(g_n_medium_count_resampled, file='output_figures/g_n_medium_count_resampled.png')

# 
# formula_prob_medium = formula(has_nonzero_n_med ~ (Ploidy_level + geneticSexID + Cos.aspect + Elevation+ STB.0 + SWE.0 + STB.1 + SWE.1 + population_density_m2))
# m_prob_medium_all = glm(formula_prob_medium, 
#                      family=binomial,
#                      data=counts_recruit_all)
# m_prob_medium_all_reduced = stepAIC(m_prob_medium_all)
# summary(m_prob_medium_all)

# 
# pdf(file='output_figures/g_ipm_prob_medium.pdf',width=10,height=10)
# plot_model(m_prob_medium_all_reduced, sort.est=TRUE) + 
#   theme_bw() + 
#   geom_hline(yintercept = 1) +
#   ggtitle('n medium, standardized effect')
# try(plot_model(m_prob_medium_all_reduced, type='int'))
# simulateResiduals(m_prob_medium_all_reduced) %>% plot
# dev.off()
# 
# 
# tab_model(m_prob_medium_all_reduced, file='output_figures/table_model_tree_level_prob_medium.html')
# 
# Anova(m_prob_medium_all_reduced,type=3)
# 
# m_prob_medium_all_resampled = lapply(counts_recruit_all_resampled, function(df) {
#   glm(formula=m_prob_medium_all_reduced$formula,
#       family=binomial,
#       data = df)
# })



# 
# # fit size residuals
# m_growth_resid2 = resid(m_growth_all)^2
# m_growth_size = m_growth_all$data[complete.cases(m_growth_all$data[,all.vars(m_growth_all$formula),with=FALSE]),]$size
# df_m_growth = data.frame(resid2=m_growth_resid2, size=m_growth_size)
# 
# # ok to ditch outliers?
# m_growth_size_variance = glm(I(sqrt(resid2))~size,data=df_m_growth %>% filter(resid2 < 30),family = Gamma)
# #plot(m_growth_size_variance)
# summary(m_growth_size_variance)
# 
# pdf(file='output_figures/g_ipm_growth_size_variance.pdf',width=10,height=10)
# plot_model(m_growth_size_variance, sort.est=TRUE) + 
#   theme_bw() + 
#   geom_hline(yintercept = 1) +
#   ggtitle('growth size variance')
# #plot_model(m_recruit_count_all, type='int')
# simulateResiduals(m_growth_size_variance) %>% plot
# dev.off()
# 
# tab_model(m_growth_size_variance, file='output_figures/table_model_tree_level_growth_size_variance.html')
# 
# Anova(m_growth_size_variance,type=3)
# 
# m_growth_size_variance_resampled = lapply(m_growth_all_resampled, function(m) {
#   m_growth_resid2_this = resid(m)^2
#   m_growth_size_this = m$data[complete.cases(m$data[,all.vars(m$formula),with=FALSE]),]$size
#   df_m_growth_this = data.frame(resid2=m_growth_resid2_this, size=m_growth_size_this)
#   
#   # ok to ditch outliers?
#   m_growth_size_variance_this = glm(I(sqrt(resid2))~size,data=df_m_growth_this %>% filter(resid2 < 30),family = Gamma)
# })

m_growth_size_variance_all = 
y_pred_all = predict(m_growth_all,
                 newdata=transitions_all_filtered_joined_no_na, 
                 na.action=na.pass)
y_obs_all = transitions_all_filtered_joined_no_na$sizeNext

df_m_growth_all = data.frame(resid2=as.numeric((y_pred_all - y_obs_all)^2), 
                              size=transitions_all_filtered_joined_no_na$size)
  
m_growth_size_variance_all = glm(I(sqrt(resid2))~size,data=df_m_growth_all,family = Gamma)


pdf(file='output_figures/g_m_growth_size_variance_all.pdf',width=10,height=10)
plot_model(m_growth_size_variance_all, sort.est=TRUE) +
  theme_bw() +
  geom_hline(yintercept = 1) +
  ggtitle('growth size variance')
simulateResiduals(m_growth_size_variance_all) %>% plot
dev.off()

tab_model(m_growth_size_variance_all, file='output_figures/table_m_growth_size_variance_all.html')




m_growth_size_variance_resampled = lapply(1:length(m_growth_resampled), function(i) {
  y_pred = predict(m_growth_resampled[[i]],
                   newdata=transitions_resampled[[i]], 
                   na.action=na.pass)
  y_obs = transitions_resampled[[i]]$sizeNext

  df_m_growth_this = data.frame(resid2=as.numeric((y_pred - y_obs)^2), 
                                size=transitions_resampled[[i]]$size)

  m_growth_size_variance_this = glm(I(sqrt(resid2))~size,data=df_m_growth_this,family = Gamma)
})

g_growth_size_variance_resampled = plot_coefficient_table(m_growth_size_variance_resampled)
ggsave(g_growth_size_variance_resampled, file='output_figures/g_growth_size_variance_resampled.pdf')
ggsave(g_growth_size_variance_resampled, file='output_figures/g_growth_size_variance_resampled.png')


pval_anova = rbind(
  Anova(m_survival_all,type=3) %>%
    as.data.frame %>%
    dplyr::select(pvalue=`Pr(>Chisq)`) %>%
    mutate(xvar=row.names(.)) %>%
    mutate(var='m_survival_all'),
  
  Anova(m_growth_all,type=3) %>%
    as.data.frame %>%
    dplyr::select(pvalue=`Pr(>Chisq)`) %>%
    mutate(xvar=row.names(.)) %>%
    mutate(var='m_growth_all'),
  
  # Anova(m_n_medium_all,type=3) %>%
  #   as.data.frame %>%
  #   dplyr::select(pvalue=`Pr(>Chisq)`) %>%
  #   mutate(xvar=row.names(.)) %>%
  #   mutate(var='m_n_medium_all'),

  Anova(m_prob_medium_all,type=3) %>%
    as.data.frame %>%
    dplyr::select(pvalue=`Pr(>Chisq)`) %>%
    mutate(xvar=row.names(.)) %>%
    mutate(var='m_prob_medium_all'),
    
  Anova(m_n_medium_count_all,type=3) %>%
    as.data.frame %>%
    dplyr::select(pvalue=`Pr(>Chisq)`) %>%
    mutate(xvar=row.names(.)) %>%
    mutate(var='m_n_medium_count_all'),
  
  Anova(m_recruit_count_all,type=3) %>%
    as.data.frame %>%
    dplyr::select(pvalue=`Pr(>Chisq)`) %>%
    mutate(xvar=row.names(.)) %>%
    mutate(var='m_recruit_count_all'),
  
  Anova(m_growth_size_variance_all,type=3) %>%
    as.data.frame %>%
    dplyr::select(pvalue=`Pr(>Chisq)`) %>%
    mutate(xvar=row.names(.)) %>%
    mutate(var='m_growth_size_variance_all')
)
row.names(pval_anova) = NULL
pval_anova$pvalue=format_p(pval_anova$pvalue,stars=TRUE)
pval_anova_wide = pivot_wider(pval_anova, names_from = 'var', values_from='pvalue')

write.csv(pval_anova_wide, file='output_figures/pval_anova_wide.csv',row.names = FALSE)


plot_tree_level_effect <- function(m,label)
{
  plot(ggpredict(m,terms=list(geneticSexID=c('female','male'),Ploidy_level=c('diploid','triploid')))) +
    labs(color='Ploidy level') + ggtitle(label) + ylab('') + xlab('Sex') +
    theme_bw() +
    scale_color_manual(values=c('blue','red'))
}

g_tree_level_effect = ggarrange(
  plot_tree_level_effect(m_growth_all, 'Size progression'),
  plot_tree_level_effect(m_survival_all, 'Probability of survival'),
  plot_tree_level_effect(m_recruit_count_all, 'Number of recruits'),
  plot_tree_level_effect(m_prob_medium_all, 'Probability of medium stem'),
  plot_tree_level_effect(m_n_medium_count_all, 'Number medium stems / population density'),

  #plot_tree_level_effect(m_growth_size_variance_all, 'Size variance'),
  align='hv',nrow=3,ncol=2,labels='AUTO', common.legend = TRUE, legend='bottom')
ggsave(g_tree_level_effect, file='output_figures/g_tree_level_effect.pdf',width=6,height=6.5)
ggsave(g_tree_level_effect, file='output_figures/g_tree_level_effect.png',width=6,height=6.5)



















update_coefficients_full <- function(coefficients, xvar, other_vars)
{
  names_all = names(coefficients) 
  indices_slope = grepl(xvar, names_all)
  
  terms_intercept = coefficients[!indices_slope]
  terms_slope = coefficients[indices_slope]
  # get rid of the xvar notation
  names(terms_slope) = gsub("^:","", gsub(xvar, "1", names(terms_slope)))
  # get rid of the intercept
  names(terms_intercept) = gsub("\\(Intercept\\)","1",names(terms_intercept))
  
  # rewrite variable names according to values
  for (i in 1:length(other_vars))
  {
    names(terms_intercept) = gsub(names(other_vars)[i], other_vars[i], names(terms_intercept))
    names(terms_slope) = gsub(names(other_vars)[i], other_vars[i], names(terms_slope))
  }
  # replace colons with multiplication
  names(terms_intercept) = gsub(":","*", names(terms_intercept))
  names(terms_slope) = gsub(":","*", names(terms_slope))
  
  # calculate combined new values
  string_intercept = paste(paste(names(terms_intercept),terms_intercept, sep="*"),collapse="+")
  string_slope = paste(paste(names(terms_slope),terms_slope, sep="*"),collapse="+")
  
  value_intercept = eval(parse(text=string_intercept))
  value_slope = eval(parse(text=string_slope))
  
  return(c(value_intercept, value_slope))
  #return(list(terms_intercept, terms_slope, string_intercept,string_slope, value_intercept, value_slope))
}







kernel_for_plotting <- function(minsize, maxsize, m, k)
{
  dbh_range = seq(minsize, maxsize, length.out=m)
  
  k_long = k %>% 
    as.matrix %>% 
    as.data.frame %>% 
    mutate(row=row_number()) %>% 
    pivot_longer(!row) %>%
    mutate(sizeTo=dbh_range[row]) %>%
    mutate(sizeFrom=dbh_range[as.numeric(gsub("V","",name))])
}

plot_kernel <- function(minsize, maxsize, m, k, name)
{
  ggplot(kernel_for_plotting(minsize, maxsize, m, k), 
         aes(x=sizeFrom,y=sizeTo,fill=value)) +
    geom_raster() +
    theme_bw() +
    scale_fill_viridis_c(option='C',name='Value') +
    xlab('Size (cm)') + ylab('Size next (cm)') +
    #geom_abline(slope=1,intercept=0,color='white',alpha=0.5) +
    coord_equal() +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(sprintf('%s',name)) 
}

# predict_n_medium <- function(current_conditions, m_n_medium_this)
# {
#   coef_n_medium = update_coefficients_full(get_coefficient_vector(m_n_medium_this), xvar='size', other_vars = current_conditions)
#   
#   # inverse link of neg bin regression
#   n_medium = exp(coef_n_medium) * current_conditions$population_density_m2
#   
#   return(n_medium)
# }


predict_n_medium <- function(current_conditions,
                             m_prob_medium_this = m_prob_medium_all,
                             m_n_medium_count_this = m_n_medium_count_all)
{
  # initial recruits set by regression
  # hurdle model (binomial family)
  coef_prob_medium = update_coefficients_full(coef(m_prob_medium_this), xvar='size', other_vars = current_conditions)
  #print(coef_prob_medium)
  prob_medium = exp(coef_prob_medium)/(1 + exp(coef_prob_medium))
  if (prob_medium >= 0.5)
  {
    # negbin family, then need to
    mu = update_coefficients_full(coef(m_n_medium_count_this), xvar='size', other_vars = current_conditions)
    # exp link
    n_S = (exp(mu))*current_conditions$population_density_m2
  }
  else
  {
    n_S = 0
  }

  return(n_S)
}
# 
# clamp <- function(x, minval, maxval)
# {
#   x[x<minval] = minval
#   x[x>maxval] = maxval
#   
#   return(x)
# }


run_model <- function(prefix_this='test',
                      minsize=5, # L
                      maxsize=60, # U
                      m=100, # mesh points
                      n_iterations=100, # time steps
                      population_density_initial=0.2, # initial density
                      size_mean_initial = mean(transitions_all_filtered_joined_no_na$size, na.rm=TRUE), # initial tagged tree size distribution
                      size_sd_initial = sd(transitions_all_filtered_joined_no_na$size, na.rm=TRUE), # initial tagged tree size distribution
                      Ploidy_leveltriploid='0',
                      #Ploidy_leveldiploid='0',
                      Ploidy_levelunknown='0',
                      geneticSexIDmale='0',
                      #geneticSexIDfemale='0',
                      geneticSexIDunknown='0',
                      Cos.aspect=-1,
                      Elevation=3000,
                      STB.sequence=rep(0,n_iterations+1), # need one extra point for the time-1 lag
                      SWE.sequence=rep(100, n_iterations+1), # need one extra point for the time-1 lag
                      n_S_lag_factor = 0.2,
                      n_S_scale_factor = 1,
                      
                      # regression models
                      m_growth_this = m_growth_all,
                      m_prob_medium_this = m_prob_medium_all,
                      m_n_medium_count_this = m_n_medium_count_all,
                      m_recruit_count_this = m_recruit_count_all,
                      m_survival_this = m_survival_all,
                      m_growth_size_variance_this = m_growth_size_variance_all,
                      progress=TRUE
)
{
  # conditions
  current_conditions = list(Ploidy_leveltriploid=Ploidy_leveltriploid,
                            #Ploidy_leveldiploid=Ploidy_leveldiploid,
                            Ploidy_levelunknown=Ploidy_levelunknown,
                            geneticSexIDmale=geneticSexIDmale,
                            #geneticSexIDfemale=geneticSexIDfemale,
                            geneticSexIDunknown=geneticSexIDunknown,
                            Cos.aspect=Cos.aspect,
                            Elevation=Elevation,
                            STB.0=STB.sequence[1+1],
                            SWE.0=SWE.sequence[1+1],
                            STB.1=STB.sequence[1],
                            SWE.1=SWE.sequence[1],
                            population_density_m2=population_density_initial)#, 
                            #delta_years=delta_years)
  
  # set up bins
  size_bins = seq(minsize, maxsize, length.out=m)
  
  n_A = sapply(size_bins, function(x) {
    size = dnorm(x, mean=size_mean_initial, sd=size_sd_initial)
    return(size)
  })
  # normalize the size to a population density of the requested value
  n_A = n_A * current_conditions$population_density_m2 / sum(n_A)
  
  # get initial number of recruits
  n_S_initial = predict_n_medium(current_conditions, 
                                 m_prob_medium_this = m_prob_medium_this, 
                                 m_n_medium_count_this = m_n_medium_count_this)
  n_S = n_S_initial
  
  # number of adult trees
  n_A_time_series = matrix(NA, nrow=n_iterations, ncol=m)
  # number of saplings
  n_S_time_series = rep(NA, n_iterations)
  n_S_component_time_series = data.frame(matrix(NA, nrow=n_iterations, ncol=2))
  names(n_S_component_time_series) = c('total','recruits')
  # kernel
  P_initial = NULL
  P_final = NULL

  ### do iteration
  if (progress==TRUE)
  {
    pb = progress_bar$new(total=n_iterations)
  }
  for (i in 1:n_iterations)
  {
    n_A_time_series[i,] = n_A
    n_S_time_series[i] = n_S
    
    # update population density
    current_conditions$population_density_m2 = sum(n_A)
    
    # update environment
    current_conditions$STB.0=STB.sequence[i+1]
    current_conditions$SWE.0=SWE.sequence[i+1]
    current_conditions$STB.1=STB.sequence[i]
    current_conditions$SWE.1=SWE.sequence[i]
    
    # get coefficients on link scale
    # binomial
    coef_survival = update_coefficients_full(get_coefficient_vector(m_survival_this), xvar='size', other_vars=current_conditions)
    # gaussian
    coef_growth = update_coefficients_full(get_coefficient_vector(m_growth_this), xvar='size', other_vars=current_conditions)
    # poisson
    coef_recruit_count = update_coefficients_full(get_coefficient_vector(m_recruit_count_this), xvar='n_medium_trees', other_vars=current_conditions)
    # gamma
    coef_size_variance = as.numeric(get_coefficient_vector(m_growth_size_variance_this))

    # survival
    sx <- function(x) {
      xbeta = coef_survival[1] + x*coef_survival[2]
      mu = exp(xbeta)/(1 + exp(xbeta)) #
      return(mu)
    }
    # growth
    gxy<-function(x,y) {
      xbeta_sd <- coef_size_variance[1] + x*coef_size_variance[2]
      mu_sd = 1/xbeta_sd # gamma inv link
      sigmax = sqrt(mu_sd);
      
      mux<-coef_growth[1] + x*coef_growth[2]
      fac1<-sqrt(2*pi)*sigmax;
      fac2<-((y-mux)^2)/(2*sigmax^2);
      return(exp(-fac2)/fac1);
    }
    
    pxy<-function(x,y) { 
      return(sx(x)*gxy(x,y)) 
    }
    
    Kyx<-function(y,x) {
      xeval<-max(x,minsize); xeval<-min(xeval,maxsize);
      yeval<-max(y,minsize); yeval<-min(yeval,maxsize);
      return(pxy(xeval,yeval))
    }
    
    bigmatrix<-function(n) {
      # upper and lower integration limits
      L<-minsize; U<-maxsize;
      
      # boundary points b and mesh points y
      b<-L+c(0:n)*(U-L)/n;
      y<-0.5*(b[1:n]+b[2:(n+1)]);
      h<-(U-L)/n
      
      # loop to construct the matrix
      M<-matrix(0,n,n);
      for (i in 1:n){
        for(j in 1:n){
          M[i,j]<-Kyx(y[i],y[j])
        }
      }
      # Kyx_vectorized<- Vectorize(Kyx)
      # M = outer(1:n,1:n,function(i,j){ Kyx_vectorized(y[i],y[j]) } )
      M<-(U-L)*M/n;
      return(list(matrix=M,meshpts=y, h=h)); 
    }
    
    # get survival growth kernel
    P = bigmatrix(m)
    
    # define recruit size distribution
    c_R <- function(x)
    {
      dunif(x, min=5,max=5.5) # uniform distribution of sizes for recruits
    }
  
    # iterate saplings
    # allow to be regulated by density dependence of adults alone
    # equilibrium value to be multiplied by a scale factor
    n_S_star = n_S_scale_factor * predict_n_medium(current_conditions,
                                 m_prob_medium_this = m_prob_medium_this,
                                 m_n_medium_count_this = m_n_medium_count_this) # n_medium_factor is to account for browsing, etc.
    
    # introduce a lag factor to allow n_S to vary more smoothly
    # experimental hack no evidence for this
    n_S_next = (1-n_S_lag_factor)*n_S + n_S_lag_factor*n_S_star 
    
    # update this with the new n_S
    # assume we can't progress more trees than we currently have
    n_recruits = min(n_S_star, exp(coef_recruit_count[1] + coef_recruit_count[2] * n_S)) # poisson link
  
    # iterate tagged trees
    n_A_next = P$matrix %*% n_A + P$h * c_R(P$meshpts) * n_recruits  
    
    # store iterated value
    n_A = n_A_next
    n_S = n_S_next
    
    # this logging has to come here based on how i organized... (awkward)
    n_S_component_time_series[i,'star'] = n_S_star
    n_S_component_time_series[i,'recruits'] = n_recruits
    
    # save kernels
    if (i==1)
    {
      P_initial = P$matrix
    }
    if (i==n_iterations)
    {
      P_final = P$matrix
    }
    
    if (progress==TRUE)
    {
      pb$tick()
    }
  }
  
  # compute statistics of final state
  size_mean_final = sum(n_A_time_series[n_iterations,] * size_bins) / sum(n_A_time_series[n_iterations,])
  #sx5_final = sx(5)
  population_density_final = sum(n_A_time_series[n_iterations,])
  
  df_n_density = data.frame(time=1:nrow(n_A_time_series), value=apply(n_A_time_series, 1, sum),type='A') %>%
    rbind(data.frame(time=1:length(n_S_time_series),value=n_S_time_series,type='S'))
  
  df_nt = n_A_time_series %>% 
    as.data.frame
  names(df_nt) = round(P$meshpts,3)
  df_nt = df_nt %>%
    mutate(time=row_number()) %>% 
    pivot_longer(!time) %>%
    mutate(size=as.numeric(name))
  
  n_A_final = df_n_density %>%
    filter(time==n_iterations & type=='A') %>%
    pull(value)  
  n_S_final = df_n_density %>%
    filter(time==n_iterations & type=='S') %>%
    pull(value)
  
  longevity_90_initial = NA
  longevity_90_final = NA
  try(longevity_90_initial <- longevity(matU = P_initial, start = 1, lx_crit = 0.1))
  try(longevity_90_final <- longevity(matU = P_final, start = 1, lx_crit = 0.1))
  
  n_A_final_ten = df_n_density %>% filter(type=='A') %>% tail(10)
  cv_n_A_final_ten = sd(n_A_final_ten$value)/mean(n_A_final_ten$value)
  
  # get density of 30 cm trees
  n_A_final_temp = n_A_time_series[n_iterations,]
  n_A_final_gt_30_cm = sum(n_A_final_temp[ which(size_bins>=30)])
  
  # get basal area density  in m2/m2
  basal_area_density = sum(pi*(size_bins/2/100)^2 * n_A_final_temp)
  
  # write out diagnostics
  params_to_output = c(
                       current_conditions, 
                       list(prefix_this=prefix_this,
                            minsize=minsize,
                         maxsize=maxsize,
                         n_iterations=n_iterations,
                         size_mean_initial=size_mean_initial,
                         size_sd_initial=size_sd_initial,
                         n_S_initial = n_S_initial,
                         n_A_initial = population_density_initial,
                         size_mean_final=size_mean_final,
                         n_A_final = population_density_final,
                         n_S_final=n_S_final,
                         longevity_90_final = longevity_90_final,
                         n_S_lag_factor=n_S_lag_factor,
                         n_S_scale_factor=n_S_scale_factor,
                         SWE.mean=mean(SWE.sequence),
                         STB.mean=mean(STB.sequence),
                         basal_area_density = basal_area_density,
                         n_A_final_gt_30_cm=n_A_final_gt_30_cm,
                         cv_n_A_final_ten=cv_n_A_final_ten)
  )
  # drop unnecessary columns for printing
  params_to_output$population_density_m2 = NULL
  params_to_output$SWE.0 = NULL
  params_to_output$SWE.1 = NULL
  params_to_output$STB.0 = NULL
  params_to_output$STB.1 = NULL
  
  # make plots
  g1 = ggplot(df_nt, aes(x=size,y=value,color=time,group=time)) +
    geom_line() +
    theme_bw() + 
    scale_x_log10() +
    xlab('DBH of tagged stems (cm)') +
    ylab('Population density (# m^-2)') +
    scale_color_viridis_c() +
    theme(legend.position = 'bottom')
  
  g2 = ggplot(df_n_density, aes(x=time,y=value,color=type)) +
    geom_line() +
    theme_bw() +
    ylab('Population density (# m^-2)') +
    scale_color_discrete(labels=c(`S`='S (medium stems)',`A`='A (adult stems)')) +
    xlab('Timestep') +
    facet_wrap(~type,scales='free') +
    theme(legend.position = 'bottom')

  # g3 = ggplot(n_S_component_time_series %>%
  #                 mutate(time=row_number()) %>%
  #                 pivot_longer(!time), aes(x=time,y=value,color=name)) +
  #   geom_line() +
  #   theme_bw() +
  #   ylab('Contribution to n_S')
  
  t_this = paste(names(params_to_output), "=", as.character(params_to_output),collapse = '\n',sep='')
  g4 = ggplot() +
    theme_void() +
    annotate('text', x=0, y=0, label=t_this,size=2)
  
  #g_pi = plot_kernel(minsize, maxsize, m=m, k=P_initial,'P initial')
  #g_pf = plot_kernel(minsize, maxsize, m=m, k=P_final,'P final')
  
  df_env = rbind(data.frame(name='SWE',value=SWE.sequence) %>% mutate(time=1:nrow(.)),data.frame(name='STB',value=STB.sequence) %>% mutate(time=1:nrow(.)))
  g_env = ggplot(df_env,
                 aes(x=time,y=value)) + 
    geom_line() +
    facet_wrap(~name,scales='free_y') +
    theme_bw()
  
  g_model = ggarrange(g1, g2, g4, g_env, nrow=2, ncol=2)
  if (!file.exists('output_figures/model_outputs'))
  {
    dir.create('output_figures/model_outputs')
  }
  ggsave(g_model, file=sprintf('output_figures/model_outputs/m_time_series_%s.pdf',prefix_this),width=12,height=8)
  

  return(list(P_initial=P_initial,
              P_final=P_final,
              size_bins=size_bins,
              n_S = n_S_time_series,
              n_A = n_A_time_series,
           results=params_to_output))
}

run_model_resampled <- function(prefix_this, cores=1, ...)
{
  mclapply(1:NUM_RESAMPLES, function(i)
  {
    print(i)
    run_model(
      prefix_this = paste(prefix_this, i, sep='_'),
      m_growth_this = m_growth_resampled[[i]],
      m_n_medium_count_this = m_n_medium_count_resampled[[i]],
      m_prob_medium_this = m_prob_medium_resampled[[i]],
      m_recruit_count_this = m_recruit_count_resampled[[i]],
      m_survival_this = m_survival_resampled[[i]],
      m_growth_size_variance_this = m_growth_size_variance_resampled[[i]],
      ...)
  }, mc.cores=cores)
}












# test run
model_output_test = run_model(geneticSexIDmale = '0',#geneticSexIDfemale = '0',
          Ploidy_leveltriploid = '0',
          Cos.aspect = 0,
          Elevation=3200,
          prefix_this = 'testfemale', 
          n_iterations = 100,
          n_S_scale_factor=1.5,
          SWE.sequence=rep(100,1000+1),
          STB.sequence=rep(0,1000+1))

 model_output_test_resampled = run_model_resampled(geneticSexIDmale = '0',#geneticSexIDfemale = '1',
                              Ploidy_leveltriploid = '0',
                              Cos.aspect = 0,
                              Elevation=3000,
                              prefix_this = 'test_resampled', 
                              n_iterations = 100,
                              SWE.sequence=rep(100,1000+1),
                              STB.sequence=rep(0,1000+1),
                              progress = TRUE,
                              cores=5)

sapply(model_output_test_resampled, function(x) { try(x$results$n_A_final) }) %>%
  as.data.frame() %>%
  rename(n_A_final =1) %>%
  ggplot(aes(x=n_A_final)) + 
  geom_histogram() +
  theme_bw()
# 
# z = do.call("rbind",lapply(1:length(model_output_test_resampled), function(i) {
#   
#   x=model_output_test_resampled[[i]]
#   data.frame(replicate=i, 
#              dbh=x$size_bins, 
#              density=as.numeric(tail(x$n_A,1)))
#   
#   }))
# 
# ggplot(z, aes(x=dbh,y=density,group=replicate)) +
#   geom_line() +
#   theme_bw()


model_output_test_stochastic = run_model(geneticSexIDmale = '1',#geneticSexIDfemale = '0',
                              Ploidy_leveltriploid = '1',
                              Cos.aspect = -1,
                              prefix_this = 'test_stochastic', 
                              n_iterations = 500,
                              n_S_lag_factor=0.2,
                              SWE.sequence=runif(500+1,0,300),
                              STB.sequence=runif(500+1,-1,1))

model_output_test_drought = run_model(geneticSexIDmale = '1',
                              Ploidy_leveltriploid = '0',
                              Cos.aspect = 0,
                              prefix_this = 'test_drought', 
                              n_iterations = 500,
                              SWE.sequence=c(rep(200,100),seq(200,0,length.out=50+1),rep(0,1000)),
                              STB.sequence=c(rep(2,100),seq(2,-2,length.out=50+1),rep(-2,1000)))


## test
# 
# matR_this = matrix(0, nrow=nrow(model_output_test$P_initial), ncol=ncol(model_output_test$P_initial))
# size_bins = seq(5,60,length.out=nrow(model_output_test$P_initial))
# matR_this[,head(which(size_bins>30),1):ncol(matR_this)] = 1
# mature_age(matU=model_output_test$P_initial, matR=matR_this,start=1) # output in # of timesteps
# 
# # compare against actual growth trajectories
# with(transitions_all_filtered_joined_no_na, (sizeNext - size)/delta_years) %>% summary
# # seems reasonable
# # envelope calculation of # of timesteps (years/3) required to reach 30 cm from 5 cm
# (30-5)/0.15/3





# run parameter sweep
params = expand.grid(geneticSexIDmale=c('0','1'), 
                     Ploidy_leveltriploid=c('0','1'), 
                     Cos.aspect=c(-1,1),
                     Elevation=c(2800,3200),
                     SWE=mean(transitions_all_filtered_joined_no_na$SWE.0),
                     STB=mean(transitions_all_filtered_joined_no_na$STB.0),
                     population_density_initial = mean(transitions_all_filtered_joined_no_na$population_density_m2),
                     n_iterations=300,
                     cores=5,
                     stringsAsFactors=FALSE
              )
print(nrow(params))
params$prefix_this=paste('param_sweep',1:nrow(params),sep='_')

model_results = lapply(1:nrow(params), function(i) {
  cat(sprintf('%d %.3f\n',i, i/nrow(params)))
  
  params_this = as.list(params[i,])
  # make environment sequence
  params_this$SWE.sequence=rep(params_this$SWE,params_this$n_iterations+1)
  params_this$STB.sequence=rep(params_this$STB,params_this$n_iterations+1)
  # remove unneeded columns
  params_this$SWE = NULL
  params_this$STB = NULL
  #print(params_this)
  
  result = do.call("run_model_resampled",params_this)
  return(result)
})
model_results_df = do.call('rbind',lapply(1:length(model_results), function(r) { do.call('rbind',lapply(model_results[[r]], function(x) { as.data.frame(x$results) %>% mutate(index=r) })) }))


model_results_df = model_results_df %>%
  rename(Ploidy_level=Ploidy_leveltriploid) %>%
  rename(geneticSexID=geneticSexIDmale) %>%
  mutate(Ploidy_level=ifelse(Ploidy_level=='1','triploid','diploid')) %>%
  mutate(geneticSexID=ifelse(geneticSexID=='1','male','female'))


write.csv(model_results_df, file='output_data/ipm_outcomes_param_sweep.csv', row.names = FALSE)



model_results_size_df = do.call('rbind',lapply(1:length(model_results), function(r) {
  model_results_size = do.call("rbind",lapply(1:length(model_results[[r]]), function(i) {
    
    x=model_results[[r]][[i]]
    data.frame(replicate=i, 
               dbh=x$size_bins, 
               density=as.numeric(tail(x$n_A,1)))
    
  }))
  
  model_results_size$index = r
  return(model_results_size)
    
}))


write.csv(model_results_size_df, file='output_data/ipm_outcomes_parameter_sweep_size_distribution.csv', row.names = FALSE)


# try some plotting
# slope changes
# ggplot(model_results_df, aes(x=coef_nt_stable)) +
#   geom_histogram(binwidth=0.0002) +
#   facet_grid(Ploidy_leveltriploid~geneticSexIDmale,labeller = label_both) +
#   theme_bw()


# ggplot(model_results_df, aes(x=n_A_final_gt_30_cm)) +
#   geom_histogram(binwidth=0.01) +
#   facet_grid(Ploidy_leveltriploid~geneticSexIDmale,labeller = label_both) +
#   theme_bw()

# look at size structure across runs
g_results_by_run = ggplot(model_results_size_df , aes(x=dbh,y=density, color=factor(index), group=paste(replicate, index))) +
  geom_line() +
  theme_bw() +
  scale_y_sqrt() +
  facet_wrap(~index) +
  theme(legend.position = 'none')
ggsave(g_results_by_run, file='output_figures/g_results_by_run.pdf',width=10,height=10)
ggsave(g_results_by_run, file='output_figures/g_results_by_run.png',width=10,height=10)






# variation in final density
g_param_sweep_n_A = ggplot(model_results_df, aes(x=geneticSexID, color=Ploidy_level,y=n_A_final)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitterdodge(),alpha=0.5) + 
  theme_bw() +
  facet_grid(Cos.aspect~Elevation,labeller = label_both) +
  ylab('Population density, N_A (# per m2)') +
  scale_color_manual(values=c('blue','red'),name='Ploidy level') +
  xlab('Sex') +
  ylim(0,1.5)

g_param_sweep_bad = ggplot(model_results_df, aes(x=geneticSexID, color=Ploidy_level,y=basal_area_density)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitterdodge(),alpha=0.5) + 
  theme_bw() +
  facet_grid(Cos.aspect~Elevation,labeller = label_both) +
  ylab('Basal area density, BAD_A (m2 / m2)') +
  scale_color_manual(values=c('blue','red'),name='Ploidy level') +
  xlab('Sex')

g_param_sweep_n_S = ggplot(model_results_df, aes(x=geneticSexID, color=Ploidy_level,y=n_S_final)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitterdodge(),alpha=0.5) + 
  theme_bw() +
  facet_grid(Cos.aspect~Elevation,labeller = label_both) +
  ylab('# medium stems, n_S') +
  scale_color_manual(values=c('blue','red'),name='Ploidy level') +
  xlab('Sex')


g_param_sweep_longevity = ggplot(model_results_df, aes(x=geneticSexID, color=Ploidy_level,y=2.5*longevity_90_final)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitterdodge(),alpha=0.5) + 
  theme_bw() +
  facet_grid(Cos.aspect~Elevation,labeller = label_both) +
  ylab('90% quantile longevity, L90 (years)') +
  scale_color_manual(values=c('blue','red'),name='Ploidy level') +
  xlab('Sex')


g_param_sweep_combined = ggarrange(g_param_sweep_n_A, g_param_sweep_bad, g_param_sweep_n_S, g_param_sweep_longevity, nrow=2, ncol=2, align='hv',common.legend = TRUE,legend='bottom',labels='auto')

ggsave(g_param_sweep_combined, file='output_figures/g_param_sweep_combined.png',width=7,height=7)
ggsave(g_param_sweep_combined, file='output_figures/g_param_sweep_combined.pdf',width=7,height=7)


g_param_sweep_pairs = ggpairs(model_results_df, mapping=aes(color=factor(Ploidy_level)), 
        columns=c('basal_area_density', 'n_A_final', 'n_S_final', 'longevity_90_final'),
        diag='blank') +
  theme_bw() +
  scale_color_manual(values=c('blue','red'),name='Ploidy level')
ggsave(g_param_sweep_pairs, file='output_figures/g_param_sweep_pairs.png',width=7,height=7)
ggsave(g_param_sweep_pairs, file='output_figures/g_param_sweep_pairs.pdf',width=7,height=7)




# compare observed and predicted outcomes
# not quite a fair fight because predicted is a parameter sweep, the observed is for plots...
df_n_comparison = rbind(rbind(data.frame(type='predicted',name='n_A',value=model_results_df$n_A_final), 
data.frame(type='observed',name='n_A',value=counts_recruit$population_density_m2)),
rbind(data.frame(type='predicted',name='n_S',value=model_results_df$n_S_final), 
                          data.frame(type='observed',name='n_S',value=counts_recruit$n_medium_trees)),
data.frame(type='predicted',name='basal_area_density',value=model_results_df$basal_area_density), 
data.frame(type='observed',name='basal_area_density',value=transitions_all_filtered_joined$basal_area_density_live))

g_param_sweep_obs_pred = ggplot(df_n_comparison, aes(x=value,color=type,fill=type)) + 
  geom_density(alpha=0.5) + 
  facet_wrap(~name,scales='free',ncol=2,
             labeller=as_labeller(c(n_S='# medium stems, n_S',
                                  n_A='Population density, N_A',
                                  basal_area_density='Basal area density, BAD_A'))) +
  theme_bw() +
  scale_color_brewer(palette='Set2') +
  scale_fill_brewer(palette='Set2') +
  xlab('Value') + 
  ylab('Density') +
  theme(legend.position='bottom')
ggsave(g_param_sweep_obs_pred, file='output_figures/g_param_sweep_obs_pred.pdf',width=6,height=6)
ggsave(g_param_sweep_obs_pred, file='output_figures/g_param_sweep_obs_pred.png',width=6,height=6)


# try other summary stats

m_param_sweep_summary_n_A_final = lm(n_A_final ~ Ploidy_level + geneticSexID + Cos.aspect + Elevation,
   data=model_results_df)
tab_model(m_param_sweep_summary_n_A_final, file='output_figures/table_param_sweep_m_summary_n_A_final.html')

m_param_sweep_summary_basal_area_density = lm(basal_area_density ~ Ploidy_level + geneticSexID + Cos.aspect + Elevation,
                                     data=model_results_df)
tab_model(m_param_sweep_summary_basal_area_density, file='output_figures/table_param_sweep_m_summary_basal_area_density.html')


m_param_sweep_summary_n_S_final = lm(n_S_final ~ Ploidy_level + geneticSexID + Cos.aspect + Elevation,
                                     data=model_results_df)
tab_model(m_param_sweep_summary_n_S_final, file='output_figures/table_param_sweep_m_summary_n_S_final.html')


m_param_sweep_summary_longevity_90_final = lm(longevity_90_final ~ Ploidy_level + geneticSexID + Cos.aspect + Elevation,
                                              data=model_results_df)
tab_model(m_param_sweep_summary_longevity_90_final, file='output_figures/table_param_sweep_m_summary_longevity_90_final.html')

















# load in climate data
df_climate = read.csv('output_data/df_climate.csv')
df_sites = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/manuscript ploidy 2019/new phyt/SI new phyt/File S1 - aspen data site-level processed 30 Mar 2020.csv') %>%
  select(site_code=Site_Code, Latitude, Longitude, X.UTM, Y.UTM, Basal.area.density)



# run across all sites
df_sites_for_ipm = transitions_all_filtered_joined_no_na %>%
  filter(site_type!='grid') %>%
  left_join(df_sites,by='site_code') %>%
  # rename a few things to standardize with the IPM code
  select(geneticSexID, Ploidy_level, Cos.aspect, Elevation,
         Latitude,Longitude,
         #STB.mean, STB.sd, SWE.meansqrt, SWE.sdsqrt,
         population_density_initial = population_density_m2, 
         site_code,
         size_mean, size_sd,
         #n_S_initial = n_medium_trees, # not currently used
         site_type,
         year) %>%
  group_by(site_code) %>%
  # get plot in most recent year
  slice_max(year) %>%
  unique %>%
  mutate(n_iterations=300, cores=5, progress=FALSE) %>%
  # rename columns for the script
  mutate(geneticSexIDmale=as.character(as.numeric((geneticSexID=='male'))),
         Ploidy_leveltriploid=as.character(as.numeric((Ploidy_level=='triploid'))),
         geneticSexIDunknown=as.character(as.numeric((geneticSexID=='unknown'))),
         Ploidy_levelunknown=as.character(as.numeric((Ploidy_level=='unknown')))) %>%
  select(-geneticSexID, -Ploidy_level) %>% 
  slice(rep(1:n(), each=3)) %>%
  group_by(site_code) %>%
  mutate(prefix_this=paste('rmbl', site_code, row_number(), sep='_'))

source('get_climate.R')

set.seed(1)
climate_all = lapply(1:nrow(df_sites_for_ipm), function(i)
{
  cat(sprintf('%d %.3f\n',i, i/nrow(df_sites_for_ipm)))
  lat_this = df_sites_for_ipm$Latitude[i]
  lon_this = df_sites_for_ipm$Longitude[i]
  climate_this = make_climate_ts_at_location(lon=lon_this, lat=lat_this, num_time_points = df_sites_for_ipm$n_iterations[i]+1)
  return(climate_this)
})

# try to map out all the sites
#1:nrow(df_sites_for_ipm) 
# code seems to get slower with more runs - maybe a virtual memory issue? consider writing each iteration to disk then re-assembling
ipms_rmbl_sites = lapply(1:nrow(df_sites_for_ipm), function(i)
{
  cat(sprintf('%d %.3f\n',i, i/nrow(df_sites_for_ipm)))
  params_this = as.list(df_sites_for_ipm[i,])
  # make environment sequence
  climate_this = climate_all[[i]]
  params_this$SWE.sequence=climate_this$SWE
  params_this$STB.sequence=climate_this$STB
  # remove unneeded columns
  params_this$year = NULL
  params_this$site_type = NULL
  params_this$site_code = NULL
  params_this$Latitude = NULL
  params_this$Longitude = NULL
  #print(params_this)
  
  result = do.call("run_model_resampled",params_this)
  
  return(result)
})



ipm_rmbl_df = do.call('rbind',lapply(1:length(ipms_rmbl_sites), function(r) { do.call('rbind',lapply(ipms_rmbl_sites[[r]], function(x) { as.data.frame(x$results) %>% mutate(index=r) })) })) %>%
  mutate(site_code=sapply(strsplit(prefix_this,split='_'),function(x) {x[2]})) %>%
  mutate(replicate_climate_resample=sapply(strsplit(prefix_this,split='_'),function(x) {x[3]})) %>%
  mutate(replicate_data_resample=sapply(strsplit(prefix_this,split='_'),function(x) {x[4]})) %>%
  left_join(df_sites,by='site_code') %>% 
  mutate(Ploidy_level=ifelse(Ploidy_leveltriploid=='1','triploid',ifelse(Ploidy_levelunknown=='1','unknown','diploid'))) %>%
  mutate(geneticSexID=ifelse(geneticSexIDmale=='1','male',ifelse(geneticSexIDunknown=='1','unknown','female')))

saveRDS(ipms_rmbl_sites, file='output_data/ipms_rmbl_sites.Rdata')

write.csv(ipm_rmbl_df, file='output_data/ipm_rmbl_df.csv',row.names = FALSE)




# do some analyses

m_varcomp_rmbl_n_A = lmer(n_A_final~ (1|site_code/replicate_data_resample),data=ipm_rmbl_df)
variance_fractions_rmbl_n_A = VarCorr(m_varcomp_rmbl_n_A,comp='Variance') %>%
  as.data.frame %>%
  select(grp, vcov) %>%
  mutate(fraction=vcov/sum(vcov)) %>%
  mutate(grp_nice = c(`replicate_data_resample:site_code`='Among dataset resamples',
                      site_code='Among sites',
                      Residual='Among climate resamples')[grp])
variance_fractions_rmbl_n_A

m_varcomp_rmbl_bad = lmer(basal_area_density~ (1|site_code/replicate_data_resample),data=ipm_rmbl_df)
variance_fractions_rmbl_bad = VarCorr(m_varcomp_rmbl_bad,comp='Variance') %>%
  as.data.frame %>%
  select(grp, vcov) %>%
  mutate(fraction=vcov/sum(vcov)) %>%
  mutate(grp_nice = c(`replicate_data_resample:site_code`='Among dataset resamples',
                      site_code='Among sites',
                      Residual='Among climate resamples')[grp])
variance_fractions_rmbl_bad

m_varcomp_rmbl_n_S = lmer(n_S_final ~ (1|site_code/replicate_data_resample),data=ipm_rmbl_df)
variance_fractions_rmbl_n_S = VarCorr(m_varcomp_rmbl_n_S,comp='Variance') %>%
  as.data.frame %>%
  select(grp, vcov) %>%
  mutate(fraction=vcov/sum(vcov)) %>%
  mutate(grp_nice = c(`replicate_data_resample:site_code`='Among dataset resamples',
                      site_code='Among sites',
                      Residual='Among climate resamples')[grp])
variance_fractions_rmbl_n_S

m_varcomp_rmbl_L_90 = lmer(longevity_90_final ~ (1|site_code/replicate_data_resample),data=ipm_rmbl_df)
variance_fractions_rmbl_L_90 = VarCorr(m_varcomp_rmbl_L_90,comp='Variance') %>%
  as.data.frame %>%
  select(grp, vcov) %>%
  mutate(fraction=vcov/sum(vcov)) %>%
  mutate(grp_nice = c(`replicate_data_resample:site_code`='Among dataset resamples',
                      site_code='Among sites',
                      Residual='Among climate resamples')[grp])
variance_fractions_rmbl_L_90


g_varcomp_rmbl_n_A = ggplot(variance_fractions_rmbl_n_A, aes(x='',fill=grp_nice,y=fraction)) + 
  geom_col() +
  coord_polar(theta='y') +
  theme_minimal() +
  scale_fill_brewer(palette='Set2',name='Scale') +
  ylab('Population density N_A') + xlab('')
g_varcomp_rmbl_bad = ggplot(variance_fractions_rmbl_bad, aes(x='',fill=grp_nice,y=fraction)) + 
  geom_col() +
  coord_polar(theta='y') +
  theme_minimal() +
  scale_fill_brewer(palette='Set2',name='Scale') +
  ylab('Basal area density BAD') + xlab('')
g_varcomp_rmbl_n_S = ggplot(variance_fractions_rmbl_n_S, aes(x='',fill=grp_nice,y=fraction)) + 
  geom_col() +
  coord_polar(theta='y') +
  theme_minimal() +
  scale_fill_brewer(palette='Set2',name='Scale') +
  ylab('Number medium stems') + xlab('')
g_varcomp_rmbl_L_90 = ggplot(variance_fractions_rmbl_L_90, aes(x='',fill=grp_nice,y=fraction)) + 
  geom_col() +
  coord_polar(theta='y') +
  theme_minimal() +
  scale_fill_brewer(palette='Set2',name='Scale') +
  ylab('90% longevity') + xlab('')

g_varcomp_rmbl = ggarrange(g_varcomp_rmbl_n_A, g_varcomp_rmbl_bad,
                           g_varcomp_rmbl_n_S, g_varcomp_rmbl_L_90,
                           labels='AUTO',align='hv',common.legend = TRUE,legend='bottom')

ggsave(g_varcomp_rmbl, file='output_figures/g_varcomp_rmbl.png',width=7,height=7)
ggsave(g_varcomp_rmbl, file='output_figures/g_varcomp_rmbl.pdf',width=7,height=7)
# 
# ggplot(ipm_rmbl_df, aes(x=Elevation,y=n_A_final,color=site_code)) +
#   geom_violin()





# 
# # do histograms
# g_rmbl_delta_n_A = ggplot(ipm_rmbl_df %>% 
#                             filter(Ploidy_level!='unknown') %>% 
#                             filter(geneticSexID!='unknown') %>%
#                             mutate(delta_n_A = n_A_final - n_A_initial), 
#        aes(x=geneticSexID,color=Ploidy_level,y=delta_n_A)) +
#   geom_violin() +
#   theme_bw() + 
#   xlab('Sex') +
#   ylab('Delta population density N_A') +
#   scale_color_manual(values=c('blue','red'),name='Ploidy level')
# 
# 
# g_rmbl_bad = ggplot(ipm_rmbl_df %>% filter(Ploidy_level!='unknown') %>% filter(geneticSexID!='unknown'), 
#                     aes(x=geneticSexID,color=Ploidy_level,y=basal_area_density)) +
#   geom_boxplot() +
#   theme_bw() + 
#   xlab('Sex') +
#   ylab('Basal area density BAD') +
#   scale_color_manual(values=c('blue','red'),name='Ploidy level')
# g_rmbl_bad = ggplot(ipm_rmbl_df %>% filter(Ploidy_level!='unknown') %>% filter(geneticSexID!='unknown'), 
#                     aes(x=geneticSexID,color=Ploidy_level,y=basal_area_density)) +
#   geom_boxplot() +
#   theme_bw() + 
#   xlab('Sex') +
#   ylab('Basal area density BAD') +
#   scale_color_manual(values=c('blue','red'),name='Ploidy level')
#   
# # g_rmbl_hist = ggplot(ipm_rmbl_df, aes(x=n_A_final)) + 
# #   geom_histogram(binwidth = 0.01) + 
# #   scale_x_sqrt() +
# #   theme_bw() +
# #   facet_grid(geneticSexID~Ploidy_level,scales='free_y')
# # ggsave(g_rmbl_hist, file='output_figures/g_rmbl_hist.pdf')
# # ggsave(g_rmbl_hist, file='output_figures/g_rmbl_hist.png')
# 
# 
# g_distribution_rmbl = ggplot(ipm_rmbl_df, aes(x=n_A_final,y=n_S_final,group=site_code)) + 
#   #geom_convexhull(alpha=0.1,color='gray',fill='lightgray') + 
#   geom_point(size=0.5,alpha=0.8) + 
#   scale_x_sqrt() +
#   scale_y_sqrt() +
#   theme_bw()
# ggsave(g_distribution_rmbl, file='output_figures/g_distribution_rmbl.pdf',width=8,height=8)
# ggsave(g_distribution_rmbl, file='output_figures/g_distribution_rmbl.png',width=8,height=8)
# 










#write.csv(age_structure_sites_all, file='output_data/sites_age_structure.csv', row.names = FALSE)



# %>%
#mutate(lambda_binned = cut(lambda, breaks=c(0,0.99,1.01,Inf),labels=c('decreasing','stable','increasing')))

# df_sites_for_ipm_joined$dbh_modal = dbh_range[apply(age_structure_sites_all, 1, function(x) { 
#   result = which.max(x)
#   if (length(result)==0)
#   {
#     result = NA
#   }
#   return(result)
#   })]
# 
# write.csv(df_sites_for_ipm_joined, file='output_data/ipm_outcomes_lambda_by_plot.csv', row.names = FALSE)
# 
# df_sites_for_ipm_joined = read.csv('output_data/ipm_outcomes_lambda_by_plot.csv')
# 

# add in initial basal area densities
df_bad_initial = transitions_all_filtered_joined %>% 
  select(site_code, year, basal_area_density_live) %>% 
  group_by(site_code) %>%
  summarize(basal_area_density_initial = mean(basal_area_density_live,na.rm=TRUE))

ipm_rmbl_df_mean = ipm_rmbl_df %>% 
  left_join(df_bad_initial, by='site_code') %>%
  group_by(site_code, X.UTM, Y.UTM, Ploidy_level, geneticSexID) %>% 
  summarize(n_A.mean=mean(n_A_final,na.rm=TRUE), 
            n_A.sd=sd(n_A_final,na.rm=TRUE),
            n_S.mean=mean(n_S_final,na.rm=TRUE), 
            n_S.sd=sd(n_S_final,na.rm=TRUE),
            bad.mean=mean(basal_area_density, na.rm=TRUE),
            ba.sd=sd(basal_area_density,na.rm=TRUE),
            delta_n_A.mean = mean(n_A_final - n_A_initial,na.rm=TRUE),
            delta_bad.mean = mean(basal_area_density - basal_area_density_initial, na.rm=TRUE),
            delta_n_S.mean = mean(n_S_final - n_S_initial, na.rm=TRUE))


(ipm_rmbl_df_mean$delta_n_A.mean > 0) %>% table
(ipm_rmbl_df_mean$delta_bad.mean > 0) %>% table
(ipm_rmbl_df_mean$delta_n_S.mean > 0) %>% table
# mean(df_sites_for_ipm_joined_summarized$lambda.mean,na.rm=TRUE)
# sd(df_sites_for_ipm_joined_summarized$lambda.mean,na.rm=TRUE)
# 
# table(df_sites_for_ipm_joined_summarized$lambda_binned)
# table(df_sites_for_ipm_joined_summarized$lambda_binned)/nrow(df_sites_for_ipm_joined_summarized)

g_rmbl_delta_n_A = ggplot(ipm_rmbl_df_mean %>% 
                            filter(Ploidy_level!='unknown') %>% 
                            filter(geneticSexID!='unknown'), 
                          aes(x=delta_n_A.mean)) +
  facet_grid(Ploidy_level~geneticSexID) +
  geom_density(fill='gray',alpha=0.5,color='gray20') +
  geom_vline(xintercept = 0) +
  #geom_histogram(binwidth=0.02) +
  theme_bw() + 
  xlab('Delta population density N_A') +
  ylab('Probability density')# +
  #scale_color_manual(values=c('blue','red'),name='Ploidy level')

g_rmbl_delta_bad = ggplot(ipm_rmbl_df_mean %>% 
                            filter(Ploidy_level!='unknown') %>% 
                            filter(geneticSexID!='unknown'), 
                          aes(x=delta_bad.mean)) +
  facet_grid(Ploidy_level~geneticSexID) +
  geom_vline(xintercept = 0) +
  geom_density(fill='gray',alpha=0.5) +
  #geom_histogram(binwidth=0.02) +
  theme_bw() + 
  xlab('Delta basal area density BAD') +
  ylab('Probability density')

g_rmbl_delta_n_S = ggplot(ipm_rmbl_df_mean %>% 
                            filter(Ploidy_level!='unknown') %>% 
                            filter(geneticSexID!='unknown'), 
                          aes(x=delta_n_S.mean)) +
  facet_grid(Ploidy_level~geneticSexID) +
  geom_vline(xintercept = 0) +
  geom_density(fill='gray',alpha=0.5) +
  #geom_histogram(binwidth=0.02) +
  theme_bw() + 
  xlab('Delta number medium n_S') +
  ylab('Probability density')



coords = ipm_rmbl_df_mean %>% st_as_sf(coords=c("X.UTM","Y.UTM"), crs = 32613)
elev = get_elev_raster(st_bbox(coords), z=12)
elev_cropped = crop(elev, st_bbox(coords)+c(-500,-500,500,500))
sl <- terrain(rast(elev_cropped), "slope", unit = "radians")
asp <- terrain(rast(elev_cropped), "aspect", unit = "radians")
hill_single <- shade(sl, asp,
                     angle = 45,
                     direction = 300,
                     normalize = TRUE
)



plot_map <- function(data_this, yvar, scale_color_this)
{
  g = ggplot() + 
    theme_void() + 
    geom_stars(data=st_as_stars(hill_single),show.legend=FALSE,alpha=0.35) +
    scale_fill_distiller(palette = "Greys") +
    new_scale_fill() +
    geom_point(data=data_this,
               aes(x=X.UTM,y=Y.UTM, color=.data[[yvar]]), alpha=0.8) +
    geom_point(data=ipm_rmbl_df_mean,
               aes(x=X.UTM,y=Y.UTM), alpha=0.25,color='black',shape=1, inherit.aes = FALSE) +
    scale_color_this +
    
    coord_equal() +
    #xlab('Easting (m)') + ylab('Northing (m)') +
    theme(legend.position='bottom')
  return(g)
}
# g_map_rmbl_n_A = plot_map(ipm_rmbl_df_mean, 'n_A.mean',scale_color_viridis_c(name='Tagged tree density (# m-2)'))
# g_map_rmbl_n_S = plot_map(ipm_rmbl_df_mean, 'n_S.mean',scale_color_viridis_c(name='# medium'))
# g_map_rmbl_bad = plot_map(ipm_rmbl_df_mean, 'bad.mean',scale_color_viridis_c('Basal area density'))
g_map_rmbl_delta_n_A = plot_map(ipm_rmbl_df_mean, 'delta_n_A.mean',scale_color_gradient2(name='Delta n_A',midpoint=0,low='darkorange1',high='olivedrab',mid='ivory2'))
g_map_rmbl_delta_bad = plot_map(ipm_rmbl_df_mean, 'delta_bad.mean',scale_color_gradient2(name='Delta basal area density',midpoint=0,low='darkorange1',high='olivedrab',mid='ivory2'))
g_map_rmbl_delta_n_S = plot_map(ipm_rmbl_df_mean, 'delta_n_S.mean',scale_color_gradient2(name='Delta number medium n_S',midpoint=0,low='darkorange1',high='olivedrab',mid='ivory2'))


g_rmbl_maps = ggarrange(ggarrange(g_rmbl_delta_n_A, g_rmbl_delta_bad, g_rmbl_delta_n_S, nrow=1, labels='AUTO'), 
          ggarrange(g_map_rmbl_delta_n_A, g_map_rmbl_delta_bad, g_map_rmbl_delta_n_S, nrow=1,align='hv', labels=c('D','E','F')),
          nrow=2,align='hv')
ggsave(g_rmbl_maps, file='output_figures/g_rmbl_maps.png',width=10,height=7)
ggsave(g_rmbl_maps, file='output_figures/g_rmbl_maps.pdf',width=10,height=7)


# ggsave(g_map_rmbl_n_A, file='output_figures/g_map_rmbl_n_A.pdf',width=8,height=5)
# ggsave(g_map_rmbl_n_S, file='output_figures/g_map_rmbl_n_S.pdf',width=8,height=5)
# ggsave(g_map_rmbl_bad, file='output_figures/g_map_rmbl_bad.pdf',width=8,height=5)
# ggsave(g_map_rmbl_delta_n_A, file='output_figures/g_map_rmbl_delta_n_A.pdf',width=8,height=5)
# ggsave(g_map_rmbl_delta_bad, file='output_figures/g_map_rmbl_delta_bad.pdf',width=8,height=5)
# ggsave(g_map_rmbl_n_A, file='output_figures/g_map_rmbl_n_A.png',width=8,height=5)
# ggsave(g_map_rmbl_n_S, file='output_figures/g_map_rmbl_n_S.png',width=8,height=5)
# ggsave(g_map_rmbl_bad, file='output_figures/g_map_rmbl_bad.png',width=8,height=5)
# ggsave(g_map_rmbl_delta_n_A, file='output_figures/g_map_rmbl_delta_n_A.png',width=8,height=5)
# ggsave(g_map_rmbl_delta_bad, file='output_figures/g_map_rmbl_delta_bad.png',width=8,height=5)



df_n_comparison_rmbl = rbind(rbind(data.frame(type='predicted',name='n_A',value=ipm_rmbl_df_mean$n_A.mean), 
                              data.frame(type='observed',name='n_A',value=counts_recruit$population_density_m2)),
                        rbind(data.frame(type='predicted',name='n_S',value=ipm_rmbl_df_mean$n_S.mean), 
                              data.frame(type='observed',name='n_S',value=counts_recruit$n_medium_trees)),
                        data.frame(type='predicted',name='basal_area_density',value=ipm_rmbl_df_mean$bad.mean), 
                        data.frame(type='observed',name='basal_area_density',value=transitions_all_filtered_joined$basal_area_density_live))

g_param_sweep_obs_pred_rmbl = ggplot(df_n_comparison_rmbl, aes(x=value,color=type)) + 
  geom_density() +
  facet_wrap(~name,scales='free',
             labeller=as_labeller(c(n_S='# medium trees',
                                    n_A='# tagged trees',
                                    basal_area_density='Basal area density'))) +
  theme_bw() +
  scale_color_brewer(palette='Set2') +
  xlab('Value') + 
  ylab('Density')
ggsave(g_param_sweep_obs_pred_rmbl, file='output_figures/g_param_sweep_obs_pred_rmbl.pdf',width=8,height=4)
ggsave(g_param_sweep_obs_pred_rmbl, file='output_figures/g_param_sweep_obs_pred_rmbl.png',width=8,height=4)







# 
# 
# # make new params with n_S_star_factor
# 
# #ipms_rmbl_sites = lapply(1:nrow(df_sites_for_ipm), function(i)
# {
#   cat(sprintf('%d %.3f\n',i, i/nrow(df_sites_for_ipm)))
#   params_this = as.list(df_sites_for_ipm[i,])
#   # make environment sequence
#   climate_this = climate_all[[i]]
#   params_this$SWE.sequence=climate_this$SWE
#   params_this$STB.sequence=climate_this$STB
#   # remove unneeded columns
#   params_this$year = NULL
#   params_this$site_type = NULL
#   params_this$site_code = NULL
#   params_this$Latitude = NULL
#   params_this$Longitude = NULL
#   #print(params_this)
#   
#   result = do.call("run_model_resampled",params_this)
#   
#   return(result)
# # })
# # 
# # ipms_rmbl_sites_dry_scenario = lapply(1:nrow(df_sites_for_ipm), function(i)
# # {
# #   cat(sprintf('%d %.3f\n',i, i/nrow(df_sites_for_ipm)))
# #   params_this = as.list(df_sites_for_ipm[i,])
# #   # make environment sequence
# #   climate_this = climate_all[[i]]
# #   # make the climate more stressful
# #   params_this$SWE.sequence=climate_this$SWE
# #   params_this$STB.sequence=climate_this$STB
# #   # remove unneeded columns
# #   params_this$year = NULL
# #   params_this$site_type = NULL
# #   params_this$site_code = NULL
# #   params_this$Latitude = NULL
# #   params_this$Longitude = NULL
# #   # change run name
# #   params_this$prefix_this = gsub('rmbl','rmbl_dry_scenario', params_this$prefix_this)
# #   #print(params_this)
# #   
# #   result = do.call("run_model_resampled",params_this)
# #   
# #   return(result)
# # })
# 
# 
# 
# ipm_rmbl_dry_scenario_df = do.call('rbind',lapply(1:length(ipms_rmbl_sites_dry_scenario), function(r) { do.call('rbind',lapply(ipms_rmbl_sites_dry_scenario[[r]], function(x) { as.data.frame(x$results) %>% mutate(index=r) })) })) %>%
#   mutate(site_code=sapply(strsplit(prefix_this,split='_'),function(x) {x[4]})) %>%
#   mutate(replicate_climate_resample=sapply(strsplit(prefix_this,split='_'),function(x) {x[5]})) %>%
#   mutate(replicate_data_resample=sapply(strsplit(prefix_this,split='_'),function(x) {x[6]})) %>%
#   left_join(df_sites,by='site_code') %>% 
#   mutate(Ploidy_level=ifelse(Ploidy_leveltriploid=='1','triploid',ifelse(Ploidy_levelunknown=='1','unknown','diploid'))) %>%
#   mutate(geneticSexID=ifelse(geneticSexIDmale=='1','male',ifelse(geneticSexIDunknown=='1','unknown','female')))
# 
# ipm_rmbl_dry_scenario_df_mean = ipm_rmbl_dry_scenario_df %>% 
#   left_join(df_bad_initial, by='site_code') %>%
#   group_by(site_code, X.UTM, Y.UTM) %>% 
#   summarize(n_A.mean=mean(n_A_final,na.rm=TRUE), 
#             n_A.sd=sd(n_A_final,na.rm=TRUE),
#             n_S.mean=mean(n_S_final,na.rm=TRUE), 
#             n_S.sd=sd(n_S_final,na.rm=TRUE),
#             bad.mean=mean(basal_area_density, na.rm=TRUE),
#             ba.sd=sd(basal_area_density,na.rm=TRUE),
#             delta_n_A.mean = mean(n_A_final - n_A_initial,na.rm=TRUE),
#             delta_bad.mean = mean(basal_area_density - basal_area_density_initial, na.rm=TRUE))
# 
# 
# 
# write.csv(ipm_rmbl_dry_scenario_df, file='output_data/ipm_rmbl_dry_scenario_df.csv',row.names = FALSE)
# 
# # analyze dry scenario
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



# save info for ipm
varnames_to_save = c('update_coefficients_full',
                     'run_model_resampled',
                     'run_model',
                     'plot_kernel',
                     'kernel_for_plotting',
                     'NUM_RESAMPLES',
                     'transitions_resampled',
                     'transitions_all_filtered_joined_no_na',
                     'coords', 'elev', 'elev_cropped', 'sl', 'asp', 'hill_single',
                     ls(pattern='^m_')) %>%
  sort
save(list=varnames_to_save,file='output_data/workspace for ipm.Rdata')













# 
# 
# 
# 
# 
# # OLD BELOW
# 
# 
# 
# # initial population structure
# hist(transitions_all_filtered_joined_no_na$population_density_m2,breaks=mesh_points)
# hist(transitions_all_filtered_joined_no_na$size)
# # make initial size distribution to be the grand mean/sd of the training data
# 
# # set the recruitment conditions
# make_ipm_for_site <- function(m_survival,
#                               m_growth,
#                               m_recruit,
#                               other_vars, 
#                               size_mean = mean(transitions_all_filtered_joined_no_na$size, na.rm=TRUE),
#                               size_sd = sd(transitions_all_filtered_joined_no_na$size, na.rm=TRUE),
#                               population_density_initial = mean(transitions_all_filtered_joined_no_na$population_density_m2))
# {
#   n_initial = sapply(1:mesh_points, function(i) {
#     size = dnorm(dbh_range[i], mean=size_mean, sd=size_sd)
#     size = ifelse(dbh_range[i] < DBH_min, 0, size)
#     return(size)
#   })
#   # normalize the size to a population density of the requested value
#   n_initial = n_initial * population_density_initial / sum(n_initial)
#   
#   plot(dbh_range, n_initial,type='h'); title(sprintf('population density=%.3f m-2',sum(n_initial)))
#   
#   
#   
#   # get size dependent coefficients, by setting population density to zero to eliminate its effect
#   # also standardize to a 3-year interval to interpret the transitions
#   coef_survival_size = update_coefficients_full(coef(m_survival), xvar='size', other_vars=c(other_vars,population_density_m2=0, delta_years=3))
#   coef_growth_size = update_coefficients_full(coef(m_growth), xvar='size', other_vars=c(other_vars,population_density_m2=0, delta_years=3))
#   coef_recruit_n_medium_trees = update_coefficients_full(coef(m_recruit), xvar='size', other_vars=c(other_vars,population_density_m2=0, delta_years=3))
#   # get density dependent coefficients, assming there are no size interactions
#   coef_growth_density = coef(m_growth)["population_density_m2"]
#   coef_survival_density = coef(m_survival)["population_density_m2"]
#   coef_recruit_density = coef(m_recruit)["population_density_m2"]
#   
#   
#   my_data_list = list(s_int            = coef_survival_size[1],
#                       s_slope_size     = coef_survival_size[2],
#                       s_slope_density  = coef_survival_density[1],
#                       g_int            = coef_growth_size[1],
#                       g_slope_size     = coef_growth_size[2],
#                       g_slope_density  = coef_growth_density[1],
#                       sd_g             = sd(resid(m_growth)),
#                       r_int            = coef_recruit_n_medium_trees[1],
#                       r_slope_density  = coef_recruit_density[1],
#                       r_d_mu           = DBH_min, # assume all come in at exactly size 5
#                       r_d_sd           = 0.5) # assume there is a little bit of fluctuation in sizes (needed so discretization works too) 
#   
#   ipm_aspen <- init_ipm(sim_gen = "simple",
#                             di_dd   = "dd",
#                             det_stoch = "det")
#   
#   ipm_aspen <- define_kernel(
#     
#     proto_ipm = ipm_aspen,
#     
#     # Name of the kernel
#     
#     name      = "P",
#     
#     # The type of transition it describes (e.g. continuous - continuous, discrete - continuous).
#     # These must be specified for all kernels!
#     
#     family    = "CC",
#     
#     # The formula for the kernel. We dont need to tack on the "z'/z"s here.  
#     
#     formula   = s * G,
#     
#     # A named set of expressions for the vital rates it includes. 
#     # note the use of user-specified functions here. Additionally, each 
#     # state variable has a stateVariable_1 and stateVariable_2, corresponding to
#     # z and z' in the equations above. We don't need to define these variables ourselves,
#     # just reference them correctly based on the way we've set up our model on paper.
#     
#     # Perform the inverse logit transformation to get survival probabilities
#     # from your model. plogis from the "stats" package does this for us. 
#     
#     s         = plogis(s_int + s_slope_size * dbh_1 + s_slope_density * sum(n_dbh_t)), 
#     
#     # The growth model requires a function to compute the mean as a function of dbh.
#     # The SD is a constant, so we don't need to define that in ... expression, 
#     # just the data_list.
#     
#     G         = dnorm(dbh_2, mu_g, sd_g), # dbh_2 corresponds to size_next
#     mu_g      = g_int + g_slope_size * dbh_1 + g_slope_density * sum(n_dbh_t),
#     
#     
#     # Specify the constant parameters in the model in the data_list. 
#     
#     data_list = my_data_list,
#     states    = list(c('dbh')),
#     
#     # If you want to correct for eviction, set evict_cor = TRUE and specify an
#     # evict_fun. ipmr provides truncated_distributions() to help. This function
#     # takes 2 arguments - the type of distribution, and the name of the parameter/
#     # vital rate that it acts on.
#     
#     evict_cor = TRUE,
#     evict_fun = truncated_distributions(fun    = 'norm',
#                                         target = 'G')
#   ) 
#   
#   ipm_aspen <- define_kernel(
#     proto_ipm = ipm_aspen,
#     name      = 'F',
#     formula   = r_r * r_d,
#     family    = 'CC',
#     
#     r_r       = plogis(r_int + r_slope_density * sum(n_dbh_t)), # overall rate of recruitment
#     r_d       = dnorm(dbh_2, r_d_mu, r_d_sd), # multiply this rate by the appropriate size factor (using size_2, the new size)
#     data_list = my_data_list,
#     states    = list(c('dbh')),
#     
#     # Again, we'll correct for eviction in new recruits by
#     # truncating the normal distribution.
#     
#     evict_cor = TRUE,
#     evict_fun = truncated_distributions(fun    = 'norm',
#                                         target = 'r_d')
#   ) 
#   
#   # Next, we have to define the implementation details for the model. 
#   # We need to tell ipmr how each kernel is integrated, what state
#   # it starts on (i.e. z from above), and what state
#   # it ends on (i.e. z' above). In simple_* models, state_start and state_end will 
#   # always be the same, because we only have a single continuous state variable. 
#   # General_* models will be more complicated.
#   
#   ipm_aspen <- define_impl(
#     proto_ipm = ipm_aspen,
#     make_impl_args_list(
#       kernel_names = c("P", "F"),
#       int_rule     = rep("midpoint", 2),
#       state_start  = rep("dbh", 2),
#       state_end    = rep("dbh", 2)
#     )
#   ) 
#   
#   ipm_aspen <- define_domains(
#     proto_ipm = ipm_aspen,
#     dbh = c(0, # the first entry is the lower bound of the domain.
#             DBH_max, # the second entry is the upper bound of the domain.
#             mesh_points # third entry is the number of meshpoints for the domain.
#     ) 
#   ) 
#   
#   # Next, we define the initial state of the population. We must do this because
#   # ipmr computes everything through simulation, and simulations require a 
#   # population state.
#   
#   ipm_aspen <- define_pop_state(
#     proto_ipm = ipm_aspen,
#     n_dbh = n_initial
#   )
#   
#   ipm_aspen <- make_ipm(proto_ipm = ipm_aspen, 
#                         iterations=MAX_ITERATIONS, 
#                         normalize_pop_size = FALSE,
#                         return_sub_kernels = TRUE)
#   
#   return(ipm_aspen)
# }
# 
# ipm_female_triploid_south = make_ipm_for_site(
#   m_survival = m_survival_all,
#   m_growth = m_growth_all,
#   m_recruit = m_recruit_all,
#   other_vars=list(geneticSexIDM="0",
#     Ploidy_levelTriploid="1",
#     n_medium_trees=5,
#     Cos.aspect=-1,
#     Elevation=3000)
#   )
# 
# # this is a density dependent deterministic model, so let's look
# # for elasticities close to equilibrium...
# K_female_triploid_south = as.matrix(ipm_female_triploid_south$sub_kernels$P_it_300 + ipm_female_triploid_south$sub_kernels$F_it_300)
# image(sqrt(elasticity(K_female_triploid_south)))
# e_female_triploid_south = elasticity(K_female_triploid_south)
# # # code below from https://compadre-db.org/Education/article/sensitivity-and-elasticity-matrices
# # w <- eigen(mat)$vectors
# # v <- Conj(solve(w))
# # senmat <- Re(v[1,] %*% t(w[,1]))
# # emat <- (1/(Re(eigen(mat)$values[1]))) * senmat * mat
# 
# 
# 
# # get lambda
# plot(ipmr::lambda(ipm_female_triploid_south),xlab='time step',ylab='lambda')
# 
# # look at trajectory of age distribution
# image(ipm_female_triploid_south$pop_state$n_dbh %>% sqrt,xlab='size',ylab='time')
# 
# # look at final age structure (sqrt transformed)
# plot(dbh_range, 
#      sqrt(ipm_male_triploid_south$pop_state$n_dbh[,MAX_ITERATIONS+1]),
#      type='h',xlab='DBH (cm)',ylab='sqrt # m-2')
# 
# 
# # # look at kernels at initial densities
# # pdf(file='output_figures/g_ipm_kernels.pdf',width=10,height=5)
# # par(mfrow=c(1,2))
# # plot(ipm_male_triploid_south$sub_kernels$P_it_1); title('P subkernel')
# # plot(ipm_all$sub_kernels$F_it_1); title('F subkernel')
# # dev.off()
# 
# k1 = plot_kernel(ipm_female_triploid_south$sub_kernels$P_it_200, 'P sub-kernel')
# k2 = plot_kernel(ipm_female_triploid_south$sub_kernels$F_it_200, 'F sub-kernel')
# k3 = plot_kernel((e_female_triploid_south)^(1/2), 'Elasticity (square root)')
# ggsave(ggarrange(k1, k2, k3, labels='AUTO', align='hv',nrow=2,ncol=2),file='output_figures/g_ipm_kernel_female_triploid_south.pdf',width=7,height=5)
# ggsave(ggarrange(k1, k2, k3, labels='AUTO', align='hv',nrow=2,ncol=2),file='output_figures/g_ipm_kernel_female_triploid_south.png',width=7,height=5)
# 
# 
# 
# 
# # now look at performance variation across parameters
# 
# params = expand.grid(replicate=1:NUM_RESAMPLES,
#                      geneticSexIDM=c("0","1"), 
#                      Ploidy_levelTriploid=c("0","1"), 
#                      n_medium_trees=c(0,5,15),
#                      Cos.aspect=c(-1,0,1),
#                      Elevation=c(2800,3000,3200))
# 
# lambdas = pbsapply(1:nrow(params), function(i) {
#   
#   cat(sprintf('%d %.3f\n',i, i/nrow(params)))
#   
#   ipm_this = NULL
#   try(ipm_this <- make_ipm_for_site(
#       m_survival = m_survival_resampled[[ params$replicate[i] ]],
#       m_growth = m_growth_resampled[[ params$replicate[i] ]],
#       m_recruit = m_recruit_resampled[[ params$replicate[i] ]],
#       other_vars = as.list(params[i,])
#     ))
#   
#   if (!is.null(ipm_this))
#   {
#     lambda = mean(tail(as.numeric(ipmr::lambda(ipm_this)),5)) # take the final 5 lambda values and average
#   }
#   else
#   {
#     lambda = NA
#   }
#   print(lambda)
#  
#   # get the mean lambda over time
#   return(lambda)
# }, cl = NUM_CORES)
# params$lambda = lambdas
# 
# write.csv(params, file='output_data/ipm_outcomes_lambda.csv', row.names = FALSE)
# 
# 
# #
# params = read.csv('output_data/ipm_outcomes_lambda.csv')
# 
# # summarize this model
# # library(glmmTMB)
# m_lambda_summary = glmmTMB(lambda ~ geneticSexIDM*Ploidy_levelTriploid + n_medium_trees + Cos.aspect + Elevation + (1|replicate), data=params)
# tab_model(m_lambda_summary, file='output_figures/table_model_lambda.html')
# 
# 
# params_for_plotting = params %>%
#   filter(Cos.aspect==-1 & Elevation==3000) %>%
#   mutate(sex=ifelse(geneticSexIDM==1,'M','F')) %>%
#   mutate(cytotype=ifelse(Ploidy_levelTriploid==1,'triploid','diploid')) %>%
#   mutate(n_medium_trees=factor(paste('# medium=',n_medium_trees,sep=''),levels=paste('# medium=',c(0,5,15),sep=''),ordered=TRUE)) %>%
#   mutate(Cos.aspect=paste('Cosine aspect=',Cos.aspect,sep=''))
# 
# 
# g_lambda =  ggplot(params_for_plotting, aes(color=factor(cytotype),x=factor(sex),y=lambda)) + 
#   geom_hline(yintercept = 1) +
#   geom_boxplot(outliers = FALSE) +
#   geom_point(alpha=0.5, position=position_jitterdodge()) +
#   # geom_point(alpha=0.5) + 
#   # geom_jitter(height=0,width=0.2) + 
#   #geom_bar(stat='identity',alpha=0.5,position='dodge') +
#   facet_grid(~n_medium_trees) +
#   theme_bw() +
#   #scale_y_log10() +
#   scale_color_manual(values=c('blue','red'),name='Cytotype') +
#   xlab('Sex') +
#   ylab(expression(paste(lambda)))
# 
# ggsave(g_lambda, file='output_figures/g_ipm_lambda.pdf',width=6,height=3)
# ggsave(g_lambda, file='output_figures/g_ipm_lambda.png',width=6,height=3)
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
# # now map lambdas for all the plots
# df_sites_for_ipm = transitions_all_filtered_joined_no_na %>%
#   select(geneticSexID, Ploidy_level, n_medium_trees, Cos.aspect, Elevation,
#          population_density_m2, 
#          size_mean, size_sd,
#          site_code, site_type,
#          year) %>%
#   group_by(site_code) %>%
#   slice_max(year) %>%
#   unique %>%
#   mutate(replicate=NUM_RESAMPLES) %>%
#   uncount(replicate) %>%
#   group_by(site_code) %>%
#   mutate(replicate=row_number())
# #df_sites_for_ipm$lambda = NA
# 
# 
# #age_structure_sites_all = matrix(NA, nrow=nrow(df_sites_for_ipm), ncol=mesh_points)
# 
# # try to map out all the sites
# lambdas_sites = pblapply(1:nrow(df_sites_for_ipm), function(i)
# {
#   ipm_this = NULL
#   
#   try(ipm_this <- make_ipm_for_site(
#     m_survival = m_survival_resampled[[ df_sites_for_ipm$replicate[i] ]],
#     m_growth = m_growth_resampled[[ df_sites_for_ipm$replicate[i] ]],
#     m_recruit = m_recruit_resampled[[ df_sites_for_ipm$replicate[i] ]],
#     other_vars = list(
#         geneticSexIDM=ifelse(df_sites_for_ipm$geneticSexID[i]=="M","1","0"),
#         Ploidy_levelTriploid=ifelse(df_sites_for_ipm$Ploidy_level[i]=="Triploid","1","0"),
#         n_medium_trees=df_sites_for_ipm$n_medium_trees[i],
#         Cos.aspect=df_sites_for_ipm$Cos.aspect[i],
#         Elevation=df_sites_for_ipm$Elevation[i]),
#       size_mean = df_sites_for_ipm$size_mean[i],
#       size_sd = df_sites_for_ipm$size_sd[i],
#       population_density_initial = df_sites_for_ipm$population_density_m2[i]
#   ))
#   if (!is.null(ipm_this))
#   {
#     lambda = mean(tail(as.numeric(ipmr::lambda(ipm_this)),5))
#     #age_structure_sites_all[i,] = ipm_this$pop_state$n_dbh[,MAX_ITERATIONS+1]
#   }
#   else
#   {
#     lambda = NA
#     #age_structure_sites_all[i,] = NA
#   }
#   
#   return(lambda)
#   #print(df_sites_for_ipm$lambda[i])
# }, cl=NUM_CORES)
# df_sites_for_ipm$lambda = as.numeric(lambdas_sites)
# 
# #write.csv(age_structure_sites_all, file='output_data/sites_age_structure.csv', row.names = FALSE)
# 
# 
# df_site_level = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv')
# df_sites_for_ipm_joined = df_sites_for_ipm %>%
#   left_join(df_site_level %>% 
#               select(site_code=Site_Code, X.UTM, Y.UTM, Watershed, 
#                      Elevation, Slope, 
#                      Canopy_openness,
#                      Summer.Insolation, Soil.type), by='site_code')# %>%
#   #mutate(lambda_binned = cut(lambda, breaks=c(0,0.99,1.01,Inf),labels=c('decreasing','stable','increasing')))
# 
# # df_sites_for_ipm_joined$dbh_modal = dbh_range[apply(age_structure_sites_all, 1, function(x) { 
# #   result = which.max(x)
# #   if (length(result)==0)
# #   {
# #     result = NA
# #   }
# #   return(result)
# #   })]
# 
# write.csv(df_sites_for_ipm_joined, file='output_data/ipm_outcomes_lambda_by_plot.csv', row.names = FALSE)
# 
# df_sites_for_ipm_joined = read.csv('output_data/ipm_outcomes_lambda_by_plot.csv')
# 
# 
# 
# df_sites_for_ipm_joined_summarized = df_sites_for_ipm_joined %>% 
#   filter(site_type=='random') %>%
#   group_by(site_code, X.UTM, Y.UTM) %>% 
#   summarize(lambda.mean=mean(lambda,na.rm=TRUE), 
#             lambda.sd=sd(lambda,na.rm=TRUE), 
#             p=ifelse(length(na.omit(lambda)) > 1, t.test(na.omit(lambda), mu=1)$p.value,NA)) %>%
#   mutate(p.adj = p.adjust(p, method="BH")) %>%
#   mutate(lambda_binned = ifelse(p<0.05, ifelse(lambda.mean < 1, 'decreasing', 'increasing'), 'stable'))
# 
# mean(df_sites_for_ipm_joined_summarized$lambda.mean,na.rm=TRUE)
# sd(df_sites_for_ipm_joined_summarized$lambda.mean,na.rm=TRUE)
# 
# table(df_sites_for_ipm_joined_summarized$lambda_binned)
# table(df_sites_for_ipm_joined_summarized$lambda_binned)/nrow(df_sites_for_ipm_joined_summarized)
# 
# g_map_lambda_binned = ggplot(df_sites_for_ipm_joined_summarized %>%
#                                filter(!is.na(lambda.mean)), 
#                            aes(x=X.UTM,y=Y.UTM, color=lambda.mean,shape=lambda_binned)) + 
#   geom_point(alpha=0.8) +
#   scale_color_gradient2(midpoint=1,
#                        low = 'darkorange',high='darkorchid1',mid = 'gray',
#                        limits=c(0.9,1.1),
#                         #limits=c(0.95,1.05),
#                         name=expression(paste(lambda, " (mean)"))) +
#   scale_shape_manual(name='Inference',values=c(4,3,16)) +
#   theme_bw() + 
#   coord_equal() +
#   xlab('Easting (m)') + ylab('Northing (m)')
# 
# g_lambda_hist = ggplot(df_sites_for_ipm_joined_summarized, aes(x=lambda.mean)) +
#   geom_histogram(binwidth = 0.005,fill='#333333') +
#   geom_vline(xintercept = 1,color='black') +
#   theme_bw() +
#   xlab(expression(paste(lambda, " (mean)"))) +
#   ylab('Count')
# 
# g_map_hist_lambda = ggarrange(g_map_lambda_binned, g_lambda_hist,
#           nrow=2,ncol=1,
#           labels='AUTO')
# 
# 
# ggsave(g_map_hist_lambda, file='output_figures/g_map_hist_lambda.png',width=6,height=7)
# ggsave(g_map_hist_lambda, file='output_figures/g_map_hist_lambda.pdf',width=6,height=7)
# 
# 
# # save info for ipm
# 
# save(list=c("make_ipm_for_site",
#           "transitions_all_filtered_joined_no_na",
#           "update_coefficients_full",
#           ls(pattern='^m_')),file='output_data/workspace for ipm.Rdata')
