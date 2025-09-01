library(MASS)
library(dplyr)
library(ipmr)
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


source('ipm parameters.R')
set.seed(1) # reproducibility of resampling

# load in data
transitions_all_filtered_joined = read.csv('output_data/transitions_all_filtered_joined.csv')
#message('could fit only the grid plots too')

# get no NA version for model selection
transitions_all_filtered_joined_no_na = transitions_all_filtered_joined %>%
  mutate(site_type = ifelse(nchar(site_code)==4,'grid','random')) %>%
  mutate(population_density_m2 = 10 / plot_area_m2) %>%
  group_by(site_code, year) %>%
  mutate(size_mean = mean(size), size_sd = sd(size)) %>%
  select(surv, sizeNext, recruit, size, 
                site_type, site_code, 
                delta_years,
                Ploidy_level, geneticSexID, 
                Cos.aspect, Elevation, 
                n_medium_trees, 
                population_density_m2, size_mean, size_sd) %>%
  # estimate a size_sd for the cases we don't have one
  ungroup %>%
  mutate(size_sd = ifelse(is.na(size_sd), mean(size_sd, na.rm=TRUE), size_sd)) %>%
  mutate(Ploidy_level=factor(Ploidy_level), geneticSexID=factor(geneticSexID))
#%>%
#  na.omit %>%
#  filter(site_type == 'random')

# add weights
counts_by_plot_year = transitions_all_filtered_joined_no_na %>% 
  group_by(year, site_code) %>% 
  tally %>%
  arrange(n) %>%
  mutate(weight = 1/n)

transitions_all_filtered_joined_no_na = transitions_all_filtered_joined_no_na %>%
  left_join(counts_by_plot_year %>% select(year, site_code, weight),
            by=c('year','site_code'))




# get resampled dataset for bootstrapping
transitions_resampled = lapply(1:NUM_RESAMPLES, function(x) { 
  transitions_all_filtered_joined_no_na[sample(1:nrow(transitions_all_filtered_joined_no_na), size=nrow(transitions_all_filtered_joined_no_na), replace=TRUE),] 
  })


# fit vital rate models
formula_survival = formula(factor(surv) ~ delta_years + size*(Ploidy_level*geneticSexID + Cos.aspect) + population_density_m2)
m_survival_all = glm(formula=formula_survival, 
                 data = transitions_all_filtered_joined_no_na, 
                 family = binomial(),
                 weights=round(10*transitions_all_filtered_joined_no_na$weight))

# m_survival_unweighted = glm(factor(surv) ~ delta_years + (size*(Ploidy_level*geneticSexID + Cos.aspect) + population_density_m2), 
#                      data = transitions_all_filtered_joined_no_na, 
#                      family = binomial())

m_survival_resampled = lapply(transitions_resampled, function(x) {
  glm(formula=formula_survival, 
      data = x, 
      family = binomial(),
      weights=round(10*transitions_all_filtered_joined_no_na$weight))
  })

pdf(file='output_figures/g_ipm_survival.pdf',width=10,height=10)
plot_model(m_survival_all, sort.est=TRUE) + 
  theme_bw() + 
  geom_hline(yintercept = 1) +
  ggtitle('survival, standardized effect')
plot_model(m_survival_all, type='int')
simulateResiduals(m_survival_all) %>% plot
dev.off()

tab_model(m_survival_all, file='output_figures/table_model_tree_level_survival.html')

Anova(m_survival_all)



formula_growth = formula(sizeNext ~ delta_years + size*(Ploidy_level*geneticSexID + Cos.aspect) + population_density_m2)
m_growth_all = glm(formula=formula_growth, 
               data = transitions_all_filtered_joined_no_na,
               weights=round(10*transitions_all_filtered_joined_no_na$weight))
m_growth_resampled = lapply(transitions_resampled, function(x) {
  lm(formula=formula_growth, 
     data = x,
     weights=round(10*transitions_all_filtered_joined_no_na$weight))
})

tab_model(m_growth_all, file='output_figures/table_model_tree_level_growth.html')


pdf(file='output_figures/g_ipm_growth.pdf',width=10,height=10)
plot_model(m_growth_all, type='est', sort.est=TRUE, rm.terms='size') + 
  theme_bw() + 
  geom_hline(yintercept = 0) +
  ggtitle('growth, standardized effect')
plot_model(m_growth_all, type='int')
simulateResiduals(m_growth_all) %>% plot
dev.off()

Anova(m_growth_all)




formula_recruit = formula(factor(recruit) ~ delta_years + n_medium_trees*(Ploidy_level*geneticSexID + Cos.aspect) + population_density_m2)
m_recruit_all = glm(formula=formula_recruit,
                data=transitions_all_filtered_joined_no_na,
                family = binomial(),
                weights=round(10*transitions_all_filtered_joined_no_na$weight))

m_recruit_resampled = lapply(transitions_resampled, function(x) {
  glm(formula=formula_recruit,
      data=x,
      family = binomial(),
      weights=round(10*transitions_all_filtered_joined_no_na$weight))
})

tab_model(m_recruit_all, file='output_figures/table_model_tree_level_recruit.html')


pdf(file='output_figures/g_ipm_recruit.pdf',width=10,height=10)
plot_model(m_recruit_all, sort.est=TRUE) + 
  theme_bw() + 
  geom_hline(yintercept = 1) +
  ggtitle('recruit, standardized effect')
plot_model(m_recruit_all, type='int')
simulateResiduals(m_recruit_all) %>% plot
dev.off()

Anova(m_recruit_all)





g_tree_level_rate_survival = plot(ggpredict(m_survival_all,terms=c('geneticSexID','Ploidy_level')),colors=c('blue','red')) + 
  labs(color='Cytotype') + ggtitle('Probability of survival') + ylab('') + xlab('Sex')
g_tree_level_rate_size = plot(ggpredict(m_growth_all,terms=c('geneticSexID','Ploidy_level')),colors=c('blue','red')) + 
  labs(color='Cytotype') + ggtitle('Size (cm)') + ylab('') + xlab('Sex')
g_tree_level_rate_recruit = plot(ggpredict(m_recruit_all,terms=c('geneticSexID','Ploidy_level')),colors=c('blue','red')) + 
  labs(color='Cytotype') + ggtitle('Probability of recruitment') + ylab('') + xlab('Sex')

g_tree_level_rates = ggarrange(g_tree_level_rate_survival, g_tree_level_rate_size, g_tree_level_rate_recruit,
          nrow=2,ncol=2,common.legend = TRUE, legend='bottom',labels='AUTO')
ggsave(g_tree_level_rates, file='output_figures/g_tree_level_rates.pdf',width=6,height=6)
ggsave(g_tree_level_rates, file='output_figures/g_tree_level_rates.png',width=6,height=6)
























# 
# # look at the distribution of resampled coefficients
# g_m_survival_resampled = sapply(m_survival_resampled, coef) %>% 
#   t %>%
#   as.data.frame %>%
#   pivot_longer(cols=everything()) %>%
#   ggplot(aes(x=name,y=value)) + 
#   geom_point() +
#   coord_flip() + 
#   theme_bw() +
#   ggtitle('Survival model bootstrap')
# ggsave(g_m_survival_resampled, file='output_figures/g_m_survival_resampled.pdf')
# 
# g_m_growth_resampled = sapply(m_growth_resampled, coef) %>% 
#   t %>%
#   as.data.frame %>%
#   pivot_longer(cols=everything()) %>%
#   ggplot(aes(x=name,y=value)) + 
#   geom_point() +
#   coord_flip() + 
#   theme_bw() +
#   ggtitle('Growth model bootstrap')
# ggsave(g_m_growth_resampled, file='output_figures/g_m_growth_resampled.pdf')
# 
# g_m_recruit_resampled = sapply(m_recruit_resampled, coef) %>% 
#   t %>%
#   as.data.frame %>%
#   pivot_longer(cols=everything()) %>%
#   ggplot(aes(x=name,y=value)) + 
#   geom_point() +
#   coord_flip() + 
#   theme_bw() +
#   ggtitle('Recruitment model bootstrap')
# ggsave(g_m_recruit_resampled, file='output_figures/g_m_recruit_resampled.pdf')
# 





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










# initial population structure
hist(transitions_all_filtered_joined_no_na$population_density_m2,breaks=mesh_points)
hist(transitions_all_filtered_joined_no_na$size)
# make initial size distribution to be the grand mean/sd of the training data

# set the recruitment conditions
make_ipm_for_site <- function(m_survival,
                              m_growth,
                              m_recruit,
                              other_vars, 
                              size_mean = mean(transitions_all_filtered_joined_no_na$size, na.rm=TRUE),
                              size_sd = sd(transitions_all_filtered_joined_no_na$size, na.rm=TRUE),
                              population_density_initial = mean(transitions_all_filtered_joined_no_na$population_density_m2))
{
  n_initial = sapply(1:mesh_points, function(i) {
    size = dnorm(dbh_range[i], mean=size_mean, sd=size_sd)
    size = ifelse(dbh_range[i] < DBH_min, 0, size)
    return(size)
  })
  # normalize the size to a population density of the requested value
  n_initial = n_initial * population_density_initial / sum(n_initial)
  
  plot(dbh_range, n_initial,type='h'); title(sprintf('population density=%.3f m-2',sum(n_initial)))
  
  
  
  # get size dependent coefficients, by setting population density to zero to eliminate its effect
  # also standardize to a 3-year interval to interpret the transitions
  coef_survival_size = update_coefficients_full(coef(m_survival), xvar='size', other_vars=c(other_vars,population_density_m2=0, delta_years=3))
  coef_growth_size = update_coefficients_full(coef(m_growth), xvar='size', other_vars=c(other_vars,population_density_m2=0, delta_years=3))
  coef_recruit_n_medium_trees = update_coefficients_full(coef(m_recruit), xvar='size', other_vars=c(other_vars,population_density_m2=0, delta_years=3))
  # get density dependent coefficients, assming there are no size interactions
  coef_growth_density = coef(m_growth)["population_density_m2"]
  coef_survival_density = coef(m_survival)["population_density_m2"]
  coef_recruit_density = coef(m_recruit)["population_density_m2"]
  
  
  my_data_list = list(s_int            = coef_survival_size[1],
                      s_slope_size     = coef_survival_size[2],
                      s_slope_density  = coef_survival_density[1],
                      g_int            = coef_growth_size[1],
                      g_slope_size     = coef_growth_size[2],
                      g_slope_density  = coef_growth_density[1],
                      sd_g             = sd(resid(m_growth)),
                      r_int            = coef_recruit_n_medium_trees[1],
                      r_slope_density  = coef_recruit_density[1],
                      r_d_mu           = DBH_min, # assume all come in at exactly size 5
                      r_d_sd           = 0.5) # assume there is a little bit of fluctuation in sizes (needed so discretization works too) 
  
  ipm_aspen <- init_ipm(sim_gen = "simple",
                            di_dd   = "dd",
                            det_stoch = "det")
  
  ipm_aspen <- define_kernel(
    
    proto_ipm = ipm_aspen,
    
    # Name of the kernel
    
    name      = "P",
    
    # The type of transition it describes (e.g. continuous - continuous, discrete - continuous).
    # These must be specified for all kernels!
    
    family    = "CC",
    
    # The formula for the kernel. We dont need to tack on the "z'/z"s here.  
    
    formula   = s * G,
    
    # A named set of expressions for the vital rates it includes. 
    # note the use of user-specified functions here. Additionally, each 
    # state variable has a stateVariable_1 and stateVariable_2, corresponding to
    # z and z' in the equations above. We don't need to define these variables ourselves,
    # just reference them correctly based on the way we've set up our model on paper.
    
    # Perform the inverse logit transformation to get survival probabilities
    # from your model. plogis from the "stats" package does this for us. 
    
    s         = plogis(s_int + s_slope_size * dbh_1 + s_slope_density * sum(n_dbh_t)), 
    
    # The growth model requires a function to compute the mean as a function of dbh.
    # The SD is a constant, so we don't need to define that in ... expression, 
    # just the data_list.
    
    G         = dnorm(dbh_2, mu_g, sd_g), # dbh_2 corresponds to size_next
    mu_g      = g_int + g_slope_size * dbh_1 + g_slope_density * sum(n_dbh_t),
    
    
    # Specify the constant parameters in the model in the data_list. 
    
    data_list = my_data_list,
    states    = list(c('dbh')),
    
    # If you want to correct for eviction, set evict_cor = TRUE and specify an
    # evict_fun. ipmr provides truncated_distributions() to help. This function
    # takes 2 arguments - the type of distribution, and the name of the parameter/
    # vital rate that it acts on.
    
    evict_cor = TRUE,
    evict_fun = truncated_distributions(fun    = 'norm',
                                        target = 'G')
  ) 
  
  ipm_aspen <- define_kernel(
    proto_ipm = ipm_aspen,
    name      = 'F',
    formula   = r_r * r_d,
    family    = 'CC',
    
    r_r       = plogis(r_int + r_slope_density * sum(n_dbh_t)), # overall rate of recruitment
    r_d       = dnorm(dbh_2, r_d_mu, r_d_sd), # multiply this rate by the appropriate size factor (using size_2, the new size)
    data_list = my_data_list,
    states    = list(c('dbh')),
    
    # Again, we'll correct for eviction in new recruits by
    # truncating the normal distribution.
    
    evict_cor = TRUE,
    evict_fun = truncated_distributions(fun    = 'norm',
                                        target = 'r_d')
  ) 
  
  # Next, we have to define the implementation details for the model. 
  # We need to tell ipmr how each kernel is integrated, what state
  # it starts on (i.e. z from above), and what state
  # it ends on (i.e. z' above). In simple_* models, state_start and state_end will 
  # always be the same, because we only have a single continuous state variable. 
  # General_* models will be more complicated.
  
  ipm_aspen <- define_impl(
    proto_ipm = ipm_aspen,
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep("midpoint", 2),
      state_start  = rep("dbh", 2),
      state_end    = rep("dbh", 2)
    )
  ) 
  
  ipm_aspen <- define_domains(
    proto_ipm = ipm_aspen,
    dbh = c(0, # the first entry is the lower bound of the domain.
            DBH_max, # the second entry is the upper bound of the domain.
            mesh_points # third entry is the number of meshpoints for the domain.
    ) 
  ) 
  
  # Next, we define the initial state of the population. We must do this because
  # ipmr computes everything through simulation, and simulations require a 
  # population state.
  
  ipm_aspen <- define_pop_state(
    proto_ipm = ipm_aspen,
    n_dbh = n_initial
  )
  
  ipm_aspen <- make_ipm(proto_ipm = ipm_aspen, 
                        iterations=MAX_ITERATIONS, 
                        normalize_pop_size = FALSE,
                        return_sub_kernels = TRUE)
  
  return(ipm_aspen)
}

ipm_male_triploid_south = make_ipm_for_site(
  m_survival = m_survival_all,
  m_growth = m_growth_all,
  m_recruit = m_recruit_all,
  other_vars=list(geneticSexIDM="0",
    Ploidy_levelTriploid="1",
    n_medium_trees=2,
    Cos.aspect=-1)
  )



# get lambda
plot(lambda(ipm_male_triploid_south),xlab='time step',ylab='lambda')

# look at trajectory of age distribution
image(ipm_male_triploid_south$pop_state$n_dbh %>% sqrt,xlab='size',ylab='time')

# look at final age structure (sqrt transformed)
plot(dbh_range, 
     sqrt(ipm_male_triploid_south$pop_state$n_dbh[,MAX_ITERATIONS+1]),
     type='h',xlab='DBH (cm)',ylab='sqrt # m-2')


# # look at kernels at initial densities
# pdf(file='output_figures/g_ipm_kernels.pdf',width=10,height=5)
# par(mfrow=c(1,2))
# plot(ipm_male_triploid_south$sub_kernels$P_it_1); title('P subkernel')
# plot(ipm_all$sub_kernels$F_it_1); title('F subkernel')
# dev.off()

kernel_for_plotting <- function(k)
{
  k_long = k %>% 
    as.matrix %>% 
    as.data.frame %>% 
    mutate(row=row_number()) %>% 
    pivot_longer(!row) %>%
    mutate(sizeTo=dbh_range[row]) %>%
    mutate(sizeFrom=dbh_range[as.numeric(gsub("V","",name))])
}

plot_kernel <- function(k, name)
{
  ggplot(kernel_for_plotting(k), 
         aes(x=sizeFrom,y=sizeTo,fill=value)) +
    geom_raster() +
    theme_bw() +
    scale_fill_viridis_c(option='C',name='Value') +
    xlab('Size (cm)') + ylab('Size next (cm)') +
    geom_abline(slope=1,intercept=0,color='white',alpha=0.5) +
    coord_equal() +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(sprintf('%s sub-kernel',name)) 
}

k1 = plot_kernel(ipm_male_triploid_south$sub_kernels$P_it_200, 'P')
k2 = plot_kernel(ipm_male_triploid_south$sub_kernels$F_it_200, 'F')
ggsave(ggarrange(k1, k2, labels='AUTO', align='hv'),file='output_figures/g_ipm_kernel_male_triploid_south.pdf',width=7,height=2.5)
ggsave(ggarrange(k1, k2, labels='AUTO', align='hv'),file='output_figures/g_ipm_kernel_male_triploid_south.png',width=7,height=2.5)




# now look at performance variation across parameters

params = expand.grid(replicate=1:NUM_RESAMPLES,
                     geneticSexIDM=c("0","1"), 
                     Ploidy_levelTriploid=c("0","1"), 
                     n_medium_trees=c(0,5,15),
                     Cos.aspect=c(-1,0,1))

lambdas = pbsapply(1:nrow(params), function(i) {
  
  cat(sprintf('%d %.3f\n',i, i/nrow(params)))
  
  ipm_this = NULL
  try(ipm_this <- make_ipm_for_site(
      m_survival = m_survival_resampled[[ params$replicate[i] ]],
      m_growth = m_growth_resampled[[ params$replicate[i] ]],
      m_recruit = m_recruit_resampled[[ params$replicate[i] ]],
      other_vars = as.list(params[i,])
    ))
  
  if (!is.null(ipm_this))
  {
    lambda = mean(tail(as.numeric(lambda(ipm_this)),5)) # take the final 5 lambda values and average
  }
  else
  {
    lambda = NA
  }
  print(lambda)
 
  # get the mean lambda over time
  return(lambda)
}, cl = NUM_CORES)
params$lambda = lambdas

write.csv(params, file='output_data/ipm_outcomes_lambda.csv', row.names = FALSE)


params_for_plotting = params %>%
  mutate(sex=ifelse(geneticSexIDM==1,'M','F')) %>%
  mutate(cytotype=ifelse(Ploidy_levelTriploid==1,'triploid','diploid')) %>%
  mutate(n_medium_trees=factor(paste('# medium=',n_medium_trees,sep=''),levels=paste('# medium=',c(0,5,15),sep=''),ordered=TRUE)) %>%
  mutate(Cos.aspect=paste('Cosine aspect=',Cos.aspect,sep=''))


g_lambda =  ggplot(params_for_plotting, aes(color=factor(cytotype),x=factor(sex),y=lambda)) + 
  geom_hline(yintercept = 1) +
  geom_boxplot(outliers = FALSE) +
  geom_point(alpha=0.5, position=position_jitterdodge()) +
  # geom_point(alpha=0.5) + 
  # geom_jitter(height=0,width=0.2) + 
  #geom_bar(stat='identity',alpha=0.5,position='dodge') +
  facet_grid(Cos.aspect~n_medium_trees) +
  theme_bw() +
  #scale_y_log10() +
  scale_color_manual(values=c('blue','red'),name='Cytotype') +
  xlab('Sex') +
  ylab(expression(paste(lambda)))

ggsave(g_lambda, file='output_figures/g_ipm_lambda.pdf',width=6,height=6)
ggsave(g_lambda, file='output_figures/g_ipm_lambda.png',width=6,height=6)










# now map lambdas for all the plots
df_sites_for_ipm = transitions_all_filtered_joined_no_na %>%
  select(geneticSexID, Ploidy_level, n_medium_trees, Cos.aspect,
         population_density_m2, 
         size_mean, size_sd,
         site_code, site_type,
         year) %>%
  group_by(site_code) %>%
  slice_max(year) %>%
  unique %>%
  mutate(replicate=NUM_RESAMPLES) %>%
  uncount(replicate) %>%
  group_by(site_code) %>%
  mutate(replicate=row_number())
#df_sites_for_ipm$lambda = NA


#age_structure_sites_all = matrix(NA, nrow=nrow(df_sites_for_ipm), ncol=mesh_points)

# try to map out all the sites
lambdas_sites = pblapply(1:nrow(df_sites_for_ipm), function(i)
{
  ipm_this = NULL
  
  try(ipm_this <- make_ipm_for_site(
    m_survival = m_survival_resampled[[ df_sites_for_ipm$replicate[i] ]],
    m_growth = m_growth_resampled[[ df_sites_for_ipm$replicate[i] ]],
    m_recruit = m_recruit_resampled[[ df_sites_for_ipm$replicate[i] ]],
    other_vars = list(
      geneticSexIDM=ifelse(df_sites_for_ipm$geneticSexID[i]=="M","1","0"),
      Ploidy_levelTriploid=ifelse(df_sites_for_ipm$Ploidy_level[i]=="Triploid","1","0"),
      n_medium_trees=df_sites_for_ipm$n_medium_trees[i],
      Cos.aspect=df_sites_for_ipm$Cos.aspect[i]),
      size_mean = df_sites_for_ipm$size_mean[i],
      size_sd = df_sites_for_ipm$size_sd[i],
      population_density_initial = df_sites_for_ipm$population_density_m2[i]
  ))
  if (!is.null(ipm_this))
  {
    lambda = mean(tail(as.numeric(lambda(ipm_this)),5))
    #age_structure_sites_all[i,] = ipm_this$pop_state$n_dbh[,MAX_ITERATIONS+1]
  }
  else
  {
    lambda = NA
    #age_structure_sites_all[i,] = NA
  }
  
  return(lambda)
  #print(df_sites_for_ipm$lambda[i])
}, cl=NUM_CORES)
df_sites_for_ipm$lambda = as.numeric(lambdas_sites)

#write.csv(age_structure_sites_all, file='output_data/sites_age_structure.csv', row.names = FALSE)


df_site_level = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv')
df_sites_for_ipm_joined = df_sites_for_ipm %>%
  left_join(df_site_level %>% 
              select(site_code=Site_Code, X.UTM, Y.UTM, Watershed, 
                     Elevation, Slope, 
                     Canopy_openness,
                     Summer.Insolation, Soil.type), by='site_code')# %>%
  #mutate(lambda_binned = cut(lambda, breaks=c(0,0.99,1.01,Inf),labels=c('decreasing','stable','increasing')))

# df_sites_for_ipm_joined$dbh_modal = dbh_range[apply(age_structure_sites_all, 1, function(x) { 
#   result = which.max(x)
#   if (length(result)==0)
#   {
#     result = NA
#   }
#   return(result)
#   })]

write.csv(df_sites_for_ipm_joined, file='output_data/ipm_outcomes_lambda_by_plot.csv', row.names = FALSE)




df_sites_for_ipm_joined_summarized = df_sites_for_ipm_joined %>% 
  filter(site_type=='random') %>%
  group_by(site_code, X.UTM, Y.UTM) %>% 
  summarize(lambda.mean=mean(lambda,na.rm=TRUE), 
            lambda.sd=sd(lambda,na.rm=TRUE), 
            p=ifelse(length(na.omit(lambda)) > 1, t.test(na.omit(lambda), mu=1)$p.value,NA)) %>%
  mutate(p.adj = p.adjust(p, method="BH")) %>%
  mutate(lambda_binned = ifelse(p<0.05, ifelse(lambda.mean < 1, 'decreasing', 'increasing'), 'stable'))

g_map_lambda_binned = ggplot(df_sites_for_ipm_joined_summarized %>%
                               filter(!is.na(lambda.mean)), 
                           aes(x=X.UTM,y=Y.UTM, color=lambda.mean,shape=lambda_binned)) + 
  geom_point(alpha=0.8) +
  scale_color_gradient2(low='red',high='blue',mid='gray',midpoint=1,
                        #limits=c(0.95,1.05),
                        name=expression(paste(lambda, " (mean)"))) +
  scale_shape_manual(name='Inference',values=c(16,3)) +
  theme_bw() + 
  coord_equal() +
  xlab('Easting (m)') + ylab('Northing (m)')

g_lambda_hist = ggplot(df_sites_for_ipm_joined_summarized, aes(x=lambda.mean)) +
  geom_histogram(binwidth = 0.01,fill='lightgray') +
  geom_vline(xintercept = 1,color='black') +
  theme_bw() +
  xlab(expression(paste(lambda, " (mean)"))) +
  ylab('Count')

g_map_hist_lambda = ggarrange(g_map_lambda_binned, g_lambda_hist,
          nrow=2,ncol=1,
          labels='AUTO')


ggsave(g_map_hist_lambda, file='output_figures/g_map_hist_lambda.png',width=6,height=7)
ggsave(g_map_hist_lambda, file='output_figures/g_map_hist_lambda.pdf',width=6,height=7)


# save info for ipm

save(list=c("make_ipm_for_site",
          "transitions_all_filtered_joined_no_na",
          "update_coefficients_full",
          ls(pattern='^m_')),file='output_data/workspace for ipm.Rdata')
