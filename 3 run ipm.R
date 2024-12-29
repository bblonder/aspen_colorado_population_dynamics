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


source('ipm parameters.R')

# load in data
transitions_all_filtered_joined = read.csv('output_data/transitions_all_filtered_joined.csv')
message('could fit only the grid plots too')

# get no NA version for model selection
transitions_all_filtered_joined_no_na = transitions_all_filtered_joined %>%
  mutate(site_type = ifelse(nchar(site_code)==4,'grid','random')) %>%
  mutate(population_density_m2 = 10 / plot_area_m2) %>%
  group_by(site_code, year) %>%
  mutate(size_mean = mean(size), size_sd = sd(size)) %>%
  select(surv, sizeNext, recruit, size, 
                site_type, site_code, 
                Ploidy_level, geneticSexID, 
                Cos.aspect, 
                n_medium_trees, 
                population_density_m2, size_mean, size_sd) %>%
  na.omit %>%
  filter(site_type == 'random')


# fit vital rate models
m_survival = glm(factor(surv) ~ size*(Ploidy_level*geneticSexID + Cos.aspect) + population_density_m2, 
                 data = transitions_all_filtered_joined_no_na, 
                 family = binomial())

pdf(file='output_figures/g_ipm_survival.pdf',width=10,height=10)
plot_model(m_survival, type='std', sort.est=TRUE) + 
  theme_bw() + 
  geom_hline(yintercept = 1) +
  ggtitle('survival, standardized effect')
plot_model(m_survival, type='int')
simulateResiduals(m_survival) %>% plot
dev.off()

coef(m_survival)
Anova(m_survival)
# simulateResiduals(m_survival) %>% plot
# #visreg(m_survival,xvar='size',by='geneticSexID',gg=TRUE, overlay=TRUE)
# visreg(m_survival,xvar='size',by='Ploidy_level',gg=TRUE, overlay=TRUE)


m_growth = lm(sizeNext ~ size*(Ploidy_level*geneticSexID + Cos.aspect) + population_density_m2, 
               data = transitions_all_filtered_joined_no_na)

pdf(file='output_figures/g_ipm_growth.pdf',width=10,height=10)
plot_model(m_growth, type='est', sort.est=TRUE, rm.terms='size') + 
  theme_bw() + 
  geom_hline(yintercept = 0) +
  ggtitle('growth, standardized effect')
plot_model(m_growth, type='int')
simulateResiduals(m_growth) %>% plot
dev.off()

coef(m_growth)
Anova(m_growth)
# simulateResiduals(m_growth) %>% plot
# visreg(m_growth,xvar='size',by='geneticSexID',gg=TRUE, overlay=TRUE)
# visreg(m_growth,xvar='size',by='Ploidy_level',gg=TRUE, overlay=TRUE)

# n small trees does not matter, basal area density might
m_recruit = glm(factor(recruit) ~ n_medium_trees*(Ploidy_level*geneticSexID + Cos.aspect) + population_density_m2,
                data=transitions_all_filtered_joined_no_na,family = binomial())

pdf(file='output_figures/g_ipm_recruit.pdf',width=10,height=10)
plot_model(m_recruit, type='std', sort.est=TRUE) + 
  theme_bw() + 
  geom_hline(yintercept = 1) +
  ggtitle('recruit, standardized effect')
plot_model(m_recruit, type='int')
simulateResiduals(m_recruit) %>% plot
dev.off()

coef(m_recruit)
Anova(m_recruit)
simulateResiduals(m_recruit) %>% plot
# visreg(m_recruit,xvar='n_medium_trees',by='geneticSexID',gg=TRUE, overlay=TRUE)
# #visreg(m_recruit,xvar='n_medium_trees',by='Ploidy_level',gg=TRUE, overlay=TRUE)

# update_coefficients <- function(coefs, xvar, is_triploid, is_male)
# {
#   intercept = coefs["(Intercept)"]
#   slope = coefs[xvar]
#   
#   if (is_triploid==TRUE)
#   {
#     intercept = intercept + coefs["Ploidy_levelTriploid"]
#     slope = slope + coefs[sprintf("%s:Ploidy_levelTriploid",xvar)]
#   }
#   if (is_male==TRUE)
#   {
#     intercept = intercept + coefs["geneticSexIDM"]
#     slope = slope + coefs[sprintf("%s:geneticSexIDM",xvar)]
#   }
#   
#   if (is_triploid==TRUE & is_male==TRUE)
#   {
#     intercept = intercept + coefs["Ploidy_levelTriploid:geneticSexIDM"]
#     slope = slope + coefs[sprintf("%s:Ploidy_levelTriploid:geneticSexIDM",xvar)]
#   }
#   
#   return(c(intercept=as.numeric(intercept), slope=as.numeric(slope)))
# }

update_coefficients_full <- function(coefficients, xvar, othervars)
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
  for (i in 1:length(othervars))
  {
    names(terms_intercept) = gsub(names(othervars)[i], othervars[i], names(terms_intercept))
    names(terms_slope) = gsub(names(othervars)[i], othervars[i], names(terms_slope))
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
make_ipm_for_site <- function(othervars, 
                              size_mean = mean(transitions_all_filtered_joined_no_na$size),
                              size_sd = sd(transitions_all_filtered_joined_no_na$size),
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
  coef_survival_size = update_coefficients_full(coef(m_survival), xvar='size', othervars=c(othervars,population_density_m2=0))
  coef_growth_size = update_coefficients_full(coef(m_growth), xvar='size', othervars=c(othervars,population_density_m2=0))
  coef_recruit_n_medium_trees = update_coefficients_full(coef(m_recruit), xvar='size', othervars=c(othervars,population_density_m2=0))
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
                      r_d_sd           = 0.5) # assume there is a little bit of fluctuation in sizes (needed so discretizaiton works too)s 
  
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

ipm_this = make_ipm_for_site(othervars=list(geneticSexIDM="0",
                                            Ploidy_levelTriploid="1",
                                            n_medium_trees=0,
                                            Cos.aspect=-1))



# get lambda
plot(lambda(ipm_this))

# look at trajectory of age distribution
image(ipm_this$pop_state$n_dbh %>% sqrt,xlab='size',ylab='time')

# look at final age structure (sqrt transformed)
plot(dbh_range, 
     sqrt(ipm_this$pop_state$n_dbh[,MAX_ITERATIONS+1]),
     type='h',xlab='DBH (cm)',ylab='sqrt # m-2')

# look at age structure     
# w_ipmr      <- right_ev(ipm_this)
# plot(dbh_range, w_ipmr$dbh_w,type='h')

# look at kernels at initial densities
pdf(file='output_figures/g_ipm_kernels.pdf',width=10,height=5)
par(mfrow=c(1,2))
plot(ipm_this$sub_kernels$P_it_1); title('P subkernel')
plot(ipm_this$sub_kernels$F_it_1); title('F subkernel')
dev.off()

rm(ipm_this)







# now look at performance variation across parameters

params = expand.grid(geneticSexIDM=c("0","1"), Ploidy_levelTriploid=c("0","1"), n_medium_trees=c(0,2,4,6),Cos.aspect=c(-1,-0.5,0,0.5,1))

lambdas = sapply(1:nrow(params), function(i) {
  ipm_this = make_ipm_for_site(as.list(params[i,]))

  cat(sprintf('%.3f\n',i/nrow(params)))
  
  # get the mean lambda over time
  return(mean(as.numeric(lambda(ipm_this))))
})
params$lambda = lambdas

write.csv(params, file='output_data/ipm_outcomes_lambda.csv', row.names = FALSE)











# now map lambdas for all the plots
df_sites_for_ipm = transitions_all_filtered_joined_no_na %>%
  select(geneticSexID, Ploidy_level, n_medium_trees, Cos.aspect, 
         population_density_m2, 
         size_mean, size_sd,
         site_code, site_type) %>%
  group_by(site_code) %>%
  slice_max(year) %>%
  unique
df_sites_for_ipm$lambda = NA

age_structure_sites_for_ipm = matrix(NA, nrow=nrow(df_sites_for_ipm), ncol=mesh_points)

# try to map out all the sites
for (i in 1:nrow(df_sites_for_ipm))
{
  print(i)
  ipm_this = make_ipm_for_site(othervars = list(
    geneticSexIDM=ifelse(df_sites_for_ipm$geneticSexID[i]=="M","1","0"),
    Ploidy_levelTriploid=ifelse(df_sites_for_ipm$Ploidy_level[i]=="Triploid","1","0"),
    n_medium_trees=df_sites_for_ipm$n_medium_trees[i],
    Cos.aspect=df_sites_for_ipm$Cos.aspect[i]),
    size_mean = df_sites_for_ipm$size_mean[i],
    size_sd = df_sites_for_ipm$size_sd[i],
    population_density_initial = df_sites_for_ipm$population_density_m2[i]
  )
  df_sites_for_ipm$lambda[i] = mean(as.numeric(lambda(ipm_this)))
  age_structure_sites_for_ipm[i,] = ipm_this$pop_state$n_dbh[,MAX_ITERATIONS+1]
  
  print(df_sites_for_ipm$lambda[i])
}

write.csv(age_structure_sites_for_ipm, file='output_data/sites_age_structure.csv', row.names = FALSE)


df_site_level = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv')
df_sites_for_ipm_joined = df_sites_for_ipm %>%
  left_join(df_site_level %>% 
              select(site_code=Site_Code, X.UTM, Y.UTM, Watershed, 
                     Elevation, Slope, 
                     Canopy_openness,
                     Cow_Damage, Summer.Insolation, Soil.type), by='site_code') %>%
  mutate(lambda_binned = cut(lambda, breaks=c(0,0.99,1.01,Inf),labels=c('decreasing','stable','increasing')))

df_sites_for_ipm_joined$dbh_modal = dbh_range[apply(age_structure_sites_for_ipm, 1, which.max)]

write.csv(df_sites_for_ipm_joined, file='output_data/sites_data_frame.csv', row.names = FALSE)