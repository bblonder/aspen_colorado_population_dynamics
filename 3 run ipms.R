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
# library(popbio)
# library(boot)
library(insight)
library(progress)

#source('ipm parameters.R')
NUM_RESAMPLES = 20

set.seed(1) # reproducibility of resampling

# load in data
transitions_all_filtered_joined = read.csv('output_data/transitions_all_filtered_joined_with_climate_year_final.csv')
#message('could fit only the grid plots too')

# get no NA version for model selection
transitions_all_filtered_joined_no_na = transitions_all_filtered_joined %>%
  mutate(site_type = ifelse(nchar(site_code)==4,'grid','random')) %>%
  mutate(population_density_m2 = 10 / plot_area_m2) %>%
  group_by(site_code, year) %>%
  mutate(size_mean = mean(size), size_sd = sd(size)) %>%
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

# add weights
counts_by_plot_year = transitions_all_filtered_joined_no_na %>% 
  group_by(year, site_code) %>% 
  tally %>%
  arrange(n) %>%
  mutate(weight = 1/n)

transitions_all_filtered_joined_no_na = transitions_all_filtered_joined_no_na %>%
  left_join(counts_by_plot_year %>% dplyr::select(year, site_code, weight),
            by=c('year','site_code'))




# get resampled dataset for bootstrapping
transitions_resampled = lapply(1:NUM_RESAMPLES, function(x) { 
  transitions_all_filtered_joined_no_na[sample(1:nrow(transitions_all_filtered_joined_no_na), size=nrow(transitions_all_filtered_joined_no_na), replace=TRUE),] 
  })


# fit vital rate models
formula_survival = formula(factor(surv) ~ delta_years + size * (Ploidy_level + Cos.aspect + Elevation + STB.0 + SWE.0) + population_density_m2)
m_survival_all = glm(formula=formula_survival, 
                 data = transitions_all_filtered_joined_no_na, 
                 family = binomial(),
                 weights=round(10*transitions_all_filtered_joined_no_na$weight))
summary(m_survival_all)
 
pdf(file='output_figures/g_ipm_survival.pdf',width=10,height=10)
plot_model(m_survival_all, sort.est=TRUE) + 
  theme_bw() + 
  geom_hline(yintercept = 1) +
  ggtitle('survival, standardized effect')
plot_model(m_survival_all, type='int')
simulateResiduals(m_survival_all) %>% plot
dev.off()

tab_model(m_survival_all, file='output_figures/table_model_tree_level_survival.html')

Anova(m_survival_all,type=3)



formula_growth = formula(sizeNext ~ delta_years + size * (geneticSexID) + population_density_m2)
m_growth_all = glm(formula=formula_growth, 
               data = transitions_all_filtered_joined_no_na,
               weights=round(10*transitions_all_filtered_joined_no_na$weight))
summary(m_growth_all)


tab_model(m_growth_all, file='output_figures/table_model_tree_level_growth.html')


pdf(file='output_figures/g_ipm_growth.pdf',width=10,height=10)
plot_model(m_growth_all, type='est', sort.est=TRUE, rm.terms='size') + 
  theme_bw() + 
  geom_hline(yintercept = 0) +
  ggtitle('growth, standardized effect')
plot_model(m_growth_all, type='int')
simulateResiduals(m_growth_all) %>% plot
dev.off()

Anova(m_growth_all,type=3)




counts_recruit_all = transitions_all_filtered_joined_no_na %>%
  group_by(site_code, weight, delta_years, year, n_medium_trees, Ploidy_level, geneticSexID, Cos.aspect, Elevation, STB.0, SWE.0, population_density_m2) %>%
  tally(recruit) %>%
  mutate(log_ratio_n_med_per_pop_dens = log(1+n_medium_trees/population_density_m2)) %>%
  as.data.frame

formula_recruit_all = formula(n ~ delta_years + n_medium_trees + population_density_m2 + geneticSexID)

m_recruit_count_all = glm(formula = formula_recruit_all,
                          data=counts_recruit_all,
                          family = poisson,
                          weights=round(10*counts_recruit_all$weight))
summary(m_recruit_count_all)

pdf(file='output_figures/g_ipm_recruit_count.pdf',width=10,height=10)
plot_model(m_recruit_count_all, sort.est=TRUE) + 
  theme_bw() + 
  geom_hline(yintercept = 1) +
  ggtitle('recruit count, standardized effect')
#plot_model(m_recruit_count_all, type='int')
simulateResiduals(m_recruit_count_all) %>% plot
dev.off()

tab_model(m_recruit_count_all, file='output_figures/table_model_tree_level_recruit_count.html')

Anova(m_recruit_count_all,type=3)



# this is effectively a no-intercept model for n_medium
formula_n_medium = formula(log_ratio_n_med_per_pop_dens ~ Cos.aspect + Ploidy_level + population_density_m2)
m_n_medium = glm(formula_n_medium, 
                 family=gaussian,
                 data=counts_recruit_all,
                weights=round(10*counts_recruit_all$weight))
 
summary(m_n_medium)

pdf(file='output_figures/g_ipm_n_medium.pdf',width=10,height=10)
plot_model(m_n_medium, sort.est=TRUE) + 
  theme_bw() + 
  geom_hline(yintercept = 1) +
  ggtitle('n medium, standardized effect')
#plot_model(m_recruit_count_all, type='int')
simulateResiduals(m_n_medium) %>% plot
dev.off()

tab_model(m_n_medium, file='output_figures/table_model_tree_level_n_medium.html')

Anova(m_n_medium,type=3)


# fit size residuals
m_growth_resid2 = resid(m_growth_all)^2
m_growth_size = m_growth_all$data[complete.cases(m_growth_all$data[,all.vars(m_growth_all$formula),with=FALSE]),]$size
df_m_growth = data.frame(resid2=m_growth_resid2, size=m_growth_size)

# ok to ditch outliers?
m_growth_size_variance = glm(I(sqrt(resid2))~size,data=df_m_growth %>% filter(resid2 < 30),family = Gamma)
#plot(m_growth_size_variance)
summary(m_growth_size_variance)

pdf(file='output_figures/g_ipm_growth_size_variance.pdf',width=10,height=10)
plot_model(m_growth_size_variance, sort.est=TRUE) + 
  theme_bw() + 
  geom_hline(yintercept = 1) +
  ggtitle('growth size variance')
#plot_model(m_recruit_count_all, type='int')
simulateResiduals(m_growth_size_variance) %>% plot
dev.off()

tab_model(m_growth_size_variance, file='output_figures/table_model_tree_level_growth_size_variance.html')

Anova(m_growth_size_variance,type=3)









Anova(m_survival_all,type=3) %>% 
  as.data.frame %>% 
  dplyr::select(pvalue=`Pr(>Chisq)`) %>% 
  mutate(pvalue=format_p(pvalue,name=NULL)) %>%
  write.csv(file='output_figures/table_anova_survival_tree_level.csv',row.names=TRUE)

Anova(m_growth_all,type=3) %>% 
  as.data.frame %>% 
  dplyr::select(pvalue=`Pr(>Chisq)`) %>% 
  mutate(pvalue=format_p(pvalue,name=NULL)) %>%
  write.csv(file='output_figures/table_anova_growth_tree_level.csv',row.names=TRUE)

Anova(m_recruit_count_all,type=3) %>% 
  as.data.frame %>% 
  dplyr::select(pvalue=`Pr(>Chisq)`) %>% 
  mutate(pvalue=format_p(pvalue,name=NULL)) %>%
  write.csv(file='output_figures/table_anova_survival_recruit_count_tree_level.csv',row.names=TRUE)

Anova(m_n_medium,type=3) %>% 
  as.data.frame %>% 
  dplyr::select(pvalue=`Pr(>Chisq)`) %>% 
  mutate(pvalue=format_p(pvalue,name=NULL)) %>%
  write.csv(file='output_figures/table_anova_n_medium_tree_level.csv',row.names=TRUE)

Anova(m_growth_size_variance,type=3) %>% 
  as.data.frame %>% 
  dplyr::select(pvalue=`Pr(>Chisq)`) %>% 
  mutate(pvalue=format_p(pvalue,name=NULL)) %>%
  write.csv(file='output_figures/table_anova_growth_size_variance_tree_level.csv',row.names=TRUE)



g_tree_level_rate_survival = plot(ggpredict(m_survival_all,terms=c('SWE.0','Ploidy_level')),colors=c('blue','red','gray')) + 
  labs(color='Cytotype') + ggtitle('Probability of survival') + ylab('') + xlab('SWE.0')
g_tree_level_rate_size = plot(ggpredict(m_growth_all,terms=c('size','geneticSexID')),colors=c('blue','red','gray')) + 
  labs(color='Cytotype') + ggtitle('Size (cm)') + ylab('') + xlab('Size')
g_tree_level_rate_recruit = plot(ggpredict(m_recruit_count_all,terms=c('n_medium_trees','geneticSexID')),colors=c('blue','red','gray')) + 
  labs(color='Cytotype') + ggtitle('# recruits') + ylab('') + xlab('# medium trees')
g_tree_level_rate_n_medium = plot(ggpredict(m_n_medium,terms=c('Cos.aspect','Ploidy_level')),colors=c('blue','red','gray')) + 
  labs(color='Cytotype') + ggtitle('# medium trees (log x+1)') + ylab('') + xlab('Cos aspect')

g_tree_level_rates = ggarrange(g_tree_level_rate_survival, g_tree_level_rate_size, g_tree_level_rate_recruit, g_tree_level_rate_n_medium,
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




run_model <- function(prefix_this='test',
                      minsize=5, # L
                      maxsize=60, # U
                      m=100, # mesh points
                      n_iterations=20, # time steps
                      n_S_initial=0, # initial # saplings
                      population_density_initial=0.1, # initial density
                      size_mean_initial = mean(transitions_all_filtered_joined_no_na$size, na.rm=TRUE), # initial tagged tree size distribution
                      size_sd_initial = sd(transitions_all_filtered_joined_no_na$size, na.rm=TRUE), # initial tagged tree size distribution
                      Ploidy_leveltriploid='1',
                      geneticSexIDmale='0',
                      Cos.aspect=-1,
                      Elevation=3000,
                      STB.0=0,
                      SWE.0=100,
                      delta_years=3,
                      ns_factor_density=1, # fudge factor to account for ns density effect
                      ns_factor_survival=1 # fudge factor to account for ns survival effect
)
{
  # conditions
  current_conditions = list(Ploidy_leveltriploid=Ploidy_leveltriploid,
                            geneticSexIDmale=geneticSexIDmale,
                            Ploidy_levelunknown='0',
                            geneticSexIDunknown='0',
                            Cos.aspect=Cos.aspect,
                            Elevation=Elevation,
                            STB.0=STB.0,
                            SWE.0=SWE.0,
                            population_density_m2=population_density_initial, 
                            delta_years=delta_years)
  
  # set up bins
  size_bins = seq(minsize, maxsize, length.out=m)
  
  n_T = sapply(size_bins, function(x) {
    size = dnorm(x, mean=size_mean_initial, sd=size_sd_initial)
    return(size)
  })
  # normalize the size to a population density of the requested value
  n_T = n_T * current_conditions$population_density_m2 / sum(n_T)
  
  # initial recruits
  n_S = n_S_initial
  
  # number of adult trees
  n_T_time_series = matrix(NA, nrow=n_iterations, ncol=m)
  # number of saplings
  n_S_time_series = rep(NA, n_iterations)
  n_S_component_time_series = data.frame(matrix(NA, nrow=n_iterations, ncol=3))
  names(n_S_component_time_series) = c('density','survival','progression')
  
  # kernel
  P_initial = NULL
  P_final = NULL

  ### do iteration
  pb = progress_bar$new(total=n_iterations)
  for (i in 1:n_iterations)
  {
    n_T_time_series[i,] = n_T
    n_S_time_series[i] = n_S
    
    # update population density
    current_conditions$population_density_m2 = sum(n_T)
    
    # get coefficients on link scale
    # binomial
    coef_survival = update_coefficients_full(coef(m_survival_all), xvar='size', other_vars=current_conditions)
    # gaussian
    coef_growth = update_coefficients_full(coef(m_growth_all), xvar='size', other_vars=current_conditions)
    # poisson
    coef_recruit_count = update_coefficients_full(coef(m_recruit_count_all), xvar='n_medium_trees', other_vars=current_conditions)
    # gamma
    coef_size_variance = as.numeric(coef(m_growth_size_variance))
    # gaussian
    coef_n_medium = update_coefficients_full(coef(m_n_medium), xvar='size', other_vars = current_conditions)
    
    #
    sx <- function(x) {
      xbeta = coef_survival[1] + x*coef_survival[2]
      mu = exp(xbeta)/(1 + exp(xbeta)) #
      return(mu)
    }
    
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
      L<-minsize; U<-maxsize; # i changed this! OK to integrate over this range only?
      
      # boundary points b and mesh points y
      b<-L+c(0:n)*(U-L)/n;
      y<-0.5*(b[1:n]+b[2:(n+1)]);
      h<-(U-L)/n
      
      # loop to construct the matrix
      M<-matrix(0,n,n);
      for (i in 1:n){
        #cat(i);
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
  
    # iterate sapling trees
    n_S_density_effect = (exp(coef_n_medium)-1)*sum(n_T)*ns_factor_density # gaussian link but need to back transform the log(x+1) and then multiply by pop density.
    # define prior year recruit survival
    # assume medium trees survive like 5 cm adults
    n_S_survival_effect = sx(5)*n_S*ns_factor_survival # no link, check the 1 factor
    n_S_progression_effect = exp(coef_recruit_count[1] + coef_recruit_count[2] * n_S) # poisson link
  
    n_S_next = n_S_density_effect + n_S_survival_effect - n_S_progression_effect
    # hack - clamp this
    if (n_S_next<0)
    {
      # hack!
      # print('clamping nS at 0')
      n_S_next = 0
    }
  
    # iterate tagged trees
    n_T_next = P$matrix %*% n_T + P$h * c_R(P$meshpts) * n_S_progression_effect  
    
    # store iterated value
    n_T = n_T_next
    n_S = n_S_next
    
    # this logging has to come here based on how i organized... (awkward)
    n_S_component_time_series[i,'density'] = n_S_density_effect
    n_S_component_time_series[i,'survival'] = n_S_survival_effect
    n_S_component_time_series[i,'progression'] = n_S_progression_effect
    
    # save kernels
    if (i==1)
    {
      P_initial = P$matrix
    }
    if (i==n_iterations)
    {
      P_final = P$matrix
    }
    
    pb$tick()
  }
  
  # write out diagnostics
  params_to_output = c(current_conditions, 
                       minsize=minsize,
                       maxsize=maxsize,
                       n_iterations=n_iterations,
                       size_mean_initial=size_mean_initial,
                       size_sd_initial=size_sd_initial,
                       n_S_initial = n_S_initial,
                       population_density_initial=population_density_initial
  )
  params_to_output$population_density_m2 = NULL
  
  
  df_n_density = data.frame(time=1:nrow(n_T_time_series), value=apply(n_T_time_series, 1, sum),type='T') %>%
    rbind(data.frame(time=1:length(n_S_time_series),value=n_S_time_series,type='S'))
  
  df_nt = n_T_time_series %>% 
    as.data.frame
  names(df_nt) = round(P$meshpts,3)
  df_nt = df_nt %>%
    mutate(time=row_number()) %>% 
    pivot_longer(!time) %>%
    mutate(size=as.numeric(name))
  
  g1 = ggplot(df_nt, aes(x=size,y=value,color=time,group=time)) +
    geom_line() +
    theme_bw() + 
    scale_x_log10() +
    xlab('DBH of tagged trees (cm)') +
    ylab('Population density (# m^-2)') +
    scale_color_viridis_c() +
    theme(legend.position = 'bottom')
  
  g2 = ggplot(df_n_density, aes(x=time,y=value,color=type)) +
    geom_line() +
    theme_bw() +
    ylab('Population density (# m^-2)') +
    scale_color_discrete(labels=c(S='S (medium saplings)',T='T (tagged trees)')) +
    xlab('Timestep') +
    facet_wrap(~type,scales='free_y')
    theme(legend.position = 'bottom')
  
  g3 = ggplot(n_S_component_time_series %>%
                  mutate(time=row_number()) %>%
                  pivot_longer(!time), aes(x=time,y=value,color=name)) +
    geom_line() +
    theme_bw() +
    ylab('Contribution to n_S')
  
  t_this = paste(names(params_to_output), "=", as.character(params_to_output),collapse = '\n',sep='')
  g4 = ggplot() +
    theme_void() +
    annotate('text', x=0, y=0, label=t_this,size=3)
  
  g_pi = plot_kernel(minsize, maxsize, m=m, k=P_initial,'P initial')
  g_pf = plot_kernel(minsize, maxsize, m=m, k=P_final,'P final')
  
  g_model = ggarrange(g1, g2, g3, g4, g_pi, g_pf, nrow=3, ncol=2)
  ggsave(g_model, file=sprintf('output_figures/m_time_series_%s.pdf',prefix_this),width=12,height=12)
  
  n_T_final = df_n_density %>%
    filter(time==n_iterations & type=='T') %>%
    pull(value)  
  n_S_final = df_n_density %>%
    filter(time==n_iterations & type=='S') %>%
    pull(value)
  
  return(c(params_to_output, 
           # add final abundances
           n_S_final=n_S_final, 
           n_T_final=n_T_final, 
           # add key parameters (hack)
           ns_factor_density=ns_factor_density, 
           ns_factor_survival=ns_factor_survival))
}

# test run
run_model(geneticSexIDmale = '0',
          Ploidy_leveltriploid = '1',
          n_S_initial = 1,
          Cos.aspect = 0,
          prefix_this = 'test', 
          ns_factor_density = 1,
          ns_factor_survival = 1,
          n_iterations = 500)


params = expand.grid(geneticSexIDmale=c("0","1"), 
                     Ploidy_leveltriploid=c("0","1"), 
                     n_S_initial=c(1),
                     Cos.aspect=c(-1,1),
                     Elevation=c(2800,3200),
                     #SWE.0=c(0,300),
                     #ns_factor_density=c(0,1),
                     n_iterations=500,
                     m=100
              )
params$prefix_this=paste('test',1:nrow(params),sep='_')

model_results = lapply(1:nrow(params), function(i) {
  cat(sprintf('%d %.3f\n',i, i/nrow(params)))
  
  result = do.call("run_model",as.list(params[i,]))
  return(as.data.frame(result))
})
model_results = do.call("rbind",model_results)

write.csv(model_results, file='output_data/ipm_outcomes_parameter_sweep.csv', row.names = FALSE)

# try some plotting

ggplot(model_results, aes(x=geneticSexIDmale, color=Ploidy_leveltriploid,y=n_T_final)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~Cos.aspect+Elevation,labeller = label_both) +
  ylab('Tagged tree density (# per m2)')

ggplot(model_results, aes(x=geneticSexIDmale, color=Ploidy_leveltriploid,y=n_S_final)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~Cos.aspect+Elevation,labeller = label_both) +
  ylab('Sapling density (# per 3 m subplot)')

ggplot(model_results, aes(x=n_S_final,y=n_T_final)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~Cos.aspect+Elevation,labeller = label_both) +
  ylab('Tagged tree density (# per m2)') +
  xlab('Sapling density (# per 3 m subplot)')




















# OLD BELOW



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

ipm_female_triploid_south = make_ipm_for_site(
  m_survival = m_survival_all,
  m_growth = m_growth_all,
  m_recruit = m_recruit_all,
  other_vars=list(geneticSexIDM="0",
    Ploidy_levelTriploid="1",
    n_medium_trees=5,
    Cos.aspect=-1,
    Elevation=3000)
  )

# this is a density dependent deterministic model, so let's look
# for elasticities close to equilibrium...
K_female_triploid_south = as.matrix(ipm_female_triploid_south$sub_kernels$P_it_300 + ipm_female_triploid_south$sub_kernels$F_it_300)
image(sqrt(elasticity(K_female_triploid_south)))
e_female_triploid_south = elasticity(K_female_triploid_south)
# # code below from https://compadre-db.org/Education/article/sensitivity-and-elasticity-matrices
# w <- eigen(mat)$vectors
# v <- Conj(solve(w))
# senmat <- Re(v[1,] %*% t(w[,1]))
# emat <- (1/(Re(eigen(mat)$values[1]))) * senmat * mat



# get lambda
plot(ipmr::lambda(ipm_female_triploid_south),xlab='time step',ylab='lambda')

# look at trajectory of age distribution
image(ipm_female_triploid_south$pop_state$n_dbh %>% sqrt,xlab='size',ylab='time')

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

k1 = plot_kernel(ipm_female_triploid_south$sub_kernels$P_it_200, 'P sub-kernel')
k2 = plot_kernel(ipm_female_triploid_south$sub_kernels$F_it_200, 'F sub-kernel')
k3 = plot_kernel((e_female_triploid_south)^(1/2), 'Elasticity (square root)')
ggsave(ggarrange(k1, k2, k3, labels='AUTO', align='hv',nrow=2,ncol=2),file='output_figures/g_ipm_kernel_female_triploid_south.pdf',width=7,height=5)
ggsave(ggarrange(k1, k2, k3, labels='AUTO', align='hv',nrow=2,ncol=2),file='output_figures/g_ipm_kernel_female_triploid_south.png',width=7,height=5)




# now look at performance variation across parameters

params = expand.grid(replicate=1:NUM_RESAMPLES,
                     geneticSexIDM=c("0","1"), 
                     Ploidy_levelTriploid=c("0","1"), 
                     n_medium_trees=c(0,5,15),
                     Cos.aspect=c(-1,0,1),
                     Elevation=c(2800,3000,3200))

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
    lambda = mean(tail(as.numeric(ipmr::lambda(ipm_this)),5)) # take the final 5 lambda values and average
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


#
params = read.csv('output_data/ipm_outcomes_lambda.csv')

# summarize this model
library(glmmTMB)
m_lambda_summary = glmmTMB(lambda ~ geneticSexIDM*Ploidy_levelTriploid + n_medium_trees + Cos.aspect + Elevation + (1|replicate), data=params)
tab_model(m_lambda_summary, file='output_figures/table_model_lambda.html')


params_for_plotting = params %>%
  filter(Cos.aspect==-1 & Elevation==3000) %>%
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
  facet_grid(~n_medium_trees) +
  theme_bw() +
  #scale_y_log10() +
  scale_color_manual(values=c('blue','red'),name='Cytotype') +
  xlab('Sex') +
  ylab(expression(paste(lambda)))

ggsave(g_lambda, file='output_figures/g_ipm_lambda.pdf',width=6,height=3)
ggsave(g_lambda, file='output_figures/g_ipm_lambda.png',width=6,height=3)










# now map lambdas for all the plots
df_sites_for_ipm = transitions_all_filtered_joined_no_na %>%
  select(geneticSexID, Ploidy_level, n_medium_trees, Cos.aspect, Elevation,
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
        Cos.aspect=df_sites_for_ipm$Cos.aspect[i],
        Elevation=df_sites_for_ipm$Elevation[i]),
      size_mean = df_sites_for_ipm$size_mean[i],
      size_sd = df_sites_for_ipm$size_sd[i],
      population_density_initial = df_sites_for_ipm$population_density_m2[i]
  ))
  if (!is.null(ipm_this))
  {
    lambda = mean(tail(as.numeric(ipmr::lambda(ipm_this)),5))
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

df_sites_for_ipm_joined = read.csv('output_data/ipm_outcomes_lambda_by_plot.csv')



df_sites_for_ipm_joined_summarized = df_sites_for_ipm_joined %>% 
  filter(site_type=='random') %>%
  group_by(site_code, X.UTM, Y.UTM) %>% 
  summarize(lambda.mean=mean(lambda,na.rm=TRUE), 
            lambda.sd=sd(lambda,na.rm=TRUE), 
            p=ifelse(length(na.omit(lambda)) > 1, t.test(na.omit(lambda), mu=1)$p.value,NA)) %>%
  mutate(p.adj = p.adjust(p, method="BH")) %>%
  mutate(lambda_binned = ifelse(p<0.05, ifelse(lambda.mean < 1, 'decreasing', 'increasing'), 'stable'))

mean(df_sites_for_ipm_joined_summarized$lambda.mean,na.rm=TRUE)
sd(df_sites_for_ipm_joined_summarized$lambda.mean,na.rm=TRUE)

table(df_sites_for_ipm_joined_summarized$lambda_binned)
table(df_sites_for_ipm_joined_summarized$lambda_binned)/nrow(df_sites_for_ipm_joined_summarized)

g_map_lambda_binned = ggplot(df_sites_for_ipm_joined_summarized %>%
                               filter(!is.na(lambda.mean)), 
                           aes(x=X.UTM,y=Y.UTM, color=lambda.mean,shape=lambda_binned)) + 
  geom_point(alpha=0.8) +
  scale_color_gradient2(midpoint=1,
                       low = 'darkorange',high='darkorchid1',mid = 'gray',
                       limits=c(0.9,1.1),
                        #limits=c(0.95,1.05),
                        name=expression(paste(lambda, " (mean)"))) +
  scale_shape_manual(name='Inference',values=c(4,3,16)) +
  theme_bw() + 
  coord_equal() +
  xlab('Easting (m)') + ylab('Northing (m)')

g_lambda_hist = ggplot(df_sites_for_ipm_joined_summarized, aes(x=lambda.mean)) +
  geom_histogram(binwidth = 0.005,fill='#333333') +
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
