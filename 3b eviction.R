library(progress)
library(dplyr)
library(tidyr)
library(Rage)
library(ggpubr)
library(ggplot2)
# assumes script 3 has been run up to creating model test
load('output_data/workspace for ipm.Rdata')


source('get_climate.R')
set.seed(4)

climate_test = make_climate_ts_at_location(lat=38.959158,lon=-106.9897676, num_time_points = 1000, years_this = 2010:2026) # Gothic, CO
climate_test$SWE.lagged=make_lagged(climate_test$SWE, 6)
climate_test$PPT.lagged=make_lagged(climate_test$PPT, 6)


# test run
model_output_test = run_model(Ploidy_leveltriploid='0',
                              Ploidy_levelunknown='0',
                              geneticSexIDmale='0',
                              geneticSexIDunknown='0',
                              Cos.aspect = -1,
                              Elevation=2800,
                              prefix_this = 'test', 
                              n_iterations = 500,
                              t_average = 300,
                              n_S_scale_factor=1,
                              survival_fraction_max = 0.99,
                              SWE.lagged=climate_test$SWE.lagged,#rep(200, 1000),
                              PPT.lagged=climate_test$PPT.lagged,#rep(100, 1000),
                              save_plot = TRUE, save_P = TRUE)
model_output_test$results$longevity_90_final_average
model_output_test$results$basal_area_density_final_average


coef_survival = model_output_test$coef_survival
coef_sizenext_size_variance = model_output_test$coef_sizenext_size_variance
coef_sizenext = model_output_test$coef_sizenext
P_final = readRDS(sprintf('output_figures/model_outputs/P_%d_test.Rdata',NUM_ITERATIONS))

# survival
sx <- function(x) {
  xbeta = coef_survival[1] + x*coef_survival[2]
  mu = exp(xbeta)/(1 + exp(xbeta)) #
  
  mu = mu * survival_fraction_max
  
  # # apply correction for larger individuals
  if (x>50)
  {
    mu = mu*0
  }
  
  #print(sprintf('survival=%.3f'))
  return(mu)
}
# growth
gxy<-function(x,y) {
  if (x>50)
  {
    return(gxy(50,y))
  } 
  else if (x<5)
  {
    return(gxy(5,y))
  } 
  else
  {
    xbeta_sd <- coef_sizenext_size_variance[1] + x*coef_sizenext_size_variance[2]
    mu_sd = 1/xbeta_sd # gamma inv link
    sigmax = sqrt(mu_sd);
    
    mux<-coef_sizenext[1] + x*coef_sizenext[2]
    fac1<-sqrt(2*pi)*sigmax;
    fac2<-((y-mux)^2)/(2*sigmax^2);
    return(exp(-fac2)/fac1);
  }
}

size_bins = seq(4,60,by=0.1)
epsilon_x = sapply(size_bins, function(x) {1-integrate(function(y) {gxy(x,y)}, head(size_bins,1), tail(size_bins,1))$value})

df_eviction = data.frame(size=size_bins, epsilon=epsilon_x)

g_eviction1 = ggplot(df_eviction, aes(x=size,y=epsilon)) + 
  geom_line() + 
  geom_point() +
  theme_bw() +
  ylim(-1e-2,1) +
  ylab(expression(epsilon)) +
  xlab(expression(paste('Size (cm)')))

g_eviction2 = plot_kernel(5, 60, 
                          100, P_final, sprintf('P after %d timesteps', NUM_ITERATIONS))


g_eviction_all = ggarrange(g_eviction1, g_eviction2, labels='auto',align='hv')
ggsave(g_eviction_all, file='output_figures/g_eviction.png',width=7,height=3)
ggsave(g_eviction_all, file='output_figures/g_eviction.pdf',width=7,height=3)


       