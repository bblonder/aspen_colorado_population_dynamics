library(progress)
library(dplyr)
library(tidyr)
library(Rage)
library(ggpubr)
library(ggplot2)
# assumes script 3 has been run up to creating model test
load('output_data/workspace for ipm.Rdata')


model_output_test = run_model(Ploidy_leveltriploid='0',
                              Ploidy_levelunknown='0',
                              geneticSexIDmale='1',
                              geneticSexIDunknown='0',
                              Cos.aspect = 0,
                              Elevation=3000,
                              prefix_this = 'test', 
                              n_iterations = NUM_ITERATIONS,
                              n_S_scale_factor=1,
                              SWE.sequence=rnormTrunc(n = NUM_ITERATIONS+1, 
                                                      mean = mean(transitions_all_filtered_joined_no_na$SWE.0), 
                                                      sd = sd(transitions_all_filtered_joined_no_na$SWE.0),
                                                      min=0,
                                                      max=1000),
                              Tmax.sequence=rnormTrunc(n = NUM_ITERATIONS+1, 
                                                       mean = mean(transitions_all_filtered_joined_no_na$Tmax.0), 
                                                       sd = sd(transitions_all_filtered_joined_no_na$Tmax.0),
                                                       min=0,
                                                       max=1000),
                              save_plot = TRUE, save_P = TRUE)


coef_survival = model_output_test$coef_survival
coef_sizenext_size_variance = model_output_test$coef_sizenext_size_variance
coef_sizenext = model_output_test$coef_sizenext
P_final = readRDS(sprintf('output_figures/model_outputs/P_%d_test.Rdata',NUM_ITERATIONS))

sx <- function(x) {
  if (x>50)
  {
    mu = 0.90 # 5% mortality rate, so 1/0.05*2.5 = 50 year expected time to death once this size class reached
  }
  else
  {
    xbeta = coef_survival[1] + x*coef_survival[2]
    mu = exp(xbeta)/(1 + exp(xbeta)) #
    return(mu)
  }
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


       