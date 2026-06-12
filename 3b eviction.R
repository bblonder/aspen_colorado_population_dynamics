coef_survival = model_output_test$coef_survival
coef_size_variance = model_output_test$coef_size_variance
coef_growth = model_output_test$coef_growth
P_final = readRDS('output_figures/model_outputs/P_test.Rdata')

sx <- function(x) {
  if (x>50)
  {
    mu = 0.95 # 5% mortality rate, so 1/0.05*2.5 = 50 year expected time to death once this size class reached
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
    xbeta_sd <- coef_size_variance[1] + x*coef_size_variance[2]
    mu_sd = 1/xbeta_sd # gamma inv link
    sigmax = sqrt(mu_sd);
    
    mux<-coef_growth[1] + x*coef_growth[2]
    fac1<-sqrt(2*pi)*sigmax;
    fac2<-((y-mux)^2)/(2*sigmax^2);
    return(exp(-fac2)/fac1);
  }
}

size_bins = 5:70
epsilon_x = sapply(size_bins, function(x) {1-integrate(function(y) {gxy(x,y)}, head(size_bins,1), tail(size_bins,1))$value})

pdf(file='output_figures/g_eviction.pdf')
plot_kernel(5, 70, 
            100, P_final, 'P after 300 timesteps')
plot(size_bins, epsilon_x,type='l',main='flattening gxy outside 5-60 cm',ylab='epsilon(x)',ylim=c(-1,1),color='red')
plot(model_output_test$bad, main='basal area density',xlab='timestep')
dev.off()


       