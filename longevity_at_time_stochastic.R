#P_this = readRDS('output_figures/model_outputs/P_400_test.Rdata')

longevity_at_time_stochastic <- function(P_list, lx_crit=0.1, n_A_initial=0.1, m=100, L=4, U=60)
{
  size_bins = seq(L,U,length.out=m)
  n_A = rep(0, length(size_bins))
  n_A[which(size_bins > 5)[1]] = n_A_initial # all saplings
  #plot(size_bins, n_A)
  
  i_max = 1000#
  i = 1
  
  survivorship = rep(NA, length(P_list))
  
  while (i < i_max)
  {
    if (i > length(P_list))
    {
      i_effective = (i %% length(P_list)) + 1
    }
    else
    {
      i_effective = i
    }
    n_A = P_list[[i_effective]] %*% n_A # loop over the P matrix if we run too far
    total_density = sum(n_A)
    survivorship[i] = (total_density/n_A_initial)
    
    i = i+1
  }
  plot(survivorship,xlab='time',ylab='survivorship',type='l',col='red')
  longevity_crit = which(survivorship < lx_crit)[1] # first time we get below threshold
  if (length(longevity_crit)==0)
  {
    longevity_crit = NA
  }
  return(longevity_crit)
  #return(n_A)
}

#longevity_at_time_stochastic(rep(list(P_this), 1000), n_A_initial = 0.1, L=4, U=60)
