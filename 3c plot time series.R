library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)

ipms_rmbl = readRDS('output_data/ipms_rmbl_sites.Rdata')

df_rmbl = read.csv('output_data/ipm_rmbl_df.csv') %>%
  filter(replicate_data_resample==1) %>%
  mutate(row=row_number())

set.seed(1)
sites_random = sample(unique(df_rmbl$site_code), 9)

data_bad = do.call('rbind',lapply(sites_random, function(site_code_this)
{
  rows_this = df_rmbl %>% 
    filter(site_code==site_code_this) %>%
    pull(row)
    
  result = NULL
  for (i in rows_this) #length(ipms_rmbl)
  {
    for (j in 1:10)
    {
      df_this = data.frame(site_code=site_code_this,
                           climate_replicate=i, 
                           dataset_replicate=j, 
                           bad=ipms_rmbl[[i]][[j]]$bad) %>%
        mutate(time=row_number())
      result = rbind(result, df_this)
    }
  }
  return(result)
}))

g_bad = ggplot(data_bad, aes(x=time,
                       y=bad,
                       color=factor(dataset_replicate),
                       group=paste(climate_replicate, dataset_replicate))) + 
  facet_wrap(~site_code,scales='free_y',labeller = label_both) +
  geom_line() +
  theme_bw() +
  #theme(legend.position='none') +
  xlab('Timestep') + 
  ylab('Basal area density (m2 m-2)') +
  scale_y_sqrt() +
  scale_color_discrete(name='Dataset resample') +
  theme(legend.position='bottom')
ggsave(g_bad, file='output_figures/g_example_time_series_bad.pdf', width=7,height=7)
ggsave(g_bad, file='output_figures/g_example_time_series_bad.png', width=7,height=7)





size_bins = ipms_rmbl[[1]][[1]]$size_bins

data_n_A = do.call('rbind',lapply(sites_random, function(site_code_this)
{
  rows_this = df_rmbl %>% 
    filter(site_code==site_code_this) %>%
    pull(row)
  
  result = NULL
  for (i in rows_this) #length(ipms_rmbl)
  {
    for (j in 1:10)
    {
      print(i)
      print(j)
      n_A_this = ipms_rmbl[[i]][[j]]$n_A %>% 
        as.data.frame
      names(n_A_this) = size_bins
      
      n_A_this = n_A_this %>% 
        mutate(time=row_number()) 
      
      
      n_A_this = n_A_this %>% 
        pivot_longer(!time)
      
      df_this = data.frame(site_code=site_code_this,
                           climate_replicate=i, 
                           dataset_replicate=j, 
                           n_A=n_A_this)
      result = rbind(result, df_this)
    }
  }
  return(result)
}))

str(data_n_A)

g_n_A = ggplot(data_n_A %>% filter(n_A.time %in% seq(1,300,by=100)), 
       aes(x=as.numeric(n_A.name),
           y=n_A.value,group=paste(climate_replicate,dataset_replicate, n_A.time),
           linetype=factor(dataset_replicate),
           color=n_A.time)) +
  geom_line() +
  facet_wrap(~site_code,scales='free_y',labeller = label_both) +
  scale_y_sqrt() +
  theme_bw() +
  scale_color_viridis_c(option='D',name='Timestep') +
  xlab('Size (cm)') + ylab(expression(paste('n'['A']))) +
  theme(legend.position='bottom') +
  scale_linetype(name='Dataset resample')

ggsave(g_n_A, file='output_figures/g_example_time_series_n_A.pdf', width=7,height=7)
ggsave(g_n_A, file='output_figures/g_example_time_series_n_A.png', width=7,height=7)

