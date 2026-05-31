ipms_rmbl_sites_dry_scenario = lapply(1:nrow(df_sites_for_ipm), function(i)
{
  cat(sprintf('%d %.3f\n',i, i/nrow(df_sites_for_ipm)))
  params_this = as.list(df_sites_for_ipm[i,])
  # make environment sequence
  climate_this = climate_all[[i]]
  # make the climate more stressful
  params_this$SWE.sequence=climate_this$SWE * 0.5
  params_this$STB.sequence=climate_this$STB - 1
  # remove unneeded columns
  params_this$year = NULL
  params_this$site_type = NULL
  params_this$site_code = NULL
  params_this$Latitude = NULL
  params_this$Longitude = NULL
  # change run name
  params_this$prefix_this = gsub('rmbl','rmbl_dry_scenario', params_this$prefix_this)
  #print(params_this)
  
  result = do.call("run_model_resampled",params_this)
  
  return(result)
})



ipm_rmbl_dry_scenario_df = do.call('rbind',lapply(1:length(ipms_rmbl_sites_dry_scenario), function(r) { do.call('rbind',lapply(ipms_rmbl_sites_dry_scenario[[r]], function(x) { as.data.frame(x$results) %>% mutate(index=r) })) })) %>%
  mutate(site_code=sapply(strsplit(prefix_this,split='_'),function(x) {x[2]})) %>%
  mutate(replicate_climate_resample=sapply(strsplit(prefix_this,split='_'),function(x) {x[3]})) %>%
  mutate(replicate_data_resample=sapply(strsplit(prefix_this,split='_'),function(x) {x[4]})) %>%
  left_join(df_sites,by='site_code') %>% 
  mutate(Ploidy_level=ifelse(Ploidy_leveltriploid=='1','triploid',ifelse(Ploidy_levelunknown=='1','unknown','diploid'))) %>%
  mutate(geneticSexID=ifelse(geneticSexIDmale=='1','male',ifelse(geneticSexIDunknown=='1','unknown','female')))

write.csv(ipm_rmbl_dry_scenario_df, file='output_data/ipm_rmbl_dry_scenario_df.csv',row.names = FALSE)

