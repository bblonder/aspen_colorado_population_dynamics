library(dplyr)

# get site-level info
data_climate = read.csv('output_data/df_lagged_climate.csv')
df_site_level = read.csv('data/aspen_data_site-level_2018-2023_v1allAspen_2024-11-27.csv')
data_sex = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_aug_11_2021.csv') %>%
  mutate(site_code = Site_Code) %>%
  dplyr::select(-Site_Code,-X.UTM,-Y.UTM) %>%
  mutate(geneticSexID=ifelse(is.na(geneticSexID),'unknown',geneticSexID))
data_site_info = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv') %>%
  dplyr::select(site_code=Site_Code,Cos.aspect,Elevation,Slope,X.UTM,Y.UTM,Watershed)
data_ploidy_new = read.csv('/Users/benjaminblonder/Documents/berkeley/roxy cruz postdoc/aspen greenhouse/data analysis - best version with no exclusion bugs/00 update ploidy rmbl/data_rmbl_ploidy_updated_2025-12-04.csv') %>%
  dplyr::rename(site_code=Site_Code)
df_site_level = df_site_level %>%
  left_join(data_sex, by='site_code') %>%
  mutate(geneticSexID=as.character(factor(geneticSexID, levels=c('M','F'), labels = c('male','female')))) %>%
  mutate(geneticSexID=ifelse(is.na(geneticSexID),'unknown',geneticSexID)) %>%
  left_join(data_site_info, by='site_code') %>%
  left_join(data_ploidy_new, by='site_code') %>%
  mutate(Point_Type = ifelse(nchar(site_code)==4,'Grid','Random'))





df_tree_level = read.csv('data/aspen_data_tree-level_2018-2023_2024-11-22.csv')
# keep only the aspens
df_tree_level = df_tree_level %>% 
  filter(species=="P. tremuloides")
# this is a hack to fix a typo in the raw data (12/9/2024)
df_tree_level$tree_num[df_tree_level$site_code=="PHRAZ01" & df_tree_level$year==2020 & df_tree_level$tree_num==11] = c(11,12)



# # augment the tree level data with the mini-census data from 2023 (code commented out - this is already being done)
# df_to_add = df_site_level %>% filter(year==2023 & is.na(plot_radius_m)) %>%
#   mutate(DBH=dbh_center_live) %>%
#   mutate(species='P. tremuloides') %>%
#   mutate(tree_num=0, found='x',dist=0,dir=NA,
#          tree_status_dead=FALSE, tree_status_damaged=FALSE) %>% #### CHECK THIS WITH ERIN!!!
#   select(names(.)[names(.) %in% names(df_tree_level)])
# 
# # add in missing columns
# names_missing = setdiff(names(df_tree_level),names(df_to_add))
# df_blank = matrix(NA, nrow=nrow(df_to_add),ncol=length(names_missing), dimnames=list(NULL, names_missing))
# df_to_add_full = cbind(df_blank, df_to_add) %>%
#   select(names(df_tree_level))
# 
# # add in the new observations
# df_tree_level = rbind(df_tree_level, df_to_add_full)





# df_tree_level %>% dplyr::select(site_code, year, tree_num, DBH) %>%
#   pivot_wider(id_cols=c(site_code, tree_num), names_from=year, values_from=DBH,values_fill=NA) %>%
#   arrange(site_code) %>%
#   View


add_transition <- function(year_start, year_end, site_code_this)
{
  stopifnot(year_start < year_end)
  
  # get unique data from each year
  df_start = df_tree_level %>% 
    filter(year==year_start & site_code==site_code_this)
  df_end = df_tree_level %>% 
    filter(year==year_end & site_code==site_code_this)
  
  # identify what tree ids we have
  tree_nums = union(df_start$tree_num, df_end$tree_num)
  # print(tree_nums)
  # 
  # print(nrow(df_start))
  # print(nrow(df_end))
  
  df_output = NULL
  for (tree_num_this in tree_nums)
  {
    print(sprintf('%s - %s', site_code_this, tree_num_this))
    
    df_start_tree = df_start %>% 
      filter(tree_num==tree_num_this)
    df_end_tree = df_end %>% 
      filter(tree_num==tree_num_this)

    # if the tree is in both censuses
    if (nrow(df_start_tree) > 0 & nrow(df_end_tree) > 0)
    {
      # check if died
      died_start = df_start_tree$tree_status_dead
      damaged_start = df_start_tree$tree_status_damaged
      if (died_start)
      {
        # don't do anything
        message('tree dropped, already dead in start census')
        temp = NULL       
      }
      else
      {
        message('adding regular tree')
        died_end = df_end_tree$tree_status_dead
        damaged_end = df_end_tree$tree_status_damaged
        # add the transition
        temp = data.frame(size=df_start_tree$DBH,
                          sizeNext=df_end_tree$DBH,
                          surv=!died_end,
                          fec=0, # no sexual reproduction
                          stage='continuous',
                          stageNext='continuous',
                          number=1,
                          damagedNext=damaged_end & !damaged_start,
                          delta_years = year_end - year_start,
                          year = year_start,
                          site_code=site_code_this)
      }
    }
    # if the tree is in the first census but not the later one
    else if (nrow(df_start_tree) > 0 & nrow(df_end_tree) == 0)
    {
      # don't do anything
      message('tree dropped, only in start census, no action')
      temp = NULL
    }
    # if the tree is the later census but not the first one
    else if (nrow(df_start_tree) == 0 & nrow(df_end_tree) > 0)
    {
      # if it is a 'new' tree because it is small-ish
      if (df_end_tree$DBH < 6)
      {
        message('adding new sapling')
        # add the transition
        temp = data.frame(size=5,
                          sizeNext=df_end_tree$DBH,
                          surv=!died_end,
                          fec=0, # no sexual reproduction
                          stage='new_5_cm',
                          stageNext='continuous',
                          number=1,
                          damagedNext=damaged_end,
                          delta_years = year_end - year_start,
                          year = year_start,
                          site_code=site_code_this)    
      }
      else
      {
        # don't do anything
        message('tree added, only in end census and too big to be sapling, no action')
        temp = NULL
      }
    }
    df_output = rbind(df_output, temp)
  }
  
  return(df_output)
}

transitions_all = NULL

site_codes = sort(unique(df_tree_level$site_code))
for (site_code_this in site_codes)
{
  years = df_tree_level %>% 
    filter(site_code==site_code_this) %>%
    pull(year) %>%
    unique %>%
    sort
  
  print(years)
  
  if (length(years) > 1)
  {
    for (i in 1:(length(years)-1))
    {
        transitions_all = rbind(transitions_all,
                                add_transition(years[i], years[i+1], site_code_this))
        message('.\n')
    }
  }
  else
  {
    message(sprintf('*** only one year of data for site %s, skipping all trees',site_code_this))
  }
}
transitions_all = transitions_all %>%
  mutate(growth_rate = (sizeNext - size)/delta_years)

hist(transitions_all$growth_rate,breaks=100)

transitions_all$delta_years %>% table
transitions_all$stage %>% table
transitions_all$damaged %>% table
transitions_all$surv %>% table


# remove dodgy size transitions that are too big
transitions_all_filtered = transitions_all %>%
  mutate(sizeNext = ifelse(abs(growth_rate) < 1,sizeNext,NA)) %>%
  mutate(growth_rate = ifelse(abs(growth_rate) < 1,growth_rate,NA)) %>%
# identify recruits as one that started at exactly 5 and did not grow much 
# (maybe there is a better way, could flag these 
  mutate(recruit = (size==5 & sizeNext < 6))

# add in the small/medium/dead/density/genetics info
transitions_all_filtered_joined = transitions_all_filtered %>%
  left_join(df_site_level %>% 
#              mutate(basal_area_density_live = basal_area_density_live / 4) # to fix a bug in the original BAD code from erin
              select(site_code, year,
                     Ploidy_level, geneticSexID, 
                     n_small_trees, n_medium_trees, 
                     Cos.aspect, Elevation,
                     basal_area_density_live, plot_area_m2), 
            by=join_by(site_code,year)) %>%
  mutate(year_final = year + delta_years) %>%
  select(year, year_final, delta_years, everything())

# now add in the climate data merging on the final year 
# (not sure how to best handle this transition)
transitions_all_filtered_joined_with_climate_year_final = transitions_all_filtered_joined %>% 
  left_join(data_climate %>% 
              rename(year_final=year), 
            by=c('site_code','year_final'))


write.csv(transitions_all_filtered_joined_with_climate_year_final, file='output_data/transitions_all_filtered_joined_with_climate_year_final.csv',row.names=FALSE)

