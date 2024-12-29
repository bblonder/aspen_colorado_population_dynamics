# classes: small, medium, large5_10, large10_15, large15_20, large20_plus

# estimate probability of new small or new medium as change in number of small per time per number of small, assume no mortality? or of the total, 2/3 new, ?




# make separate transition probabilities for males, females
# 
# 
# df_for_regression_full_surveys = df_for_regression %>% 
#   filter(n_aspen_trees_surveyed==11)


#library(IPMpack)
library(dplyr)

# get site-level info
df_site_level = read.csv('data/aspen_data_site-level_2018-2023_v1allAspen_2024-11-27.csv')
data_sex = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_aug_11_2021.csv') %>%
  mutate(site_code = Site_Code) %>%
  select(-Site_Code,-X.UTM,-Y.UTM)
data_ploidy = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv') %>%
  select(site_code=Site_Code,Ploidy_level,Cos.aspect,Elevation,Slope,X.UTM,Y.UTM,Watershed)
df_site_level = df_site_level %>%
  left_join(data_sex, by='site_code') %>%
  left_join(data_ploidy, by='site_code') %>%
  mutate(Point_Type = ifelse(nchar(site_code)==4,'Grid','Random'))


df_tree_level = read.csv('data/aspen_data_tree-level_2018-2023_2024-11-22.csv')
# keep only the aspens
df_tree_level = df_tree_level %>% 
  filter(species=="P. tremuloides")
# this is a hack to fix a typo in the raw data (12/9/2024)
df_tree_level$tree_num[df_tree_level$site_code=="PHRAZ01" & df_tree_level$year==2020 & df_tree_level$tree_num==11] = c(11,12)



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
                          stage='continuous',
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


# remove dodgy size transitions
transitions_all_filtered = transitions_all %>%
  mutate(sizeNext = ifelse(abs(growth_rate) < 1.5,sizeNext,NA)) %>%
# identify recruits as one that started at exactly 5 and did not grow much 
# (maybe there is a better way, could flag these 
  mutate(recruit = (size==5 & sizeNext < 5.5))

# add in the small/medium/dead/density/genetics info
transitions_all_filtered_joined = transitions_all_filtered %>%
  left_join(df_site_level %>% 
              select(site_code, year, 
                     Ploidy_level, geneticSexID, 
                     n_small_trees, n_medium_trees, 
                     Cos.aspect, Elevation,
                     basal_area_density_live, plot_area_m2), 
            by=join_by(site_code,year)) 


write.csv(transitions_all_filtered_joined, file='output_data/transitions_all_filtered_joined.csv')


  
# 
# 
# growObj = makeGrowthObj(transitions_all_filtered)
# survObj = makeSurvObj(transitions_all_filtered)
# fecObj = makeFecObj(transitions_all_filtered)
# 
# Pmatrix = makeIPMPmatrix(minSize=5, maxSize=100, growObj = growObj, survObj = survObj)
# Fmatrix = makeIPMFmatrix(minSize=5, maxSize=100, fecObj = fecObj)
# 
# fullMatrix = largeMatrixCalc(Pmatrix, Fmatrix)


# # really no recruits??
# df_tree_level %>% filter(DBH < 5.5 & tree_num>10) %>% View


# 
# # figure out which sites only had a 2018 and 2023 census
# years_censused = df_tree_level %>% 
#   group_by(site_code, year) %>% 
#   tally %>% 
#   pivot_wider(id_cols=site_code,names_from=year,values_from=n) %>%
#   dplyr::select(site_code, `2018`,`2020`,`2023`) %>%
#   mutate(only_five = is.na(`2020`) & !is.na(`2023`) & !is.na(`2018`))
# 
# 

# add_transition(2018, 2023, "PHRAZ01")

# add_transition(2018, 2020, "ERNAK")

# only need to do the tagged trees, can use the small/medium to estimate the fecundity