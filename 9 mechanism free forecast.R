PREFIX_TYPE = 'v2closest11'
source('0 analyze 2018-2013.R')

make_delta <- function(df, year1, year2, response)
{
  stopifnot(year1 < year2)
  
  df_year1 = df %>% 
    filter(year %in% year1) %>%
    mutate(stem_density_live = n_standing / plot_area_m2)
  
  df_year2 = df %>% 
    filter(year %in% year2) %>%
    mutate(stem_density_live = n_standing / plot_area_m2)

  
  df_response = df_year1 %>%
    select(site_code, !!response) %>%
    left_join(df_year2 %>% select(site_code, !!response), by='site_code')
  
  # year1 before year2
  names(df_response) = c("site_code",paste(response, "prev", sep='_'),response)
    
  #df_response$delta_response = df_response[,3] - df_response[,2]
  df_response$delta_t = year2 - year1
  df_response$year_prev = year1
  df_response = df_response[,-2] # remove the year1 ('prev') value because it will be joined in below in a more general way
  
  df_predictors_year1 = df_year1 %>%
    select(n_small_trees, 
           n_medium_trees, 
           dbh_center_live, 
           n_adult_dead_w_background_mortality_estimate_per_m2, 
           n_adult_damaged_per_m2, 
           basal_area_density_live, 
           stem_density_live)

  names(df_predictors_year1) = paste(names(df_predictors_year1),"prev",sep='_')
  df_predictors_year1$site_code = df_year1$site_code # only works if we don't filter() in this block
  
  df_static = df_year1 %>%
    select(site_code,
           Elevation,
           Cos.aspect,
           Slope,
           Ploidy_level,
           geneticSexID,
           X.UTM, Y.UTM)
  
  df_response_joined = df_response %>% 
    left_join(df_predictors_year1, by='site_code') %>%
    left_join(df_static, by='site_code')
  
  return(df_response_joined)
}

make_demography_data <- function(response)
{
  rbind(
    make_delta(df_for_regression %>% filter(outlier_growth==FALSE), 2018, 2020, response),
    make_delta(df_for_regression %>% filter(outlier_growth==FALSE), 2020, 2023, response)
  )
}

fit_demography_model <- function(response, data, family, zi)
{
  m = glmmTMB(
    formula= formula(sprintf("%s ~ delta_t * (
        n_small_trees_prev +
        n_medium_trees_prev + 
        dbh_center_live_prev + 
        n_adult_dead_w_background_mortality_estimate_per_m2_prev +
        n_adult_damaged_per_m2_prev +
        basal_area_density_live_prev +
        stem_density_live_prev + 
        Elevation +
        Cos.aspect +
        Slope +
        Ploidy_level * geneticSexID)",response)),
    data=data,
    family=family,
    ziformula = formula(ifelse(zi==TRUE, "~1", "~0"))) 

  return(m)
}

# 



m_small_trees = fit_demography_model(response="n_small_trees",
                                     data=make_demography_data("n_small_trees"),
                                     family=lognormal,zi=TRUE)

m_medium_trees = fit_demography_model(response="n_medium_trees",
                                     data=make_demography_data("n_medium_trees"),
                                     family=lognormal,zi=TRUE)

m_dbh_center_live = fit_demography_model(response="dbh_center_live",
                                      data=make_demography_data("dbh_center_live"),
                                      family=gaussian,zi=FALSE)
m_n_adult_dead_w_background_mortality_estimate_per_m2 = fit_demography_model(response="n_adult_dead_w_background_mortality_estimate_per_m2",
                                      data=make_demography_data("n_adult_dead_w_background_mortality_estimate_per_m2"),
                                      family=lognormal,zi=TRUE)
m_n_adult_damaged_per_m2 = fit_demography_model(response="frac_adult_damaged_per_m2",
                                                                                data=make_demography_data("frac_adult_damaged_per_m2"),
                                                                                family=lognormal,zi=TRUE)
m_basal_area_density_live = fit_demography_model(response="basal_area_density_live",
                                                data=make_demography_data("basal_area_density_live"),
                                                family=lognormal,zi=FALSE)
m_stem_density_live = fit_demography_model(response="stem_density_live",
                                                 data=make_demography_data("stem_density_live"),
                                                 family=lognormal,zi=FALSE)







# initialize the current state as the 2018 data
# this is a hack to get the key predictors in the right tabular form/names, with a 1-year delta_t and 2023 as starting year
data_current = make_delta(data_site %>% filter(Point_Type=='Random'), 2018, 2020, "cow") %>% 
  select(-cow) %>% 
  na.omit

# store a copy of the state
predictions_demography = data_current
for (i in 1:100)
{
  print(i)
  # do an iteration loop advancing all the state variables
  n_small_trees_current = predict(m_small_trees, newdata=data_current, type='response')
  n_medium_trees_current = predict(m_medium_trees, newdata=data_current, type='response')
  dbh_center_live_current = predict(m_dbh_center_live, newdata=data_current, type='response')
  n_adult_dead_w_background_mortality_estimate_per_m2_current = predict(m_n_adult_dead_w_background_mortality_estimate_per_m2, newdata=data_current, type='response')
  n_adult_damaged_per_m2_current = predict(m_n_adult_damaged_per_m2, newdata=data_current, type='response')
  basal_area_density_live_current = predict(m_basal_area_density_live, newdata=data_current, type='response')
  stem_density_live_current = predict(m_stem_density_live, newdata=data_current, type='response')
  
  data_predicted = data_current %>%
    # advance the year
    mutate(year_prev = year_prev+1) %>%
    # advance the demography variables based on predictions
    mutate(n_small_trees_prev = n_small_trees_current,
           n_medium_trees_prev = n_medium_trees_current,
           dbh_center_live_prev = dbh_center_live_current,
           n_adult_dead_w_background_mortality_estimate_per_m2_prev = n_adult_dead_w_background_mortality_estimate_per_m2_current,
           n_adult_damaged_per_m2_prev = n_adult_damaged_per_m2_current,
           basal_area_density_live_prev = basal_area_density_live_current, 
           stem_density_live_prev = stem_density_live_current)
  
  # # cap predictions that are too large to try and regularize the model
  # data_predicted = data_predicted %>%
  #   mutate(n_small_trees_prev = min(100,n_small_trees_prev),
  #          n_medium_trees_prev = min(50,n_medium_trees_prev),
  #          dbh_center_live_prev = min(80,dbh_center_live_prev),
  #          n_adult_dead_w_background_mortality_estimate_per_m2_prev = min(1,n_adult_dead_w_background_mortality_estimate_per_m2_prev),
  #          n_adult_damaged_per_m2_prev = min(1,n_adult_damaged_per_m2_prev),
  #          basal_area_density_live_prev = min(0.2,basal_area_density_live_prev),
  #          stem_density_live_prev = min(10,min(stem_density_live_prev)))
  # copy the result
  predictions_demography = rbind(predictions_demography, data_predicted)
  
  # replace the current state with the predicted state
  data_current = data_predicted
}
# calculate overall population size
predictions_demography = predictions_demography %>%
  mutate(population_density_m2 = n_small_trees_prev / (pi*9^2) + n_medium_trees_prev / (pi*9^2) + stem_density_live_prev)



int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}

plot_predicted_demography <- function(data, yvar, ymax)
{
  ggplot(data, aes(x=year_prev,y=.data[[yvar]],group=site_code,color=Ploidy_level)) +
    geom_line() +
    geom_point() + 
    theme_bw() +
    facet_wrap(~Ploidy_level+geneticSexID) +
    ggtitle(yvar) +
    scale_color_manual(values=c('blue','red')) +
    scale_x_continuous(breaks = int_breaks) +
    scale_y_sqrt(limits=c(0,ymax))
}

plot_predicted_demography(predictions_demography, "n_small_trees_prev", 100)
plot_predicted_demography(predictions_demography, "n_medium_trees_prev", 40)
plot_predicted_demography(predictions_demography, "dbh_center_live_prev", 60)
plot_predicted_demography(predictions_demography, "n_adult_dead_w_background_mortality_estimate_per_m2_prev", 0.5)
plot_predicted_demography(predictions_demography, "n_adult_damaged_per_m2_prev", 0.2)
plot_predicted_demography(predictions_demography, "basal_area_density_live_prev", 0.1)
plot_predicted_demography(predictions_demography, "stem_density_live_prev", 1)
g_density = plot_predicted_demography(predictions_demography, "population_density_m2", 1)
ggsave(g_density, file='figures/g_black_box_population_density.pdf',width=10,height=10)

predictions_demography %>%
  #filter(population_density_m2 < 10000) %>%
  slice_max(year_prev) %>%
  ggplot(aes(x=X.UTM,y=Y.UTM,color=sqrt(population_density_m2))) +
  geom_point() +
  scale_color_viridis_c() +
  theme_bw()


# is the issue something with how we are dealing with absolute abundance vs fractional abundance vs abundance density?
# are there really no sites that are going extinct? check some examples like the ones with high 2018 mortality
# maybe we need to restrict our data to within the original fixed radius to ensure we are not biased against dead trees
# should we restrict the analysis to only the full demograhpy sites and get more information on 

# i think the standing alive density appears to increase because of how we have defined the variable radius plots. maybe we need
# to have a sampling model because when the plot radius is just to the next tree we bias towards increased density
# could make a plot of tree number density vs distance from the full census to understand this... add trees in increasing order?
# maybe we actually need to do the analysis on a fractional basis 

ggplot(df_for_regression, aes(x=year,y=I(n_standing/plot_area_m2),group=site_code)) + geom_line() + geom_point()
