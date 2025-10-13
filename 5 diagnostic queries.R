library(ggplot2)
library(dplyr)
library(glmmTMB)
library(DHARMa)

transitions_all_filtered_joined = read.csv('output_data/transitions_all_filtered_joined.csv')

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
                n_medium_trees, 
                population_density_m2, size_mean, size_sd) %>%
  # estimate a size_sd for the cases we don't have one
  ungroup %>%
  mutate(size_sd = ifelse(is.na(size_sd), mean(size_sd, na.rm=TRUE), size_sd)) %>%
  mutate(Ploidy_level=factor(Ploidy_level), geneticSexID=factor(geneticSexID))


# add weights
counts_by_plot_year = transitions_all_filtered_joined_no_na %>% 
  group_by(year, site_code) %>% 
  tally %>%
  arrange(n) %>%
  mutate(weight = 1/n)

transitions_all_filtered_joined_no_na = transitions_all_filtered_joined_no_na %>%
  left_join(counts_by_plot_year %>% dplyr::select(year, site_code, weight),
            by=c('year','site_code'))

# fit vital rate models
pdf(file='output_figures/residuals_vital_rate.pdf')
formula_survival = formula(factor(surv) ~ delta_years + size * (Ploidy_level + geneticSexID + Cos.aspect + Elevation) + population_density_m2)
m_survival_all = glm(formula=formula_survival, 
                     data = transitions_all_filtered_joined_no_na, 
                     family = binomial(),
                     weights=round(10*transitions_all_filtered_joined_no_na$weight))
plot(m_survival_all)

formula_growth = formula(sizeNext ~ delta_years + I(size^2)*delta_years + size * (Ploidy_level + geneticSexID + Cos.aspect + Elevation) + population_density_m2)
m_growth_all = glm(formula=formula_growth, 
                   #family=quasi(variance="mu^2"),
                   data = transitions_all_filtered_joined_no_na,
                   weights=round(10*transitions_all_filtered_joined_no_na$weight))
plot(m_growth_all)
#plot(resid(m_growth_all)~predict(m_growth_all))
#plot(simulateResiduals(m_growth_all))

formula_recruit = formula(factor(recruit) ~ delta_years + n_medium_trees * (Ploidy_level + geneticSexID + Cos.aspect + Elevation) + population_density_m2)
m_recruit_all = glm(formula=formula_recruit,
                    data=transitions_all_filtered_joined_no_na,
                    family = binomial(),
                    weights=round(10*transitions_all_filtered_joined_no_na$weight))
plot(m_recruit_all)
dev.off()


# look at outliers
transitions_all_filtered_joined_no_na[669,]
# ECGX has a big jump in size from 2018 to 2020 but the annualized rate is OK


pdf(file='output_figures/size_distribution.pdf')
hist(transitions_all_filtered_joined_no_na$size,breaks=100,xlim=c(0,70))
abline(v=5,col='red')
abline(v=60,col='blue')
dev.off()

pdf(file='output_figures/growth.pdf')
plot(transitions_all_filtered_joined_no_na$size, transitions_all_filtered_joined_no_na$sizeNext)
dev.off()





# look at density issue
ipm_outcome_lambda_by_plot = read.csv('output_data/ipm_outcomes_lambda_by_plot.csv')
data_site = read.csv('data/aspen_data_site-level_2018-2023_v2closest11_2024-11-27.csv')
ipm_outcome_lambda_by_plot_joined = ipm_outcome_lambda_by_plot %>%
  left_join(data_site %>% select(site_code, basal_area_density_live, year) %>% filter(year==2018))

pdf(file='output_figures/basal_area_lambda.pdf')
plot(ipm_outcome_lambda_by_plot_joined$basal_area_density_live, ipm_outcome_lambda_by_plot_joined$lambda,
     xlab='basal area density 2018 (m2/m2)', ylab='lambda')
dev.off()