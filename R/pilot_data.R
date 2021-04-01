
# Load packages -----------------------------------------------------------

library(tidyverse)
library(psych)
# library(summarytools)
library(lme4)
library(car)
library(performance)
library(BayesFactor)
library(bayestestR)
# library(brms)  # only needed when running Bayesian analysis
# library(rstan)  # only needed when running Bayesian analysis
# library(papaja)
library(here)

# Source functions --------------------------------------------------------

# source(here("R/functions/glmm_stability.R"))
# source(here("R/functions/boot_glmm.R"))


# Input and process data --------------------------------------------------

# Read data
prereg_data <- read_csv(here("data/md1_data_vienna_long.csv")) %>%
  filter(preregistered_data == "after_prereg") %>% # only data that were collected after preregistration
  na.exclude() %>%
  droplevels()

# summary(prereg_data)
# view(dfSummary(prereg_data))


# Aggregating data
agg_data <- prereg_data %>%
  group_by(subject_ID, sex, breed, age, condition) %>%
  summarise(mean.resp = mean(response)) %>%
  ungroup()

# demographics:
agg_data %>%
  filter(condition == "ost") %>%
  summarise(mean.age = mean(age), min.age = min(age), max.age = max(age), sd.age = sd(age), females = sum(sex == "F"), males = sum(sex == "M"))

summary_data <- agg_data %>%
  group_by(condition) %>%
  summarise(mean.resp2 = mean(mean.resp), median.resp = median(mean.resp), sd = sd(mean.resp), se = sd(mean.resp) / sqrt(length(mean.resp)))

num_subjects <- length(unique(prereg_data$subject_ID))
# summary_data


# One-sample t-test to compare against chance level -----------------------

# Frequentist
tt.ost <- t.test(agg_data$mean.resp[agg_data$condition == "ost"], mu = 0.5, alternative = "two.sided")
tt.non <- t.test(agg_data$mean.resp[agg_data$condition == "non"], mu = 0.5, alternative = "two.sided")
tt.oc <- t.test(agg_data$mean.resp[agg_data$condition == "odour"], mu = 0.5, alternative = "two.sided")
# tt.ost
# tt.non
# tt.oc

# Bayes factors
ttbf.ost <- ttestBF(agg_data$mean.resp[agg_data$condition == "ost"], mu = 0.5, alternative = "two.sided")
ttbf.non <- ttestBF(agg_data$mean.resp[agg_data$condition == "non"], mu = 0.5, alternative = "two.sided")
ttbf.oc <- ttestBF(agg_data$mean.resp[agg_data$condition == "odour"], mu = 0.5, alternative = "two.sided")
# extractBF(ttbf.ost)$bf
# extractBF(ttbf.non)$bf
# extractBF(ttbf.oc)$bf



# Generalized linear mixed model ------------------------------------------

# Centering variables for modeling

model_data <- prereg_data %>%
  filter(condition == "ost" | condition == "non") %>%
  mutate(
    z.age = as.numeric(scale(age, scale = T, center = T)),
    sex.c = as.numeric(scale(as.numeric(as.factor(sex)), scale = F, center = T)),
    condition.c = as.numeric(scale(as.numeric(as.factor(condition)), scale = F, center = T)),
    condition_order.c = as.numeric(scale(as.numeric(as.factor(condition_order)), scale = F, center = T)),
    z.trial_num = as.numeric(scale(trial_num, scale = T, center = T)),
    z.training_experience = as.numeric(scale(dog_training_experience, scale = T, center = T))
  )

# view(dfSummary(model_data))


contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))

full_model <- glmer(response ~ condition + condition_order + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID),
  family = binomial, data = model_data, control = contr
)


#* Likelihood ratio test*
drop1_full_model <- drop1(full_model, test = "Chisq", control = contr)
# round(drop1_full_model, 3)


#* Estimates of the fixed effects*
# round(summary(full_model)$coefficients, 3)


#* Random effects*
summary(full_model)$varcor


## Bootstrap 95% CIs ------------------------------------------------------
# #*Calculating confidence intervals with 1000 bootstraps*
# #perform bootstraps for all predictors
# full_model_boot_ci <- boot.glmm.pred(model.res = full_model, excl.warnings = T, nboots = 1000, para = T, n.cores="all-1", level = 0.95)
# write_csv(full_model_boot_ci$ci.estimates, here("R/full_model_boot_ci.csv"))
full_model_boot_ci <- read_csv(here("R/full_model_boot_ci.csv"))

# round(full_model_boot_ci$ci.estimates, 3)


## Bayes factors ----------------------------------------------------------

# ncores = parallel::detectCores()
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# 
# bayes_full <- brm(response ~ condition + condition_order + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)
# 
# bayes_nocondition <- brm(response ~ condition_order + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)
# 
# (condition_bf <- bayes_factor(bayes_full, bayes_nocondition, repetitions = 10, silent = TRUE))
# 
# bayes_noconditionorder <- brm(response ~ condition + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)
# 
# (order_bf <- bayes_factor(bayes_full, bayes_noconditionorder, repetitions = 10, silent = TRUE))
# 
# bayes_notrialnum <- brm(response ~ condition + condition_order + sex + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)
# 
# (trialnum_bf <- bayes_factor(bayes_full, bayes_notrialnum, repetitions = 10, silent = TRUE))
# 
# bayes_nosex <- brm(response ~ condition + condition_order + z.trial_num + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)
# 
# (sex_bf <- bayes_factor(bayes_full, bayes_nosex, repetitions = 10, silent = TRUE))
# 
# bayes_noage <- brm(response ~ condition + condition_order + z.trial_num + sex + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)
# 
# (age_bf <- bayes_factor(bayes_full, bayes_noage, repetitions = 10, silent = TRUE))
# 
# bayes_notraining <- brm(response ~ condition + condition_order + z.trial_num + sex + z.age + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)
# 
# (training_bf <- bayes_factor(bayes_full, bayes_notraining, repetitions = 10, silent = TRUE))
# 
# model_bfs <- tibble(effect = c("(Intercept)", "Condition", "Order of condition", "Trial number", "Sex", "Age", "C-BARQ trainability score"), bf = c(NA, condition_bf$bf_median_based, order_bf$bf_median_based, trialnum_bf$bf_median_based, sex_bf$bf_median_based, age_bf$bf_median_based, training_bf$bf_median_based))
# model_bfs
# write_csv(model_bfs, here("R/model_bfs.csv"))
model_bfs <- read_csv(here("R/model_bfs.csv"))

## Build table ------------------------------------------------------------

model_table <- bind_cols(as.data.frame(summary(full_model)$coefficients),
                         drop1_full_model,
                         full_model_boot_ci,
                         model_bfs) %>%
  select(effect, Estimate, SE = `Std. Error`, LowerCI = X2.5., UpperCI = X97.5., Chi2 = LRT, df = npar, p = `Pr(Chi)`, bf) %>%
  mutate(across(.cols = c(Chi2, p, bf), ~ round(.x, 2))) %>% 
  mutate(across(Chi2:bf, ~replace_na(.x, ""))) %>% 
  remove_rownames()

# ## Check assumptions -------------------------------------------------------
# 
# # Plot visualizations of model checks
# check_model(full_model)
# 
# # check for collinearity
# xx <- lm(response ~ condition + condition_order + z.trial_num + sex + z.age,
#   data = model_data
# )
# vif(xx)
# # Collinearity was no issue (maximum variance inflation factor: 1.02)


## # Model stability
## # One subject at a time excluded to assess the impact of outliers.
##
## m.stab.b <- glmm.model.stab(model.res = full_model, use = c("subject_ID"))
## m.stab.b$detailed$warnings
## xx <- as.data.frame(round(m.stab.b$summary[, -1], 3))
##
## png("graphs/full_model_stability_plot.png")
## m.stab.plot(round(m.stab.b$summary[, -1], 3))
## dev.off()
# m.stab.plot(round(m.stab.b$summary[, -1], 3))
# The model appeared to be stable with respect to the fixed effects (see full_model_stability_plot).


## # Calculate CIs for plot
## # Model with all predictor variables centered except for condition:
## full_model.CI.con <- glmer(response ~ condition + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID),
##   family = binomial, data = model_data, control = contr
## )
##
## pred.con.ci=boot.glmm.pred(model.res=full_model.CI.con,  level=0.95, use="condition", n.cores="all-1", para=T)
##
## pred.con.ci$ci.predicted
##
## write.csv(pred.con.ci$ci.predicted, file = "data/full_model_predicted_ci_for_conditionMD1_Vienna.csv")


# Breed analysis ----------------------------------------------------------

breed_data <- prereg_data %>%
  filter(condition == "ost" | condition == "non") %>%
  filter(breed!="MIX", breed!="Samoyed", breed!="White Swiss Shepherd Dog", breed!="Siberian_Husky")%>% #removing breeds with N<8
  mutate(
    z.age = as.numeric(scale(age, scale = T, center = T)),
    sex.c = as.numeric(scale(as.numeric(as.factor(sex)), scale = F, center = T)),
    condition.c = as.numeric(scale(as.numeric(as.factor(condition)), scale = F, center = T)),
    condition_order.c = as.numeric(scale(as.numeric(as.factor(condition_order)), scale = F, center = T)),
    z.trial_num = as.numeric(scale(trial_num, scale = T, center = T)),
    z.training_experience = as.numeric(scale(dog_training_experience, scale = T, center = T))
  )

contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))

breed_model <- glmer(response ~ condition + condition_order + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num |subject_ID) + (condition.c + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience|breed),
             family = binomial, data = breed_data, control = contr
)

drop1_breed_model<-drop1(breed_model, test = "Chisq", control = contr) 
round(drop1_breed_model, 3)
round(summary(breed_model)$coefficients, 3)

summary(breed_model)$varcor


## Bootstrap 95% CIs ------------------------------------------------------
# perform bootstraps for all predictors
# breed_model_boot_ci <- boot.glmm.pred(model.res = breed_model, excl.warnings = T, nboots = 1000, para = T, level = 0.95, n.cores="all-1")
# write_csv(breed_model_boot_ci$ci.estimates, here("R/breed_model_boot_ci.csv"))
breed_model_boot_ci <- read_csv(here("R/breed_model_boot_ci.csv"))

round(boot.res_breed_model$ci.estimates, 3)


## Bayes factors ----------------------------------------------------------

ncores = parallel::detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

bayes_breed_full <- brm(response ~ condition + condition_order + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num |subject_ID) + (condition.c + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience|breed), family = bernoulli, data = breed_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

bayes_breed_nocondition <- brm(response ~ condition_order + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num |subject_ID) + (condition.c + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience|breed), family = bernoulli, data = breed_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(condition_breed_bf <- bayes_factor(bayes_breed_full, bayes_breed_nocondition, repetitions = 10, silent = TRUE))

bayes_breed_noconditionorder <- brm(response ~ condition + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num |subject_ID) + (condition.c + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience|breed), family = bernoulli, data = breed_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(order_breed_bf <- bayes_factor(bayes_breed_full, bayes_breed_noconditionorder, repetitions = 10, silent = TRUE))

bayes_breed_notrialnum <- brm(response ~ condition + condition_order + sex + z.age + z.training_experience + (condition.c + z.trial_num |subject_ID) + (condition.c + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience|breed), family = bernoulli, data = breed_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(trialnum_breed_bf <- bayes_factor(bayes_breed_full, bayes_breed_notrialnum, repetitions = 10, silent = TRUE))

bayes_breed_nosex <- brm(response ~ condition + condition_order + z.trial_num + z.age + z.training_experience + (condition.c + z.trial_num |subject_ID) + (condition.c + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience|breed), family = bernoulli, data = breed_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(sex_breed_bf <- bayes_factor(bayes_breed_full, bayes_breed_nosex, repetitions = 10, silent = TRUE))

bayes_breed_noage <- brm(response ~ condition + condition_order + z.trial_num + sex + z.training_experience + (condition.c + z.trial_num |subject_ID) + (condition.c + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience|breed), family = bernoulli, data = breed_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(age_breed_bf <- bayes_factor(bayes_breed_full, bayes_breed_noage, repetitions = 10, silent = TRUE))

bayes_breed_notraining <- brm(response ~ condition + condition_order + z.trial_num + sex + z.age + (condition.c + z.trial_num |subject_ID) + (condition.c + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience|breed), family = bernoulli, data = breed_data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(training_breed_bf <- bayes_factor(bayes_breed_full, bayes_breed_notraining, repetitions = 10, silent = TRUE))

model_breed_bfs <- tibble(effect = c("(Intercept)", "Condition", "Order of condition", "Trial number", "Sex", "Age", "C-BARQ trainability score"), bf = c(NA, condition_breed_bf$bf_median_based, order_breed_bf$bf_median_based, trialnum_breed_bf$bf_median_based, sex_breed_bf$bf_median_based, age_breed_bf$bf_median_based, training_breed_bf$bf_median_based))
model_breed_bfs
write_csv(model_breed_bfs, here("R/model_breed_bfs.csv"))
model_breed_bfs <- read_csv(here("R/model_breed_bfs.csv"))

## Build table ------------------------------------------------------------

breed_table <- bind_cols(as.data.frame(summary(breed_model)$coefficients),
                         drop1_breed_model,
                         breed_model_boot_ci,
                         model_breed_bfs) %>%
  select(effect, Estimate, SE = `Std. Error`, LowerCI = X2.5., UpperCI = X97.5., Chi2 = LRT, df = npar, p = `Pr(Chi)`, bf) %>%
  mutate(across(.cols = c(Chi2, p, bf), ~ round(.x, 2))) %>% 
  mutate(across(Chi2:bf, ~replace_na(.x, ""))) %>% 
  remove_rownames()




# pred.names = tibble(effect=c("(Intercept)", "Condition", "Order of condition", "Trial number", "Sex", "Age", "C-BARQ trainability score"))
# 
# model_table_breed_model <- bind_cols(as.data.frame(summary(breed_model)$coefficients),
#                                      drop1_breed_model,
#                                      boot.res_breed_model$ci.estimates,
#                                      pred.names) %>% #,model_bfs
#   select(effect, Estimate, SE = `Std. Error`, LowerCI = X2.5., UpperCI = X97.5., Chi2 = LRT, df = npar, p = `Pr(Chi)`) %>%#, bf
#   mutate(across(.cols = c(p), ~ round(.x, 3))) %>% #bf
#   mutate(across(.cols = c(Estimate:Chi2), ~ round(.x, 2))) %>% 
#   mutate(across(Chi2:p, ~replace_na(.x, ""))) %>% #:bf
#   remove_rownames()
# 
# model_table_breed_model.re <-as.data.frame(VarCorr(breed_model))%>%
#   filter(is.na(var2))%>%
#   select("Random.intercept" = grp, "Random.slope" = var1, "SD" = sdcor)%>%
#   mutate(across(.cols = c(SD), ~ round(.x, 2)))%>%
#   mutate(Random.intercept=fct_recode(Random.intercept,  "Subject"="subject_ID", "Breed"="breed"), Random.slope=fct_recode(Random.slope,  "Condition"="condition.c", "Trial number"="z.trial_num",  "Order of conditions"="condition_order.c",  "Sex"="sex.c",  "Age"="z.age",  "C-BARQ trainability score"="z.training_experience"))
# # check for colliniarity
# xx3 <- lm(response ~ condition + condition_order + z.trial_num + sex + z.age,
#           data = breed_data
# )
# vif(xx3)
# 
# ### Model stability
# m.stab.b_breed_model <- glmm.model.stab(model.res = breed_model, use = c("subject_ID", "breed"))
# m.stab.b_breed_model$detailed$warnings
# xx2 <- as.data.frame(round(m.stab.b$summary[, -1], 3))
# # png("graphs/breed_model_stability_plot.png")
# # m.stab.plot(round(m.stab.b_breed_model$summary[, -1], 3))
# # dev.off()
# m.stab.plot(round(m.stab.b_breed_model$summary[, -1], 3))
