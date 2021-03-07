library(tidyverse)
library(psych)
library(summarytools)
library(lme4)
library(car)
library(BayesFactor)
library(bayestestR)
library(papaja)
library(brms)

source("./functions/glmm_stability.r")
source("./functions/boot_glmm.r")

load("md1_vienna_glmm_workspace.RData")



model.data <- prereg.data  %>%
  filter(condition == "ost" | condition == "non") %>%
  mutate(
    z.age = as.numeric(scale(age, scale = T, center = T)),
    sex.c = as.numeric(scale(as.numeric(as.factor(sex)), scale = F, center = T)),
    condition.c = as.numeric(scale(as.numeric(as.factor(condition)), scale = F, center = T)),
    condition_order.c = as.numeric(scale(as.numeric(as.factor(condition_order)), scale = F, center = T)),
    z.trial_num = as.numeric(scale(trial_num, scale = T, center = T)),
    z.training_experience = as.numeric(scale(dog_training_experience, scale = T, center = T))
  )

library(rstan)
ncores = parallel::detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

bayes_full <- brm(response ~ condition + condition_order + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model.data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

bayes_nocondition <- brm(response ~ condition_order + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model.data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(condition_bf <- bayes_factor(bayes_full, bayes_nocondition, repetitions = 10, silent = TRUE))

bayes_noconditionorder <- brm(response ~ condition + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model.data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(order_bf <- bayes_factor(bayes_full, bayes_noconditionorder, repetitions = 10, silent = TRUE))

bayes_notrialnum <- brm(response ~ condition + condition_order + sex + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model.data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(trialnum_bf <- bayes_factor(bayes_full, bayes_notrialnum, repetitions = 10, silent = TRUE))

bayes_nosex <- brm(response ~ condition + condition_order + z.trial_num + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model.data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(sex_bf <- bayes_factor(bayes_full, bayes_nosex, repetitions = 10, silent = TRUE))

bayes_noage <- brm(response ~ condition + condition_order + z.trial_num + sex + z.training_experience + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model.data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(age_bf <- bayes_factor(bayes_full, bayes_noage, repetitions = 10, silent = TRUE))

bayes_notraining <- brm(response ~ condition + condition_order + z.trial_num + sex + z.age + (condition.c + z.trial_num | subject_ID), family = bernoulli, data = model.data, save_pars = save_pars(all = TRUE), iter = 12000, warmup = 2000)

(training_bf <- bayes_factor(bayes_full, bayes_notraining, repetitions = 10, silent = TRUE))

model_bfs <- tibble(fixedeffect = c("Condition", "Order of condition", "Trial number", "Sex", "Age", "C-BARQ trainability score"), bf = c(condition_bf$bf_median_based, order_bf$bf_median_based, trialnum_bf$bf_median_based, sex_bf$bf_median_based, age_bf$bf_median_based, training_bf$bf_median_based))
model_bfs
readr::write_csv(model_bfs, "mm1_bfs.csv")