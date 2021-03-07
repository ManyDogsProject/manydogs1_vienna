## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------

knitr::opts_chunk$set(echo = TRUE)

rm(list = ls())

library(tidyverse)
library(psych)
library(summarytools)
library(lme4)
library(car)
library(BayesFactor)
library(bayestestR)

source("./functions/glmm_stability.r")
source("./functions/boot_glmm.r")

load("md1_vienna_glmm_workspace.RData")


## ----include = FALSE---------------------------------------------------------------------------------------------------------------------
# Read data
prereg.data <- read.csv("data/md1_data_vienna_long.csv") %>%
  filter(preregistered_data=="after_prereg")%>% #only data that were collected after preregistration
  na.exclude() %>%
  droplevels()

summary(prereg.data)
view(dfSummary(prereg.data))


## ----include = FALSE---------------------------------------------------------------------------------------------------------------------
# Aggregating data
agg.data <- prereg.data  %>%
  group_by(subject_ID, sex, breed, age, condition) %>%
  summarise(mean.resp = mean(response)) %>%
  ungroup()

#demographics:
agg.data%>%filter(condition=="ost")%>%
  summarise(mean.age = mean(age), min.age = min(age), max.age = max(age), sd.age = sd(age), females=sum(sex=="F"), males=sum(sex=="M"))

summary.data <- agg.data %>%
  group_by(condition) %>%
  summarise(mean.resp2 = mean(mean.resp), median.resp = median(mean.resp), sd = sd(mean.resp), se = sd(mean.resp) / sqrt(length(mean.resp)))


summary.data


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------

tt.ost <- t.test(agg.data$mean.resp[agg.data$condition == "ost"], mu = 0.5, alternative = "two.sided")
tt.non <- t.test(agg.data$mean.resp[agg.data$condition == "non"], mu = 0.5, alternative = "two.sided")
tt.oc <- t.test(agg.data$mean.resp[agg.data$condition == "odour"], mu = 0.5, alternative = "two.sided")

tt.ost
tt.non
tt.oc


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------

ttbf.ost <- ttestBF(agg.data$mean.resp[agg.data$condition == "ost"], mu = 0.5, alternative = "two.sided")
ttbf.non <- ttestBF(agg.data$mean.resp[agg.data$condition == "non"], mu = 0.5, alternative = "two.sided")
ttbf.oc <- ttestBF(agg.data$mean.resp[agg.data$condition == "odour"], mu = 0.5, alternative = "two.sided")

extractBF(ttbf.ost)$bf
extractBF(ttbf.non)$bf
extractBF(ttbf.oc)$bf


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------
# centering variables for modeling

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

#view(dfSummary(model.data))


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------
contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))

mm1 <- glmer(response ~ condition + condition_order + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID),
  family = binomial, data = model.data, control = contr
)


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------
#*Likelihood ratio test*
drop1.mm1<-drop1(mm1, test = "Chisq", control = contr) 
round(drop1.mm1, 3)


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------
#*Estimates of the fixed effects*
round(summary(mm1)$coefficients, 3)


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------
#*Random effects* 
summary(mm1)$varcor


## ----eval=FALSE, include=FALSE-----------------------------------------------------------------------------------------------------------
## #*Calculating confidence intervals with 1000 bootstraps*
## #perform bootstraps for all predictors
## mm1.boot.ci <- boot.glmm.pred(model.res = mm1, excl.warnings = T, nboots = 1000, para = T, n.cores="all-1", level = 0.95)


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------
round(mm1.boot.ci$ci.estimates, 3)


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------
### Check assumptions
# check for colliniarity
xx <- lm(response ~ condition + condition_order + z.trial_num + sex + z.age,
  data = model.data
)
vif(xx)
#Collinearity was no issue (maximum variance inflation factor: 1.02)


## ----eval=FALSE, include=FALSE-----------------------------------------------------------------------------------------------------------
## # Model stability
## # One subject at a time excluded to assess the impact of outliers.
## 
## m.stab.b <- glmm.model.stab(model.res = mm1, use = c("subject_ID"))
## m.stab.b$detailed$warnings
## xx <- as.data.frame(round(m.stab.b$summary[, -1], 3))
## 
## png("graphs/mm1_stability_plot.png")
## m.stab.plot(round(m.stab.b$summary[, -1], 3))
## dev.off()
## 


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------
m.stab.plot(round(m.stab.b$summary[, -1], 3))
#The model appeared to be stable with respect to the fixed effects (see mm1_stability_plot).


## ----eval=FALSE, include=FALSE-----------------------------------------------------------------------------------------------------------
## # Calculate CIs for plot
## # Model with all predictor variables centered except for condition:
## mm1.CI.con <- glmer(response ~ condition + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID),
##   family = binomial, data = model.data, control = contr
## )


## ----eval=FALSE, include=FALSE-----------------------------------------------------------------------------------------------------------
## pred.con.ci=boot.glmm.pred(model.res=mm1.CI.con,  level=0.95, use="condition", n.cores="all-1", para=T)
## 
## pred.con.ci$ci.predicted
## 
## write.csv(pred.con.ci$ci.predicted, file = "data/mm1_predicted_ci_for_conditionMD1_Vienna.csv")


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------
#save.image("md1_vienna_glmm_workspace.RData")


## ----include=FALSE-----------------------------------------------------------------------------------------------------------------------
knitr::knit_exit()


## ----include = FALSE---------------------------------------------------------------------------------------------------------------------
all.data <- read.csv("data/md1_data_vienna_long.csv") %>%
  na.exclude() %>%
  droplevels()

summary(all.data)
#view(dfSummary(all.data))


## ----include = FALSE---------------------------------------------------------------------------------------------------------------------

agg.data.all <- all.data %>%
  group_by(subject_ID, sex, breed, age, condition) %>%
  summarise(mean.resp = mean(response)) %>%
  ungroup()

#demographics:
agg.data.all%>%filter(condition=="ost")%>%
  summarise(mean.age = mean(age), min.age = min(age), max.age = max(age), sd.age = sd(age), females=sum(sex=="F"), males=sum(sex=="M"))

summary.data.all <- agg.data.all %>%
  group_by(condition) %>%
  summarise(mean.resp2 = mean(mean.resp), median.resp = median(mean.resp), sd = sd(mean.resp), se = sd(mean.resp) / sqrt(length(mean.resp)))


summary.data.all


## ----------------------------------------------------------------------------------------------------------------------------------------

tt.ost.all <- t.test(agg.data.all$mean.resp[agg.data.all$condition == "ost"], mu = 0.5, alternative = "two.sided")
tt.non.all <- t.test(agg.data.all$mean.resp[agg.data.all$condition == "non"], mu = 0.5, alternative = "two.sided")
tt.oc.all <- t.test(agg.data.all$mean.resp[agg.data.all$condition == "odour"], mu = 0.5, alternative = "two.sided")

tt.ost.all
tt.non.all
tt.oc.all


## ----------------------------------------------------------------------------------------------------------------------------------------
# centering variables for modeling

model.data.all <- all.data %>%
  filter(condition == "ost" | condition == "non") %>%
  mutate(
    z.age = as.numeric(scale(age, scale = T, center = T)),
    sex.c = as.numeric(scale(as.numeric(as.factor(sex)), scale = F, center = T)),
    condition.c = as.numeric(scale(as.numeric(as.factor(condition)), scale = F, center = T)),
    condition_order.c = as.numeric(scale(as.numeric(as.factor(condition_order)), scale = F, center = T)),
    z.trial_num = as.numeric(scale(trial_num, scale = T, center = T)),
    z.training_experience = as.numeric(scale(dog_training_experience, scale = T, center = T))
  )

#view(dfSummary(model.data.all))


## ----------------------------------------------------------------------------------------------------------------------------------------
contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))

mm2 <- glmer(response ~ condition + condition_order + z.trial_num + sex + z.age + z.training_experience + (condition.c + z.trial_num |subject_ID),
  family = binomial, data = model.data.all, control = contr
)


## ----------------------------------------------------------------------------------------------------------------------------------------
drop1.mm2<-drop1(mm2, test = "Chisq", control = contr) 
round(drop1.mm2, 3)


## ----------------------------------------------------------------------------------------------------------------------------------------
round(summary(mm2)$coefficients, 3)


## ----------------------------------------------------------------------------------------------------------------------------------------
summary(mm2)$varcor


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------
## # perform bootstraps for all predictors
## boot.res.mm2 <- boot.glmm.pred(model.res = mm2, excl.warnings = T, nboots = 1000, para = T, level = 0.95, n.cores="all-1")

## ----------------------------------------------------------------------------------------------------------------------------------------
round(boot.res.mm2$ci.estimates, 3)


## ----------------------------------------------------------------------------------------------------------------------------------------
# check for colliniarity
xx2 <- lm(response ~ condition + condition_order + z.trial_num + sex + z.age,
  data = model.data
)
vif(xx2)


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------
## m.stab.b.mm2 <- glmm.model.stab(model.res = mm1, use = c("subject_ID"))
## m.stab.b.mm2$detailed$warnings
## xx2 <- as.data.frame(round(m.stab.b$summary[, -1], 3))
## png("graphs/mm2_stability_plot.png")
## m.stab.plot(round(m.stab.b.mm2$summary[, -1], 3))
## dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------------
m.stab.plot(round(m.stab.b.mm2$summary[, -1], 3))


## ----------------------------------------------------------------------------------------------------------------------------------------
# Model with all predictor variables centered except for condition:
mm2.CI.con <- glmer(response ~ condition + condition_order.c + z.trial_num + sex.c + z.age + z.training_experience + (condition.c + z.trial_num | subject_ID),
  family = binomial, data = model.data, control = contr
)


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------
## mm2.pred.con.ci=boot.glmm.pred(model.res=mm2.CI.con,  level=0.95, use="condition", n.cores="all-1", para=T)
## 
## mm2.pred.con.ci$ci.predicted
## 
## write.csv(mm2.pred.con.ci$ci.predicted, file = "data/mm2_predicted_ci_for_condition_MD1_Vienna.csv")


## ----------------------------------------------------------------------------------------------------------------------------------------
#save.image("md1_vienna_glmm_workspace.RData")

