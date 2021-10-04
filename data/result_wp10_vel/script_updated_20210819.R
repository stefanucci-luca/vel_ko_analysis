library('lm.beta')
library("questionr")
library('esc')
library('effectsize')
library('tidyverse')
library('see')
library('GGally')
library(RColorBrewer)
library(RNOmni)

dt_WP10 = xlsx::read.xlsx("/Volumes/GoogleDrive-105684671225474146017/.shortcut-targets-by-id/1zpTtL0t1Um6lyloW4ZxkWYpH4dMkbiBD/Manuscript/Tables/Table S4.xlsx", 1)

str(dt_WP10)
colnames(dt_WP10)

tt = c("LAR", "GLUC", "CHOL", "TG", "HDL", "ALT", "LEPT", "ADPN", "AST" ,"Inorganic_Phosphate",
       "Ferritin", "hsCRP", "FFA", "Total_T3", "Total_T4", "TSH", "LDL")
# Not gender stratified and BMI as covariate
feature = c()
effect = c()
cohort = c()
CI_low = c()
CI_high = c()
p.val = c()

for (variable in tt) {
  df_tmp = dt_WP10[which(!is.na(dt_WP10[,variable])),]
  normal_rINV = RankNorm(as.vector(df_tmp[,variable]))
  model = lm(data = df_tmp, 
              formula = normal_rINV ~ COHORT + GENDER + BMI + AGE)
  param = parameters::parameters(model)
  feature = c(feature, variable)
  effect = c(effect, param$Coefficient[2])
  cohort = c(cohort, param$Parameter[2])
  CI_low = c(CI_low, param$CI_low[2])
  CI_high = c(CI_high, param$CI_high[2])
  p.val = c(p.val, param$p[2])
}

stat_df = data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  'CI_low' = CI_low, 
  'CI_high' = CI_high,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr")
) %>%  arrange(p.val)


# female and BMI as covariate
feature = c()
effect = c()
cohort = c()
CI_low = c()
CI_high = c()
p.val = c()

for (variable in tt) {
  df_tmp = dt_WP10[which(!is.na(dt_WP10[,variable]) & dt_WP10[,"GENDER"] == "Female"),]
  normal_rINV = RankNorm(as.vector(df_tmp[,variable]))
  model = lm(data = df_tmp, 
             formula = normal_rINV ~ COHORT + BMI + AGE)
  param = parameters::parameters(model)
  feature = c(feature, variable)
  effect = c(effect, param$Coefficient[2])
  cohort = c(cohort, param$Parameter[2])
  CI_low = c(CI_low, param$CI_low[2])
  CI_high = c(CI_high, param$CI_high[2])
  p.val = c(p.val, param$p[2])
}

stat_df = data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  'CI_low' = CI_low, 
  'CI_high' = CI_high,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr")
) %>%  arrange(p.val)

# male and BMI as covariate
feature = c()
effect = c()
cohort = c()
CI_low = c()
CI_high = c()
p.val = c()

for (variable in tt) {
  df_tmp = dt_WP10[which(!is.na(dt_WP10[,variable]) & dt_WP10[,"GENDER"] == "Male"),]
  normal_rINV = RankNorm(as.vector(df_tmp[,variable]))
  model = lm(data = df_tmp, 
             formula = normal_rINV ~ COHORT + BMI + AGE)
  param = parameters::parameters(model)
  feature = c(feature, variable)
  effect = c(effect, param$Coefficient[2])
  cohort = c(cohort, param$Parameter[2])
  CI_low = c(CI_low, param$CI_low[2])
  CI_high = c(CI_high, param$CI_high[2])
  p.val = c(p.val, param$p[2])
}

stat_df = data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  'CI_low' = CI_low, 
  'CI_high' = CI_high,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr")
) %>%  arrange(p.val)








# Not gender stratified and NO BMI as covariate
feature = c()
effect = c()
cohort = c()
CI_low = c()
CI_high = c()
p.val = c()

for (variable in tt) {
  df_tmp = dt_WP10[which(!is.na(dt_WP10[,variable])),]
  normal_rINV = RankNorm(as.vector(df_tmp[,variable]))
  model = lm(data = df_tmp, 
             formula = normal_rINV ~ COHORT + GENDER + AGE)
  param = parameters::parameters(model)
  feature = c(feature, variable)
  effect = c(effect, param$Coefficient[2])
  cohort = c(cohort, param$Parameter[2])
  CI_low = c(CI_low, param$CI_low[2])
  CI_high = c(CI_high, param$CI_high[2])
  p.val = c(p.val, param$p[2])
}

stat_df = data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  'CI_low' = CI_low, 
  'CI_high' = CI_high,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr")
) %>%  arrange(p.val)


# female and NO BMI as covariate
feature = c()
effect = c()
cohort = c()
CI_low = c()
CI_high = c()
p.val = c()

for (variable in tt) {
  df_tmp = dt_WP10[which(!is.na(dt_WP10[,variable]) & dt_WP10[,"GENDER"] == "Female"),]
  normal_rINV = RankNorm(as.vector(df_tmp[,variable]))
  model = lm(data = df_tmp, 
             formula = normal_rINV ~ COHORT + AGE)
  param = parameters::parameters(model)
  feature = c(feature, variable)
  effect = c(effect, param$Coefficient[2])
  cohort = c(cohort, param$Parameter[2])
  CI_low = c(CI_low, param$CI_low[2])
  CI_high = c(CI_high, param$CI_high[2])
  p.val = c(p.val, param$p[2])
}

stat_df = data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  'CI_low' = CI_low, 
  'CI_high' = CI_high,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr")
) %>%  arrange(p.val)

# male and NO BMI as covariate
feature = c()
effect = c()
cohort = c()
CI_low = c()
CI_high = c()
p.val = c()

for (variable in tt) {
  df_tmp = dt_WP10[which(!is.na(dt_WP10[,variable]) & dt_WP10[,"GENDER"] == "Male"),]
  normal_rINV = RankNorm(as.vector(df_tmp[,variable]))
  model = lm(data = df_tmp, 
             formula = normal_rINV ~ COHORT + AGE)
  param = parameters::parameters(model)
  feature = c(feature, variable)
  effect = c(effect, param$Coefficient[2])
  cohort = c(cohort, param$Parameter[2])
  CI_low = c(CI_low, param$CI_low[2])
  CI_high = c(CI_high, param$CI_high[2])
  p.val = c(p.val, param$p[2])
}

stat_df = data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  'CI_low' = CI_low, 
  'CI_high' = CI_high,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr")
) %>%  arrange(p.val)















