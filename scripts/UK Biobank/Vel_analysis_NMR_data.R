
# Read in the NMR data
nmr_df = readRDS("NMR_data.RData")
# relabel the genotype column
colnames(nmr_df)[25] = "genotype"
headers = unique(stringr::str_split(colnames(nmr_df), "_f", simplify=T, n=2)[,1])[26:length(unique(stringr::str_split(colnames(nmr_df), "_f", simplify=T, n=2)[,1]))]



for (feat in headers){
  print(feat)
  cols = grepl(paste0("^", feat), colnames(nmr_df))
  nmr_df[,feat] = apply(nmr_df[,cols], 1, function(x) mean(x, na.rm = T))
}

saveRDS(object = nmr_df, file = "/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/VEL/NMR_data.RData" )








df_to_use_pheno = nmr_df

gender = c()
feature = c()
effect = c()
cohort = c()
sd = c()
CI_low = c()
CI_high = c()
p.val = c()
bmi = c()

for(i in 362:dim(df_to_use_pheno)[2]){
  print(i)
  normalise_tmp = df_to_use_pheno[,i] # no normalisation
  model1 = lm(normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages + df_to_use_pheno$sex + df_to_use_pheno$bmi) 
  param = parameters::model_parameters(model1, df_method="wald")
  for (num in 1:length(unique(df_to_use_pheno$genotype))) {
    feature = c(feature, colnames(df_to_use_pheno)[i])
    gender = c(gender, 'not_gender_stratified but as covariates')
    cohort = c(cohort,param$Parameter[num])
    effect = c(effect, param$Coefficient[num])
    sd = c(sd, param$SE[num])
    CI_low = c(CI_low, param$CI_low[num] )
    CI_high = c(CI_high, param$CI_high[num] )
    p.val = c(p.val, param$p[num])
    bmi = c(bmi, 'BMI_as_covariates')
  }
}

for(i in 362:dim(df_to_use_pheno)[2]){
  print(i)
  normalise_tmp = df_to_use_pheno[,i] # no normalisation
  model1 = lm(normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages + df_to_use_pheno$sex, family = "binomial") 
  param = parameters::model_parameters(model1, df_method="wald")
  for (num in 1:length(unique(df_to_use_pheno$genotype))) {
    feature = c(feature, colnames(df_to_use_pheno)[i])
    gender = c(gender, 'not_gender_stratified but as covariates')
    cohort = c(cohort,param$Parameter[num])
    effect = c(effect, param$Coefficient[num])
    sd = c(sd, param$SE[num])
    CI_low = c(CI_low, param$CI_low[num] )
    CI_high = c(CI_high, param$CI_high[num] )
    p.val = c(p.val, param$p[num])
    bmi = c(bmi, 'NOT_BMI_as_covariates')
  }
}


#Male

df_to_use_pheno = nmr_df[which(nmr_df$sex == "Male"),]

for(i in 362:dim(df_to_use_pheno)[2]){
  print(i)
  normalise_tmp = df_to_use_pheno[,i] # no normalisation
  model1 = lm(normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages + df_to_use_pheno$bmi) 
  param = parameters::model_parameters(model1, df_method="wald")
  for (num in 1:length(unique(df_to_use_pheno$genotype))) {
    feature = c(feature, colnames(df_to_use_pheno)[i])
    gender = c(gender, 'Male')
    cohort = c(cohort,param$Parameter[num])
    effect = c(effect, param$Coefficient[num])
    sd = c(sd, param$SE[num])
    CI_low = c(CI_low, param$CI_low[num] )
    CI_high = c(CI_high, param$CI_high[num] )
    p.val = c(p.val, param$p[num])
    bmi = c(bmi, 'BMI_as_covariates')
  }
}

stat_ukb_nmr =  data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  'sd' = sd,
  'CI_low' = CI_low,
  'CI_high' = CI_high ,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr"),
  'gender' = gender,
  'bmi' = bmi
) %>%  
  arrange(p.val) %>%
  filter(cohort != "(Intercept)") 

for(i in 362:dim(df_to_use_pheno)[2]){
  print(i)
  normalise_tmp = df_to_use_pheno[,i] # no normalisation
  model1 = lm(normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages) 
  param = parameters::model_parameters(model1, df_method="wald")
  for (num in 1:length(unique(df_to_use_pheno$genotype))) {
    feature = c(feature, colnames(df_to_use_pheno)[i])
    gender = c(gender, 'Male')
    cohort = c(cohort,param$Parameter[num])
    effect = c(effect, param$Coefficient[num])
    sd = c(sd, param$SE[num])
    CI_low = c(CI_low, param$CI_low[num] )
    CI_high = c(CI_high, param$CI_high[num] )
    p.val = c(p.val, param$p[num])
    bmi = c(bmi, 'NOT_BMI_as_covariates')
  }
}

stat_ukb_nmr =  data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  'sd' = sd,
  'CI_low' = CI_low,
  'CI_high' = CI_high ,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr"),
  'gender' = gender,
  'bmi' = bmi
) %>%  
  arrange(p.val) %>%
  filter(cohort != "(Intercept)") 

# Female

df_to_use_pheno = nmr_df[which(nmr_df$sex == "Female"),]

for(i in 362:dim(df_to_use_pheno)[2]){
  print(i)
  normalise_tmp = df_to_use_pheno[,i] # no normalisation
  model1 = lm(normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages + df_to_use_pheno$bmi) 
  param = parameters::model_parameters(model1, df_method="wald")
  for (num in 1:length(unique(df_to_use_pheno$genotype))) {
    feature = c(feature, colnames(df_to_use_pheno)[i])
    gender = c(gender, 'Female')
    cohort = c(cohort,param$Parameter[num])
    effect = c(effect, param$Coefficient[num])
    sd = c(sd, param$SE[num])
    CI_low = c(CI_low, param$CI_low[num] )
    CI_high = c(CI_high, param$CI_high[num] )
    p.val = c(p.val, param$p[num])
    bmi = c(bmi, 'BMI_as_covariates')
  }
}

for(i in 362:dim(df_to_use_pheno)[2]){
  print(i)
  normalise_tmp = df_to_use_pheno[,i] # no normalisation
  model1 = lm(normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages) 
  param = parameters::model_parameters(model1, df_method="wald")
  for (num in 1:length(unique(df_to_use_pheno$genotype))) {
    feature = c(feature, colnames(df_to_use_pheno)[i])
    gender = c(gender, 'Female')
    cohort = c(cohort,param$Parameter[num])
    effect = c(effect, param$Coefficient[num])
    sd = c(sd, param$SE[num])
    CI_low = c(CI_low, param$CI_low[num] )
    CI_high = c(CI_high, param$CI_high[num] )
    p.val = c(p.val, param$p[num])
    bmi = c(bmi, 'NOT_BMI_as_covariates')
  }
}

library(tidyverse)
# Build the df for the statistical analysis
stat_ukb_nmr =  data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  'sd' = sd,
  'CI_low' = CI_low,
  'CI_high' = CI_high ,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr"),
  'gender' = gender,
  'bmi' = bmi
) %>%  
  arrange(p.val) %>%
  filter(cohort != "(Intercept)") 

