library(tidyverse)

df_large = data.table::fread("/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/UKB_VEL_relevant_fields_20210308.csv")
df_large = as.data.frame(df_large)
# remove un-assigned genotypes
df_large = filter(df_large, genotype != "./." )
# Rename genotype as:
# 0/0 = 0
# 0/1 = 1
# 1/1 = 2
df_large$genotype = str_replace_all(string = df_large$genotype, pattern = "0/0", replacement = "0")
df_large$genotype = str_replace_all(string = df_large$genotype, pattern = "0/1", replacement = "1")
df_large$genotype = str_replace_all(string = df_large$genotype, pattern = "1/1", replacement = "2")
# set WT allele (i.e. "0/0" aka "0") as the one to use for reference
df_large$genotype <- relevel(as.factor(df_large$genotype), ref = "0")

# import ancestry information 
eur_unre_500k = data.table::fread('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/british-ancestry-MSUP_unrel_europ_500K_from_will_20210317.tsv') %>% 
  select(ID)

#Covariates
load('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/Luanluan/data/ukb_dvt_pe_combined.rdata')
# times in the hospital
df_final = merge(df, df_large, by = 'eid' )

df_times_hosp = data.table::fread('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/ukbcvo/times_in_hopital_otUKBcnv.tab', 
                                  col.names = c('eid','times_in_hospital'))

df_final = merge(df_final, df_times_hosp, by = 'eid')
# n_hom = dim(df_final[df_final$genotype == 2,])[1]
# n_het = dim(df_final[df_final$genotype == 2,])[1]

# DF that has:
# Unrelated european
# WT
# Het
# Hom
df_final_un_eu = df_final[df_final$eid %in% eur_unre_500k$ID ,]
# count the number of VEL cases:
# df_final_un_eu = df_final[df_final$eid %in% eur_unre_500k$ID ,] # there are 90
# df_final_un_eu = df_final[df_final$eid %in% eur_unre_500k$ID ,] # there are 11849

eqtl_variant_genotype = read.delim('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/VEL/eQTL_snps_in_UKB.raw', stringsAsFactors = F)
eqtl_variant_genotype=as.data.frame(eqtl_variant_genotype)

notLD_var = as.data.frame(data.table::fread('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/VEL/SNP_not_in_LD.tsv', header = F))

df_final = merge(df_final_un_eu, eqtl_variant_genotype, by.x  = 'eid', by.y = "IID")
df_to_use = df_final

for (snp in notLD_var$V1) {

  # vector for stats
  gender = c()
  feature = c()
  effect = c()
  cohort = c()
  sd = c()
  CI_low = c()
  CI_high = c()
  p.val = c()
  bmi = c()
  
  for (round_bmi in 1:2) {
    if (round_bmi == 1 ) {
      for (gen in 0:1) {
        tmp_df = df_to_use[which(df_to_use$sex_f31_0_0 == gen),]
        if (gen == 0){
          for(i in c(23,27:77)){
            # normalise_tmp = tmp_df[,i] # no normalisation
            normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
            model1 = lm(normalise_tmp ~ factor(tmp_df[,pmatch(snp, colnames(tmp_df))]) + tmp_df$ages ) 
            param = parameters::model_parameters(model1)
            for (num in 1:2) {
              feature = c(feature, colnames(tmp_df)[i])
              gender = c(gender, gen)
              cohort = c(cohort,param$Parameter[num])
              effect = c(effect, param$Coefficient[num])
              sd = c(sd, param$SE[num])
              CI_low = c(CI_low, param$CI_low[num] )
              CI_high = c(CI_high, param$CI_high[num] )
              p.val = c(p.val, param$p[num])
              bmi = c(bmi, 'NOT_BMI_as_covariates')
            }
          }
        } else {
          for(i in c(23,27,30:77)){
            # normalise_tmp = tmp_df[,i] # no normalisation
            normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
            model1 = lm(normalise_tmp ~ factor(tmp_df[,pmatch(snp, colnames(tmp_df))]) + tmp_df$ages )
            param = parameters::model_parameters(model1)
            for (num in 1:2) {
              feature = c(feature, colnames(tmp_df)[i])
              gender = c(gender, gen)
              cohort = c(cohort,param$Parameter[num])
              effect = c(effect, param$Coefficient[num])
              sd = c(sd, param$SE[num])
              CI_low = c(CI_low, param$CI_low[num] )
              CI_high = c(CI_high, param$CI_high[num] )
              p.val = c(p.val, param$p[num])
              bmi = c(bmi, 'NOT_BMI_as_covariates')
            }
          }
        }
      }
      for(i in c(23,27:77)){
        normalise_tmp = ( df_to_use[,i] - mean(as.vector(df_to_use[,i]), na.rm = TRUE) ) / sd(df_to_use[,i], na.rm = TRUE) 
        model1 = lm(normalise_tmp ~ factor(df_to_use[,pmatch(snp, colnames(tmp_df))]) + df_to_use$ages + df_to_use$sex_f31_0_0 ) 
        param = parameters::model_parameters(model1)
        for (num in 1:2) {
          feature = c(feature, colnames(df_to_use)[i])
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
    } else {
      for (gen in 0:1) {
        tmp_df = df_to_use[which(df_to_use$sex_f31_0_0 == gen),]
        if (gen == 0){
          for(i in c(23,27:77)){
            # normalise_tmp = tmp_df[,i] # no normalisation
            normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
            model1 = lm(normalise_tmp ~ factor(tmp_df[,pmatch(snp, colnames(tmp_df))]) + tmp_df$ages + tmp_df$bmi )
            param = parameters::model_parameters(model1)
            for (num in 1:2) {
              feature = c(feature, colnames(tmp_df)[i])
              gender = c(gender, gen)
              cohort = c(cohort,param$Parameter[num])
              effect = c(effect, param$Coefficient[num])
              sd = c(sd, param$SE[num])
              CI_low = c(CI_low, param$CI_low[num] )
              CI_high = c(CI_high, param$CI_high[num] )
              p.val = c(p.val, param$p[num])
              bmi = c(bmi, 'BMI_as_covariates')
            }
          }
        } else {
          for(i in c(23,27,30:77)){
            # normalise_tmp = tmp_df[,i] # no normalisation
            normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
            model1 = lm(normalise_tmp ~ factor(tmp_df[,pmatch(snp, colnames(tmp_df))]) + tmp_df$ages + tmp_df$bmi )
            param = parameters::model_parameters(model1)
            for (num in 1:2) {
              feature = c(feature, colnames(tmp_df)[i])
              gender = c(gender, gen)
              cohort = c(cohort,param$Parameter[num])
              effect = c(effect, param$Coefficient[num])
              sd = c(sd, param$SE[num])
              CI_low = c(CI_low, param$CI_low[num] )
              CI_high = c(CI_high, param$CI_high[num] )
              p.val = c(p.val, param$p[num])
              bmi = c(bmi, 'BMI_as_covariates')
            }
          }
        }
      }
      for(i in c(23,27:77)){
        normalise_tmp = ( df_to_use[,i] - mean(as.vector(df_to_use[,i]), na.rm = TRUE) ) / sd(df_to_use[,i], na.rm = TRUE) 
        model1 = lm(normalise_tmp ~ factor(df_to_use[,pmatch(snp, colnames(tmp_df))]) + df_to_use$ages + df_to_use$sex_f31_0_0 + df_to_use$bmi ) 
        param = parameters::model_parameters(model1)
        for (num in 1:2) {
          feature = c(feature, colnames(df_to_use)[i])
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
    }
  }
  
  # Build the df for the statistical analysis
  stat_ukb = data.frame(
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
    filter(cohort != "(Intercept)") %>% 
    filter(feature != "statin" ) %>% filter( feature != "vte"  )
  
  filename=paste0(snp,'stats_table_20210325.csv')
  
  write.table(stat_ukb, file = paste0('Desktop/VEL/result_ukb_vel/', filename), row.names = F, append = F, quote = TRUE, sep = ",")
  
  
}

