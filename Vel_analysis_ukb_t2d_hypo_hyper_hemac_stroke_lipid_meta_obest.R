
library(tidyverse)
library(skimr)

load("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/Luca/returned_datasets/ukb_first_events.rdata")

ukb_vel_ICD = ukb_first_events %>% 
  select( c("identifier", "cal_ps_d_079", "cal_ps_d_084", "cal_ps_d_085", "cal_ps_d_086",
            "cal_ps_d_087", "cal_ps_d_088", "cal_ps_d_089"))
# Set the name to something more interpretable
# Swap vector
switch_name = c( 
  "identifier" = "eid",
  "cal_ps_d_079"= "VEL_T2D",
  "cal_ps_d_084"= "VEL_haemochromatosis",
  "cal_ps_d_085"= "VEL_hypothrodism",
  "cal_ps_d_086"= "VEL_hyperthyroidism",
  "cal_ps_d_087"= "VEL_lipid_metabolism", 
  "cal_ps_d_088"= "VEL_obesity", 
  "cal_ps_d_089"= "VEL_ischemic_heart_disease"
)
# Change names
for (repl in 1:length(switch_name)) {
  # position of the match
  position = grep(pattern = names(switch_name)[repl], 
                  x = colnames(ukb_vel_ICD), )
  # change
  valUK = gsub(pattern = names(switch_name)[repl], 
               replacement = switch_name[repl], 
               x = colnames(ukb_vel_ICD), )
  # replace colnames vector
  colnames(ukb_vel_ICD)[position] = valUK[position]
}
# explore the data
skim(ukb_vel_ICD)

# Add basic data from UKB
# original UKB df shared with Mattia
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

# DF that has:
# Unrelated european
# WT
# Het
# Hom
df_final_un_eu = df_final[df_final$eid %in% eur_unre_500k$ID ,]
# count the number of VEL cases:
# df_final_un_eu = df_final[df_final$eid %in% eur_unre_500k$ID ,] # there are 90
# df_final_un_eu = df_final[df_final$eid %in% eur_unre_500k$ID ,] # there are 11849

# DF that has:
# Unrelated european
# WT
# Hom
df_final_un_eu_only_hom = df_final_un_eu[which(df_final_un_eu$genotype == '0' | df_final_un_eu$genotype == '2'),]

# Select the basic infromation to transfer to the new df
df_additional_info = df_final_un_eu %>% 
  select( "eid", "sex", "bmi", "ages", "smallbin", "genotype")

# Create the df for stat analysis
df_analysis = merge(df_additional_info, ukb_vel_ICD, by = "eid")

#convert all data columns to 0/1
value_cols = c("VEL_T2D", "VEL_haemochromatosis", "VEL_hypothrodism", "VEL_hyperthyroidism", 
               "VEL_lipid_metabolism", "VEL_obesity", "VEL_ischemic_heart_disease")
for (vect in value_cols) {
  df_analysis[,vect] = ifelse(is.na(df_analysis[,vect]), 0, 1)
}

#Convert male and female to 0/1
df_analysis$sex = ifelse( df_analysis$sex  == "Female", 0, 1) 
# Female == 0 
# Male == 1

# df that has to be used in the statistical analysis below
df_to_use = df_analysis
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
      tmp_df = df_to_use[which(df_to_use$sex == gen),]
      if (gen == 0){
        for(i in c(7:length(colnames(tmp_df)))){
          normalise_tmp = tmp_df[,i] # no normalisation
          # normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
          model1 = glm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages ) 
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
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
        for(i in c(7:length(colnames(tmp_df)))){
          normalise_tmp = tmp_df[,i] # no normalisation
          # normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
          model1 = glm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages )
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
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
    for(i in c(7:length(colnames(df_to_use)))){
      normalise_tmp = df_to_use[,i] # no normalisation
      # normalise_tmp = ( df_to_use[,i] - mean(as.vector(df_to_use[,i]), na.rm = TRUE) ) / sd(df_to_use[,i], na.rm = TRUE) 
      model1 = glm(normalise_tmp ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex ) 
      param = parameters::model_parameters(model1)
      for (num in 1:length(unique(df_to_use$genotype))) {
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
      tmp_df = df_to_use[which(df_to_use$sex == gen),]
      if (gen == 0){
        for(i in c(7:length(colnames(tmp_df)))){
          normalise_tmp = tmp_df[,i] # no normalisation
          # normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
          model1 = glm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi )
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
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
        for(i in c(7:length(colnames(tmp_df)))){
          normalise_tmp = tmp_df[,i] # no normalisation
          # normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
          model1 = glm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi )
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
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
    for(i in c(7:length(colnames(df_to_use)))){
      model1 = glm(df_to_use[,i] ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex + df_to_use$bmi ) 
      param = parameters::model_parameters(model1)
      for (num in 1:length(unique(df_to_use$genotype))) {
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



height_df = data.table::fread('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/ukbcvo/height_otUKBcnv.tab')

height_df = height_df %>% 
  select(c("f.eid", "f.50.0.0", "f.20015.0.0"))

colnames(height_df) = c("eid", "standing_height", "sitting height")

df2 = merge(
  df_to_use,
  height_df,
  by="eid"
)

df_to_use = df2

df_to_use$genotype = relevel(df_to_use$genotype, ref = "0")

for (round_bmi in 1:2) {
  if (round_bmi == 1 ) {
    for (gen in 0:1) {
      tmp_df = df_to_use[which(df_to_use$sex == gen),]
      if (gen == 0){
        for(i in c(14:length(colnames(tmp_df)))){
          normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
          model1 = lm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages ) 
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
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
        for(i in c(14:length(colnames(tmp_df)))){
          normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
          model1 = lm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages )
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
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
    for(i in c(14:length(colnames(df_to_use)))){
      normalise_tmp = ( df_to_use[,i] - mean(as.vector(df_to_use[,i]), na.rm = TRUE) ) / sd(df_to_use[,i], na.rm = TRUE) 
      model1 = lm(normalise_tmp ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex ) 
      param = parameters::model_parameters(model1)
      for (num in 1:length(unique(df_to_use$genotype))) {
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
      tmp_df = df_to_use[which(df_to_use$sex == gen),]
      if (gen == 0){
        for(i in c(14:length(colnames(tmp_df)))){
          normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
          model1 = lm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi )
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
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
        for(i in c(14:length(colnames(tmp_df)))){
          normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE) 
          model1 = lm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi )
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
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
    for(i in c(14:length(colnames(df_to_use)))){
      model1 = lm(df_to_use[,i] ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex + df_to_use$bmi ) 
      param = parameters::model_parameters(model1)
      for (num in 1:length(unique(df_to_use$genotype))) {
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
  filter(cohort != "(Intercept)") 

write.table(stat_ukb, file = '~/Desktop/VEL/result_ukb_vel/UKB_VEL_analysis_ICD_T2D_haemochromatosis_hypothrodism_hyperthyroidism_lipid_met_obesity_ischemic.csv', row.names = F, append = F, quote = F, sep = ",")


ggplot(df_to_use,
       aes(x=genotype, y=statin)) +
  geom_col() + 
  facet_wrap('sex_f31_0_0') +
  ggtitle(label = feat)
ggsave(filename = paste0("boxplot_", feat, ".svg" ),
       path = 'Desktop/VEL/result_ukb_vel/' ,
       device = 'svg', plot = p2) 

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

raincloud_theme <- theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.position = "right",
  plot.title = element_text(lineheight = .8, face = "bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
  axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"))

for (feata in colnames(df_to_use)[14:length(colnames(df_to_use))]) {
  if (is.numeric(df_to_use[,feata])){
    png(filename = paste0("Desktop/VEL/result_ukb_vel/rainplot_", 'feata', "_20210514.png" ), 
        bg = "transparent",
        width = 1024, 
        height = 768 )
    p1 = ggplot(data = df_to_use, aes( x = genotype, y = df_to_use[,feata], fill = genotype)) +
      ggtitle(label = feata) +
      geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, bw = 1) +
      geom_point(aes(color = genotype), position = position_jitter(width = .15), size = .3, alpha = 0.5) +
      geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
      #expand_limits(x = 3.25) +
      guides(fill = FALSE) +
      guides(color = FALSE) +
      scale_color_manual(values = ggsci::pal_npg("nrc")(3) ) +
      scale_fill_manual(values = ggsci::pal_npg("nrc")(3) ) +
      theme_minimal() +
      facet_wrap('sex') 
    ggsave(filename = paste0("rainplot_", feata, "_20210514.png.png" ),
           path = 'Desktop/VEL/result_ukb_vel/' ,
           device = 'png', plot = p1)
    dev.off()
  } else {
    p2 = ggplot(as.data.frame(table(df_to_use[,feata])),
                aes(x=Var1, y=Freq)) +
      geom_col()  + 
      ggtitle(label = feata)
    ggsave(filename = paste0("boxplot_", feata, ".svg" ),
           path = 'Desktop/VEL/result_ukb_vel/' ,
           device = 'svg', plot = p2)
    dev.off()
  }
}
