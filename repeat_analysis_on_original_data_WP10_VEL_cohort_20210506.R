library('lm.beta')
library("questionr")
library('esc')
library('effectsize')
library('tidyverse')
library('see')
library('GGally')
library(skimr)
library(RColorBrewer)

# Import VEL phenotype data
wp10_VEL=xlsx::read.xlsx("/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/VEL/Original files CBAL/C19-078 Manifest Assays March 2019 - results (2).xlsx", 
                sheetIndex = 1, endRow = 261)

wp10_VEL_wrong=xlsx::read.xlsx("Desktop/VEL/clinicalValues_WP10_bariatric_lipo_Vel.xlsx", 
                         sheetIndex = 1)

df1 = merge(wp10_VEL, wp10_VEL_wrong, by.x = "Barcode", by.y = "ID")

# SNP = rs566629828

vel_het = c('S01HRB','S019YT','S001GV', 'S00278', 'S00PWE','S00U3F','S00WRX') # these 2 samples are heterozygous in the vel cohort

wdf=df1

set.seed(1)
# convert class to numeric where possible
for(i in 1:length(colnames(wdf))){
  if ( sum(is.na(as.numeric(as.character(wdf[,i])))) == dim(wdf)[1] ) {
    message(colnames(wdf)[i], ' not converted')
  } else {
    wdf[,i] = as.numeric(as.character(wdf[,i]))
  }
}

str(wdf)
summary(wdf)
skim(wdf)

wdf_sbt = wdf

wdf_sbt$gt = 'ref' 
wdf_sbt$gt[which(wdf_sbt$cohort == "VEL")] = "hom"
wdf_sbt$gt[which(wdf_sbt$Barcode %in% vel_het)] = "het"
wdf_sbt$gt[which(wdf_sbt$Barcode %in% 'S005VM')] = "hom"
wdf_sbt$cohort[which(wdf_sbt$cohort %in% 'S005VM')] = "VEL"

# wdf_sbt = wdf_sbt %>% select("gt", colnames(wdf_sbt)[-55]) 

# FBC counts are not present in the VEL data.
columns_to_remove = c("WBC", "RBC", "HGB", "MCV", "MCH",     
                      "MCHC", "PLT", "RDW_SD", "PDW",     
                      "MPV", "NEUT", "LYMPH", "MONO",    
                      "EO",  "BASO", "IG",  "RET_He",  
                      "PLT_I", "WBC_D",
                      'LDL', 'DLK1', 'WEIGHT',
                      'Free_T4')

wdf_sbt_2 = wdf_sbt[, !(colnames(wdf_sbt) %in% columns_to_remove)]

wdf_sbt_2$LDL = ( wdf_sbt_2$CHOL - wdf_sbt_2$HDL - wdf_sbt_2$TG ) / 5 # calculates LDL

ggpairs(wdf_sbt_2, 
        columns = colnames(wdf_sbt_2[,-2]), )

wdf_sbt_2$gt <- relevel(as.factor(wdf_sbt_2$gt), ref = "ref")

wdf_sbt_2 = wdf_sbt_2 %>% 
  filter(gt != "het") %>% 
  filter(COHORT == "DonorWP10" | COHORT == "VEL")

gender = c()
feature = c()
effect = c()
cohort = c()
CI_low = c()
CI_high = c()
p.val = c()
bmi = c()

for (round_bmi in 1:2) {
  if (round_bmi == 1 ) {
    for (gen in c('Male','Female')) {
      tmp_df = wdf_sbt_2[wdf_sbt_2$GENDER == gen,]
      for(i in c(22,24,25)){
      #for(i in 9:length(colnames(tmp_df))){
        model1 = lm( data = tmp_df, tmp_df[,i] ~ gt + AGE + BMI) # + Smoker + Smoker_past + Alcohol )
        param = standardize_parameters(model1, method = "posthoc", robust = TRUE)
        feature = c(feature, colnames(wdf_sbt_2)[i])
        effect = c(effect, param$Std_Coefficient[2])
        cohort = c(cohort, param$Parameter[2])
        CI_low = c(CI_low, param$CI_low[2])
        CI_high = c(CI_high, param$CI_high[2])
        param = parameters::model_parameters(model1)
        p.val = c(p.val, param$p[2])
        gender = c(gender, gen)
        bmi = c(bmi, 'BMI_as_covariates')
      }
    }
    for(i in c(22,24,25)){
#    for(i in 9:length(colnames(tmp_df))){
      tmp_var = ( wdf_sbt_2[,i] - mean(wdf_sbt_2[,i], na.rm = T) ) / sd(wdf_sbt_2[,i], na.rm = T)
      model1 = lm( data = wdf_sbt_2, tmp_var ~ gt + AGE + BMI + GENDER) # + Smoker + Smoker_past + Alcohol )
      param = parameters::parameters(model1)
      feature = c(feature, colnames(wdf_sbt_2)[i])
      effect = c(effect, param$Coefficient[2])
      cohort = c(cohort, param$Parameter[2])
      CI_low = c(CI_low, param$CI_low[2])
      CI_high = c(CI_high, param$CI_high[2])
      p.val = c(p.val, param$p[2])
      gender = c(gender, 'not_gender_stratified, but as covariates')
      bmi = c(bmi, 'BMI_as_covariates')
    }
  } else {
    for (gen in c('Male','Female')) {
      tmp_df = wdf_sbt_2[wdf_sbt_2$GENDER == gen,]
      for(i in c(22,24,25)){
#      for(i in 9:length(colnames(tmp_df))){
        model1 = lm( data = tmp_df, tmp_df[,i] ~ gt + AGE)
        param = standardize_parameters(model1, method = "posthoc", robust = TRUE)
        feature = c(feature, colnames(tmp_df)[i])
        effect = c(effect, param$Std_Coefficient[2])
        cohort = c(cohort, param$Parameter[2])
        CI_low = c(CI_low, param$CI_low[2])
        CI_high = c(CI_high, param$CI_high[2])
        param = parameters::model_parameters(model1)
        p.val = c(p.val, param$p[2])
        gender = c(gender, gen)
        bmi = c(bmi, 'NOT_BMI_as_covariates')
      }
    }
    for(i in c(22,24,25)){
#    for(i in 9:length(colnames(tmp_df))){
      tmp_var = ( wdf_sbt_2[,i] - mean(wdf_sbt_2[,i], na.rm = T) ) / sd(wdf_sbt_2[,i], na.rm = T)
      model1 = lm( data = wdf_sbt_2, tmp_var ~ gt + AGE + GENDER) 
      param = parameters::parameters(model1)
      feature = c(feature, colnames(wdf_sbt_2)[i])
      effect = c(effect, param$Coefficient[2])
      cohort = c(cohort, param$Parameter[2])
      CI_low = c(CI_low, param$CI_low[2])
      CI_high = c(CI_high, param$CI_high[2])
      p.val = c(p.val, param$p[2])
      gender = c(gender, 'not_gender_stratified, but as covariates')
      bmi = c(bmi, 'NOT_BMI_as_covariates')
    }
  }
}

stat_df = data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  'CI_low' = CI_low, 
  'CI_high' = CI_high,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr"),
  'gender' = gender,
  'bmi' = bmi
) %>%  arrange(p.val)
xlsx::write.xlsx(x = stat_df,
                 file = 'Desktop/VEL/result_wp10_vel/summary_stats_repeated_lm_wgs10_vel_with_het_20210513.xls')


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

df_to_use = wdf_sbt_2
for (feata in colnames(df_to_use)[9:length(colnames(df_to_use))]) {
  if (is.numeric(df_to_use[,feata])){
    png(filename = paste0("Desktop/VEL/result_wp10_vel/rainplot_", 'feata', "20210513.png" ), 
        bg = "transparent",
        width = 1024, 
        height = 768 )
    p1 = ggplot(data = df_to_use, aes( x = gt, y = df_to_use[,feata], fill = gt)) +
      ggtitle(label = feata) +
      geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, bw = 1) +
      geom_point(aes(color = gt), position = position_jitter(width = .15), size = .3, alpha = 0.5) +
      geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
      #expand_limits(x = 3.25) +
      guides(fill = FALSE) +
      guides(color = FALSE) +
      scale_color_manual(values = ggsci::pal_npg("nrc")(3) ) +
      scale_fill_manual(values = ggsci::pal_npg("nrc")(3) ) +
      theme_minimal() +
      facet_wrap('GENDER') 
    ggsave(filename = paste0("rainplot_", feata, "20210513.png" ),
           path = 'Desktop/VEL/result_wp10_vel/' ,
           device = 'png', plot = p1)
    dev.off()
  } else {
    next()
    ggsave(filename = paste0("boxplot_", feata, "20210513.svg" ),
           path = 'Desktop/VEL/result_wp10_vel/' ,
           device = 'svg', plot = p2)
    dev.off()
  }
}


