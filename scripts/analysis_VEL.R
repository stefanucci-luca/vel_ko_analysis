# lod the libraries used for the analysis
library('lm.beta')
library("questionr")
library('esc')
library('effectsize')
library('tidyverse')
library('see')
library('GGally')
library('RColorBrewer')

# set seed for reproducibility
set.seed(1)

######################################################
#                   _______     ___    __________ 
#    _      ______ <  / __ \   /_/ |  / / ____/ / 
#   | | /| / / __ \/ / / / / /_/ | | / / __/ / /  
#   | |/ |/ / /_/ / / /_/ //_/   | |/ / /___/ /___
#   |__/|__/ .___/_/\____/_/     |___/_____/_____/
#         /_/                                     
#####################################################

# SNP = rs566629828

# Import VEL phenotype data from the donors cohort

wdf = xlsx::read.xlsx('../data/result_wp10_vel/data_NIHR.xlsx', sheetIndex = 1)

#summary(df)
# skimr::skim(df)

# When importing the class of certain numberic valued became factor
# convert class back to numeric
for(i in 1:length(colnames(wdf))){
  # check if the column can coerce to numeric
  if ( sum(is.na(as.numeric(as.character(wdf[,i])))) == dim(wdf)[1] ) {
    # if not possible don't convert and print the message
    message(colnames(wdf)[i], ' not converted')
  } else {
    # convert to numeric
    wdf[,i] = as.numeric(as.character(wdf[,i]))
  }
}

# explore the df
# str(wdf)
# summary(wdf)

#### Start to control from here

wdf_sbt$gt = 'ref' 
wdf_sbt$gt[which(wdf_sbt$COHORT == "VEL")] = "hom"
wdf_sbt$gt[which(wdf_sbt$ID %in% vel_het)] = "het"
wdf_sbt$gt[which(wdf_sbt$ID %in% 'S005VM')] = "hom"
wdf_sbt$COHORT[which(wdf_sbt$ID %in% 'S005VM')] = "VEL"

wdf_sbt = wdf_sbt %>% select("gt", colnames(wdf_sbt)[-55]) 
  
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

xlsx::write.xlsx(wdf_sbt_2, "/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/Tables/To_CONTROL_tableS5.xlsx")

ggpairs(wdf_sbt_2, 
        columns = colnames(wdf_sbt_2[,-2]), )

gender = c()
feature = c()
effect = c()
cohort = c()
CI_low = c()
CI_high = c()
p.val = c()
bmi = c()

wdf_sbt_2$gt <- relevel(as.factor(wdf_sbt_2$gt), ref = "ref")

for (round_bmi in 1:2) {
    if (round_bmi == 1 ) {
        for (gen in c('Male','Female')) {
          tmp_df = wdf_sbt_2[wdf_sbt_2$GENDER == gen,]
          for(i in 9:length(colnames(tmp_df))){
            model1 = lm( data = tmp_df, tmp_df[,i] ~ gt + AGE + BMI) # + Smoker + Smoker_past + Alcohol )
            param = standardize_parameters(model1, method = "posthoc", robust = TRUE)
            feature = c(feature, colnames(tmp_df)[i])
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
        for(i in 9:length(colnames(tmp_df))){
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
        for(i in 9:length(colnames(tmp_df))){
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
      for(i in 9:length(colnames(tmp_df))){
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
                 file = 'ls760@platgenbio.net - Google Drive/My Drive/Desktop_Macbook_PhD/Desktop/VEL/result_wp10_vel/summary_stats_lm_wgs10_vel_with_het_20210809.xls')

# Forest plot
p <- ggplot(data=stat_df, aes(y= reorder (feature, -p.val) , x=effect, xmin=CI_low, xmax=CI_high))+ 
  #this adds the effect sizes to the plot
  geom_point()+ 
  #adds the CIs
  geom_errorbarh(height=.1)+
  #adding a vertical line at the effect = 0 mark
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)+
  #faceting based on my subgroups
  facet_grid(gender~., scales= "free", space="free")+
  #thematic stuff
  ggtitle("Target Effects")+
  theme_minimal()+
  theme(text=element_text(family="Arial",size=18, color="black"))+
  theme(panel.spacing = unit(1, "lines"))
ggsave('Desktop/VEL/result_wp10_vel/OR_picture_wgs_vel.svg' ,
       device = 'svg', plot = p)

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


plot_value = c( "LEPT", 'LAR', 'AT_FFA', 'Total_T3', 'Total_T4','ALT', 'AST', 'hsCRP', 'Ferritin' )

for (feata in colnames(wdf_sbt)[which(colnames(wdf_sbt) %in% plot_value)]) {
  p1 = ggplot(data = wdf_sbt, aes( x = COHORT, y = wdf_sbt[,feata], fill = COHORT)) +
    ggtitle(label = feata) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, bw = 0.1) +
    geom_point(aes(color = COHORT), position = position_jitter(width = .15), size = .3, alpha = 0.8) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
    #expand_limits(x = 3.25) +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    scale_color_manual(values = c("#4DBBD5FF" ,"#E64B35FF")) +
    scale_fill_manual(values = c("#4DBBD5FF" ,"#E64B35FF")) +
    theme_minimal() # +
  # facet_wrap('sex') 
  ggsave(filename = paste0("rainplot_paper_", feata, "_paper_20210809.png" ),
         path = 'ls760@platgenbio.net - Google Drive/My Drive/Desktop_Macbook_PhD/Desktop/VEL/result_wp10_vel/' , dpi = "retina", width = 10, height = 20, units = "cm",
         device = 'png', plot = p1)
} 

  


############################################################
#  ____       _          _            _   
# |__ /_ _ __| |  __ ___| |_  ___ _ _| |_ 
# |_  \ '_/ _` | / _/ _ \ ' \/ _ \ '_|  _|
# |___/_| \__,_| \__\___/_||_\___/_|  \__|
#############################################################

#control df 
df_ctrl3 = xlsx::read.xlsx("Desktop/VEL/3rd_cohort/All individuals on 149 (control data).xlsx", sheetIndex = 1)
df_ctrl3$status = 'control'

df_ctrl3 = df_ctrl3 %>% filter(Ethnicity=="White")


#case df 
df_case3 = xlsx::read.xlsx("Desktop/VEL/3rd_cohort/VEL NEG Protocol 7 subjects.xlsx", sheetIndex = 4) 
df_case3 = df_case3[1:13,]
df_case3$status = 'vel_negative'
colnames(df_case3) = str_replace(string = colnames(df_case3), pattern = 'Sex', replacement = 'Gender' )
colnames(df_case3) = str_replace(string = colnames(df_case3), pattern = 'Age.at.study', replacement = 'Age' )
colnames(df_case3) = str_replace(string = colnames(df_case3), pattern = 'NEFA.umol.L', replacement = 'NEFA' )
colnames(df_case3) = str_replace(string = colnames(df_case3), pattern = 'Leptin.ng.ml', replacement = 'Leptin' )
colnames(df_case3) = str_replace(string = colnames(df_case3), pattern = 'Adiponectin.ug.ml', replacement = 'Adiponectin' )
colnames(df_case3) = str_replace(string = colnames(df_case3), pattern = 'Insulin.pmol.L', replacement = 'Insulin' )
colnames(df_case3) = str_replace(string = colnames(df_case3), pattern = 'Plasma.Glucose.mmol.l.', replacement = 'Plasma.Glucose.mmol.l' )
df_case3$Age = as.numeric(as.character(df_case3$Age))

df_case3_2 = df_case3[,6:dim(df_case3)[2]-1] %>% mutate_all(as.character) %>% mutate_all(as.double)
df_case3_2 = cbind(df_case3[,c(1:4,44)], df_case3_2)
df_case3_2$Ethnicity = 'White'

df_ctrl3_2 = df_ctrl3[,5:dim(df_ctrl3)[2]-1] %>% mutate_all(as.character) %>% mutate_all(as.numeric) 
df_ctrl3_2 = cbind(df_ctrl3[,c(1:3,dim(df_ctrl3)[2])], df_ctrl3_2)

# merge the 2 reasources
df_3rd = df_case3_2 %>%
  full_join(df_ctrl3_2, by = intersect(colnames(df_case3_2), colnames(df_ctrl3_2)))

# vector for stats
gender = c()
feature = c()
effect = c()
cohort = c()
sd = c()
p.val = c()
CI_low = c()
CI_high = c()

# Not gender stratified
tmp_df = df_3rd
  for(i in 6:length(colnames(tmp_df))){
    if ( (sum(is.na(tmp_df[which(tmp_df$status == 'vel_negative'),i])) != length( tmp_df[which(tmp_df$status == 'vel_negative'),i] ) &&
          sum(is.na(tmp_df[which(tmp_df$status == 'control'),i])) != length( tmp_df[which(tmp_df$status == 'control'),i] ) ) &&
         class(tmp_df[,i]) != "character" 
    ) {
      normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE)
      model1 = lm( data = tmp_df, normalise_tmp ~ status + Age) # + BMI + Smoker + Smoker_past + Alcohol )
      param = parameters::parameters(model1, method = "posthoc", robust = TRUE)
      feature = c(feature, colnames(tmp_df)[i])
      effect = c(effect, param$Coefficient[2])
      cohort = c(cohort, param$Parameter[2])
      CI_low = c(CI_low, param$CI_low[2])
      CI_high = c(CI_high, param$CI_high[2])
      param = parameters::model_parameters(model1)
      p.val = c(p.val, param$p[2])
      gender = c(gender, 'Not stratified')
    }
  }


# Gender stratified
for (gen in c('Male','Female')) {
  tmp_df = df_3rd[df_3rd$Gender == gen,]
  for(i in 6:length(colnames(tmp_df))){
    if ( (sum(is.na(tmp_df[which(tmp_df$status == 'vel_negative'),i])) != length( tmp_df[which(tmp_df$status == 'vel_negative'),i] ) &&
         sum(is.na(tmp_df[which(tmp_df$status == 'control'),i])) != length( tmp_df[which(tmp_df$status == 'control'),i] ) ) &&
      class(tmp_df[,i]) != "character" 
         ) {
    normalise_tmp = ( tmp_df[,i] - mean(as.vector(tmp_df[,i]), na.rm = TRUE) ) / sd(tmp_df[,i], na.rm = TRUE)
    model1 = lm( data = tmp_df, normalise_tmp ~ status + Age) # + BMI + Smoker + Smoker_past + Alcohol )
    param = parameters::parameters(model1, method = "posthoc", robust = TRUE)
    feature = c(feature, colnames(tmp_df)[i])
    effect = c(effect, param$Coefficient[2])
    cohort = c(cohort, param$Parameter[2])
    CI_low = c(CI_low, param$CI_low[2])
    CI_high = c(CI_high, param$CI_high[2])
    param = parameters::model_parameters(model1)
    p.val = c(p.val, param$p[2])
    gender = c(gender, gen)
    }
  }
}

# Build the df for the statistical analysis
stat_ukb = data.frame(
  'feature' = feature,
  'effect' = effect, 
  'cohort' = cohort, 
  #'sd' = sd,
  'p.val' = p.val,
  'p.val.FDR' = p.adjust(p.val, method = "fdr"),
  'gender' = gender
) %>%  
  arrange(p.val) %>%
  filter(cohort != "(Intercept)") 

write.csv(stat_ukb,file = 'Desktop/VEL/result_3rd_cohort_vel/stats_table.csv')

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


df_to_use = df_3rd
for (feata in colnames(df_to_use)[6:length(colnames(df_to_use))]) {
  if (is.numeric(df_to_use[,feata])){
    p1 = ggplot(data = df_to_use, aes( x = status, y = df_to_use[,feata], fill = status)) +
      ggtitle(label = feata) +
      geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, bw = 1) +
      geom_point(aes(color = status), position = position_jitter(width = .15), size = .3, alpha = 0.5) +
      geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
      #expand_limits(x = 3.25) +
      guides(fill = FALSE) +
      guides(color = FALSE) +
      scale_color_manual(values = ggsci::pal_npg("nrc")(3) ) +
      scale_fill_manual(values = ggsci::pal_npg("nrc")(3) ) +
      theme_minimal() +
      facet_wrap('Gender') 
    ggsave(filename = paste0("rainplot_", feata, ".png" ),
           path = 'Desktop/VEL/result_3rd_cohort_vel/' ,
           device = 'png', plot = p1)
  } else {
    p2 = ggplot(as.data.frame(table(df_to_use[,feata])),
                aes(x=Var1, y=Freq)) +
      geom_col()  + 
      ggtitle(label = feata)
    ggsave(filename = paste0("boxplot_", feata, ".svg" ),
           path = 'Desktop/VEL/result_3rd_cohort_vel/' ,
           device = 'svg', plot = p2)
  }
}


###########################################################
#          __      __    _       __                __  
#   __  __/ /__   / /_  (_)___  / /_  ____ _____  / /__
#  / / / / //_/  / __ \/ / __ \/ __ \/ __ `/ __ \/ //_/
# / /_/ / ,<    / /_/ / / /_/ / /_/ / /_/ / / / / ,<   
# \__,_/_/|_|  /_.___/_/\____/_.___/\__,_/_/ /_/_/|_|  
# 
############################################################

# analysis of UKB data - extracted with the script commented below
df_large = data.table::fread("/Volumes/GoogleDrive-105684671225474146017/My Drive/Phd/Shared Luca Mattia/UKB_VEL_relevant_fields_20210308.csv")
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
eur_unre_500k = data.table::fread('/Users/ls31/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/british-ancestry-MSUP_unrel_europ_500K_from_will_20210317.tsv') %>% 
  select(ID)

#Covariates
load('/Users/ls31/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luanluan/data/ukb_dvt_pe_combined.rdata')
# times in the hospital
df_final = merge(df, df_large, by = 'eid' )

df_times_hosp = data.table::fread('/Users/ls31/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/ukbcvo/times_in_hopital_otUKBcnv.tab', 
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

# DF that has:
# Unrelated european
# WT
# Hom
df_final_un_eu_only_hom = df_final_un_eu[which(df_final_un_eu$genotype == '0' | df_final_un_eu$genotype == '2'),]

df_final_un_eu$ldl_indirect_calculation = ( df_final_un_eu$cholesterol_f30690_0_0 - df_final_un_eu$hdl_cholesterol_f30760_0_0 - df_final_un_eu$triglycerides_f30870_0_0 ) / 5

std_direct_ldl = (df_final_un_eu$ldl_direct_f30780_0_0 - mean(df_final_un_eu$ldl_direct_f30780_0_0, na.rm = T))/ sd (df_final_un_eu$ldl_direct_f30780_0_0, na.rm = T)

# plot correlation between direct and indirect LDL
plot(df_final_un_eu$ldl_indirect_calculation, std_direct_ldl)

# df that has to be used in the statistical analysis below
df_to_use = df_final_un_eu
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
          for(i in c(23,27:length(colnames(tmp_df)))){
            # normalise_tmp = tmp_df[,i] # no normalisation
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
            for(i in c(23,27,30:length(colnames(tmp_df)))){
              # normalise_tmp = tmp_df[,i] # no normalisation
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
    for(i in c(23,27:length(colnames(df_to_use)))){
              normalise_tmp = ( df_to_use[,i] - mean(as.vector(df_to_use[,i]), na.rm = TRUE) ) / sd(df_to_use[,i], na.rm = TRUE) 
              model1 = lm(normalise_tmp ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex_f31_0_0 ) 
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
      tmp_df = df_to_use[which(df_to_use$sex_f31_0_0 == gen),]
      if (gen == 0){
        for(i in c(23,27:length(colnames(tmp_df)))){
          # normalise_tmp = tmp_df[,i] # no normalisation
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
        for(i in c(23,27,30:length(colnames(tmp_df)))){
          # normalise_tmp = tmp_df[,i] # no normalisation
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
    for(i in c(23,27:length(colnames(tmp_df)))){
      normalise_tmp = ( df_to_use[,i] - mean(as.vector(df_to_use[,i]), na.rm = TRUE) ) / sd(df_to_use[,i], na.rm = TRUE) 
      model1 = lm(normalise_tmp ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex_f31_0_0 + df_to_use$bmi ) 
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
  filter(cohort != "(Intercept)") %>% 
  filter(feature != "statin" ) %>% filter( feature != "vte"  )

write.table(stat_ukb, file = 'Desktop/VEL/result_ukb_vel/ukb_stats_table_20210325.csv', row.names = F, append = F, quote = F, sep = ",")

################################################
#          __          __   __        
#  _______/  |______ _/  |_|__| ____  
# /  ___/\   __\__  \\   __\  |/    \ 
# \___ \  |  |  / __ \|  | |  |   |  \
# /____  > |__| (____  /__| |__|___|  /
#      \/            \/             \/ 
#################################################

df_to_use$statin = ifelse(df_to_use$statin != 0, yes = 1, no = 0)

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
          model1 = glm(tmp_df$statin ~ factor(tmp_df$genotype) + tmp_df$ages ) 
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
            feature = c(feature,  'statin' )
            gender = c(gender, gen)
            cohort = c(cohort,param$Parameter[num])
            effect = c(effect, param$Coefficient[num])
            sd = c(sd, param$SE[num])
            CI_low = c(CI_low, param$CI_low[num] )
            CI_high = c(CI_high, param$CI_high[num] )
            p.val = c(p.val, param$p[num])
            bmi = c(bmi, 'NOT_BMI_as_covariates')
          }
      } else {
          model1 = glm(tmp_df$statin ~ factor(tmp_df$genotype) + tmp_df$ages ) 
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
            feature = c(feature,  'statin' )
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
      model1 = glm(df_to_use$statin ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex_f31_0_0 ) 
      param = parameters::model_parameters(model1)
      for (num in 1:length(unique(df_to_use$genotype))) {
        feature = c(feature, 'statin' )
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
   else {
    for (gen in 0:1) {
      tmp_df = df_to_use[which(df_to_use$sex_f31_0_0 == gen),]
      if (gen == 0){
          model1 = glm(tmp_df$statin ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi )
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
            feature = c(feature, 'statin' )
            gender = c(gender, gen)
            cohort = c(cohort,param$Parameter[num])
            effect = c(effect, param$Coefficient[num])
            sd = c(sd, param$SE[num])
            CI_low = c(CI_low, param$CI_low[num] )
            CI_high = c(CI_high, param$CI_high[num] )
            p.val = c(p.val, param$p[num])
            bmi = c(bmi, 'BMI_as_covariates')
          }
      } else {
        model1 = glm(tmp_df$statin ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi )
        param = parameters::model_parameters(model1)
        for (num in 1:length(unique(tmp_df$genotype))) {
          feature = c(feature, 'statin' )
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
      model1 = glm(df_to_use$statin ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex_f31_0_0 + df_to_use$bmi ) 
      param = parameters::model_parameters(model1)
      for (num in 1:length(unique(df_to_use$genotype))) {
        feature = c(feature, 'statin')
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

# Build the df for the statistical analysis
stat_ukb_statin =  data.frame(
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

write.table(stat_ukb_statin, file = 'Desktop/VEL/result_ukb_vel/ukb_stats_table_20210325.csv', row.names = F, append = T, quote = F, sep = ",")

##########__#######################################     
#  _   __/ /____ 
# | | / / __/ _ \
# | |/ / /_/  __/
# |_|___/\__/\___/ 
# #################################################

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
        model1 = glm(tmp_df$vte ~ factor(tmp_df$genotype) + tmp_df$ages ) 
        param = parameters::model_parameters(model1)
        for (num in 1:length(unique(tmp_df$genotype))) {
          feature = c(feature,  'vte' )
          gender = c(gender, gen)
          cohort = c(cohort,param$Parameter[num])
          effect = c(effect, param$Coefficient[num])
          sd = c(sd, param$SE[num])
          CI_low = c(CI_low, param$CI_low[num] )
          CI_high = c(CI_high, param$CI_high[num] )
          p.val = c(p.val, param$p[num])
          bmi = c(bmi, 'NOT_BMI_as_covariates')
        }
      } else {
        model1 = glm(tmp_df$vte ~ factor(tmp_df$genotype) + tmp_df$ages ) 
        param = parameters::model_parameters(model1)
        for (num in 1:length(unique(tmp_df$genotype))) {
          feature = c(feature,  'vte' )
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
    model1 = glm(df_to_use$vte ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex_f31_0_0 ) 
    param = parameters::model_parameters(model1)
    for (num in 1:length(unique(df_to_use$genotype))) {
      feature = c(feature, 'vte' )
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
  else {
    for (gen in 0:1) {
      tmp_df = df_to_use[which(df_to_use$sex_f31_0_0 == gen),]
      if (gen == 0){
        model1 = glm(tmp_df$vte ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi )
        param = parameters::model_parameters(model1)
        for (num in 1:length(unique(tmp_df$genotype))) {
          feature = c(feature, 'vte' )
          gender = c(gender, gen)
          cohort = c(cohort,param$Parameter[num])
          effect = c(effect, param$Coefficient[num])
          sd = c(sd, param$SE[num])
          CI_low = c(CI_low, param$CI_low[num] )
          CI_high = c(CI_high, param$CI_high[num] )
          p.val = c(p.val, param$p[num])
          bmi = c(bmi, 'BMI_as_covariates')
        }
      } else {
        model1 = glm(tmp_df$vte ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi )
        param = parameters::model_parameters(model1)
        for (num in 1:length(unique(tmp_df$genotype))) {
          feature = c(feature, 'vte' )
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
    model1 = glm(df_to_use$vte ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex_f31_0_0 + df_to_use$bmi ) 
    param = parameters::model_parameters(model1)
    for (num in 1:length(unique(df_to_use$genotype))) {
      feature = c(feature, 'vte')
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

# Build the df for the statistical analysis
stat_ukb_vte =  data.frame(
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

#################################################################
#  _                 _   _                         _            
# | | _____   _____ | |_| |__  _   _ _ __ _____  _(_)_ __   ___ 
# | |/ _ \ \ / / _ \| __| '_ \| | | | '__/ _ \ \/ / | '_ \ / _ \
# | |  __/\ V / (_) | |_| | | | |_| | | | (_) >  <| | | | |  __/
# |_|\___| \_/ \___/ \__|_| |_|\__, |_|  \___/_/\_\_|_| |_|\___|
#                              |___/                            
#################################################################

df_levo = data.table::fread('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/VEL/UKB_VEL_relevant_fields_20210324.csv') %>% 
  select(- 'V1' ) %>% 
  select('V1','levotitoxine_yes_no')
df_levo = as.data.frame(df_levo)
df_final_un_eu = merge(df_final_un_eu, df_levo, by.x = 'eid', by.y = "V1", all.x = T)

# df that has to be used in the statistical analysis below
df_to_use = df_final_un_eu

df_to_use$genotype <- relevel(as.factor(df_to_use$genotype), ref = '0')

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
      tmp_df = df_to_use[which(df_to_use$levotitoxine_yes_no == gen),]
      if (gen == 0){
        model1 = glm(tmp_df$levotitoxine_yes_no ~ factor(tmp_df$genotype) + tmp_df$ages ) 
        param = parameters::model_parameters(model1)
        for (num in 1:length(unique(tmp_df$genotype))) {
          feature = c(feature,  'levotitoxine_yes_no' )
          gender = c(gender, gen)
          cohort = c(cohort,param$Parameter[num])
          effect = c(effect, param$Coefficient[num])
          sd = c(sd, param$SE[num])
          CI_low = c(CI_low, param$CI_low[num] )
          CI_high = c(CI_high, param$CI_high[num] )
          p.val = c(p.val, param$p[num])
          bmi = c(bmi, 'NOT_BMI_as_covariates')
        }
      } else {
        model1 = glm(tmp_df$levotitoxine_yes_no ~ factor(tmp_df$genotype) + tmp_df$ages ) 
        param = parameters::model_parameters(model1)
        for (num in 1:length(unique(tmp_df$genotype))) {
          feature = c(feature,  'levotitoxine_yes_no' )
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
    model1 = glm(df_to_use$levotitoxine_yes_no ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex_f31_0_0 ) 
    param = parameters::model_parameters(model1)
    for (num in 1:length(unique(df_to_use$genotype))) {
      feature = c(feature, 'levotitoxine_yes_no' )
      gender = c(gender, 'not_gender_stratified, but as covariates')
      cohort = c(cohort,param$Parameter[num])
      effect = c(effect, param$Coefficient[num])
      sd = c(sd, param$SE[num])
      CI_low = c(CI_low, param$CI_low[num] )
      CI_high = c(CI_high, param$CI_high[num] )
      p.val = c(p.val, param$p[num])
      bmi = c(bmi, 'NOT_BMI_as_covariates')
    }
  }
  else {
    for (gen in 0:1) {
      tmp_df = df_to_use[which(df_to_use$sex_f31_0_0 == gen),]
      if (gen == 0){
        model1 = glm(tmp_df$levotitoxine_yes_no ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi )
        param = parameters::model_parameters(model1)
        for (num in 1:length(unique(tmp_df$genotype))) {
          feature = c(feature, 'levotitoxine_yes_no' )
          gender = c(gender, gen)
          cohort = c(cohort,param$Parameter[num])
          effect = c(effect, param$Coefficient[num])
          sd = c(sd, param$SE[num])
          CI_low = c(CI_low, param$CI_low[num] )
          CI_high = c(CI_high, param$CI_high[num] )
          p.val = c(p.val, param$p[num])
          bmi = c(bmi, 'BMI_as_covariates')
        }
      } else {
        model1 = glm(tmp_df$levotitoxine_yes_no ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi )
        param = parameters::model_parameters(model1)
        for (num in 1:length(unique(tmp_df$genotype))) {
          feature = c(feature, 'levotitoxine_yes_no' )
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
    model1 = glm(df_to_use$levotitoxine_yes_no ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex_f31_0_0 + df_to_use$bmi ) 
    param = parameters::model_parameters(model1)
    for (num in 1:length(unique(df_to_use$genotype))) {
      feature = c(feature, 'levotitoxine_yes_no')
      gender = c(gender, 'not_gender_stratified, but as covariates')
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

# Build the df for the statistical analysis
stat_ukb_levo =  data.frame(
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

write.table(stat_ukb_levo, file = 'Desktop/VEL/result_ukb_vel/ukb_stats_table_20210325.csv', row.names = F, append = T, quote = F, sep = ",")


# All phenotypes 

load('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luca/returned_datasets/ukb_first_events.rdata')
conv_tab = read.csv('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luca/returned_datasets/phenotype_event_names.csv')

df_to_use_pheno = merge(df_to_use, ukb_first_events, by.x = "eid", "identifier")

pheno_interest = which(startsWith(colnames(df_to_use_pheno), "cal_p"))

for (var in colnames(df_to_use_pheno[,pheno_interest])) {
  df_to_use_pheno[,var] = ifelse(is.na(df_to_use_pheno[,var]), 0, 1)  
}

gender = c()
feature = c()
effect = c()
cohort = c()
sd = c()
CI_low = c()
CI_high = c()
p.val = c()
bmi = c()

df_to_use = df_to_use_pheno


for(i in which(startsWith(colnames(df_to_use_pheno), "cal_p")) ){
  print(i)
  normalise_tmp = df_to_use_pheno[,i] # no normalisation
  model1 = glm(normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages + df_to_use_pheno$sex_f31_0_0 + df_to_use_pheno$bmi) 
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

# Build the df for the statistical analysis
stat_ukb_pheno =  data.frame(
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


# All phenotypes GP

load('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luca/returned_datasets/ukb_event_clinics.RData')

df_to_use_pheno = merge(df_to_use, gp_pheno, by = "eid")

pheno_interest = which(startsWith(colnames(df_to_use_pheno), "cal_p"))

for (var in colnames(df_to_use_pheno[,pheno_interest])) {
  df_to_use_pheno[,var] = ifelse(is.na(df_to_use_pheno[,var]), 0, 1)  
}

gender = c()
feature = c()
effect = c()
cohort = c()
sd = c()
CI_low = c()
CI_high = c()
p.val = c()
bmi = c()

df_to_use = df_to_use_pheno


for(i in grep(x = colnames(df_to_use_pheno), pattern = "*to_RC$") ){
  print(i)
  normalise_tmp = df_to_use_pheno[,i] # no normalisation
  model1 = glm(normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages + df_to_use_pheno$sex_f31_0_0 )# + df_to_use_pheno$bmi) 
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

# Build the df for the statistical analysis
stat_ukb_pheno =  data.frame(
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

#######
#######
####### All phenotypes GP + HES
#######
#######

conv_vec = conv_tab$phenotype_short
names(conv_vec) = paste0("cal_ps_d_", sprintf("%03d", conv_tab$event_code))

for (value in grep(pattern = "cal_ps_d_", x = colnames(ukb_first_events), value = T)) {
  position = grep(pattern = value, names(conv_vec) )
  replace = grep(pattern = value, colnames(ukb_first_events) )
  colnames(ukb_first_events)[replace] = as.character(conv_vec[position])
}

for (col_name in colnames(ukb_first_events)[2:length(colnames(ukb_first_events))]) {
  ukb_first_events[,col_name] = ifelse(is.na(ukb_first_events[,col_name]) == TRUE, FALSE, TRUE)
}

setting_df_1 = merge(eur_unre_500k, gp_pheno, by.x = "ID", by.y = "eid", all.x = T)
setting_df_2 = merge(eur_unre_500k, ukb_first_events, by.y = "identifier", by.x = "ID", all.x = T)
setting_df = merge(setting_df_1, setting_df_2, by = "ID")

setting_df = as.data.frame(setting_df[20:dim(setting_df)[1],])

for (column in colnames(setting_df)) {
  setting_df[,column] = setting_df[,column] %>% replace_na(FALSE) 
}

df_tmp = data.frame("eid" = setting_df$ID)
for (phen in conv_tab$phenotype_short) {
  columns = grep(pattern = phen, x = colnames(setting_df))
  pheno_vec = apply(setting_df[,columns], 
                    1,
                    function(x)
                      sum(x))
  df_tmp[,phen] = pheno_vec 
}

df_pheno_hos_clin = df_tmp

for (column in colnames(df_pheno_hos_clin)[2:dim(df_pheno_hos_clin)[2]]) {
  df_pheno_hos_clin[,column] = ifelse(df_pheno_hos_clin[,column] == 0, 0, 1)
}

df_to_use_pheno = merge(df, df_pheno_hos_clin, by = "eid" )
df_to_use_pheno = merge(df_large, df_to_use_pheno, by = "eid" )

gender = c()
feature = c()
effect = c()
cohort = c()
sd = c()
CI_low = c()
CI_high = c()
p.val = c()
bmi = c()

for(i in 74:dim(df_to_use_pheno)[2]){
  print(i)
  normalise_tmp = df_to_use_pheno[,i] # no normalisation
  model1 = glm(normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages + df_to_use_pheno$sex_f31_0_0 + df_to_use_pheno$bmi, family = "binomial") 
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

for(i in 74:dim(df_to_use_pheno)[2]){
  print(i)
  normalise_tmp = df_to_use_pheno[,i] # no normalisation
  model1 = glm(normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages + df_to_use_pheno$sex_f31_0_0, family = "binomial") 
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

# Build the df for the statistical analysis
stat_ukb_pheno =  data.frame(
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

xlsx::write.xlsx(stat_ukb_pheno, file = 'Desktop/VEL/result_ukb_vel/ukb_stats_pheno_clinic_hosp_20210522.xlsx', sheetName = "VEL_ICD_UKB", col.names = T, row.names = F, append = F)

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

plot_value = c("triglycerides_f30870_0_0", "alanine_aminotransferase_f30620_0_0", "aspartate_aminotransferase_f30650_0_0", "gamma_glutamyltransferase_f30730_0_0",
               "leg_fat_mass_right_f23112_0_0", "leg_fatfree_mass_right_f23113_0_0", "leg_fat_mass_left_f23116_0_0",                     
               "leg_fatfree_mass_left_f23117_0_0", "arm_fat_mass_right_f23120_0_0", "arm_fatfree_mass_right_f23121_0_0",                    
               "arm_predicted_mass_right_f23122_0_0", "arm_fat_mass_left_f23124_0_0", "arm_fatfree_mass_left_f23125_0_0", "arm_predicted_mass_left_f23126_0_0" ,
               "urate_f30880_0_0")


plot_value = c("triglycerides_f30870_0_0", "alanine_aminotransferase_f30620_0_0", "aspartate_aminotransferase_f30650_0_0", "gamma_glutamyltransferase_f30730_0_0",
               "leg_fat_mass_right_f23112_0_0", "leg_fatfree_mass_right_f23113_0_0", "leg_fat_mass_left_f23116_0_0",                     
               "leg_fatfree_mass_left_f23117_0_0", "arm_fat_mass_right_f23120_0_0", "arm_fatfree_mass_right_f23121_0_0",                    
               "arm_predicted_mass_right_f23122_0_0", "arm_fat_mass_left_f23124_0_0", "arm_fatfree_mass_left_f23125_0_0", "arm_predicted_mass_left_f23126_0_0" ,
               "urate_f30880_0_0", "weight_f21002_0_0")

for (feata in colnames(df_to_use)[which(colnames(df_to_use) %in% plot_value)]) {
    p1 = ggplot(data = df_to_use[df_to_use$genotype %in% c(0,2),], aes( x = genotype, y = df_to_use[df_to_use$genotype %in% c(0,2),feata], fill = genotype)) +
      ggtitle(label = feata) +
      # geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, bw = 5) +
      geom_jitter(aes(color = genotype), position = position_jitter(width = .3), size = .28, alpha = 0.3) +
      geom_violin(width = .40, outlier.shape = NA, alpha = 0.5, position = position_nudge(x=0.5), draw_quantiles = c(0.25, 0.5, 0.75), trim = F) +
      #expand_limits(x = 3.25) +
      guides(fill = FALSE) +
      guides(color = FALSE) +
      scale_color_manual(values = c("#3C5488FF","#990000")) +
      scale_fill_manual(values = c("#3C5488FF","#990000")) +
      theme_minimal() + 
      #coord_cartesian(
      #ylim = c(0, quantile(df_to_use[df_to_use$genotype %in% c(0,2),feata], na.rm = T)[[4]])
      #)+
      theme(panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA))# +
      # facet_wrap('sex') 
    ggsave(filename = paste0("rainplot_paper_", feata, "_bw5_zoomed_20210609.png" ),
           path = 'Desktop/VEL/result_ukb_vel/', dpi = "retina", width = 10, height = 20, units = "cm",
           device = 'png', plot = p1)
  } 


  #######################################################################################################################################
#  _____   _       _          _  _             _    _                     _  _   __   __                                        
# |  __ \ (_)     | |        (_)| |           | |  (_)                   | |(_) / _| / _|                                       
# | |  | | _  ___ | |_  _ __  _ | |__   _   _ | |_  _   ___   _ __     __| | _ | |_ | |_  ___  _ __  ___  _ __    ___  ___  ___ 
# | |  | || |/ __|| __|| '__|| || '_ \ | | | || __|| | / _ \ | '_ \   / _` || ||  _||  _|/ _ \| '__|/ _ \| '_ \  / __|/ _ \/ __|
# | |__| || |\__ \| |_ | |   | || |_) || |_| || |_ | || (_) || | | | | (_| || || |  | | |  __/| |  |  __/| | | || (__|  __/\__ \
# |_____/ |_||___/ \__||_|   |_||_.__/  \__,_| \__||_| \___/ |_| |_|  \__,_||_||_|  |_|  \___||_|   \___||_| |_| \___|\___||___/
########################################################################################################################################


parameters::parameters(t.test(wdf_sbt_2$AGE[which(wdf_sbt_2$COHORT == "VEL")],
                              wdf_sbt_2$AGE[which(wdf_sbt_2$COHORT == "DonorWP10")]))

parameters::parameters(t.test(wdf_sbt_2$AGE[which(wdf_sbt_2$COHORT == "VEL" & wdf_sbt_2$GENDER == "Female")],
                              wdf_sbt_2$AGE[which(wdf_sbt_2$COHORT == "DonorWP10" & wdf_sbt_2$GENDER == "Female")]))

parameters::parameters(t.test(wdf_sbt_2$AGE[which(wdf_sbt_2$COHORT == "VEL" & wdf_sbt_2$GENDER == "Male")],
                              wdf_sbt_2$AGE[which(wdf_sbt_2$COHORT == "DonorWP10" & wdf_sbt_2$GENDER == "Male")])
)

parameters::parameters(t.test(df_to_use$ages[which(df_to_use$genotype == 2)],
                              df_to_use$ages[which(df_to_use$genotype == 0)]))

parameters::parameters(t.test(df_to_use$ages[which(df_to_use$genotype == 2 & df_to_use$sex == "Female")],
                              df_to_use$ages[which(df_to_use$genotype == 0 & df_to_use$sex == "Female")]))

parameters::parameters(t.test(df_to_use$ages[which(df_to_use$genotype == 2 & df_to_use$sex == "Male")],
                              df_to_use$ages[which(df_to_use$genotype == 0 & df_to_use$sex == "Male")])
)

# Age VEL
ks.test(wdf_sbt_2$AGE[which(wdf_sbt_2$COHORT == "VEL")],
                              df_to_use$ages[which(df_to_use$genotype == 2)])
# Age - not VELs
ks.test(wdf_sbt_2$AGE[which(wdf_sbt_2$COHORT == "DonorWP10")],
                              df_to_use$ages[which(df_to_use$genotype == 0)])

# BMI VEL
ks.test(wdf_sbt_2$BMI[which(wdf_sbt_2$COHORT == "VEL")],
                              df_to_use$bmi[which(df_to_use$genotype == 2)])
# BMI - not VEL
ks.test(wdf_sbt_2$BMI[which(wdf_sbt_2$COHORT == "DonorWP10")],
                              df_to_use$bmi[which(df_to_use$genotype == 0)])



# # Forest plot
# stat_ukb = na.omit(stat_ukb)
# ggplot(data=stat_ukb[which(stat_ukb$p.val.FDR < 0.05),], aes(y= reorder (feature, -p.val) , x=effect, xmin=effect-sd, xmax=effect+sd))+ 
#   #this adds the effect sizes to the plot
#   geom_point()+ 
#   #xlim(c(-1,1)) +
#   #adds the CIs
#   geom_errorbarh(height=.1)+
#   #adding a vertical line at the effect = 0 mark
#   geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)+
#   #faceting based on my subgroups
#   facet_grid(gender+cohort~., scales= "free", space="free")+
#   #thematic stuff
#   ggtitle("Target Effects")+
#   theme_minimal()+
#   theme(text=element_text(family="Arial",size=18, color="black"))+
#   theme(panel.spacing = unit(1, "lines"))

# #############
# # Analysis on UKB
# library(ukbtools)
# library(stringr)
# # On the HPC
# df_ukb = ukb_df('ukb44092', n_threads = "max")
# df_field_ukb <- ukb_df_field("ukb44092")
# 
# fields = c('f50_', 'f2714_', 'f3581_', 'f23099_', 'f23112_', 'f23113_','f23116_', 'f23117_', 'f23120_',
#             'f23121_', 'f23122_', 'f23122_', 'f23124_', 'f23125_', 'f23125_', 'f23126_', 'f30600_', 'f30610_',
#             'f30620_', 'f30630_', 'f30640_', 'f30650_', 'f30660_', 'f30670_', 'f30680_', 'f30690_', 'f30700_', 'f30710_',
#             'f30720_', 'f30730_', 'f30740_', 'f30750_', 'f30760_', 'f30770_', 'f30780_', 'f30790_', 'f30800_',
#             'f30810_', 'f30820_', 'f30830_', 'f30840_', 'f30850_', 'f30860_', 'f30870_', 'f30880_',
#             'f30890_', 'f48_', 'f21001_', 'f21002_', 'f2976_')
# 
# df_to_use = as.data.frame(as.character(df_ukb[,1]))
# for (field in fields) {
#   if(is.numeric(df_ukb[,grepl(pattern = paste0(field, '[0-9]'), x = colnames(df_ukb) )][,1]))
#     {
#   mean_v = apply(X = df_ukb[,grepl(pattern = paste0(field, '[0-9]'), x = colnames(df_ukb) )],
#         MARGIN = 1,
#         function(x)
#           mean(x, na.rm = TRUE)
#         )
#   df_to_use[,colnames(df_ukb)[grepl(pattern = paste0(field, '[0-9]'), x = colnames(df_ukb) )][1]] = mean_v
#   } else {
#     mean_v = apply(X = df_ukb[,grepl(pattern = paste0(field, '[0-9]'), x = colnames(df_ukb) )],
#                    MARGIN = 1,
#                    function(x)
#                      paste0(x, collapse = ';')
#     )
#     df_to_use[,field] = mean_v
#   }
#   }
# colnames(df_to_use)[1] = 'eid'
# 
# # Need to append the gender in a second time (field 31, 'sex_f31_0_0')
# df_to_use = merge( df_ukb[, c("eid","sex_f31_0_0")], df_to_use, by = 'eid' )
# 
# # Need to find a way to put in f20003_
# medication_columns <- grep(pattern = 'f20003_', x = colnames(df_ukb) )
# # the drugs of interest are coded in UKB as:
# # rosuvastatin = 1141192410
# # simvastatin = 1140861958
# # atorvastatin = 1141146234
# # pravastatin = 1140888648
# statin.info.vector = function(df) {
#     apply(df,
#     1,
#     function(x)
#     sum(grepl(pattern="1141192410|1140861958|1141146234|1140888648",
#     x= x))
#     )
# }
# 
# df_to_use$statin = statin.info.vector(df_ukb[,medication_columns])
# ifelse(df_to_use$statin == 0, print(0), print(1))
# # field 41234 has to go in individually
# df_to_use = merge( df_ukb[, c("eid","records_in_hes_inpatient_diagnoses_dataset_f41234_0_0")], df_to_use, by = 'eid' )
# 
# 
#This has goes to in 41235 (times in hospita)
# 
# # Age diabetes diagnose goes in individually (field 2976)
# age_diagnose_diabetes_col = c(1477, 1478, 1479, 1480)
# 
# age.diabetes = function(df) {
#     apply(df,
#     1,
#     function(x)
#     mean(x, na.rm=TRUE)
#     )
# }
# 
# df_to_use$age_diagno_diabetes = age.diabetes(df_ukb[,age_diagnose_diabetes_col])
# 
# 
# # Genotype 
# gt = data.table::fread('../VEL/UKBB_vel_genotypes.txt', sep = ' ', header = F
# )
# gt$V1 = str_split(gt$V1, pattern = "_", n = 2, simplify = TRUE)[,1]
# df_to_use$eid = as.character(df_to_use$eid)
# df_to_use = merge( gt, df_to_use, by.x = 'V1', by.y = "eid"  )
# 
# 
# To calculate the LD table I used this https://ldlink.nci.nih.gov/?tab=home. 
# Firs tool LDmatrix Tool with Population GBR
# Then to remove the one in LD I used SNPClip, R2>0.7
# thne on the HPC I extracted the list of genotypes and partecipants with plink:
# plink --bed /rds/project/wja24/rds-wja24-uk-biobank-gen/plink/ukb_cal_chr1_v2.bed 
# --bim /rds/project/wja24/rds-wja24-uk-biobank-gen/plink/ukb_snp_chr1_v2.bim 
# --fam /rds/project/wja24/rds-wja24-uk-biobank-gen/13745_specific/fam/ukb13745_cal_chr1_v2_s488292.fam 
# --extract ~/shared_luanluan_luca_UKB/VEL/VEL_snps_eqtl_after_pruning_20210317.txt \
# --recode A include-alt --out ./shared_luanluan_luca_UKB/VEL/genotype_VEL_LD_FILTER_ukb13745_20210315
# 
# the henotype id 0  fo hom reference, 1 for het and 2 for hom alt.
# I extracted the pruned list and the one I got from Mattia.

