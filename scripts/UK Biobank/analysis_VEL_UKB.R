# load the libraries used for the analysis
library('lm.beta')
library("questionr")
library('esc')
library('effectsize')
library('tidyverse')
library('see')
library('GGally')
library('RColorBrewer')
library('RNOmni')

# set seed for reproducibility
set.seed(1)

# analysis of UKB data - extracted with the script commented below
df_large = data.table::fread("UKB_VEL_relevant_fields.csv")
df_large = as.data.frame(df_large)
# remove un-assigned genotypes for rs566629828
df_large = filter(df_large, genotype != "./.")

# Rename genotypes as:
# 0/0 = 0
# 0/1 = 1
# 1/1 = 2
df_large$genotype = str_replace_all(string = df_large$genotype,
                                    pattern = "0/0",
                                    replacement = "0")
df_large$genotype = str_replace_all(string = df_large$genotype,
                                    pattern = "0/1",
                                    replacement = "1")
df_large$genotype = str_replace_all(string = df_large$genotype,
                                    pattern = "1/1",
                                    replacement = "2")

# set WT allele (i.e. "0/0" aka "0") as the one to use for reference
# this will be used in the mathematical model as reference. Having hom ref as reference for the effect sizes.
df_large$genotype <-
  relevel(as.factor(df_large$genotype), ref = "0")

# Import ancestry information
# This file contains the list of id of unrelated europen
eur_unre_500k = data.table::fread('british-ancestry-MSUP_unrel_europ_500K.tsv') %>%
  select(ID)

# Import covariates information
load('unb_covariates.rdata')

# times in the hospital
df_final = merge(df, df_large, by = 'eid')

# DF used in the analyses has:
# Unrelated european
# WT/Hom alt genotype
df_final_un_eu = df_final[df_final$eid %in% eur_unre_500k$ID , ]
df_final_un_eu_only_hom = df_final_un_eu[which(df_final_un_eu$genotype == '0' |
                                                 df_final_un_eu$genotype == '2'),]

# df that has to be used in the statistical analysis below
df_to_use = df_final_un_eu_only_hom

# Initialise vector for stats
gender = c()
feature = c()
effect = c()
cohort = c()
sd = c()
CI_low = c()
CI_high = c()
p.val = c()
bmi = c()

# continuous variables
# The for-loop loops across the experimental variables we set. i.e. W/ and W/O BMI and gender as covariates. If gender is not a covariate, the db is gender strata.
for (round_bmi in 1:2) {
  if (round_bmi == 1) {
    for (gen in 0:1) {
      tmp_df = df_to_use[which(df_to_use$sex_f31_0_0 == gen), ]
      if (gen == 0) {
        for (i in c(23, 27:length(colnames(tmp_df)))) {
          # standardise the data
          normalise_tmp = (tmp_df[, i] - mean(as.vector(tmp_df[, i]), na.rm = TRUE)) / sd(tmp_df[, i], na.rm = TRUE)
          normalise_tmp = RankNorm(normalise_tmp)
          model1 = lm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages)
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
            feature = c(feature, colnames(tmp_df)[i])
            gender = c(gender, gen)
            cohort = c(cohort, param$Parameter[num])
            effect = c(effect, param$Coefficient[num])
            sd = c(sd, param$SE[num])
            CI_low = c(CI_low, param$CI_low[num])
            CI_high = c(CI_high, param$CI_high[num])
            p.val = c(p.val, param$p[num])
            bmi = c(bmi, 'NOT_BMI_as_covariates')
          }
        }
      } else {
        for (i in c(23, 27, 30:length(colnames(tmp_df)))) {
          # standardise the data
          normalise_tmp = (tmp_df[, i] - mean(as.vector(tmp_df[, i]), na.rm = TRUE)) / sd(tmp_df[, i], na.rm = TRUE)
          normalise_tmp = RankNorm(normalise_tmp)
          model1 = lm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages)
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
            feature = c(feature, colnames(tmp_df)[i])
            gender = c(gender, gen)
            cohort = c(cohort, param$Parameter[num])
            effect = c(effect, param$Coefficient[num])
            sd = c(sd, param$SE[num])
            CI_low = c(CI_low, param$CI_low[num])
            CI_high = c(CI_high, param$CI_high[num])
            p.val = c(p.val, param$p[num])
            bmi = c(bmi, 'NOT_BMI_as_covariates')
          }
        }
      }
    }
    for (i in c(23, 27:length(colnames(df_to_use)))) {
      # standardise the data
      normalise_tmp = (df_to_use[, i] - mean(as.vector(df_to_use[, i]), na.rm = TRUE)) / sd(df_to_use[, i], na.rm = TRUE)
      normalise_tmp = RankNorm(normalise_tmp)
      model1 = lm(
        normalise_tmp ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex
      )
      param = parameters::model_parameters(model1)
      for (num in 1:length(unique(df_to_use$genotype))) {
        feature = c(feature, colnames(df_to_use)[i])
        gender = c(gender, 'not_gender_stratified but as covariates')
        cohort = c(cohort, param$Parameter[num])
        effect = c(effect, param$Coefficient[num])
        sd = c(sd, param$SE[num])
        CI_low = c(CI_low, param$CI_low[num])
        CI_high = c(CI_high, param$CI_high[num])
        p.val = c(p.val, param$p[num])
        bmi = c(bmi, 'NOT_BMI_as_covariates')
      }
    }
  } else {
    for (gen in 0:1) {
      tmp_df = df_to_use[which(df_to_use$sex_f31_0_0 == gen), ]
      if (gen == 0) {
        for (i in c(23, 27:length(colnames(tmp_df)))) {
          # standardise the data
          normalise_tmp = (tmp_df[, i] - mean(as.vector(tmp_df[, i]), na.rm = TRUE)) / sd(tmp_df[, i], na.rm = TRUE)
          normalise_tmp = RankNorm(normalise_tmp)
          model1 = lm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi)
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
            feature = c(feature, colnames(tmp_df)[i])
            gender = c(gender, gen)
            cohort = c(cohort, param$Parameter[num])
            effect = c(effect, param$Coefficient[num])
            sd = c(sd, param$SE[num])
            CI_low = c(CI_low, param$CI_low[num])
            CI_high = c(CI_high, param$CI_high[num])
            p.val = c(p.val, param$p[num])
            bmi = c(bmi, 'BMI_as_covariates')
          }
        }
      } else {
        for (i in c(23, 27, 30:length(colnames(tmp_df)))) {
          # standardise the data
          normalise_tmp = (tmp_df[, i] - mean(as.vector(tmp_df[, i]), na.rm = TRUE)) / sd(tmp_df[, i], na.rm = TRUE)
          normalise_tmp = RankNorm(normalise_tmp)
          model1 = lm(normalise_tmp ~ factor(tmp_df$genotype) + tmp_df$ages + tmp_df$bmi)
          param = parameters::model_parameters(model1)
          for (num in 1:length(unique(tmp_df$genotype))) {
            feature = c(feature, colnames(tmp_df)[i])
            gender = c(gender, gen)
            cohort = c(cohort, param$Parameter[num])
            effect = c(effect, param$Coefficient[num])
            sd = c(sd, param$SE[num])
            CI_low = c(CI_low, param$CI_low[num])
            CI_high = c(CI_high, param$CI_high[num])
            p.val = c(p.val, param$p[num])
            bmi = c(bmi, 'BMI_as_covariates')
          }
        }
      }
    }
    for (i in c(23, 27:length(colnames(tmp_df)))) {
      # standardise the data
      normalise_tmp = (df_to_use[, i] - mean(as.vector(df_to_use[, i]), na.rm = TRUE)) / sd(df_to_use[, i], na.rm = TRUE)
      normalise_tmp = RankNorm(normalise_tmp)
      model1 = lm(
        normalise_tmp ~ factor(df_to_use$genotype) + df_to_use$ages + df_to_use$sex + df_to_use$bmi
      )
      param = parameters::model_parameters(model1)
      for (num in 1:length(unique(df_to_use$genotype))) {
        feature = c(feature, colnames(df_to_use)[i])
        gender = c(gender, 'not_gender_stratified but as covariates')
        cohort = c(cohort, param$Parameter[num])
        effect = c(effect, param$Coefficient[num])
        sd = c(sd, param$SE[num])
        CI_low = c(CI_low, param$CI_low[num])
        CI_high = c(CI_high, param$CI_high[num])
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

# export table
write.table(
  stat_ukb,
  file = 'ukb_stats.csv',
  row.names = F,
  append = F,
  quote = F,
  sep = ","
)

# All disease phenotypes - Categorical data analyses

# load the disease dataframe (i.e. ukb_events)
load('first_events.rdata')

# add the phenotype to the genotype df
df_to_use_pheno = merge(df_to_use, ukb_events, by.x = "eid", "identifier")

# select the columns which have the disease of interest (i.e. labelled with "cal_p")
pheno_interest = which(startsWith(colnames(df_to_use_pheno), "cal_p"))

# convert entries to 0 (i.e. phenotype not available or absent) and 1 (phenotype present)
for (var in colnames(df_to_use_pheno[, pheno_interest])) {
  df_to_use_pheno[, var] = ifelse(is.na(df_to_use_pheno[, var]), 0, 1)
}

# initialise vector
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

for(i in which(startsWith(colnames(df_to_use_pheno), "cal_p"))) {
  print(i)
  normalise_tmp = df_to_use_pheno[, i] # no normalisation # simply to use the previous structure of the script
  model1 = glm(
    normalise_tmp ~ factor(df_to_use_pheno$genotype) + df_to_use_pheno$ages + df_to_use_pheno$sex_f31_0_0 + df_to_use_pheno$bmi, 
    family = "binomial"
  )
  param = parameters::model_parameters(model1, df_method = "wald")
  for (num in 1:length(unique(df_to_use_pheno$genotype))) {
    feature = c(feature, colnames(df_to_use_pheno)[i])
    gender = c(gender, 'not_gender_stratified but as covariates')
    cohort = c(cohort, param$Parameter[num])
    effect = c(effect, param$Coefficient[num])
    sd = c(sd, param$SE[num])
    CI_low = c(CI_low, param$CI_low[num])
    CI_high = c(CI_high, param$CI_high[num])
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

# export the df
xlsx::write.xlsx(
  stat_ukb_pheno,
  file = 'ukb_stats_pheno_clinic_hosp.xlsx',
  col.names = T,
  row.names = F,
  append = F
)

# source the script to draw the rain/cloud plots
source(
  "https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R"
)

# set the theme for the rain/could plots
raincloud_theme <- theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.position = "right",
  plot.title = element_text(
    lineheight = .8,
    face = "bold",
    size = 16
  ),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(
    colour = "black",
    size = 0.5,
    linetype = "solid"
  ),
  axis.line.y = element_line(
    colour = "black",
    size = 0.5,
    linetype = "solid"
  )
)

plot_value = c(
  "triglycerides",
  "alanine_aminotransferase",
  "aspartate_aminotransferase",
  "gamma_glutamyltransferase",
  "leg_fat_mass_right",
  "leg_fatfree_mass_right",
  "leg_fat_mass_left",
  "leg_fatfree_mass_left",
  "arm_fat_mass_right",
  "arm_fatfree_mass_right",
  "arm_predicted_mass_right",
  "arm_fat_mass_left",
  "arm_fatfree_mass_left",
  "arm_predicted_mass",
  "urate"
)

# the for loop goes across the plot_value vector and plot each entry 
for (feata in colnames(df_to_use)[which(colnames(df_to_use) %in% plot_value)]) {
  p1 = ggplot(data = df_to_use[df_to_use$genotype %in% c(0, 2), ], aes(x = genotype, y = df_to_use[df_to_use$genotype %in% c(0, 2), feata], fill = genotype)) +
    ggtitle(label = feata) +
    # geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, bw = 5) +
    geom_jitter(
      aes(color = genotype),
      position = position_jitter(width = .3),
      size = .28,
      alpha = 0.3
    ) +
    geom_violin(
      width = .40,
      outlier.shape = NA,
      alpha = 0.5,
      position = position_nudge(x = 0.5),
      draw_quantiles = c(0.25, 0.5, 0.75),
      trim = F
    ) +
    #expand_limits(x = 3.25) +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    scale_color_manual(values = c("#3C5488FF", "#990000")) +
    scale_fill_manual(values = c("#3C5488FF", "#990000")) +
    theme_minimal() +
    #coord_cartesian(
    #ylim = c(0, quantile(df_to_use[df_to_use$genotype %in% c(0,2),feata], na.rm = T)[[4]])
    #)+
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA)
    )# +
  # facet_wrap('sex')
  ggsave(
    filename = paste0("rainplot_paper_", feata, ".png"),
    path = 'result_ukb_vel/',
    dpi = "retina",
    width = 10,
    height = 20,
    units = "cm",
    device = 'png',
    plot = p1
  )
}


# Check the differences in the distribution across the UKB and NIHR cohorts 

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
