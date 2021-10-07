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

## explore the df
# str(wdf)
# summary(wdf)
# ggpairs(wdf, 
#        columns = colnames(wdf[,-2]), )

# Statistical analyses
# Initialise the empty vectors
gender = c()
feature = c()
effect = c()
cohort = c()
CI_low = c()
CI_high = c()
p.val = c()
bmi = c()

# For loop creates 2 steps:
# The first one use BMI as covariate, the second does not.
for (round_bmi in 1:2) {
    if (round_bmi == 1 ) {
        for (gen in c('Male','Female')) {
          tmp_df = wdf[wdf$GENDER == gen,]
          for(i in 9:length(colnames(tmp_df))){
            tmp_var = ( tmp_df[,i] - mean(tmp_df[,i], na.rm = T) ) / sd(tmp_df[,i], na.rm = T)
            model1 = lm( data = tmp_df, tmp_var ~ COHORT + AGE + BMI) 
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
          tmp_var = ( wdf[,i] - mean(wdf[,i], na.rm = T) ) / sd(wdf[,i], na.rm = T)
          model1 = lm( data = wdf, tmp_var ~ COHORT + AGE + BMI + GENDER) 
          param = parameters::parameters(model1)
          feature = c(feature, colnames(wdf)[i])
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
        tmp_df = wdf[wdf$GENDER == gen,]
        for(i in 9:length(colnames(tmp_df))){
          tmp_var = ( tmp_df[,i] - mean(tmp_df[,i], na.rm = T) ) / sd(tmp_df[,i], na.rm = T)
          model1 = lm( data = tmp_df, tmp_var ~ COHORT + AGE)
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
        tmp_var = ( wdf[,i] - mean(wdf[,i], na.rm = T) ) / sd(wdf[,i], na.rm = T)
        model1 = lm( data = wdf, tmp_var ~ COHORT + AGE + GENDER) 
        param = parameters::parameters(model1)
        feature = c(feature, colnames(wdf)[i])
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

# assemble the vector in a data frame
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

# export the data frame in a txt
xlsx::write.xlsx(x = stat_df,
                 file = 'summary_stats_lm_wgs10_vel.xls')

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

# import the aesthetics for the violin/cloud plots
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
# set the theme for the plots
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

# used in the for loop
plot_value = c( "LEPT", 'LAR', 'AT_FFA', 'Total_T3', 'Total_T4','ALT', 'AST', 'hsCRP', 'Ferritin' )
# loop goes across all the entries in the plot_value and plot them
for (feata in colnames(wdf)[which(colnames(wdf) %in% plot_value)]) {
  p1 = ggplot(data = wdf, aes( x = COHORT, y = wdf[,feata], fill = COHORT)) +
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
  ggsave(filename = paste0("rainplot_paper_", feata, "_paper.png" ),
         path = 'result_wp10_vel/' , dpi = "retina", width = 10, height = 20, units = "cm",
         device = 'png', plot = p1)
} 

#######################################################################################################################################
#  _____   _       _          _  _             _    _                     _  _   __   __                                        
# |  __ \ (_)     | |        (_)| |           | |  (_)                   | |(_) / _| / _|                                       
# | |  | | _  ___ | |_  _ __  _ | |__   _   _ | |_  _   ___   _ __     __| | _ | |_ | |_  ___  _ __  ___  _ __    ___  ___  ___ 
# | |  | || |/ __|| __|| '__|| || '_ \ | | | || __|| | / _ \ | '_ \   / _` || ||  _||  _|/ _ \| '__|/ _ \| '_ \  / __|/ _ \/ __|
# | |__| || |\__ \| |_ | |   | || |_) || |_| || |_ | || (_) || | | | | (_| || || |  | | |  __/| |  |  __/| | | || (__|  __/\__ \
# |_____/ |_||___/ \__||_|   |_||_.__/  \__,_| \__||_| \___/ |_| |_|  \__,_||_||_|  |_|  \___||_|   \___||_| |_| \___|\___||___/
#######################################################################################################################################


parameters::parameters(t.test(wdf$AGE[which(wdf$COHORT == "VEL")],
                              wdf$AGE[which(wdf$COHORT == "DonorWP10")]))

parameters::parameters(t.test(wdf$AGE[which(wdf$COHORT == "VEL" & wdf$GENDER == "Female")],
                              wdf$AGE[which(wdf$COHORT == "DonorWP10" & wdf$GENDER == "Female")]))

parameters::parameters(t.test(wdf$AGE[which(wdf$COHORT == "VEL" & wdf$GENDER == "Male")],
                              wdf$AGE[which(wdf$COHORT == "DonorWP10" & wdf$GENDER == "Male")])
)


