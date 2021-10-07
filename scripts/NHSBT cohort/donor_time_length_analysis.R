library(tidyverse)
library(ggpubr)
library(rstatix)
library(RColorBrewer)

df_VEL = xlsx::read.xlsx('~/Desktop/VEL/Rare_red_cell_pheno_donors V2.xlsx', sheetIndex = 1, header = T, colClasses = c("character", "character", "numeric", "Date", "character", "character", "Date", "character"))
df_VEL$cohort = "VEL"
df_K_neg = xlsx::read.xlsx('~/Desktop/VEL/Rare_red_cell_pheno_donors V2.xlsx', sheetIndex = 2, header = T, colClasses = c("character", "character", "numeric", "Date", "character", "character", "Date", "character"))
df_K_neg$cohort = "K_negative"
df_Kpb_neg = xlsx::read.xlsx('~/Desktop/VEL/Rare_red_cell_pheno_donors V2.xlsx', sheetIndex = 3, header = T, colClasses = c("character", "character", "numeric", "Date", "character", "character", "Date", "character"))
df_Kpb_neg$cohort = "Kpb_negative"

df_final = Reduce(rbind, list(df_VEL, df_Kpb_neg, df_K_neg))

#I am assuming that if referra date is missing, the person is still in the active donor list
df_final$Deferral_date = replace_na(df_final$Deferral_date, Sys.time()) 

df_final$diffr = ifelse(!is.na(df_final$Deferral_date),
              difftime(df_final$Deferral_date,df_final$Registration_date,units="days"),
              difftime(Sys.time(),df_final$Registration_date,units="days"))

df_final = df_final[which(df_final$diffr > 0 ),]

# check time enrolled as donor 
######################################################################################
#                       _                      _     _   _                
#                      | |                    | |   | | (_)               
#   ___ _ __  _ __ ___ | |_ __ ___   ___ _ __ | |_  | |_ _ _ __ ___   ___ 
#  / _ \ '_ \| '__/ _ \| | '_ ` _ \ / _ \ '_ \| __| | __| | '_ ` _ \ / _ \
# |  __/ | | | | | (_) | | | | | | |  __/ | | | |_  | |_| | | | | | |  __/
#  \___|_| |_|_|  \___/|_|_| |_| |_|\___|_| |_|\__|  \__|_|_| |_| |_|\___|
#######################################################################################

# Length (in days) that a participants have been in the NHSBT cohort
# quick preview to the data
ggboxplot(df_final, x = "cohort", y = "diffr")
# QQplot
ggqqplot(df_final, "diffr", facet.by = "cohort")
# As expected QQplot show that data are not normally distibuted, but are right skewed.
# Shapito test for normality
df_final %>%
  group_by(cohort) %>%
  shapiro_test(diffr)

# shapiro is showing that data are not normally distributed. I'll go for a non-parametric ANOVA
# i.e Kruskal-wallis anova
# check p-val
res.aov.kw <- df_final %>% kruskal_test(diffr ~ cohort)
res.aov.kw
# effect size
res.aov.kw.eff <- df_final %>% kruskal_effsize(diffr ~ cohort)
res.aov.kw.eff
# pairwise comparison with post hoc bonferroni correction
pwc <- df_final %>% dunn_test(diffr ~ cohort, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "cohort")
pwc
# plot the data
ggboxplot(df_final, x = "cohort", y = "diffr") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov.kw, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )
#
#
#
#
#################################################################
#              _ _             _                   _           _ 
#             (_) |           | |                 | |         | |
#  _   _ _ __  _| |_ ___    __| | ___  _ __   __ _| |_ ___  __| |
# | | | | '_ \| | __/ __|  / _` |/ _ \| '_ \ / _` | __/ _ \/ _` |
# | |_| | | | | | |_\__ \ | (_| | (_) | | | | (_| | ||  __/ (_| |
#  \__,_|_| |_|_|\__|___/  \__,_|\___/|_| |_|\__,_|\__\___|\__,_|
#################################################################
# number of donations as proxy for units donated
# quick look at the data
ggboxplot(df_final, x = "cohort", y = "Number_of_donations")
# QQplot
ggqqplot(df_final, "Number_of_donations", facet.by = "cohort")
# Again not normally distributed
# try to correct with og transformation
df_final$log_Num_don = log10(df_final$Number_of_donations)
ggqqplot(df_final, "log_Num_don", facet.by = "cohort")
# log doesn't improve a lot. I'll go for a non-parametric again
# test with shapio if data are normal. i.e Kruskal-wallis anova
df_final %>%
  group_by(cohort) %>%
  shapiro_test(Number_of_donations)
# Kruskal-wallis anova
res.aov.kw1 <- df_final %>% kruskal_test(Number_of_donations ~ cohort)
res.aov.kw1
# Kruskal-wallis effect size
res.aov.kw.eff1 <- df_final %>% kruskal_effsize(Number_of_donations ~ cohort)
res.aov.kw.eff1
# pairwise comparison with post hoc bonferroni correction
pwc1 <- df_final %>% dunn_test(Number_of_donations ~ cohort, p.adjust.method = "bonferroni") 
pwc1 <- pwc1 %>% add_xy_position(x = "cohort")
pwc1
# plot the results
ggboxplot(df_final, x = "cohort", y = "Number_of_donations") +
  stat_pvalue_manual(pwc1, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov.kw1, detailed = TRUE),
    caption = get_pwc_label(pwc1)
  )
#
#
#
#
#
###############################################################################
#                                   _ _              _               _ _       
#                                  | (_)            | |             (_) |      
#  _ __   ___  _ __ _ __ ___   __ _| |_ ___  ___  __| |  _   _ _ __  _| |_ ___ 
# | '_ \ / _ \| '__| '_ ` _ \ / _` | | / __|/ _ \/ _` | | | | | '_ \| | __/ __|
# | | | | (_) | |  | | | | | | (_| | | \__ \  __/ (_| | | |_| | | | | | |_\__ \
# |_| |_|\___/|_|  |_| |_| |_|\__,_|_|_|___/\___|\__,_|  \__,_|_| |_|_|\__|___/
###############################################################################
# Number of units donated by the time spent as a donor (in days)
# normalise units by the number of days
df_final$donation_norm = df_final$Number_of_donations / ifelse(df_final$SEX == "M", (df_final$diffr/4)/365, (df_final$diffr/3)/365)

# Quick look at the data
ggboxplot(df_final, x = "cohort", y = "donation_norm")
# QQplot
ggqqplot(df_final, "donation_norm", facet.by = "cohort")
# they looks better
df_final %>%
  group_by(cohort) %>%
  shapiro_test(donation_norm)
# these data are not normal again
# Kruskal-wallis anova
res.aov.kw2 <- df_final %>% kruskal_test(donation_norm ~ cohort)
res.aov.kw2
# Effect sizes
res.aov.kw.eff2 <- df_final %>% kruskal_effsize(donation_norm ~ cohort)
res.aov.kw.eff2
# pairwise comparison with post hoc bonferroni correction
pwc2 <- df_final %>% dunn_test(donation_norm ~ cohort, p.adjust.method = "bonferroni") 
pwc2 <- pwc2 %>% add_xy_position(x = "cohort")
# Plot the results 
ggboxplot(df_final, x = "cohort", y = "donation_norm") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov.kw2, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )
# zoom in graph
ggboxplot(df_final, x = "cohort", y = "donation_norm") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  ylim(c(0,50))+
  labs(
    subtitle = get_test_label(res.aov.kw2, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )

############################################################################################
#           _ _   _         _                                                               
#          (_) | | |       | |                                                              
# __      ___| |_| |__   __| |_ __ __ ___      ___ __    _ __ ___  __ _ ___  ___  _ __  ___ 
# \ \ /\ / / | __| '_ \ / _` | '__/ _` \ \ /\ / / '_ \  | '__/ _ \/ _` / __|/ _ \| '_ \/ __|
#  \ V  V /| | |_| | | | (_| | | | (_| |\ V  V /| | | | | | |  __/ (_| \__ \ (_) | | | \__ \
#   \_/\_/ |_|\__|_| |_|\__,_|_|  \__,_| \_/\_/ |_| |_| |_|  \___|\__,_|___/\___/|_| |_|___/
############################################################################################  

df_final = df_final %>% 
  filter(Activity_code != "T") %>% 
  filter(Activity_code != "D") %>% 
  filter(Activity_code != "B")
# preview data
ggplot(df_final %>% 
         group_by(cohort,Activity_code) %>% 
         tally(), 
       aes(x = cohort, y = n, fill = Activity_code)) +
  geom_bar(stat = 'identity', position = position_dodge2()) +
  scale_fill_manual(values = ggsci::pal_startrek("uniform")(5))+
  theme_bw()
# create contingency matrix
don_tab <- as.table(rbind(c(1653, 440, 100, 97, 756), c(74, 17, 5, 8, 30)))
dimnames(don_tab) <- list(genotype = c("K_negative","VEL"),
                          activity =c("A","M","R","S","W")  
                          )
# Chi Squared test
chi_donrs = chisq.test(don_tab)
chi_donrs
# Expected values
message("Expected")
round(chi_donrs$expected)
# Observed values in out cohort
message("Observed")
chi_donrs$observed
# plot the values for each activity 
corrplot::corrplot(don_tab, is.cor = FALSE)
# normalise the value of each activity by the number of people enrolled
# i.e. proportion of activity
don_tab_norm = matrix(rbind(don_tab[1,] / sum(don_tab[1,]), 
                            don_tab[2,] / sum(don_tab[2,])
                            ), 
                            ncol = 5)
dimnames(don_tab_norm) <- list(genotype = c("K_negative","VEL"),
                          activity =c("A","M","R","S","W")  
)
corrplot::corrplot(don_tab_norm, is.cor = FALSE)
#
#
#
#
#
#
##########################
#           _____ ______ 
#     /\   / ____|  ____|
#    /  \ | |  __| |__   
#   / /\ \| | |_ |  __|  
#  / ____ \ |__| | |____ 
# /_/    \_\_____|______|
##########################

df_final$AGEGROUP <- factor(df_final$AGEGROUP, levels = c(
  '18 to 22',
  '23 to 27',
  '28 to 32',
  '33 to 37',
  '38 to 42',
  '43 to 47',
  '48 to 52',
  '53 to 57', 
  '58 to 62', 
  '63 to 67', 
  '68 to 72', 
  '73 to 77',
  '78 to 82',
  '83 to 87',
  '88 to 92', 
  '93 to 97'))

# preview data
ggplot(df_final %>% 
         group_by(cohort,AGEGROUP) %>% 
         tally(), 
       aes(x = cohort, y = n, fill = AGEGROUP)) +
  geom_bar(stat = 'identity', position = position_dodge2()) +
  scale_fill_manual(values = rev(colorRampPalette(brewer.pal(9, "Blues"))(16)))+
  theme_bw() +
  theme(legend.position = "bottom")

# create contingency matrix
don_tab <- as.table(rbind(c(47, 177, 225, 297, 
                            292, 287,337,349,
                            315,253,229,160,
                            54,25,3,2), 
                          c(1,0,14,20,
                            14,10,8,21,
                            19,11,8,5,
                            1,2,0,0)))
dimnames(don_tab) <- list(genotype = c("K_negative","VEL"),
                          activity =c("18 to 22","23 to 27","28 to 32","33 to 37",
                                      "38 to 42","43 to 47","48 to 52","53 to 57",
                                      "58 to 62","63 to 67","68 to 72","73 to 77",
                                      "78 to 82","83 to 87","88 to 92","93 to 97")  
)
# Chi Squared test
chi_donrs = chisq.test(don_tab)
chi_donrs
# Expected values
message("Expected")
round(chi_donrs$expected)
# Observed values in out cohort
message("Observed")
chi_donrs$observed
# plot the values for each activity 
corrplot::corrplot(don_tab, is.cor = FALSE)
# normalise the value of each activity by the number of people enrolled
# i.e. proportion of activity
don_tab_norm = matrix(rbind(don_tab[1,] / sum(don_tab[1,]), 
                            don_tab[2,] / sum(don_tab[2,])
), 
ncol = 16)
dimnames(don_tab_norm) <- list(genotype = c("K_negative","VEL"),
                          activity =c("18 to 22","23 to 27","28 to 32","33 to 37",
                                      "38 to 42","43 to 47","48 to 52","53 to 57",
                                      "58 to 62","63 to 67","68 to 72","73 to 77",
                                      "78 to 82","83 to 87","88 to 92","93 to 97") 
)
corrplot::corrplot(don_tab_norm, is.cor = FALSE)
