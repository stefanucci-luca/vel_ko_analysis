# load data
load("ukb_data.RData")
# initialise vectors
measure = c()
mean_WT = c()
median_WT = c()
sd_WT = c()
mean_KO = c()
meadian_KO = c()
sd_KO = c()
# explore and check
skimr::skim(df_final_un_eu)
# loop acroos columns and calculate data
for (column in colnames(df_final_un_eu)) {
  if (is.numeric(df_final_un_eu[,column]) == TRUE){
  measure = c(column, measure)
  mean_WT = c(mean(df_final_un_eu[which(df_final_un_eu$genotype == "0"),column], na.rm = T), mean_WT)
  median_WT = c(median(df_final_un_eu[which(df_final_un_eu$genotype == "0"),column], na.rm = T), median_WT)
  sd_WT = c(sd(df_final_un_eu[which(df_final_un_eu$genotype == "0"),column], na.rm = T), sd_WT)
  mean_KO = c(mean(df_final_un_eu[which(df_final_un_eu$genotype == "2"),column], na.rm = T), mean_KO)
  meadian_KO = c(median(df_final_un_eu[which(df_final_un_eu$genotype == "2"),column], na.rm = T), meadian_KO)
  sd_KO = c(sd(df_final_un_eu[which(df_final_un_eu$genotype == "2"),column], na.rm = T), sd_KO)
  }
}
# export data
write.csv(x = 
data.frame(
  "measure" = measure,
  'mean_WT' = mean_WT,
  'median_WT' = median_WT,
  'sd_WT' = sd_WT,
  'mean_KO' = mean_KO,
  'meadian_KO' = meadian_KO,
  'sd_KO' = sd_KO
), file = "summary_stat_values.csv",
quote = F, row.names = F)
