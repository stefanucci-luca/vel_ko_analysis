library(skimr)
library(ggplot2)
library(gghighlight)
library(tidyverse)
library(reshape2)
library(plyr)
library(ggExtra)
library(ggpubr)
library(rstatix)

df=xlsx::read.xlsx("/Users/luca/Desktop/VEL/laura_data/VEL NEG additional data_working_copy.xlsx", sheetIndex = 1, colIndex = c(1:5))
# Statistical test
test = t_test(REE.Z.Score ~ cohort, data = df)
test
# Effect
tesst_eff = esc::esc_t(test$statistic, test$p, totaln = length(df$cohort), es.type = "d")
tesst_eff
# The standardised mean difference og 0.5 is reporting a medium effect (0.8 is considered large, see https://www-ncbi-nlm-nih-gov.ezp.lib.cam.ac.uk/pmc/articles/PMC2730804/#b5-ptj33_12p700)
dat = data.frame(x=runif(1), y=runif(1))
p = ggplot(df, aes(REE.Z.Score, LM.Z.score, color = cohort, alpha=c(1.8))) + 
  geom_point(size = as.numeric(df$cohort)^2) +
  scale_color_manual(values = c("#4DBBD5FF" ,"#E64B35FF")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    title = get_test_label(test, detailed = TRUE),
    subtitle = paste( "standardized mean difference =", round(tesst_eff$es,3)),
    caption = paste("N.B.: statistal test is referring only to REE score.\n Analysis performed on:" , Sys.Date())
  )# +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1))+
  geom_point(aes(x=x, y=y), data=dat, size=50, shape=1, color="gold4")
# Add distribution on the side
ggMarginal(p, type = 'densigram', groupFill = T)

# differences lean mass are not significant
test_lean = t.test(LM.Z.score ~ cohort, data = df)
test_lean
