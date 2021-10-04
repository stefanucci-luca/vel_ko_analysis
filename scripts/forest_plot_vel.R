library(ggplot2)

ukb_bmi = readxl::read_xlsx("/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/Tables/ukb_stats_table_paper_formatted.xlsx", sheet = 1)
ukb_no_bmi = readxl::read_xlsx("/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/Tables/ukb_stats_table_paper_formatted.xlsx", sheet = 4)

ukb_no_bmi$feature <- factor(ukb_no_bmi$feature,levels = c("bmi", "weight", "waist circumference",
                                                           "triglycerides", "ldl", "arm fat mass left", "arm fat mass right",
                                                           "alanine aminotransferase", "aspartate aminotransferase", "gamma-glutamyl transferase", "urate"
                                                     ))

p1 = ggplot(data=ukb_no_bmi,
           aes(x = "VEL",y = effect, ymin = CI_low, ymax = CI_high))+
  geom_pointrange(color = ggsci::pal_startrek("uniform")(7)[1])+
  geom_hline(yintercept =0, linetype=2)+
  xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high),width=0.05,cex=0.21, color = ggsci::pal_startrek("uniform")(7)[1])+ 
  scale_color_continuous(ggsci::pal_startrek("uniform")(7)[1]) +
  facet_wrap(~feature, nrow=15, scales = "free_y", strip.position = "left" )+
  theme_bw(base_line_size = 0, base_rect_size = 0) +
  coord_flip() + 
  scale_y_continuous(breaks = c(-1.3,-0.9,-.6,-.3,0,.5,1,1.5)) 

p1 +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), 
        panel.background = element_rect("transparent"), 
        panel.spacing=unit(-5, units = "cm"),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(angle = 0) )


ukb_bmi$feature <- factor(ukb_bmi$feature,levels = c("alanine aminotransferase", 
                                                     "aspartate aminotransferase",
                                                     "gamma-glutamyl transferase",
                                                     "urate",
                                                     "ldl"))

p = ggplot(data=ukb_bmi,
           aes(x = "VEL",y = effect, ymin = CI_low, ymax = CI_high))+
  geom_pointrange(color = ggsci::pal_startrek("uniform")(7)[1])+
  geom_hline(yintercept =0, linetype=2)+
  xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high),width=0.05,cex=0.21, color = ggsci::pal_startrek("uniform")(7)[1])+ 
  scale_color_continuous(ggsci::pal_startrek("uniform")(7)[1]) +
  facet_wrap(~feature, nrow=9, scales = "free_y", strip.position = "left" )+
  theme_bw(base_line_size = 0, base_rect_size = 0) +
  coord_flip() +
  scale_y_continuous(breaks = c(-1.3,-0.9,-.6,-.3,0,.5,1,1.5))

p+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), 
        panel.background = element_rect("transparent"), 
        panel.spacing=unit(-5, units = "cm"),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(angle = 0) )



wgs10_bmi = readxl::read_xls("/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/Tables/WGS10_stats_table_paper_formatted.xls", sheet = 1)
wgs10_no_bmi = readxl::read_xls("/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/Tables/WGS10_stats_table_paper_formatted.xls", sheet = 4)

p3 = ggplot(data=wgs10_no_bmi,
            aes(x = "VEL",y = effect, ymin = CI_low, ymax = CI_high))+
  geom_pointrange(color = ggsci::pal_startrek("uniform")(7)[1])+
  geom_hline(yintercept =0, linetype=2)+
  xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high),width=0.05,cex=0.21, color = ggsci::pal_startrek("uniform")(7)[1])+ 
  scale_color_continuous(ggsci::pal_startrek("uniform")(7)[1]) +
  facet_wrap(~feature, nrow=12, scales = "free_y", strip.position = "top" )+
  theme_bw(base_line_size = 0, base_rect_size = 0) +
  coord_flip() +
  scale_y_continuous(breaks = c(-1.3,-0.9,-.6,-.3,0,.5,1,1.5))

p3+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), 
        panel.background = element_rect("transparent"), 
        panel.spacing=unit(-5, units = "cm"),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(angle = 0) )

p2 = ggplot(data=wgs10_bmi,
           aes(x = "VEL",y = effect, ymin = CI_low, ymax = CI_high))+
  geom_pointrange(color = ggsci::pal_startrek("uniform")(7)[1])+
  geom_hline(yintercept =0, linetype=2)+
  xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high),width=0.05,cex=0.21, color = ggsci::pal_startrek("uniform")(7)[1])+ 
  scale_color_continuous(ggsci::pal_startrek("uniform")(7)[1]) +
  facet_wrap(~feature, nrow=12, strip.position = "left" )+
  theme_bw(base_line_size = 0, base_rect_size = 0) +
  coord_flip() +
  scale_y_continuous(breaks = c(-1.3,-0.9,-.6,-.3,0,.5,1,1.5))

  
p2+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), 
        panel.background = element_rect("transparent"), 
        panel.spacing=unit(-5, units = "cm"),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(angle = 0) )


ukb_pheno = readxl::read_xlsx("/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/Tables/ukb_stats_pheno_clinic_paper_format.xlsx", sheet = 1)
ukb_pheno_bmi_cov = readxl::read_xlsx("/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/Tables/ukb_stats_pheno_clinic_paper_format.xlsx", sheet = 2)

p4 = ggplot(data=ukb_pheno,
            aes(x = "VEL",y = `effect log(OR)`, ymin = CI_low, ymax = CI_high))+
  geom_pointrange(color = ggsci::pal_startrek("uniform")(7)[1])+
  geom_hline(yintercept =0, linetype=2)+
  xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high),width=0.05,cex=0.21, color = ggsci::pal_startrek("uniform")(7)[1])+ 
  scale_color_continuous(ggsci::pal_startrek("uniform")(7)[1]) +
  facet_wrap(~Disease, nrow=12, strip.position = "left" )+
  theme_bw(base_line_size = 0, base_rect_size = 0) +
  coord_flip() +
  scale_y_continuous(breaks = c(-1.3,-0.9,-.6,-.3,0,.5,1,1.5))


p4+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), 
        panel.background = element_rect("transparent"), 
        panel.spacing=unit(-5, units = "cm"),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(angle = 0) )


p5 = ggplot(data=ukb_pheno_bmi_cov,
            aes(x = "VEL",y = `effect log(OR)`, ymin = CI_low, ymax = CI_high))+
  geom_pointrange(color = ggsci::pal_startrek("uniform")(7)[1])+
  geom_hline(yintercept =0, linetype=2)+
  xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high),width=0.05,cex=0.21, color = ggsci::pal_startrek("uniform")(7)[1])+ 
  scale_color_continuous(ggsci::pal_startrek("uniform")(7)[1]) +
  facet_wrap(~Disease, nrow=12, strip.position = "left" )+
  theme_bw(base_line_size = 0, base_rect_size = 0) +
  coord_flip() +
  scale_y_continuous(breaks = c(-1.3,-0.9,-.6,-.3,0,.5,1,1.5))


p5+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), 
        panel.background = element_rect("transparent"), 
        panel.spacing=unit(-5, units = "cm"),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(angle = 0) )


gridExtra::grid.arrange(p, p1,p2, p3, p4, p5, ncol = 1)

