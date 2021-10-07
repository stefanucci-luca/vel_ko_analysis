library(tidyverse)
col_classes = c("character", "character","character","character",rep("numeric", 4), rep("factor", 34) )
df_mouse_m = read.csv("Downloads/mouse Tissue and Body Weight_SMIM1_190327.xlsx - Male - body and tissue weight.csv", header = T, colClasses = col_classes)
df_mouse_f = read.csv("Downloads/mouse Tissue and Body Weight_SMIM1_190327.xlsx - Female -body and tissues weight.csv", header = T, colClasses = col_classes)

df_mouse = dplyr::bind_rows(df_mouse_m, df_mouse_f)
# remove the columns with the old stats
colnames(df_mouse)
df_row_data = 
  df_mouse %>% select(1:13) %>% 
  filter(Genotype != 'Het')

#check remover all the Het
unique(df_row_data$Genotype)


df_row_data$Genotype = relevel(as.factor(df_row_data$Genotype), ref = "WT")
df_row_data$normalised_weight = ( df_row_data$`Body..g.` - mean(df_row_data$`Body..g.`) ) / sd(df_row_data$`Body..g.`)

plot1 = ggplot(df_row_data, aes(x = `Age..weeks.`, y = `Body..g.`, color = Genotype)) +
  geom_point() +
  scale_color_manual(values = c("#3C5488FF","#990000")) +
  scale_fill_manual(values = c("#3C5488FF","#990000")) +
  geom_smooth(method = "lm") +
  ggtitle("Not gender stratified",
        paste(
          "effect", format(parameters::parameters(lm(`normalised_weight` ~ Genotype, 
                                                      data = df_row_data))$Coefficient[2], digits = 3 ),
          "p-val", format(parameters::parameters(lm(`normalised_weight` ~ Genotype, data = df_row_data))$p[2],  digits = 3)
              )
  )

plot2 = ggplot(df_row_data, aes(x = `Age..weeks.`, y = `Body..g.`, color = Genotype)) +
  geom_point() +
  scale_color_manual(values = c("#3C5488FF","#990000")) +
  scale_fill_manual(values = c("#3C5488FF","#990000")) +
  geom_smooth(method = "lm") +
  facet_wrap("Gender", scales = "free_x")+
  ggtitle("Gender stratified",
          paste(
            "effect_Male", format(parameters::parameters(lm(`normalised_weight` ~ Genotype, 
                                                 data = df_row_data[which(df_row_data$Gender == "M"),]))$Coefficient[2],  digits = 3),
            "p-val_Male", format(parameters::parameters(lm(`normalised_weight` ~ Genotype, 
                                                data = df_row_data[which(df_row_data$Gender == "M"),]))$p[2],  digits = 3),
            "\neffect_Female", format(parameters::parameters(lm(`normalised_weight` ~ Genotype, 
                                                      data = df_row_data[which(df_row_data$Gender == "F"),]))$Coefficient[2],  digits = 3),
            "p-val_Female", format(parameters::parameters(lm(`normalised_weight` ~ Genotype, 
                                                     data = df_row_data[which(df_row_data$Gender == "F"),]))$p[2],  digits = 3)
          )
  )


library(magrittr)
library(multipanelfigure)
figure1 <- multi_panel_figure(columns = 2, rows = 2, panel_label_type = "lower-alpha")
figure1 %<>%
  fill_panel(plot1, column = 1:2, row = 1) %<>%
  fill_panel(plot2, column = 1:2, row = 2)
ggsave(filename = "ls760@platgenbio.net - Google Drive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/Figures/mouse_weight.png", plot = figure1, dpi = "retina", device = "png")
         
ggsave(filename = "plot_mouse_weight.png", 
       path = "ls760@platgenbio.net - Google Drive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/mouse/mouse_body_weight_from_Rita copy/",
       plot = gridExtra::grid.arrange(plot1, plot2, nrow = 1), width = 16,
       device = "png")



cols_to_loop = c("Body..g.","Spleen.Height..cm.", 
                 "Spleen.width..cm.", "Spleen..g.", 
                 "Heart..g.", "Liver..g.", "Kidneys..g.")

for (var_col in cols_to_loop) {

  df_row_data$normalised_var = ( as.numeric(as.character(df_row_data[,var_col])) -
                                  mean(as.numeric(as.character(df_row_data[,var_col])), 
                                       na.rm = T) ) /
                                  sd(as.numeric(as.character(df_row_data[,var_col])), 
                                     na.rm = T)
  
  plot_tmp1 = ggplot(df_row_data, aes(x = Genotype, y = as.numeric(as.character(df_row_data[,var_col])), color = Genotype)) +
    geom_boxplot() +
    ylab(var_col) +
    scale_color_manual(values = c("#3C5488FF","#990000")) +
    scale_fill_manual(values = c("#3C5488FF","#990000")) +
    ggtitle("Not gender stratified",
            paste(
              "effect", format(parameters::parameters(lm(`normalised_var` ~ Genotype + `Age..weeks.`, 
                                                         data = df_row_data))$Coefficient[2], digits = 3 ),
              "p-val", format(parameters::parameters(lm(`normalised_var` ~ Genotype + `Age..weeks.`, data = df_row_data))$p[2],  digits = 3)
            )
    )
  
  plot_tmp2 = ggplot(df_row_data, aes(x = Genotype, y = as.numeric(as.character(df_row_data[,var_col])), color = Genotype)) +
    geom_boxplot() +
    ylab(var_col) +
    facet_wrap("Gender", scales = "free_x")+
    scale_color_manual(values = c("#3C5488FF","#990000")) +
    scale_fill_manual(values = c("#3C5488FF","#990000")) +
    ggtitle("Gender stratified",
            paste(
              "effect_Male", format(parameters::parameters(lm(`normalised_var` ~ Genotype + `Age..weeks.`, 
                                                              data = df_row_data[which(df_row_data$Gender == "M"),]))$Coefficient[2],  digits = 3),
              "p-val_Male", format(parameters::parameters(lm(`normalised_var` ~ Genotype + `Age..weeks.`, 
                                                             data = df_row_data[which(df_row_data$Gender == "M"),]))$p[2],  digits = 3),
              "\neffect_Female", format(parameters::parameters(lm(`normalised_var` ~ Genotype + `Age..weeks.`, 
                                                                  data = df_row_data[which(df_row_data$Gender == "F"),]))$Coefficient[2],  digits = 3),
              "p-val_Female", format(parameters::parameters(lm(`normalised_var` ~ Genotype + `Age..weeks.`, 
                                                               data = df_row_data[which(df_row_data$Gender == "F"),]))$p[2],  digits = 3)
            )
    )
  
  name1=paste0("plot_not_gender_stratified", var_col, "png")
  ggsave(filename = name1,
         path = "ls760@platgenbio.net - Google Drive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/mouse/mouse_organs_weight_from_Rita", 
         plot = plot_tmp1, 
         device = "png")
  name2=paste0("plot_gender_stratified", var_col, "png")
  ggsave(filename = name2, 
         path = "ls760@platgenbio.net - Google Drive/My Drive/Phd/Shared Luca Mattia/VEL/Manuscript/mouse/mouse_organs_weight_from_Rita",
         plot = plot_tmp2, 
         device = "png")
  }

