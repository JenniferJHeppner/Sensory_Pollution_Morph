#####################
# 2018-2025 Adult Morphology & Condition - Figures

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(lubridate)
library(ggpubr)
library(ggnewscale)
library(ggpubr)

dat <- read.csv("Data/18-26_Adult_Morphology_Cavity.csv")

#---- Clean Columns ----
dat <- dat %>%
  mutate(
    species = as.factor(species),
    sex = as.factor(sex),
    ID = as.factor(ID),
    color_band = as.character(color_band),
    trt = as.factor(trt),
    section = as.factor(section),
    site = as.factor(site),
    date = as.Date(date, "%m/%d/%y"),
    julian = yday(date),
    year = as.factor(year),
    nest = as.factor(nest),
    age = as.factor(age),
    stage = as.factor(stage),
    stage = factor(stage, levels = c("TE", "NB", "Eggs", "Chx", "Unk")),
    lux = as.numeric(lux),
    light = as.factor(light),
    LAeq = as.numeric(LAeq),
    LZeq = as.numeric(LZeq),
    LAF90 = as.numeric(LAF90),
    mass = as.numeric(mass),
    tarsus = as.numeric(tarsus),
    wing = as.numeric(wing),
    tail = as.numeric(tail),
    eye = as.numeric(eye),
    bill.l = as.numeric(bill.l),
    bill.tn = as.numeric(bill.tn),
    bill.d = as.numeric(bill.d),
    bill.w = as.numeric(bill.w),
    muscle = as.integer(muscle),
    fat = as.integer(fat),
    cp_score = as.numeric(cp_score),
    measurer = as.factor(measurer),
    nb_date = as.Date(nb_date, "%m/%d/%y"),
    nb_julian = yday(nb_date),
    egg_date = as.Date(egg_date, "%m/%d/%y"),
    egg_julian = yday(egg_date)
  )

# Rename groups
levels(dat$trt)
levels(dat$trt) <- c("Control", "Light", "Combined", "Noise")
levels(dat$trt)
levels(dat$sex)
levels(dat$sex) <- c("Female", "Male", "Unknown")
levels(dat$sex)

# Reorder levels 
dat <- dat %>%
  mutate(trt = fct_relevel(trt, "Control", "Light", "Noise", "Combined"))

# Reorder the levels of the age factor
dat <- dat %>%
  mutate(age = factor(age, levels = c("SY", "AHY", "ASY")))

# Define colors
trt_cols <- c("Control" = "steelblue", "Noise" = "firebrick",
              "Light" = "#FFE082", "Combined" = "coral")


#---- Body Condition Function ----
# Body Condition Function with mass and length
smi  = function (M, L, plot = TRUE, ...)
{
  if (plot) plot(log(M)~log(L), ...)
  {
    if (require(smatr)) {
      SMA = sma(log(M)~log(L))
      bSMA = coef(SMA)[2]}
    else {
      OLS = lm(log(M)~log(L))
      bOLS = coef(OLS)[2]
      r = cor.test(~log(M)+log(L), method = "pearson")$estimate
      #outliers = which(abs(rstandard(ols))>3)
      bSMA = bOLS/r }
  }
  L0 = mean(L, na.rm = T)
  SMi = M*((L0/L)^bSMA)
  return(SMi)
}

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Figure 1 - MOBL Morphology -----
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Figure 1A - MOBL Wing Length -----
# Subset the Data to only include MOBL
mobl <- dat %>% filter(species == "MOBL") %>% droplevels()

# Linear Model
mobl_wing <- lm(wing ~ trt + sex, data = mobl)
summary(mobl_wing)
confint(mobl_wing)
confint(mobl_wing, level = 0.85)

# MOBL WING FIGURE
mobl_wing <- ggplot(mobl, aes(x = trt, y = wing, color = trt)) +
  geom_point(aes(shape = sex), position = position_jitterdodge(jitter.width = 0.3, 
            dodge.width = 0.5), size = 4, alpha = 0.8) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               aes(group = sex), position = position_dodge(width = 0.5),
               width = 0.15, linewidth = 1, color = "black", alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", 
               aes(group = sex, shape = sex), position = position_dodge(width = 0.5),
               size = 4, color = "black", stroke = 1.5, show.legend = FALSE) +
  annotate("text", x = 2, y = 117.5, label = "*", size = 10, color = "grey25") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16)) + 
   scale_color_manual(values = trt_cols) +
   labs(x = "Treatment Group",
       y = "Wing Length (mm)",
       #title = "Mountain Bluebird", 
       color = "Treatment", shape = "Sex") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(color = "black") ,
        #legend.position = "none",
        axis.line = element_line(linewidth = 0.8),
        plot.title = element_text(hjust = 0.5))
mobl_wing



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Figure 1B - MOBL Body Mass -----

# LMER
mobl_mass <- lmer(mass ~ scale(lux) + sex + (1 | section) + (1 | year) + (1 | stage) +  (1 | measurer), data = mobl, REML = FALSE)
summary(mobl_mass)
confint(mobl_mass)
confint(mobl_mass, level = 0.85)

# Mass across lux values using model residuals
predict_lux_sex <- predict_response(mobl_mass, 
                                    terms = c("lux [all]", "sex"), 
                                    margin = "marginalmeans", 
                                    ci.lvl = 0.95)


# MOBL Mass vs Lux
# Marginal effect (the predicted mean)
mobl_mass_lux <- ggplot() +
  geom_ribbon(data = predict_lux_sex, 
              aes(x = x, ymin = conf.low, ymax = conf.high, group = group), 
              fill = "#FFE082", alpha = 0.25) +
  geom_line(data = predict_lux_sex, 
            aes(x = x, y = predicted, group = group, linetype = group), 
            color = "#FFE082", linewidth = 1.5) +
  geom_point(data = mobl, 
             aes(x = lux, y = mass, shape = sex), 
             color = "black", size = 4, alpha = 0.8) +
  annotate("text", x = 1.9, y = 36, label = "*", size = 10, color = "grey25") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(
   x = "Light Intensity (Lux)",
    y = "Body Mass (g)"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    #legend.position = "right",
    legend.position = "none")

mobl_mass_lux


#------------------------------------------------------------------------------
#---- Figure 1C - MOBL BC -----
# Calculate Body Condition using mass and wing
BCwing <- smi(mobl$mass, mobl$wing)
mobl$bc_wing <- BCwing

# Linear Model
mobl_bc <- lm(bc_wing ~ trt, data = mobl)
summary(mobl_bc)
confint(mobl_bc)
confint(mobl_bc, level = 0.85)

# MOBL BC FIGURE
mobl_BC <- ggplot(mobl, aes(x = trt, y = bc_wing, color = trt)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.5), 
             size = 4, alpha = 0.8) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               position = position_dodge(width = 0.7),
               width = 0.15, linewidth = 1, color = "black", alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", 
               position = position_dodge(width = 0.5),
               size = 4, color = "black", stroke = 1.5, show.legend = FALSE) +
  scale_color_manual(values = trt_cols) +
  annotate("text", x = 2, y = 38.5, label = "*", size = 10, color = "grey25") +
  annotate("text", x = 3, y = 39, label = "†", size = 5, color = "grey25") +
  labs(x = "Treatment Group",
       y = "Body Condition",
       #title = "Mountain Bluebird",
       color = "Treatment") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(color = "black") ,
        legend.position = "none",
        #legend.position = "left",
        axis.line = element_line(linewidth = 0.8),
        plot.title = element_text(hjust = 0.5))
mobl_BC



#------------------------------------------------------------------------------
#---- Figure 1 Multi-Panel -----
# Combining figures
figure1 <- ggarrange(mobl_wing, mobl_mass_lux, mobl_BC,
                     labels = c("A", "B", "C"),
                     ncol = 3, nrow = 1)
figure1

# Add Multi-panel Figure Title
figure1_final <- annotate_figure(figure1,
                                 top = text_grob("Mountain Bluebird Morphology", 
                                                 color = "black", 
                                                 #face = "bold", 
                                                 size = 20,
                                                 vjust = -0.1))
figure1_final

# Add extra margin at the top for Title
figure1_final <- figure1_final + 
  theme(plot.margin = margin(t = 8, r = 0, b = 0, l = 0)) 
figure1_final

ggsave(
  filename = "Figure_1_MOBL_Morph_26.png", 
  plot = last_plot(),         
  width = 330, 
  height = 100, 
  units = "mm", 
  dpi = 300,                  
  bg = "white"                
)

ggsave(
  filename = "Figure_1_Legend.png", 
  plot = last_plot(),         
  width = 330, 
  height = 100, 
  units = "mm", 
  dpi = 300,                  
  bg = "white"                
)








#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Figure 2 - WEBL Morphology -----
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Figure 2A - WEBL Mass Light -----
# LMER
webl_mass <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + age + (1 | section) + (1 | stage), data = webl, REML = FALSE)
summary(webl_mass)
confint(webl_mass)
confint(webl_mass, level = 0.85)

# Mass across lux values using model residuals
predict_webl_lux_sex <- predict_response(webl_mass, 
                                         terms = c("lux [all]", "sex"), 
                                         margin = "marginalmeans", 
                                         ci.lvl = 0.85)


# WEBL Mass vs Lux
# Marginal effect (the predicted mean)
webl_mass_lux <- ggplot() +
  geom_ribbon(data = predict_webl_lux_sex, 
              aes(x = x, ymin = conf.low, ymax = conf.high, group = group), 
              fill = "#FFE082", alpha = 0.3) +
  geom_line(data = predict_webl_lux_sex, 
            aes(x = x, y = predicted, group = group, linetype = group), 
            color = "#FFE082", linewidth = 1.5) +
  geom_point(data = webl, 
             aes(x = lux, y = mass, shape = sex), 
             color = "black", size = 4, alpha = 0.8) +
  annotate("text", x = 1.5, y = 32, label = "†", size = 5, color = "grey25") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(x = "Light Intensity (Lux)",
       y = "Body Mass (g)"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    #legend.position = "right")
    legend.position = "none")

webl_mass_lux


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Figure 2B - WEBL Mass Noise -----

# Same model as Mass Lux

# Mass across noise values using model residuals
predict_webl_noise_sex <- predict_response(webl_mass, 
                                           terms = c("LAeq [all]", "sex"), 
                                           margin = "marginalmeans", 
                                           ci.lvl = 0.85)


# WEBL Mass vs Noise
# Marginal effect (the predicted mean)
webl_mass_noise <- ggplot() +
  geom_ribbon(data = predict_webl_noise_sex, 
              aes(x = x, ymin = conf.low, ymax = conf.high, group = group), 
              fill = "firebrick", alpha = 0.2) +
  geom_line(data = predict_webl_noise_sex, 
            aes(x = x, y = predicted, group = group, linetype = group), 
            color = "firebrick", linewidth = 1.5) +
  geom_point(data = webl, 
             aes(x = LAeq, y = mass, shape = sex), 
             color = "black", size = 4, alpha = 0.8) +
  annotate("text", x = 53, y = 32, label = "†", size = 5, color = "grey25") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(x = "Noise Level (LAeq)",
       y = "Body Mass (g)"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    #legend.position = "right")
    legend.position = "none")

webl_mass_noise




#------------------------------------------------------------------------------
#---- Figure 2C - WEBL BC Lux -----
# Calculate Body Condition using mass and wing
BCwing <- smi(webl$mass, webl$wing)
webl$bc_wing <- BCwing

# LMER
webl_BC_int <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + (1 | section) + (1 | stage), data = webl, REML = FALSE)
summary(webl_BC_int)
confint(webl_BC_int)
confint(webl_BC_int, level = 0.85)

# Body Condition across lux values using model residuals
predict_webl_lux_sex <- predict_response(webl_BC_int, 
                                         terms = c("lux [all]", "sex"), 
                                         margin = "marginalmeans", 
                                         ci.lvl = 0.95)


# WEBL BC vs Lux
# Marginal effect (the predicted mean)
webl_BC_lux <- ggplot() +
  geom_ribbon(data = predict_webl_lux_sex, 
              aes(x = x, ymin = conf.low, ymax = conf.high, group = group), 
              fill = "#FFE082", alpha = 0.3) +
  geom_line(data = predict_webl_lux_sex, 
            aes(x = x, y = predicted, group = group, linetype = group), 
            color = "#FFE082", linewidth = 1.5) +
  geom_point(data = webl, 
             aes(x = lux, y = bc_wing, shape = sex), 
             color = "black", size = 4, alpha = 0.8) +
  annotate("text", x = 1.5, y = 32.2, label = "*", size = 10, color = "grey25") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(x = "Light Intensity (Lux)",
       y = "Body Condition"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    #legend.position = "right")
    legend.position = "none")

webl_BC_lux




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Figure 2E - WEBL BC Noise -----
# Body Condition across noise values using model residuals
predict_webk_noise_sex <- predict_response(webl_BC_int, 
                                           terms = c("LAeq [all]", "sex"), 
                                           margin = "marginalmeans", 
                                           ci.lvl = 0.85)


# WEBL BC vs Noise
# Marginal effect (the predicted mean)
webl_BC_noise <- ggplot() +
  geom_ribbon(data = predict_webl_noise_sex, 
              aes(x = x, ymin = conf.low, ymax = conf.high, group = group), 
              fill = "firebrick", alpha = 0.2) +
  geom_line(data = predict_webl_noise_sex, 
            aes(x = x, y = predicted, group = group, linetype = group), 
            color = "firebrick", linewidth = 1.5) +
  geom_point(data = webl, 
             aes(x = LAeq, y = bc_wing, shape = sex), 
             color = "black", size = 4, alpha = 0.8) +
  annotate("text", x = 52, y = 32.5, label = "†", size = 5, color = "grey25") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(x = "Noise Level (LAeq)",
       y = "Body Condition"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    #legend.position = "right")
    legend.position = "none")

webl_BC_noise






#------------------------------------------------------------------------------
#---- Figure 2E - WEBL BC Int -----

# Base off of same model in Figure 2B
# LMER
webl_BC_int <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + (1 | section) + (1 | stage), data = webl, REML = FALSE)
summary(webl_BC_int)
confint(webl_BC_int)
confint(webl_BC_int, level = 0.85)


# Using model results, split noise values between lux SD -1, 0, 1
inter_data_webl <- predict_response(webl_BC_int, terms = c("LAeq [all]", "lux [-1, 0, 1]"))

plot(inter_data_webl)

# Convert the prediction to a data frame
plot_df_webl <- as.data.frame(inter_data_webl)

# Ensure your plot has the labels
plot_df_webl$lux_label <- factor(plot_df_webl$group, 
                                 levels = c("-1", "0", "1"),
                                 labels = c("Low Lux (-1 SD)", "Mean Lux", "High Lux (+1 SD)"))

# Find mean lux to coordinate in figure
scales::rescale(webl$lux)
mean(scales::rescale(webl$lux), na.rm = T) # 0.2


# Define your consistent anchor colors
low_lux_col <- "midnightblue"  
mean_lux_col <- "magenta4" 
high_lux_col <- "#FFE082"

# Discrete palette for model lines
int_cols <- c(
  "Low Lux (-1 SD)"  = low_lux_col, 
  "Mean Lux"         = mean_lux_col, 
  "High Lux (+1 SD)" = high_lux_col
)

webl_BC_int <- ggplot() +
  geom_ribbon(data = plot_df_webl, 
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = lux_label), 
              alpha = 0.3, show.legend = FALSE) + 
  geom_line(data = plot_df_webl, 
            aes(x = x, y = predicted, color = lux_label), 
            linewidth = 1.5, alpha = 0.8) +
  annotate("text", x = 52, y = 35.7, label = "†", size = 5, color = "grey25") +
  scale_color_manual(values = int_cols, name = "Model Slices") +
  scale_fill_manual(values = int_cols) +
  new_scale_fill() +
  new_scale_color() +
  geom_point(data = filter(webl, !is.na(lux)),
             aes(x = LAeq, y = bc_wing, color = lux, shape = sex),   
             size = 3, alpha = 0.95) +     
  scale_color_gradientn(
    colours = c("midnightblue", "magenta4", "#FFE082"),
    values = c(0, 0.23, 1),
    guide = guide_colorbar(title = "Lux")
  ) +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") +
  labs(
    subtitle = "Interaction of Light and Noise on Body Condition",
    x = "Noise Level (LAeq)",
    y = "Body Condition" ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
    #legend.position = "right",
    legend.position = "none",
    legend.title = element_text(size = 10, face = "bold") ) +
  guides(
    fill = guide_colorbar(order = 2),
    shape = guide_legend(order = 3) )
webl_BC_int











#------------------------------------------------------------------------------
#---- Figure 2 Multi-Panel -----
# Arrange Top Row 
row1 <- ggarrange(webl_mass_lux, webl_mass_noise,
                  labels = c("A", "B"),
                  ncol = 2)

# Arrange Bottom Row 
row2 <- ggarrange(webl_BC_lux, webl_BC_noise,
                  labels = c("C", "D"),
                  ncol = 2)

# Arrange Bottom Row 
row3 <- ggarrange(webl_BC_int,
                  labels = c("E"),
                  ncol = 2)

# Add a NULL spacer row in between
figure2 <- ggarrange(row1, NULL, row2, NULL, row3,
                     nrow = 5,
                     heights = c(1, 0.05, 1, 0.05, 1))

# Add Multi-panel Figure Title
figure2_final <- annotate_figure(figure2,
                                 top = text_grob("Western Bluebird Morphology", 
                                                 color = "black", 
                                                 size = 20,
                                                 vjust = -0.1))

# Add extra margin at the top for Title
figure2_final <- figure2_final + 
  theme(plot.margin = margin(t = 8, r = 0, b = 0, l = 0)) 

# Display plot
figure2_final

# 6. Save high-res output
ggsave(
  filename = "Figure_2_WEBL_Morph_26.png", 
  plot = figure2_final,        
  width = 220,
  height = 300, 
  units = "mm", 
  dpi = 300,                  
  bg = "white"                
)




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Figure 3 - ATFL Morphology -----
#------------------------------------------------------------------------------
#---- Figure 3A - ATFL Wing Length & Lux -----
# Subset the Data to only include ATFL
atfl <- dat %>% filter(species == "ATFL") %>% droplevels()

# Does not change results, helps with residual normality
atfl <- atfl %>%
  mutate(
    mass = ifelse(ID %in% c("2921-17874", "3251-02916"), NA, mass),
    tarsus = ifelse(ID %in% c("2921-49602"), NA, tarsus)
  )

# LMER
atfl_wing_lux <- lmer(wing ~ scale(lux) * sex + (1 | measurer), data = atfl, REML = FALSE)
summary(atfl_wing_lux)
confint(atfl_wing_lux)
confint(atfl_wing_lux, level = 0.85)

#Body Condition across lux values using model residuals
predict_atfl_wing_lux <- predict_response(atfl_wing_lux, 
                                          terms = c("lux [all]", "sex"), 
                                          margin = "marginalmeans", 
                                          ci.lvl = 0.85)


# ATFL Tarsus vs Lux
# Marginal effect (the predicted mean)
atfl_wing_lux <- ggplot() +
  geom_ribbon(data = predict_atfl_wing_lux, 
              aes(x = x, ymin = conf.low, ymax = conf.high, group = group), 
              fill = "#FFE082", alpha = 0.3) +
  geom_line(data = predict_atfl_wing_lux, 
            aes(x = x, y = predicted, group = group, linetype = group), 
            color = "#FFE082", linewidth = 1.5) +
  geom_point(data = atfl, 
             aes(x = lux, y = wing, shape = sex), 
             color = "black", size = 4, alpha = 0.8) +
  annotate("text", x = 2, y = 110, label = "†", size = 5, color = "grey25") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(x = "Light Intensity (Lux)",
       y = "Wing Length (mm)"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    #legend.position = "right")
    legend.position = "none")
atfl_wing_lux








#------------------------------------------------------------------------------
#---- Figure 3B - ATFL Mass -----

# LMER
atfl_mass_lux <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + (1|stage), data = atfl, REML = FALSE)
summary(atfl_mass_lux)
confint(atfl_mass_lux)
confint(atfl_mass_lux, level = 0.85)

# Mass across lux values using model residuals
predict_atlf_mass_lux <- predict_response(atfl_mass_lux, terms = c("lux [all]", "sex"), margin = "marginalmeans", ci.lvl = 0.85)

# ATFL Mass vs Lux
# Marginal effect (the predicted mean)
atfl_mass_lux <- ggplot() +
  geom_ribbon(data = predict_atlf_mass_lux, 
              aes(x = x, ymin = conf.low, ymax = conf.high, group = group), 
              fill = "#FFE082", alpha = 0.3) +
  geom_line(data = predict_atlf_mass_lux, 
            aes(x = x, y = predicted, group = group, linetype = group), 
            color = "#FFE082", linewidth = 1.5) +
  geom_point(data = atfl, 
             aes(x = lux, y = mass, shape = sex), 
             color = "black", size = 4, alpha = 0.8) +
  annotate("text", x = 2, y = 30.5, label = "†", size = 5, color = "grey25") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(
    x = "Light Intensity (Lux)",
    y = "Mass (g)"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    #legend.position = "right")
    legend.position = "none")

atfl_mass_lux







#------------------------------------------------------------------------------
#---- Figure 3D - ATFL BC Lux -----
# Calculate Body Condition using mass and wing
BCwing <- smi(atfl$mass, atfl$wing)
atfl$bc_wing <- BCwing

# Does not change results, helps with residual normality
atfl <- atfl %>%
  mutate(
    bc_wing = ifelse(ID %in% c("2921-49611"), NA, bc_wing)
  )

# LMER
atfl_bc_lux <- lmer(bc_wing ~ scale(lux) * sex + (1 | section) + (1 | stage) + (1 | measurer), data = atfl, REML = FALSE)
summary(atfl_bc_lux)
confint(atfl_bc_lux)
confint(atfl_bc_lux, level = 0.85)

# Body Condition across lux values using model residuals
predict_atfl_lux_sex <- predict_response(atfl_bc_lux, 
                                         terms = c("lux [all]", "sex"), 
                                         margin = "marginalmeans", 
                                         ci.lvl = 0.85)


# ATFL BC vs Lux
# Marginal effect (the predicted mean)
atfl_BC_lux <- ggplot() +
  geom_ribbon(data = predict_atfl_lux_sex, 
              aes(x = x, ymin = conf.low, ymax = conf.high, group = group), 
              fill = "#FFE082", alpha = 0.3) +
  geom_line(data = predict_atfl_lux_sex, 
            aes(x = x, y = predicted, group = group, linetype = group), 
            color = "#FFE082", linewidth = 1.5) +
  geom_point(data = atfl, 
             aes(x = lux, y = bc_wing, shape = sex), 
             color = "black", size = 4, alpha = 0.8) +
  annotate("text", x = 2, y = 32.7, label = "*", size = 10, color = "grey25") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(
    x = "Light Intensity (Lux)",
    y = "Body Condition"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    #legend.position = "right")
    legend.position = "none")
atfl_BC_lux










#------------------------------------------------------------------------------
#---- Figure 3 Multi-Panel -----
# Combining figures
figure3 <- ggarrange(
  atfl_wing_lux, atfl_mass_lux, atfl_BC_lux, 
  labels = c("A", "B", "C"), 
  ncol = 3, 
  nrow = 1
)

figure3

# Add Multi-panel Figure Title
figure3_final <- annotate_figure(figure3,
                                 top = text_grob("Ash-throated Flycatcher Morphology", 
                                                 color = "black", 
                                                 #face = "bold", 
                                                 size = 20,
                                                 vjust = -0.1))
figure3_final

# Add extra margin at the top for Title
figure3_final <- figure3_final + 
  theme(plot.margin = margin(t = 8, r = 0, b = 0, l = 0)) 
figure3_final



ggsave(
  filename = "Figure_3_ATFL_Morph_26.png", 
  plot = last_plot(),         
  width = 330,
  height = 100, 
  units = "mm", 
  dpi = 300,                  
  bg = "white"                
)

ggsave(
  filename = "Figure_3_Legend.png", 
  plot = last_plot(),         
  width = 225, 
  height = 300, 
  units = "mm", 
  dpi = 300,                  
  bg = "white"                
)






#------------------------------------------------------------------------------
#---- Figure S1 - WEBL BC Int SEPARATED BY SEX -----

# Base off of same model in Figure 2B
# LMER
webl_BC_int <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + (1 | section) + (1 | stage), data = webl, REML = FALSE)
summary(webl_BC_int)
confint(webl_BC_int)
confint(webl_BC_int, level = 0.85)


# Add sex to the terms argument
inter_data_webl_sex <- predict_response(webl_BC_int, terms = c("LAeq [all]", "lux [-1, 0, 1]", "sex"))

# Convert to data frame
plot_df_webl_sex <- as.data.frame(inter_data_webl_sex)

# Format the Lux labels (from the 'group' column)
plot_df_webl_sex$lux_label <- factor(plot_df_webl_sex$group, 
                                 levels = c("-1", "0", "1"),
                                 labels = c("Low Lux (-1 SD)", "Mean Lux", "High Lux (+1 SD)"))

# Format the Sex labels (from the 'facet' column)
plot_df_webl_sex$sex_label <- factor(plot_df_webl_sex$facet, levels = c("Female", "Male"))


# Anchor colors
low_lux_col  <- "midnightblue"  
mean_lux_col <- "magenta4" 
high_lux_col <- "#FFE082"

int_cols <- c(
  "Low Lux (-1 SD)"  = low_lux_col, 
  "Mean Lux"         = mean_lux_col, 
  "High Lux (+1 SD)" = high_lux_col
)

webl_BC_int <- ggplot() +
  geom_ribbon(data = plot_df_webl_sex, 
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = lux_label, 
                  group = interaction(lux_label, sex_label)), 
              alpha = 0.2, show.legend = FALSE) + 
  geom_line(data = plot_df_webl_sex, 
            aes(x = x, y = predicted, color = lux_label, linetype = sex_label), 
            linewidth = 1.5, alpha = 0.8) +
    scale_color_manual(values = int_cols, name = "Model Slices") +
  scale_fill_manual(values = int_cols) +
  scale_linetype_manual(values = c("Female" = "solid", "Male" = "dashed"), guide = "none") +
  new_scale_fill() +
  new_scale_color() +
    geom_point(data = webl %>% 
               filter(!is.na(lux)) %>% 
               mutate(sex_label = sex),
             aes(x = LAeq, y = bc_wing, color = lux, shape = sex),   
             size = 3, alpha = 0.95) +     
    scale_color_gradientn(
    colours = c("midnightblue", "magenta4", "#FFE082"),
    values = c(0, 0.23, 1),
    guide = guide_colorbar(title = "Lux")
  ) +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") +
    labs(
    subtitle = "Interaction of Light and Noise on Body Condition by Sex",
    x = "Noise Level (LAeq)",
    y = "Body Condition" ) +
    theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
    legend.position = "right", 
    legend.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"), 
    strip.background = element_blank()                   
  ) +
  guides(
    fill = guide_colorbar(order = 2),
    shape = guide_legend(order = 3)
  ) +
  facet_wrap(~ sex_label)

webl_BC_int



ggsave(
  filename = "Figure_S1_WEBL_BC_sex.png", 
  plot = last_plot(),         
  width = 250, 
  height = 125,  
  units = "mm", 
  dpi = 300,                  
  bg = "white"                
)



