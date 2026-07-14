#####################
# 2018-2025 Adult Morphology on Clutch Initiation Phenology 

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(lubridate)    # to transform dates
library(AICcmodavg)   # model selection
library(emmeans)
library(ggeffects)
library(lubridate)
library(ggpubr)
library(ggnewscale)


#---- Nest Building and Adult Size
dat <- read.csv("Data/18-25_Adult_Morphology_Cavity.csv")

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
    LAeq_N = as.numeric(LAeq_N),
    mass = as.numeric(mass),
    tarsus = as.numeric(tarsus),
    wing = as.numeric(wing),
    tail = as.numeric(tail),
    eye = as.numeric(eye),
    measurer = as.factor(measurer),
    nb_date = as.Date(nb_date, "%m/%d/%y"),
    nb_julian = yday(nb_date),
    egg_date = as.Date(egg_date, "%m/%d/%y"),
    egg_julian = yday(egg_date)
  )

# Reorder the levels of the age factor
dat <- dat %>%
  mutate(age = factor(age, levels = c("SY", "AHY", "ASY")))


#---- Body Condition Function ----
# Body Condition Function with mass and tarsus
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
#---- Mountain Bluebird ----

# Subset the Data to only include MOBL
mobl <- dat %>% filter(species == "MOBL") %>% droplevels()

# Determine # of nests that have egg laying date and adult
mobl <- mobl %>% filter(!is.na(egg_date))
  # 32 nests have an egg laying date

# Remove nests that are clearly 2nd clutches
mobl <- mobl[!(mobl$nest %in% c("PC13B1-nest2")), ]
  # 31 nests wtih eggs
  # Only 1 individual without a wing measurement


#------------------------------------------------------------------------------
#---- MOBL Egg Laying with Fixed effects for Adult dataset
# Egg Laying across years
ggplot(mobl, aes(x = egg_julian, y = year, colour = year)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "MOBL Egg Laying Dates across Years",x = "Egg Laying Dates per Nest", y = "Year")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

egg_year <- lm(egg_julian ~ year, data = mobl)
summary(egg_year)
anova(egg_year) # sig
post_hoc <- emmeans(egg_year, ~ year)
pairs(post_hoc)
plot(post_hoc)

# Egg Laying dates across sections
custom_colors <- c("green3", "orange", "red", "yellow2")
ggplot(mobl, aes(x = egg_julian, y = section, colour = section)) +
  geom_jitter(size = 3, position = position_jitter(width = 0, height = 0.2)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  scale_color_manual(values = custom_colors) +
  labs(title = "MOBL Egg Laying across Sections",x = "Egg Laying Dates per Nest", y = "Section")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

egg_section <- lm(egg_julian ~ section, data = mobl)
summary(egg_section)
anova(egg_section) # Differs at 85% CI


# Egg Laying dates across age
ggplot(mobl, aes(x = egg_julian, y = age, colour = age)) +
  geom_jitter(size = 2, position = position_jitter(height = 0.3)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "Egg Laying Dates across Age",x = "Egg Laying Dates per Age", y = "Age")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(egg_julian ~ age, data = mobl)
summary(age)
anova(age) # NS

# Egg Laying dates across sex
ggplot(mobl, aes(x = egg_julian, y = sex, colour = sex)) +
  geom_jitter(size = 2, position = position_jitter(height = 0.3)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "Egg Laying Dates across Sex",x = "Egg Laying Dates per Sex", y = "Sex")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(egg_julian ~ sex, data = mobl)
summary(sex)
anova(sex) # NS

# Egg Laying dates across measurere of bird
ggplot(mobl, aes(x = egg_julian, y = measurer, colour = measurer)) +
  geom_jitter(size = 2, position = position_jitter(height = 0.3)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "MOBL Egg Laying Dates across Measurer ID",x = "Egg Laying Dates", y = "Measurer ID")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

measurer <- lm(egg_julian ~ measurer, data = mobl)
summary(measurer)
anova(measurer) # NS



#------------------------------------------------------------------------------
#---- MOBL TARSUS Length Relationship with Egg Laying Date ----

# Determine # of individuals with nest building date and tarsus length
mobl %>% filter(!is.na(egg_date) & !is.na(tarsus)) %>% nrow()
  # 31 nests have a nest building date AND tarsus

# Tarsus length & Egg Laying
ggplot(mobl, aes(x = egg_julian, y = tarsus)) +
  geom_jitter(position = position_jitter(width = 0.2)) +  
  geom_smooth(method = "lm", se = FALSE, color = "blue3") +
  labs(title = "MOBL Tarsus across Egg Laying",x = "Egg Laying Dates per Nest", y = "Tarsus Length")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


# Tarsus length across sex
ggplot(mobl, aes(x = sex, y = tarsus, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "MOBL Tarsus Across Sex", x = "Sex", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(tarsus ~ sex, data = mobl)
summary(sex)
# NS

# Tarsus length across age
ggplot(mobl, aes(x = age, y = tarsus, fill = age)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "MOBL Tarsus Across Sex", x = "Sex", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(tarsus ~ age, data = mobl)
summary(age) # NS



#---- Test Random Effects for Egg Laying TARSUS Model ----
# Model with (1|year), (1|section)
model_3 <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|year) + (1|section), data = mobl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|year), data = mobl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|section), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(egg_julian ~ scale(tarsus) + sex + age, data = mobl)

# Compare all models
anova(model_3, model_2, model_1, model_0)
# Model 0 is best, 2 next

# Move Forward: LMER Year



#---- Model Selection TARSUS vs Egg Laying -----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

#### Interaction
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + sex + age + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + sex + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + age + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + sex + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + age + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + (1|year), data = mobl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(egg_julian ~ 1 + (1|year), data = mobl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[9]], cand.models[[8]])  
top_models

# Top was Null

# 2nd Top model: egg_julian ~ scale(tarsus) + (1 | year)
top_mod2 <- lmer(egg_julian ~ scale(tarsus) + (1 | year), data = mobl)
summary(top_mod2)
anova(top_mod2)
confint(top_mod2)
confint(top_mod2, level = 0.85)

# No effect of Tarsus on MOBL Egg Laying date



#------------------------------------------------------------------------------
#---- MOBL WING Length Relationship with Egg Laying Date ----

# Determine # of individuals with nest building date and tarsus length
mobl %>% filter(!is.na(egg_date) & !is.na(wing)) %>% nrow()
  # 31 nests have a nest building date AND tarsus

# Wing length & Egg Laying
ggplot(mobl, aes(x = egg_julian, y = wing)) +
  geom_jitter(position = position_jitter(width = 0.2)) +  
  geom_smooth(method = "lm", se = FALSE, color = "blue3") +
  labs(title = "MOBL Wing across Egg Laying",x = "Egg Laying Dates per Nest", y = "Wing Length")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing length across sex
ggplot(mobl, aes(x = sex, y = wing, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "MOBL Wing Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(wing ~ sex, data = mobl)
summary(sex)


# Wing length across age
ggplot(mobl, aes(x = age, y = wing, fill = age)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "MOBL Winh Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(wing ~ age, data = mobl)
summary(age) # NS



#---- Test Random Effects for Egg Laying WING Model ----
# Model with (1|year), (1|section)
model_3 <- lmer(egg_julian ~ scale(wing) + sex + age + (1|year) + (1|section), data = mobl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(egg_julian ~ scale(wing) + sex + age + (1|year), data = mobl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(egg_julian ~ scale(wing) + sex + age + (1|section), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(egg_julian ~ scale(wing) + sex + age, data = mobl)

# Compare all models
anova(model_3, model_2, model_1, model_0)
# Model 2 is best

# Move Forward: LMER Year




#---- Model Selection WING vs Egg Laying -----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

#### Interaction
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + sex + age + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + sex + age + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + sex + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + age + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + sex + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + age + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + (1|year), data = mobl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(egg_julian ~ 1 + (1|year), data = mobl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[8]], cand.models[[9]], cand.models[[6]])  
top_models

# Top model: egg_julian ~ scale(wing) + (1 | year)
top_mod <- lmer(egg_julian ~ scale(wing) + (1 | year), data = mobl)
summary(top_mod)
anova(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
  # wing length has a negative effect at 85% CI

# 2nd Top model: NULL

# 3rd Top model: egg_julian ~ scale(wing) + age + (1 | year)
top_mod3 <- lmer(egg_julian ~ scale(wing) + age + (1 | year), data = mobl)
summary(top_mod3)
anova(top_mod3)
confint(top_mod3)
confint(top_mod3, level = 0.85)
  # wing has a neg effect at 85%
  # no effect of age





#------------------------------------------------------------------------------
#---- MOBL BODY CONDITION Length Relationship with Egg Laying Date ----

# Calculate Body Condition using mass and wing
BCwing <- smi(mobl$mass, mobl$wing)
mobl$bc_wing <- BCwing

# Determine # of individuals with nest building date and tarsus length
mobl %>% filter(!is.na(egg_date) & !is.na(bc_wing)) %>% nrow()
  # 31 nests have a nest building date AND tarsus

# BC length & Egg Laying
ggplot(mobl, aes(x = egg_julian, y = bc_wing)) +
  geom_jitter(position = position_jitter(width = 0.2)) +  
  geom_smooth(method = "lm", se = FALSE, color = "blue3") +
  labs(title = "MOBL BC across Egg Laying",x = "Egg Laying Dates per Nest", y = "BC")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# BC length across sex
ggplot(mobl, aes(x = sex, y = bc_wing, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "MOBL BC Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(bc_wing ~ sex, data = mobl)
summary(sex)


# BC length across age
ggplot(mobl, aes(x = age, y = bc_wing, fill = age)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "MOBL BC Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(bc_wing ~ age, data = mobl)
summary(age) # NS



#---- Test Random Effects for Egg Laying WING Model ----
# Model with (1|year), (1|section)
model_3 <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|year) + (1|section), data = mobl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|year), data = mobl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|section), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(egg_julian ~ scale(bc_wing) + sex + age, data = mobl)

# Compare all models
anova(model_3, model_2, model_1, model_0)
# Model 2 is best

# Move Forward: LMER Year




#---- Model Selection WING vs Egg Laying -----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

#### Interaction
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + sex + age + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + sex + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + age + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + sex + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + age + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + section + (1|year), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + (1|year), data = mobl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(egg_julian ~ 1 + (1|year), data = mobl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[8]], cand.models[[9]], cand.models[[5]])  
top_models

# Top model: egg_julian ~ scale(bc_wing) + (1 | year)
top_mod <- lmer(egg_julian ~ scale(bc_wing) + (1 | year), data = mobl)
summary(top_mod)
anova(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
  # wing length has a negative effect at 85% CI

# 2nd Top model: NULL

# 3rd Top model: egg_julian ~ scale(bc_wing) + sex + (1 | year)
top_mod3 <- lmer(egg_julian ~ scale(bc_wing) + sex + (1 | year), data = mobl)
summary(top_mod3)
anova(top_mod3)
confint(top_mod3)
confint(top_mod3, level = 0.85)
  # wing has a SIG neg effect
# no effect of age





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Western Bluebird ----

# Subset the Data to only include MOBL
webl <- dat %>% filter(species == "WEBL") %>% droplevels()

# Determine # of nests that have egg laying date and adult
webl <- webl %>% filter(!is.na(egg_date))
  # 44 nests have an egg laying date

# Remove nests that are clearly 2nd clutches
#webl <- webl[!(webl$nest %in% c("TP01B06")), ]
webl <- webl[!(webl$nest %in% c("TP01B06") | webl$ID == "3011-18045"), ]
  # 42 nests wtih eggs
  # Only 1 individual without a wing measurement



#------------------------------------------------------------------------------
#---- WEBL Egg Laying with Fixed effects for Adult dataset
# Egg Laying across years
ggplot(webl, aes(x = egg_julian, y = year, colour = year)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "WEBL Egg Laying Dates across Years",x = "Egg Laying Dates per Nest", y = "Year")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

egg_year <- lm(egg_julian ~ year, data = webl)
summary(egg_year)
anova(egg_year)
post_hoc <- emmeans(egg_year, ~ year)
pairs(post_hoc)
plot(post_hoc)

# Egg Laying dates across sections
custom_colors <- c("green3", "blue", "orange", "pink", "red", "yellow2")
ggplot(webl, aes(x = egg_julian, y = section, colour = section)) +
  geom_jitter(size = 3, position = position_jitter(width = 0, height = 0.2)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  scale_color_manual(values = custom_colors) +
  labs(title = "WEBL Egg Laying across Sections",x = "Egg Laying Dates per Nest", y = "Section")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

egg_section <- lm(egg_julian ~ section, data = webl)
summary(egg_section)
anova(egg_section) # Differs at 85% CI


# Egg Laying dates across age
ggplot(webl, aes(x = egg_julian, y = age, colour = age)) +
  geom_jitter(size = 2, position = position_jitter(height = 0.3)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "WEBL Egg Laying Dates across Age",x = "Egg Laying Dates per Age", y = "Age")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(egg_julian ~ age, data = webl)
summary(age)
anova(age) # 85% Effect

# Egg Laying dates across sex
ggplot(webl, aes(x = egg_julian, y = sex, colour = sex)) +
  geom_jitter(size = 2, position = position_jitter(height = 0.3)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "WEBL Egg Laying Dates across Sex",x = "Egg Laying Dates per Sex", y = "Sex")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(egg_julian ~ sex, data = webl)
summary(sex)
anova(sex) # ASY and SY are SIG diff

# Egg Laying dates across measurere of bird
ggplot(webl, aes(x = egg_julian, y = measurer, colour = measurer)) +
  geom_jitter(size = 2, position = position_jitter(height = 0.3)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "WEBL Egg Laying Dates across Measurer ID",x = "Egg Laying Dates", y = "Measurer ID")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

measurer <- lm(egg_julian ~ measurer, data = webl)
summary(measurer)
anova(measurer) # NS



#------------------------------------------------------------------------------
#---- WEBL TARSUS Length Relationship with Egg Laying Date ----

# Determine # of individuals with nest building date and tarsus length
webl %>% filter(!is.na(egg_date) & !is.na(tarsus)) %>% nrow()
  # 39 nests have a nest building date AND tarsus

# Tarsus length & Egg Laying
ggplot(webl, aes(x = egg_julian, y = tarsus)) +
  geom_jitter(position = position_jitter(width = 0.2)) +  
  geom_smooth(method = "lm", se = FALSE, color = "blue3") +
  labs(title = "WEBL Tarsus across Egg Laying",x = "Egg Laying Dates per Nest", y = "Tarsus Length")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


# Tarsus length across sex
ggplot(webl, aes(x = sex, y = tarsus, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "WEBL Tarsus Across Sex", x = "Sex", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(tarsus ~ sex, data = webl)
summary(sex)
# NS

# Tarsus length across age
ggplot(webl, aes(x = age, y = tarsus, fill = age)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "WEBL Tarsus Across Sex", x = "Sex", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(tarsus ~ age, data = webl)
summary(age) # NS



#---- Test Random Effects for Egg Laying TARSUS Model ----
# Model with (1|year), (1|section)
model_3 <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|year) + (1|section), data = webl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|year), data = webl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|section), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(egg_julian ~ scale(tarsus) + sex + age, data = webl)

# Compare all models
anova(model_3, model_2, model_1, model_0)
# Model 0 is best

# Move Forward: LM



#---- Model Selection TARSUS vs Egg Laying -----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

#### Interaction
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + sex + age + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + sex + age + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + sex + age + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + age + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + sex + age, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + sex + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + sex + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + age + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + age + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + sex, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + age, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus) + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(tarsus), data = webl); y <- y + 1
# Null
cand.models[[y]] <- lm(egg_julian ~ 1, data = webl); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[12]], cand.models[[15]])  
top_models


# Top model: egg_julian ~ scale(tarsus) + age
top_mod <- lm(egg_julian ~ scale(tarsus) + age, data = webl)
summary(top_mod)
anova(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
  # No effect of tarsus
  # SIG effect of ASY 

# 2nd Top model: egg_julian ~ scale(tarsus)
top_mod2 <- lm(egg_julian ~ scale(tarsus), data = webl)
summary(top_mod2)
anova(top_mod2)
confint(top_mod2)
confint(top_mod2, level = 0.85)
  # No effect



#------------------------------------------------------------------------------
#---- WING Length Relationship with Egg Laying Date ----

# Determine # of individuals with nest building date and tarsus length
webl %>% filter(!is.na(egg_date) & !is.na(wing)) %>% nrow()
  # 39 nests have a nest building date AND tarsus

# Wing length & Egg Laying
ggplot(webl, aes(x = egg_julian, y = wing)) +
  geom_jitter(position = position_jitter(width = 0.2)) +  
  geom_smooth(method = "lm", se = FALSE, color = "blue3") +
  labs(title = "WEBL Wing across Egg Laying",x = "Egg Laying Dates per Nest", y = "Wing Length")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing length across sex
ggplot(webl, aes(x = sex, y = wing, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "WEBL Wing Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(wing ~ sex, data = webl)
summary(sex) # SIG


# Wing length across age
ggplot(webl, aes(x = age, y = wing, fill = age)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "WEBL Winh Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(wing ~ age, data = webl)
summary(age) 
  # ASY has an effect at 85%



#---- Test Random Effects for Egg Laying WING Model ----
# Model with (1|year), (1|section)
model_3 <- lmer(egg_julian ~ scale(wing) + sex + age + (1|year) + (1|section), data = webl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(egg_julian ~ scale(wing) + sex + age + (1|year), data = webl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(egg_julian ~ scale(wing) + sex + age + (1|section), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(egg_julian ~ scale(wing) + sex + age, data = webl)

# Compare all models
anova(model_3, model_2, model_1, model_0)
# Model 0

# Move Forward: LM



#---- Model Selection WING vs Egg Laying -----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

#### Interaction
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + sex + age + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + sex + age + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + sex + age + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + sex + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + age + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + sex + age, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + sex + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + sex + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + age + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + age + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + sex, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + age, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing) + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(wing), data = webl); y <- y + 1
# Null
cand.models[[y]] <- lm(egg_julian ~ 1, data = webl); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[13]], cand.models[[16]], cand.models[[12]], cand.models[[6]])  
top_models


# Top model: egg_julian ~ scale(wing) + age
top_mod <- lm(egg_julian ~ scale(wing) + age, data = webl)
summary(top_mod)
anova(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
  # No effect of wing 
  # ASY differs from SY at 85% 

# 2nd Top model: egg_julian ~ scale(wing)
top_mod2 <- lm(egg_julian ~ scale(wing), data = webl)
summary(top_mod2)
anova(top_mod2)
confint(top_mod2)
confint(top_mod2, level = 0.85)
  # no effect of wing

# 3rd Top model: egg_julian ~ scale(wing) + sex
top_mod3 <- lm(egg_julian ~ scale(wing) + sex , data = webl)
summary(top_mod3)
anova(top_mod3)
confint(top_mod3)
confint(top_mod3, level = 0.85)
  # wing length has a negative effect at 85% CI
  # no effect of sex

# 4th Top model: egg_julian ~ scale(wing) + sex
top_mod4 <- lm(egg_julian ~ scale(wing) + sex + age , data = webl)
summary(top_mod4)
anova(top_mod4)
confint(top_mod4)
confint(top_mod4, level = 0.85)
# No effect of wing or sex
#   # effect of age





#------------------------------------------------------------------------------
#---- WEBL BODY CONDITION Length Relationship with Egg Laying Date ----

# Calculate Body Condition using mass and wing
BCwing <- smi(webl$mass, webl$wing)
webl$bc_wing <- BCwing

# Determine # of individuals with nest building date and tarsus length
webl %>% filter(!is.na(egg_date) & !is.na(bc_wing)) %>% nrow()
  # 37 nests have a nest building date AND tarsus

# BC length & Egg Laying
ggplot(webl, aes(x = egg_julian, y = bc_wing)) +
  geom_jitter(position = position_jitter(width = 0.2)) +  
  geom_smooth(method = "lm", se = FALSE, color = "blue3") +
  labs(title = "WEBL BC across Egg Laying",x = "Egg Laying Dates per Nest", y = "BC")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# BC length across sex
ggplot(webl, aes(x = sex, y = bc_wing, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "WEBL BC Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(bc_wing ~ sex, data = webl)
summary(sex)

# BC length across age
ggplot(webl, aes(x = age, y = bc_wing, fill = age)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "WEBL BC Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(bc_wing ~ age, data = webl)
summary(age) # ASY is SIG diff from SY



#---- Test Random Effects for Egg Laying WING Model ----
# Model with (1|year), (1|section)
model_3 <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|year) + (1|section), data = webl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|year), data = webl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|section), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(egg_julian ~ scale(bc_wing) + sex + age, data = webl)

# Compare all models
anova(model_3, model_2, model_1, model_0)
  # Model 0 is best, then 1

# Move Forward: LM



#---- Model Selection WING vs Egg Laying -----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

#### Interaction
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + sex + age + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + sex + age + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + sex + age + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + sex + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + age + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + sex + age, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + sex + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + sex + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + age + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + age + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + section + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + sex, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + age, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + section, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing) + year, data = webl); y <- y + 1
cand.models[[y]] <- lm(egg_julian ~ scale(bc_wing), data = webl); y <- y + 1
# Null
cand.models[[y]] <- lm(egg_julian ~ 1, data = webl); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[6]], cand.models[[12]], cand.models[[16]])  
top_models

# Top model: egg_julian ~ scale(bc_wing) + sex + age
top_mod <- lm(egg_julian ~ scale(bc_wing) + sex + age, data = webl)
summary(top_mod)
anova(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
  # Condition has a SIG neg effect
  # Sex was SIG positive
  # ASY was SIG neg

# 2nd Top model: egg_julian ~ scale(bc_wing) + sex
top_mod2 <- lm(egg_julian ~ scale(bc_wing) + sex, data = webl)
summary(top_mod2)
anova(top_mod2)
confint(top_mod2)
confint(top_mod2, level = 0.85)
  # Condition has a SIG neg effect
  # sex had a pos effect at 85%

# 3rd Top model: egg_julian ~ scale(bc_wing)
top_mod3 <- lm(egg_julian ~ scale(bc_wing), data = webl)
summary(top_mod3)
anova(top_mod3)
confint(top_mod3)
confint(top_mod3, level = 0.85)
  # Condition has a SIG neg effect




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Ash-throated Flycatcher ----

# Subset the Data to only include MOBL
atfl <- dat %>% filter(species == "ATFL") %>% droplevels()

# Determine # of nests that have egg laying date and adult
atfl <- atfl %>% filter(!is.na(egg_date))
  # 32 nests have an egg laying date

# Remove mass outlier
# Does not change results, helps with residual normality
atfl <- atfl %>%
  mutate(mass = ifelse(ID %in% c("2921-17874"), NA, mass))




#------------------------------------------------------------------------------
#---- ATFL Egg Laying with Fixed effects for Adult dataset
# Egg Laying across years
ggplot(atfl, aes(x = egg_julian, y = year, colour = year)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "ATFL Egg Laying Dates across Years",x = "Egg Laying Dates per Nest", y = "Year")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

egg_year <- lm(egg_julian ~ year, data = atfl)
summary(egg_year)
anova(egg_year)
post_hoc <- emmeans(egg_year, ~ year)
pairs(post_hoc)
plot(post_hoc)

# Egg Laying dates across sections
custom_colors <- c("green3", "blue", "orange", "pink", "red", "yellow2")
ggplot(atfl, aes(x = egg_julian, y = section, colour = section)) +
  geom_jitter(size = 3, position = position_jitter(width = 0, height = 0.2)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  scale_color_manual(values = custom_colors) +
  labs(title = "ATFL Egg Laying across Sections",x = "Egg Laying Dates per Nest", y = "Section")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

egg_section <- lm(egg_julian ~ section, data = atfl)
summary(egg_section)
anova(egg_section) # SIG


# Egg Laying dates across age
ggplot(atfl, aes(x = egg_julian, y = age, colour = age)) +
  geom_jitter(size = 2, position = position_jitter(height = 0.3)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "ATFL Egg Laying Dates across Age",x = "Egg Laying Dates per Age", y = "Age")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(egg_julian ~ age, data = atfl)
summary(age)
anova(age) 

# Egg Laying dates across sex
ggplot(atfl, aes(x = egg_julian, y = sex, colour = sex)) +
  geom_jitter(size = 2, position = position_jitter(height = 0.3)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "ATFL Egg Laying Dates across Sex",x = "Egg Laying Dates per Sex", y = "Sex")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(egg_julian ~ sex, data = atfl)
summary(sex)
anova(sex) 

# Egg Laying dates across measurere of bird
ggplot(atfl, aes(x = egg_julian, y = measurer, colour = measurer)) +
  geom_jitter(size = 2, position = position_jitter(height = 0.3)) +  
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.3,
               position = position_dodge(width = 0.5)) +
  labs(title = "ATFL Egg Laying Dates across Measurer ID",x = "Egg Laying Dates", y = "Measurer ID")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

measurer <- lm(egg_julian ~ measurer, data = atfl)
summary(measurer)
anova(measurer) # NS



#------------------------------------------------------------------------------
#---- ATFL TARSUS Length Relationship with Egg Laying Date ----

# Determine # of individuals with nest building date and tarsus length
atfl %>% filter(!is.na(egg_date) & !is.na(tarsus)) %>% nrow()
  # 32 nests have a nest building date AND tarsus

# Tarsus length & Egg Laying
ggplot(atfl, aes(x = egg_julian, y = tarsus)) +
  geom_jitter(position = position_jitter(width = 0.2)) +  
  geom_smooth(method = "lm", se = FALSE, color = "blue3") +
  labs(title = "ATFL Tarsus across Egg Laying",x = "Egg Laying Dates per Nest", y = "Tarsus Length")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus length across sex
ggplot(atfl, aes(x = sex, y = tarsus, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "ATFL Tarsus Across Sex", x = "Sex", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(tarsus ~ sex, data = atfl)
summary(sex)

# Tarsus length across age
ggplot(atfl, aes(x = age, y = tarsus, fill = age)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "ATFL Tarsus Across Sex", x = "Sex", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(tarsus ~ age, data = atfl)
summary(age) # NS



#---- Test Random Effects for Egg Laying TARSUS Model ----
# Model with (1|year), (1|section)
model_3 <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|year) + (1|section), data = atfl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|year), data = atfl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|section), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(egg_julian ~ scale(tarsus) + sex + age, data = atfl)

# Compare all models
anova(model_3, model_2, model_1, model_0)
# Model 1

# Move Forward: LMER - section



#---- Model Selection TARSUS vs Egg Laying -----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

#### Interaction
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + sex + age + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + sex + age + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + sex + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + age + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + sex + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + age + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(tarsus) + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(egg_julian ~ 1 + (1|section), data = atfl, REML = FALSE); y <- y + 1


# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[8]])  
top_models

# Top model: Null




#------------------------------------------------------------------------------
#---- WING Length Relationship with Egg Laying Date ----

# Determine # of individuals with nest building date and tarsus length
atfl %>% filter(!is.na(egg_date) & !is.na(wing)) %>% nrow()
  # 32 nests have a nest building date AND tarsus

# Wing length & Egg Laying
ggplot(atfl, aes(x = egg_julian, y = wing)) +
  geom_jitter(position = position_jitter(width = 0.2)) +  
  geom_smooth(method = "lm", se = FALSE, color = "blue3") +
  labs(title = "ATFL Wing across Egg Laying",x = "Egg Laying Dates per Nest", y = "Wing Length")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing length across sex
ggplot(atfl, aes(x = sex, y = wing, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "ATFL Wing Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(wing ~ sex, data = atfl)
summary(sex) # SIG


# Wing length across age
ggplot(atfl, aes(x = age, y = wing, fill = age)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "ATFL Winh Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(wing ~ age, data = atfl)
summary(age) 



#---- Test Random Effects for Egg Laying WING Model ----
# Model with (1|year), (1|section)
model_3 <- lmer(egg_julian ~ scale(wing) + sex + age + (1|year) + (1|section), data = atfl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(egg_julian ~ scale(wing) + sex + age + (1|year), data = atfl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(egg_julian ~ scale(wing) + sex + age + (1|section), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(egg_julian ~ scale(wing) + sex + age, data = atfl)

# Compare all models
anova(model_3, model_2, model_1, model_0)
# Model 1

# Move Forward: LMER - section




#---- Model Selection WING vs Egg Laying -----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

#### Interaction
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + sex + age + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + sex + age + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + sex + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + age + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + sex + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + age + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(wing) + (1|section), data = atfl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(egg_julian ~ 1 + (1|section), data = atfl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[9]])  
top_models

# Top model: Null






#------------------------------------------------------------------------------
#---- WEBL BODY CONDITION Length Relationship with Egg Laying Date ----

# Calculate Body Condition using mass and wing
BCwing <- smi(atfl$mass, atfl$wing)
atfl$bc_wing <- BCwing

# Determine # of individuals with nest building date and tarsus length
atfl %>% filter(!is.na(egg_date) & !is.na(bc_wing)) %>% nrow()
  # 31 nests have a nest building date AND tarsus

# BC length & Egg Laying
ggplot(atfl, aes(x = egg_julian, y = bc_wing)) +
  geom_jitter(position = position_jitter(width = 0.2)) +  
  geom_smooth(method = "lm", se = FALSE, color = "blue3") +
  labs(title = "ATFL BC across Egg Laying",x = "Egg Laying Dates per Nest", y = "BC")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# BC length across sex
ggplot(atfl, aes(x = sex, y = bc_wing, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "ATFL BC Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

sex <- lm(bc_wing ~ sex, data = atfl)
summary(sex)

# BC length across age
ggplot(atfl, aes(x = age, y = bc_wing, fill = age)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "ATFL BC Across Sex", x = "Sex", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

age <- lm(bc_wing ~ age, data = atfl)
summary(age) # ASY is SIG diff from SY



#---- Test Random Effects for Egg Laying WING Model ----
# Model with (1|year), (1|section)
model_3 <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|year) + (1|section), data = atfl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|year), data = atfl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|section), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(egg_julian ~ scale(bc_wing) + sex + age, data = atfl)

# Compare all models
anova(model_3, model_2, model_1, model_0)
# Model 1 - SECTION

# Move Forward: LM




#---- Model Selection WING vs Egg Laying -----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

#### Interaction
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + sex + age + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + sex + age + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + sex + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + age + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + sex + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + age + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + year + (1|section), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(egg_julian ~ scale(bc_wing) + (1|section), data = atfl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(egg_julian ~ 1 + (1|section), data = atfl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[8]])  
top_models

# Top model: egg_julian ~ scale(bc_wing) + (1 | section)
top_mod <- lmer(egg_julian ~ scale(bc_wing) + (1 | section), data = atfl)
summary(top_mod)
anova(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
  # No effect of BC or sex



