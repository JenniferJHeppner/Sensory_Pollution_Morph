
#####################
# Analysis of Adult MOBL (Mountain Bluebird) Morphology & Condition 
# Wing length, tarsus length, body mass, body condition

# Load Packages ----
library(tidyverse)
library(MuMIn)        # model selection
library(lme4)
library(lmerTest)
library(ggplot2)
library(lubridate)    # to transform dates
library(AICcmodavg)   # model selection
library(emmeans)
library(car)          # testing normality
library(performance)  # heteroscedasticity check
library(DHARMa)       # Check Dispersion
library(ggeffects)    # Predicts Model Trends

# Load "Adult_Morphology_Cavity.csv" 
# Hereon, data file is named "dat"
dat <- read.csv("Data/Adult_Morphology_Cavity.csv")

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
    LAeq = as.numeric(LAeq),
    mass = as.numeric(mass),
    tarsus = as.numeric(tarsus),
    wing = as.numeric(wing),
    measurer = as.factor(measurer))


#---- Rename groups ----
levels(dat$trt)
levels(dat$trt) <- c("Control", "Light", "Combined", "Noise")
levels(dat$trt)
levels(dat$sex)
levels(dat$sex) <- c("Female", "Male")
levels(dat$sex)


#---- Reorder levels ----
# Reorder Treatments
dat <- dat %>%
  mutate(trt = fct_relevel(trt, "Control", "Light", "Noise", "Combined"))

# Reorder Ages
dat <- dat %>%
  mutate(age = factor(age, levels = c("SY", "AHY", "ASY")))


#---- Define colors ----
trt_cols <- c("Control" = "steelblue", "Noise" = "firebrick",
              "Light" = "#FFE082", "Combined" = "coral")


#---- Subset Data to only include MOBL ----
mobl <- dat %>% filter(species == "MOBL") %>% droplevels()


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

# Calculate Body Condition using mass and wing
BCwing <- smi(mobl$mass, mobl$wing)
mobl$bc_wing <- BCwing






#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- WING Exploratory Graphs and Sample Sizes -----

#------------------------------------------------------------------------------
#---- WING - Random Effects Model Selection ----

##  Macro-habitat Treatment Models
# Full model with random effects (1|year) (1|measurer) and (1|section)
model_7 <- lmer(wing ~ trt + sex + julian + (1|section) + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) and (1|measurer)
model_6 <- lmer(wing ~ trt + sex + julian + (1|section) + (1|year), data = mobl, REML = FALSE)
# Model with (1|year) and (1|section)
model_5 <- lmer(wing ~ trt + sex + julian + (1|section) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|measurer) and (1|section)
model_4 <- lmer(wing ~ trt + sex + julian + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with only (1|year)
model_3 <- lmer(wing ~ trt + sex + julian + (1|section), data = mobl, REML = FALSE)
# Model with only (1|measurer)
model_2 <- lmer(wing ~ trt + sex + julian + (1|year), data = mobl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(wing ~ trt + sex + julian + (1|measurer), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(wing ~ trt + sex + julian, data = mobl)

# Compare all model AICc
cand_rand <- list(
  model_7 = model_7, 
  model_6 = model_6, 
  model_5 = model_5, 
  model_4 = model_4, 
  model_3 = model_3, 
  model_2 = model_2, 
  model_1 = model_1, 
  model_0 = model_0)
model.sel(cand_rand) 


## Micro-habitat Lux and Noise Models
# Full model with random effects (1|year) (1|measurer) and (1|section)
model_7 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) and (1|measurer)
model_6 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year), data = mobl, REML = FALSE)
# Model with (1|year) and (1|section)
model_5 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|measurer) and (1|section)
model_4 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with only (1|year)
model_3 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section), data = mobl, REML = FALSE)
# Model with only (1|measurer)
model_2 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year), data = mobl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|measurer), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian, data = mobl)

# Compare all model AICc
cand_rand <- list(
  model_7 = model_7, 
  model_6 = model_6, 
  model_5 = model_5, 
  model_4 = model_4, 
  model_3 = model_3, 
  model_2 = model_2, 
  model_1 = model_1, 
  model_0 = model_0)
model.sel(cand_rand) 


##-- Optimized Random Effect Structure: Model 0 Linear model



#------------------------------------------------------------------------------
#---- WING - Model Selection on Fixed Effects ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lm(wing ~ trt + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt, data = mobl); y <- y + 1
# lux * noise
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq), data = mobl); y <- y + 1
# lux + noise
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq), data = mobl); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lm(wing ~ scale(lux), data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(LAeq), data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(LAeq) * sex, data = mobl); y <- y + 1
# Null
cand.models[[y]] <- lm(wing ~ 1, data = mobl); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models (models with ΔAIC < 2)
top_models <- list(cand.models[[2]])  
top_models

# Top model:  wing ~ trt + sex
top_mod <- lm(wing ~ trt + sex, data = mobl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)

# Make Reference group Light to see differences
mobl$trt <- relevel(factor(mobl$trt), ref = "Noise")

# Top model with Light as Reference group:  wing ~ trt + sex
top_mod <- lm(wing ~ trt + sex, data = mobl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)

# Remake Control the Reference group
mobl$trt <- relevel(factor(mobl$trt), ref = "Control")



#---- Assumption Tests
#### Normality Test
qqPlot(residuals(top_mod)) 
plot(density(resid(top_mod)))
shapiro.test(residuals(top_mod));length(residuals(top_mod)) # p-value > 0.05 indicates normal residuals

#### Homogeneity of variance 
# Test for Heteroscedasticity (performance)
check_heteroscedasticity(top_mod) # p > 0.05 indicates homoscedasticity (constant variance)

#### Residual Diagnostics & Dispersion (DHARMa)
# Test for over/under-dispersion and residual patterns
simulation_output <- simulateResiduals(top_mod)
plot(simulation_output)
testDispersion(simulation_output)  # p-value > 0.05 indicates no significant over- or under-dispersion problems




#------------------------------------------------------------------------------
#---- WING - MAIN  MANUSCRIPT FIGURE ----
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
  scale_color_manual(values = trt_cols) +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16)) + 
  labs(x = "Treatment Group",
       y = "Wing Length (mm)",
       #title = "Mountain Bluebird", 
       shape = "Sex",
       color = "Treatment") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(color = "black") ,
        #legend.position = "none",
        axis.line = element_line(linewidth = 0.8),
        plot.title = element_text(hjust = 0.5))
mobl_wing






#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Tarsus Length -----

#------------------------------------------------------------------------------
#---- TARSUS - Random Effects Model Selection ----

##  Macro-habitat Treatment Models
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(tarsus ~ trt + sex + (1|section) + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(tarsus ~ trt + sex + (1|section) + (1|year), data = mobl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(tarsus ~ trt + sex + (1|section) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(tarsus ~ trt + sex + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(tarsus ~ trt + sex + (1|section), data = mobl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(tarsus ~ trt + sex + (1|year), data = mobl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(tarsus ~ trt + sex + (1|measurer), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(tarsus ~ trt + sex , data = mobl)

# Compare all model AICc
cand_rand <- list(
  model_7 = model_7, 
  model_6 = model_6, 
  model_5 = model_5, 
  model_4 = model_4, 
  model_3 = model_3, 
  model_2 = model_2, 
  model_1 = model_1, 
  model_0 = model_0)
model.sel(cand_rand)

## Micro-habitat Lux and Noise Models
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|year), data = mobl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section), data = mobl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|year), data = mobl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|measurer), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex, data = mobl)

# Compare all model AICc
cand_rand <- list(
  model_7 = model_7, 
  model_6 = model_6, 
  model_5 = model_5, 
  model_4 = model_4, 
  model_3 = model_3, 
  model_2 = model_2, 
  model_1 = model_1, 
  model_0 = model_0)
model.sel(cand_rand) 


##-- Optimized Random Effect Structure: Model 0 Linear Model




#------------------------------------------------------------------------------
#---- TARSUS - Model Selection on Fixed Effects ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lm(tarsus ~ trt + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt, data = mobl); y <- y + 1
# lux * noise
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq), data = mobl); y <- y + 1
# lux + noise
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq), data = mobl); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lm(tarsus ~ scale(lux), data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(LAeq), data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(LAeq) * sex, data = mobl); y <- y + 1
# Null
cand.models[[y]] <- lm(tarsus ~ 1, data = mobl); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models (models with ΔAIC < 2)
top_models <- list(cand.models[[13]], cand.models[[10]])  
top_models

#  Top model: NULL

# 2nd Top model: tarsus ~ scale(LAeq) + (1|measurer)
top_mod_2 <- lm(tarsus ~ scale(LAeq), data = mobl)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)



#---- Assumption Tests
#### Normality Test
qqPlot(residuals(top_mod_2)) 
plot(density(resid(top_mod_2)))
shapiro.test(residuals(top_mod_2));length(residuals(top_mod_2)) # p-value > 0.05 indicates normal residuals

#### Homogeneity of variance 
# Test for Heteroscedasticity (performance)
check_heteroscedasticity(top_mod_2) # p > 0.05 indicates homoscedasticity (constant variance)

#### Residual Diagnostics & Dispersion (DHARMa)
# Test for over/under-dispersion and residual patterns
simulation_output <- simulateResiduals(top_mod_2)
plot(simulation_output)
testDispersion(simulation_output)  # p-value > 0.05 indicates no significant over- or under-dispersion problems







#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Body Mass -----

#------------------------------------------------------------------------------
#---- MASS - Random Effects Model Selection ----

##  Macro-habitat Treatment Models
# Full model with random effects (1|section) (1|year) (1|stage) and (1|measurer)
model_15 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|year) (1|stage)
model_14 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|section) (1|year) (1|measurer)
model_13 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|stage) (1|measurer)
model_12 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) (1|stage) (1|measurer)
model_11 <- lmer(mass ~ trt + sex + julian + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|year)
model_10 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|year), data = mobl, REML = FALSE)
# Model with (1|section) (1|stage)
model_9 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|section) (1|measurer)
model_8 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) (1|stage)
model_7 <- lmer(mass ~ trt + sex + julian + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|year) (1|measurer)
model_6 <- lmer(mass ~ trt + sex + julian + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|stage) (1|measurer)
model_5 <- lmer(mass ~ trt + sex + julian + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section)
model_4 <- lmer(mass ~ trt + sex + julian + (1|section), data = mobl, REML = FALSE)
# Model with (1|year)
model_3 <- lmer(mass ~ trt + sex + julian + (1|year), data = mobl, REML = FALSE)
# Model with (1|stage)
model_2 <- lmer(mass ~ trt + sex + julian + (1|stage), data = mobl, REML = FALSE)
# Model with (1|measurer)
model_1 <- lmer(mass ~ trt + sex + julian + (1|measurer), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(mass ~ trt + sex + julian, data = mobl)

# Compare all model AICc
cand_rand <- list(
  model_15 = model_15, 
  model_14 = model_14, 
  model_13 = model_13, 
  model_12 = model_12, 
  model_11 = model_11, 
  model_10 = model_10,
  model_9 = model_9,
  model_8 = model_8,
  model_7 = model_7, 
  model_6 = model_6, 
  model_5 = model_5, 
  model_4 = model_4, 
  model_3 = model_3, 
  model_2 = model_2, 
  model_1 = model_1, 
  model_0 = model_0)
model.sel(cand_rand) 


## Micro-habitat Lux and Noise Models
# Full model with random effects (1|section) (1|year) (1|stage) and (1|measurer)
model_15 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|year) (1|stage)
model_14 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|section) (1|year) (1|measurer)
model_13 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|stage) (1|measurer)
model_12 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) (1|stage) (1|measurer)
model_11 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|year)
model_10 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year), data = mobl, REML = FALSE)
# Model with (1|section) (1|stage)
model_9 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|section) (1|measurer)
model_8 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) (1|stage)
model_7 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|year) (1|measurer)
model_6 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|stage) (1|measurer)
model_5 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section)
model_4 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section), data = mobl, REML = FALSE)
# Model with (1|year)
model_3 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year), data = mobl, REML = FALSE)
# Model with (1|stage)
model_2 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|stage), data = mobl, REML = FALSE)
# Model with (1|measurer)
model_1 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|measurer), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(mass ~ scale(lux) * scale(LAeq) + sex + julian, data = mobl)

# Compare all model AICc
cand_rand <- list(
  model_15 = model_15, 
  model_14 = model_14, 
  model_13 = model_13, 
  model_12 = model_12, 
  model_11 = model_11, 
  model_10 = model_10,
  model_9 = model_9,
  model_8 = model_8,
  model_7 = model_7, 
  model_6 = model_6, 
  model_5 = model_5, 
  model_4 = model_4, 
  model_3 = model_3, 
  model_2 = model_2, 
  model_1 = model_1, 
  model_0 = model_0)
model.sel(cand_rand) 


##-- Optimized Random Effect Structure: Model 11 (1|year) + (1|stage) + (1|measurer)





#------------------------------------------------------------------------------
#---- MASS - Model Selection on Fixed Effects ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(mass ~ trt + sex + julian  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + julian + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + julian  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq)  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq)  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lmer(mass ~ scale(lux)  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + sex  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * sex  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq)  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + sex  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) * sex  + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(mass ~ 1 + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models (models with ΔAIC < 2)
top_models <- list(cand.models[[14]])  
top_models

# Top model
top_mod <- lmer(mass ~ scale(lux) + sex + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)



#---- Assumption Tests
#### Normality Test
qqPlot(residuals(top_mod)) 
plot(density(resid(top_mod)))
shapiro.test(residuals(top_mod));length(residuals(top_mod)) # p-value > 0.05 indicates normal residuals

#### Homogeneity of variance 
# Test for Heteroscedasticity (performance)
check_heteroscedasticity(top_mod) # p > 0.05 indicates homoscedasticity (constant variance)

#### Residual Diagnostics & Dispersion (DHARMa)
# Test for over/under-dispersion and residual patterns
simulation_output <- simulateResiduals(top_mod)
plot(simulation_output)
testDispersion(simulation_output)  # p-value > 0.05 indicates no significant over- or under-dispersion problems





#------------------------------------------------------------------------------
#---- MASS - MAIN  MANUSCRIPT FIGURE ----
# LMER
mobl_mass <- lmer(mass ~ scale(lux) + sex + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
summary(mobl_mass)
confint(mobl_mass)

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
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(
    #title = "Mountain Bluebird",
    x = "Light Intensity (Lux)",
    y = "Body Mass (g)"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right")

mobl_mass_lux






##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
#---- Wing Body Condition ----
# Calculate Body Condition using mass and wing
BCwing <- smi(mobl$mass, mobl$wing)
mobl$bc_wing <- BCwing


#------------------------------------------------------------------------------
#---- WING BODY CONDITION - Random Effects Model Selection ----

##  Macro-habitat Treatment Models
# Full model with random effects (1|section) (1|year) (1|stage) and (1|measurer)
model_15 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|year) (1|stage)
model_14 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|section) (1|year) (1|measurer)
model_13 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|stage) (1|measurer)
model_12 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) (1|stage) (1|measurer)
model_11 <- lmer(bc_wing ~ trt + sex + julian + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|year)
model_10 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|year), data = mobl, REML = FALSE)
# Model with (1|section) (1|stage)
model_9 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|section) (1|measurer)
model_8 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) (1|stage)
model_7 <- lmer(bc_wing ~ trt + sex + julian + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|year) (1|measurer)
model_6 <- lmer(bc_wing ~ trt + sex + julian + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|stage) (1|measurer)
model_5 <- lmer(bc_wing ~ trt + sex + julian + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section)
model_4 <- lmer(bc_wing ~ trt + sex + julian + (1|section), data = mobl, REML = FALSE)
# Model with (1|year)
model_3 <- lmer(bc_wing ~ trt + sex + julian + (1|year), data = mobl, REML = FALSE)
# Model with (1|stage)
model_2 <- lmer(bc_wing ~ trt + sex + julian + (1|stage), data = mobl, REML = FALSE)
# Model with (1|measurer)
model_1 <- lmer(bc_wing ~ trt + sex + julian + (1|measurer), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(bc_wing ~ trt + sex + julian, data = mobl)

# Compare all model AICc
cand_rand <- list(
  model_15 = model_15, 
  model_14 = model_14, 
  model_13 = model_13, 
  model_12 = model_12, 
  model_11 = model_11, 
  model_10 = model_10,
  model_9 = model_9,
  model_8 = model_8,
  model_7 = model_7, 
  model_6 = model_6, 
  model_5 = model_5, 
  model_4 = model_4, 
  model_3 = model_3, 
  model_2 = model_2, 
  model_1 = model_1, 
  model_0 = model_0)
model.sel(cand_rand) 

## Micro-habitat Lux and Noise Models
# Full model with random effects (1|section) (1|year) (1|stage) and (1|measurer)
model_15 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|year) (1|stage)
model_14 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|section) (1|year) (1|measurer)
model_13 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|stage) (1|measurer)
model_12 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) (1|stage) (1|measurer)
model_11 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section) (1|year)
model_10 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year), data = mobl, REML = FALSE)
# Model with (1|section) (1|stage)
model_9 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|section) (1|measurer)
model_8 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|year) (1|stage)
model_7 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|year) (1|measurer)
model_6 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|stage) (1|measurer)
model_5 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|stage) + (1|measurer), data = mobl, REML = FALSE)
# Model with (1|section)
model_4 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section), data = mobl, REML = FALSE)
# Model with (1|year)
model_3 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year), data = mobl, REML = FALSE)
# Model with (1|stage)
model_2 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|stage), data = mobl, REML = FALSE)
# Model with (1|measurer)
model_1 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|measurer), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian, data = mobl)

# Compare all model AICc
cand_rand <- list(
  model_15 = model_15, 
  model_14 = model_14, 
  model_13 = model_13, 
  model_12 = model_12, 
  model_11 = model_11, 
  model_10 = model_10,
  model_9 = model_9,
  model_8 = model_8,
  model_7 = model_7, 
  model_6 = model_6, 
  model_5 = model_5, 
  model_4 = model_4, 
  model_3 = model_3, 
  model_2 = model_2, 
  model_1 = model_1, 
  model_0 = model_0)
model.sel(cand_rand) 


##-- Optimized Random Effect Structure: Model 0 Linear Model



#------------------------------------------------------------------------------
#---- WING BODY CONDITION - Model Selection on Fixed Effects ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1

# trt models
cand.models[[y]] <- lm(bc_wing ~ trt + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt, data = mobl); y <- y + 1
# lux * noise
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq), data = mobl); y <- y + 1
# lux + noise
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq), data = mobl); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lm(bc_wing ~ scale(lux), data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(LAeq), data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(LAeq) * sex, data = mobl); y <- y + 1
# Null
cand.models[[y]] <- lm(bc_wing ~ 1, data = mobl); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models (models with ΔAIC < 2)
top_models <- list(cand.models[[4]], cand.models[[19]], cand.models[[2]])  
top_models

# Top model
top_mod <- lm(bc_wing ~ trt, data = mobl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)

# 2nd Top model: NULL

# 3rd Top model
top_mod3 <- lm(bc_wing ~ trt + sex, data = mobl)
summary(top_mod3)
confint(top_mod3)
confint(top_mod3, level = 0.85)


#---- Assumption Tests
#### Normality Test
qqPlot(residuals(top_mod)) 
plot(density(resid(top_mod)))
shapiro.test(residuals(top_mod));length(residuals(top_mod)) # p-value > 0.05 indicates normal residuals

#### Homogeneity of variance 
# Test for Heteroscedasticity (performance)
check_heteroscedasticity(top_mod) # p > 0.05 indicates homoscedasticity (constant variance)

#### Residual Diagnostics & Dispersion (DHARMa)
# Test for over/under-dispersion and residual patterns
simulation_output <- simulateResiduals(top_mod)
plot(simulation_output)
testDispersion(simulation_output)  # p-value > 0.05 indicates no significant over- or under-dispersion problems





#------------------------------------------------------------------------------
#---- WING BODY CONDITION - MAIN  MANUSCRIPT FIGURE ----
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
  labs(x = "Treatment Group",
       y = "Body Condition",
       #title = "Mountain Bluebird",
      color = "Treatment") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(color = "black") ,
        #legend.position = "none",
        legend.position = "right",
        axis.line = element_line(linewidth = 0.8),
        plot.title = element_text(hjust = 0.5))
mobl_BC
