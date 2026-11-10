
#####################
# Analysis of Adult WEBL (Western Bluebird) Morphology & Condition 
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

age_cols <- c("SY" = "#88CCEE", "AHY" = "#3B7A9E", "ASY" = "#103048")


#---- Subset Data to only include WEBL ----
webl <- dat %>% filter(species == "WEBL") %>% droplevels()

# Remove outlier -  Does not change results
webl <- webl %>%
  mutate(tarsus = ifelse(ID %in% c("3011-18015"), NA, tarsus)) 


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

# Calculate Body Condition using mass and wing
BCwing <- smi(webl$mass, webl$wing)
webl$bc_wing <- BCwing





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- WING Length -----

#------------------------------------------------------------------------------
#---- WING - Random Effects Model Selection ----

##  Macro-habitat Treatment Models
# Full model with random effects (1|year) (1|measurer) and (1|section)
model_7 <- lmer(wing ~ trt + sex + age + julian + (1|section) + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) and (1|measurer)
model_6 <- lmer(wing ~ trt + sex + age + julian + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|year) and (1|section)
model_5 <- lmer(wing ~ trt + sex + age + julian + (1|section) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|measurer) and (1|section)
model_4 <- lmer(wing ~ trt + sex + age + julian + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with only (1|year)
model_3 <- lmer(wing ~ trt + sex + age + julian + (1|section), data = webl, REML = FALSE)
# Model with only (1|measurer)
model_2 <- lmer(wing ~ trt + sex + age + julian + (1|year), data = webl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(wing ~ trt + sex + age + julian + (1|measurer), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(wing ~ trt + sex + age + julian, data = webl)

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
model_7 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) and (1|measurer)
model_6 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|year) and (1|section)
model_5 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|measurer) and (1|section)
model_4 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with only (1|year)
model_3 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section), data = webl, REML = FALSE)
# Model with only (1|measurer)
model_2 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year), data = webl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|measurer), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(wing ~ scale(lux) * scale(LAeq) + sex + age + julian, data = webl)

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


##-- Optimized Random Effect Structure: Model 1 (1|measurer)



#------------------------------------------------------------------------------
#---- WING - Model Selection on Fixed Effects ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(wing ~ trt + sex + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + age  + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + age  + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + age  + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# Individual 
cand.models[[y]] <- lmer(wing ~ scale(lux) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(LAeq) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(LAeq) + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(LAeq) * sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(LAeq) + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(LAeq) * age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt * sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt * age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(wing ~ 1 + (1|measurer), data = webl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models (models with ΔAIC < 2)
top_models <- list(cand.models[[18]], cand.models[[2]], cand.models[[10]])  
top_models


# Top model:  
top_mod <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + (1 | measurer), data = webl, REML = FALSE)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)

# Find means of wing by Age
adj_means <- emmeans(top_mod, ~ age)
adj_means
pairs(adj_means)

# 2nd Top model: 
top_mod_2 <- lmer( wing ~ trt + sex + age + (1 | measurer), data = webl, REML = FALSE)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)

# 3rd Top model: 
top_mod_3 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + (1 | measurer), data = webl, REML = FALSE)
summary(top_mod_3)
confint(top_mod_3)
confint(top_mod_3, level = 0.85)

# Find means of wing by Age
adj_means <- emmeans(top_mod_3, ~ age)
adj_means
pairs(adj_means)



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
#------------------------------------------------------------------------------
#---- Tarsus Length -----

#------------------------------------------------------------------------------
#---- TARSUS - Random Effects Model Selection ----

##  Macro-habitat Treatment Models
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(tarsus ~ trt + sex + (1|section) + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(tarsus ~ trt + sex + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(tarsus ~ trt + sex + (1|section) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(tarsus ~ trt + sex + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(tarsus ~ trt + sex + (1|section), data = webl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(tarsus ~ trt + sex + (1|year), data = webl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(tarsus ~ trt + sex + (1|measurer), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(tarsus ~ trt + sex, data = webl)

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
model_7 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section), data = webl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|year), data = webl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|measurer), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex, data = webl)

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


##-- Optimized Random Effect Structure: Model 4 (1|year) + (1|measurer)



#------------------------------------------------------------------------------
#---- TARSUS - Model Selection on Fixed Effects ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(tarsus ~ trt + sex + age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + sex + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * sex + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) + sex + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) * sex + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) + age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) * age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt * sex + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt * age + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(tarsus ~ 1 + (1|year) + (1|measurer), data = webl, REML = FALSE); y <- y + 1


# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models (models with ΔAIC < 2)
top_models <- list(cand.models[[25]], cand.models[[18]], cand.models[[20]]) 
top_models


# Top model: Null

# 2nd Top model:
top_mod2 <- lmer(tarsus ~ scale(LAeq) + (1 | year) + (1 | measurer), data = webl, REML = FALSE)
summary(top_mod2)
confint(top_mod2)
confint(top_mod2, level = 0.85)

# 3rd Top model:
top_mod3 <- lmer(tarsus ~ scale(LAeq) * sex + (1 | year) + (1 | measurer), data = webl, REML = FALSE)
summary(top_mod3)
confint(top_mod3)
confint(top_mod3, level = 0.85)



#---- Assumption Tests
#### Normality Test
qqPlot(residuals(top_mod3)) 
plot(density(resid(top_mod3)))
shapiro.test(residuals(top_mod3));length(residuals(top_mod3)) # p-value > 0.05 indicates normal residuals

#### Homogeneity of variance 
# Test for Heteroscedasticity (performance)
check_heteroscedasticity(top_mod3) # p > 0.05 indicates homoscedasticity (constant variance)

#### Residual Diagnostics & Dispersion (DHARMa)
# Test for over/under-dispersion and residual patterns
simulation_output <- simulateResiduals(top_mod3)
plot(simulation_output)
testDispersion(simulation_output)  # p-value > 0.05 indicates no significant over- or under-dispersion problems






#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Body Mass -----

#------------------------------------------------------------------------------
#---- MASS - Random Effects Model Selection ----

##  Macro-habitat Treatment Models
# Full model with random effects (1|section) (1|year) (1|stage) and (1|measurer)
model_15 <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|year) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|year) (1|stage)
model_14 <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with (1|section) (1|year) (1|measurer)
model_13 <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|stage) (1|measurer)
model_12 <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) (1|stage) (1|measurer)
model_11 <- lmer(mass ~ trt + sex + age + julian + (1|year) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|year)
model_10 <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|section) (1|stage)
model_9 <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE)
# Model with (1|section) (1|measurer)
model_8 <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) (1|stage)
model_7 <- lmer(mass ~ trt + sex + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with (1|year) (1|measurer)
model_6 <- lmer(mass ~ trt + sex + age + julian + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|stage) (1|measurer)
model_5 <- lmer(mass ~ trt + sex + age + julian + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section)
model_4 <- lmer(mass ~ trt + sex + age + julian + (1|section), data = webl, REML = FALSE)
# Model with (1|year)
model_3 <- lmer(mass ~ trt + sex + age + julian + (1|year), data = webl, REML = FALSE)
# Model with (1|stage)
model_2 <- lmer(mass ~ trt + sex + age + julian + (1|stage), data = webl, REML = FALSE)
# Model with (1|measurer)
model_1 <- lmer(mass ~ trt + sex + age + julian + (1|measurer), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(mass ~ trt + sex + age + julian, data = webl)

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
model_15 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|year) (1|stage)
model_14 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with (1|section) (1|year) (1|measurer)
model_13 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|stage) (1|measurer)
model_12 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) (1|stage) (1|measurer)
model_11 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|year)
model_10 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|section) (1|stage)
model_9 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE)
# Model with (1|section) (1|measurer)
model_8 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) (1|stage)
model_7 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with (1|year) (1|measurer)
model_6 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|stage) (1|measurer)
model_5 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section)
model_4 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section), data = webl, REML = FALSE)
# Model with (1|year)
model_3 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year), data = webl, REML = FALSE)
# Model with (1|stage)
model_2 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|stage), data = webl, REML = FALSE)
# Model with (1|measurer)
model_1 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|measurer), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(mass ~ scale(lux) * scale(LAeq) + sex + age + julian, data = webl)

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


##-- Optimized Random Effect Structure: Model 9 (1|section) + (1|stage)



#------------------------------------------------------------------------------
#---- MASS - Model Selection on Fixed Effects ----

# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lmer(mass ~ scale(lux) + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) * sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) * age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt * sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt * age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(mass ~ 1 + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models (models with ΔAIC < 2)
top_models <- list(cand.models[[18]], cand.models[[31]], cand.models[[33]], cand.models[[17]]
                   , cand.models[[26]] , cand.models[[21]])  
top_models


# Top model:
top_mod <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + age + (1 | section) + (1 | stage), data = webl, REML = FALSE)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)

# Find means and SD of BC by Age
adj_means <- emmeans(top_mod, ~ age)
adj_means
pairs(adj_means)

# 2nd Top model: 
top_mod_2 <- lmer(mass ~ scale(LAeq) + sex + (1 | section) + (1 | stage), data = webl, REML = FALSE)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)

# 3rd Top model: 
top_mod_3 <- lmer(mass ~ scale(LAeq) + age + (1 | section) + (1 | stage), data = webl, REML = FALSE)
summary(top_mod_3)
confint(top_mod_3)
confint(top_mod_3, level = 0.85)

# 4th Top model:
top_mod_4 <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + age + julian + (1 | section) +      (1 | stage), data = webl, REML = FALSE)
summary(top_mod_4)
confint(top_mod_4)
confint(top_mod_4, level = 0.85)

# 5th Top model: 
top_mod_5 <- lmer(mass ~ scale(lux) + sex + (1 | section) + (1 | stage), data = webl, REML = FALSE)
summary(top_mod_5)
confint(top_mod_5)
confint(top_mod_5, level = 0.85)

# 6th Top model: 
top_mod_6 <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + (1 | section) + (1 |      stage), data = webl, REML = FALSE)
summary(top_mod_6)
confint(top_mod_6)
confint(top_mod_6, level = 0.85)



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
ggplot() +
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
    legend.position = "right")
    #legend.position = "none")



# Mass across noise values using model residuals
predict_webl_noise_sex <- predict_response(webl_mass, 
                                         terms = c("LAeq [all]", "sex"), 
                                         margin = "marginalmeans", 
                                         ci.lvl = 0.85)

# WEBL Mass vs Noise
# Marginal effect (the predicted mean)
ggplot() +
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
    legend.position = "right")
    #legend.position = "none")






##------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Wing Body Condition ----
# Calculate Body Condition using mass and wing
BCwing <- smi(webl$mass, webl$wing)
webl$bc_wing <- BCwing

#------------------------------------------------------------------------------
#---- WING BODY CONDITION - Random Effects Model Selection ----

##  Macro-habitat Treatment Models
# Full model with random effects (1|section) (1|year) (1|stage) and (1|measurer)
model_15 <- lmer(bc_wing ~ trt + sex + age + julian + (1|section) + (1|year) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|year) (1|stage)
model_14 <- lmer(bc_wing ~ trt + sex + age + julian + (1|section) + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with (1|section) (1|year) (1|measurer)
model_13 <- lmer(bc_wing ~ trt + sex + age + julian + (1|section) + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|stage) (1|measurer)
model_12 <- lmer(bc_wing ~ trt + sex + age + julian + (1|section) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) (1|stage) (1|measurer)
model_11 <- lmer(bc_wing ~ trt + sex + age + julian + (1|year) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|year)
model_10 <- lmer(bc_wing ~ trt + sex + age + julian + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|section) (1|stage)
model_9 <- lmer(bc_wing ~ trt + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE)
# Model with (1|section) (1|measurer)
model_8 <- lmer(bc_wing ~ trt + sex + age + julian + (1|section) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) (1|stage)
model_7 <- lmer(bc_wing ~ trt + sex + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with (1|year) (1|measurer)
model_6 <- lmer(bc_wing ~ trt + sex + age + julian + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|stage) (1|measurer)
model_5 <- lmer(bc_wing ~ trt + sex + age + julian + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section)
model_4 <- lmer(bc_wing ~ trt + sex + age + julian + (1|section), data = webl, REML = FALSE)
# Model with (1|year)
model_3 <- lmer(bc_wing ~ trt + sex + age + julian + (1|year), data = webl, REML = FALSE)
# Model with (1|stage)
model_2 <- lmer(bc_wing ~ trt + sex + age + julian + (1|stage), data = webl, REML = FALSE)
# Model with (1|measurer)
model_1 <- lmer(bc_wing ~ trt + sex + age + julian + (1|measurer), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(bc_wing ~ trt + sex + age + julian, data = webl)

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
model_15 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|year) (1|stage)
model_14 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with (1|section) (1|year) (1|measurer)
model_13 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|stage) (1|measurer)
model_12 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) (1|stage) (1|measurer)
model_11 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year) + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) (1|year)
model_10 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|section) (1|stage)
model_9 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE)
# Model with (1|section) (1|measurer)
model_8 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) (1|stage)
model_7 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with (1|year) (1|measurer)
model_6 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|stage) (1|measurer)
model_5 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|stage) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section)
model_4 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section), data = webl, REML = FALSE)
# Model with (1|year)
model_3 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year), data = webl, REML = FALSE)
# Model with (1|stage)
model_2 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|stage), data = webl, REML = FALSE)
# Model with (1|measurer)
model_1 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|measurer), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian, data = webl)

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


##-- Optimized Random Effect Structure: Model 9 (1|section) + (1|stage)



#------------------------------------------------------------------------------
#---- WING BODY CONDITION - Model Selection on Fixed Effects ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1

# trt models
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) * sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) * sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(LAeq) + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(LAeq) + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(LAeq) * sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(LAeq) + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(LAeq) * age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt * sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt * age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(bc_wing ~ 1 + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models (models with ΔAIC < 2)
top_models <- list(cand.models[[10]], cand.models[[19]])  
top_models

# Top model: 
top_mod <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + (1 | section) + (1 | stage), data = webl, REML = FALSE)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)

# Find means and SD of BC by Age
adj_means <- emmeans(top_mod, ~ age)
adj_means
pairs(adj_means)

# 2nd Top model: 
top_mod_2 <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + (1 | section) +  (1 | stage), data = webl, REML = FALSE)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)



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
ggplot() +
  geom_ribbon(data = predict_webl_lux_sex, 
              aes(x = x, ymin = conf.low, ymax = conf.high, group = group), 
              fill = "#FFE082", alpha = 0.3) +
  geom_line(data = predict_webl_lux_sex, 
            aes(x = x, y = predicted, group = group, linetype = group), 
            color = "#FFE082", linewidth = 1.5) +
  geom_point(data = webl, 
             aes(x = lux, y = bc_wing, shape = sex), 
             color = "black", size = 4, alpha = 0.8) +
  annotate("text", x = 1.5, y = 33, label = "*", size = 10, color = "grey25") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(x = "Light Intensity (Lux)",
       y = "Body Condition"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right")
    #legend.position = "none")



# Body Condition across noise values using model residuals
predict_webl_noise_sex <- predict_response(webl_BC_int, 
                                           terms = c("LAeq [all]", "sex"), 
                                           margin = "marginalmeans", 
                                           ci.lvl = 0.85)

# WEBL BC vs Noise
# Marginal effect (the predicted mean)
ggplot() +
  geom_ribbon(data = predict_webl_noise_sex, 
              aes(x = x, ymin = conf.low, ymax = conf.high, group = group), 
              fill = "firebrick", alpha = 0.2) +
  geom_line(data = predict_webl_noise_sex, 
            aes(x = x, y = predicted, group = group, linetype = group), 
            color = "firebrick", linewidth = 1.5) +
  geom_point(data = webl, 
             aes(x = LAeq, y = bc_wing, shape = sex), 
             color = "black", size = 4, alpha = 0.8) +
  annotate("text", x = 52, y = 33, label = "†", size = 5, color = "grey25") +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  labs(x = "Noise Level (LAeq)",
       y = "Body Condition"  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right")
    #legend.position = "none")




# WEBL BC Interaction between Lux & Noise
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

ggplot() +
  geom_ribbon(data = plot_df_webl, 
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = lux_label), 
              alpha = 0.3, show.legend = FALSE) + 
  geom_line(data = plot_df_webl, 
            aes(x = x, y = predicted, color = lux_label), 
            linewidth = 1.5, alpha = 0.8) +
  annotate("text", x = 52, y = 35, label = "†", size = 5, color = "grey25") +
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
    legend.position = "right",
    #legend.position = "none",
    legend.title = element_text(size = 10, face = "bold") ) +
  guides(
    fill = guide_colorbar(order = 2),
    shape = guide_legend(order = 3) )



# WEBL AGE BC FIGURE
webl_age_bc <- ggplot(webl, aes(x = age, y = bc_wing, color = age)) +
  geom_point(aes(shape = sex), position = position_jitterdodge(jitter.width = 0.3, 
                                                               dodge.width = 0.5), size = 4, alpha = 0.8) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               aes(group = sex), position = position_dodge(width = 0.5),
               width = 0.15, linewidth = 1, color = "black", alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", 
               aes(group = sex, shape = sex), position = position_dodge(width = 0.5),
               size = 4, color = "black", stroke = 1.5, show.legend = FALSE) +
  annotate("text", x = 3, y = 32.9, label = "*", size = 10, color = "grey25") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16)) + 
  scale_color_manual(values = age_cols) +
  scale_y_continuous(breaks = seq(22, 32, by = 2)) +
  coord_cartesian(ylim = c(21.3, 33)) +
  labs(x = "Age Group",
       y = "Body Condition",
       color = "Age", shape = "Sex") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(color = "black") ,
        legend.position = "right",
        axis.line = element_line(linewidth = 0.8),
        plot.title = element_text(hjust = 0.5))
webl_age_bc
