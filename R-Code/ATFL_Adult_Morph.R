#####################
# 2018-2025 Adult ATFL BC - Model Selection

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(lubridate)    # to transform dates
library(AICcmodavg)   # model selection

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
    measurer = as.factor(measurer)
  )


# Subset the Data to only include ATFL
atfl <- dat %>% filter(species == "ATFL") %>% droplevels()

# Remove mass outlier
# Does not change results, helps with residual normality
atfl <- atfl %>%
  mutate(mass = ifelse(ID %in% c("2921-17874"), NA, mass))


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
BCwing <- smi(atfl$mass, atfl$wing)
atfl$bc_wing <- BCwing


# Print the # of individuals in which we have a lux value of 
sum(!is.na(atfl$lux))
# 48 of the 61 individuals have a lux value



#------------------------------------------------------------------------------
#---- Exploratory Graphs and Sample Sizes -----
# Look at Lux Values within each trt 
ggplot(atfl, aes(x= trt, y = lux, fill = trt)) +
  geom_jitter(width = 0.25, size = 2, shape = 21) +
  stat_summary(colour = "black") +
  labs(title = "ATFL Lux Values per Trt Group", x = "Trt", y = "Lux") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# Look at Noise values within each trt 
ggplot(atfl, aes(x= trt, y = LAeq, fill = trt)) +
  geom_jitter(width = 0.25, size = 2, shape = 21) +
  stat_summary(colour = "black") +
  labs(title = "ATFL LAeq Values per Trt Group", x = "Trt", y = "LAeq")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# Individuals across lux values
ggplot(atfl, aes(x = 1, y = lux)) +
  geom_jitter(width = 0.2, height = 0, color = "yellow4", alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Individuals across LAeq values
ggplot(atfl, aes(x = 1, y = LAeq)) +
  geom_jitter(width = 0.2, height = 0, color = "darkgreen", alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- WING LENGTH -----

# Remove individuals that have NA for mass
wing <- atfl %>% filter(!is.na(wing))
# 61 individuals from 24-25
wing_lux <- wing %>% filter(!is.na(lux))
# 48 individuals have lux of those that have mass

# Sample sizes of individuals per Trt group and that have mass 
ggplot(wing, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
wing %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
wing %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()

# Sample sizes of individuals per Trt group and that have BOTH mass and Tarsus AND a lux value
ggplot(wing_lux, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
wing_lux %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
wing_lux %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()



#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for WING Models ####
# Wing across sexes
ggplot(atfl, aes(x = sex, y = wing, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Wing Across Sexes",x = "Sex", y = "Wing Length")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# test for Mass across sexes
sex <- lm(wing ~ sex, data = atfl)
summary(sex)
anova(sex) # SIG M have longer wings

# Wing across ages
ggplot(atfl, aes(x = age, y = wing, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Wing Across Ages",x = "Age", y = "Wing")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Ages
age <- lm(wing ~ age, data = atfl)
summary(age)
anova(age) # NS

# Wing across years
ggplot(atfl, aes(x = year, y = wing, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Wing Across Years", x = "Year", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Years
year <- lm(wing ~ year, data = atfl)
summary(year)
anova(year) # NS

# Wing across Section
ggplot(atfl, aes(x = section, y = wing, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Wing Across Sections", x = "Section", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Sections
section <- lm(wing ~ section, data = atfl)
summary(section)
anova(section) # NS

# Wing across Julian Date/Season
ggplot(atfl, aes(x = julian, y = wing)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "ATFL Wing Across Julian Date",
       x = "Julian Date", y = "Wing Length") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Model to see tarsus across the season
julian <- lm(wing ~ julian, data = atfl)
summary(julian) # NS 

# Wing across Reproductive Stages
ggplot(atfl, aes(x = stage, y = wing, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "ATFL Wing Across Reproductive Stages", x = "Stage", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Reproductive Stages
stage <- lm(wing ~ stage, data = atfl)
summary(stage)
anova(stage) # NS

# Wing across Measure ID 
ggplot(atfl, aes(x = measurer, y = wing, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Wing Across Measurers", x = "Individual Measurer", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Individual Measureres
measurer <- lm(wing ~ measurer, data = atfl)
summary(measurer)
anova(measurer) 



#------------------------------------------------------------------------------
#---- WING - MAIN FIGURES ----
# Wing across treatments 
ggplot(atfl, aes(x = trt, y = wing)) +
  geom_jitter(width = 0.2, aes(color = trt), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "ATFL Wing", x = "Treatment Group", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing across treatments BY SEX
ggplot(atfl, aes(x = trt, y = wing, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "ATFL Wing by Sex", x = "Treatment Group", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing across Lux light values 
ggplot(atfl, aes(x = lux, y = wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7) + 
  labs(title = "ATFL Wing across Light", x = "Lux", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing across Lux light values  BY SEX
ggplot(atfl, aes(x = lux, y = wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "ATFL Wing across Light", x = "Lux", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing across LAeq noise values 
ggplot(atfl, aes(x = LAeq, y = wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7) + 
  labs(title = "ATFL Wing across Noise", x = "LAeq", y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing across LAeq noise values BY SEX 
ggplot(atfl, aes(x = LAeq, y = wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "ATFL Wing across Noise", x = "LAeq", y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass between Light and non light sites
ggplot(atfl, aes(x = light, y = wing)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "ATFL Wing between Light Sites", x = "Explosed to Light?", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing between Light and non light sites BY SEX
ggplot(atfl, aes(x = light, y = wing, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "ATFL Wing by Light Site and Sex", x = "Exposed to Light?",  y = "Wing Length",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


#------------------------------------------------------------------------------
#---- WING - Random Effects Model Selection ----
##  Treatment Model 
# Full model with random effects (1|year) (1|measurer) and (1|section)
model_7 <- lmer(wing ~ trt + sex + julian + (1|section) + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|year) and (1|measurer)
model_6 <- lmer(wing ~ trt + sex + julian + (1|section) + (1|year), data = atfl, REML = FALSE)
# Model with (1|year) and (1|section)
model_5 <- lmer(wing ~ trt + sex + julian + (1|section) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|measurer) and (1|section)
model_4 <- lmer(wing ~ trt + sex + julian + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with only (1|year)
model_3 <- lmer(wing ~ trt + sex + julian + (1|section), data = atfl, REML = FALSE)
# Model with only (1|measurer)
model_2 <- lmer(wing ~ trt + sex + julian + (1|year), data = atfl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(wing ~ trt + sex + julian + (1|measurer), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(wing ~ trt + sex + julian, data = atfl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best Models: 1, 5 (4 is two away)


## Lux and Noise Model 
# Full model with random effects (1|year) (1|measurer) and (1|section)
model_7 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|year) and (1|measurer)
model_6 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year), data = atfl, REML = FALSE)
# Model with (1|year) and (1|section)
model_5 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|measurer) and (1|section)
model_4 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with only (1|year)
model_3 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section), data = atfl, REML = FALSE)
# Model with only (1|measurer)
model_2 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year), data = atfl, REML = FALSE)
# Model with only (1|section)
model_1 <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|measurer), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian, data = atfl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best Models: 1 (4, 5 two away)


# Moving forward - LMER (1|measurer) 


#------------------------------------------------------------------------------
#---- WING - Model Selection across all model types ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(wing ~ trt + sex + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lmer(wing ~ scale(lux) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(LAeq) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(LAeq) * sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(LAeq) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(wing ~ 1 + (1|measurer), data = atfl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[51]])  
top_models

# Top model:  wing ~ scale(lux) + sex + (1 | measurer)
top_mod <- lmer(wing ~ scale(lux) + sex + (1 | measurer), data = atfl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
# No effect of lux
# SIG effect of Sex




#### Normality Test
library(car)
shapiro.test(residuals(top_mod));length(residuals(top_mod)) # residuals are normal
qqPlot(residuals(top_mod)) 
plot(density(resid(top_mod)))
# residuals are normal if p is > 0.05.


#### Homogeneity of variance 
library(DHARMa)
simulation_output <- simulateResiduals(top_mod)
plot(simulation_output)
testDispersion(simulation_output)  # This will test for over/underdispersion
# p-value > 0.05 means residuals are homogeneous and that there are no dispersion problems

# Manual BP (Breusch-Pagan) Test  
resids_sq <- (residuals(top_mod, type = "pearson"))^2 # Extract residuals and fitted values
fitted_vals <- fitted(top_mod)
bp_manual <- lm(resids_sq ~ fitted_vals) # Run a simple linear model to check for a relationship
summary(bp_manual) # variance is constant
# p-value > 0.05 means variance is constant (homoscedasitity)



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Tarsus LENGTH -----
# Remove individuals that have NA for tarsus
tarsus <- atfl %>% filter(!is.na(tarsus))
# 61 individuals from 24-25
tarsus_lux <- tarsus %>% filter(!is.na(lux))
# 48 individuals have lux of those that have tarsus

# Sample sizes of individuals per Trt group and that have mass 
ggplot(tarsus, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
tarsus %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
tarsus %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()

# Sample sizes of individuals per Trt group and that have BOTH mass and Tarsus AND a lux value
ggplot(tarsus_lux, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
tarsus_lux %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
tarsus_lux %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()




#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for Tarsus Models ####
# Tarsus across sexes
ggplot(atfl, aes(x = sex, y = tarsus, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Tarsus Across Sexes",x = "Sex", y = "Tarsus")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# test for Mass across sexes
sex <- lm(tarsus ~ sex, data = atfl)
summary(sex)
anova(sex)  

# Tarsus across ages
ggplot(atfl, aes(x = age, y = tarsus, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Tarsus Across Ages",x = "Age", y = "Tarsus")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Tarsus across Ages
age <- lm(tarsus ~ age, data = atfl)
summary(age)
anova(age) # NS

# Tarsus across years
ggplot(atfl, aes(x = year, y = tarsus, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Tarsus Across Years", x = "Year", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Tarsus across Years
year <- lm(tarsus ~ year, data = atfl) # Weak effect of year
summary(year)
anova(year) 

# Tarsus across Section
ggplot(atfl, aes(x = section, y = tarsus, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Tarsus Across Sections", x = "Section", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Tarsus across Sections
section <- lm(tarsus ~ section, data = atfl)
summary(section)
anova(section) # NS

# Tarsus across Julian Date/Season
ggplot(atfl, aes(x = julian, y = tarsus)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "ATFL Tarsus Across Julian Date",
       x = "Julian Date", y = "Tarsus") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Model to see tarsus across the season
julian <- lm(tarsus ~ julian, data = atfl)
summary(julian) # NS effect of julian

# Tarsus across Reproductive Stages
ggplot(atfl, aes(x = stage, y = tarsus, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "ATFL Tarsus Across Reproductive Stages", x = "Stage", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Tarsus across Reproductive Stages
stage <- lm(tarsus ~ stage, data = atfl)
summary(stage)
anova(stage) # No effect of stage

# Tarsus across Measure ID 
ggplot(atfl, aes(x = measurer, y = tarsus, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Tarsus Across Measurers", x = "Individual Measurer", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing Body Condition across Individual Measureres
measurer <- lm(tarsus ~ measurer, data = atfl)
summary(measurer)
anova(measurer) # SIG effect of measurer




#------------------------------------------------------------------------------
#---- TARSUS - MAIN FIGURES ----
# Tarsus across treatments 
ggplot(atfl, aes(x = trt, y = tarsus)) +
  geom_jitter(width = 0.2, aes(color = trt), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "ATFL Tarsus", x = "Treatment Group", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")


# Tarsus across treatments BY SEX
ggplot(atfl, aes(x = trt, y = tarsus, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "ATFL Tarsus by Sex", x = "Treatment Group", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus across Lux light values 
ggplot(atfl, aes(x = lux, y = tarsus)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7) + 
  labs(title = "ATFL Tarsus across Lux", x = "Lux", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus across Lux light values BY SEX
ggplot(atfl, aes(x = lux, y = tarsus, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "ATFL Tarsus across Lux", x = "Lux", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus across LAeq noise values 
ggplot(atfl, aes(x = LAeq, y = tarsus)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7) + 
  labs(title = "ATFL Tarsus across Noise", x = "LAeq", y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus across LAeq noise values BY SEX 
ggplot(atfl, aes(x = LAeq, y = tarsus, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "ATFL Tarsus across Noise", x = "LAeq", y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus Light and non light sites
ggplot(atfl, aes(x = light, y = tarsus)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "ATFL Tarsus between Light Sites", x = "Explosed to Light?", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus between Light and non light sites BY SEX
ggplot(atfl, aes(x = light, y = tarsus, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "ATFL Tarsus by Light Site and Sex",
       x = "Exposed to Light?", 
       y = "Tarsus Length",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))




#------------------------------------------------------------------------------
#---- TARSUS - Random Effects Model Selection ----
##  Treatment Model 
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(tarsus ~ trt + sex + (1|section) + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(tarsus ~ trt + sex + (1|section) + (1|year), data = atfl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(tarsus ~ trt + sex + (1|section) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(tarsus ~ trt + sex + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(tarsus ~ trt + sex + (1|section), data = atfl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(tarsus ~ trt + sex + (1|year), data = atfl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(tarsus ~ trt + sex + (1|measurer), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(tarsus ~ trt + sex , data = atfl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 1 then 4


## Lux and Noise Model 
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|year), data = atfl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|section), data = atfl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|year), data = atfl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|measurer), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex, data = atfl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 0, 1


# Moving forward - LMER (1|measurer)




#------------------------------------------------------------------------------
#---- TARSUS - Model Selection across all model types ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(tarsus ~ trt + sex + section + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + section + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + section + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + section + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + section + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + section + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# Individual lux, noise, light
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) * sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(tarsus ~ 1 + (1|measurer), data = atfl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[25]], cand.models[[24]], cand.models[[23]], cand.models[[16]], cand.models[[15]])  
top_models

# Top model: tarsus ~ scale(lux) + (1 | measurer) 
top_mod <- lmer(tarsus ~ scale(lux) + (1 | measurer), data = atfl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
# No effect of lux

# 2nd Top model:  tarsus ~ scale(lux) + scale(LAeq) + (1 | measurer)
top_mod_2 <- lmer(tarsus ~ scale(lux) + scale(LAeq) + (1 | measurer), data = atfl)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)
# No effect of lux or noise

# 3rd Top model: tarsus ~ scale(lux) * scale(LAeq) + year + (1 | measurer)
top_mod_3 <- lmer(tarsus ~ scale(lux) + scale(LAeq) + year + (1 | measurer), data = atfl)
summary(top_mod_3)
confint(top_mod_3)
confint(top_mod_3, level = 0.85)
# No effect of lux or noise
# effect of year at 85% CI

# 4th Top model: tarsus ~ scale(lux) * scale(LAeq) + (1 | measurer)
top_mod_4 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + (1 | measurer), data = atfl)
summary(top_mod_4)
confint(top_mod_4)
confint(top_mod_4, level = 0.85)
# No effect of lux, noise, or int

# 5th Top model: tarsus ~ scale(lux) * scale(LAeq) + year + (1 | measurer)
top_mod_5 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + year + (1 | measurer), data = atfl)
summary(top_mod_5)
confint(top_mod_5)
confint(top_mod_5, level = 0.85)



#### Normality Test
library(car)
shapiro.test(residuals(top_mod_5));length(residuals(top_mod_5)) # residuals are normal
qqPlot(residuals(top_mod_5)) 
plot(density(resid(top_mod_5)))
# residuals are normal if p is > 0.05.


#### Homogeneity of variance 
# Performance Package to test for Heteroscedasticity
library(performance)
check_heteroscedasticity(top_mod_5) # model variance is homoscedastic

library(DHARMa)
simulation_output <- simulateResiduals(top_mod_5)
plot(simulation_output)
testQuantiles(simulation_output)  # This will test for homoscedasticity
testDispersion(simulation_output)  # This will test for over/underdispersion
# p-value > 0.05 means residuals are homogeneous and that there are no dispersion problems

# Manual BP (Breusch-Pagan) Test  
resids_sq <- (residuals(top_mod_5, type = "pearson"))^2 # Extract residuals and fitted values
fitted_vals <- fitted(top_mod_5)
bp_manual <- lm(resids_sq ~ fitted_vals) # Run a simple linear model to check for a relationship
summary(bp_manual) # variance is constant
# p-value > 0.05 means variance is constant (homoscedasitity)



#### 2 categories Tarsus across LAeq noise AND Lux values
atfl <- atfl %>%
  mutate(lux_binary = if_else(lux == 0, "No Lux", "Lux"),
         lux_binary = factor(lux_binary, levels = c("No Lux", "Lux")))

atfl |> 
  filter(!is.na(lux_binary)) |> 
  ggplot(aes(x = lux_binary, y = lux)) +
  geom_point()
atfl |> 
  filter(!is.na(lux_binary)) |> 
  ggplot(aes(x = LAeq, y = tarsus, group = lux_binary, color = lux_binary)) +
  geom_point(size = 3) +  
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Lux" = "red", "No Lux" = "grey50")) +
  labs(title = "ATFL Tarsus across Noise & Lux", x = "LAeq", y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



### 3 categories - Tarsus across LAeq noise AND Lux values
atfl <- atfl %>%
  mutate(pos_lux_mean = mean(lux[lux > 0], na.rm = TRUE)) %>% # Calculate the mean of ONLY the positive values
  mutate(lux_cat3 = case_when(  # Assign categories based on that threshold
    lux == 0 ~ "No Lux",
    lux > 0 & lux <= pos_lux_mean ~ "Low Lux",
    lux > pos_lux_mean ~ "High Lux",
    TRUE ~ NA_character_ # Handles any existing NAs in the lux data
  )) %>%
  mutate(lux_cat3 = factor(lux_cat3, levels = c("No Lux", "Low Lux", "High Lux")))

atfl |> 
  filter(!is.na(lux_cat3)) |> 
  ggplot(aes(x = lux_cat3, y = lux)) +
  geom_point()


atfl %>%
  filter(!is.na(lux_cat3)) %>%
  ggplot(aes(x = LAeq, y = tarsus, color = lux_cat3)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  scale_color_manual(values = c("No Lux" = "black", 
                                "Low Lux" = "blue", 
                                "High Lux" = "red")) +
  labs(
    title = "ATFL Tarsus: Interaction of Noise and Light Intensity",
    x = "Noise Level (LAeq)",
    y = "Tarsus Length",
    color = "Lux Category"
  ) +
  theme_minimal()





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Body Mass -----

# Remove individuals that have NA for mass
mass <- atfl %>% filter(!is.na(mass))
# 59 individuals from 24-25
mass_lux <- mass %>% filter(!is.na(lux))
# 47 individuals have lux of those that have mass

# Sample sizes of individuals per Trt group and that have mass 
ggplot(mass, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
mass %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
mass %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()

# Sample sizes of individuals per Trt group and that have mass AND a lux value
ggplot(mass_lux, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
mass_lux %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
mass_lux %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()




#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for Mass Models ####
# Mass across sexes
ggplot(atfl, aes(x = sex, y = mass, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Mass Across Sexes",x = "Sex", y = "Mass (g)")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# test for Mass across sexes
sex <- lm(mass ~ sex, data = atfl)
summary(sex)
anova(sex) # NS

# Mass across ages
ggplot(atfl, aes(x = age, y = mass, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Mass Across Ages",x = "Age", y = "Mass (g)")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Ages
age <- lm(mass ~ age, data = atfl)
summary(age)
anova(age) # NS

# Mass across years
ggplot(atfl, aes(x = year, y = mass, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Mass Across Years", x = "Year", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Years
year <- lm(mass ~ year, data = atfl)
summary(year)
anova(year) # NS

# Mass across Section
ggplot(atfl, aes(x = section, y = mass, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Mass Across Sections", x = "Section", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Sections
section <- lm(mass ~ section, data = atfl)
summary(section)
anova(section) # NS

# Mass across Julian Date/Season
ggplot(atfl, aes(x = julian, y = mass)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "ATFL Mass Across Julian Date",
       x = "Julian Date", y = "Mass (g)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Model to see mass across the season
julian <- lm(mass ~ julian, data = atfl)
summary(julian) # Julian date SIG affects mass
anova(julian)

# Mass across Reproductive Stages
ggplot(atfl, aes(x = stage, y = mass, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "ATFL Mass Across Reproductive Stages", x = "Stage", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Reproductive Stages
stage <- lm(mass ~ stage, data = atfl)
summary(stage)
anova(stage) # SIG effect of stage

# Reproductive stages across Julian Date
ggplot(atfl, aes(x = julian, y = stage, fill = stage)) +
  geom_boxplot(alpha = 0.5, width = 0.3) +
  geom_jitter(width = 0.1, alpha = 1) +
  labs(title = "ATFL Reproductive Stage by Julian Day",
       x = "Julian Date",
       y = "Stage") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Correlation test between julain date and stage
corr_test <- aov(julian ~ stage, data = atfl)
summary(corr_test) # Julian date and stage are highly correlated

# Mass across Measure ID
ggplot(atfl, aes(x = measurer, y = mass, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL Mass Across Measurers", x = "Individual Measurer", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Individual Measureres
measurer <- lm(mass ~ measurer, data = atfl)
summary(measurer)
anova(measurer) # NS


#------------------------------------------------------------------------------
#---- MASS - MAIN FIGURES ----
# Mass across treatments 
ggplot(atfl, aes(x = trt, y = mass)) +
  geom_jitter(width = 0.2, aes(color = trt), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "ATFL Mass", x = "Treatment Group", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Mass across treatments BY SEX
ggplot(atfl, aes(x = trt, y = mass, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "ATFL Mass by Sex", x = "Treatment Group", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across Lux light values 
ggplot(atfl, aes(x = lux, y = mass)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, size = 2, alpha = 1) + 
  labs(title = "ATFL Mass across Lux", x = "Lux", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across Lux light values BY SEX
ggplot(atfl, aes(x = lux, y = mass, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "ATFL Mass across Lux", x = "Lux", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across LAeq noise values 
ggplot(atfl, aes(x = LAeq, y = mass)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, size = 2, alpha = 1) + 
  labs(title = "ATFL Mass across Noise", x = "LAeq", y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across LAeq noise values BY SEX 
ggplot(atfl, aes(x = LAeq, y = mass, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "ATFL Mass across Noise", x = "LAeq", y = "Mass(g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass between Light and non light sites
ggplot(atfl, aes(x = light, y = mass)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "ATFL Mass between Light Sites", x = "Explosed to Light?", 
       y = "Mass") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass between Light and non light sites BY SEX
ggplot(atfl, aes(x = light, y = mass, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "ATFL Mass by Light Site and Sex",
       x = "Exposed to Light?", 
       y = "Mass (g)",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across LAeq noise values in regard to Light or No Light
ggplot(atfl, aes(x = LAeq, y = mass, group = light)) +
  geom_smooth(method = "lm", aes(color = light), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = light), size = 2, alpha = 0.7) + 
  labs(title = "ATFL Mass across Noise and Light (Y/N)", x = "LAeq", y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


#------------------------------------------------------------------------------
#---- MASS - Random Effects Model Selection ----
##  Treatment Model 
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|year) + (1|stage), data = atfl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|year), data = atfl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|stage), data = atfl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(mass ~ trt + sex + julian + (1|year) + (1|stage), data = atfl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(mass ~ trt + sex + julian + (1|section), data = atfl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(mass ~ trt + sex + julian + (1|year), data = atfl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(mass ~ trt + sex + julian + (1|stage), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(mass ~ trt + sex + julian, data = atfl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 0, 1


## Lux and Noise Model 
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|stage), data = atfl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year), data = atfl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage), data = atfl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage), data = atfl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section), data = atfl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year), data = atfl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|stage), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(mass ~ scale(lux) * scale(LAeq) + sex + julian, data = atfl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best model: 1, 0


# Moving forward - LMER (1|stage)

#------------------------------------------------------------------------------
#---- MASS - Model Selection across all model types ----

# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(mass ~ trt + sex + julian + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + julian + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + julian + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + julian + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + julian + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + julian + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + julian + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + julian + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + (1|stage), data = atfl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + julian + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + julian + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + julian + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + julian + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + (1|stage), data = atfl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + year + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + year + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + section + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + (1|stage), data = atfl, REML = FALSE); y <- y + 1
# Individual lux, noise, light
cand.models[[y]] <- lmer(mass ~ scale(lux) + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + sex + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * sex + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + sex + (1|stage), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) * sex + (1|stage), data = atfl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(mass ~ 1 + (1|stage), data = atfl, REML = FALSE); y <- y + 1

# Generate names for models
modnames1 <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table1 <- aictab(cand.models, modnames = modnames1) 
print(aic_table1)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[49]], cand.models[[48]]) 
top_models

# Top model: mass ~ scale(lux) + (1 | stage)
top_mod <- lmer(mass ~ scale(lux) + (1 | stage), data = atfl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
# lux has neg effect at 85% CI

# *** 2nd Top model: mass ~ scale(lux) + scale(LAeq) + (1 | stage)
top_mod_2 <- lmer(mass ~ scale(lux) + scale(LAeq) + (1 | stage), data = atfl)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)
# Lux has a negative effect at 85% CI
# No effect of noise



#### Normality Test
library(car)
shapiro.test(residuals(top_mod));length(residuals(top_mod))
qqPlot(residuals(top_mod)) 
plot(density(resid(top_mod)))
# residuals are normal if p is > 0.05.
# residuals are normally distributed. 


#### Homogeneity of variance 
# Performance Package to test for Heteroscedasticity
library(performance)
check_heteroscedasticity(top_mod) # model variance is homoscedastic

library(DHARMa)
simulation_output <- simulateResiduals(top_mod)
plot(simulation_output)
testQuantiles(simulation_output)  # This will test for homoscedasticity
testDispersion(simulation_output)  # This will test for over/underdispersion
# p-value > 0.05 means residuals are homogeneous and that there are no dispersion problems

# Manual BP (Breusch-Pagan) Test  
resids_sq <- (residuals(top_mod, type = "pearson"))^2 # Extract residuals and fitted values
fitted_vals <- fitted(top_mod)
bp_manual <- lm(resids_sq ~ fitted_vals) # Run a simple linear model to check for a relationship
summary(bp_manual) # variance is constant
# p-value > 0.05 means variance is constant (homoscedasitity)




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Wing Body Condition ----
# Calculate Body Condition using mass and wing
BCwing <- smi(atfl$mass, atfl$wing)
atfl$bc_wing <- BCwing

# Mass by Wing
ggplot(atfl, aes(x = mass, y = wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "ATFL Mass by Wing", x = "Mass", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")


# Mass by Wing by Sex
ggplot(atfl, aes(x = mass, y = wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "ATFL Mass by Wing", x = "Mass", y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))






#------------------------------------------------------------------------------
#---- Exploratory Graphs and Sample Sizes -----
# Remove individuals that have NA for mass or wing
bc_w <- atfl %>% filter(!is.na(mass), !is.na(wing))
# 60 individuals
bc_w_lux <- bc_w %>% filter(!is.na(lux))
# 48 individuals have lux of those that have mass and wing

# Sample sizes of individuals per Trt group that have BOTH mass and Wing
ggplot(bc_w, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
bc_w %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
bc_w %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()

# Sample sizes of individuals per Trt group and that have BOTH mass and Wing AND a lux value
ggplot(bc_w_lux, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
bc_w_lux %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
bc_w_lux %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()






#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for Wing Body Condition Models ####
# Body condition across sexes
ggplot(atfl, aes(x = sex, y = bc_wing, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL W-BC Across Sexes",x = "Sex", y = "Body Condition")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model for Body Condition across sexes
sex <- lm(bc_wing ~ sex, data = atfl)
summary(sex)
anova(sex) # SIG

# Wing BC across ages
ggplot(atfl, aes(x = age, y = bc_wing, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL W-BC Across Ages",x = "Age", y = "Body Condition")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing BC across Ages
age <- lm(bc_wing ~ age, data = atfl)
summary(age)
anova(age) # NS

# Body Condition across years
ggplot(atfl, aes(x = year, y = bc_wing, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL W-BC Across Years", x = "Year", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Body Condition across Years
year <- lm(bc_wing ~ year, data = atfl)
summary(year)
anova(year) # No effect of year

# Body Condition across Section
ggplot(atfl, aes(x = section, y = bc_wing, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL W-BC Across Sections", x = "Section", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Sections
section <- lm(bc_wing ~ section, data = atfl)
summary(section) # Pink section trending to reduce BC
anova(section) # No effect of section

# Body Condition across Julian Date/Season
ggplot(atfl, aes(x = julian, y = bc_wing)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "ATFL W-BC Across Julian Date",
       x = "Julian Date", y = "Body Condition") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Model to see body condition across the season
julian <- lm(bc_wing ~ julian, data = atfl)
summary(julian) 
anova(julian) # SIG effect of Julian

# Wing Body Condition across Reproductive Stages
ggplot(atfl, aes(x = stage, y = bc_wing, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "ATFL W-BC Across Reproductive Stages", x = "Stage", 
       y = "Wing Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing Body Condition across Reproductive Stages
stage <- lm(bc_wing ~ stage, data = atfl)
summary(stage)
anova(stage) # SIG effect of stage

# Post Hoc Tests
post_hoc <- emmeans(stage, ~ stage)
pairs(post_hoc, simple = "each")
plot(post_hoc) # chx stage is sig lower

# Check Correlation between Julian date and Reproductive stage
ggplot(atfl, aes(x = stage, y = julian, fill = stage)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  geom_line(aes(group = nest)) +
  labs(title = "Julian Day by Reproductive Stage",
       x = "Reproductive Stage",
       y = "Julian Day") +
  theme_minimal()

# Correlation test between julain date and stage
corr_test <- aov(julian ~ stage, data = atfl)
summary(corr_test) # Julian date and stage are highly correlated

# Wing Body Condition across Measure ID
ggplot(atfl, aes(x = measurer, y = bc_wing, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "ATFL W-BC Across Measurers", x = "Individual Measurer", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model ofWing  Body Condition across Individual Measureres
measurer <- lm(bc_wing ~ measurer, data = atfl)
summary(measurer)
anova(measurer) # No effect of individual measurer



#------------------------------------------------------------------------------
#---- WING BODY CONDITION - MAIN FIGURES ----
# Wing Body condition across treatments 
ggplot(atfl, aes(x = trt, y = bc_wing)) +
  geom_jitter(width = 0.2, aes(color = trt), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "ATFL W-BC", x = "Treatment Group", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing Body condition across treatments BY SEX
ggplot(atfl, aes(x = trt, y = bc_wing, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "ATFL W-BC by Sex", x = "Treatment Group", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition across Lux light values 
ggplot(atfl, aes(x = lux, y = bc_wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7) + 
  labs(title = "ATFL W-BC across Lux", x = "Lux", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition across Lux light values BY SEX
ggplot(atfl, aes(x = lux, y = bc_wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "ATFL W-BC across Lux by Sex", x = "Lux", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition across LAeq noise values 
ggplot(atfl, aes(x = LAeq, y = bc_wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7) + 
  labs(title = "ATFL W-BC across Noise", x = "LAeq", y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing body condition across LAeq noise values BY SEX 
ggplot(atfl, aes(x = LAeq, y = bc_wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "ATFL W-BC across Noise by Sex", x = "LAeq", y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition between Light and non light sites
ggplot(atfl, aes(x = light, y = bc_wing)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "ATFL W-BC between Light Sites", x = "Explosed to Light?", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition between Light and non light sites BY SEX
ggplot(atfl, aes(x = light, y = bc_wing, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "ATFL W-BC by Light Site and Sex", x = "Exposed to Light?", y = "Body Condition",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition across LAeq noise values in regard to Light or No Light
ggplot(atfl, aes(x = LAeq, y = bc_wing, group = light)) +
  geom_smooth(method = "lm", aes(color = light), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = light), size = 2, alpha = 0.7) + 
  labs(title = "ATFL Wing Body Condition across Noise and Light (Y/N)", x = "LAeq", y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))




#------------------------------------------------------------------------------
#---- WING BODY CONDITION - Random Effects Model Selection ----
##  Treatment Model 
# Full model with random effects (1|section) (1|year) (1|stage) and (1|measurer)
model_15 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|year) + (1|stage) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|section) (1|year) (1|stage)
model_14 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|year) + (1|stage), data = atfl, REML = FALSE)
# Model with (1|section) (1|year) (1|measurer)
model_13 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|section) (1|stage) (1|measurer)
model_12 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|stage) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|year) (1|stage) (1|measurer)
model_11 <- lmer(bc_wing ~ trt + sex + julian + (1|year) + (1|stage) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|section) (1|year)
model_10 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|year), data = atfl, REML = FALSE)
# Model with (1|section) (1|stage)
model_9 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|stage), data = atfl, REML = FALSE)
# Model with (1|section) (1|measurer)
model_8 <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|year) (1|stage)
model_7 <- lmer(bc_wing ~ trt + sex + julian + (1|year) + (1|stage), data = atfl, REML = FALSE)
# Model with (1|year) (1|measurer)
model_6 <- lmer(bc_wing ~ trt + sex + julian + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|stage) (1|measurer)
model_5 <- lmer(bc_wing ~ trt + sex + julian + (1|stage) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|section)
model_4 <- lmer(bc_wing ~ trt + sex + julian + (1|section), data = atfl, REML = FALSE)
# Model with (1|year)
model_3 <- lmer(bc_wing ~ trt + sex + julian + (1|year), data = atfl, REML = FALSE)
# Model with (1|stage)
model_2 <- lmer(bc_wing ~ trt + sex + julian + (1|stage), data = atfl, REML = FALSE)
# Model with (1|measurer)
model_1 <- lmer(bc_wing ~ trt + sex + julian + (1|measurer), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(bc_wing ~ trt + sex + julian, data = atfl)

# Compare all models
anova(model_15, model_14, model_13, model_12, model_11, model_10, model_9, model_8, model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 1, 0, (8,6, two away)


## Lux and Noise Model 
# Full model with random effects (1|section) (1|year) (1|stage) and (1|measurer)
model_15 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|stage) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|section) (1|year) (1|stage)
model_14 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|stage), data = atfl, REML = FALSE)
# Model with (1|section) (1|year) (1|measurer)
model_13 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|section) (1|stage) (1|measurer)
model_12 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|year) (1|stage) (1|measurer)
model_11 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|section) (1|year)
model_10 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year), data = atfl, REML = FALSE)
# Model with (1|section) (1|stage)
model_9 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage), data = atfl, REML = FALSE)
# Model with (1|section) (1|measurer)
model_8 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|year) (1|stage)
model_7 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage), data = atfl, REML = FALSE)
# Model with (1|year) (1|measurer)
model_6 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|stage) (1|measurer)
model_5 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|stage) + (1|measurer), data = atfl, REML = FALSE)
# Model with (1|section)
model_4 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section), data = atfl, REML = FALSE)
# Model with (1|year)
model_3 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|year), data = atfl, REML = FALSE)
# Model with (1|stage)
model_2 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|stage), data = atfl, REML = FALSE)
# Model with (1|measurer)
model_1 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|measurer), data = atfl, REML = FALSE)
# No random effects
model_0 <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian, data = atfl)

# Compare all models
anova(model_15, model_14, model_13, model_12, model_11, model_10, model_9, model_8, model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best Models: 0, 1, (2, 3, 4 two away)



# Move forward with LMER (1|measurer)



#------------------------------------------------------------------------------
#---- WING BODY CONDITION - Model Selection across all model types ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1

# trt models
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + year + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + year + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + year + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + section + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + year + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + section + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + stage + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# Individual lux, noise, light
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(LAeq) + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(LAeq) * sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(LAeq) + sex + (1|measurer), data = atfl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(bc_wing ~ 1 + (1|measurer), data = atfl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[84]], cand.models[[75]])  
top_models

# Top model: bc_wing ~ scale(lux) + scale(LAeq) + sex + stage + (1 | measurer)
top_mod <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + stage + (1 | measurer), data = atfl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
# Lux has a SIG neg effect 
# no effect of noise
# Sex and Stage SIG effects

# 2nd Top model: bc_wing ~ scale(lux) + scale(LAeq) + sex + year + stage + (1 |measurer)
top_mod_2 <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + year + stage + (1 | measurer), data = atfl)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)
# Lux has a neg effect at 85% CI
# no effect of noise
# Sex and stage have SIG effects
# no effect of year




#### Normality Test
library(car)
qqPlot(residuals(top_mod))
shapiro.test(residuals(top_mod));length(residuals(top_mod)) # normally distributed
plot(density(resid(top_mod)))
# residuals are normal if p is > 0.05.


#### Homogeneity of variance 
# Performance Package to test for Heteroscedasticity
library(performance)
check_heteroscedasticity(top_mod_2) # model variance is homoscedastic

library(DHARMa)
simulation_output <- simulateResiduals(top_mod_2)
plot(simulation_output)
testQuantiles(simulation_output)  # This will test for homoscedasticity
testDispersion(simulation_output)  # This will test for over/underdispersion
# p-value > 0.05 means residuals are homogeneous and that there are no dispersion problems

# Manual BP (Breusch-Pagan) Test  
resids_sq <- (residuals(top_mod_2, type = "pearson"))^2 # Extract residuals and fitted values
fitted_vals <- fitted(top_mod_2)
bp_manual <- lm(resids_sq ~ fitted_vals) # Run a simple linear model to check for a relationship
summary(bp_manual) # variance is constant
# p-value > 0.05 means variance is constant (homoscedasitity)





