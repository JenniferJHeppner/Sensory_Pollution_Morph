#####################
# 2018-2025 Adult MOBL Morphology & Condition - Model Selection

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(lubridate)    # to transform dates
library(AICcmodavg)   # model selection
library(emmeans)

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

# Subset the Data to only include MOBL
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
#---- Exploratory Graphs and Sample Sizes -----
# Individuals across lux values
ggplot(mobl, aes(x = 1, y = lux)) +
  geom_jitter(width = 0.2, height = 0, color = "yellow4", alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
# 11 samples experience light - 32% of data

# Individuals across LAeq values
ggplot(mobl, aes(x = 1, y = LAeq)) +
  geom_jitter(width = 0.2, height = 0, color = "darkgreen", alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Look at Lux Values within each trt 
ggplot(mobl, aes(x= trt, y = lux, fill = trt)) +
  geom_jitter(width = 0.25, size = 2, shape = 21) +
  stat_summary(colour = "black") +
  labs(title = "Lux Values per Trt Group", x = "Trt", y = "Lux") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
# Really lacking samples that experience both light and noise

# Look at Noise values within each trt 
ggplot(mobl, aes(x= trt, y = LAeq, fill = trt)) +
  geom_jitter(width = 0.25, size = 2, shape = 21) +
  stat_summary(colour = "black") +
  labs(title = "LAeq Values per Trt Group", x = "Trt", y = "LAeq")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- WING Exploratory Graphs and Sample Sizes -----

# Remove individuals that have NA for wing
wing <- mobl %>% filter(!is.na(wing))
# 33 individuals from 18-25

# Sample sizes of individuals per Trt group and sex for ONLY individuals who have wing
ggplot(wing, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
wing %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
wing %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()



#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for WING Models ####
# Wing across sexes
ggplot(mobl, aes(x = sex, y = wing, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Wing Across Sexes",x = "Sex", y = "Wing Length")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# test for Mass across sexes
sex <- lm(wing ~ sex, data = mobl)
summary(sex)
anova(sex) # Weak effect of sex 

# Wing across ages
ggplot(mobl, aes(x = age, y = wing, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Wing Across Ages",x = "Age", y = "Wing")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Ages
age <- lm(wing ~ age, data = mobl)
summary(age)
anova(age) # NS

# Wing across years
ggplot(mobl, aes(x = year, y = wing, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Wing Across Years", x = "Year", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Years
year <- lm(wing ~ year, data = mobl)
summary(year)
anova(year) # NS difference across years

# Wing across Section
ggplot(mobl, aes(x = section, y = wing, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Wing Across Sections", x = "Section", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Sections
section <- lm(wing ~ section, data = mobl)
summary(section)
anova(section) # NS section

# Wing across Julian Date/Season
ggplot(mobl, aes(x = julian, y = wing)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "Wing Across Julian Date",
       x = "Julian Date", y = "Wing Length") +
  scale_x_continuous(limits = c(100, 160)) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Model to see tarsus across the season
julian <- lm(wing ~ julian, data = mobl)
summary(julian) # NS effect of julian
anova(julian)

# Wing across Reproductive Stages
ggplot(mobl, aes(x = stage, y = wing, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "MOBL Wing Across Reproductive Stages", x = "Stage", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

ggplot(mobl, aes(x = stage, y = wing, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              aes(color = sex), alpha = 1) +
  labs(title = "MOBL Wing Across Reproductive Stages", x = "Stage", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Linear Model of Tarsus across Reproductive Stages
stage <- lm(wing ~ stage*sex, data = mobl)
summary(stage)
anova(stage) # no effect of stage or sex

# Post Hoc Tests
post_hoc <- emmeans(stage, ~ stage | sex)
pairs(post_hoc, simple = "each")
plot(post_hoc)

# Wing across Measure ID 
ggplot(mobl, aes(x = measurer, y = wing, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Wing Across Measurers", x = "Individual Measurer", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Individual Measureres
measurer <- lm(wing ~ measurer, data = mobl)
summary(measurer)
anova(measurer) # NS effect of measurer



#------------------------------------------------------------------------------
#---- WING - MAIN FIGURES ----
# Wing across treatments 
ggplot(mobl, aes(x = trt, y = wing)) +
  geom_jitter(width = 0.2, aes(color = sex), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "MOBL Wing", x = "Treatment Group", 
       y = "Wing Length") +
  #facet_wrap(~sex) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing across treatments BY SEX
ggplot(mobl, aes(x = trt, y = wing, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "MOBL Wing by Sex", x = "Treatment Group", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing across Lux light values 
ggplot(mobl, aes(x = lux, y = wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Wing across Lux", x = "Lux", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing across Lux light values BY SEX
ggplot(mobl, aes(x = lux, y = wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Wing across Lux", x = "Lux", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing across LAeq noise values 
ggplot(mobl, aes(x = LAeq, y = wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Wing across Noise", x = "LAeq", y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing across LAeq noise values BY SEX 
ggplot(mobl, aes(x = LAeq, y = wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "MOBL Wing across Noise", x = "LAeq", y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing between Light and non light sites
ggplot(mobl, aes(x = light, y = wing)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "MOBL Wing between Light Sites", x = "Explosed to Light?", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing between Light and non light sites BY SEX
ggplot(mobl, aes(x = light, y = wing, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "MOBL Wing by Light Site and Sex",
       x = "Exposed to Light?", 
       y = "Wing Length",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#------------------------------------------------------------------------------
#---- WING ALL INDIVIDUALS - Random Effect Model Selection ----
##  Treatment Model 
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

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best Models: 0 (1,2,3 all two away)


## Lux and Noise Model 
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

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best Models: 0 (1,2,3 all two away)


# Moving forward - Linear Model




#------------------------------------------------------------------------------
#---- WING ALL INDIVIDUALS - Model Selection across all model types ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lm(wing ~ trt + sex + julian + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + julian + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + julian + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + julian + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + julian + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + julian + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + julian + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + julian + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + julian + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + julian + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ trt, data = mobl); y <- y + 1
# lux * noise
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + julian + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + julian + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + julian + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + julian + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + julian + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq) + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) * scale(LAeq), data = mobl); y <- y + 1
# lux + noise
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + julian + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + julian + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + julian + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + julian + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + julian + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + julian + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + julian + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + julian + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + julian + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + julian + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(wing ~ scale(lux) + scale(LAeq) + measurer, data = mobl); y <- y + 1
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

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[27]], cand.models[[101]], cand.models[[103]], cand.models[[32]], 
                   cand.models[[97]], cand.models[[98]], cand.models[[100]])  
top_models

# Top model:  wing ~ trt + sex
top_mod <- lm(wing ~ trt + sex, data = mobl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
# Noise has a SIG neg effect
# No effect of Light or Light+Noise
# Sex has an effect at 85% CI

# Make Reference group Light to see differences
mobl$trt <- relevel(factor(mobl$trt), ref = "Noise")

# Top model with Light as Reference group:  wing ~ trt + sex
top_mod <- lm(wing ~ trt + sex, data = mobl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
# Noise has SIG shorter wings than Light

# Remake Control the Reference group
mobl$trt <- relevel(factor(mobl$trt), ref = "Control")

# 2nd Top model: wing ~ scale(LAeq) + sex
top_mod_2 <- lm(wing ~ scale(LAeq) + sex, data = mobl)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)
# No effect of noise
# Sex has an effect at 85% CI

# 3rd Top model:  NULL

# 4th Top model:  wing ~ trt
top_mod_4 <- lm(wing ~ trt, data = mobl)
summary(top_mod_4)
confint(top_mod_4)
confint(top_mod_4, level = 0.85)
# Noise has a SIG neg effect

# Post Hoc Tests
post_hoc_mod_4 <- emmeans(top_mod_4, list(pairwise ~ trt))
pairs(post_hoc_mod_4, simple = "each")
plot(post_hoc_mod_4)


# 5th Top model: wing ~ scale(lux)
top_mod_5 <- lm(wing ~ scale(lux), data = mobl)
summary(top_mod_5)
confint(top_mod_5)
confint(top_mod_5, level = 0.85)
# No effect of lux

# 6th Top model: wing ~ scale(lux) + sex
top_mod_6 <- lm(wing ~ scale(lux) + sex, data = mobl)
summary(top_mod_6)
confint(top_mod_6)
confint(top_mod_6, level = 0.85)
# No effect of lux
# Sex has effect at 85% CI

# 7th Top model: wing ~ scale(LAeq)
top_mod_7 <- lm(wing ~ scale(LAeq), data = mobl)
summary(top_mod_7)
confint(top_mod_7)
confint(top_mod_7, level = 0.85)
# No effect of noise


#### Normality Test
library(car)
shapiro.test(residuals(top_mod_4));length(residuals(top_mod_4)) # residuals are normal
qqPlot(residuals(top_mod_4)) 
plot(density(resid(top_mod_4)))
# residuals are top_mod_4 if p is > 0.05.


#### Homogeneity of variance 
# Performance Package to test for Heteroscedasticity
library(performance)
check_heteroscedasticity(top_mod_4) # model variance is homoscedastic

library(DHARMa)
simulation_output <- simulateResiduals(top_mod_2)
plot(simulation_output)
testDispersion(simulation_output)  # This will test for over/underdispersion
# p-value > 0.05 means residuals are homogeneous and that there are no dispersion problems

# Manual BP (Breusch-Pagan) Test  
resids_sq <- (residuals(top_mod_2, type = "pearson"))^2 # Extract residuals and fitted values
fitted_vals <- fitted(top_mod_2)
bp_manual <- lm(resids_sq ~ fitted_vals) # Run a simple linear model to check for a relationship
summary(bp_manual) # variance is constant
# p-value > 0.05 means variance is constant (homoscedasitity)





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Tarsus Length -----

# Remove individuals that have NA for tarsus
tarsus <- mobl %>% filter(!is.na(tarsus))
# 33 individuals from 18-25

# Sample sizes of individuals per Trt group and sex for ONLY individuals who have tarsus
ggplot(tarsus, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
tarsus %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
tarsus %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()



#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for Tarsus Models ####
# Tarsus across sexes
ggplot(mobl, aes(x = sex, y = tarsus, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Tarsus Across Sexes",x = "Sex", y = "Tarsus")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear model of tarsus across sex
sex <- lm(tarsus ~ sex, data = mobl)
summary(sex)
anova(sex)  

# Tarsus across ages
ggplot(mobl, aes(x = age, y = tarsus, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Tarsus Across Ages",x = "Age", y = "Tarsus")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Tarsus across Ages
age <- lm(tarsus ~ age, data = mobl)
summary(age)
anova(age) # NS

# Tarsus across years
ggplot(mobl, aes(x = year, y = tarsus, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Tarsus Across Years", x = "Year", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Tarsus across Years
year <- lm(tarsus ~ year, data = mobl)
summary(year)
anova(year) # NS difference across years

# Tarsus across Section
ggplot(mobl, aes(x = section, y = tarsus, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Tarsus Across Sections", x = "Section", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Tarsus across Sections
section <- lm(tarsus ~ section, data = mobl)
summary(section)
anova(section) # No effect of section

# Tarsus across Reproductive Stages
ggplot(mobl, aes(x = stage, y = tarsus, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "MOBL Tarsus Across Reproductive Stages", x = "Stage", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

ggplot(mobl, aes(x = stage, y = tarsus, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              aes(color = sex), alpha = 1) +
  labs(title = "MOBL Tarsus Across Reproductive Stages", x = "Stage", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Linear Model of Tarsus across Reproductive Stages
stage <- lm(tarsus ~ stage*sex, data = mobl)
summary(stage)
anova(stage) # no effect of stage

# Post Hoc Tests
post_hoc <- emmeans(stage, ~ stage| sex)
pairs(post_hoc, simple = "each")
plot(post_hoc)

# Tarsus across Measure ID 
ggplot(mobl, aes(x = measurer, y = tarsus, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Tarsus Across Measurers", x = "Individual Measurer", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing Body Condition across Individual Measureres
measurer <- lm(tarsus ~ measurer, data = mobl)
summary(measurer)
anova(measurer) # No effect of measurer
confint(measurer)
confint(measurer, level = 0.85)



#------------------------------------------------------------------------------
#---- TARSUS - MAIN FIGURES ----
# Tarsus across treatments 
ggplot(mobl, aes(x = trt, y = tarsus)) +
  geom_jitter(width = 0.2, aes(color = sex), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "MOBL Tarsus", x = "Treatment Group", 
       y = "Tarsus Length") +
  #facet_wrap(~sex) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Tarsus across treatments BY SEX
ggplot(mobl, aes(x = trt, y = tarsus, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "MOBL Tarsus by Sex", x = "Treatment Group", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus across Lux light values 
ggplot(mobl, aes(x = lux, y = tarsus)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Tarsus across Light", x = "Lux", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Tarsus across Lux light values BY SEX
ggplot(mobl, aes(x = lux, y = tarsus, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Tarsus across Light", x = "Lux", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus across LAeq noise values 
ggplot(mobl, aes(x = LAeq, y = tarsus)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Tarsus across Noise", x = "LAeq", y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Tarsus across LAeq noise values BY SEX 
ggplot(mobl, aes(x = LAeq, y = tarsus, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "MOBL Tarsus across Noise", x = "LAeq", y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus between Light and non light sites
ggplot(mobl, aes(x = light, y = tarsus)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "MOBL Tarsus between Light Sites", x = "Explosed to Light?", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus between Light and non light sites BY SEX
ggplot(mobl, aes(x = light, y = tarsus, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "MOBL Tarsus by Light Site and Sex",
       x = "Exposed to Light?", 
       y = "Tarsus Length",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


#------------------------------------------------------------------------------
#---- TARSUS - Random Effects Model Selection ----
##  Treatment Model 
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

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 0, (1,2,3 two away)


## Lux and Noise Model 
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

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 0, (1,2,3 all two away)

# Moving forward - LM 




#------------------------------------------------------------------------------
#---- TARSUS - Model Selection across all model types ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lm(tarsus ~ trt + sex + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + sex + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + sex + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + sex + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + sex + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + sex + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + sex + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ trt, data = mobl); y <- y + 1
# lux * noise
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq) + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) * scale(LAeq), data = mobl); y <- y + 1
# lux + noise
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + sex + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + sex + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + sex + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + sex + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + section + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + sex + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + sex + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + sex + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + section + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + section + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + year + measurer, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(tarsus ~ scale(lux) + scale(LAeq) + measurer, data = mobl); y <- y + 1
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

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[55]], cand.models[[52]])  
top_models

# Top model: NUll

# 2nd Top model: tarsus ~ scale(LAeq)
top_mod_2 <- lm(tarsus ~ scale(LAeq), data = mobl)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)
# No effect of noise


#### Normality Test
library(car)
shapiro.test(residuals(top_mod_2));length(residuals(top_mod_2)) # residuals are normal
qqPlot(residuals(top_mod_2)) 
plot(density(resid(top_mod_2)))
# residuals are normal if p is > 0.05.


#### Homogeneity of variance 
# Performance Package to test for Heteroscedasticity
library(performance)
check_heteroscedasticity(top_mod_2) # model variance is homoscedastic

library(DHARMa)
simulation_output <- simulateResiduals(top_mod_2)
plot(simulation_output)
testDispersion(simulation_output)  # This will test for over/underdispersion
# p-value > 0.05 means residuals are homogeneous and that there are no dispersion problems

# Manual BP (Breusch-Pagan) Test  
resids_sq <- (residuals(top_mod_2, type = "pearson"))^2 # Extract residuals and fitted values
fitted_vals <- fitted(top_mod_2)
bp_manual <- lm(resids_sq ~ fitted_vals) # Run a simple linear model to check for a relationship
summary(bp_manual) # variance is constant
# p-value > 0.05 means variance is constant (homoscedasitity)




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Body Mass -----

# Remove individuals that have NA for mass
mass <- mobl %>% filter(!is.na(mass))
# 33 individuals from 18-25

# Sample sizes of individuals per Trt group and sex for ONLY individuals who have mass
ggplot(mass, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
mass %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
mass %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()



#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for Mass Models ####
# Mass across sexes
ggplot(mobl, aes(x = sex, y = mass, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Mass Across Sexes",x = "Sex", y = "Mass (g)")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# test for Mass across sexes
sex <- lm(mass ~ sex, data = mobl)
summary(sex)
anova(sex) # SIG

# Mass across ages
ggplot(mobl, aes(x = age, y = mass, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Mass Across Ages",x = "Age", y = "Mass (g)")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Ages
age <- lm(mass ~ age, data = mobl)
summary(age)
anova(age) # NS

# Mass across years
ggplot(mobl, aes(x = year, y = mass, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Mass Across Years", x = "Year", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Years
year <- lm(mass ~ year, data = mobl)
summary(year)
anova(year) # NS difference across years

# Mass across Section
ggplot(mobl, aes(x = section, y = mass, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Mass Across Sections", x = "Section", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Sections
section <- lm(mass ~ section, data = mobl)
summary(section)
anova(section) # Weak effect of Section

# Mass across Julian Date/Season
ggplot(mobl, aes(x = julian, y = mass)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "Mass Across Julian Date",
       x = "Julian Date", y = "Mass (g)") +
  scale_x_continuous(limits = c(100, 160)) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Model to see mass across the season
julian <- lm(mass ~ julian, data = mobl)
summary(julian) # NS
anova(julian)

# Mass Across Reproductive Stages
ggplot(mobl, aes(x = stage, y = mass, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "MOBL Mass Across Reproductive Stages", x = "Stage", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Mass across Reproductive Stages by sex
ggplot(mobl, aes(x = stage, y = mass, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              aes(color = sex), alpha = 1) +
  labs(title = "MOBL Mass Across Reproductive Stages", x = "Stage", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Reproductive Stages
stage <- lm(mass ~ stage * sex, data = mobl)
summary(stage)
anova(stage) # SIG effect of stage, NS sex
# egg stage is SIG higher

# Post Hoc Tests
post_hoc <- emmeans(stage, ~ stage | sex)
pairs(post_hoc, simple = "each")
plot(post_hoc)
# F - SIG diff bet TE and Eggs and bet Eggs and Chicks

# Reproductive stages across Julian Date
ggplot(mobl, aes(x = julian, y = stage, fill = stage)) +
  geom_boxplot(alpha = 0.5, width = 0.3) +
  geom_jitter(width = 0.1, alpha = 1) +
  labs(title = "MOBL Reproductive Stage by Julian Day",
       x = "Julian Date",
       y = "Stage") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Correlation test between julain date and stage
corr_test <- aov(julian ~ stage, data = mobl)
summary(corr_test) # Julian date and stage are highly correlated

# Mass across Measure ID
ggplot(mobl, aes(x = measurer, y = mass, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL Mass Across Measurers", x = "Individual Measurer", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Individual Measureres
measurer <- lm(mass ~ measurer, data = mobl)
summary(measurer)
anova(measurer) # No effect of individual measurer



#------------------------------------------------------------------------------
#---- MASS - MAIN FIGURES ----
# Mass across treatments 
ggplot(mobl, aes(x = trt, y = mass)) +
  geom_jitter(width = 0.2, aes(color = sex), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "MOBL Mass", x = "Treatment Group", 
       y = "Mass (g)") +
  #facet_wrap(~sex) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Mass across treatments BY SEX
ggplot(mobl, aes(x = trt, y = mass, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "MOBL Mass by Sex", x = "Treatment Group", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across Lux light values 
ggplot(mobl, aes(x = lux, y = mass)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Mass across Light", x = "Lux", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across Lux light values BY SEX
ggplot(mobl, aes(x = lux, y = mass, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Mass across Light", x = "Lux", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across LAeq noise values 
ggplot(mobl, aes(x = LAeq, y = mass)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Mass across Noise", x = "LAeq", y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across LAeq noise values BY SEX 
ggplot(mobl, aes(x = LAeq, y = mass, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "MOBL Mass across Noise", x = "LAeq", y = "Mass(g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


# Mass between Light and non light sites
ggplot(mobl, aes(x = light, y = mass)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5, aes(color = sex)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "MOBL Mass between Light Sites", x = "Explosed to Light?", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass between Light and non light sites BY SEX
ggplot(mobl, aes(x = light, y = mass, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "MOBL Mass by Light Site and Sex",
       x = "Exposed to Light?", 
       y = "Mass (g)",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#------------------------------------------------------------------------------
#---- MASS - Test Random Effects ----
##  Treatment Model 
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|year), data = mobl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(mass ~ trt + sex + julian + (1|section) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(mass ~ trt + sex + julian + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(mass ~ trt + sex + julian + (1|section), data = mobl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(mass ~ trt + sex + julian + (1|year), data = mobl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(mass ~ trt + sex + julian + (1|stage), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(mass ~ trt + sex + julian, data = mobl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 4, 5, 1, 0, 7


## Lux and Noise Model 
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|year), data = mobl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage), data = mobl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage), data = mobl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|section), data = mobl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year), data = mobl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|stage), data = mobl, REML = FALSE)
# No random effects
model_0 <- lm(mass ~ scale(lux) * scale(LAeq) + sex + julian, data = mobl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best model: 4, 7


# Moving forward - LMER (1|year) + (1|stage)





#------------------------------------------------------------------------------
#---- MASS - Model Selection across all model types ----

# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(mass ~ trt + sex + julian + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + julian + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + julian + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + julian + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + julian + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + julian + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + section + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lmer(mass ~ scale(lux) + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + sex + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * sex + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + sex + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) * sex + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(mass ~ 1 + (1|year) + (1|stage), data = mobl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table1 <- aictab(cand.models, modnames = modnames) 
print(aic_table1)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[26]], cand.models[[25]], cand.models[[21]])  
top_models

# Top model: mass ~ scale(lux) + sex + (1 | year) + (1 | stage)
top_mod <- lmer(mass ~ scale(lux) + sex + (1 | year) + (1 | stage), data = mobl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
# Lux has a neg effect at 85% CI
# Sex has a SIG effect

# 2nd Top model: mass ~ scale(lux) + (1 | year) + (1 | stage)
top_mod_2 <- lmer(mass ~ scale(lux) + (1 | year) + (1 | stage), data = mobl)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)
# Lux has SIG neg effect

# 3rd Top model: mass ~ scale(lux) + scale(LAeq) + sex + (1 | year) + (1 | stage)
top_mod_3 <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + (1 | year) + (1 | stage), data = mobl)
summary(top_mod_3)
confint(top_mod_3)
confint(top_mod_3, level = 0.85)
# Lux has a neg effect at 85% CI
# No effect of noise
# Sex has an effect at 85% CI



#### Normality Test
library(car)
shapiro.test(residuals(top_mod_3));length(residuals(top_mod_3))
qqPlot(residuals(top_mod_3)) 
plot(density(resid(top_mod_3)))
# residuals are normal if p is > 0.05.
# residuals are normally distributed. 


#### Homogeneity of variance 
# Performance Package to test for Heteroscedasticity
library(performance)
check_heteroscedasticity(top_mod_3) # model variance is homoscedastic

library(DHARMa)
simulation_output <- simulateResiduals(top_mod_3)
plot(simulation_output)
testQuantiles(simulation_output)  # This will test for homoscedasticity
testDispersion(simulation_output)  # This will test for over/underdispersion
# p-value > 0.05 means residuals are homogeneous and that there are no dispersion problems

# Manual BP (Breusch-Pagan) Test  
resids_sq <- (residuals(top_mod_3, type = "pearson"))^2 # Extract residuals and fitted values
fitted_vals <- fitted(top_mod_3)
bp_manual <- lm(resids_sq ~ fitted_vals) # Run a simple linear model to check for a relationship
summary(bp_manual) # variance is constant
# p-value > 0.05 means variance is constant (homoscedasitity)







##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
#---- Wing Body Condition ----
# Calculate Body Condition using mass and wing
BCwing <- smi(mobl$mass, mobl$wing)
mobl$bc_wing <- BCwing

# Mass by Wing
ggplot(mobl, aes(x = mass, y = wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Mass by Wing", x = "Mass", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")


# Mass by Wing by Sex
ggplot(mobl, aes(x = mass, y = wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Mass by Wing", x = "Mass", y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#------------------------------------------------------------------------------
#---- Exploratory Graphs and Sample Sizes -----
# Remove individuals that have NA for mass or wing
bc_w <- mobl %>% filter(!is.na(mass), !is.na(wing))
# 32 individuals from 2024-2025

# Sample sizes of individuals per Trt group and sex that have BOTH mass and Wing
ggplot(bc_w, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
bc_w %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
bc_w %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()




#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for Wing Body Condition Models ####
# Body condition across sexes
ggplot(mobl, aes(x = sex, y = bc_wing, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL W-BC Across Sexes",x = "Sex", y = "Body Condition")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model for Body Condition across sexes
sex <- lm(bc_wing ~ sex, data = mobl)
summary(sex)
anova(sex) # NS

# Wing BC across ages
ggplot(mobl, aes(x = age, y = bc_wing, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL W-BC Across Ages",x = "Age", y = "Body Condition")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing BC across Ages
age <- lm(bc_wing ~ age, data = mobl)
summary(age)
anova(age) # NS

# Body Condition across years
ggplot(mobl, aes(x = year, y = bc_wing, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL W-BC Across Years", x = "Year", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Body Condition across Years
year <- lm(bc_wing ~ year, data = mobl)
summary(year)
anova(year) # No effect of year

# Body Condition across Section
ggplot(mobl, aes(x = section, y = bc_wing, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL W-BC Across Sections", x = "Section", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Sections
section <- lm(bc_wing ~ section, data = mobl)
summary(section)
anova(section) # No effect of section

# Body Condition across Julian Date/Season
ggplot(mobl, aes(x = julian, y = bc_wing)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "MOBL W-BC Across Julian Date",
       x = "Julian Date", y = "Body Condition") +
  scale_x_continuous(limits = c(100, 160)) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Model to see body condition across the season
julian <- lm(bc_wing ~ julian, data = mobl)
summary(julian) # no effect of date
anova(julian)

# Wing Body Condition across Reproductive Stages
ggplot(mobl, aes(x = stage, y = bc_wing, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "MOBL W-BC Across Reproductive Stages", x = "Stage", 
       y = "Wing Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

ggplot(mobl, aes(x = stage, y = bc_wing, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              aes(color = sex), alpha = 1) +
  labs(title = "MOBL W-BC Across Reproductive Stages", x = "Stage", 
       y = "Wing Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing Body Condition across Reproductive Stages
stage <- lm(bc_wing ~ stage*sex, data = mobl)
summary(stage)
anova(stage) # NS

# Post Hoc Tests
post_hoc <- emmeans(stage, ~ stage | sex)
pairs(post_hoc, simple = "each")
plot(post_hoc) # No differences

# Wing Body Condition across Measure ID
ggplot(mobl, aes(x = measurer, y = bc_wing, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "MOBL W-BC Across Measurers", x = "Individual Measurer", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model ofWing  Body Condition across Individual Measureres
measurer <- lm(bc_wing ~ measurer, data = mobl)
summary(measurer)
anova(measurer) # No effect of individual measurer


#------------------------------------------------------------------------------
#---- MAIN FIGURES ----
# Wing Body condition across treatments 
ggplot(mobl, aes(x = trt, y = bc_wing, color = trt)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "MOBL W-BC", x = "Treatment Group", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing Body condition across treatments BY SEX
ggplot(mobl, aes(x = trt, y = bc_wing, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "MOBL Wing Body Condition by Sex", x = "Treatment Group", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition across Lux light values 
ggplot(mobl, aes(x = lux, y = bc_wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7) + 
  labs(title = "MOBL W-BC across Lux", x = "Lux", y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition across Lux light values BY SEX
ggplot(mobl, aes(x = lux, y = bc_wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "MOBL Wing Body Condition across Lux", x = "Lux", y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition across LAeq noise values 
ggplot(mobl, aes(x = LAeq, y = bc_wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7) + 
  labs(title = "MOBL W-BC across Noise", x = "LAeq", y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing body condition across LAeq noise values BY SEX 
ggplot(mobl, aes(x = LAeq, y = bc_wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "MOBL Wing Body Condition across Noise", x = "LAeq", y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition between Light and non light sites
ggplot(mobl, aes(x = light, y = bc_wing)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "MOBL Wing Body Condition between Light Sites", x = "Explosed to Light?", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#------------------------------------------------------------------------------
#---- WING BODY CONDITION - Test Random Effects ----
##  Treatment Model 
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

# Compare all models
anova(model_15, model_14, model_13, model_12, model_11, model_10, model_9, model_8, model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 0 (1,2,3,4 all two away)


## Lux and Noise Model 
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

# Compare all models
anova(model_15, model_14, model_13, model_12, model_11, model_10, model_9, model_8, model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best Models: 0, 4, (1,2,3 two away)

# Move forward with LM



#------------------------------------------------------------------------------
#---- WING BODY CONDITION - Model Selection across all model types ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1

# trt models
cand.models[[y]] <- lm(bc_wing ~ trt + sex + julian + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + sex + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + sex + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + sex + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + julian + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + sex + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + sex + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ trt * sex , data = mobl); y <- y + 1
# lux * noise
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + julian + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq) + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) * scale(LAeq), data = mobl); y <- y + 1
# lux + noise
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + sex + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + julian + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + sex + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + sex + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + julian + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + julian + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + year + section, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + sex, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + julian, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + year, data = mobl); y <- y + 1
cand.models[[y]] <- lm(bc_wing ~ scale(lux) + scale(LAeq) + section, data = mobl); y <- y + 1
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

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[12]], cand.models[[16]])  
top_models

# Top model: bc_wing ~ trt + sex
top_mod_1 <- lm(bc_wing ~ trt + sex, data = mobl)
summary(top_mod_1)
confint(top_mod_1)
confint(top_mod_1, level = 0.85)
# Light only SIG neg effects
# Noise only SIG neg effects
# Light and Noise has a neg effect at 85% CI
# Sex has an effect at 85% CI

# 2nd Top model: bc_wing ~ trt
top_mod_2 <- lm(bc_wing ~ trt, data = mobl)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)
# Light only SIG neg effects
# Noise only SIG neg effects


#### Normality Test
library(car)
qqPlot(residuals(top_mod_2))
shapiro.test(residuals(top_mod_2));length(residuals(top_mod_2)) # normally distributed
plot(density(resid(top_mod_2)))
# residuals are normal if p is > 0.05.


#### Homogeneity of variance 
# Performance Package to test for Heteroscedasticity
library(performance)
check_heteroscedasticity(top_mod_2) # model variance is homoscedastic

library(DHARMa)
simulation_output <- simulateResiduals(top_mod_1)
plot(simulation_output)
testQuantiles(simulation_output)  # This will test for homoscedasticity
testDispersion(simulation_output)  # This will test for over/underdispersion
# p-value > 0.05 means residuals are homogeneous and that there are no dispersion problems

# Manual BP (Breusch-Pagan) Test  
resids_sq <- (residuals(top_mod_1, type = "pearson"))^2 # Extract residuals and fitted values
fitted_vals <- fitted(top_mod_1)
bp_manual <- lm(resids_sq ~ fitted_vals) # Run a simple linear model to check for a relationship
summary(bp_manual) # variance is constant
# p-value > 0.05 means variance is constant (homoscedasitity)






