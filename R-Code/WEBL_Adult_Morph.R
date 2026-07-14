#####################
# 2018-2025 Adult WEBL Morphology & Condition - Model Selection

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(lubridate)    # to transform dates
library(AICcmodavg)   # model selection
library(emmeans)
library(ggeffects)

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

# Reorder the levels of the age factor
dat <- dat %>%
  mutate(age = factor(age, levels = c("SY", "AHY", "ASY")))


# Subset the Data to only include WEBL
webl <- dat %>% filter(species == "WEBL") %>% droplevels()


# Remove outlier -  Does not change results
webl <- webl %>%
  mutate(tarsus = ifelse(ID %in% c("3011-18015"), NA, tarsus)) %>%
  mutate(bill.l = ifelse(ID %in% c("2791-95473", "2941-73548", "2791-95475"), NA, bill.l))%>%
  mutate(bill.w = ifelse(ID %in% c("2941-73548"), NA, bill.w))


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
#---- Exploratory Graphs and Sample Sizes -----
# Look at Lux Values within each trt 
ggplot(webl, aes(x= trt, y = lux, fill = trt)) +
  geom_jitter(width = 0.25, size = 2, shape = 21) +
  stat_summary(colour = "black") +
  labs(title = "WEBL Lux Values per Trt Group", x = "Trt", y = "Lux") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# Look at Noise values within each trt 
ggplot(webl, aes(x= trt, y = LAeq, fill = trt)) +
  geom_jitter(width = 0.25, size = 2, shape = 21) +
  stat_summary(colour = "black") +
  labs(title = "WEBL LAeq Values per Trt Group", x = "Trt", y = "LAeq")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# Individuals across lux values
ggplot(webl, aes(x = 1, y = lux)) +
  geom_jitter(width = 0.2, height = 0, color = "yellow4", alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Individuals across LAeq values
ggplot(webl, aes(x = 1, y = LAeq)) +
  geom_jitter(width = 0.2, height = 0, color = "darkgreen", alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- WING Length -----

# Remove individuals that have NA for wing
wing <- webl %>% filter(!is.na(wing))
# 47 individuals from 18-25

# Sample sizes of individuals per Trt group and sex for ONLY individuals who have wing
ggplot(wing, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
wing %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
wing %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()




#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for WING Models ####
# Wing across sexes
ggplot(webl, aes(x = sex, y = wing, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Wing Across Sexes",x = "Sex", y = "Wing Length")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# test for Mass across sexes
sex <- lm(wing ~ sex, data = webl)
summary(sex)
anova(sex) # Weak effect of sex 

# Wing across ages
ggplot(webl, aes(x = age, y = wing, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Wing Across Ages",x = "Age", y = "Wing")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Ages
age <- lm(wing ~ age, data = webl)
summary(age)
anova(age) # NS

# Wing across years
ggplot(webl, aes(x = year, y = wing, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Wing Across Years", x = "Year", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Years
year <- lm(wing ~ year, data = webl)
summary(year)
anova(year) # NS difference across years

# Wing across Section
ggplot(webl, aes(x = section, y = wing, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Wing Across Sections", x = "Section", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Sections
section <- lm(wing ~ section, data = webl)
summary(section)
anova(section) # NS difference across sections

# Wing across Julian Date/Season
ggplot(webl, aes(x = julian, y = wing)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "WEBL Wing Across Julian Date",
       x = "Julian Date", y = "Wing Length") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Model to see tarsus across the season
julian <- lm(wing ~ julian, data = webl)
summary(julian) # NS effect of julian
anova(julian)

# Wing across Reproductive Stages
ggplot(webl, aes(x = stage, y = wing, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "WEBL Wing Across Reproductive Stages", x = "Stage", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

ggplot(webl, aes(x = stage, y = wing, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              aes(color = sex), alpha = 1) +
  labs(title = "WEBL Wing Across Reproductive Stages", x = "Stage", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Linear Model of Tarsus across Reproductive Stages
stage <- lm(wing ~ stage, data = webl)
summary(stage)
anova(stage) # No stage effect

stage <- lm(wing ~ stage*sex, data = webl)
summary(stage)
anova(stage) # no effect of stage or sex

# Post Hoc Tests
post_hoc <- emmeans(stage, ~ stage | sex)
pairs(post_hoc, simple = "each")
plot(post_hoc)

# Wing across Measure ID 
ggplot(webl, aes(x = measurer, y = wing, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Wing Across Measurers", x = "Individual Measurer", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing across Individual Measureres
measurer <- lm(wing ~ measurer, data = webl)
summary(measurer)
anova(measurer) # SIG effect of measurer




#------------------------------------------------------------------------------
#---- WING - MAIN FIGURES ----
# Wing across treatments 
ggplot(webl, aes(x = trt, y = wing)) +
  geom_jitter(width = 0.2, aes(color = sex), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "WEBL Wing", x = "Treatment Group", 
       y = "Wing Length") +
  #facet_wrap(~sex) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing across treatments BY SEX
ggplot(webl, aes(x = trt, y = wing, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "WEBL Wing by Sex", x = "Treatment Group", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing across treatments by AGE
ggplot(webl, aes(x = trt, y = wing)) +
  geom_jitter(width = 0.2, aes(color = age), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "WEBL Wing", x = "Treatment Group", 
       y = "Wing Length") +
  #facet_wrap(~sex) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing across treatments BY AGE
ggplot(webl, aes(x = trt, y = wing, group = age)) +
  geom_jitter(aes(color = age), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "WEBL Wing by Age", x = "Treatment Group", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing across Lux light values 
ggplot(webl, aes(x = lux, y = wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Wing across Lux", x = "Lux", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing across Lux light values  BY SEX
ggplot(webl, aes(x = lux, y = wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Wing across Lux", x = "Lux", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing across LAeq noise values 
ggplot(webl, aes(x = LAeq, y = wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Wing across Noise", x = "LAeq", y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing across LAeq noise values BY SEX 
ggplot(webl, aes(x = LAeq, y = wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "WEBL Wing across Noise", x = "LAeq", y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing between Light and non light sites
ggplot(webl, aes(x = light, y = wing)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "WEBL Wing between Light Sites", x = "Explosed to Light?", 
       y = "Wing Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing between Light and non light sites BY SEX
ggplot(webl, aes(x = light, y = wing, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "WEBL Wing by Light Site and Sex",
       x = "Exposed to Light?", 
       y = "Wing Length",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#------------------------------------------------------------------------------
#---- WING ALL INDIVIDUALS - Test Random Effects ----
##  Treatment Model 
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

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best Models: 1, (4,5 all two away)


## Lux and Noise Model 
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

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best Models: 1,0, (4,5 two away)


# Moving forward - LMER (1|measurer)



#------------------------------------------------------------------------------
#---- WING ALL INDIVIDUALS LMER (1|measurer) - Model Selection across all model types ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(wing ~ trt + sex + age + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + age + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + age + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + age + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + age + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + age + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + age + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + age + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ trt + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + age + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + age + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + age + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + age + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) * scale(LAeq) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + age + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + age + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + age + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + age + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + julian + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + age + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + julian + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + julian + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + year + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + julian + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(wing ~ scale(lux) + scale(LAeq) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# Individual lux, noise
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

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[103]], cand.models[[17]], cand.models[[98]], cand.models[[81]], cand.models[[27]])  
top_models

# Top model:  wing ~ scale(LAeq) + sex + (1 | measurer)
top_mod <- lmer(wing ~ scale(LAeq) + sex + (1 | measurer), data = webl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
# No effect of noise
# Sex has SIG effect

# 2nd Top model: wing ~ trt + sex + age + (1 | measurer)
top_mod_2 <- lmer(wing ~ trt + sex + age + (1 | measurer), data = webl)
summary(top_mod_2)
anova(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)
# Light has a SIG neg effect 
# Noise has a SIG neg effect 
# No effect of Light+Noise 
# Sex is SIG
# Age is SIG 

# Find means of wing by Age
adj_means <- emmeans(top_mod_2, ~ age)
adj_means
pairs(adj_means)

# 3rd Top model: wing ~ scale(lux) + sex + (1 | measurer)
top_mod_3 <- lmer(wing ~ scale(lux) + sex + (1 | measurer), data = webl)
summary(top_mod_3)
confint(top_mod_3)
confint(top_mod_3, level = 0.95)
# No effect of lux
# Sex has SIG effect

# 4th Top model: wing ~ scale(lux) + scale(LAeq) + sex + age + (1 | measurer)
top_mod_4 <- lmer(wing ~ scale(lux) + scale(LAeq) + sex + age + (1 | measurer), data = webl)
summary(top_mod_4)
confint(top_mod_4)
confint(top_mod_4, level = 0.85)
# no effect of lux or noise
# Sex is SIG
# Age at ASY is SIG

# 5th Top model: wing ~ trt + sex + (1 | measurer)
top_mod_5 <- lmer(wing ~ trt + sex + (1 | measurer), data = webl)
summary(top_mod_5)
confint(top_mod_5)
confint(top_mod_5, level = 0.85)
# Light has neg effect at 85%
# Noise has neg effect at 95%
# Light+Noise has neg effect at 85%
# Sex is SIG


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
simulation_output <- simulateResiduals(top_mod_4)
plot(simulation_output)
testQuantiles(simulation_output)  # This will test for homoscedasticity
testDispersion(simulation_output)  # This will test for over/underdispersion
# p-value > 0.05 means residuals are homogeneous and that there are no dispersion problems

# Manual BP (Breusch-Pagan) Test  
resids_sq <- (residuals(top_mod_4, type = "pearson"))^2 # Extract residuals and fitted values
fitted_vals <- fitted(top_mod_4)
bp_manual <- lm(resids_sq ~ fitted_vals) # Run a simple linear model to check for a relationship
summary(bp_manual) # variance is constant
# p-value > 0.05 means variance is constant (homoscedasitity)






#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#---- Tarsus Length -----

# Remove individuals that have NA for tarsus
tarsus <- webl %>% filter(!is.na(tarsus))
# 46 individuals from 18-25

# Sample sizes of individuals per Trt group and sex for ONLY individuals who have tarsus
ggplot(tarsus, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
tarsus %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
tarsus %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()



#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for Tarsus Models ####
# Tarsus across sexes
ggplot(webl, aes(x = sex, y = tarsus, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Tarsus Across Sexes",x = "Sex", y = "Tarsus")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear model of tarsus across sex
sex <- lm(tarsus ~ sex, data = webl)
summary(sex)
anova(sex)  # NS

# Tarsus across ages
ggplot(webl, aes(x = age, y = tarsus, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Tarsus Across Ages",x = "Age", y = "Tarsus")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Tarsus across Ages
age <- lm(tarsus ~ age, data = webl)
summary(age)
anova(age) # NS

# Tarsus across years
ggplot(webl, aes(x = year, y = tarsus, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Tarsus Across Years", x = "Year", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Tarsus across Years
year <- lm(tarsus ~ year, data = webl)
summary(year)
anova(year) # NS difference across years

# Tarsus across Section
ggplot(webl, aes(x = section, y = tarsus, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Tarsus Across Sections", x = "Section", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Tarsus across Sections
section <- lm(tarsus ~ section, data = webl)
summary(section)
anova(section) # No effect of section

# Tarsus across Julian Date/Season
ggplot(webl, aes(x = julian, y = tarsus)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "Tarsus Across Julian Date",
       x = "Julian Date", y = "Tarsus") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Model to see tarsus across the season
julian <- lm(tarsus ~ julian, data = webl)
summary(julian) # NS effect of julian
anova(julian)

# Tarsus across Reproductive Stages
ggplot(webl, aes(x = stage, y = tarsus, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "WEBL Tarsus Across Reproductive Stages", x = "Stage", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

ggplot(webl, aes(x = stage, y = tarsus, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              aes(color = sex), alpha = 1) +
  labs(title = "WEBL Tarsus Across Reproductive Stages", x = "Stage", 
       y = "Tarsus") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Linear Model of Tarsus across Reproductive Stages
stage <- lm(tarsus ~ stage * sex, data = webl)
summary(stage)
anova(stage) # No effect of stage or sex

# Post Hoc Tests
post_hoc <- emmeans(stage, ~ stage | sex)
pairs(post_hoc, simple = "each")
plot(post_hoc) # Tarsus does not differ by stage or sex

# Tarsus across Measure ID 
ggplot(webl, aes(x = measurer, y = tarsus, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Tarsus Across Measurers", x = "Individual Measurer", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing Body Condition across Individual Measureres
measurer <- lm(tarsus ~ measurer, data = webl)
summary(measurer)
anova(measurer) # SIG effect of measurer



#------------------------------------------------------------------------------
#---- TARSUS - MAIN FIGURES ----
# Tarsus across treatments 
ggplot(webl, aes(x = trt, y = tarsus)) +
  geom_jitter(width = 0.2, aes(color = sex), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "WEBL Tarsus", x = "Treatment Group", 
       y = "Tarsus Length") +
  #facet_wrap(~sex) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Tarsus across treatments BY SEX
ggplot(webl, aes(x = trt, y = tarsus, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "WEBL Tarsus by Sex", x = "Treatment Group", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus across Lux light values 
ggplot(webl, aes(x = lux, y = tarsus)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Tarsus across Lux", x = "Lux", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Tarsus across Lux light values BY SEX
ggplot(webl, aes(x = lux, y = tarsus, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Tarsus across Lu", x = "Lux", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus across LAeq noise values 
ggplot(webl, aes(x = LAeq, y = tarsus)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Tarsus across Noise", x = "LAeq", y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Tarsus across LAeq noise values BY SEX 
ggplot(webl, aes(x = LAeq, y = tarsus, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "WEBL Tarsus across Noise", x = "LAeq", y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus between Light and non light sites
ggplot(webl, aes(x = light, y = tarsus)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "WEBL Tarsus between Light Sites", x = "Explosed to Light?", 
       y = "Tarsus Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Tarsus between Light and non light sites BY SEX
ggplot(webl, aes(x = light, y = tarsus, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "WEBL Tarsus by Light Site and Sex",
       x = "Exposed to Light?", 
       y = "Tarsus Length",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#------------------------------------------------------------------------------
#---- TARSUS - Random Effects Model Selection ----
##  Treatment Model 
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(tarsus ~ trt + sex + age + (1|section) + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(tarsus ~ trt + sex + age + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(tarsus ~ trt + sex + age + (1|section) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(tarsus ~ trt + sex + age + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(tarsus ~ trt + sex + age + (1|section), data = webl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(tarsus ~ trt + sex + age + (1|year), data = webl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(tarsus ~ trt + sex + age + (1|measurer), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(tarsus ~ trt + sex + age, data = webl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 1, 4, (two away: 5) 


## Lux and Noise Model 
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + (1|section) + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + (1|section) + (1|measurer), data = webl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + (1|year) + (1|measurer), data = webl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + (1|section), data = webl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + (1|year), data = webl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + (1|measurer), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(tarsus ~ scale(lux) * scale(LAeq) + sex + age, data = webl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 1, 5, (two away: 4) 

# Moving forward - LMER (1|meaasurer)




#------------------------------------------------------------------------------
#---- TARSUS - Model Selection across all model types ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(tarsus ~ trt + sex + age + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + age + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + age + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * scale(LAeq) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + age + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + age + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + age + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + age + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + section + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + section + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + year + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + scale(LAeq) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(lux) * age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) + sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) * sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) + age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ scale(LAeq) * age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt * sex + (1|measurer), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(tarsus ~ trt * age + (1|measurer), data = webl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(tarsus ~ 1 + (1|measurer), data = webl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models_T <- list(cand.models[[61]], cand.models[[54]])  
top_models_T

# Top model is Null

# Top model: tarsus ~ scale(LAeq) + (1 | measurer)
top_mod <- lmer(tarsus ~ scale(LAeq) + (1 | measurer), data = webl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
# No effect of lux


#### Normality Test
library(car)
shapiro.test(residuals(top_mod));length(residuals(top_mod)) # residuals are normal
qqPlot(residuals(top_mod)) 
plot(density(resid(top_mod)))
# residuals are normal if p is > 0.05.


#### Homogeneity of variance 
# Performance Package to test for Heteroscedasticity
library(performance)
check_heteroscedasticity(top_mod) # model variance is homoscedastic

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
#---- Body Mass -----

# Remove individuals that have NA for mass
mass <- webl %>% filter(!is.na(mass))
# 48 individuals from 18-25

# Sample sizes of individuals per Trt group and sex for ONLY individuals who have mass
ggplot(mass, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
mass %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
mass %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()



#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for Mass Models ####
# Mass across sexes
ggplot(webl, aes(x = sex, y = mass, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Mass Across Sexes",x = "Sex", y = "Mass (g)")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# test for Mass across sexes
sex <- lm(mass ~ sex, data = webl)
summary(sex)
anova(sex) # SIG

# Mass across ages
ggplot(webl, aes(x = age, y = mass, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Mass Across Ages",x = "Age", y = "Mass (g)")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Ages
age <- lm(mass ~ age, data = webl)
summary(age)
anova(age) # NS

# Mass across years
ggplot(webl, aes(x = year, y = mass, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Mass Across Years", x = "Year", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Years
year <- lm(mass ~ year, data = webl)
summary(year)
anova(year) # Year SIG affects mass

# Mass across Section
ggplot(webl, aes(x = section, y = mass, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Mass Across Sections", x = "Section", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Sections
section <- lm(mass ~ section, data = webl)
summary(section)
anova(section) # No effect of section

# Mass across Julian Date/Season
ggplot(webl, aes(x = julian, y = mass)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "WEBL Mass Across Julian Date",
       x = "Julian Date", y = "Mass (g)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Model to see mass across the season
julian <- lm(mass ~ julian, data = webl)
summary(julian) # Weak effect of julian date
anova(julian)

# Mass across Reproductive Stages
ggplot(webl, aes(x = stage, y = mass, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "WEBL Mass Across Reproductive Stages", x = "Stage", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

ggplot(webl, aes(x = stage, y = mass, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              aes(color = sex), alpha = 1) +
  labs(title = "WEBL Mass Across Reproductive Stages", x = "Stage", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Linear Model of Mass across Reproductive Stages
stage <- lm(mass ~ stage, data = webl)
summary(stage)
anova(stage) # SIG effect of stage and sex

# Post Hoc Tests
post_hoc <- emmeans(stage, ~ stage | sex)
pairs(post_hoc, simple = "each")
plot(post_hoc)

# Reproductive stages across Julian Date
ggplot(webl, aes(x = julian, y = stage, fill = stage)) +
  geom_boxplot(alpha = 0.5, width = 0.3) +
  geom_jitter(width = 0.1, alpha = 1) +
  labs(title = "WEBL Reproductive Stage by Julian Day",
       x = "Julian Date",
       y = "Stage") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Correlation test between julain date and stage
corr_test <- aov(julian ~ stage, data = webl)
summary(corr_test) # Julian date and stage are highly correlated

# Mass across Measure ID
ggplot(webl, aes(x = measurer, y = mass, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL Mass Across Measurers", x = "Individual Measurer", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Individual Measureres
measurer <- lm(mass ~ measurer, data = webl)
summary(measurer)
anova(measurer) # No effect of individual measurer



#------------------------------------------------------------------------------
#---- MASS - MAIN FIGURES ----
# Mass across treatments 
ggplot(webl, aes(x = trt, y = mass)) +
  geom_jitter(width = 0.2, aes(color = sex), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "WEBL Mass", x = "Treatment Group", 
       y = "Mass (g)") +
  #facet_wrap(~sex) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Mass across treatments BY SEX
ggplot(webl, aes(x = trt, y = mass, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "WEBL Mass by Sex", x = "Treatment Group", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across Lux light values 
ggplot(webl, aes(x = lux, y = mass)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Mass across Lux", x = "Lux", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Mass across Lux light values BY SEX
ggplot(webl, aes(x = lux, y = mass, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Mass across Lux", x = "Lux", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass across LAeq noise values 
ggplot(webl, aes(x = LAeq, y = mass)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Mass across Noise", x = "LAeq", y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Mass across LAeq noise values BY SEX 
ggplot(webl, aes(x = LAeq, y = mass, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "WEBL Mass across Noise", x = "LAeq", y = "Mass(g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass between Light and non light sites
ggplot(webl, aes(x = light, y = mass)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "WEBL Mass between Light Sites", x = "Explosed to Light?", 
       y = "Mass (g)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Mass between Light and non light sites BY SEX
ggplot(webl, aes(x = light, y = mass, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "WEBL Mass by Light Site and Sex",
       x = "Exposed to Light?", 
       y = "Mass (g)",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#------------------------------------------------------------------------------
#---- MASS - Random Effects Model Selection ----
##  Treatment Model 
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(mass ~ trt + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(mass ~ trt + sex + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(mass ~ trt + sex + age + julian + (1|section), data = webl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(mass ~ trt + sex + age + julian + (1|year), data = webl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(mass ~ trt + sex + age + julian + (1|stage), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(mass ~ trt + sex + age + julian, data = webl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 4, 1, 7


## Lux and Noise Model 
# Full model with random effects (1|section) (1|year) and (1|stage)
model_7 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with (1|section) and (1|year)
model_6 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|year), data = webl, REML = FALSE)
# Model with (1|section) and (1|stage)
model_5 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE)
# Model with (1|year) and (1|stage)
model_4 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE)
# Model with only (1|section)
model_3 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section), data = webl, REML = FALSE)
# Model with only (1|year)
model_2 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year), data = webl, REML = FALSE)
# Model with only (1|stage)
model_1 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|stage), data = webl, REML = FALSE)
# No random effects
model_0 <- lm(mass ~ scale(lux) * scale(LAeq) + sex + age + julian, data = webl)

# Compare all models
anova(model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best model: 4, 7, 1, 5


# Moving forward - LMER (1|year) + (1|stage)




#------------------------------------------------------------------------------
#---- MASS - Model Selection across all model types ----

# Create an empty list to store candidate models 
cand.models <- list()
y <- 1 # Counter for model indexing

# trt models
cand.models[[y]] <- lmer(mass ~ trt + sex + age + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + age + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + age + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + age + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + age + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + sex + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + age + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ trt + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + age + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + age + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + age + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + age + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * scale(LAeq) + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + age + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + age + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + age + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + age + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + age + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + age + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + age + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + section + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + scale(LAeq) + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# Individual lux, noise
cand.models[[y]] <- lmer(mass ~ scale(lux) + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + sex + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * sex + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) + age + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(lux) * age + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + sex + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) * sex + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) + age + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(mass ~ scale(LAeq) * age + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# Null
cand.models[[y]] <- lmer(mass ~ 1 + (1|year) + (1|stage), data = webl, REML = FALSE); y <- y + 1

# Generate names for models
modnames <- paste0("Model ", seq_along(cand.models))

# Compare models using AIC
aic_table <- aictab(cand.models, modnames = modnames) 
print(aic_table)

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[39]], cand.models[[50]], cand.models[[23]], cand.models[[55]], 
                   cand.models[[46]], cand.models[[59]], cand.models[[7]])  
top_models

# Top model: mass ~ scale(lux) + scale(LAeq) + sex + julian + (1 | year) + (1 | stage)
top_mod <- lmer(mass ~ scale(lux) + scale(LAeq) + sex + julian + (1 | year) + (1 | stage), data = webl)
summary(top_mod)
confint(top_mod)
confint(top_mod, level = 0.85)
# No effect of lux or noise
# Sex and julian have a SIG effect

# 2nd Top model: mass ~ scale(lux) + sex + (1 | year) + (1 | stage)
top_mod_2 <- lmer(mass ~ scale(lux) + sex + (1 | year) + (1 | stage), data = webl)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)
# No effect of lux
# SIG effect of sex

# 3rd Top model: mass ~ scale(lux) * scale(LAeq) + sex + julian + (1 | year) + (1 | stage)
top_mod_3 <- lmer(mass ~ scale(lux) * scale(LAeq) + sex + julian + (1 | year) + (1 | stage), data = webl)
summary(top_mod_3)
confint(top_mod_3)
confint(top_mod_3, level = 0.85)
# No effect of noise
# Sex has SIG effect

# 4th Top model: mass ~ scale(LAeq) + sex + (1 | year) + (1 | stage)
top_mod_4 <- lmer(mass ~ scale(LAeq) + sex + (1 | year) + (1 | stage), data = webl)
summary(top_mod_4)
confint(top_mod_4)
confint(top_mod_4, level = 0.85)
# No effect of trt
# SIG effect of Sex and Julian

# 5th Top model: mass ~ scale(lux) + scale(LAeq) + julian + (1 | year) + (1 | stage)
top_mod_5 <- lmer(mass ~ scale(lux) + scale(LAeq) + julian + (1 | year) + (1 | stage), data = webl)
summary(top_mod_5)
confint(top_mod_5)
confint(top_mod_5, level = 0.85)
# No effect of lux or noise
# Sig effect of julian

# 6th is null

# 7th Top model: mass ~ trt + sex + julian + (1 | year) + (1 | stage)
top_mod_7 <- lmer(mass ~ trt + sex + julian + (1 | year) + (1 | stage), data = webl)
summary(top_mod_7)
confint(top_mod_7)
confint(top_mod_7, level = 0.85)
# No effect of lux, noise, or int
# SIG effects of sex and julian




#### Normality Test
library(car)
shapiro.test(residuals(top_mod_7));length(residuals(top_mod_7))
qqPlot(residuals(top_mod_7)) 
plot(density(resid(top_mod_7)))
# residuals are normal if p is > 0.05.


#### Homogeneity of variance 
# Performance Package to test for Heteroscedasticity
library(performance)
check_heteroscedasticity(top_mod_7) # model variance is homoscedastic

library(DHARMa)
simulation_output <- simulateResiduals(top_mod_7)
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







##------------------------------------------------------------------------------
#---- Wing Body Condition ----
# Calculate Body Condition using mass and wing
BCwing <- smi(webl$mass, webl$wing)
webl$bc_wing <- BCwing

# Mass by Wing
ggplot(webl, aes(x = bc_wing, y = wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Mass by Wing", x = "bc_wing", 
       y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")


# Mass by Wing by Sex
ggplot(webl, aes(x = mass, y = wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL Mass by Wing", x = "Mass", y = "Wing") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))




#------------------------------------------------------------------------------
#---- Exploratory Graphs and Sample Sizes -----
# Remove individuals that have NA for mass or wing
bc_w <- webl %>% filter(!is.na(mass), !is.na(wing))
# 45 individuals

# Sample sizes of individuals per Trt group and sex that have BOTH mass and Tarsus
ggplot(bc_w, aes(x = trt)) + geom_bar() +
  labs(title = "Individuals per Treatment", x = "Treatment Group", y = "Count") +
  theme_minimal()
bc_w %>%  group_by(trt) %>%  summarise(sample_size = n()) %>% print()
bc_w %>%  group_by(trt, sex) %>%  summarise(sample_size = n()) %>% print()



#------------------------------------------------------------------------------
#---- Visualize Fixed Effects for Wing Body Condition Models ####
# Body condition across sexes
ggplot(webl, aes(x = sex, y = bc_wing, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL W-BC Across Sexes",x = "Sex", y = "Body Condition")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model for Body Condition across sexes
sex <- lm(bc_wing ~ sex, data = webl)
summary(sex)
anova(sex) # NS

# Wing BC across ages
ggplot(webl, aes(x = age, y = bc_wing, fill = age)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL W-BC Across Ages",x = "Age", y = "Body Condition")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Wing BC across Ages
age <- lm(bc_wing ~ age, data = webl)
summary(age)
anova(age) 

# Body Condition across years
ggplot(webl, aes(x = year, y = bc_wing, fill = year)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL W-BC Across Years", x = "Year", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Body Condition across Years
year <- lm(bc_wing ~ year, data = webl)
summary(year)
anova(year) # No effect of year

# Body Condition across Section
ggplot(webl, aes(x = section, y = bc_wing, fill = section)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = "WEBL W-BC Across Sections", x = "Section", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model of Mass across Sections
section <- lm(bc_wing ~ section, data = webl)
summary(section) # Pink section trending to reduce BC
anova(section) # No effect of section

# Body Condition across Julian Date/Season
ggplot(webl, aes(x = julian, y = bc_wing)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  
  geom_smooth(method = "loess", color = "red", se = TRUE) +  
  labs(title = "WEBL W-BC Across Julian Date",
       x = "Julian Date", y = "Body Condition") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Model to see body condition across the season
julian <- lm(bc_wing ~ julian, data = webl)
summary(julian) # no effect of date
anova(julian)

# Wing Body Condition across Reproductive Stages
ggplot(webl, aes(x = stage, y = bc_wing, fill = stage)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              alpha = 1) +
  labs(title = "WEBL W-BC Across Reproductive Stages", x = "Stage", 
       y = "Wing Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

ggplot(webl, aes(x = stage, y = bc_wing, fill = sex)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
              aes(color = sex), alpha = 1) +
  labs(title = "WEBL Wing Body Condition Across Reproductive Stages", x = "Stage", 
       y = "Wing Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Linear Model of Wing Body Condition across Reproductive Stages
stage <- lm(bc_wing ~ stage, data = webl)
summary(stage)
anova(stage) # SIG across stages and sexes alone

# Post Hoc Tests
post_hoc <- emmeans(stage, ~ stage | sex)
pairs(post_hoc, simple = "each")
plot(post_hoc)

# Wing Body Condition across Measure ID
ggplot(webl, aes(x = measurer, y = bc_wing, fill = measurer)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), color = "black",
              alpha = 0.7) +  
  labs(title = " WEBL W-BC Across Measurers", x = "Individual Measurer", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Linear Model ofWing  Body Condition across Individual Measureres
measurer <- lm(bc_wing ~ measurer, data = webl)
summary(measurer)
anova(measurer) # No effect of individual measurer



#------------------------------------------------------------------------------
#---- MAIN FIGURES ----
# Wing Body condition across treatments 
ggplot(webl, aes(x = trt, y = bc_wing)) +
  geom_jitter(width = 0.2, aes(color = sex), size = 2, alpha = 0.7) + 
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  labs(title = "WEBL W-BC", x = "Treatment Group", 
       y = "Body Condition") +
  #facet_wrap(~sex) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing Body condition across treatments BY SEX
ggplot(webl, aes(x = trt, y = bc_wing, group = sex)) +
  geom_jitter(aes(color = sex), size = 2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.7)) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.7),
               color = "black", alpha = 0.8, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.7)) +
  labs(title = "WEBL Wing Body Condition by Sex", x = "Treatment Group", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition across Lux light values 
ggplot(webl, aes(x = lux, y = bc_wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "yellow4") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL W-BC across Lux", x = "Lux", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing Body condition across Lux light values BY SEX
ggplot(webl, aes(x = lux, y = bc_wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL W-BC across Lux", x = "Lux", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition across LAeq noise values 
ggplot(webl, aes(x = LAeq, y = bc_wing)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  geom_jitter(width = 0.01, aes(color = sex), size = 2, alpha = 0.7) + 
  labs(title = "WEBL W-BC across Noise", x = "LAeq", y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Wing body condition across LAeq noise values BY SEX 
ggplot(webl, aes(x = LAeq, y = bc_wing, group = sex)) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.7, aes(color = sex)) + 
  labs(title = "WEBL Wing Body Condition across Noise", x = "LAeq", y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition between Light and non light sites
ggplot(webl, aes(x = light, y = bc_wing)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, 
               position = position_dodge(width = 0.9),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, 
               linewidth = 1, color = "black", alpha = 0.7, 
               position = position_dodge(width = 0.9)) +
  labs(title = "WEBL Wing Body Condition between Light Sites", x = "Explosed to Light?", 
       y = "Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Wing Body condition between Light and non light sites BY SEX
ggplot(webl, aes(x = light, y = bc_wing, group = sex, color = sex)) +
  geom_jitter(size = 2, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 5,
               position = position_dodge(width = 0.5),
               color = "black", alpha = 1, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1,
               color = "black", linewidth = 1, alpha = 0.7,
               position = position_dodge(width = 0.5)) +
  labs(title = "WEBL Wing Body Condition by Light Site and Sex",
       x = "Exposed to Light?", 
       y = "Body Condition",
       color = "Sex") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#------------------------------------------------------------------------------
#---- WING BODY CONDITION - Test Random Effects ----
##  Treatment Model 
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

# Compare all models
anova(model_15, model_14, model_13, model_12, model_11, model_10, model_9, model_8, model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best models: 2, 9, 12, 5,(7 are two away)



## Lux and Noise Model 
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

# Compare all models
anova(model_15, model_14, model_13, model_12, model_11, model_10, model_9, model_8, model_7, model_6, model_5, model_4, model_3, model_2, model_1, model_0)
# Best Models: 2, 9, 7, (5 and  are two away)

# Move forward with LMER (Section) (Stage) 



#------------------------------------------------------------------------------
#---- WING BODY CONDITION - Model Selection across all model types SECTION AND STAGE ----
# Create an empty list to store candidate models 
cand.models <- list()
y <- 1

# trt models
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + age + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + age + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + age + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + age + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + age + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + age + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + age + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + age + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + age + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + age + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + age + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + age + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ trt + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# lux * noise
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + age + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + age + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + age + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + age + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + age + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + age + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
# lux + noise
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + age + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + age + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + age + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + age + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + age + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + age + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + age + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + year + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + age + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + julian + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + year + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
cand.models[[y]] <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + measurer + (1|section) + (1|stage), data = webl, REML = FALSE); y <- y + 1
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

# Select top models manually (models with ΔAIC < 2)
top_models <- list(cand.models[[98]], cand.models[[49]], cand.models[[100]], cand.models[[81]], 
                   cand.models[[103]], cand.models[[39]])  
top_models

# Top model: bc_wing ~ scale(lux) + sex + (1 | section) + (1 | stage)
top_mod_1 <- lmer(bc_wing ~ scale(lux) + sex + (1 | section) + (1 | stage), data = webl)
summary(top_mod_1)
confint(top_mod_1)
confint(top_mod_1, level = 0.85)
# No effect of lux
# SIG effect of sex- M have higher BNC

# 2nd Top model: bc_wing ~ scale(lux) * scale(LAeq) + sex + age + (1 | section) + (1 | stage)
top_mod_2 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + (1 | section) + (1 | stage), data = webl)
summary(top_mod_2)
confint(top_mod_2)
confint(top_mod_2, level = 0.85)
# Lux has a SIG positive effect
# Noise has a pos effect at 85% CI
# Int has a SIG pos effect
# Sex have SIG effect - M higher
# Age has SIG effect with AHY and ASY having higher BC 

# Find means and SD of BC by Age
adj_means <- emmeans(top_mod_2, ~ age)
adj_means
pairs(adj_means)

# 3rd Top model: bc_wing ~ scale(lux) + age + (1 | section) + (1 | stage)
top_mod_3 <- lmer(bc_wing ~ scale(lux) + age + (1 | section) + (1 | stage), data = webl)
summary(top_mod_3)
confint(top_mod_3)
confint(top_mod_3, level = 0.85)
# Lux has a pos effect at 85% CI
# Age has SIG effect with AHY and ASY having higher BC 

# Find means and SD of BC by Age
adj_means <- emmeans(top_mod_3, ~ age)
adj_means
pairs(adj_means)

# 4th Top model: bc_wing ~ scale(lux) + scale(LAeq) + sex + age + (1 | section) + (1 | stage)
top_mod_4 <- lmer(bc_wing ~ scale(lux) + scale(LAeq) + sex + age + (1 | section) + (1 | stage), data = webl)
summary(top_mod_4)
confint(top_mod_4)
confint(top_mod_4, level = 0.85)
# No effect of lux or noise
# Males have higher BC

# Find means and SD of BC by Age
adj_means <- emmeans(top_mod_4, ~ age)
adj_means
pairs(adj_means)

# 5th Top model: bc_wing ~ scale(LAeq) + sex + (1 | section) + (1 | stage)
top_mod_5 <- lmer(bc_wing ~ scale(LAeq) + sex + (1 | section) + (1 | stage), data = webl)
summary(top_mod_5)
confint(top_mod_5)
confint(top_mod_5, level = 0.85)
# No effect of noise
# Males have higher BC

# 6th Top model: bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1 | section) + (1 | stage)
top_mod_6 <- lmer(bc_wing ~ scale(lux) * scale(LAeq) + sex + age + julian + (1 | section) + (1 | stage), data = webl)
summary(top_mod_6)
confint(top_mod_6)
confint(top_mod_6, level = 0.85)
# Lux has a positive effect at 85% CI
# Noise has a pos effect at 85% CI
# Int has a SIG pos effect
# Sex have SIG effect - M higher
# Age has SIG effect with AHY and ASY having higher BC 
# No effect of julian

# Find means and SD of BC by Age
adj_means <- emmeans(top_mod_6, ~ age)
adj_means
pairs(adj_means)





#### Normality Test
library(car)
qqPlot(residuals(top_mod_1))
shapiro.test(residuals(top_mod_1));length(residuals(top_mod_1)) # normally distributed
plot(density(resid(top_mod_1)))
# residuals are normal if p is > 0.05.


#### Homogeneity of variance 
# Performance Package to test for Heteroscedasticity
library(performance)
check_heteroscedasticity(top_mod_1) # model variance is homoscedastic

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




#### 2 categories W-BC across LAeq noise AND Lux values
webl <- webl %>%
  mutate(lux_binary = if_else(lux == 0, "No Lux", "Lux"),
         lux_binary = factor(lux_binary, levels = c("No Lux", "Lux")))

webl |> 
  filter(!is.na(lux_binary)) |> 
  ggplot(aes(x = lux_binary, y = lux)) +
  geom_point()
webl |> 
  filter(!is.na(lux_binary)) |> 
  ggplot(aes(x = LAeq, y = bc_wing, group = lux_binary, color = lux_binary)) +
  geom_point(size = 3) +  
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Lux" = "red", "No Lux" = "grey50")) +
  labs(title = "WEBL W-BC across Noise & Lux", x = "LAeq", y = "Wing Body Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



### 3 categories - W-BC across LAeq noise AND Lux values
webl <- webl %>%
  mutate(pos_lux_mean = mean(lux[lux > 0], na.rm = TRUE)) %>% # Calculate the mean of ONLY the positive values
  mutate(lux_cat3 = case_when(  # Assign categories based on that threshold
    lux == 0 ~ "No Lux",
    lux > 0 & lux <= pos_lux_mean ~ "Low Lux",
    lux > pos_lux_mean ~ "High Lux",
    TRUE ~ NA_character_ # Handles any existing NAs in the lux data
  )) %>%
  mutate(lux_cat3 = factor(lux_cat3, levels = c("No Lux", "Low Lux", "High Lux")))

webl |> 
  filter(!is.na(lux_cat3)) |> 
  ggplot(aes(x = lux_cat3, y = lux)) +
  geom_point()


webl %>%
  filter(!is.na(lux_cat3)) %>%
  ggplot(aes(x = LAeq, y = bc_wing, color = lux_cat3)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  scale_color_manual(values = c("No Lux" = "black", 
                                "Low Lux" = "blue", 
                                "High Lux" = "red")) +
  labs(
    title = "ATFL W-BC: Interaction of Noise and Light Intensity",
    x = "Noise Level (LAeq)",
    y = "Wing Body Condition",
    color = "Lux Category"
  ) +
  theme_minimal()







##################################

# Using specific values that match the 'scaled' distribution
inter_data <- predict_response(top_mod_2, terms = c("LAeq [all]", "lux [-1, 0, 1]"))

plot(inter_data)

# Convert the prediction to a data frame for more control
plot_df <- as.data.frame(inter_data)

# 1. Ensure your plot_df has the labels
plot_df$lux_label <- factor(plot_df$group, 
                            levels = c("-1", "0", "1"),
                            labels = c("Low Lux (-1 SD)", "Mean Lux", "High Lux (+1 SD)"))



library(ggnewscale)

ggplot() +
  geom_ribbon(data = plot_df, 
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = lux_label), 
              alpha = 0.12, show.legend = FALSE) + 
  scale_fill_manual(values = c("firebrick", "orange2","#FDDC5C")) +
  geom_line(data = plot_df, 
            aes(x = x, y = predicted, color = lux_label), 
            linewidth = 1.5) +
  scale_color_manual(values = c("firebrick", "orange2","#FDDC5C")) +
  new_scale_fill() +
  geom_point(data = webl, 
             aes(x = LAeq, y = bc_wing, fill = lux), 
             shape = 21, color = "grey60", size = 3, alpha = 0.9) +
  scale_fill_gradient(low = "red2", high = "yellow2") + 
  labs(
    title = "Noise & Lux Interacting on Wing-Body Condition",
    subtitle = "Lines = Model Predictions; Points = Raw Data colored by Lux",
    x = "Noise Level (LAeq)",
    y = "Body Condition",
    color = "Model 'Slices' (±1 SD)",
    fill = "Observed Lux Value"
  ) +
  theme_minimal() +
  theme(legend.position = "right")





