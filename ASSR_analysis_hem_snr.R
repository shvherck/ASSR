# ----------------------------------------------------------------- #
# name              ASSR_analysis_hem_snr.R                         #
# created by        Shauni Van Herck                                #
# description       used to analyse ASSR data (snr)                 #
# version.string    R version 3.5.1 (2018-07-02)                    #
# platform          x86_64-w64-mingw32                              #
# date created      20/01/2020                                      #
# ----------------------------------------------------------------- #


# load libraries -------------------
# reading
library(readxl)
# reshaping
library(stringr)      # split data into numbers/letters
library(reshape)
library(dplyr)
# plotting
library(car) #qqplot
library(ggplot2)
library(gridExtra)
# mixed models
library(lme4)
library(nlme)
# permutations
library(lmPerm)
# post-hoc comparisons
library(emmeans)

# prepare data ---------------------

setwd("C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/R")

d              = read.csv("ASSRdata.csv", header=TRUE, sep=",")

# only keep rows for electrodes 'Right' and 'Left'
# only keep rows for frequencies 4 & 20
d              = subset(d, Channel%in%c('Right', 'Left')) 
d              = subset(d, StimulationFrequency%in%c('4', '20')) 

#strsplit only works with character arguments
d$Recording    = as.character(d$Recording)

# create column subject and group
d$subject      = sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[1]))
#x[1]: keep the part before '_', 2: keep the part after '_'

groups         = read_excel("c1groups.xlsx")
d$group        = groups$Groep[match(unlist(d$subject), groups$I_code)]
d              = subset(d, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
#use groups dataframe to find correct group for each i-code

# create column stimulus type (AM/PULS)
d$stimtype     = (str_extract(sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[2])), "[aA-zZ]+"))
#strsplit to keep AM4, ... + str_extract to only keep letter part

# create column frequency
d$frequency    = d$StimulationFrequency

# create column hemisphere
d$hem          = d$Channel

# create column condition (stimtype + freq)
d$condition    = sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[2]))

# create column testing phase (pre/post)
d$prepost      = sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[3]))

# only keep relevant columns
d              = d[c("subject", "group", "condition", "stimtype", "frequency", "prepost", "hem", "AssrHt2BiasedRecordingSnrDb")]
colnames(d)   = c("subject", "group", "cond", "stimtype", "freq", "prepost", "hem", "snr")
head(d)

d$subject   = factor(d$subject)
d$group     = factor(d$group, levels=c("GG_EE", "GG_NE", "ActiveControl"))    #reorder the levels
d$cond      = factor(d$cond, levels=c("AM4", "AM20", "PULS4", "PULS20"))
d$stimtype  = factor(d$stimtype, levels=c("AM", "PULS"))
d$freq      = factor(d$freq, levels=c("4", "20"))
d$prepost   = factor(d$prepost, levels=c("pre", "post")) #reorder the levels
d$hem       = factor(d$hem)


d_4              = subset(d, freq%in%c('4'))
d_4$logsnr       = log(d_4$snr)

d_4$row          = 1:nrow(d_4)     #know row number of Cook's data points (see 'check assumptions')
d_4_cook         = d_4[ -c(168, 172, 467, 437, 468, 424, 10, 627, 47, 534, 518, 48, 163, 51, 501, 261, 423, 300, 101, 361, 214, 264, 547, 3), ]   #remove data points based on Cook's distance


d_20             = subset(d, freq%in%c('20')) 
d_20$logsnr      = log(d_20$snr)

d_20$row         = 1:nrow(d_20)                                
d_20_cook        = d_20[ -c(576, 262), ]   


d_AM4       = subset(d, cond%in%c('AM4'))
d_AM20      = subset(d, cond%in%c('AM20'))
d_PULS4       = subset(d, cond%in%c('PULS4'))
d_PULS20       = subset(d, cond%in%c('PULS20'))


# plot data ------------------

setwd("C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/R/Plots")


# AM4
boxplot1         = ggplot(d_AM4, aes(x=interaction(hem, group), y=snr, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                                    #makes sure that boxplots belonging to the same group are close to each other                                                                            
  labs(title="AM4",x="group : hemisphere", y = "snr (dB)", fill="test phase") + theme_classic() +  
  coord_cartesian(ylim=c(-10,25)) 
print(boxplot1)

row.names(d_AM4) <- as.character(paste(d_AM4$subject,"_", d_AM4$cond, "_" , d_AM4$prepost, "_", d_AM4$hem)) 
Boxplot(data=d_AM4, snr ~ prepost*hem*group)

# AM20
boxplot2         = ggplot(d_AM20, aes(x=interaction(hem,group), y=snr, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
  labs(title="AM20",x="group : hemisphere", y = "snr (dB)", fill="test phase") + theme_classic() +
  coord_cartesian(ylim=c(-10,25)) 
print(boxplot2)

row.names(d_AM20) <- as.character(paste(d_AM20$subject,"_", d_AM20$cond, "_" , d_AM20$prepost, "_", d_AM20$hem)) 
Boxplot(data=d_AM20, snr ~ prepost*hem*group)


# PULS4
boxplot3         = ggplot(d_PULS4, aes(x=interaction(hem, group), y=snr, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
  labs(title="PULS4",x="group : hemisphere", y = "snr (dB)", fill="test phase") + theme_classic() +
  coord_cartesian(ylim=c(-10,25)) 
print(boxplot3)

row.names(d_PULS4) <- as.character(paste(d_PULS4$subject,"_", d_PULS4$cond, "_" , d_PULS4$prepost, "_", d_PULS4$hem)) 
Boxplot(data=d_PULS4, snr ~ prepost*hem*group)

# PULS20
boxplot4         = ggplot(d_PULS20, aes(x=interaction(hem, group), y=snr, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
  labs(title="PULS20",x="group : hemisphere", y = "snr (dB)", fill="test phase") + theme_classic() +
  coord_cartesian(ylim=c(-10,25)) 
print(boxplot4)

row.names(d_PULS20) <- as.character(paste(d_PULS20$subject,"_", d_PULS20$cond, "_" , d_PULS20$prepost, "_", d_PULS20$hem)) 
Boxplot(data=d_PULS20, snr ~ prepost*hem*group)


tiff("snr_4Hz_hem_withnegativesnr.tiff",width=1600,height=800)                 #save figure in wd
print(grid.arrange(boxplot1, boxplot3, ncol = 2, top="4 Hz"))
dev.off()

tiff("snr_20Hz_hem_withnegativesnr.tiff",width=1600,height=800)                 #save figure in wd
print(grid.arrange(boxplot2, boxplot4, ncol = 2, top="20 Hz"))
dev.off()

# check assumptions general -------------------

#
# normality assumption #
#

## 4 Hz 
fitQQ4 <- lm(snr ~ prepost*stimtype*hem*group, data=d_4)
qqPlot(fitQQ4, main="QQ Plot")
plot(fitQQ4)
# Extract the residuals
fit_residuals4 <- residuals(object = fitQQ4)
# Shapiro-Wilk test -violated
shapiro.test(fit_residuals4)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), sd=sd(fit_residuals4))

# 4 Hz with data points based on Cook's distance removed (24 points removed)
fitQQ4cook <- lm(snr ~ prepost*stimtype*hem*group, data=d_4_cook)
qqPlot(fitQQ4cook, main="QQ Plot")
fit_residuals4cook <- residuals(object = fitQQ4cook)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals4cook, "pnorm", mean=mean(fit_residuals4cook), sd=sd(fit_residuals4cook))

## 4 Hz (done on log) -> doesn't work!
fitQQ4 <- lm(logsnr ~ prepost*stimtype*hem*group, data=d_4)
qqPlot(fitQQ4, main="QQ Plot")
plot(fitQQ4)
# Extract the residuals
fit_residuals4 <- residuals(object = fitQQ4)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), d=sd(fit_residuals4))


## 20 Hz
fitQQ20 <- lm(snr ~ prepost*stimtype*hem*group, data=d_20)
qqPlot(fitQQ20, main="QQ Plot")
# Extract the residuals
fit_residuals20 <- residuals(object = fitQQ20)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals20, "pnorm", mean=mean(fit_residuals20), sd=sd(fit_residuals20))

# 20 Hz with data points based on Cook's distance removed (2 points removed)
fitQQ20cook <- lm(snr ~ prepost*stimtype*group, data=d_20_cook)
qqPlot(fitQQ20cook, main="QQ Plot")
fit_residuals20cook <- residuals(object = fitQQ20cook)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals20cook, "pnorm", mean=mean(fit_residuals20cook), sd=sd(fit_residuals20cook))

## 20 Hz (done on log) -> doesn't work!
fitQQ20 <- lm(logsnr ~ prepost*stimtype*hem*group, data=d_20)
qqPlot(fitQQ20, main="QQ Plot")
# Extract the residuals
fit_residuals20 <- residuals(object = fitQQ20)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals20, "pnorm", mean=mean(fit_residuals20), sd=sd(fit_residuals20))

#
# sphericity assumption #
#

# this only needs to be tested when there are 3 levels or more of a repeated measures factor
# for now sphericity is not an issue yet
# it will be after consolidation test though!!!

#
# homogeneity of variances assumption #
#

## 4 Hz
leveneTest(snr ~ prepost*stimtype*group, data = d_4_cook)

plot(fitQQ4cook, 1)
# violated

## 20 Hz
leveneTest(snr ~ prepost*stimtype*group, data = d_20_cook)
plot(fitQQ20, 1)
# not violated


# analyses general -------------------

#
# 'modern' mixed model approach #
#

## 4 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit4             = lmer(snr  ~ group*stimtype*prepost*hem + (1|subject), data=d_4_cook)

summary(fit4)

Anova(fit4, type = "III", test.statistic = "F")


# stimtype: PULS > AM
# group:prepost: post-hoc emmeans: nothing significant
# prepost:hem: post-hoc emmeans: nothing significant
# stimtype:hem:
  # AM left < PULS left
  # AM left < PULS right
  # AM right < PULS left
  # AM right < PULS right
  # AM left < AM right

# post-hoc
ph41              = emmeans(fit4, specs = pairwise  ~ group:prepost, adjust = "holm")              
ph41$contrasts    # no significances even though there is a significant interaction

ph42              = emmeans(fit4, specs = pairwise  ~ stimtype:hem, adjust = "holm")              
ph42$contrasts

d_4_AML          = subset(d_4_cook, stimtype%in%c("AM") & hem%in%c("Left"))
d_4_AMR          = subset(d_4_cook, stimtype%in%c("AM") & hem%in%c("Right"))
mean(d_4_AML$snr)
mean(d_4_AMR$snr)
sd(d_4_AML$snr)
sd(d_4_AMR$snr)

d_4_PULSR        = subset(d_4_cook, stimtype%in%c("PULS") & hem%in%c("Right"))
d_4_PULSL          = subset(d_4_cook, stimtype%in%c("PULS") & hem%in%c("Left"))
mean(d_4_PULSR$snr)
mean(d_4_PULSL$snr)
sd(d_4_AML$snr)
sd(d_4_AMR$snr)

ph43              = emmeans(fit4, specs = pairwise  ~ prepost:hem, adjust = "holm")              
ph43$contrasts    # no significances even though there is a significant interaction

## model building try out

library(MuMIn)

fit4             = lmer(snr  ~ group*stimtype*prepost*hem + (1|subject), data=d_4_cook)
fit              = lmer(snr  ~ group*stimtype*prepost*hem + (prepost|subject), data=d_4_cook)

r.squaredGLMM(fit4) # this model has the best marginal R squared (or at least doesn't differ from fit)
r.squaredGLMM(fit)

# m marginal fixed effects - c conditional marginal effects (usually interested in marginal)


# 20 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit20             = lmer(snr  ~ group*stimtype*prepost*hem + (1|subject), data=d_20_cook)

summary(fit20)

Anova(fit20, type = "III", test.statistic = "F")

# stimtype: AM > PULS
# prepost: POST > PRE
# hem: right > left
# group:hem
  # GGEE left < GGEE right


# post-hoc
d_20_AM              = subset(d_20_cook, stimtype%in%c('AM'))
mean(d_20_AM$snr)
sd(d_20_AM$snr)
d_20_PULS            = subset(d_20_cook, stimtype%in%c('PULS'))
mean(d_20_PULS$snr)
sd(d_20_PULS$snr)

d_20_pre             = subset(d_20_cook, prepost%in%c('pre'))
mean(d_20_pre$snr)
sd(d_20_pre$snr)
d_20_post            = subset(d_20_cook, prepost%in%c('post'))
mean(d_20_post$snr)
sd(d_20_post$snr)

d_20_left            = subset(d_20_cook, hem%in%c('Left'))
mean(d_20_left$snr)
sd(d_20_left$snr)
d_20_right            = subset(d_20_cook, hem%in%c('Right'))
mean(d_20_right$snr)
sd(d_20_right$snr)


ph20            = emmeans(fit20, specs = pairwise  ~ group:hem, adjust = "holm")              
ph20$contrasts 

d_20_GGEEL      = subset(d_20_cook, group%in%c("GG_EE") & hem%in%c("Left"))
mean(d_20_GGEEL$snr)
sd(d_20_GGEEL$snr)
d_20_GGEER      = subset(d_20_cook, group%in%c("GG_EE") & hem%in%c("Right"))
mean(d_20_GGEER$snr)
sd(d_20_GGEER$snr)

