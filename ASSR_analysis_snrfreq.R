# ----------------------------------------------------------------- #
# name              ASSR_analysis_snrfreq.R                         #
# created by        Shauni Van Herck                                #
# description       used to analyse ASSR data (snr)                 #
# version.string    R version 3.5.1 (2018-07-02)                    #
# platform          x86_64-w64-mingw32                              #
# date created      07/10/2019                                      #
# ----------------------------------------------------------------- #


# load libraries -------------------
# reading
library(readxl)
# reshaping
library(stringr)      # split data into numbers/letters
library(reshape)
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

ASSRdata              = read.csv("ASSRdata.csv", header=TRUE, sep=",")

# only keep rows for electrodes 'Right', 'Left' and 'All'
# only keep rows for frequencies 4 & 20
ASSRdata              = subset(ASSRdata, Channel%in%c('Right', 'Left', 'All')) 
ASSRdata              = subset(ASSRdata, StimulationFrequency%in%c('4', '20')) 

#strsplit only works with character arguments
ASSRdata$Recording    = as.character(ASSRdata$Recording)

# create column subject and group
ASSRdata$subject      = sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[1]))
  #x[1]: keep the part before '_', 2: keep the part after '_'

groups                = read_excel("c1groups.xlsx")
ASSRdata$group        = groups$Groep[match(unlist(ASSRdata$subject), groups$I_code)]
  #use groups dataframe to find correct group for each i-code

# create column stimulus type (AM/PULS)
ASSRdata$stimtype     = (str_extract(sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[2])), "[aA-zZ]+"))
  #strsplit to keep AM4, ... + str_extract to only keep letter part

# create column frequency
ASSRdata$frequency    = ASSRdata$StimulationFrequency

# create column condition (stimtype + freq)
ASSRdata$condition    = sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[2]))

# create column testing phase (pre/post)
ASSRdata$prepost      = sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[3]))

# only keep relevant columns
ASSRdata              = ASSRdata[c("subject", "group", "condition", "stimtype", "frequency", "prepost", "Channel", "AssrHt2BiasedRecordingSnrDb")]

# subset data
ASSRdata_all          = subset(ASSRdata, Channel%in%c('All'))
colnames(ASSRdata_all)= c("subject", "group", "cond", "stimtype", "freq", "prepost", "channel", "snr")
head(ASSRdata_all)

ASSRdata_all$subject   = factor(ASSRdata_all$subject)
ASSRdata_all$group     = factor(ASSRdata_all$group, levels=c("GG_EE", "GG_NE", "ActiveControl", "PassiveControl"))    #reorder the levels
ASSRdata_all$cond      = factor(ASSRdata_all$cond, levels=c("AM4", "AM20", "PULS4", "PULS20"))
ASSRdata_all$stimtype  = factor(ASSRdata_all$stimtype, levels=c("AM", "PULS"))
ASSRdata_all$freq      = factor(ASSRdata_all$freq, levels=c("4", "20"))
ASSRdata_all$prepost   = factor(ASSRdata_all$prepost, levels=c("pre", "post")) #reorder the levels
ASSRdata_all$channel   = factor(ASSRdata_all$channel)


ASSRdata_all_4              = subset(ASSRdata_all, freq%in%c('4'))
ASSRdata_all_4              = subset(ASSRdata_all_4, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_4$row          = 1:nrow(ASSRdata_all_4)                                 #know row number of Cook's data points (see 'check assumptions')
ASSRdata_all_4_cook         = ASSRdata_all_4[ -c(219, 234), ]                   #remove data points based on Cook's distance
ASSRdata_all_4_final_cook   = subset(ASSRdata_all_4_cook, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_4$logsnr       = log(ASSRdata_all_4$snr)
ASSRdata_all_4_final        = ASSRdata_all_4
ASSRdata_all_20             = subset(ASSRdata_all, freq%in%c('20'))                  
ASSRdata_all_20             = subset(ASSRdata_all_20, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_20$row         = 1:nrow(ASSRdata_all_20)                                
ASSRdata_all_20_cook        = ASSRdata_all_20[ -c(131, 247, 48), ]
ASSRdata_all_20_final       = subset(ASSRdata_all_20, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_20_final_cook  = subset(ASSRdata_all_20_cook, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))


ASSRdata_all_AM4       = subset(ASSRdata_all, cond%in%c('AM4'))
ASSRdata_all_AM4_final = subset(ASSRdata_all_AM4, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_AM4_GGEE  = subset(ASSRdata_all_AM4, group%in%c('GG_EE'))
ASSRdata_all_AM4_GGNE  = subset(ASSRdata_all_AM4, group%in%c('GG_NE'))
ASSRdata_all_AM4_AC    = subset(ASSRdata_all_AM4, group%in%c('ActiveControl'))

ASSRdata_all_AM20      = subset(ASSRdata_all, cond%in%c('AM20'))
ASSRdata_all_AM20_final= subset(ASSRdata_all_AM20, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_AM20_GGEE = subset(ASSRdata_all_AM20, group%in%c('GG_EE'))
ASSRdata_all_AM20_GGNE = subset(ASSRdata_all_AM20, group%in%c('GG_NE'))
ASSRdata_all_AM20_AC   = subset(ASSRdata_all_AM20, group%in%c('ActiveControl'))

ASSRdata_all_PULS4       = subset(ASSRdata_all, cond%in%c('PULS4'))
ASSRdata_all_PULS4_final = subset(ASSRdata_all_PULS4, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_PULS4_GGEE  = subset(ASSRdata_all_PULS4, group%in%c('GG_EE'))
ASSRdata_all_PULS4_GGNE  = subset(ASSRdata_all_PULS4, group%in%c('GG_NE'))
ASSRdata_all_PULS4_AC    = subset(ASSRdata_all_PULS4, group%in%c('ActiveControl'))

ASSRdata_all_PULS20       = subset(ASSRdata_all, cond%in%c('PULS20'))
ASSRdata_all_PULS20_final = subset(ASSRdata_all_PULS20, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_PULS20_GGEE  = subset(ASSRdata_all_PULS20, group%in%c('GG_EE'))
ASSRdata_all_PULS20_GGNE  = subset(ASSRdata_all_PULS20, group%in%c('GG_NE'))
ASSRdata_all_PULS20_AC    = subset(ASSRdata_all_PULS20, group%in%c('ActiveControl'))

ASSRdata_all_AM4_melt          = melt(ASSRdata_all_AM4, id.var=c("subject", "group", "prepost"), measure.var=c("snr"))
ASSRdata_all_AM4_wide          = cast(subject + group ~ prepost , data=ASSRdata_all_AM4_melt, mean)
ASSRdata_all_AM4_wide          = ASSRdata_all_AM4_wide[c("subject", "group", "pre", "post")]
head(ASSRdata_all_AM4_wide)

ASSRdata_all_AM20_melt          = melt(ASSRdata_all_AM20, id.var=c("subject", "group", "prepost"), measure.var=c("snr"))
ASSRdata_all_AM20_wide          = cast(subject + group ~ prepost , data=ASSRdata_all_AM20_melt, mean)
ASSRdata_all_AM20_wide          = ASSRdata_all_AM20_wide[c("subject", "group", "pre", "post")]
head(ASSRdata_all_AM20_wide)

ASSRdata_all_PULS4_melt          = melt(ASSRdata_all_PULS4, id.var=c("subject", "group", "prepost"), measure.var=c("snr"))
ASSRdata_all_PULS4_wide          = cast(subject + group ~ prepost , data=ASSRdata_all_PULS4_melt, mean)
ASSRdata_all_PULS4_wide          = ASSRdata_all_PULS4_wide[c("subject", "group", "pre", "post")]
head(ASSRdata_all_PULS4_wide)

ASSRdata_all_PULS20_melt          = melt(ASSRdata_all_PULS20, id.var=c("subject", "group", "prepost"), measure.var=c("snr"))
ASSRdata_all_PULS20_wide          = cast(subject + group ~ prepost , data=ASSRdata_all_PULS20_melt, mean)
ASSRdata_all_PULS20_wide          = ASSRdata_all_PULS20_wide[c("subject", "group", "pre", "post")]
head(ASSRdata_all_PULS20_wide)

# plot data ------------------

setwd("C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/R/Plots")


# AM4
boxplot1         = ggplot(ASSRdata_all_AM4, aes(x=group, y=snr, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        # position dodge + geom_jitter(position_dodge) is for
#  geom_jitter(position=position_dodge(0.8)) +                                        # getting dots with multiple groups (so that dots for)
  labs(title="AM4",x="group", y = "snr (dB)", fill="test phase") + theme_classic() +  # pre and post are also separated per group
  coord_cartesian(ylim=c(-10,25))
# print(boxplot1)

# AM20
boxplot2         = ggplot(ASSRdata_all_AM20, aes(x=group, y=snr, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
#  geom_jitter(position=position_dodge(0.8)) +                                         
  labs(title="AM20",x="group", y = "snr (dB)", fill="test phase") + theme_classic() +
  coord_cartesian(ylim=c(-10,25))
#print(boxplot2)

# PULS4
boxplot3         = ggplot(ASSRdata_all_PULS4, aes(x=group, y=snr, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
#  geom_jitter(position=position_dodge(0.8)) +                                         
  labs(title="PULS4",x="group", y = "snr (dB)", fill="test phase") + theme_classic() +
  coord_cartesian(ylim=c(-10,25))
#print(boxplot3)

# PULS20
boxplot4         = ggplot(ASSRdata_all_PULS20, aes(x=group, y=snr, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
#  geom_jitter(position=position_dodge(0.8)) +                                         
  labs(title="PULS20",x="group", y = "snr (dB)", fill="test phase") + theme_classic() +
  coord_cartesian(ylim=c(-10,25))
#print(boxplot4)


tiff("snrfreq_4Hz.tiff",width=1360,height=1360)                 #save figure in wd
print(grid.arrange(boxplot1, boxplot3, ncol = 2, top="4 Hz"))
dev.off()

tiff("snrfreq_20Hz.tiff",width=1360,height=1360)                 #save figure in wd
print(grid.arrange(boxplot2, boxplot4, ncol = 2, top="20 Hz"))
dev.off()


# spaghettiplot AM4
ASSRdata_all_AM4_1 = subset(ASSRdata_all_AM4, group=="GG_EE")
ASSRdata_all_AM4_2 = subset(ASSRdata_all_AM4, group=="GG_NE")
ASSRdata_all_AM4_3 = subset(ASSRdata_all_AM4, group=="ActiveControl")

p1            = ggplot(data = ASSRdata_all_AM4_1, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_EE", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")
p2            = ggplot(data = ASSRdata_all_AM4_2, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_NE", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")
p3            = ggplot(data = ASSRdata_all_AM4_3, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="Active control", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")

tiff("snrfreq_AM4_spaghetti_preliminary.tiff",width=1360,height=1360)
print(grid.arrange(p1, p2, p3, ncol = 2, top="Spaghetti plot for AM4"))
dev.off()

# spaghettiplot AM20
ASSRdata_all_AM20_1 = subset(ASSRdata_all_AM20, group=="GG_EE")
ASSRdata_all_AM20_2 = subset(ASSRdata_all_AM20, group=="GG_NE")
ASSRdata_all_AM20_3 = subset(ASSRdata_all_AM20, group=="ActiveControl")

p1            = ggplot(data = ASSRdata_all_AM20_1, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_EE", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")
p2            = ggplot(data = ASSRdata_all_AM20_2, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_NE", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")
p3            = ggplot(data = ASSRdata_all_AM20_3, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="Active control", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")

tiff("snrfreq_AM20_spaghetti_preliminary.tiff",width=1360,height=1360)
print(grid.arrange(p1, p2, p3, ncol = 2, top="Spaghetti plot for AM20"))
dev.off()

# spaghettiplot PULS4
ASSRdata_all_PULS4_1 = subset(ASSRdata_all_PULS4, group=="GG_EE")
ASSRdata_all_PULS4_2 = subset(ASSRdata_all_PULS4, group=="GG_NE")
ASSRdata_all_PULS4_3 = subset(ASSRdata_all_PULS4, group=="ActiveControl")

p1            = ggplot(data = ASSRdata_all_PULS4_1, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_EE", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")
p2            = ggplot(data = ASSRdata_all_PULS4_2, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_NE", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")
p3            = ggplot(data = ASSRdata_all_PULS4_3, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="Active control", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")

tiff("snrfreq_PULS4_spaghetti_preliminary.tiff",width=1360,height=1360)
print(grid.arrange(p1, p2, p3, ncol = 2, top="Spaghetti plot for PULS4"))
dev.off()

# spaghettiplot PULS20
ASSRdata_all_PULS20_1 = subset(ASSRdata_all_PULS20, group=="GG_EE")
ASSRdata_all_PULS20_2 = subset(ASSRdata_all_PULS20, group=="GG_NE")
ASSRdata_all_PULS20_3 = subset(ASSRdata_all_PULS20, group=="ActiveControl")

p1            = ggplot(data = ASSRdata_all_PULS20_1, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_EE", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")
p2            = ggplot(data = ASSRdata_all_PULS20_2, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_NE", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")
p3            = ggplot(data = ASSRdata_all_PULS20_3, aes(x = prepost, y = snr, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="Active control", x="test phase", y="snr (dB)") + theme_bw() +
  theme(legend.position="none")

tiff("snrfreq_PULS20_spaghetti_preliminary.tiff",width=1360,height=1360)
print(grid.arrange(p1, p2, p3, ncol = 2, top="Spaghetti plot for PULS20"))
dev.off()


# descriptives --------------

# AM4
descriptives            = cast(ASSRdata_all_AM4_melt, group*prepost ~ ., fun.aggregate=c(median,IQR))
names(descriptives)[3]  = "median"
names(descriptives)[4]  = "IQR"
descriptives
AM4descr                = as.data.frame(descriptives)

# AM20
descriptives            = cast(ASSRdata_all_AM20_melt, group*prepost ~ ., fun.aggregate=c(median,IQR))
names(descriptives)[3]  = "median"
names(descriptives)[4]  = "IQR"
descriptives
AM20descr               = as.data.frame(descriptives)

# PULS4
descriptives            = cast(ASSRdata_all_PULS4_melt, group*prepost ~ ., fun.aggregate=c(median,IQR))
names(descriptives)[3]  = "median"
names(descriptives)[4]  = "IQR"
descriptives
PULS4descr              = as.data.frame(descriptives)


# PULS20
descriptives            = cast(ASSRdata_all_PULS20_melt, group*prepost ~ ., fun.aggregate=c(median,IQR))
names(descriptives)[3]  = "median"
names(descriptives)[4]  = "IQR"
descriptives
PULS20descr             = as.data.frame(descriptives)








# check assumptions general -------------------

#
# normality assumption #
#

## 4 Hz 
fitQQ4 <- lm(snr ~ prepost*stimtype*group, data=ASSRdata_all_4)
qqPlot(fitQQ4, main="QQ Plot")
plot(fitQQ4)
library(MASS)
sresid <- studres(fitQQ4) 
hist(sresid, freq=FALSE, main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40) 
yfit<-dnorm(xfit) 
lines(xfit, yfit)
# Extract the residuals
fit_residuals4 <- residuals(object = fitQQ4)
# Shapiro-Wilk test -violated
shapiro.test(fit_residuals4)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), sd=sd(fit_residuals4))

# 4 Hz with data points based on Cook's distance removed (2 points removed)
fitQQ4cook <- lm(snr ~ prepost*stimtype*group, data=ASSRdata_all_4_cook)
qqPlot(fitQQ4cook, main="QQ Plot")
plot(fitQQ4cook)
fit_residuals4cook <- residuals(object = fitQQ4cook)
# Shapiro-Wilk test -violated
shapiro.test(fit_residuals4cook)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals4cook, "pnorm", mean=mean(fit_residuals4cook), sd=sd(fit_residuals4cook))

## 4 Hz (done on log) 
fitQQ4 <- lm(logsnr ~ prepost*stimtype*group, data=ASSRdata_all_4)
qqPlot(fitQQ4, main="QQ Plot")
plot(fitQQ4)
# Extract the residuals
fit_residuals4 <- residuals(object = fitQQ4)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), sd=sd(fit_residuals4))


## 20 Hz
fitQQ20 <- lm(snr ~ prepost*stimtype*group, data=ASSRdata_all_20)
qqPlot(fitQQ20, main="QQ Plot")
plot(fitQQ20)
library(MASS)
sresid <- studres(fitQQ20) 
hist(sresid, freq=FALSE, main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40) 
yfit<-dnorm(xfit) 
lines(xfit, yfit)
# Extract the residuals
fit_residuals20 <- residuals(object = fitQQ20)
# Shapiro-Wilk test - violated
shapiro.test(x = fit_residuals20)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals20, "pnorm", mean=mean(fit_residuals20), sd=sd(fit_residuals20))


# 20 Hz with data points based on Cook's distance removed
fitQQ20cook <- lm(snr ~ prepost*stimtype*group, data=ASSRdata_all_20_cook)
qqPlot(fitQQ20cook, main="QQ Plot")
plot(fitQQ20cook)
fit_residuals20cook <- residuals(object = fitQQ20cook)
# Shapiro-Wilk test -violated
shapiro.test(fit_residuals20cook)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals20cook, "pnorm", mean=mean(fit_residuals20cook), sd=sd(fit_residuals20cook))


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
leveneTest(snr ~ prepost*stimtype*group, data = ASSRdata_all_4_cook)
plot(fitQQ4cook, 1)
# violated

## 20 Hz
leveneTest(snr ~ prepost*stimtype*group, data = ASSRdata_all_20_cook)
leveneTest(snr ~ prepost*stimtype*group, data = ASSRdata_all_20)
plot(fitQQ20, 1)
# not violated

# NOT DONE YET
# check equality of distributions for permutation testing
sm.density.compare(RT_post$threshold, RT_post$fgroup)
# distributions appear to be quite equal


# analyses general -------------------

#
# 'modern' mixed model approach #
#

## 4 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit4             = lmer(snr  ~ group*stimtype*prepost + (1|subject), data=ASSRdata_all_4_final_cook)
summary(fit4)

Anova(fit4, type = "III", test.statistic = "F")

# post-hoc
AM4              = subset(ASSRdata_all_4_final_cook, stimtype%in%c('AM'))
mean(AM4$snr)
sd(AM4$snr)

PULS4            = subset(ASSRdata_all_4_final_cook, stimtype%in%c('PULS'))
mean(PULS4$snr)
sd(PULS4$snr)

# 20 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit20             = lmer(snr  ~ group*stimtype*prepost + (1|subject), data=ASSRdata_all_20_final)
summary(fit20)

Anova(fit20, type = "III", test.statistic = "F")

# post-hoc
AM20              = subset(ASSRdata_all_20_final, stimtype%in%c('AM'))
mean(AM20$snr)
sd(AM20$snr)

PULS20            = subset(ASSRdata_all_20_final, stimtype%in%c('PULS'))
mean(PULS20$snr)
sd(PULS20$snr)

PRE20             = subset(ASSRdata_all_20_final, prepost%in%c('pre'))
mean(PRE20$snr)
sd(PRE20$snr)

POST20            = subset(ASSRdata_all_20_final, prepost%in%c('post'))
mean(POST20$snr)
sd(POST20$snr)



options(contrasts=c("contr.sum", "contr.poly"))
fit20cook             = lmer(snr  ~ group*stimtype*prepost + (1|subject) + (1|stimtype:subject) + (1|prepost:subject), data=ASSRdata_all_20_final_cook)
summary(fit20cook)

Anova(fit20cook, type = "III", test.statistic = "F")


#
# nlme package #
#

## 4 Hz
lme4  <-lme(snr ~ group*stimtype*prepost, random =~1|subject, data=ASSRdata_all_4_final_cook, method="ML")
summary(lme4)

lme4  <-lme(snr ~ group*stimtype*prepost, random =~1|subject, data=ASSRdata_all_4_final_cook, method="ML", weights=varIdent(form=~1|group)) # the weights factor accounts for different variances per group
summary(lme4)


# 
# non-parametric approach using aovp #
# 

## 4 Hz
perm4 = aovp(snr ~ group*stimtype*prepost + Error(subject), data=ASSRdata_all_4_final_cook, perm="Exact")
summary(perm4)

## 20 Hz
perm20 = aovp(snr ~ group*stimtype*prepost + Error(subject), data=ASSRdata_all_20_final, perm="Exact")
summary(perm20)

# check assumptions AM -----------------

#
# normality assumption #
#

## 4 Hz 
fitQQ4 <- lm(snr ~ prepost*group, data=ASSRdata_all_AM4_final)
qqPlot(fitQQ4, main="QQ Plot")
plot(fitQQ4)
# Extract the residuals
fit_residuals4 <- residuals(object = fitQQ4)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), sd=sd(fit_residuals4))


## 20 Hz
fitQQ20 <- lm(snr ~ prepost*group, data=ASSRdata_all_AM20_final)
qqPlot(fitQQ20, main="QQ Plot")
plot(fitQQ20)
# Extract the residuals
fit_residuals20 <- residuals(object = fitQQ20)
# Kolmogorov-Smirnov test - not violated
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
leveneTest(snr ~ prepost*group, data = ASSRdata_all_AM4_final)
plot(fitQQ4, 1)
# not violated

## 20 Hz
leveneTest(snr ~ prepost*group, data = ASSRdata_all_AM20_final)
plot(fitQQ20,1)
# not violated

# analyses AM ------------------


## 4 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit4             = lmer(snr  ~ group*prepost + (1|subject), data=ASSRdata_all_AM4_final)
summary(fit4)

Anova(fit4, type = "III", test.statistic = "F")

## 20 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit20             = lmer(snr  ~ group*prepost + (1|subject), data=ASSRdata_all_AM20_final)
summary(fit20)

Anova(fit20, type = "III", test.statistic = "F")




# check assumptions PULS ---------------

#
# normality assumption #
#

## 4 Hz 
fitQQ4 <- lm(snr ~ prepost*group, data=ASSRdata_all_PULS4_final)
qqPlot(fitQQ4, main="QQ Plot")
plot(fitQQ4)
# Extract the residuals
fit_residuals4 <- residuals(object = fitQQ4)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), sd=sd(fit_residuals4))


## 20 Hz
fitQQ20 <- lm(snr ~ prepost*group, data=ASSRdata_all_PULS20_final)
qqPlot(fitQQ20, main="QQ Plot")
plot(fitQQ20)
# Extract the residuals
fit_residuals20 <- residuals(object = fitQQ20)
# Kolmogorov-Smirnov test - not violated
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
leveneTest(snr ~ prepost*group, data = ASSRdata_all_PULS4_final)
plot(fitQQ4, 1)
# not violated

## 20 Hz
leveneTest(snr ~ prepost*group, data = ASSRdata_all_PULS20_final)
plot(fitQQ20,1)
# not violated



# analyses PULS --------------

## 4 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit4             = lmer(snr  ~ group*prepost + (1|subject), data=ASSRdata_all_PULS4_final)
summary(fit4)

Anova(fit4, type = "III", test.statistic = "F")

## 20 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit20             = lmer(snr  ~ group*prepost + (1|subject), data=ASSRdata_all_PULS20_final)
summary(fit20)

Anova(fit20, type = "III", test.statistic = "F")

# post-hoc
ph20              = emmeans(fit20, specs = pairwise  ~ group:prepost, adjust = "holm")              
ph20$contrasts




