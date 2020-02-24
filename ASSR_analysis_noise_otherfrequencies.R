# --------------------------------------------------------- #
# name              ASSR_analysis_noise_otherfrequencies.R  #
# created by        Shauni Van Herck                        #
# description       used to analyse ASSR data (noise amp)   #
# version.string    R version 3.5.1 (2018-07-02)            #
# platform          x86_64-w64-mingw32                      #
# date created      28/01/2020                              #
# --------------------------------------------------------- #


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

setwd(dir             = "C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/R")

d              = read.csv("ASSRdata.csv", header=TRUE, sep=",")

# only keep rows for electrodes 'Right', 'Left' and 'All'
# only keep rows for frequencies 4 & 20
d              = subset(d, Channel%in%c('Right', 'Left', 'All')) 
d              = subset(d, StimulationFrequency%in%c('3', '5', '19', '21')) 

#strsplit only works with character arguments
d$Recording    = as.character(d$Recording)

# create column subject and group
d$subject      = sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[1]))
#x[1]: keep the part before '_', 2: keep the part after '_'

groups                = read_excel("c1groups.xlsx")
d$group        = groups$Groep[match(unlist(d$subject), groups$I_code)]
#use groups dataframe to find correct group for each i-code

# create column stimulus type (AM/PULS)
d$stimtype    = (str_extract(sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[2])), "[aA-zZ]+"))
#strsplit to keep AM4, ... + str_extract to only keep letter part

# create column frequency
d$frequency    = d$StimulationFrequency

# create column condition (stimtype + freq)
d$condition    = sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[2]))

# create column testing phase (pre/post)
d$prepost      = sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[3]))

# only keep relevant columns
d              = d[c("subject", "group", "condition", "stimtype", "frequency", "prepost", "Channel", "AssrHt2RecordingNoiseAmplitude")]

# subset data
d_all          = subset(d, Channel%in%c('All'))
colnames(d_all)= c("subject", "group", "cond", "stimtype", "freq", "prepost", "channel", "noiseamp")
head(d_all)

d_all$subject  = factor(d_all$subject)
d_all$group    = factor(d_all$group, levels=c("GG_EE", "GG_NE", "ActiveControl", "PassiveControl"))    #reorder the levels
d_all$cond     = factor(d_all$cond, levels=c("AM4", "AM20", "PULS4", "PULS20"))
d_all$stimtype = factor(d_all$stimtype, levels=c("AM", "PULS"))
d_all$freq     = factor(d_all$freq, levels=c("3", "5", "19", "21"))
d_all$prepost  = factor(d_all$prepost, levels=c("pre", "post")) #reorder the levels
d_all$channel  = factor(d_all$channel)

d_all_35              = subset(d_all, freq%in%c('3', '5'))  
d_all_35              = subset(d_all_35, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
d_all_35$lognoiseamp  = log(d_all_35$noiseamp)

d_all_1921             = subset(d_all, freq%in%c('19', '21')) 
d_all_1921             = subset(d_all_1921, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
d_all_1921$lognoiseamp  = log(d_all_1921$noiseamp)

d_all_AM4                      = subset(d_all, cond%in%c('AM4'))
d_all_AM20                     = subset(d_all, cond%in%c('AM20'))
d_all_PULS4                    = subset(d_all, cond%in%c('PULS4'))
d_all_PULS20                   = subset(d_all, cond%in%c('PULS20'))


d_all_AM4_melt          = melt(d_all_AM4, id.var=c("subject", "group", "prepost"), measure.var=c("noiseamp"))
d_all_AM4_wide          = cast(subject + group ~ prepost , data=d_all_AM4_melt, mean)
d_all_AM4_wide          = d_all_AM4_wide[c("subject", "group", "pre", "post")]
head(d_all_AM4_wide)

d_all_AM20_melt          = melt(d_all_AM20, id.var=c("subject", "group", "prepost"), measure.var=c("noiseamp"))
d_all_AM20_wide          = cast(subject + group ~ prepost , data=d_all_AM20_melt, mean)
d_all_AM20_wide          = d_all_AM20_wide[c("subject", "group", "pre", "post")]
head(d_all_AM20_wide)

d_all_PULS4_melt          = melt(d_all_PULS4, id.var=c("subject", "group", "prepost"), measure.var=c("noiseamp"))
d_all_PULS4_wide          = cast(subject + group ~ prepost , data=d_all_PULS4_melt, mean)
d_all_PULS4_wide          = d_all_PULS4_wide[c("subject", "group", "pre", "post")]
head(d_all_PULS4_wide)

d_all_PULS20_melt          = melt(d_all_PULS20, id.var=c("subject", "group", "prepost"), measure.var=c("noiseamp"))
d_all_PULS20_wide          = cast(subject + group ~ prepost , data=d_all_PULS20_melt, mean)
d_all_PULS20_wide          = d_all_PULS20_wide[c("subject", "group", "pre", "post")]
head(d_all_PULS20_wide)

# plot data ------------------

setwd("C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/R/Plots")


# AM4
boxplot1         = ggplot(d_all_AM4, aes(x=group, y=noiseamp, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                                        # #makes sure that boxplots belonging to the same group are close to each other
  #  geom_jitter(position=position_dodge(0.8)) +                                                      
  labs(title="AM4",x="group", y = "noise amplitude (µV)", fill="test phase") + theme_classic()   
#print(boxplot1)

row.names(d_all_AM4) <- as.character(paste(d_all_AM4$subject,"_", d_all_AM4$cond, "_" ,d_all_AM4$prepost)) 
Boxplot(data=d_all_AM4, noiseamp ~ prepost*group)

# AM20
boxplot2         = ggplot(d_all_AM20, aes(x=group, y=noiseamp, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
  #  geom_jitter(position=position_dodge(0.8)) +                                         
  labs(title="AM20",x="group", y = "noise amplitude (µV)", fill="test phase") + theme_classic() 
#print(boxplot2)

row.names(d_all_AM20) <- as.character(paste(d_all_AM20$subject,"_", d_all_AM20$cond, "_" ,d_all_AM20$prepost)) 
Boxplot(data=d_all_AM20, noiseamp ~ prepost*group)


# PULS4
boxplot3         = ggplot(d_all_PULS4, aes(x=group, y=noiseamp, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
  #  geom_jitter(position=position_dodge(0.8)) +                                         
  labs(title="PULS4",x="group", y = "noise amplitude (µV)", fill="test phase") + theme_classic() 
#print(boxplot3)

row.names(d_all_PULS4) <- as.character(paste(d_all_PULS4$subject,"_", d_all_PULS4$cond, "_" ,d_all_PULS4$prepost)) 
Boxplot(data=d_all_PULS4, noiseamp ~ prepost*group)

# PULS20
boxplot4         = ggplot(d_all_PULS20, aes(x=group, y=noiseamp, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
  #  geom_jitter(position=position_dodge(0.8)) +                                         
  labs(title="PULS20",x="group", y = "noise amplitude (µV)p", fill="test phase") + theme_classic() 
#print(boxplot4)

row.names(d_all_PULS20) <- as.character(paste(d_all_PULS20$subject,"_", d_all_PULS20$cond, "_" ,d_all_PULS20$prepost)) 
Boxplot(data=d_all_PULS20, noiseamp ~ prepost*group)


tiff("noiseamp_otherfrequencies_box.tiff",width=1360,height=1360)                #save figure in wd
print(grid.arrange(boxplot1, boxplot2, boxplot3, boxplot4, ncol = 2, top="Boxplots"))
dev.off()


# descriptives --------------

# AM4
descriptives            = cast(d_all_AM4_melt, group*prepost ~ ., fun.aggregate=c(median,IQR))
names(descriptives)[3]  = "median"
names(descriptives)[4]  = "IQR"
descriptives
AM4descr                = as.data.frame(descriptives)

# AM20
descriptives            = cast(d_all_AM20_melt, group*prepost ~ ., fun.aggregate=c(median,IQR))
names(descriptives)[3]  = "median"
names(descriptives)[4]  = "IQR"
descriptives
AM20descr               = as.data.frame(descriptives)

# PULS4
descriptives            = cast(d_all_PULS4_melt, group*prepost ~ ., fun.aggregate=c(median,IQR))
names(descriptives)[3]  = "median"
names(descriptives)[4]  = "IQR"
descriptives
PULS4descr              = as.data.frame(descriptives)


# PULS20
descriptives            = cast(d_all_PULS20_melt, group*prepost ~ ., fun.aggregate=c(median,IQR))
names(descriptives)[3]  = "median"
names(descriptives)[4]  = "IQR"
descriptives
PULS20descr             = as.data.frame(descriptives)


# check assumptions general -------------------


#
# normality assumption #
#

## 4 (3-5) Hz 
fitQQ4 <- lm(noiseamp ~ prepost*stimtype*group, data=d_all_35)
qqPlot(fitQQ4, main="QQ Plot")
plot(fitQQ4)
# Extract the residuals
fit_residuals4 <- residuals(object = fitQQ4)
# Shapiro-Wilk test -violated
shapiro.test(fit_residuals4)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), sd=sd(fit_residuals4))

# 4 (3-5) Hz with log transformed data
fitQQ4log <- lm(lognoiseamp ~ prepost*stimtype*group, data=d_all_35)
qqPlot(fitQQ4log , main="QQ Plot")
plot(fitQQ4log)
# Extract the residuals
fit_residuals4log <- residuals(object = fitQQ4log)
# Shapiro-Wilk test - violated
shapiro.test(x = fit_residuals4log)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals4log, "pnorm", mean=mean(fit_residuals4log), sd=sd(fit_residuals4log))


## 20 (19-21) Hz
fitQQ20 <- lm(noiseamp ~ prepost*stimtype*group, data=d_all_1921)
qqPlot(fitQQ20, main="QQ Plot")
plot(fitQQ20)
# Extract the residuals
fit_residuals20 <- residuals(object = fitQQ20)
# Shapiro-Wilk test - violated
shapiro.test(x = fit_residuals20)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals20, "pnorm", mean=mean(fit_residuals20), sd=sd(fit_residuals20))


# 20 (19-21) Hz with log transformed data
fitQQ20log <- lm(lognoiseamp ~ prepost*stimtype*group, data=d_all_1921)
qqPlot(fitQQ20log , main="QQ Plot")
plot(fitQQ20log)
# Extract the residuals
fit_residuals20log <- residuals(object = fitQQ20log)
# Shapiro-Wilk test - violated
shapiro.test(x = fit_residuals20log)
# Kolmogorov-Smirnov test - sight violation (.043)
ks.test(fit_residuals20log, "pnorm", mean=mean(fit_residuals20log), sd=sd(fit_residuals20log))


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
leveneTest(lognoiseamp ~ prepost*stimtype*group, data = d_all_35)
plot(fitQQ4log, 1)
# not violated

## 20 Hz
leveneTest(lognoiseamp ~ prepost*stimtype*group, data = d_all_1921)
plot(fitQQ20log, 1)
# violated



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
fit4             = lmer(lognoiseamp  ~ group*stimtype*prepost + (1|subject), data=d_all_35)
summary(fit4)

Anova(fit4, type = "III", test.statistic = "F")


# 20 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit20             = lmer(lognoiseamp  ~ group*stimtype*prepost + (1|subject), data=d_all_1921)
summary(fit20)

Anova(fit20, type = "III", test.statistic = "F")

# without outlying subjects
options(contrasts=c("contr.sum", "contr.poly"))
fit20             = lmer(lognoiseamp  ~ group*stimtype*prepost + (1|subject), data=d_all_20_test)
summary(fit20)

Anova(fit20, type = "III", test.statistic = "F")

ph20              = emmeans(fit20, specs = pairwise  ~ group, adjust = "bonferroni")              
ph20$contrasts

# post-hoc
GGEE20            = subset(d_all_20_final, group%in%c('GG_EE'))
median(GGEE20$noiseamp)
IQR(GGEE20$noiseamp)
mean(GGEE20$noiseamp)
sd(GGEE20$noiseamp)

AC20            = subset(d_all_20_final, group%in%c('ActiveControl'))
median(AC20$noiseamp)
IQR(AC20$noiseamp)
mean(AC20$noiseamp)
sd(AC20$noiseamp)


#
# nlme package #
#

lme20  <-lme(lognoiseamp ~ group*stimtype*prepost, random =~1|subject, data=d_all_20_final, method="ML", weights=varIdent(form=~1|group))
summary(lme20)


# 
# non-parametric approach using aovp #
# 

## 4 Hz
perm4 = aovp(lognoiseamp ~ group*stimtype*prepost + Error(subject), data=d_all_4_final, perm="Exact")
summary(perm4)

## 20 Hz
perm20 = aovp(lognoiseamp ~ group*stimtype*prepost + Error(subject), data=d_all_20_final, perm="Exact")
summary(perm20)

