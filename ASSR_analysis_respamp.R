# --------------------------------------------------------- #
# name              ASSR_analysis_respamp.R                 #
# created by        Shauni Van Herck                        #
# description       used to analyse ASSR data (resp amp)    #
# version.string    R version 3.5.1 (2018-07-02)            #
# platform          x86_64-w64-mingw32                      #
# date created      07/10/2019                              #
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
ASSRdata$stimtype    = (str_extract(sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[2])), "[aA-zZ]+"))
#strsplit to keep AM4, ... + str_extract to only keep letter part

# create column frequency
ASSRdata$frequency    = ASSRdata$StimulationFrequency

# create column condition (stimtype + freq)
ASSRdata$condition    = sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[2]))

# create column testing phase (pre/post)
ASSRdata$prepost      = sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[3]))

# only keep relevant columns
ASSRdata              = ASSRdata[c("subject", "group", "condition", "stimtype", "frequency", "prepost", "Channel", "AssrHt2BiasedResponseAmplitude", "AssrHt2RecordingNoiseAmplitude")]

# subset data
ASSRdata_all          = subset(ASSRdata, Channel%in%c('All'))
colnames(ASSRdata_all)= c("subject", "group", "cond", "stimtype", "freq", "prepost", "channel", "respamp", "noiseamp")
head(ASSRdata_all)

ASSRdata_all$subject  = factor(ASSRdata_all$subject)
ASSRdata_all$group    = factor(ASSRdata_all$group, levels=c("GG_EE", "GG_NE", "ActiveControl", "PassiveControl"))    #reorder the levels
ASSRdata_all$cond     = factor(ASSRdata_all$cond, levels=c("AM4", "AM20", "PULS4", "PULS20"))
ASSRdata_all$stimtype = factor(ASSRdata_all$stimtype, levels=c("AM", "PULS"))
ASSRdata_all$freq     = factor(ASSRdata_all$freq, levels=c("4", "20"))
ASSRdata_all$prepost  = factor(ASSRdata_all$prepost, levels=c("pre", "post")) #reorder the levels
ASSRdata_all$channel  = factor(ASSRdata_all$channel)


ASSRdata_all_4             = subset(ASSRdata_all, freq%in%c('4'))  
ASSRdata_all_4             = subset(ASSRdata_all_4, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_4$row         = 1:nrow(ASSRdata_all_4)                                 #know row number of Cook's data points (see 'check assumptions')
ASSRdata_all_4_final       = subset(ASSRdata_all_4, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_20            = subset(ASSRdata_all, freq%in%c('20'))                  
ASSRdata_all_20            = subset(ASSRdata_all_20, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_20$row        = 1:nrow(ASSRdata_all_20)                                
ASSRdata_all_20_cook       = ASSRdata_all_20[ -c(292, 55, 56), ] 
ASSRdata_all_20$logrespamp = log(ASSRdata_all_20$respamp)
ASSRdata_all_20_final      = subset(ASSRdata_all_20, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))



ASSRdata_all_AM4                     = subset(ASSRdata_all, cond%in%c('AM4'))
ASSRdata_all_AM4_final               = subset(ASSRdata_all_AM4, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_AM20                    = subset(ASSRdata_all, cond%in%c('AM20'))
ASSRdata_all_AM20_final              = subset(ASSRdata_all_AM20, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_AM20_final$logrespamp   = log(ASSRdata_all_AM20_final$respamp)
ASSRdata_all_PULS4                   = subset(ASSRdata_all, cond%in%c('PULS4'))
ASSRdata_all_PULS4_final             = subset(ASSRdata_all_PULS4, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_PULS20                  = subset(ASSRdata_all, cond%in%c('PULS20'))
ASSRdata_all_PULS20_final            = subset(ASSRdata_all_PULS20, group%in%c('GG_EE', 'GG_NE', 'ActiveControl'))
ASSRdata_all_PULS20_final$logrespamp = log(ASSRdata_all_PULS20_final$respamp)

ASSRdata_all_AM4_melt          = melt(ASSRdata_all_AM4, id.var=c("subject", "group", "prepost"), measure.var=c("respamp"))
ASSRdata_all_AM4_wide          = cast(subject + group ~ prepost , data=ASSRdata_all_AM4_melt, mean)
ASSRdata_all_AM4_wide          = ASSRdata_all_AM4_wide[c("subject", "group", "pre", "post")]
head(ASSRdata_all_AM4_wide)

ASSRdata_all_AM20_melt          = melt(ASSRdata_all_AM20, id.var=c("subject", "group", "prepost"), measure.var=c("respamp"))
ASSRdata_all_AM20_wide          = cast(subject + group ~ prepost , data=ASSRdata_all_AM20_melt, mean)
ASSRdata_all_AM20_wide          = ASSRdata_all_AM20_wide[c("subject", "group", "pre", "post")]
head(ASSRdata_all_AM20_wide)

ASSRdata_all_PULS4_melt          = melt(ASSRdata_all_PULS4, id.var=c("subject", "group", "prepost"), measure.var=c("respamp"))
ASSRdata_all_PULS4_wide          = cast(subject + group ~ prepost , data=ASSRdata_all_PULS4_melt, mean)
ASSRdata_all_PULS4_wide          = ASSRdata_all_PULS4_wide[c("subject", "group", "pre", "post")]
head(ASSRdata_all_PULS4_wide)

ASSRdata_all_PULS20_melt          = melt(ASSRdata_all_PULS20, id.var=c("subject", "group", "prepost"), measure.var=c("respamp"))
ASSRdata_all_PULS20_wide          = cast(subject + group ~ prepost , data=ASSRdata_all_PULS20_melt, mean)
ASSRdata_all_PULS20_wide          = ASSRdata_all_PULS20_wide[c("subject", "group", "pre", "post")]
head(ASSRdata_all_PULS20_wide)

# plot data ------------------

setwd("C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/R/Plots")

#
## Response amplitude boxplots ##
#

# AM4
boxplot1         = ggplot(ASSRdata_all_AM4, aes(x=group, y=respamp, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                                        # position dodge + geom_jitter(position_dodge) is for
  #  geom_jitter(position=position_dodge(0.8)) +                                                      # getting dots with multiple groups (so that dots for)
  labs(title="AM4",x="group", y = "response amplitude (µV)", fill="test phase") + theme_classic() +   # pre and post are also separated per group
  coord_cartesian(ylim=c(0,3)) 
print(boxplot1)

# AM20
boxplot2         = ggplot(ASSRdata_all_AM20, aes(x=group, y=respamp, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
  #  geom_jitter(position=position_dodge(0.8)) +                                         
  labs(title="AM20",x="group", y = "response amplitude (µV)", fill="test phase") + theme_classic() + 
  coord_cartesian(ylim=c(0,1))
print(boxplot2)

# PULS4
boxplot3         = ggplot(ASSRdata_all_PULS4, aes(x=group, y=respamp, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
  #  geom_jitter(position=position_dodge(0.8)) +                                         
  labs(title="PULS4",x="group", y = "response amplitude (µV)", fill="test phase") + theme_classic() +
  coord_cartesian(ylim=c(0,5))
#print(boxplot3)

# PULS20
boxplot4         = ggplot(ASSRdata_all_PULS20, aes(x=group, y=respamp, fill=prepost)) + 
  geom_boxplot(position=position_dodge(0.8)) +                                        
  #  geom_jitter(position=position_dodge(0.8)) +                                         
  labs(title="PULS20",x="group", y = "response amplitude (µV)p", fill="test phase") + theme_classic() +
  coord_cartesian(ylim=c(0,1))
print(boxplot4)


tiff("respamp_allconditions_box_preliminary.tiff",width=1360,height=1360)                 #save figure in wd
print(grid.arrange(boxplot1, boxplot2, boxplot3, boxplot4, ncol = 2, top="Boxplots"))
dev.off()

#
## response and noise amplitude on one plot ##
## (greyscale)                              ##
#

# AM4
boxplot5         = ggplot(ASSRdata_all_AM4, aes(x=group)) + 
  geom_boxplot(aes(y=respamp, fill=prepost, colour="Response amplitude"), position=position_dodge(0.8)) +                                                       
  labs(title="AM4",x="group", y = "amplitude (µV)", fill="test phase", colour="measure") + 
  theme_classic() + scale_fill_grey() + scale_color_grey() + coord_cartesian(ylim=c(0,5))
print(boxplot5)

boxplot5         = boxplot5 + geom_boxplot(aes(y=noiseamp, fill=prepost, colour="Noise amplitude"))
print(boxplot5)

# AM20
boxplot6         = ggplot(ASSRdata_all_AM20, aes(x=group)) + 
  geom_boxplot(aes(y=respamp, fill=prepost, colour="Response amplitude"), position=position_dodge(0.8)) +                                                       
  labs(title="AM20",x="group", y = "amplitude (µV)", fill="test phase", colour="measure") + 
  theme_classic() + scale_fill_grey() + scale_color_grey() + coord_cartesian(ylim=c(0,1))
print(boxplot6)

boxplot6         = boxplot6 + geom_boxplot(aes(y=noiseamp, fill=prepost, colour="Noise amplitude"))
print(boxplot6)

# PULS4
boxplot7         = ggplot(ASSRdata_all_PULS4, aes(x=group)) + 
  geom_boxplot(aes(y=respamp, fill=prepost, colour="Response amplitude"), position=position_dodge(0.8)) +                                                       
  labs(title="PULS4",x="group", y = "amplitude (µV)", fill="test phase", colour="measure") + 
  theme_classic() + scale_fill_grey() + scale_color_grey() + coord_cartesian(ylim=c(0,5))
print(boxplot7)

boxplot7         = boxplot7 + geom_boxplot(aes(y=noiseamp, fill=prepost, colour="Noise amplitude"))
print(boxplot7)

# PULS20
boxplot8         = ggplot(ASSRdata_all_PULS20, aes(x=group)) + 
  geom_boxplot(aes(y=respamp, fill=prepost, colour="Response amplitude"), position=position_dodge(0.8)) +                                                       
  labs(title="PULS20",x="group", y = "amplitude (µV)", fill="test phase", colour="measure") + 
  theme_classic() + scale_fill_grey() + scale_color_grey() + coord_cartesian(ylim=c(0,1))
print(boxplot8)

boxplot8         = boxplot8 + geom_boxplot(aes(y=noiseamp, fill=prepost, colour="Noise amplitude"))
print(boxplot8)


tiff("4Hz_grey.tiff",width=1360,height=680)                 #save figure in wd
print(grid.arrange(boxplot5, boxplot7, ncol = 2, top="4 Hz"))
dev.off()

tiff("20Hz_grey.tiff",width=1360,height=680)                 #save figure in wd
print(grid.arrange(boxplot6, boxplot8, ncol = 2, top="20 Hz"))
dev.off()


#
## response and noise amplitude on one plot ##
## (colour scale)                           ##
#

# AM4
boxplot9         = ggplot(ASSRdata_all_AM4, aes(x=group)) + 
  geom_boxplot(aes(y=respamp, fill=prepost, colour="Response amplitude"), position=position_dodge(0.8)) +                                                       
  labs(title="AM4",x="group", y = "amplitude (µV)", fill="test phase", colour="measure") + 
  theme_classic() + scale_color_grey() + coord_cartesian(ylim=c(0,5))
print(boxplot9)

boxplot9         = boxplot9 + geom_boxplot(aes(y=noiseamp, fill=prepost, colour="Noise amplitude"))
print(boxplot9)

# AM20
boxplot10         = ggplot(ASSRdata_all_AM20, aes(x=group)) + 
  geom_boxplot(aes(y=respamp, fill=prepost, colour="Response amplitude"), position=position_dodge(0.8)) +                                                       
  labs(title="AM20",x="group", y = "amplitude (µV)", fill="test phase", colour="measure") + 
  theme_classic() + scale_color_grey() + coord_cartesian(ylim=c(0,1))
print(boxplot10)

boxplot10         = boxplot10 + geom_boxplot(aes(y=noiseamp, fill=prepost, colour="Noise amplitude"))
print(boxplot10)

# PULS4
boxplot11         = ggplot(ASSRdata_all_PULS4, aes(x=group)) + 
  geom_boxplot(aes(y=respamp, fill=prepost, colour="Response amplitude"), position=position_dodge(0.8)) +                                                       
  labs(title="PULS4",x="group", y = "amplitude (µV)", fill="test phase", colour="measure") + 
  theme_classic() + scale_color_grey() + coord_cartesian(ylim=c(0,5))
print(boxplot11)

boxplot11         = boxplot11 + geom_boxplot(aes(y=noiseamp, fill=prepost, colour="Noise amplitude"))
print(boxplot11)

# PULS20
boxplot12         = ggplot(ASSRdata_all_PULS20, aes(x=group)) + 
  geom_boxplot(aes(y=respamp, fill=prepost, colour="Response amplitude"), position=position_dodge(0.8)) +                                                     
  labs(title="PULS20",x="group", y = "amplitude (µV)", fill="test phase", colour="measure") + 
  theme_classic() + scale_color_grey() + coord_cartesian(ylim=c(0,1))
print(boxplot12)

boxplot12         = boxplot12 + geom_boxplot(aes(y=noiseamp, fill=prepost, colour="Noise amplitude"))
print(boxplot12)


tiff("4Hz_colour.tiff",width=1360,height=680)                 #save figure in wd
print(grid.arrange(boxplot9, boxplot11, ncol = 2, top="4 Hz"))
dev.off()

tiff("20Hz_colour.tiff",width=1360,height=680)                 #save figure in wd
print(grid.arrange(boxplot10, boxplot12, ncol = 2, top="20 Hz"))
dev.off()

#
## spaghetti plots ##
#

# spaghettiplot AM4
ASSRdata_all_AM4_1 = subset(ASSRdata_all_AM4, group=="GG_EE")
ASSRdata_all_AM4_2 = subset(ASSRdata_all_AM4, group=="GG_NE")
ASSRdata_all_AM4_3 = subset(ASSRdata_all_AM4, group=="ActiveControl")

p1            = ggplot(data = ASSRdata_all_AM4_1, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_EE", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")
p2            = ggplot(data = ASSRdata_all_AM4_2, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_NE", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")
p3            = ggplot(data = ASSRdata_all_AM4_3, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="Active control", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")

tiff("respamp_AM4_spaghetti_preliminary.tiff",width=1360,height=1360)
print(grid.arrange(p1, p2, p3, ncol = 2, top="Spaghetti plot for AM4"))
dev.off()


# spaghettiplot AM20
ASSRdata_all_AM20_1 = subset(ASSRdata_all_AM20, group=="GG_EE")
ASSRdata_all_AM20_2 = subset(ASSRdata_all_AM20, group=="GG_NE")
ASSRdata_all_AM20_3 = subset(ASSRdata_all_AM20, group=="ActiveControl")

p1            = ggplot(data = ASSRdata_all_AM20_1, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_EE", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")
p2            = ggplot(data = ASSRdata_all_AM20_2, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_NE", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")
p3            = ggplot(data = ASSRdata_all_AM20_3, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="Active control", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")

tiff("respamp_AM20_spaghetti_preliminary.tiff",width=1360,height=1360)
print(grid.arrange(p1, p2, p3, ncol = 2, top="Spaghetti plot for AM20"))
dev.off()
  
# spaghettiplot PULS4
ASSRdata_all_PULS4_1 = subset(ASSRdata_all_PULS4, group=="GG_EE")
ASSRdata_all_PULS4_2 = subset(ASSRdata_all_PULS4, group=="GG_NE")
ASSRdata_all_PULS4_3 = subset(ASSRdata_all_PULS4, group=="ActiveControl")
  
p1            = ggplot(data = ASSRdata_all_PULS4_1, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_EE", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")
p2            = ggplot(data = ASSRdata_all_PULS4_2, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_NE", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")
p3            = ggplot(data = ASSRdata_all_PULS4_3, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="Active control", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")
 
tiff("respamp_PULS4_spaghetti_preliminary.tiff",width=1360,height=1360) 
print(grid.arrange(p1, p2, p3, ncol = 2, top="Spaghetti plot for PULS4"))
dev.off()

# spaghettiplot PULS20
ASSRdata_all_PULS20_1 = subset(ASSRdata_all_PULS20, group=="GG_EE")
ASSRdata_all_PULS20_2 = subset(ASSRdata_all_PULS20, group=="GG_NE")
ASSRdata_all_PULS20_3 = subset(ASSRdata_all_PULS20, group=="ActiveControl")

p1            = ggplot(data = ASSRdata_all_PULS20_1, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_EE", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")
p2            = ggplot(data = ASSRdata_all_PULS20_2, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="GG_NE", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")
p3            = ggplot(data = ASSRdata_all_PULS20_3, aes(x = prepost, y = respamp, group = subject, colour=subject)) +
  geom_line(size=1) + labs(title="Active control", x="test phase", y="response amplitude (µV)") + theme_bw() +
  theme(legend.position="none")

tiff("respamp_PULS20_spaghetti_preliminary.tiff",width=1360,height=1360)
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
fitQQ4 <- lm(respamp ~ prepost*stimtype*group, data=ASSRdata_all_4)
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
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), sd=sd(fit_residuals4))


## 20 Hz
fitQQ20 <- lm(respamp ~ prepost*stimtype*group, data=ASSRdata_all_20)
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
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals20, "pnorm", mean=mean(fit_residuals20), sd=sd(fit_residuals20))


# 20 Hz with data points based on Cook's distance removed
fitQQ20cook <- lm(respamp ~ prepost*stimtype*group, data=ASSRdata_all_20_cook)
qqPlot(fitQQ20cook, main="QQ Plot")
plot(fitQQ20cook)
fit_residuals20cook <- residuals(object = fitQQ20cook)
# Shapiro-Wilk test -violated
shapiro.test(fit_residuals20cook)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals20cook, "pnorm", mean=mean(fit_residuals20cook), sd=sd(fit_residuals20cook))

# 20 Hz with log transformed data
fitQQ20log <- lm(logrespamp ~ prepost*stimtype*group, data=ASSRdata_all_20)
qqPlot(fitQQ20log , main="QQ Plot")
plot(fitQQ20log)
library(MASS)
sresid <- studres(fitQQ20log) 
hist(sresid, freq=FALSE, main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40) 
yfit<-dnorm(xfit) 
lines(xfit, yfit)
# Extract the residuals
fit_residuals20log <- residuals(object = fitQQ20log)
# Shapiro-Wilk test - violated
shapiro.test(x = fit_residuals20log)
# Kolmogorov-Smirnov test - not violated
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
leveneTest(respamp ~ prepost*stimtype*group, data = ASSRdata_all_4)
# violated

## 20 Hz
leveneTest(logrespamp ~ prepost*stimtype*group, data = ASSRdata_all_20)
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
fit4             = lmer(respamp  ~ group*stimtype*prepost + (1|subject), data=ASSRdata_all_4_final)
summary(fit4)

Anova(fit4, type = "III", test.statistic = "F")

# post-hoc
AM4              = subset(ASSRdata_all_4_final, stimtype%in%c('AM'))
median(AM4$respamp)
IQR(AM4$respamp)
mean(AM4$respamp)
sd(AM4$respamp)

PULS4            = subset(ASSRdata_all_4_final, stimtype%in%c('PULS'))
median(PULS4$respamp)
IQR(PULS4$respamp)
mean(PULS4$respamp)
sd(PULS4$respamp)


# 20 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit20             = lmer(logrespamp  ~ group*stimtype*prepost + (1|subject), data=ASSRdata_all_20_final)
summary(fit20)

Anova(fit20, type = "III", test.statistic = "F")

# post-hoc
AM20              = subset(ASSRdata_all_20_final, stimtype%in%c('AM'))
median(AM20$respamp)
IQR(AM20$respamp)
mean(AM20$respamp)
sd(AM20$respamp)

PULS20            = subset(ASSRdata_all_20_final, stimtype%in%c('PULS'))
median(PULS20$respamp)
IQR(PULS20$respamp)
mean(PULS20$respamp)
sd(PULS20$respamp)

PRE20             = subset(ASSRdata_all_20_final, prepost%in%c('pre'))
median(PRE20$respamp)
IQR(PRE20$respamp)
mean(PRE20$respamp)
sd(PRE20$respamp)

POST20            = subset(ASSRdata_all_20_final, prepost%in%c('post'))
median(POST20$respamp)
IQR(POST20$respamp)
mean(POST20$respamp)
sd(POST20$respamp)



# 
# non-parametric approach using aovp #
# 

## 4 Hz
perm4 = aovp(respamp ~ group*stimtype*prepost + Error(subject), data=ASSRdata_all_4_final, perm="Exact")
summary(perm4)
# different results than with lmer

## 20 Hz
perm20 = aovp(logrespamp ~ group*stimtype*prepost + Error(subject), data=ASSRdata_all_20_final, perm="Exact")
summary(perm20)

# check assumptions AM -----------------

#
# normality assumption #
#

## 4 Hz 
fitQQ4 <- lm(respamp ~ prepost*group, data=ASSRdata_all_AM4_final)
qqPlot(fitQQ4, main="QQ Plot")
plot(fitQQ4)
# Extract the residuals
fit_residuals4 <- residuals(object = fitQQ4)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), sd=sd(fit_residuals4))


## 20 Hz
fitQQ20 <- lm(respamp ~ prepost*group, data=ASSRdata_all_AM20_final)
qqPlot(fitQQ20, main="QQ Plot")
plot(fitQQ20)
# Extract the residuals
fit_residuals20 <- residuals(object = fitQQ20)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals20, "pnorm", mean=mean(fit_residuals20), sd=sd(fit_residuals20))

## 20 Hz (done on log)
fitQQ20 <- lm(logrespamp ~ prepost*group, data=ASSRdata_all_AM20_final)
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
leveneTest(respamp ~ prepost*group, data = ASSRdata_all_AM4_final)
plot(fitQQ4, 1)
# not violated

## 20 Hz
leveneTest(logrespamp ~ prepost*group, data = ASSRdata_all_AM20_final)
plot(fitQQ20,1)
# not violated

# analyses AM ------------------


## 4 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit4             = lmer(respamp  ~ group*prepost + (1|subject), data=ASSRdata_all_AM4_final)
summary(fit4)

Anova(fit4, type = "III", test.statistic = "F")

## 20 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit20             = lmer(logrespamp  ~ group*prepost + (1|subject), data=ASSRdata_all_AM20_final)
summary(fit20)

Anova(fit20, type = "III", test.statistic = "F")

# post-hoc

PRE20             = subset(ASSRdata_all_AM20_final, prepost%in%c('pre'))
mean(PRE20$respamp)
sd(PRE20$respamp)

POST20            = subset(ASSRdata_all_AM20_final, prepost%in%c('post'))
mean(POST20$respamp)
sd(POST20$respamp)




# check assumptions PULS ---------------

#
# normality assumption #
#

## 4 Hz 
fitQQ4 <- lm(respamp ~ prepost*group, data=ASSRdata_all_PULS4_final)
qqPlot(fitQQ4, main="QQ Plot")
plot(fitQQ4)
# Extract the residuals
fit_residuals4 <- residuals(object = fitQQ4)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), sd=sd(fit_residuals4))


## 20 Hz
fitQQ20 <- lm(respamp ~ prepost*group, data=ASSRdata_all_PULS20_final)
qqPlot(fitQQ20, main="QQ Plot")
plot(fitQQ20)
# Extract the residuals
fit_residuals20 <- residuals(object = fitQQ20)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals20, "pnorm", mean=mean(fit_residuals20), sd=sd(fit_residuals20))

## 20 Hz (done on log)
fitQQ20 <- lm(logrespamp ~ prepost*group, data=ASSRdata_all_PULS20_final)
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
leveneTest(respamp ~ prepost*group, data = ASSRdata_all_PULS4_final)
plot(fitQQ4, 1)
# not violated

## 20 Hz
leveneTest(logrespamp ~ prepost*group, data = ASSRdata_all_PULS20_final)
plot(fitQQ20,1)
# not violated



# analyses PULS --------------

## 4 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit4             = lmer(respamp  ~ group*prepost + (1|subject), data=ASSRdata_all_PULS4_final)
summary(fit4)

Anova(fit4, type = "III", test.statistic = "F")

## 20 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit20             = lmer(logrespamp  ~ group*prepost + (1|subject), data=ASSRdata_all_PULS20_final)
summary(fit20)

Anova(fit20, type = "III", test.statistic = "F")

# post-hoc

PRE20             = subset(ASSRdata_all_PULS20_final, prepost%in%c('pre'))
mean(PRE20$respamp)
sd(PRE20$respamp)

POST20            = subset(ASSRdata_all_PULS20_final, prepost%in%c('post'))
mean(POST20$respamp)
sd(POST20$respamp)

ph20              = emmeans(fit20, specs = pairwise  ~ group:prepost, adjust = "holm")              
ph20$contrasts





