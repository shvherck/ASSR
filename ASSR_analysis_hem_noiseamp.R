# ----------------------------------------------------------------------------------------- #
# name              ASSR_analysis_hem_noiseamp.R                                            #
# created by        Shauni Van Herck                                                        #
# description       used to analyse ASSR data for hemispheres (noiseamp)                    #
# version.string    R version 3.5.1 (2018-07-02)                                            #
# platform          x86_64-w64-mingw32                                                      #
# date created      11/05/2020                                                              #
# ----------------------------------------------------------------------------------------- #



# load libraries -------------------
# reading
library(readxl)
# reshaping
library(stringr)      # split data into numbers/letters
library(reshape)
library(reshape2)
library(dplyr)
# plotting
library(car) #qqplot
library(ggplot2)
library(grid)
library(gridExtra)
# load library for robust lmm
library(robustlmm)
library(pbkrtest)
library(WRS2)
# mixed models
library(lme4)
library(nlme)
# permutations
library(lmPerm)
# post-hoc comparisons
library(emmeans)
library(phia)

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
d              = d[c("subject", "group", "condition", "stimtype", "frequency", "prepost", "hem", "AssrHt2RecordingNoiseAmplitude")]
colnames(d)    = c("subject", "group", "cond", "stimtype", "freq", "prepost", "hem", "noiseamp")
head(d)

d$subject   = factor(d$subject)
d$group     = factor(d$group, levels=c("GG_EE", "GG_NE", "ActiveControl"))    #reorder the levels
d$cond      = factor(d$cond, levels=c("AM4", "AM20", "PULS4", "PULS20"))
d$stimtype  = factor(d$stimtype, levels=c("AM", "PULS"))
d$freq      = factor(d$freq, levels=c("4", "20"))
d$prepost   = factor(d$prepost, levels=c("pre", "post")) #reorder the levels
d$hem       = factor(d$hem)


d_4              = subset(d, freq%in%c('4'))

d_20             = subset(d, freq%in%c('20')) 


d_AM4       = subset(d, cond%in%c('AM4'))
d_AM20      = subset(d, cond%in%c('AM20'))
d_PULS4     = subset(d, cond%in%c('PULS4'))
d_PULS20    = subset(d, cond%in%c('PULS20'))


# plot data: plots can be found together with respamp plots in ASSR_analysis_hem_respamp.R -------------------------


# check assumptions general -------------------

#
# normality assumption #
#

## 4 Hz 
fitQQ4 <- lm(noiseamp ~ prepost*stimtype*hem*group, data=d_4)
qqPlot(fitQQ4, main="QQ Plot")
# Extract the residuals
fit_residuals4 <- residuals(object = fitQQ4)
# Kolmogorov-Smirnov test - violated
ks.test(fit_residuals4, "pnorm", mean=mean(fit_residuals4), sd=sd(fit_residuals4))
shapiro.test(x = fit_residuals4 )


## 20 Hz
fitQQ20 <- lm(noiseamp ~ prepost*stimtype*hem*group, data=d_20)
qqPlot(fitQQ20, main="QQ Plot")
# Extract the residuals
fit_residuals20 <- residuals(object = fitQQ20)
# Kolmogorov-Smirnov test - not violated
ks.test(fit_residuals20, "pnorm", mean=mean(fit_residuals20), sd=sd(fit_residuals20))
shapiro.test(x = fit_residuals20 )



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
leveneTest(noiseamp ~ prepost*stimtype*hem*group, data = d_4)
plot(fitQQ4, 1)
# violated

## 20 Hz
leveneTest(noiseamp ~ prepost*stimtype*hem*group, data = d_20)
plot(fitQQ20, 1)
# not violated


# analyses general -------------------


#
# robust lmm #
#

## 
anova_merMod<-function(model,rand,w=NULL,seed=round(runif(1,0,100),0),nsim=10000){
  data<-model@frame
  if(!is.null(w)){
    data<-data[,-grep("(weights)",names(data))]
  }
  
  resp<-names(model.frame(model))[1]
  #generate a list of reduced model formula
  fs<-list()
  fs[[1]]<-as.formula(paste(resp,"~ 1 +",rand))
  nb_terms<-length(attr(terms(model),"term.labels"))
  if(nb_terms>1){
    for(i in 1:nb_terms){
      tmp<-c(attr(terms(model),"term.labels")[1:i],rand)
      fs[[i+1]]<-reformulate(tmp,response=resp)
    }      
  }
  
  #fit the reduced model to the data
  
  fam<-family(model)[1]$family
  if(fam=="gaussian"){
    m_fit<-lapply(fs,function(x) lmer(x,data,REML=FALSE))
  } else if(fam=="binomial"){
    m_fit<-lapply(fs,function(x) glmer(x,data,family=fam,weights=w))
  }  else{
    m_fit<-lapply(fs,function(x) glmer(x,data,family=fam))
  }
  
  #compare nested model with one another and get LRT values (ie increase in the likelihood of the models as parameters are added)
  tab_out<-NULL
  
  for(i in 1:(length(m_fit)-1)){
    comp<-PBmodcomp(m_fit[[i+1]],m_fit[[i]],seed=seed,nsim=nsim)    
    term_added<-attr(terms(m_fit[[i+1]]),"term.labels")[length(attr(terms(m_fit[[i+1]]),"term.labels"))]
    #here are reported the bootstrapped p-values, ie not assuming any parametric distribution like chi-square to the LRT values generated under the null model
    #these p-values represent the number of time the simulated LRT value (under null model) are larger than the observe one
    tmp<-data.frame(term=term_added,LRT=comp$test$stat[1],p_value=comp$test$p.value[2])
    tab_out<-rbind(tab_out,tmp)
    print(paste("Variable ",term_added," tested",sep=""))
  }  
  print(paste("Seed set to:",seed))
  return(tab_out)  
}
##

## 4 Hz
rfit4 <- rlmer(noiseamp  ~ group*stimtype*prepost*hem + (1|subject), data = d_4) # robust mixed model with rlmer (robustlmm package)

summary(rfit4)

getME(rfit4, "w_e") # all robustness weights for the residuals
getME(rfit4, "w_b") # all robustness weights for the random effects

plot(rfit4) # QQ plot etc. that also indicate weights

set.seed(123)
anova_merMod(model=rfit4,rand="(1|subject)") # uses below function that wraps around the PBmodcomp function to compute bootstrapped p-values 
# for each term in a model by sequentially adding them
q
# post-hoc on main effect stimtype
d_4_PULS  = subset(d_4, stimtype%in%c("PULS"))
mean(d_4_PULS$noiseamp, 0.2)
d_4_AM    = subset(d_4, stimtype%in%c("AM"))
mean(d_4_AM$noiseamp, 0.2)

# post-hoc on interaction effect stimtype x hem
d_ph <- dcast(d_4, subject ~ stimtype + hem, value.var="noiseamp", fun.aggregate = mean)
head(d_ph)
# this format is needed for post-hoc test
# post hoc test by performing wilcox test / yuen dependent samples test (with trimmed means) on all conditions of interest

AMleft_vs_AMright     = wilcox.test(d_ph$AM_Left,                            # .002 (holm)
                                    d_ph$AM_Right,                           
                                    paired=TRUE)$p.value                       
PULSleft_vs_PULSright = wilcox.test(d_ph$PULS_Left,                          # .00014 (holm)                       
                                    d_ph$PULS_Right,                          
                                    paired=TRUE)$p.value              

p.adjust(c(AMleft_vs_AMright, PULSleft_vs_PULSright), method="holm")

# post-hoc with yuen (uses trimmed means)
AMleft_vs_AMright_yuen     = yuend(x = d_ph$AM_Left, y = d_ph$AM_Right)$p.value
PULSleft_vs_PULSright_yuen = yuend(x = d_ph$PULS_Left, y = d_ph$PULS_Right)$p.value

p.adjust(c(AMleft_vs_AMright_yuen, PULSleft_vs_PULSright_yuen), method="holm")

mean(d_ph$AM_Left, 0.2)
mean(d_ph$AM_Right, 0.2)
mean(d_ph$PULS_Left, 0.2, na.rm=TRUE)
mean(d_ph$PULS_Right, 0.2, na.rm=TRUE)


AMleft_vs_PULSleft    = wilcox.test(d_ph$AM_Left,
                                    d_ph$PULS_Left,
                                    paired=TRUE)$p.value
AMright_vs_PULSright  = wilcox.test(d_ph$AM_Right,
                                    d_ph$PULS_Right,
                                    paired=TRUE)$p.value

p.adjust(c(AMleft_vs_PULSleft, AMright_vs_PULSright), method="holm") # both significant but explained by the stimtype main effect

# post-hoc with yuen (uses trimmed means)
AMleft_vs_PULSleft_yuen    = yuend(x = d_ph$AM_Left, y = d_ph$PULS_Left)$p.value
AMright_vs_PULSright_yuen  = yuend(x = d_ph$AM_Right, y = d_ph$PULS_Right)$p.value

p.adjust(c(AMleft_vs_PULSleft_yuen, AMright_vs_PULSright_yuen), method="holm")




## 20 Hz
rfit20 <- rlmer(noiseamp  ~ group*stimtype*prepost*hem + (1|subject), data = d_20) # robust mixed model with rlmer (robustlmm package)

summary(rfit20)

getME(rfit20, "w_e") # all robustness weights for the residuals
getME(rfit20, "w_b") # all robustness weights for the random effects

plot(rfit20) # QQ plot etc. that also indicate weights

set.seed(123)
anova_merMod(model=rfit20,rand="(1|subject)") # uses below function that wraps around the PBmodcomp function to compute bootstrapped p-values 
# for each term in a model by sequentially adding them

# post-hoc on main effect stimtype
d_20_AM    = subset(d_20, stimtype%in%c("AM"))
mean(d_20_AM$noiseamp, 0.2)
d_20_PULS  = subset(d_20, stimtype%in%c("PULS"))
mean(d_20_PULS$noiseamp, 0.2)


# post-hoc on interaction effect group x hem
d_ph2 <- dcast(d_20, subject ~ group + hem, value.var="noiseamp", fun.aggregate = mean)
head(d_ph2)
# this format is needed for post-hoc test
# post hoc test by performing wilcox test on all conditions of interest

GGEEleft_vs_GGEEright     = wilcox.test(d_ph2$GG_EE_Left,                           
                                        d_ph2$GG_EE_Right,                           
                                        paired=TRUE,)$p.value                                                
GGNEleft_vs_GGNEright     = wilcox.test(d_ph2$GG_NE_Left,                                                   
                                        d_ph2$GG_NE_Right,                          
                                        paired=TRUE)$p.value 
ACleft_vs_ACright         = wilcox.test(d_ph2$ActiveControl_Left,
                                        d_ph2$ActiveControl_Right,
                                        paired=TRUE)$p.value

p.adjust(c(GGEEleft_vs_GGEEright, GGNEleft_vs_GGNEright, ACleft_vs_ACright), method="holm") # nothing signifcant 

GGEEleft_vs_GGEEright_yuen     = yuend(x = d_ph2$GG_EE_Left, d_ph2$GG_EE_Right)$p.value
GGNEleft_vs_GGNEright_yuen     = yuend(x = d_ph2$GG_NE_Left, y = d_ph2$GG_NE_Right)$p.value
ACleft_vs_ACright_yuen         = yuend(x = d_ph2$ActiveControl_Left, d_ph2$ActiveControl_Right)$p.value

p.adjust(c(GGEEleft_vs_GGEEright_yuen, GGNEleft_vs_GGNEright_yuen, ACleft_vs_ACright_yuen), method="holm") # nothing significant

d_20_ordered = d_20[order(d_20$group),]                     # data need to be ordered this way for Kruskal

d_20_left           = subset(d_20_ordered, (hem%in%c("Left")))
groupdiff_left      = kruskal.test(noiseamp ~ group, data=d_20_left)$p.value          

d_20_right          = subset(d_20_ordered, (hem%in%c("Right")))
groupdiff_right     = kruskal.test(noiseamp ~ group, data=d_20_right)$p.value         

p.adjust(c(groupdiff_left, groupdiff_right), method="holm") # nothing significant 

groupdiff_left_yuen      = t1way(noiseamp ~ group, data = d_20_left)$p.value
groupdiff_right_yuen     = t1way(noiseamp ~ group, data = d_20_left)$p.value

p.adjust(c(groupdiff_left_yuen, groupdiff_right_yuen), method="holm") # nothing significant 


#
# 'modern' mixed model approach #
#

## 4 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit4             = lmer(noiseamp  ~ group*stimtype*prepost*hem + (1|subject), data=d_4_cook)

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
mean(d_4_AML$noiseamp)
mean(d_4_AMR$noiseamp)

ph43              = emmeans(fit4, specs = pairwise  ~ prepost:hem, adjust = "holm")              
ph43$contrasts    # no significances even though there is a significant interaction

ph44              = emmeans(fit4, specs = pairwise  ~ group:stimtype:hem, adjust = "holm")              
ph44$contrasts 

## model building try out

library(MuMIn)

fit4             = lmer(noiseamp  ~ group*stimtype*prepost*hem + (1|subject), data=d_4_cook)
fit              = lmer(noiseamp  ~ group*stimtype*prepost*hem + (prepost|subject), data=d_4_cook)

r.squaredGLMM(fit4) # this model has the best marginal R squared (or at least doesn't differ from fit)
r.squaredGLMM(fit)

# m marginal fixed effects - c conditional marginal effects (usually interested in marginal)


# 20 Hz
options(contrasts=c("contr.sum", "contr.poly"))
fit20             = lmer(noiseamp  ~ group*stimtype*prepost*hem + (1|subject), data=d_20)

summary(fit20)

Anova(fit20, type = "III", test.statistic = "F")

# stimtype: AM > PULS
# prepost: POST > PRE
# hem: right > left
# group:hem
# GGEE left < GGEE right

anova_merMod(model=fit20,rand="(1|subject)") # uses below function that wraps around the PBmodcomp function to compute bootstrapped p-values 
# for each term in a model by sequentially adding them


# post-hoc
d_20_AM              = subset(d_20, stimtype%in%c('AM'))
mean(d_20_AM$noiseamp)
d_20_PULS            = subset(d_20, stimtype%in%c('PULS'))
mean(d_20_PULS$noiseamp)

d_20_pre             = subset(d_20, prepost%in%c('pre'))
mean(d_20_pre$noiseamp)
d_20_post            = subset(d_20, prepost%in%c('post'))
mean(d_20_post$noiseamp)

d_20_left            = subset(d_20, hem%in%c('Left'))
mean(d_20_left$noiseamp)
d_20_right            = subset(d_20, hem%in%c('Right'))
mean(d_20_right$noiseamp)


ph20            = emmeans(fit20, specs = pairwise  ~ group:hem, adjust = "holm")              
ph20$contrasts 

d_20_GGEEL      = subset(d_20, group%in%c("GG_EE") & hem%in%c("Left"))
mean(d_20_GGEEL$noiseamp)
d_20_GGEER      = subset(d_20, group%in%c("GG_EE") & hem%in%c("Right"))
mean(d_20_GGEER$noiseamp)




