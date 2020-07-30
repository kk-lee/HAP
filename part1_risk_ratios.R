
# load packages####
library(tableone)
library(caret)
library(reshape2)
library(tidyr)
library(dplyr)
library(binom)
library(survival)
library(gridExtra)
library(grid)
library(survMisc)
library(scales)
library(ggplot)
library(tidyr)
library(metafor)
library(xlsx)
library(xlsxjars)
library(rJava)
library(readxl)
library(stringr)
library(reshape)
library(data.table)
library(plyr)
library(RColorBrewer)
library(rworldmap)
library(classInt)
library(gdata)
library(forestplot)
library(Hmisc)
library(data.table)
library(tidyverse)
library(fmsb)
library(glue)
library(cowplot)

rm(list=ls())

###====================================
### Part 1 -  estimate pooled rate and risk ratio
###====================================

# load files

dat <- read.csv("data_raw/iap.data.v8.csv",header = T,na.strings = "")

names(dat) <- tolower(names(dat))
dat <- subset(dat,!is.na(dat$outcome))
dat$author <- paste0(dat$author," ","et al")


# sort out different analysis.codes
# ----------------------------------
# 1- polluting fuels vs clean fuels only
# 2- polluting fuels vs clean fuel + some polluting fuels (ie. Kerosene/ coal/ not biomass)
# 3- non-lpg vs lpg
# 4- lpg vs no lpg
# 5- polluting + clean fuels vs clean fuels
# 6- clean fuel vs polluting fuel
# 7- gas vs electricity
# 8- coal vs biomass
# 9- improved stove vs traditional stove
# 10- biomass vs kerosene/charcoal
# 11- smoky coal vs smokeless coal or wood
# 12- wood vs charcoal
# 13- wood vs coal
# 14- wood vs fuel oil
# 15- kerosene vs biomass
# 16- biomass vs coal
# 17- traditional stove vs improved stove
# 88- direct pollutant measurement
# 999- * not useful comparison

# reverse analysis codes 3, 6, 17
reverse <- function(x) {
  1/x
}

dat[dat$analysis.code%in% c(3,6,17), c('ratio.ce', 'ratio.ll', 'ratio.ul')] <- lapply(dat[dat$analysis.code%in% c(3,6,17), c('ratio.ce', 'ratio.ul', 'ratio.ll')],reverse)


dat <- dat %>% 
  mutate(analysis.code = case_when(
    analysis.code==3 ~ 4,
    analysis.code==6 ~ 1,
    analysis.code==17 ~ 9,
    TRUE ~ as.numeric(analysis.code)
  ))


#calculate risk ratio from 2X2 tables
dat <- dat %>% 
  mutate_at(vars(exposed.event, 
                 exposed.health, 
                 control.event,
                 control.health), 
            ~as.numeric(as.character(na_if(.,"NA"))))

dat2 <- dat %>% filter(risk.measure=="2x2")

dat2 <- escalc(measure="OR", ai=exposed.event, bi=exposed.health, ci=control.event, di=control.health, data=dat2, var.names=c("lgrr", "variance"))
dat2 <- dat2 %>% 
  mutate(se=sqrt(variance),
         ratio.ce=exp(lgrr),
         ratio.ll=exp(lgrr-1.96*se),
         ratio.ul=exp(lgrr+1.96*se)) %>% 
  select(-c(variance, lgrr, se)) 
  
dat <- rbind(dat[dat$risk.measure!="2x2",], as.data.frame(dat2))  


dat$lgrr <- logb(dat$ratio.ce)
dat$lgll <- logb(dat$ratio.ll)
dat$lgul <- logb(dat$ratio.ul)
dat$se <- (dat$lgul-dat$lgll)/(2*1.96)

num_iter <- 10000


dat <- dat %>% mutate(se=case_when(
  ratio.ll==0 ~ (lgul-lgrr)/1.96,
  TRUE~se
)) 


saveRDS(dat, "data_derived.v2/dat.rds")


##============================================
## calculate meta-estimates for risk ratio of main outcomes
##============================================

dat.ratio <- dat

dat.ratio <- subset(dat.ratio,!is.na(dat.ratio$ratio.ll))

dat.ratio <- dat.ratio %>% 
  select("study.id", "author", "pub.year", "overall.estimate", "analysis.code", "subgroup.id", "subgroup.gender", "subgroup.exposure", "subgroup.outcome", "subgroup.other", "interventional.study", "intervention.type", "exposed.group", "control.group", "outcome.recode", "outcome.type", "adult.paediatric.all","paed.under5", "risk.measure", "ratio.ce", "ratio.ll", "ratio.ul", "lgrr", "se", "exposed.event", "exposed.health", "control.event","control.health")

# diagnosis and mortality only
dat.ratio <- dat.ratio %>% 
  filter(overall.estimate==1 & analysis.code %in% c(1,2,5) & outcome.type!="symptoms")

outcomes<-as.data.frame(table(dat.ratio$outcome.recode))

#### Pool ratios by endpoint ####
dat.ratio.resp <- dat.ratio %>% 
  filter(outcome.recode %in% c("asthma", "copd", "ari", "lung cancer", "tuberculosis", "respiratory disease"))

dat.ratio.cvd <- dat.ratio %>% 
  filter(outcome.recode %in% c("stroke", "ihd", "cardiovascular disease"))

dat.ratio.lbw <- dat.ratio %>% 
  filter(outcome.recode == "lbw")

dat.ratio.under5 <- dat.ratio %>% 
  filter(outcome.recode == "under five death")

dat.ratio.stillbirth <- dat.ratio %>% 
  filter(outcome.recode == "stillbirth")

dat.ratio.alldeath <- dat.ratio %>% 
  filter(outcome.recode == "all cause death")

# mortality only
dat.ratio.mort <- dat.ratio %>% 
  filter(overall.estimate==1 & analysis.code %in% c(1,2,5) & outcome.type=="mortality")

dat.ratio.mort <- dat.ratio.mort %>% 
  mutate(outcome.recode = case_when(
    outcome.recode %in% c("asthma", "copd", "ari", "lung cancer", "tuberculosis", "respiratory disease") ~ "resp death",
    outcome.recode %in% c("stroke", "ihd", "cardiovascular disease") ~ "cv death",
    TRUE ~ as.character(outcome.recode)
  ))

outcomes.mort<-as.data.frame(table(dat.ratio.mort$outcome.recode))


#### respiratory disease ####

rr.res <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("asthma", "copd", "ari", "lung cancer", "tuberculosis", "respiratory disease"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("asthma", "copd", "ari", "lung cancer", "tuberculosis", "respiratory disease"),]$author) 
rr.resp <- predict(rr.res, transf = exp, digit=2)
rr.resp <- as.data.frame(rr.resp)
rr.resp$subgroup <- "respiratory overall"
rr.resp$i2 <- round(rr.res$I2, digits=1)
rr.resp$estimates <- rr.res$k

rr.asth <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("asthma"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("asthma"),]$author) 
rr.asthma <- predict(rr.asth, transf = exp, digit=2)
rr.asthma <- as.data.frame(rr.asthma)
rr.asthma$subgroup <- "asthma"
rr.asthma$i2 <- round(rr.asth$I2, digits=1)
rr.asthma$estimates <- rr.asth$k

rr.cop <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("copd"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("copd"),]$author) 
rr.copd <- predict(rr.cop, transf = exp, digits=2)
rr.copd <- as.data.frame(rr.copd)
rr.copd$subgroup <- "copd"
rr.copd$i2 <- round(rr.cop$I2, digits=1)
rr.copd$estimates <- rr.cop$k

rr.ar <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("ari"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("ari"),]$author) 
rr.ari <- predict(rr.ar, transf = exp, digit=2)
rr.ari <- as.data.frame(rr.ari)
rr.ari$subgroup <- "ari- overall"
rr.ari$i2 <- round(rr.ar$I2, digits=1)
rr.ari$estimates <- rr.ar$k

rr.adult.ar <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("ari") & dat.ratio$adult.paediatric.all %in% c("adult", "both"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("ari") & dat.ratio$adult.paediatric.all %in% c("adult", "both"),]$author) 
rr.adult.ari <- predict(rr.adult.ar, transf = exp, digit=2)
rr.adult.ari <- as.data.frame(rr.adult.ari)
rr.adult.ari$subgroup <- "ari- adult"
rr.adult.ari$i2 <- round(rr.adult.ar$I2, digits=1)
rr.adult.ari$estimates <- rr.adult.ar$k

rr.paed.ar <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("ari") & dat.ratio$adult.paediatric.all %in% c("paed"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("ari") & dat.ratio$adult.paediatric.all %in% c("paed"),]$author) 
rr.paed.ari <- predict(rr.paed.ar, transf = exp, digit=2)
rr.paed.ari <- as.data.frame(rr.paed.ari)
rr.paed.ari$subgroup <- "ari- paed"
rr.paed.ari$i2 <- round(rr.paed.ar$I2, digits=1)
rr.paed.ari$estimates <- rr.paed.ar$k

rr.paed.5t16.ar <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("ari") & dat.ratio$adult.paediatric.all %in% c("paed") & dat.ratio$paed.under5==0,],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("ari") & dat.ratio$adult.paediatric.all %in% c("paed") & dat.ratio$paed.under5==0,]$author) 
rr.paed.5t16.ari <- predict(rr.paed.5t16.ar, transf = exp, digit=2)
rr.paed.5t16.ari <- as.data.frame(rr.paed.5t16.ari)
rr.paed.5t16.ari$subgroup <- "ari- paed (5-16)"
rr.paed.5t16.ari$i2 <- round(rr.paed.5t16.ar$I2, digits=1)
rr.paed.5t16.ari$estimates <- rr.paed.5t16.ar$k

rr.paed.u5.ar <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("ari") & dat.ratio$adult.paediatric.all %in% c("paed") & dat.ratio$paed.under5==1,],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("ari") & dat.ratio$adult.paediatric.all %in% c("paed") & dat.ratio$paed.under5==1,]$author) 
rr.paed.u5.ari <- predict(rr.paed.u5.ar, transf = exp, digit=2)
rr.paed.u5.ari <- as.data.frame(rr.paed.u5.ari)
rr.paed.u5.ari$subgroup <- "ari- paed (under5)"
rr.paed.u5.ari$i2 <- round(rr.paed.u5.ar$I2, digits=1)
rr.paed.u5.ari$estimates <- rr.paed.u5.ar$k

rr.lc <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("lung cancer"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("lung cancer"),]$author) 
rr.lca <- predict(rr.lc, transf = exp, digit=2)
rr.lca <- as.data.frame(rr.lca)
rr.lca$subgroup <- "lung cancer"
rr.lca$i2 <- round(rr.lc$I2, digits=1)
rr.lca$estimates <- rr.lc$k

rr.t <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("tuberculosis"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("tuberculosis"),]$author) 
rr.tb <- predict(rr.t, transf = exp, digit=2)
rr.tb <- as.data.frame(rr.tb)
rr.tb$subgroup <- "tb"
rr.tb$i2 <- round(rr.t$I2, digits=1)
rr.tb$estimates <- rr.t$k

rr.resd <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("respiratory disease"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("respiratory disease"),]$author) 
rr.resdis <- predict(rr.resd, transf = exp, digit=2)
rr.resdis <- as.data.frame(rr.resdis)
rr.resdis$subgroup <- "respiratory disease"
rr.resdis$i2 <- round(rr.resd$I2, digits=1)
rr.resdis$estimates <- rr.resd$k

resp.pooled <- rbind (rr.asthma, rr.copd, rr.ari, rr.adult.ari, rr.paed.ari,rr.paed.5t16.ari, rr.paed.u5.ari, rr.lca, rr.tb, rr.resdis, rr.resp)

#### cardiovascular disease ####

rr.cv <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("stroke", "ihd", "cardiovascular disease"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("stroke", "ihd", "cardiovascular disease"),]$author) 
rr.cvd <- predict(rr.cv, transf = exp, digit=2)
rr.cvd <- as.data.frame(rr.cvd)
rr.cvd$subgroup <- "cvd overall"
rr.cvd$i2 <- round(rr.cv$I2, digits=1)
rr.cvd$estimates <- rr.cv$k

rr.strok <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("stroke"),],method = "ML", slab=) 
rr.stroke <- predict(rr.strok, transf = exp, digit=2)
rr.stroke <- as.data.frame(rr.stroke)
rr.stroke$subgroup <- "stroke"
rr.stroke$i2 <- round(rr.strok$I2, digits=1)
rr.stroke$estimates <- rr.strok$k

rr.ih <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("ihd"),],method = "ML", slab=) 
rr.ihd <- predict(rr.ih, transf = exp, digit=2)
rr.ihd <- as.data.frame(rr.ihd)
rr.ihd$subgroup <- "ihd"
rr.ihd$i2 <- round(rr.ih$I2, digits=1)
rr.ihd$estimates <- rr.ih$k

rr.cardi <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("cardiovascular disease"),],method = "ML", slab=) 
rr.cardis <- predict(rr.cardi, transf = exp, digit=2)
rr.cardis <- as.data.frame(rr.cardis)
rr.cardis$subgroup <- "cardiovascular disease"
rr.cardis$i2 <- round(rr.cardi$I2, digits=1)
rr.cardis$estimates <- rr.cardi$k

cvd.pooled <- rbind (rr.stroke, rr.ihd, rr.cardis, rr.cvd)

#### low birth weight ####
rr.lb <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("lbw"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("lbw"),]$author) 
rr.lbw <- predict(rr.lb, transf = exp, digit=2)
rr.lbw <- as.data.frame(rr.lbw)
rr.lbw$subgroup <- "lbw"
rr.lbw$i2 <- round(rr.lb$I2, digits=1)
rr.lbw$estimates <- rr.lb$k

#### stillbirth ####
rr.stil <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("stillbirth"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("stillbirth"),]$author) 
rr.still <- predict(rr.stil, transf = exp, digit=2)
rr.still <- as.data.frame(rr.still)
rr.still$subgroup <- "stillbirth"
rr.still$i2 <- round(rr.stil$I2, digits=1)
rr.still$estimates <- rr.stil$k

#### under five mortality ####

rr.ufmor <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("under five death", "infant mortality", "neonatal death", "post neonatal mortality"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("under five death", "infant mortality", "neonatal death", "post neonatal mortality"),]$author) 
rr.ufmort <- predict(rr.ufmor, transf = exp, digit=2)
rr.ufmort <- as.data.frame(rr.ufmort)
rr.ufmort$subgroup <- "under five death"
rr.ufmort$i2 <- round(rr.ufmor$I2, digits=1)
rr.ufmort$estimates <- rr.ufmor$k

#### cvd death ####
rr.cvdmor <- rma(yi = lgrr,sei = se, data = dat.ratio.mort[dat.ratio.mort$outcome.recode %in% c("cv death"),],method = "ML", slab=dat.ratio.mort[dat.ratio.mort$outcome.recode %in% c("cv death"),]$author) 
rr.cvdmort <- predict(rr.cvdmor, transf = exp, digit=2)
rr.cvdmort <- as.data.frame(rr.cvdmort)
rr.cvdmort$subgroup <- "cardio death"
rr.cvdmort$i2 <- round(rr.cvdmor$I2, digits=1)
rr.cvdmort$estimates <- rr.cvdmor$k

#### resp death ####
rr.respmor <- rma(yi = lgrr,sei = se, data = dat.ratio.mort[dat.ratio.mort$outcome.recode %in% c("resp death"),],method = "ML", slab=dat.ratio.mort[dat.ratio.mort$outcome.recode %in% c("resp death"),]$author) 
rr.respmort <- predict(rr.respmor, transf = exp, digit=2)
rr.respmort <- as.data.frame(rr.respmort)
rr.respmort$subgroup <- "resp death"
rr.respmort$i2 <- round(rr.respmor$I2, digits=1)
rr.respmort$estimates <- rr.respmor$k

##### all cause death ####
rr.allmor <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("all cause death"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("all cause death"),]$author) 
rr.allmort <- predict(rr.allmor, transf = exp, digit=2)
rr.allmort <- as.data.frame(rr.allmort)
rr.allmort$subgroup <- "all cause death"
rr.allmort$i2 <- round(rr.allmor$I2, digits=1)
rr.allmort$estimates <- rr.allmor$k

mort.pooled <- rbind (rr.still, rr.ufmort, rr.cvdmort, rr.respmort, rr.allmort)

all.pooled <- rbind(resp.pooled, cvd.pooled, rr.lbw, mort.pooled)
all.pooled <- all.pooled %>% mutate_at(vars(pred, ci.lb, ci.ub), funs(round(., 2)))
all.pooled$rr <-paste0(all.pooled$pred," ","(",all.pooled$ci.lb,"-",all.pooled$ci.ub,")") 
all.pooled <- all.pooled %>% 
  select(subgroup, estimates,pred, ci.lb, ci.ub,rr,i2)

write.csv(all.pooled, "results.v2/all.pooled.csv")


#### place beta and se into dataframe for burden analysis ####
endpoints <- c("asthma", "copd", "ari", "lungca", "tb", "ihd", "stroke")

rr.sim.para <- data.frame(endpoints = endpoints, 
                 rr.beta = c(rr.asth$b, rr.cop$b, rr.ar$b, rr.lc$b, rr.t$b, rr.ih$b, rr.strok$b),
                 rr.se = c(rr.asth$se, rr.cop$se, rr.ar$se, rr.lc$se, rr.t$se, rr.ih$se, rr.strok$se))

write.csv(rr.sim.para, "data_derived.v2/rr.sim.para.csv")


#### place beta and se into dataframe for death analysis ####
endpoints <- c("resp death", "cv death", "ufive")

rr.sim.para.death <- data.frame(endpoints = endpoints, 
                          rr.beta = c(rr.respmor$b, rr.cvdmor$b, rr.ufmor$b),
                          rr.se = c(rr.respmor$se, rr.cvdmor$se, rr.ufmor$se))

write.csv(rr.sim.para.death, "data_derived.v2/rr.sim.para.death.csv")


#### stratified by gender ####
dat.adult <- dat %>% 
  filter(adult.paediatric.all=="adult" & overall.estimate %in% c(1,2) & analysis.code %in% c(1,2,5) & outcome.type!="symptoms")

dat.women <- dat.adult %>% filter(subgroup.gender=="women")
dat.women <- subset(dat.women,!is.na(dat.women$ratio.ll))

dat.women <- dat.women %>% 
  select("study.id", "author", "pub.year", "overall.estimate", "analysis.code", "subgroup.id", "subgroup.gender", "subgroup.exposure", "subgroup.outcome", "subgroup.other", "interventional.study", "intervention.type", "exposed.group", "control.group", "outcome.recode", "outcome.type", "risk.measure", "ratio.ce", "ratio.ll", "ratio.ul", "lgrr", "lgll", "lgul", "se", "exposed.event", "exposed.health", "control.event","control.health")


dat.men <- dat.adult %>% filter(subgroup.gender=="men")
dat.men <- subset(dat.men,!is.na(dat.men$ratio.ll))

dat.men <- dat.men %>% 
  select("study.id", "author", "pub.year", "overall.estimate", "analysis.code", "subgroup.id", "subgroup.gender", "subgroup.exposure", "subgroup.outcome", "subgroup.other", "interventional.study", "intervention.type", "exposed.group", "control.group", "outcome.recode", "outcome.type", "risk.measure", "ratio.ce", "ratio.ll", "ratio.ul", "lgrr", "lgll", "lgul", "se", "exposed.event", "exposed.health", "control.event","control.health")

outcomes.women<-as.data.frame(table(dat.women$outcome.recode))
outcomes.men<-as.data.frame(table(dat.men$outcome.recode))
outcomes.gender <- left_join (outcomes.women, outcomes.men, by = "Var1")

#### women ####
# women- respiratory disease #
rr.women.res <- rma(yi = lgrr,sei = se, data = dat.women[dat.women$outcome.recode %in% c("asthma", "copd", "ari", "lung cancer", "tuberculosis", "respiratory disease"),],method = "ML")
rr.women.resp <- predict(rr.women.res, transf = exp, digit=2)
rr.women.resp <- as.data.frame(rr.women.resp)
rr.women.resp$subgroup <- "respiratory overall"
rr.women.resp$i2 <- round(rr.women.res$I2, digits=1)
rr.women.resp$estimates <- rr.women.res$k

rr.women.asth <- rma(yi = lgrr,sei = se, data = dat.women[dat.women$outcome.recode %in% c("asthma"),],method = "ML")
rr.women.asthma <- predict(rr.women.asth, transf = exp, digit=2)
rr.women.asthma <- as.data.frame(rr.women.asthma)
rr.women.asthma$subgroup <- "asthma"
rr.women.asthma$i2 <- round(rr.women.asth$I2, digits=1)
rr.women.asthma$estimates <- rr.women.asth$k

rr.women.cop <- rma(yi = lgrr,sei = se, data = dat.women[dat.women$outcome.recode %in% c("copd"),],method = "ML")
rr.women.copd <- predict(rr.women.cop, transf = exp, digit=2)
rr.women.copd <- as.data.frame(rr.women.copd)
rr.women.copd$subgroup <- "copd"
rr.women.copd$i2 <- round(rr.women.cop$I2, digits=1)
rr.women.copd$estimates <- rr.women.cop$k

rr.women.ar <- rma(yi = lgrr,sei = se, data = dat.women[dat.women$outcome.recode %in% c("ari"),],method = "ML")
rr.women.ari <- predict(rr.women.ar, transf = exp, digit=2)
rr.women.ari <- as.data.frame(rr.women.ari)
rr.women.ari$subgroup <- "ari"
rr.women.ari$i2 <- round(rr.women.ar$I2, digits=1)
rr.women.ari$estimates <- rr.women.ar$k

rr.women.lc <- rma(yi = lgrr,sei = se, data = dat.women[dat.women$outcome.recode %in% c("lung cancer"),],method = "ML")
rr.women.lca <- predict(rr.women.lc, transf = exp, digit=2)
rr.women.lca <- as.data.frame(rr.women.lca)
rr.women.lca$subgroup <- "lung cancer"
rr.women.lca$i2 <- round(rr.women.lc$I2, digits=1)
rr.women.lca$estimates <- rr.women.lc$k

rr.women.t <- rma(yi = lgrr,sei = se, data = dat.women[dat.women$outcome.recode %in% c("tuberculosis"),],method = "ML")
rr.women.tb <- predict(rr.women.t, transf = exp, digit=2)
rr.women.tb <- as.data.frame(rr.women.tb)
rr.women.tb$subgroup <- "tb"
rr.women.tb$i2 <- round(rr.women.t$I2, digits=1)
rr.women.tb$estimates <- rr.women.t$k

rr.women.resd <- rma(yi = lgrr,sei = se, data = dat.women[dat.women$outcome.recode %in% c("respiratory disease"),],method = "ML")
rr.women.resdis <- predict(rr.women.resd, transf = exp, digit=2)
rr.women.resdis <- as.data.frame(rr.women.resdis)
rr.women.resdis$subgroup <- "respiratory disease"
rr.women.resdis$i2 <- round(rr.women.resd$I2, digits=1)
rr.women.resdis$estimates <- rr.women.resd$k

resp.women.pooled <- rbind (rr.women.asthma, rr.women.copd, rr.women.ari, rr.women.lca, rr.women.tb, rr.women.resdis, rr.women.resp)

# women- cardiovascular disease #
rr.women.cv <- rma(yi = lgrr,sei = se, data = dat.women[dat.women$outcome.recode %in% c("stroke", "ihd", "cardiovascular disease"),],method = "ML")
rr.women.cvd <- predict(rr.women.cv, transf = exp, digit=2)
rr.women.cvd <- as.data.frame(rr.women.cvd)
rr.women.cvd$subgroup <- "cvd overall"
rr.women.cvd$i2 <- round(rr.women.cv$I2, digits=1)
rr.women.cvd$estimates <- rr.women.cv$k

rr.women.strok <- rma(yi = lgrr,sei = se, data = dat.women[dat.women$outcome.recode %in% c("stroke"),],method = "ML")
rr.women.stroke <- predict(rr.women.strok, transf = exp, digit=2)
rr.women.stroke <- as.data.frame(rr.women.stroke)
rr.women.stroke$subgroup <- "stroke"
rr.women.stroke$i2 <- round(rr.women.strok$I2, digits=1)
rr.women.stroke$estimates <- rr.women.strok$k

rr.women.ih <- rma(yi = lgrr,sei = se, data = dat.women[dat.women$outcome.recode %in% c("ihd"),],method = "ML")
rr.women.ihd <- predict(rr.women.ih, transf = exp, digit=2)
rr.women.ihd <- as.data.frame(rr.women.ihd)
rr.women.ihd$subgroup <- "ihd"
rr.women.ihd$i2 <- round(rr.women.ih$I2, digits=1)
rr.women.ihd$estimates <- rr.women.ih$k

cvd.women.pooled <- rbind (rr.women.stroke, rr.women.ihd, rr.women.cvd)

women.pooled <- rbind(resp.women.pooled, cvd.women.pooled)
women.pooled <- women.pooled %>% mutate_at(vars(pred, ci.lb, ci.ub), funs(round(., 2)))
women.pooled$rr <-paste0(women.pooled$pred," ","(",women.pooled$ci.lb,"-",women.pooled$ci.ub,")") 
women.pooled <- women.pooled %>% 
  select(subgroup, estimates,pred, ci.lb, ci.ub,rr,i2)

write.csv(women.pooled, "results.v2/women.pooled.csv")



#### men- insufficient estimates for TB and ARI ####
# men- respiratory disease #
rr.men.res <- rma(yi = lgrr,sei = se, data = dat.men[dat.men$outcome.recode %in% c("asthma", "copd", "ari", "lung cancer", "tuberculosis", "respiratory disease"),],method = "ML")
rr.men.resp <- predict(rr.men.res, transf = exp, digit=2)
rr.men.resp <- as.data.frame(rr.men.resp)
rr.men.resp$subgroup <- "respiratory overall"
rr.men.resp$i2 <- round(rr.men.res$I2, digits=1)
rr.men.resp$estimates <- rr.men.res$k

rr.men.asth <- rma(yi = lgrr,sei = se, data = dat.men[dat.men$outcome.recode %in% c("asthma"),],method = "ML")
rr.men.asthma <- predict(rr.men.asth, transf = exp, digit=2)
rr.men.asthma <- as.data.frame(rr.men.asthma)
rr.men.asthma$subgroup <- "asthma"
rr.men.asthma$i2 <- round(rr.men.asth$I2, digits=1)
rr.men.asthma$estimates <- rr.men.asth$k

rr.men.cop <- rma(yi = lgrr,sei = se, data = dat.men[dat.men$outcome.recode %in% c("copd"),],method = "ML")
rr.men.copd <- predict(rr.men.cop, transf = exp, digit=2)
rr.men.copd <- as.data.frame(rr.men.copd)
rr.men.copd$subgroup <- "copd"
rr.men.copd$i2 <- round(rr.men.cop$I2, digits=1)
rr.men.copd$estimates <- rr.men.cop$k

rr.men.lc <- rma(yi = lgrr,sei = se, data = dat.men[dat.men$outcome.recode %in% c("lung cancer"),],method = "ML")
rr.men.lca <- predict(rr.men.lc, transf = exp, digit=2)
rr.men.lca <- as.data.frame(rr.men.lca)
rr.men.lca$subgroup <- "lung cancer"
rr.men.lca$i2 <- round(rr.men.lc$I2, digits=1)
rr.men.lca$estimates <- rr.men.lc$k

rr.men.resd <- rma(yi = lgrr,sei = se, data = dat.men[dat.men$outcome.recode %in% c("respiratory disease"),],method = "ML")
rr.men.resdis <- predict(rr.men.resd, transf = exp, digit=2)
rr.men.resdis <- as.data.frame(rr.men.resdis)
rr.men.resdis$subgroup <- "respiratory disease"
rr.men.resdis$i2 <- round(rr.men.resd$I2, digits=1)
rr.men.resdis$estimates <- rr.men.resd$k

resp.men.pooled <- rbind (rr.men.asthma, rr.men.copd, rr.men.lca, rr.men.resdis, rr.men.resp)

# men- cardiovascular disease #
rr.men.cv <- rma(yi = lgrr,sei = se, data = dat.men[dat.men$outcome.recode %in% c("stroke", "ihd", "cardiovascular disease"),],method = "ML")
rr.men.cvd <- predict(rr.men.cv, transf = exp, digit=2)
rr.men.cvd <- as.data.frame(rr.men.cvd)
rr.men.cvd$subgroup <- "cvd overall"
rr.men.cvd$i2 <- round(rr.men.cv$I2, digits=1)
rr.men.cvd$estimates <- rr.men.cv$k

rr.men.strok <- rma(yi = lgrr,sei = se, data = dat.men[dat.men$outcome.recode %in% c("stroke"),],method = "ML")
rr.men.stroke <- predict(rr.men.strok, transf = exp, digit=2)
rr.men.stroke <- as.data.frame(rr.men.stroke)
rr.men.stroke$subgroup <- "stroke"
rr.men.stroke$i2 <- round(rr.men.strok$I2, digits=1)
rr.men.stroke$estimates <- rr.men.strok$k

rr.men.ih <- rma(yi = lgrr,sei = se, data = dat.men[dat.men$outcome.recode %in% c("ihd"),],method = "ML")
rr.men.ihd <- predict(rr.men.ih, transf = exp, digit=2)
rr.men.ihd <- as.data.frame(rr.men.ihd)
rr.men.ihd$subgroup <- "ihd"
rr.men.ihd$i2 <- round(rr.men.ih$I2, digits=1)
rr.men.ihd$estimates <- rr.men.ih$k

cvd.men.pooled <- rbind (rr.men.stroke, rr.men.ihd, rr.men.cvd)

men.pooled <- rbind(resp.men.pooled, cvd.men.pooled)
men.pooled <- men.pooled %>% mutate_at(vars(pred, ci.lb, ci.ub), funs(round(., 2)))
men.pooled$rr <-paste0(men.pooled$pred," ","(",men.pooled$ci.lb,"-",men.pooled$ci.ub,")") 
men.pooled <- men.pooled %>% 
  select(subgroup, estimates,pred, ci.lb, ci.ub,rr,i2)

write.csv(men.pooled, "results.v2/men.pooled.csv")


#### paediatric only ####
dat.paed <- dat %>% 
  filter(adult.paediatric.all=="paed" & overall.estimate %in% c(1,2) & analysis.code %in% c(1,2,5) & outcome.type!="symptoms")

outcomes.paed<-as.data.frame(table(dat.paed$outcome.recode))

rr.paed.res <- rma(yi = lgrr,sei = se, data = dat.paed[dat.paed$outcome.recode %in% c("asthma", "copd", "ari", "lung cancer", "tuberculosis", "respiratory disease"),],method = "ML")
rr.paed.resp <- predict(rr.paed.res, transf = exp, digit=2)
rr.paed.resp <- as.data.frame(rr.paed.resp)
rr.paed.resp$subgroup <- "respiratory overall"
rr.paed.resp$i2 <- round(rr.paed.res$I2, digits=1)
rr.paed.resp$estimates <- rr.paed.res$k

rr.paed.asth <- rma(yi = lgrr,sei = se, data = dat.paed[dat.paed$outcome.recode %in% c("asthma"),],method = "ML")
rr.paed.asthma <- predict(rr.paed.asth, transf = exp, digit=2)
rr.paed.asthma <- as.data.frame(rr.paed.asthma)
rr.paed.asthma$subgroup <- "asthma"
rr.paed.asthma$i2 <- round(rr.paed.asth$I2, digits=1)
rr.paed.asthma$estimates <- rr.paed.asth$k

rr.paed.resd <- rma(yi = lgrr,sei = se, data = dat.paed[dat.paed$outcome.recode %in% c("respiratory disease"),],method = "ML")
rr.paed.resdis <- predict(rr.paed.resd, transf = exp, digit=2)
rr.paed.resdis <- as.data.frame(rr.paed.resdis)
rr.paed.resdis$subgroup <- "respiratory disease"
rr.paed.resdis$i2 <- round(rr.paed.resd$I2, digits=1)
rr.paed.resdis$estimates <- rr.paed.resd$k

paed.pooled <- rbind (rr.paed.asthma, rr.paed.ari, rr.paed.resp)
paed.pooled <- paed.pooled %>% mutate_at(vars(pred, ci.lb, ci.ub), funs(round(., 2)))
paed.pooled$rr <-paste0(paed.pooled$pred," ","(",paed.pooled$ci.lb,"-",paed.pooled$ci.ub,")") 
paed.pooled <- paed.pooled %>% 
  select(subgroup,estimates,rr,i2)

write.csv(paed.pooled, "results.v2/paed.pooled.csv")


#### place beta and se into dataframe for paed analysis ####
endpoints <- c("asthma", "ari", "ufive")

rr.sim.para.paed <- data.frame(endpoints = endpoints, 
                                rr.beta = c(rr.paed.asth$b, rr.paed.ar$b, rr.ufmor$b),
                                rr.se = c(rr.paed.asth$se, rr.paed.ar$se, rr.ufmor$se))

write.csv(rr.sim.para.paed, "data_derived.v2/rr.sim.para.paed.csv")


# QQ plots ----------------------------------------------------------------

pdf("results.v2/qqplot_diagnosis.pdf", height=8, width=12)

par(mfrow=c(3,3), mar=c(4,4,4,4))
qqnorm(rr.asth, main = "Asthma", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.cop, main = "COPD", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.adult.ar, main = "ARI (adults)", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.paed.ar, main = "ARI (paediatric)", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.lc, main = "Lung cancer", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)
       
qqnorm(rr.t, main = "Tuberculosis", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.strok, main = "Cerebrovascular disease", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.ih, main = "Ischaemic heart disease", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.lb, main = "Low birth weight", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

par(op)

dev.off()


pdf("results.v2/qqplot_mortality.pdf", height=8, width=12)

par(mfrow=c(2,3), mar=c(4,4,4,4))

qqnorm(rr.allmor, main = "All-cause mortality", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.cvdmor, main = "Cardiovascular mortality", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.respmor, main = "Respiratory mortality", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.ufmor, main = "Under-five mortality", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

qqnorm(rr.stil, main = "Stillbirth", xlab="", ylab="")
mtext(text = "Theoretical quantiles",side = 1, line = 2, cex=0.7)
mtext(text = "Sample quantiles",side = 2, line = 2, cex=0.7)

par(op)

dev.off()


##============================================
## plot forestplots
##============================================

dat.ratio <- dat.ratio %>% 
  mutate_at(vars(ratio.ce, ratio.ll, ratio.ul),
  funs(formatC(as.numeric(.), digits = 2, format = "f")))

dat.ratio$rr <- paste0(dat.ratio$ratio.ce," [",dat.ratio$ratio.ll,"-",dat.ratio$ratio.ul,"]")

dat.ratio <- dat.ratio %>% 
  mutate_at(vars("exposed.event", "exposed.health", "control.event","control.health"),
            funs(round(.,0)))

dat.ratio$events <- paste0(dat.ratio$exposed.event, "/",
                           (dat.ratio$exposed.event+dat.ratio$exposed.health), " vs ", dat.ratio$control.event, "/",
                           (dat.ratio$control.event+dat.ratio$control.health))

dat.ratio$events <- ifelse(is.na(dat.ratio$exposed.event), 
                           "not reported", dat.ratio$events)


## forest plot function ##

forest_plot <- function(data, endpoint, age_group, subtitle) {
  
  # data=dat.ratio
  # endpoint="copd"
  # age_group=c("adult", "both", "paed") 
  # subtitle= "all"
  
  # set population
  data_temp <- data %>% 
    filter(outcome.recode %in% endpoint & 
             adult.paediatric.all %in% age_group)

  # perform meta-analysis
  rr <- rma(yi = lgrr,
            sei = se, 
            data = data_temp,
            method = "ML", 
            slab=data_temp$author)
  rr_temp <- predict(rr, transf = exp, digits=2)
  rr_temp <- as.data.frame(rr_temp)
  
  cairo_pdf(glue("results.v2/{endpoint}_{subtitle}.pdf"), 
            height=15.5, 
            width=11, 
            onefile = FALSE)
  
  # plot forest plot
  par(mar=c(4,4,1,4), font=1)
  
  forest(rr, 
         xlim=c(-4, 4.5), 
         at=log(c(0.25,0.50, 0.75, 1,2,3,4,5,6)), 
         atransf=exp,
         ilab=cbind(data_temp$pub.year, data_temp$events ,data_temp$rr),
         ilab.xpos=c(-2, 2.5, 3.8),
         cex=0.75, 
         ylim=c(-1, nrow(data_temp)+3),
         annotate=FALSE,
         order=rev(order(data_temp$pub.year)),
         rows=c(1:nrow(data_temp)),
         xlab="Risk Ratio", mlab="",
         col = "blue", 
         border="blue")
  
  text(-4, -1, pos=4, cex=0.70, 
       bquote(paste("RE Model (df = ", .(rr$k - rr$p),"; ", 
                    I^2, " = ",.(formatC(rr$I2, digits=1, format="f")), "%),  ", 
                    "Pooled risk ratio = ", .(formatC(rr_temp$pred, digits=2, format="f")), 
                    " [95% CI ",.(formatC(rr_temp$ci.lb, digits=2, format="f"))," - ",.(formatC(rr_temp$ci.ub, digits=2, format="f")),
                    "]")))
  
  op <- par(cex=0.70, font=4)
  
  par(font=2)
  
  text(-4,nrow(data_temp)+2, "Author(s)",  pos=4)
  text(-2,nrow(data_temp)+2,"Year")
  text(2.5,nrow(data_temp)+2.2, "Number of events \n in exposed vs unexposed")
  text(3.8,nrow(data_temp)+2, "Risk Ratio [95% CI]")
  
  par(op)
  
  dev.off()
  
  rm(data_temp, rr_temp)
  
}


# plot forest plots
forest_plot(data=dat.ratio, 
            endpoint="copd", 
            age_group=c("adult", "both", "paed"), 
            subtitle= "all")

forest_plot(data=dat.ratio, 
            endpoint="asthma", 
            age_group=c("adult", "both", "paed"), 
            subtitle= "all")

forest_plot(data=dat.ratio, 
            endpoint="ari", 
            age_group=c("paed"), 
            subtitle= "paed")

forest_plot(data=dat.ratio, 
            endpoint="ari", 
            age_group=c("adult", "both"), 
            subtitle= "adult")

forest_plot(data=dat.ratio, 
            endpoint="lung cancer", 
            age_group=c("adult", "both", "paed"), 
            subtitle= "all")

forest_plot(data=dat.ratio, 
            endpoint="tuberculosis", 
            age_group=c("adult", "both", "paed"), 
            subtitle= "all")

forest_plot(data=dat.ratio, 
            endpoint="lbw", 
            age_group=c("adult", "both", "paed"), 
            subtitle= "all")

#### Plot for Cardiovascular disease ####

dat.ratio.cvd <- dat.ratio[dat.ratio$outcome.recode %in% c("stroke", "ihd", "cardiovascular disease"),]

pdf("results.v2/cvd.pdf",
    height=8,
    width=12,
    onefile=FALSE)

par(mar=c(4,4,1,2), font=1)

forest(rr.cv, xlim=c(-4, 4.5), at=log(c(0.25,0.5, 0.75, 1,2,3,4,5,6)), atransf=exp,
       ilab=cbind(dat.ratio.cvd$pub.year,dat.ratio.cvd$events,dat.ratio.cvd$rr),
       ilab.xpos=c(-2,2.5, 3.8),cex=0.75, ylim=c(-1, 53),
       order=rev(order(dat.ratio.cvd$outcome.recode,dat.ratio.cvd$pub.year)),
       rows=c(3:15,21:33,39:49),
       annotate=FALSE,
       xlab="Risk Ratio", mlab="",
       col = "red", border="red")

### add text with Q-value, dfs, p-value, and I^2 statistic
text(-4, -1, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.cv$k - rr.cv$p),
                                           "; ", I^2, " = ",.(formatC(rr.cv$I2, digits=1, format="f")), "%),  ", "Overall risk ratio = ", .(formatC(rr.cvd$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.cvd$ci.lb, digits=2, format="f")),
                                           " - ",
                                           .(formatC(rr.cvd$ci.ub, digits=2, format="f")),
                                           "]")))

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)

### add text for the subgroups
text(-4, c(50,34,16), pos=4, c("Cardiovascular events",
                               "Ischemic heart disease",
                               "Cerebrovascular disease"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-4,52, "Author(s)",  pos=4)
text(-2,52,"Year")
text(2.5,52.2, "Number of events \n in exposed vs unexposed")
text(3.8,52, "Risk Ratio [95% CI]")

### set par back to the original settings
par(op)


### add summary polygons for the three subgroups
addpoly(rr.cardi, row=37.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate=FALSE)
addpoly(rr.ih, row= 19.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate=FALSE)
addpoly(rr.strok, row= 1.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate=FALSE)

### add text with Q-value, dfs, p-value, and I^2 statistic for subgroups

text(-4, 37.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.cardi$k - rr.cardi$p),
                                             "; ", I^2, " = ",.(formatC(rr.cardi$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.cardis$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.cardis$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.cardis$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 19.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.ih$k - rr.ih$p),
                                             "; ", I^2, " = ",.(formatC(rr.ih$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.ihd$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.ihd$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.ihd$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 1.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.strok$k - rr.strok$p),
                                            "; ", I^2, " = ",.(formatC(rr.strok$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.stroke$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.stroke$ci.lb, digits=2, format="f")),
                                            " - ",
                                            .(formatC(rr.stroke$ci.ub, digits=2, format="f")),
                                            "]")))

dev.off()


#### Plot for mortality ####
dat.ratio.death <- dat.ratio.mort %>% 
  mutate(outcome.recode = case_when(
    outcome.recode %in% c("stillbirth") ~ "a.stillbirth",
    outcome.recode %in% c("under five death", "infant mortality", "neonatal death", "post neonatal mortality") ~ "b.under five death",
    outcome.recode %in% c("cv death") ~ "c.cv death",
    outcome.recode %in% c("resp death") ~ "d.resp death",
    outcome.recode %in% c("all cause death") ~ "e.all cause death",
    TRUE ~ as.character(outcome.recode)
  ))


outcomes.mort<-as.data.frame(table(dat.ratio.death$outcome.recode))


dat.ratio.death <- dat.ratio.death %>% 
  mutate_at(vars(ratio.ce, ratio.ll, ratio.ul),
            funs(round(.,2)))

dat.ratio.death$rr <- dat.ratio.death$rr <- paste0(dat.ratio.death$ratio.ce," [",dat.ratio.death$ratio.ll,"-",dat.ratio.death$ratio.ul,"]")


dat.ratio.death <- dat.ratio.death %>% 
  mutate_at(vars("exposed.event", "exposed.health", "control.event","control.health"),
            funs(round(.,0)))

dat.ratio.death$events <- paste0(dat.ratio.death$exposed.event, "/",
                           (dat.ratio.death$exposed.event+dat.ratio.death$exposed.health), " vs ", dat.ratio.death$control.event, "/",
                           (dat.ratio.death$control.event+dat.ratio.death$control.health))

dat.ratio.death$events <- ifelse(is.na(dat.ratio.death$exposed.event), 
                           "not reported", dat.ratio.death$events)


## summary plot  for death##

death.pooled <- all.pooled %>%
  mutate(analysis = case_when(
    subgroup %in% c("under five death") ~ "a.under-five death",
    subgroup %in% c("cardio death") ~ "b.cvd death",
    subgroup %in% c("resp death") ~ "c.resp death",
    subgroup %in% c("all cause death") ~ "d.all death",
    TRUE ~ "exclude"
  )) %>%
  filter(analysis!="exclude")


death.pooled <- death.pooled %>%
  mutate(disease = case_when(
    subgroup %in% c("under five death") ~ "Under-five mortality",
    subgroup %in% c("cardio death") ~ "Cardiovascular mortality",
    subgroup %in% c("resp death") ~ "Respiratory mortality",
    subgroup %in% c("all cause death") ~ "All-cause mortality",
    TRUE ~ as.character(subgroup)
  ))

death.pooled$lgrr <- logb(death.pooled$pred)
death.pooled$lgll <- logb(death.pooled$ci.lb)
death.pooled$lgul <- logb(death.pooled$ci.ub)
death.pooled$se <- (death.pooled$lgul-death.pooled$lgll)/(2*1.96)

rr.death.pool <- rma(yi = lgrr,sei = se, data = death.pooled ,method = "ML", slab=death.pooled$disease)

pdf("results.v2/Figure_2.mortality.summary.pdf", 
    height=5,
    width=10)

par(mfrow=c(1,1),mar=c(4,4,1,2), font=1)

forest(rr.death.pool, xlim=c(-4, 1), at=log(c(0.5,1,2,3,4,5)), atransf=exp,
       ilab=cbind(death.pooled$estimates, paste0(death.pooled$i2, "%"), death.pooled$rr),
       ilab.xpos=c(-2.5,-1.9, -1),cex=1, ylim=c(0, 7), psize=1, col="white",
       order=rev(order(death.pooled$analysis)),
       annotate=FALSE,
       addfit=FALSE,
       rows=c(1:4),
       xlab="Risk Ratio", mlab="")

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)

### switch to bold font
par(font=2)

### add column headings to the plot
text(-4,5.5, "Mortality",  pos=4)
text(-2.5,5.5,"No. of \n Estimates")
text(-1.9,5.5, expression("I"^2))
text(-1,5.5, "Pooled risk ratios \n [95% CI]")

### set par back to the original settings
par(op)

dev.off()


# Plot mortality forestplots ----------------------------------------------

forest_plot(data=dat.ratio, 
            endpoint=c("under five death", "infant mortality", "neonatal death", "post neonatal mortality"), 
            age_group=c("paed"), 
            subtitle= "paed")

## plot mortality
dat.ratio.mort <- dat.ratio.mort %>% 
  mutate_at(vars(ratio.ce, ratio.ll, ratio.ul),
            funs(formatC(as.numeric(.), digits = 2, format = "f")))

dat.ratio.mort$rr <- dat.ratio.mort$rr <- paste0(dat.ratio.mort$ratio.ce," [",dat.ratio.mort$ratio.ll,"-",dat.ratio.mort$ratio.ul,"]")

dat.ratio.mort <- dat.ratio.mort %>% 
  mutate_at(vars("exposed.event", "exposed.health", "control.event","control.health"),
            funs(round(.,0)))

dat.ratio.mort$events <- paste0(dat.ratio.mort$exposed.event, "/",
                                 (dat.ratio.mort$exposed.event+dat.ratio.mort$exposed.health), " vs ", dat.ratio.mort$control.event, "/",
                                 (dat.ratio.mort$control.event+dat.ratio.mort$control.health))

dat.ratio.mort$events <- ifelse(is.na(dat.ratio.mort$exposed.event), 
                                 "not reported", dat.ratio.mort$events)


forest_plot(data=dat.ratio.mort, 
            endpoint=c("cv death"), 
            age_group=c("paed", "adult", "both"), 
            subtitle= "all")

forest_plot(data=dat.ratio.mort, 
            endpoint=c("resp death"), 
            age_group=c("paed", "adult", "both"), 
            subtitle= "all")

forest_plot(data=dat.ratio.mort, 
            endpoint=c("all cause death"), 
            age_group=c("paed", "adult", "both"), 
            subtitle= "all")

## summary plot ##

all.pooled <- all.pooled %>% 
  mutate(analysis = case_when(
    subgroup %in% c("asthma", "copd", "ari- adult", "ari- paed", "lung cancer", "tb", "respiratory disease") ~ "a.resp",
    subgroup %in% c("stroke", "ihd", "cardiovascular disease") ~ "b.cvd",
    subgroup %in% c("lbw") ~ "c.lbw",
    subgroup %in% c("stillbirth") ~ "d.stillbirth",
    TRUE ~ "exclude"
  )) %>% 
  filter(analysis!="exclude")
  


all.pooled <- all.pooled %>%
  mutate(disease = case_when(
    subgroup %in% c("asthma")~ "Asthma",
    subgroup %in% c("copd")~ "COPD",
    subgroup %in% c("ari- adult")~ "ARI (adults)",
    subgroup %in% c("ari- paed")~ "ARI (paediatric)",
    subgroup %in% c("lung cancer")~ "Lung cancer",
    subgroup %in% c("tb")~ "Tuberculosis",
    subgroup %in% c("respiratory disease")~ "Respiratory disease",
    subgroup %in% c("stroke")~ "Cerebrovascular disease",
    subgroup %in% c("ihd")~ "Ischemic heart disease",
    subgroup %in% c("cardiovascular disease")~ "Cardiovascular events",
    subgroup %in% c("lbw")~ "Low birth weight",
    subgroup %in% c("stillbirth") ~ "Stillbirth",
    TRUE ~ as.character(subgroup)
  ))

all.pooled$lgrr <- logb(all.pooled$pred)
all.pooled$lgll <- logb(all.pooled$ci.lb)
all.pooled$lgul <- logb(all.pooled$ci.ub)
all.pooled$se <- (all.pooled$lgul-all.pooled$lgll)/(2*1.96)

rr.pool <- rma(yi = lgrr,sei = se, data = all.pooled ,method = "ML", slab=all.pooled$disease)

pdf("results.v2/Figure_2.overall.summary.pdf",
    height=8,
    width=12)

par(mar=c(4,4,1,2), font=1)

forest(rr.pool, xlim=c(-4, 1), at=log(c(0.5,1, 2,3,4,5)), atransf=exp,
       ilab=cbind(all.pooled$estimates, paste0(all.pooled$i2, "%"), all.pooled$rr),
       ilab.xpos=c(-2.5,-1.9, -1),cex=1, ylim=c(0, 20), psize=1, 
       order=rev(order(all.pooled$analysis)),
       rows=c(1:2, 5:7, 10:16),
       addfit=FALSE,
       annotate=FALSE,
       xlab="Risk Ratio", mlab="")

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)

### add text for the subgroups
text(-4, c(17,8,3), pos=4, c("Respiratory diseases",
                                     "Cardiovascular diseases",
                                     "Adverse pregnancy outcomes"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-4,19, "Disease",  pos=4)
text(-2.5,19,"No. of \n Estimates")
text(-1.9,19, expression("I"^2))
text(-1,19, "Pooled risk ratios \n [95% CI]")

### set par back to the original settings
par(op)

dev.off()

## Gender forest plot ##

men.pooled$gender <- "men"
women.pooled$gender <- "women"

men.pooled <- men.pooled %>% filter(!(subgroup %in% c("ari", "tb", "respiratory disease", "respiratory overall", "cvd overall"))) %>% arrange(subgroup)
women.pooled <- women.pooled %>% filter(!(subgroup %in% c("ari", "tb", "respiratory disease", "respiratory overall", "cvd overall"))) %>% arrange(subgroup)

gender.pooled <- rbind(men.pooled, women.pooled)
gender.pooled <- gender.pooled %>% 
  filter(!(subgroup %in% c("ari", "tb", "respiratory disease", "respiratory overall", "cvd overall"))) %>% 
  arrange(subgroup) %>% 
  select(subgroup, gender, estimates, i2, rr)

write.csv(gender.pooled, "results.v2/gender.pooled.csv")

women.pooled <- women.pooled %>%
  mutate(disease = case_when(
    subgroup %in% c("asthma")~ "Asthma",
    subgroup %in% c("copd")~ "COPD",
    subgroup %in% c("ari")~ "ARI (adults)",
    subgroup %in% c("lung cancer")~ "Lung cancer",
    subgroup %in% c("tb")~ "Tuberculosis",
    subgroup %in% c("respiratory disease")~ "Respiratory disease",
    subgroup %in% c("stroke")~ "Cerebrovascular disease",
    subgroup %in% c("ihd")~ "Ischemic heart disease",
    subgroup %in% c("cardiovascular disease")~ "Cardiovascular events",
    subgroup %in% c("lbw")~ "Low birth weight",
    subgroup %in% c("stillbirth") ~ "Stillbirth",
    subgroup %in% c("under five death") ~ "Under five death",
    subgroup %in% c("cardio death") ~ "Cardiovascular death",
    subgroup %in% c("resp death") ~ "Respiratory death",
    subgroup %in% c("all cause death") ~ "All cause death",
    TRUE ~ as.character(subgroup)
  ))


tabletext<- women.pooled$disease

pdf("results.v2/gender.pdf",
    height=8,
    width=12,
    onefile=FALSE)

forestplot(tabletext, 
           legend = c("Men", "Women"),
           legend_args = fpLegend(pos = list(x=.9, y=0.1), 
                                  gp=gpar(col="#CCCCCC", fill="#F9F9F9")),
           boxsize = .1, # We set the box size to better visualize the type
           mean = cbind(men.pooled$pred, women.pooled$pred),
           lower = cbind(men.pooled$ci.lb, women.pooled$ci.lb),
           upper = cbind(men.pooled$ci.ub, women.pooled$ci.ub),
           clip =c(0.5, 3.55),
           col=fpColors(box=c("blue", "darkred")),
           xlog= TRUE,
           xticks = c(0.75, 1, 2, 3, 3.5),
           ci.vertices=TRUE, ci.vertices.height = 0.03,
           lwd.ci=2, zero=1,
           txt_gp= fpTxtGp(xlab= gpar(cex=1), ticks = gpar(cex=1)),
           xlab="Risk Ratio") 
          
dev.off()


