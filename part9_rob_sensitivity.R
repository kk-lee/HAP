
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

rm(list=ls())

# load files
baseline <- read.csv("results.v2/baseline.table.csv")
dat <- readRDS("data_derived.v2/dat.rds")

dat.rob <- baseline %>% 
  select(study.id, overall.rob) %>% 
  left_join(dat,.)

dat.rob$rr <- paste0(round(dat.rob$ratio.ce,2)," [",round(dat.rob$ratio.ll,2),"-",round(dat.rob$ratio.ul,2),"]")

dat.rob <- dat.rob %>% 
  filter(overall.rob %in% c("moderate RoB", "low RoB"))

#==============================================================
# sensitivity analysis in studies at low/moderate risk of bias
#==============================================================

dat.ratio <- dat.rob

dat.ratio <- subset(dat.ratio,!is.na(dat.ratio$ratio.ll))

dat.ratio <- dat.ratio %>% 
  select("study.id", "author", "pub.year", "overall.estimate", "analysis.code", "subgroup.id", "subgroup.gender", "subgroup.exposure", "subgroup.outcome", "subgroup.other", "interventional.study", "intervention.type", "exposed.group", "control.group", "outcome.recode", "outcome.type", "adult.paediatric.all","paed.under5", "risk.measure", "ratio.ce", "ratio.ll", "ratio.ul", "lgrr", "se")

# diagnosis and mortality only
dat.ratio <- dat.ratio %>% 
  filter(overall.estimate==1 & analysis.code %in% c(1,2,5) & outcome.type!="symptoms")

outcomes<-as.data.frame(table(dat.ratio$outcome.recode))


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

resp.pooled <- rbind (rr.asthma, rr.copd, rr.ari, rr.adult.ari, rr.paed.ari, rr.lca, rr.tb, rr.resp)

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

write.csv(all.pooled, "results.v2/rob.all.pooled.csv")


#### place beta and se into dataframe for burden analysis ####
endpoints <- c("asthma", "copd", "ari", "lungca", "tb", "ihd", "stroke")

rr.sim.para <- data.frame(endpoints = endpoints, 
                          rr.beta = c(rr.asth$b, rr.cop$b, rr.ar$b, rr.lc$b, rr.t$b, rr.ih$b, rr.strok$b),
                          rr.se = c(rr.asth$se, rr.cop$se, rr.ar$se, rr.lc$se, rr.t$se, rr.ih$se, rr.strok$se))

write.csv(rr.sim.para, "data_derived.v2/rr.sim.para.rob.csv")


## summary plot ##

all.pooled.burden <- all.pooled %>% 
  mutate(analysis = case_when(
    subgroup %in% c("asthma", "copd", "ari- overall", "lung cancer", "tb", "respiratory disease") ~ "a.resp",
    subgroup %in% c("stroke", "ihd") ~ "b.cvd",
    subgroup %in% c("lbw") ~ "c.lbw",
    subgroup %in% c("stillbirth") ~ "d.stillbirth",
    TRUE ~ "exclude"
  )) %>% 
  filter(analysis!="exclude")



all.pooled.burden <- all.pooled.burden %>%
  mutate(disease = case_when(
    subgroup %in% c("asthma")~ "Asthma",
    subgroup %in% c("copd")~ "COPD",
    subgroup %in% c("ari- overall")~ "ARI",
    subgroup %in% c("lung cancer")~ "Lung cancer",
    subgroup %in% c("tb")~ "Tuberculosis",
    subgroup %in% c("stroke")~ "Cerebrovascular disease",
    subgroup %in% c("ihd")~ "Ischemic heart disease",
    subgroup %in% c("lbw")~ "Low birth weight",
    subgroup %in% c("stillbirth") ~ "Stillbirth",
    TRUE ~ as.character(subgroup)
  ))

all.pooled.burden$lgrr <- logb(all.pooled.burden$pred)
all.pooled.burden$lgll <- logb(all.pooled.burden$ci.lb)
all.pooled.burden$lgul <- logb(all.pooled.burden$ci.ub)
all.pooled.burden$se <- (all.pooled.burden$lgul-all.pooled.burden$lgll)/(2*1.96)

rr.pool <- rma(yi = lgrr,sei = se, data = all.pooled.burden ,method = "ML", slab=all.pooled.burden$disease)


pdf(file="results.v2/RoB_sensitivity.pdf", height=8, width=12)

par(mar=c(4,4,1,2), font=1)

forest(rr.pool, xlim=c(-4, 1), at=log(c(0.5,1,2,3,4,5)), atransf=exp,
       ilab=cbind(all.pooled.burden$estimates, paste0(all.pooled.burden$i2, "%"), all.pooled.burden$rr),
       ilab.xpos=c(-2.5,-1.9, -1),cex=1, ylim=c(0, 17), psize=1, 
       order=rev(order(all.pooled.burden$analysis)),
       rows=c(1:2, 5:6, 9:13),
       annotate=FALSE,
       addfit=FALSE,
       xlab="Risk Ratio", mlab="")

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)

### add text for the subgroups
text(-4, c(14,7,3), pos=4, c("Respiratory diseases",
                             "Cardiovascular diseases",
                             "Adverse pregnancy outcomes"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-4,16, "Disease",  pos=4)
text(-2.5,16,"No. of \n Estimates")
text(-1.9,16, expression("I"^2))
text(-1,16, "Pooled risk ratios \n [95% CI]")

### set par back to the original settings
par(op)

dev.off()



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


pdf(file="results.v2/RoB_death_sensitivity.pdf", height=8, width=12)

par(mar=c(4,4,1,2), font=1)

forest(rr.death.pool, xlim=c(-4, 1), at=log(c(0.5,1,2,3)), atransf=exp,
       ilab=cbind(death.pooled$estimates, paste0(death.pooled$i2, "%"), death.pooled$rr),
       ilab.xpos=c(-2.5,-1.9, -1),cex=1, ylim=c(0, 7), psize=1, col="white",
       order=rev(order(death.pooled$analysis)),
       rows=c(1:4),
       annotate=FALSE,
       addfit=FALSE,
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







