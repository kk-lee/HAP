
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

###====================================
### Sensitivity analysis for prospective studies only
###====================================


dat <- readRDS("data_derived.v2/dat.rds")
dat$rr <- paste0(round(dat$ratio.ce,2)," [",round(dat$ratio.ll,2),"-",round(dat$ratio.ul,2),"]")

dat.baseline.design <- read.csv("data_raw/iap.baseline.v8.csv") %>% 
  select(study.id, study.design)

dat <- left_join(dat, dat.baseline.design) %>% 
  filter(study.design %in% c("cohort study", "cohort study ", "prospective cohort study", "rct", "retrospective cohort study"))


##============================================
## calculate meta-estimates for risk ratio of main outcomes
##============================================

dat.ratio <- dat

dat.ratio <- subset(dat.ratio,!is.na(dat.ratio$ratio.ll))

dat.ratio <- dat.ratio %>% 
  select("study.id", "author", "pub.year", "overall.estimate", "analysis.code", "subgroup.id", "subgroup.gender", "subgroup.exposure", "subgroup.outcome", "subgroup.other", "interventional.study", "intervention.type", "exposed.group", "control.group", "outcome.recode", "outcome.type", "adult.paediatric.all","paed.under5", "risk.measure", "ratio.ce", "ratio.ll", "ratio.ul", "lgrr", "se")

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


rr.lc <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("lung cancer"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("lung cancer"),]$author) 
rr.lca <- predict(rr.lc, transf = exp, digit=2)
rr.lca <- as.data.frame(rr.lca)
rr.lca$subgroup <- "lung cancer"
rr.lca$i2 <- round(rr.lc$I2, digits=1)
rr.lca$estimates <- rr.lc$k


rr.resd <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome.recode %in% c("respiratory disease"),],method = "ML", slab=dat.ratio[dat.ratio$outcome.recode %in% c("respiratory disease"),]$author) 
rr.resdis <- predict(rr.resd, transf = exp, digit=2)
rr.resdis <- as.data.frame(rr.resdis)
rr.resdis$subgroup <- "respiratory disease"
rr.resdis$i2 <- round(rr.resd$I2, digits=1)
rr.resdis$estimates <- rr.resd$k

resp.pooled <- rbind (rr.asthma, rr.copd, rr.ari, rr.adult.ari, rr.paed.ari, rr.lca, rr.resdis, rr.resp)

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

mort.pooled <- rbind (rr.still, rr.cvdmort, rr.respmort, rr.allmort)

all.pooled <- rbind(resp.pooled, cvd.pooled, rr.lbw, mort.pooled)
all.pooled <- all.pooled %>% mutate_at(vars(pred, ci.lb, ci.ub), funs(round(., 2)))
all.pooled$rr <-paste0(all.pooled$pred," ","(",all.pooled$ci.lb,"-",all.pooled$ci.ub,")") 
all.pooled <- all.pooled %>% 
  select(subgroup, estimates,pred, ci.lb, ci.ub,rr,i2)


##============================================
## plot forestplots
##============================================

dat.ratio <- dat.ratio %>% 
  mutate_at(vars(ratio.ce, ratio.ll, ratio.ul),
  funs(round(.,2)))

dat.ratio$rr <- dat.ratio$rr <- paste0(dat.ratio$ratio.ce," [",dat.ratio$ratio.ll,"-",dat.ratio$ratio.ul,"]")


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


## all deaths ##
rr.alldeat <- rma(yi = lgrr,sei = se, data = dat.ratio.death,method = "ML", slab=dat.ratio.death$author) 
rr.alldeath <- predict(rr.alldeat, transf = exp, digit=2)
rr.alldeath <- as.data.frame(rr.alldeath)
rr.alldeath$subgroup <- "all death"
rr.alldeath$i2 <- round(rr.alldeat$I2, digits=1)
rr.alldeath$estimates <- rr.alldeat$k


par(mar=c(4,4,1,2), font=1)

forest(rr.alldeat, xlim=c(-4, 4), at=log(c(0.5,1,2,3,4)), atransf=exp,
       ilab=cbind(dat.ratio.death$pub.year,dat.ratio.death$rr),
       ilab.xpos=c(-3,-2),cex=0.75, ylim=c(-1, 107),
       order=rev(order(dat.ratio.death$outcome.recode,dat.ratio.death$pub.year)),
       rows=c(3:8,14:19,25:39, 45:82, 88:103),
       xlab="Risk Ratio", mlab="",
       col = "red", border="red")

### add text with Q-value, dfs, p-value, and I^2 statistic
text(-4, -1, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (df = ", .(rr.alldeat$k - rr.alldeat$p),
                                           "; ", I^2, " = ",.(formatC(rr.alldeat$I2, digits=1, format="f")), "%),  ", "Overall risk ratio = ", .(formatC(rr.alldeath$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.alldeath$ci.lb, digits=2, format="f")),
                                           " - ",
                                           .(formatC(rr.alldeath$ci.ub, digits=2, format="f")),
                                           "]")))

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)

### add text for the subgroups
text(-4, c(104,83,40,20,9), pos=4, c("Stillbirth",
                               "Under five mortality",
                               "Cardiovascular mortality",
                               "Respiratory mortality",
                               "All cause mortality"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-4,106, "Author(s)",  pos=4)
text(-3,106,"Year")
text(-2,106, "Risk Ratio [95% CI]")

### set par back to the original settings
par(op)


### add summary polygons for the three subgroups
addpoly(rr.stil, row=86.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")
addpoly(rr.ufmor, row=43.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")
addpoly(rr.cvdmor, row= 23.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")
addpoly(rr.respmor, row=12.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")
addpoly(rr.allmor, row= 1.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")

### add text with Q-value, dfs, p-value, and I^2 statistic for subgroups

text(-4, 86.5, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (df = ", .(rr.stil$k - rr.stil$p),
                                             "; ", I^2, " = ",.(formatC(rr.stil$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.still$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.still$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.still$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 43.5, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (df = ", .(rr.ufmor$k - rr.ufmor$p),
                                             "; ", I^2, " = ",.(formatC(rr.ufmor$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.ufmort$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.ufmort$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.ufmort$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 23.5, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (df = ", .(rr.cvdmor$k - rr.cvdmor$p),
                                             "; ", I^2, " = ",.(formatC(rr.cvdmor$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.cvdmort$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.cvdmort$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.cvdmort$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 12.5, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (df = ", .(rr.respmor$k - rr.respmor$p),
                                             "; ", I^2, " = ",.(formatC(rr.respmor$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.respmort$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.respmort$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.respmort$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 1.5, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (df = ", .(rr.allmor$k - rr.allmor$p),
                                             "; ", I^2, " = ",.(formatC(rr.allmor$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.allmort$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.allmort$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.allmort$ci.ub, digits=2, format="f")),
                                             "]")))


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


par(mar=c(4,4,1,2), font=1)

forest(rr.death.pool, xlim=c(-4, 4), at=log(c(0.5,1,2,3)), atransf=exp,
       ilab=cbind(death.pooled$estimates, paste0(death.pooled$i2, "%"), death.pooled$rr),
       ilab.xpos=c(-2.5,-1.9, -1),cex=1, ylim=c(-2.5, 7), psize=1, col="white",
       order=rev(order(death.pooled$analysis)),
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

all.pooled <- all.pooled %>% 
  filter(estimates>=3)

all.pooled$lgrr <- logb(all.pooled$pred)
all.pooled$lgll <- logb(all.pooled$ci.lb)
all.pooled$lgul <- logb(all.pooled$ci.ub)
all.pooled$se <- (all.pooled$lgul-all.pooled$lgll)/(2*1.96)

rr.pool <- rma(yi = lgrr,sei = se, data = all.pooled ,method = "ML", slab=all.pooled$disease)


pdf("results.v2/Prospective.summary.pdf",
    height=8,
    width=12)

par(mar=c(4,4,1,2), font=1)

forest(rr.pool, xlim=c(-4, 1), at=log(c(0.5,1, 2,3,4,5)), atransf=exp,
       ilab=cbind(all.pooled$estimates, paste0(all.pooled$i2, "%"), all.pooled$rr),
       ilab.xpos=c(-2.5,-1.9, -1),cex=1, ylim=c(0, 19), psize=1, 
       order=rev(order(all.pooled$analysis)),
       rows=c(1:2, 5:7, 10:15),
       addfit=FALSE,
       annotate=FALSE,
       xlab="Risk Ratio", mlab="")

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)

### add text for the subgroups
text(-4, c(16,8,3), pos=4, c("Respiratory diseases",
                             "Cardiovascular diseases",
                             "Adverse pregnancy outcomes"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-4,18, "Disease",  pos=4)
text(-2.5,18,"No. of \n Estimates")
text(-1.9,18, expression("I"^2))
text(-1,18, "Pooled risk ratios \n [95% CI]")

### set par back to the original settings
par(op)

dev.off()
