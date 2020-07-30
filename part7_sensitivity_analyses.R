
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

###====================================
### Part 8 -  estimate pooled rate and risk ratio
###====================================

# load files

dat <- readRDS("data_derived.v2/dat.rds")
dat$rr <- paste0(round(dat$ratio.ce,2)," [",round(dat$ratio.ll,2),"-",round(dat$ratio.ul,2),"]")

dat$events <- paste0(dat$exposed.event, "/",
                                 (dat$exposed.event+dat$exposed.health), " vs ", dat$control.event, "/",
                                 (dat$control.event+dat$control.health))

dat$events <- ifelse(is.na(dat$exposed.event), 
                                 "not reported", dat$events)


#### LPG forest plot ####

dat.lpg <- dat
dat.lpg <- subset(dat.lpg,!is.na(dat.lpg$ratio.ll))

dat.lpg <- dat.lpg %>% 
  select("study.id", "author", "pub.year", "overall.estimate", "analysis.code", "subgroup.id", "subgroup.gender", "subgroup.exposure", "subgroup.outcome", "subgroup.other", "interventional.study", "intervention.type", "exposed.group", "control.group", "outcome.recode", "outcome.type", "risk.measure", "ratio.ce", "ratio.ll", "ratio.ul", "lgrr", "lgll", "lgul", "se", "rr", "events")%>% 
  filter(overall.estimate==1 & analysis.code %in% c(4) & outcome.type!="symptoms")


outcomes.lpg<-as.data.frame(table(dat.lpg$outcome.recode))

# calculate pooled rr

rr.lpg.res <- rma(yi = lgrr,sei = se, data = dat.lpg[dat.lpg$outcome.recode %in% c("asthma", "ari", "lung cancer", "respiratory disease"),],method = "ML", slab=dat.lpg[dat.lpg$outcome.recode %in% c("asthma", "ari", "lung cancer", "respiratory disease"),]$author) 
rr.lpg.resp <- predict(rr.lpg.res, transf = exp, digit=2)
rr.lpg.resp <- as.data.frame(rr.lpg.resp)
rr.lpg.resp$subgroup <- "respiratory overall"
rr.lpg.resp$i2 <- round(rr.lpg.res$I2, digits=1)
rr.lpg.resp$estimates <- rr.lpg.res$k

rr.lpg.asth <- rma(yi = lgrr,sei = se, data = dat.lpg[dat.lpg$outcome.recode %in% c("asthma"),],method = "ML", slab=dat.lpg[dat.lpg$outcome.recode %in% c("asthma"),]$author) 
rr.lpg.asthma <- predict(rr.lpg.asth, transf = exp, digit=2)
rr.lpg.asthma <- as.data.frame(rr.lpg.asthma)
rr.lpg.asthma$subgroup <- "asthma"
rr.lpg.asthma$i2 <- round(rr.lpg.asth$I2, digits=1)
rr.lpg.asthma$estimates <- rr.lpg.asth$k

rr.lpg.ar <- rma(yi = lgrr,sei = se, data = dat.lpg[dat.lpg$outcome.recode %in% c("ari"),],method = "ML", slab=dat.lpg[dat.lpg$outcome.recode %in% c("ari"),]$author) 
rr.lpg.ari <- predict(rr.lpg.ar, transf = exp, digit=2)
rr.lpg.ari <- as.data.frame(rr.lpg.ari)
rr.lpg.ari$subgroup <- "ari"
rr.lpg.ari$i2 <- round(rr.lpg.ar$I2, digits=1)
rr.lpg.ari$estimates <- rr.lpg.ar$k

rr.lpg.lc <- rma(yi = lgrr,sei = se, data = dat.lpg[dat.lpg$outcome.recode %in% c("lung cancer"),],method = "ML", slab=dat.lpg[dat.lpg$outcome.recode %in% c("lung cancer"),]$author) 
rr.lpg.lca <- predict(rr.lpg.lc, transf = exp, digit=2)
rr.lpg.lca <- as.data.frame(rr.lpg.lca)
rr.lpg.lca$subgroup <- "lung cancer"
rr.lpg.lca$i2 <- round(rr.lpg.lc$I2, digits=1)
rr.lpg.lca$estimates <- rr.lpg.lc$k

rr.lpg.resd <- rma(yi = lgrr,sei = se, data = dat.lpg[dat.lpg$outcome.recode %in% c("respiratory disease"),],method = "ML", slab=dat.lpg[dat.lpg$outcome.recode %in% c("respiratory disease"),]$author) 
rr.lpg.resdis <- predict(rr.lpg.resd, transf = exp, digit=2)
rr.lpg.resdis <- as.data.frame(rr.lpg.resdis)
rr.lpg.resdis$subgroup <- "respiratory disease"
rr.lpg.resdis$i2 <- round(rr.lpg.resd$I2, digits=1)
rr.lpg.resdis$estimates <- rr.lpg.resd$k

lpg.pooled <- rbind (rr.lpg.asthma, rr.lpg.ari, rr.lpg.lca, rr.lpg.resdis, rr.lpg.resp)


# forest plot
dat.lpg <- dat.lpg[dat.lpg$outcome.recode %in% c("asthma", "ari", "lung cancer", "respiratory disease"),]

dat.lpg <- dat.lpg %>% 
  mutate(analysis = case_when(
    outcome.recode %in% c("asthma") ~ "a.asthma",
    outcome.recode %in% c("ari") ~ "b.ari",
    outcome.recode %in% c("lung cancer") ~ "c.lung cancer",
    outcome.recode %in% c("respiratory disease") ~ "d.respiratory disease",
    TRUE ~ as.character(outcome.recode)
  )) 


pdf(file="results.v2/LPG forest plot.pdf", height=13, width=11)

par(mar=c(4,4,1,2), font=1)

forest(rr.lpg.res, xlim=c(-4, 4.5), 
       at=log(c(0.25,0.5,0.75,1,2,3,4,5,6)), atransf=exp,
       ilab=cbind(dat.lpg$pub.year,dat.lpg$events,dat.lpg$rr),
       ilab.xpos=c(-2,2.5, 3.8),cex=0.75, ylim=c(1, 90),
       order=rev(order(dat.lpg$analysis,dat.lpg$pub.year)),
       rows=c(2:6,12:17,23:39, 45:86),
       annotate=FALSE,
       addfit=FALSE,
       xlab="Risk Ratio", mlab="",
       col = "red", border="red")

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.8, font=4)

### add text for the subgroups
text(-4, c(7,18,40, 87), pos=4, c("Respiratory disease", "Lung cancer", "Acute respiratory infection", "Asthma"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-4,89, "Author(s)",  pos=4)
text(-2,89,"Year")
text(2.5,89.2, "Number of events \n in exposed vs unexposed")
text(3.8,89, "Risk Ratio [95% CI]")

### set par back to the original settings
par(op)


### add summary polygons for the three subgroups
addpoly(rr.lpg.asth, row=43.5, cex=0.75, atransf=exp, mlab="",col = "red", border="red", annotate=FALSE)
addpoly(rr.lpg.ar, row=21.5, cex=0.75, atransf=exp, mlab="",col = "red", border="red", annotate=FALSE)
addpoly(rr.lpg.lc, row= 10.5, cex=0.75, atransf=exp, mlab="",col = "red", border="red", annotate=FALSE)
addpoly(rr.lpg.resd, row= 0.5, cex=0.75, atransf=exp, mlab="",col = "red", border="red", annotate=FALSE)


### add text with Q-value, dfs, p-value, and I^2 statistic for subgroups

text(-4, 43.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.lpg.asth$k - rr.lpg.asth$p),
                                             "; ", I^2, " = ",.(formatC(rr.lpg.asth$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.lpg.asthma$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.lpg.asthma$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.lpg.asthma$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 21.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.lpg.ar$k - rr.lpg.ar$p),
                                             "; ", I^2, " = ",.(formatC(rr.lpg.ar$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.lpg.ari$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.lpg.ari$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.lpg.ari$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 10.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.lpg.lc$k - rr.lpg.lc$p),
                                             "; ", I^2, " = ",.(formatC(rr.lpg.lc$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.lpg.lca$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.lpg.lca$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.lpg.lca$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 0.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.lpg.resd$k - rr.lpg.resd$p),
                                            "; ", I^2, " = ",.(formatC(rr.lpg.resd$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.lpg.resdis$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.lpg.resdis$ci.lb, digits=2, format="f")),
                                            " - ",
                                            .(formatC(rr.lpg.resdis$ci.ub, digits=2, format="f")),
                                            "]")))

dev.off()


##### diagnosis and mortality for clean fuel control only ####

dat.ratio <- dat
dat.ratio <- subset(dat.ratio,!is.na(dat.ratio$ratio.ll))

dat.ratio <- dat.ratio %>% 
  select("study.id", "author", "pub.year", "overall.estimate", "analysis.code", "subgroup.id", "subgroup.gender", "subgroup.exposure", "subgroup.outcome", "subgroup.other", "interventional.study", "intervention.type", "exposed.group", "control.group", "outcome.recode", "outcome.type", "adult.paediatric.all","paed.under5", "risk.measure", "ratio.ce", "ratio.ll", "ratio.ul", "lgrr", "se", "events")

# diagnosis and mortality only
dat.ratio <- dat.ratio %>% 
  filter(overall.estimate==1 & analysis.code %in% c(1,5) & outcome.type!="symptoms")

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

mort.pooled <- rbind (rr.still, rr.ufmort)

all.pooled <- rbind(resp.pooled, cvd.pooled, rr.lbw, mort.pooled)
all.pooled <- all.pooled %>% mutate_at(vars(pred, ci.lb, ci.ub), funs(round(., 2)))
all.pooled$rr <-paste0(all.pooled$pred," ","(",all.pooled$ci.lb,"-",all.pooled$ci.ub,")") 
all.pooled <- all.pooled %>% 
  select(subgroup, estimates,pred, ci.lb, ci.ub,rr,i2)

write.csv(all.pooled, "results.v2/all.clean.pooled.csv")

#### place beta and se into dataframe for burden analysis ####
endpoints <- c("asthma", "copd", "ari", "lungca", "tb", "ihd", "stroke")

rr.sim.para <- data.frame(endpoints = endpoints, 
                          rr.beta = c(rr.asth$b, rr.cop$b, rr.ar$b, rr.lc$b, rr.t$b, rr.ih$b, rr.strok$b),
                          rr.se = c(rr.asth$se, rr.cop$se, rr.ar$se, rr.lc$se, rr.t$se, rr.ih$se, rr.strok$se))

write.csv(rr.sim.para, "data_derived.v2/rr.sim.para.clean.csv")


## summary clean plot ##

all.pooled <- all.pooled %>% 
  mutate(analysis = case_when(
    subgroup %in% c("asthma", "copd", "ari- adult", "ari- paed", "lung cancer", "tb", "respiratory disease") ~ "a.resp",
    subgroup %in% c("stroke", "ihd") ~ "b.cvd",
    subgroup %in% c("lbw") ~ "c.lbw",
    subgroup %in% c("stillbirth") ~ "d.stillbirth",
    subgroup %in% c("under five death") ~ "e.under five death",
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
    subgroup %in% c("lbw")~ "Low birth weight",
    subgroup %in% c("stillbirth") ~ "Stillbirth",
    subgroup %in% c("under five death") ~ "Under-five death",
    TRUE ~ as.character(subgroup)
  ))


all.pooled$lgrr <- logb(all.pooled$pred)
all.pooled$lgll <- logb(all.pooled$ci.lb)
all.pooled$lgul <- logb(all.pooled$ci.ub)
all.pooled$se <- (all.pooled$lgul-all.pooled$lgll)/(2*1.96)

rr.pool <- rma(yi = lgrr,sei = se, data = all.pooled ,method = "ML", slab=all.pooled$disease)


pdf(file="results.v2/All pooled. clean comparator.pdf", height=8, width=14)

par(mar=c(4,4,1,2), font=1)

forest(rr.pool, xlim=c(-4, 1), at=log(c(0.5,1,2,3,4,5)), atransf=exp,
       ilab=cbind(all.pooled$estimates, paste0(all.pooled$i2, "%"), all.pooled$rr),
       ilab.xpos=c(-2.5,-1.9, -1),cex=1, ylim=c(-1, 20), psize=1, 
       order=rev(order(all.pooled$analysis)),
       annotate=FALSE,
       addfit=FALSE,
       rows=c(1:3, 6:7, 10:16),
       xlab="Risk Ratio", mlab="")

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)

### add text for the subgroups
text(-4, c(17,8,4), pos=4, c("Respiratory diseases",
                             "Cardiovascular diseases",
                             "Adverse pregnancy and pediatric outcomes"))

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


#### intervention forest plot ####

dat.int <- dat
dat.int <- subset(dat.int,!is.na(dat.int$ratio.ll))

dat.int <- dat.int %>% 
  select("study.id", "author", "pub.year", "overall.estimate", "analysis.code", "subgroup.id", "subgroup.gender", "subgroup.exposure", "subgroup.outcome", "subgroup.other", "interventional.study", "intervention.type", "exposed.group", "control.group", "outcome.recode", "outcome.type", "risk.measure", "ratio.ce", "ratio.ll", "ratio.ul", "lgrr", "lgll", "lgul", "se", "rr", "events")%>% 
  filter(overall.estimate %in% c(1,3) & analysis.code %in% c(9))

outcomes.int<-as.data.frame(table(dat.int$outcome.recode))


dat.int <- dat.int %>% 
  mutate(analysis = case_when(
    outcome.recode %in% c("ari") ~ "a.ari",
    outcome.recode %in% c("asthma") ~ "b.asthma",
    outcome.recode %in% c("lung cancer") ~ "d.lung cancer",
    outcome.recode %in% c("wheeze") ~ "e.wheeze",
    outcome.recode %in% c("cough") ~ "f.cough", 
    outcome.recode %in% c("sob") ~ "g.sob", 
    TRUE ~ "exclude"
  )) %>% 
  filter(analysis!="exclude")



rr.int.res <- rma(yi = lgrr,sei = se, data = dat.int[dat.int$outcome.recode %in% c("asthma", "ari", "lung cancer", "tuberculosis", "respiratory disease"),],method = "ML", slab=dat.int[dat.int$outcome.recode %in% c("asthma", "ari", "lung cancer", "tuberculosis", "respiratory disease"),]$author) 
rr.int.resp <- predict(rr.int.res, transf = exp, digit=2)
rr.int.resp <- as.data.frame(rr.int.resp)
rr.int.resp$subgroup <- "respiratory overall"
rr.int.resp$i2 <- round(rr.int.res$I2, digits=1)
rr.int.resp$estimates <- rr.int.res$k

rr.int.asth <- rma(yi = lgrr,sei = se, data = dat.int[dat.int$outcome.recode %in% c("asthma"),],method = "ML", slab=dat.int[dat.int$outcome.recode %in% c("asthma"),]$author) 
rr.int.asthma <- predict(rr.int.asth, transf = exp, digit=2)
rr.int.asthma <- as.data.frame(rr.int.asthma)
rr.int.asthma$subgroup <- "asthma"
rr.int.asthma$i2 <- round(rr.int.asth$I2, digits=1)
rr.int.asthma$estimates <- rr.int.asth$k

rr.int.ar <- rma(yi = lgrr,sei = se, data = dat.int[dat.int$outcome.recode %in% c("ari"),],method = "ML", slab=dat.int[dat.int$outcome.recode %in% c("ari"),]$author) 
rr.int.ari <- predict(rr.int.ar, transf = exp, digit=2)
rr.int.ari <- as.data.frame(rr.int.ari)
rr.int.ari$subgroup <- "ari"
rr.int.ari$i2 <- round(rr.int.ar$I2, digits=1)
rr.int.ari$estimates <- rr.int.ar$k

rr.int.sym <- rma(yi = lgrr,sei = se, data = dat.int[dat.int$outcome.recode %in% c("wheeze", "cough", "sob"),],method = "ML", slab=dat.int[dat.int$outcome.recode %in% c("wheeze", "cough", "sob"),]$author) 
rr.int.symp <- predict(rr.int.sym, transf = exp, digit=2)
rr.int.symp <- as.data.frame(rr.int.symp)
rr.int.symp$subgroup <- "respiratory symptoms"
rr.int.symp$i2 <- round(rr.int.sym$I2, digits=1)
rr.int.symp$estimates <- rr.int.sym$k

rr.int.lc <- rma(yi = lgrr,sei = se, data = dat.int[dat.int$outcome.recode %in% c("lung cancer"),],method = "ML", slab=dat.int[dat.int$outcome.recode %in% c("lung cancer"),]$author) 
rr.int.lca <- predict(rr.int.lc, transf = exp, digit=2)
rr.int.lca <- as.data.frame(rr.int.lca)
rr.int.lca$subgroup <- "lung cancer"
rr.int.lca$i2 <- round(rr.int.lc$I2, digits=1)
rr.int.lca$estimates <- rr.int.lc$k

resp.int.pooled <- rbind (rr.int.asthma, rr.int.ari, rr.int.symp, rr.int.lca, rr.int.resp)


#### Plot for interventional studies ####

rr.int <- rma(yi = lgrr,sei = se, data = dat.int,method = "ML", slab=dat.int$author)
rr.inter <- predict(rr.int, transf = exp, digit=2)
rr.inter <- as.data.frame(rr.inter)
rr.inter$subgroup <- "overall"
rr.inter$i2 <- round(rr.int$I2, digits=1)
rr.inter$estimates <- rr.int$k

#######

pdf(file="results.v2/Interventional studies.pdf", height=10, width=16)

par(mar=c(4,4,1,2), font=1)

forest(rr.int, xlim=c(-4, 3), at=log(c(0.1,0.2,0.4, 0.6, 0.8, 1,2,3,4)), atransf=exp,
       ilab=cbind(dat.int$pub.year,dat.int$events,dat.int$rr),
       ilab.xpos=c(-2.8,1.9, 2.7),cex=0.75, ylim=c(-1, 73),
       order=rev(order(dat.int$outcome.recode,dat.int$pub.year)),
       rows=c(3:34, 39:41, 46:55, 60:69),
       annotate=FALSE,
       addfit=FALSE,
       xlab="Risk Ratio", mlab="",
       col = "red", border="red")

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)

### add text for the subgroups
text(-4, c(35, 42, 56, 70), pos=4, c("Respiratory symptoms", "Lung cancer", "Asthma", "Acute respiratory infection"))

### switch to bold font
par(font=2)


### add column headings to the plot
text(-4,72, "Author(s)",  pos=4)
text(-2.8,72,"Year")
text(1.9,72.2, "Number of events \n in exposed vs unexposed")
text(2.7,72, "Risk Ratio [95% CI]")


### set par back to the original settings
par(op)


### add summary polygons for the three subgroups
addpoly(rr.int.ar, row=58.5, cex=0.5, atransf=exp, mlab="",col = "blue", border="blue", annotate = FALSE)
addpoly(rr.int.asth, row=44.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate = FALSE)
addpoly(rr.int.lc, row= 37.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate = FALSE)
addpoly(rr.int.sym, row= 1.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate = FALSE)

### add text with Q-value, dfs, p-value, and I^2 statistic for subgroups

text(-4, 58.5, pos=4, cex=0.70, bquote(paste("RE Model for acute respiratory infection (df = ", .(rr.int.ar$k - rr.int.ar$p),
                                             "; ", I^2, " = ",.(formatC(rr.int.ar$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.int.ari$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.int.ari$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.int.ari$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 44.5, pos=4, cex=0.70, bquote(paste("RE Model for asthma (df = ", .(rr.int.asth$k - rr.int.asth$p),
                                             "; ", I^2, " = ",.(formatC(rr.int.asth$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.int.asthma$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.int.asthma$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.int.asthma$ci.ub, digits=2, format="f")),
                                             "]")))


text(-4, 37.5, pos=4, cex=0.70, bquote(paste("RE Model for lung cancer (df = ", .(rr.int.lc$k - rr.int.lc$p),
                                             "; ", I^2, " = ",.(formatC(rr.int.lc$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.int.lca$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.int.lca$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.int.lca$ci.ub, digits=2, format="f")),
                                             "]")))

text(-4, 1.5, pos=4, cex=0.70, bquote(paste("RE Model for respiratory symptoms (df = ", .(rr.int.sym$k - rr.int.sym$p),
                                             "; ", I^2, " = ",.(formatC(rr.int.sym$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.int.symp$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.int.symp$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.int.symp$ci.ub, digits=2, format="f")),
                                             "]")))

dev.off()


#### intervention forest plot-RCT ####

dat.baseline.design <- read.csv("data_raw/iap.baseline.v8.csv") %>% 
  select(study.id, study.design)
  
dat.rct <- left_join(dat, dat.baseline.design) %>% 
  filter(study.design=="rct")

dat.rct <- subset(dat.rct,!is.na(dat.rct$ratio.ll))

dat.rct <- dat.rct %>% 
  select("study.id", "author", "pub.year", "overall.estimate", "analysis.code", "subgroup.id", "subgroup.gender", "subgroup.exposure", "subgroup.outcome", "subgroup.other", "interventional.study", "intervention.type", "exposed.group", "control.group", "outcome.recode", "outcome.type", "risk.measure", "ratio.ce", "ratio.ll", "ratio.ul", "lgrr", "lgll", "lgul", "se", "rr", "events")%>% 
  filter(overall.estimate %in% c(1,3) & analysis.code %in% c(9))

outcomes.rct<-as.data.frame(table(dat.rct$outcome.recode))


dat.rct <- dat.rct %>% 
  mutate(analysis = case_when(
    outcome.recode %in% c("ari") ~ "Acute respiratory infection",
    outcome.recode %in% c("asthma") ~ "Asthma",
    outcome.recode %in% c("wheeze") ~ "Wheeze",
    outcome.recode %in% c("cough") ~ "Cough", 
    outcome.recode %in% c("sob") ~ "Shortness of breath", 
    TRUE ~ "exclude"
  )) %>% 
  filter(analysis!="exclude")

rr.rc <- rma(yi = lgrr,sei = se, data = dat.rct ,method = "ML", slab=dat.rct$author) 
rr.rct <- predict(rr.rc, transf = exp, digit=2)
rr.rct <- as.data.frame(rr.rct)
rr.rct$subgroup <- "rct- overall"
rr.rct$i2 <- round(rr.rc$I2, digits=1)
rr.rct$estimates <- rr.rc$k


pdf("results.v2/rct.pdf", height=6, width=12)

par(mar=c(4,4,1,4), font=1)

forest(rr.rc, xlim=c(-6, 4), at=log(c(0.2, 0.4, 0.6, 0.8,1,2,3,4)), atransf=exp,
       ilab=cbind(dat.rct$pub.year,dat.rct$rr),
       ilab.xpos=c(-4,-2.5),cex=0.75, ylim=c(-1, 17),
       order=rev(order(dat.rct$analysis, dat.rct$pub.year)),
       rows=c(1:6, 9, 12:13),
       xlab="Risk Ratio", mlab="",
       col = "blue", border="blue")

### add text with Q-value, dfs, p-value, and I^2 statistic
text(-6, -1, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (df = ", .(rr.rc$k - rr.rc$p),
                                           "; ", I^2, " = ",.(formatC(rr.rc$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.rct$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.rct$ci.lb, digits=2, format="f")),
                                           " - ",
                                           .(formatC(rr.rct$ci.ub, digits=2, format="f")),
                                           "]")))

### set font expansion factor (as in forest() above) and use bold italic

### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)


### add text for the subgroups
text(-6, c(7, 10, 14), pos=4, c("Respiratory symptoms", "Asthma", "Acute respiratory infection"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-6,16, "Author(s)",  pos=4)
text(-4,16,"Year")
text(-2.5,16, "Risk Ratio [95% CI]")

### set par back to the original settings
par(op)


dev.off()



#### symptoms forest plot ####

dat.symp <- dat
dat.symp <- subset(dat.symp,!is.na(dat.symp$ratio.ll))

dat.symp <- dat.symp %>% 
  select("study.id", "author", "pub.year", "overall.estimate", "analysis.code", "subgroup.id", "subgroup.gender", "subgroup.exposure", "subgroup.outcome", "subgroup.other", "interventional.study", "intervention.type", "exposed.group", "control.group", "outcome.recode", "outcome.type", "adult.paediatric.all", "risk.measure", "ratio.ce", "ratio.ll", "ratio.ul", "lgrr", "lgll", "lgul", "se", "rr", "events")%>% 
  filter(overall.estimate %in% c(1, 3) & analysis.code %in% c(1,2,5) & outcome.type=="symptoms")

outcomes.symp<-as.data.frame(table(dat.symp$outcome.recode))


rr.wheez <- rma(yi = lgrr,sei = se, data = dat.symp[dat.symp$outcome.recode %in% c("wheeze"),],method = "ML", slab=dat.symp[dat.symp$outcome.recode %in% c("wheeze"),]$author) 
rr.wheeze <- predict(rr.wheez, transf = exp, digit=2)
rr.wheeze <- as.data.frame(rr.wheeze)
rr.wheeze$subgroup <- "wheeze"
rr.wheeze$i2 <- round(rr.wheez$I2, digits=1)
rr.wheeze$estimates <- rr.wheez$k


rr.coug <- rma(yi = lgrr,sei = se, data = dat.symp[dat.symp$outcome.recode %in% c("cough"),],method = "ML", slab=dat.symp[dat.symp$outcome.recode %in% c("cough"),]$author) 
rr.cough <- predict(rr.coug, transf = exp, digit=2)
rr.cough <- as.data.frame(rr.cough)
rr.cough$subgroup <- "cough"
rr.cough$i2 <- round(rr.coug$I2, digits=1)
rr.cough$estimates <- rr.coug$k


rr.so <- rma(yi = lgrr,sei = se, data = dat.symp[dat.symp$outcome.recode %in% c("sob"),],method = "ML", slab=dat.symp[dat.symp$outcome.recode %in% c("sob"),]$author) 
rr.sob <- predict(rr.so, transf = exp, digit=2)
rr.sob <- as.data.frame(rr.sob)
rr.sob$subgroup <- "sob"
rr.sob$i2 <- round(rr.so$I2, digits=1)
rr.sob$estimates <- rr.so$k


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
forest_plot(data=dat.symp, 
            endpoint="wheeze", 
            age_group=c("adult", "both", "paed"), 
            subtitle= "all")

forest_plot(data=dat.symp, 
            endpoint="cough", 
            age_group=c("adult", "both", "paed"), 
            subtitle= "all")

forest_plot(data=dat.symp, 
            endpoint="sob", 
            age_group=c("adult", "both", "paed"), 
            subtitle= "all")


#### birth weight forest plot ####

dat.weight <- dat
dat.weight <- subset(dat.weight,!is.na(dat.weight$bw.exp))

dat.weight <- escalc(measure="MD",  m1i=bw.exp, sd1i=sd.exp, n1i=n.exp, m2i=bw.con, sd2i=sd.con, n2i=n.exp, data=dat.weight)

rr.weigh <- rma(yi=yi, vi=vi, data=dat.weight, slab=dat.weight$author)

rr.weight <- predict(rr.weigh, digit=2)
rr.weight <- as.data.frame(rr.weight)
rr.weight$subgroup <- "weight"
rr.weight$i2 <- round(rr.weigh$I2, digits=1)
rr.weight$estimates <- rr.weigh$k


## Plot for weight ##
### decrease margins so the full space is used

pdf("results.v2/birthweight.pdf",
    height=8,
    width=12)

par(mar=c(4,4,1,4), font=1)

forest(rr.weigh, xlim=c(-1200, 550), 
       slab=dat.weight$author,
       ilab.xpos=c(-850), ilab=cbind(dat.weight$pub.year), 
       ylim=c(-2, 26),
       rows=c(1:23),
       order=rev(order(dat.weight$pub.year)),
       xlab="Mean difference in birth weight (grams)", mlab="",
       col = "blue", border="blue")

### add text with Q-value, dfs, p-value, and I^2 statistic
text(-1200, -1, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (df = ", .(rr.weigh$k - rr.weigh$p),
                                           "; ", I^2, " = ",.(formatC(rr.weigh$I2, digits=1, format="f")), "%)")))  

### set font expansion factor (as in forest() above) and use bold italic

### font and save original settings in object 'op'
op <- par(cex=0.90, font=4)


### switch to bold font
par(font=2)

### add column headings to the plot
text(-1200,25, "Author(s)",  pos=4)
text(-850,25,"Year")
text(380,25, "Mean difference in \n birth weight, grams [95% CI]")

### set par back to the original settings
par(op)

dev.off()

##############

#### Direct pollutant forest plot ####
dat.direct <- dat %>% filter(analysis.code==88 & overall.estimate %in% c(1,3))

dat.direct[,"pollutant.increment"] <- as.numeric(as.character(dat.direct[,"pollutant.increment"]))

dat.direct <- dat.direct %>% filter(pollutant %in% c("PM2.5", "NO2") & outcome.recode %in% c("asthma","wheeze", "sob", "cough"))

outcomes.direct<-as.data.frame(table(dat.direct$pollutant, dat.direct$outcome.recode)) %>% filter(Freq!=0)

dat.direct <- dat.direct %>% 
  mutate(analysis = case_when(
    outcome.recode %in% c("asthma") ~ "a.asthma",
    outcome.recode %in% c("wheeze") ~ "b.wheeze",
    outcome.recode %in% c("cough") ~ "c.cough", 
    outcome.recode %in% c("sob") ~ "d.sob", 
    TRUE ~ as.character(outcome.recode)
  ))



# NO2 #
dat.direct.no2 <- dat.direct %>% filter(pollutant=="NO2")

dat.direct.no2 <- dat.direct.no2 %>% 
  mutate_at(vars(pollutant.increment),
            funs(case_when(pollutant.units!="ppb" ~ ./1.88,
                 TRUE~.
                 )))

dat.direct.no2 <- dat.direct.no2 %>% 
  mutate_at(vars(ratio.ce, ratio.ll, ratio.ul),
            funs(.^(10/pollutant.increment)))

dat.direct.no2$lgrr <- logb(dat.direct.no2$ratio.ce)
dat.direct.no2$lgll <- logb(dat.direct.no2$ratio.ll)
dat.direct.no2$lgul <- logb(dat.direct.no2$ratio.ul)
dat.direct.no2$se <- (dat.direct.no2$lgul-dat.direct.no2$lgll)/(2*1.96)

num_iter <- 10000

dat.direct.no2$rr <- paste0(round(dat.direct.no2$ratio.ce,2)," [",round(dat.direct.no2$ratio.ll,2),"-",round(dat.direct.no2$ratio.ul,2),"]")

rr.no <- rma(yi = lgrr,sei = se, data = dat.direct.no2 ,method = "ML", slab=dat.direct.no2$author) 
rr.no2 <- predict(rr.no, transf = exp, digit=2)
rr.no2 <- as.data.frame(rr.no2)
rr.no2$subgroup <- "no2"
rr.no2$i2 <- round(rr.no$I2, digits=1)
rr.no2$estimates <- rr.no$k

rr.no2.asth <- rma(yi = lgrr,sei = se, data = dat.direct.no2[dat.direct.no2$outcome.recode %in% c("asthma"),] ,method = "ML", slab=dat.direct.no2[dat.direct.no2$outcome.recode %in% c("asthma"),]$author) 
rr.no2.asthma <- predict(rr.no2.asth, transf = exp, digit=2)
rr.no2.asthma <- as.data.frame(rr.no2.asthma)
rr.no2.asthma$subgroup <- "no2.asthma"
rr.no2.asthma$i2 <- round(rr.no2.asth$I2, digits=1)
rr.no2.asthma$estimates <- rr.no2.asth$k


rr.no2.wheez <- rma(yi = lgrr,sei = se, data = dat.direct.no2[dat.direct.no2$outcome.recode %in% c("wheeze"),] ,method = "ML", slab=dat.direct.no2[dat.direct.no2$outcome.recode %in% c("wheeze"),]$author) 
rr.no2.wheeze <- predict(rr.no2.wheez, transf = exp, digit=2)
rr.no2.wheeze <- as.data.frame(rr.no2.wheeze)
rr.no2.wheeze$subgroup <- "no2.wheeze"
rr.no2.wheeze$i2 <- round(rr.no2.wheez$I2, digits=1)
rr.no2.wheeze$estimates <- rr.no2.wheez$k


rr.no2.coug <- rma(yi = lgrr,sei = se, data = dat.direct.no2[dat.direct.no2$outcome.recode %in% c("cough"),] ,method = "ML", slab=dat.direct.no2[dat.direct.no2$outcome.recode %in% c("cough"),]$author) 
rr.no2.cough <- predict(rr.no2.coug, transf = exp, digit=2)
rr.no2.cough <- as.data.frame(rr.no2.cough)
rr.no2.cough$subgroup <- "no2.cough"
rr.no2.cough$i2 <- round(rr.no2.coug$I2, digits=1)
rr.no2.cough$estimates <- rr.no2.coug$k


rr.no2.so <- rma(yi = lgrr,sei = se, data = dat.direct.no2[dat.direct.no2$outcome.recode %in% c("sob"),] ,method = "ML", slab=dat.direct.no2[dat.direct.no2$outcome.recode %in% c("sob"),]$author) 
rr.no2.sob <- predict(rr.no2.so, transf = exp, digit=2)
rr.no2.sob <- as.data.frame(rr.no2.sob)
rr.no2.sob$subgroup <- "no2.sob"
rr.no2.sob$i2 <- round(rr.no2.so$I2, digits=1)
rr.no2.sob$estimates <- rr.no2.so$k

rr.pooled.no2 <- rbind(rr.no2, rr.no2.asthma, rr.no2.wheeze, rr.no2.cough, rr.no2.sob)


pdf(file="results.v2/NO2.pdf", height=10, width=12)

par(mar=c(4,4,1,2), font=1)

forest(rr.no, xlim=c(-3.5, 1), at=log(c(0.25,0.5,1,2,3,4,5,6)), atransf=exp,
       ilab=cbind(dat.direct.no2$pub.year,dat.direct.no2$rr),
       ilab.xpos=c(-2.5,-1.5),cex=0.75, ylim=c(0, 62),
       order=rev(order(dat.direct.no2$analysis,dat.direct.no2$pub.year)),
       rows=c(3:10,16:28,34:48, 54:58),
       annotate=FALSE,
       addfit=FALSE,
       xlab="Risk Ratio", mlab="",
       col = "red", border="red")

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.8, font=4)

### add text for the subgroups
text(-3.5, c(11,29,49, 59), pos=4, c("Shortness of breath", "Cough", "Wheeze", "Asthma"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-3.5,62, "Author(s)",  pos=4)
text(-2.5,62,"Year")
text(-1.5,62, "Risk Ratio [95% CI]\n(per 10 ppb increment)")

### set par back to the original settings
par(op)


### add summary polygons for the three subgroups
addpoly(rr.no2.asth, row=52.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate=FALSE)
addpoly(rr.no2.wheez, row=32.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate=FALSE)
addpoly(rr.no2.coug, row= 14.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate=FALSE)
addpoly(rr.no2.so, row= 1.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate=FALSE)


### add text with Q-value, dfs, p-value, and I^2 statistic for subgroups

text(-3.5, 52.5, pos=4, cex=0.70, bquote(paste("RE Model for asthma (df = ", .(rr.no2.asth$k - rr.no2.asth$p),
                                             "; ", I^2, " = ",.(formatC(rr.no2.asth$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.no2.asthma$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.no2.asthma$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.no2.asthma$ci.ub, digits=2, format="f")),
                                             "]")))


text(-3.5, 32.5, pos=4, cex=0.70, bquote(paste("RE Model for wheeze (df = ", .(rr.no2.wheez$k - rr.no2.wheez$p),
                                             "; ", I^2, " = ",.(formatC(rr.no2.wheez$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.no2.wheeze$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.no2.wheeze$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.no2.wheeze$ci.ub, digits=2, format="f")),
                                             "]")))


text(-3.5, 14.5, pos=4, cex=0.70, bquote(paste("RE Model for cough (df = ", .(rr.no2.coug$k - rr.no2.coug$p),
                                             "; ", I^2, " = ",.(formatC(rr.no2.coug$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.no2.cough$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.no2.cough$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.no2.cough$ci.ub, digits=2, format="f")),
                                             "]")))


text(-3.5, 1.5, pos=4, cex=0.70, bquote(paste("RE Model for shortness of breath (df = ", .(rr.no2.so$k - rr.no2.so$p),
                                            "; ", I^2, " = ",.(formatC(rr.no2.so$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.no2.sob$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.no2.sob$ci.lb, digits=2, format="f")),
                                            " - ",
                                            .(formatC(rr.no2.sob$ci.ub, digits=2, format="f")),
                                            "]")))

dev.off()


# PM2.5 #
dat.direct.pm2.5 <- dat.direct %>% filter(pollutant=="PM2.5")

dat.direct.pm2.5 <- dat.direct.pm2.5 %>% 
  mutate_at(vars(ratio.ce, ratio.ll, ratio.ul),
            funs(.^(10/pollutant.increment)))


dat.direct.pm2.5$lgrr <- logb(dat.direct.pm2.5$ratio.ce)
dat.direct.pm2.5$lgll <- logb(dat.direct.pm2.5$ratio.ll)
dat.direct.pm2.5$lgul <- logb(dat.direct.pm2.5$ratio.ul)
dat.direct.pm2.5$se <- (dat.direct.pm2.5$lgul-dat.direct.pm2.5$lgll)/(2*1.96)

num_iter <- 10000

dat.direct.pm2.5$rr <- paste0(round(dat.direct.pm2.5$ratio.ce,2)," [",round(dat.direct.pm2.5$ratio.ll,2),"-",round(dat.direct.pm2.5$ratio.ul,2),"]")

rr.pm <- rma(yi = lgrr,sei = se, data = dat.direct.pm2.5 ,method = "ML", slab=dat.direct.pm2.5$author) 
rr.pm2.5 <- predict(rr.pm, transf = exp, digit=2)
rr.pm2.5 <- as.data.frame(rr.pm2.5)
rr.pm2.5$subgroup <- "pm2.5"
rr.pm2.5$i2 <- round(rr.pm$I2, digits=1)
rr.pm2.5$estimates <- rr.pm$k

rr.pm2.5.asth <- rma(yi = lgrr,sei = se, data = dat.direct.pm2.5[dat.direct.pm2.5$outcome.recode %in% c("asthma", "wheeze"),] ,method = "ML", slab=dat.direct.pm2.5[dat.direct.pm2.5$outcome.recode %in% c("asthma", "wheeze"),]$author) 
rr.pm2.5.asthma <- predict(rr.pm2.5.asth, transf = exp, digit=2)
rr.pm2.5.asthma <- as.data.frame(rr.pm2.5.asthma)
rr.pm2.5.asthma$subgroup <- "pm2.5.asthma.wheeze"
rr.pm2.5.asthma$i2 <- round(rr.pm2.5.asth$I2, digits=1)
rr.pm2.5.asthma$estimates <- rr.pm2.5.asth$k


rr.pm2.5.coug <- rma(yi = lgrr,sei = se, data = dat.direct.pm2.5[dat.direct.pm2.5$outcome.recode %in% c("cough"),] ,method = "ML", slab=dat.direct.pm2.5[dat.direct.pm2.5$outcome.recode %in% c("cough"),]$author) 
rr.pm2.5.cough <- predict(rr.pm2.5.coug, transf = exp, digit=2)
rr.pm2.5.cough <- as.data.frame(rr.pm2.5.cough)
rr.pm2.5.cough$subgroup <- "pm2.5.cough"
rr.pm2.5.cough$i2 <- round(rr.pm2.5.coug$I2, digits=1)
rr.pm2.5.cough$estimates <- rr.pm2.5.coug$k


rr.pm2.5.so <- rma(yi = lgrr,sei = se, data = dat.direct.pm2.5[dat.direct.pm2.5$outcome.recode %in% c("sob"),] ,method = "ML", slab=dat.direct.pm2.5[dat.direct.pm2.5$outcome.recode %in% c("sob"),]$author) 
rr.pm2.5.sob <- predict(rr.pm2.5.so, transf = exp, digit=2)
rr.pm2.5.sob <- as.data.frame(rr.pm2.5.sob)
rr.pm2.5.sob$subgroup <- "pm2.5.sob"
rr.pm2.5.sob$i2 <- round(rr.pm2.5.so$I2, digits=1)
rr.pm2.5.sob$estimates <- rr.pm2.5.so$k

rr.pooled.pm2.5 <- rbind(rr.pm2.5, rr.pm2.5.asthma, rr.pm2.5.cough, rr.pm2.5.sob)


pdf("results.v2/PM2.5.pdf",
    height=8,
    width=12)

par(mar=c(4,4,1,2), font=1)

forest(rr.pm, xlim=c(-3.5, 1), at=log(c(0.25, 0.5,1,2,3,4,5,6)), atransf=exp,
       ilab=cbind(dat.direct.pm2.5$pub.year,dat.direct.pm2.5$rr),
       ilab.xpos=c(-2.5,-1.5),cex=0.75, ylim=c(0, 33),
       order=rev(order(dat.direct.pm2.5$analysis,dat.direct.pm2.5$pub.year)),
       rows=c(3:6,12:18,24:29),
       annotate=FALSE,
       addfit=FALSE,
       xlab="Risk Ratio", mlab="",
       col = "red", border="red")

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.8, font=4)

### add text for the subgroups
text(-3.5, c(7,19,30), pos=4, c("Shortness of breath", "Cough", "Asthma or Wheeze"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-3.5,32, "Author(s)",  pos=4)
text(-2.5,32,"Year")
text(-1.5, 32.5, "Risk Ratio [95% CI]")
text(-1.5, 31.8, expression(paste("(per 10 ", mu, "g/m"^"3"," increment)")))

### set par back to the original settings
par(op)


### add summary polygons for the three subgroups
addpoly(rr.pm2.5.asth, row=22.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate=FALSE)
addpoly(rr.pm2.5.coug, row= 10.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate=FALSE)
addpoly(rr.pm2.5.so, row= 1.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue", annotate=FALSE)


### add text with Q-value, dfs, p-value, and I^2 statistic for subgroups

text(-3.5, 22.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.pm2.5.asth$k - rr.pm2.5.asth$p),
                                             "; ", I^2, " = ",.(formatC(rr.pm2.5.asth$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.pm2.5.asthma$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.pm2.5.asthma$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.pm2.5.asthma$ci.ub, digits=2, format="f")),
                                             "]")))


text(-3.5, 10.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.pm2.5.coug$k - rr.pm2.5.coug$p),
                                             "; ", I^2, " = ",.(formatC(rr.pm2.5.coug$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.pm2.5.cough$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.pm2.5.cough$ci.lb, digits=2, format="f")),
                                             " - ",
                                             .(formatC(rr.pm2.5.cough$ci.ub, digits=2, format="f")),
                                             "]")))


text(-3.5, 1.5, pos=4, cex=0.70, bquote(paste("RE Model (df = ", .(rr.pm2.5.so$k - rr.pm2.5.so$p),
                                            "; ", I^2, " = ",.(formatC(rr.pm2.5.so$I2, digits=1, format="f")), "%),  ", "Pooled risk ratio = ", .(formatC(rr.pm2.5.sob$pred, digits=2, format="f")), " [95% CI ",.(formatC(rr.pm2.5.sob$ci.lb, digits=2, format="f")),
                                            " - ",
                                            .(formatC(rr.pm2.5.sob$ci.ub, digits=2, format="f")),
                                            "]")))

dev.off()