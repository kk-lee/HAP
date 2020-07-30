
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

###====================================
### Part 1 -  estimate pooled rate and risk ratio
###====================================

# load files

dat <- read.csv("data_raw/iap.baseline.v8.csv",header = T,na.strings = "")

names(dat) <- tolower(names(dat))
dat$author <- paste0(dat$author," ","et al")


## high rob==99, low rob==1, moderate rob==2

dat2 <- dat %>% 
  mutate(endpoint.rob=case_when(
    grepl('self reported', outcome.measurement)~99,
    is.na(outcome.measurement)~99,
    TRUE~1
  )) %>% 
  select("study.id", "author", "pub.year", "location.country", "study.population", "study.design", "sample.size", "male.percentage", "adult.paediatric.all", "exposed.group", "control.group", "outcome", "endpoint.rob")

dat2 <- dat %>% 
  select(matches(".adjustment"))%>% 
  mutate_if(is.factor, as.numeric) %>% 
  mutate_all(funs(case_when(
    .==1~0,
    .==2~1,
    .==3~NA_real_
  ))) %>% 
  mutate(num.adjust=rowSums(.)) %>% 
  mutate(age.sex=age.adjustment+sex.adjustment) %>% 
  mutate(confounder.rob=case_when(
    age.sex==2 & num.adjust>3~1,
    age.sex==2 & num.adjust<=3|age.sex==1~2,
    TRUE~99
  )) %>% 
  select(confounder.rob) %>% 
  cbind(dat2,.) %>% 
  mutate(sum.rob=confounder.rob+endpoint.rob) %>% 
  mutate(overall.rob=case_when(
    sum.rob==2~"low RoB",
    sum.rob>99~"high RoB",
    TRUE~"moderate RoB"
  )) %>% 
  select(-c(confounder.rob, endpoint.rob, sum.rob))

write.csv(dat2, "results.v2/baseline.table.csv")


## list of countries table ##
final <- read.csv("data_derived.v2/final.csv") %>% 
  select(country, iso3, region) %>% 
  distinct(country, .keep_all=TRUE)

income <- read.csv("data_raw/worldbank.csv") %>% 
  select(Code, Income.group) %>% 
  rename_at("Code",~"iso3") %>% 
  mutate(Income.group=factor(Income.group, levels=c("Low income", "Lower middle income", "Upper middle income", "High income")))

country.list <- inner_join(final, income) %>% 
  dplyr::arrange(region, Income.group)

write.csv(country.list, "results.v2/country.list.csv")


## number of studies evaluating exposure to kerosene

kerosene <- dat %>% 
  filter(str_detect(exposed.group, "kerosene")| str_detect(exposed.definition, "kerosene"))



