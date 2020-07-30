
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
library(ggplot2)
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


final_death <- read.csv("data_derived.v2/final_death.csv")


##==============================================
##                  PLOTS
##==============================================

# plot of deaths 

dat <- final_death %>% 
  group_by(year) %>% 
  dplyr::summarise(deathce=sum(deathce), deathll=sum(deathll), deathul=sum(deathul)) %>% 
  ungroup()

p1 <- ggplot(data= dat) 
p1 <- p1+
  geom_line(aes(y=dat$deathce, x=dat$year), colour = "blue")+
  geom_line(aes(y=dat$deathll, x=dat$year), colour = "red",linetype="dashed")+
  geom_line(aes(y=dat$deathul, x=dat$year), colour = "red",linetype="dashed")+
  scale_x_continuous(name="\nTime (year)")+
  scale_y_continuous(name="Deaths in thousands\n", limits=c(0, 6000), breaks=c(seq(0,6000, by=1000)))+
  theme(legend.position="none")+
  theme_light()

p1

# plot deaths by income and region #

income <- read.csv("data_raw/worldbank.csv")
income <- income %>% select(Code, Income.group) %>% rename_at("Code",~"iso3")

final_death <- inner_join(final_death, income)


dat2 <- final_death %>% 
  group_by(year, Income.group) %>% 
  dplyr::summarise(deathce=sum(deathce), deathll=sum(deathll), deathul=sum(deathul)) %>% 
  ungroup()

p2 <- ggplot(data= dat2) 
p2 <- p2+
  geom_line(aes(y=dat2$deathce, x=dat2$year, colour = dat2$Income.group))+
  scale_x_continuous(name="\nTime (year)")+
  scale_y_continuous(name="Deaths in thousands\n")+
  theme(legend.position="none")+
  guides(colour=guide_legend(title="Income group"))+
  theme_light()

p2



dat2$Income.group <- factor(dat2$Income.group,levels=c("High income", "Upper middle income", "Low income","Lower middle income"))


p2b <- ggplot(arrange(dat2,Income.group), aes(x=as.numeric(as.character(year)), y=deathce,order=Income.group,fill=Income.group))
p2b <- p2b + geom_area(alpha=0.9)+
  scale_fill_manual(values=c("#0c2c84", "#7fcdbb", "#c7e9b4","#ffffd9")) +
  scale_x_continuous(name="\nTime (year)", breaks = c(2000, 2005, 2010, 2015, 2017))+
  scale_y_continuous(name="Deaths in thousands\n", limits=c(0, 3000), breaks=c(seq(0,4000, by=1000)))+
  theme_light()+
  labs(fill="Income group")+
  theme(legend.position = c(0.8, 0.83),legend.text = element_text(size=10))


p2b 

pdf(file="results.v2/Death by income stack plot.pdf", height=8, width=12)
print(p2b)
dev.off()




dat3 <- final_death %>% 
  group_by(year, region) %>% 
  dplyr::summarise(deathce=sum(deathce), deathll=sum(deathll), deathul=sum(deathul)) %>% 
  ungroup() %>% 
  mutate(region= recode(region,
                        "1_Afr"= "African region", 
                        "2_Amr"= "Region of the Americas",
                        "3_Sear"= "South-East Asia region", 
                        "4_Eur"= "European region", 
                        "5_Emr"= "Eastern Mediterranean region", 
                        "6_Wpr"= "Western Pacific region"))      

p3 <- ggplot(data= dat3) 
p3 <- p3+
  geom_line(aes(y=dat3$deathce, x=dat3$year, colour = dat3$region))+
  scale_x_continuous(name="\nTime (year)")+
  scale_y_continuous(name="Deaths in thousands\n", limits=c(0, 1500))+
  theme(legend.position="none")+
  guides(colour=guide_legend(title="WHO region"))+
  theme_light()

p3



# Stack plot by region

dat3$region <- factor(dat3$region,levels=c("European region", "Region of the Americas", "Eastern Mediterranean region", "Western Pacific region","African region", "South-East Asia region"))


p4 <- ggplot(arrange(dat3,region), aes(x=as.numeric(as.character(year)), y=deathce,order=region,fill=region))
p4 <- p4 + geom_area(alpha=0.8)+scale_fill_brewer(palette="Reds", direction = 1,labels=c("European region", "Region of the Americas", "Eastern Mediterranean region", "Western Pacific region","African region", "South-East Asia region")) +
  scale_x_continuous(name="\nTime (year)", breaks = c(2000, 2005, 2010, 2015, 2017))+
  scale_y_continuous(name="Deaths in thousands\n", limits=c(0,3000))+
  theme_light()+
  labs(fill="WHO region")+
  theme(legend.position = c(0.8, 0.83),legend.text = element_text(size=10))


p4 


pdf(file="results.v2/Death by region stack plot.pdf", height=8, width=12)
print(p4)
dev.off()


# plot of deaths by cause
dat5 <- final_death %>% 
  group_by(year, subgroup) %>% 
  dplyr::summarise(deathce=sum(deathce), deathll=sum(deathll), deathul=sum(deathul)) %>% 
  ungroup() %>% 
  mutate(subgroup= recode(subgroup,
                        "cv death"= "Cardiovascular mortality", 
                        "resp death"= "Respiratory mortality",
                        "ufive"= "Under five mortality"))      

p5 <- ggplot(data= dat5) 
p5 <- p5+
  geom_line(aes(y=dat5$deathce, x=dat5$year, colour = dat5$subgroup))+
  scale_x_continuous(name="\nTime (year)")+
  scale_y_continuous(name="Deaths in thousands\n", limits=c(0, 2000))+
  theme(legend.position="none")+
  guides(colour=guide_legend(title="Cause of death"))+
  theme_light()

p5

# 
# p6 <- ggplot(data <- dat4, aes(y=dat4$deathce, x=dat4$year, ymin=dat4$deathll, ymax=dat4$deathul))+
#   geom_ribbon(aes(fill=dat4$subgroup), alpha = 0.2)+
#   geom_line(aes(color=dat4$subgroup))+
#   scale_x_continuous(name="\nTime (year)")+
#   scale_y_continuous(name="Deaths in thousands\n", limits=c(0,3000),breaks=seq(0,3000,500),labels=c("0","500","1,000","1,500","2,000","2,500", "3000"))+
#   theme(legend.position="none")+
#   guides(colour=guide_legend(title="Cause of death"), fill=FALSE)+
#   theme_light()
# 
# p6


# Stack plot by cause of death
dat4 <- final_death %>% 
  group_by(year, subgroup) %>% 
  dplyr::summarise(deathce=sum(deathce), deathll=sum(deathll), deathul=sum(deathul)) %>% 
  ungroup() %>% 
  mutate(subgroup= recode(subgroup,
                          "cv death"= "Cardiovascular mortality", 
                          "resp death"= "Respiratory mortality",
                          "ufive"= "Under five mortality"))      

dat4$subgroup <- factor(dat4$subgroup,levels=c("Cardiovascular mortality", "Under five mortality", "Respiratory mortality"))


p7 <- ggplot(arrange(dat4,subgroup), aes(x=as.numeric(as.character(year)), y=deathce,order=subgroup,fill=subgroup))
p7 <- p7 + geom_area(alpha=0.8)+scale_fill_brewer(palette="RdPu", direction = 1, labels=c("Cardiovascular mortality", "Under five mortality", "Respiratory mortality")) +
  scale_x_continuous(name="\nTime (year)", breaks = c(2000, 2005, 2010, 2015, 2017))+
  scale_y_continuous(name="Deaths in thousands\n", limits=c(0,3000))+
  theme_light()+
  labs(fill="Cause of death")+
  theme(legend.position = c(0.8, 0.83),legend.text = element_text(size=10))


p7 


pdf(file="results.v2/Death by cause stack plot.pdf", height=8, width=12)
print(p7)
dev.off()



#descriptive calculations by country
dat.nought <- final_death %>% 
  group_by(year, iso3) %>% 
  dplyr::summarise(deathce1=sum(deathce), deathll1=sum(deathll), deathul1=sum(deathul)) %>% 
  filter(year==2000) %>% 
  ungroup()

dat.seventeen <- final_death %>% 
  group_by(year, iso3) %>% 
  dplyr::summarise(deathce=sum(deathce), deathll=sum(deathll), deathul=sum(deathul)) %>% 
  filter(year==2017) %>% 
  ungroup()

x <- cbind(dat.nought, dat.seventeen) %>% 
  mutate(percentage=((deathce1-deathce)/deathce1)*100) %>% 
  filter(iso3 %in% c("CHN", "IND", "USA", "IDN", "BRA", "PAK", "NGA", "BGD", "RUS", "MEX")) %>% 
  mutate(diff=(deathce1-deathce))


# calculations by cause
dat5.nought <- dat5 %>% filter(year==2000)

dat5.seventeen <- dat5 %>% 
  group_by(year, subgroup) %>% 
  dplyr::summarise(deathce1=sum(deathce), deathll1=sum(deathll), deathul1=sum(deathul)) %>% 
  ungroup() %>% 
  filter(year==2017)

dat5.diff <- left_join(dat5.nought, dat5.seventeen, by="subgroup") %>% 
  mutate(percentage=((deathce1-deathce)/deathce)*100)

# calculations by year
dat.nought <- dat %>% filter(year==2000)

dat.seventeen1 <- dat %>% 
  filter(year==2017) %>% 
  rename.vars(from=c("deathce", "deathll", "deathul"),
              to=c("deathce1", "deathll1", "deathul1"))

dat.diff <- cbind(dat.nought, dat.seventeen1) %>% 
  mutate(percentagece=((deathce1-deathce)/deathce)*100) %>% 
  mutate(percentagell=((deathll1-deathll)/deathll)*100) %>% 
  mutate(percentageul=((deathul1-deathul)/deathul)*100) 



iso <- readRDS("data_derived/iso.rds")

table.dat.seventeen <- left_join(dat.seventeen, iso) %>% 
  mutate(death=paste0(round(deathce, digits=1)," [", round(deathll, digits=1), " - ", round(deathul, digits=1), "]")) %>% 
  select(country, death)

write.csv(table.dat.seventeen, "results.v2/country_death.table.csv")

