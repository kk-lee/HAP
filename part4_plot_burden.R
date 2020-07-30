
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


final <- read.csv("data_derived.v2/final.csv")


##==============================================
##                  PLOTS
##==============================================

# plot of all daly

dat <- final %>% 
  group_by(year) %>% 
  dplyr::summarise(burdence=sum(burdence), burdenll=sum(burdenll), burdenul=sum(burdenul)) %>% 
  ungroup()

p1 <- ggplot(data= dat) 
p1 <- p1+
  geom_line(aes(y=dat$burdence, x=dat$year), colour = "blue")+
  geom_line(aes(y=dat$burdenll, x=dat$year), colour = "red",linetype="dashed")+
  geom_line(aes(y=dat$burdenul, x=dat$year), colour = "red",linetype="dashed")+
  scale_x_continuous(name="\nTime (year)")+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,250000))+
  theme(legend.position="none")+
  theme_light()

p1

# plot all daly by income and region #

income <- read.csv("data_raw/worldbank.csv")
income <- income %>% select(Code, Income.group) %>% rename_at("Code",~"iso3")

final <- inner_join(final, income)


dat2 <- final %>% 
  group_by(year, Income.group) %>% 
  dplyr::summarise(burdence=sum(burdence), burdenll=sum(burdenll), burdenul=sum(burdenul)) %>% 
  ungroup()

p2 <- ggplot(data= dat2) 
p2 <- p2+
  geom_line(aes(y=dat2$burdence, x=dat2$year, colour = dat2$Income.group))+
  scale_x_continuous(name="\nTime (year)")+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,100000))+
  theme(legend.position="none")+
  guides(colour=guide_legend(title="Income group"))+
  theme_light()

p2

# stackplot by income

dat2$Income.group <- factor(dat2$Income.group,levels=c("High income", "Upper middle income", "Low income","Lower middle income"))

pdf(file="results.v2/Figure_5A.Burden by income stack plot.pdf", height=8, width=12)

p2b <- ggplot(arrange(dat2,Income.group), aes(x=as.numeric(as.character(year)), y=burdence,order=Income.group,fill=Income.group))
p2b <- p2b + geom_area(alpha=0.9)+
  scale_fill_manual(values=c("#0c2c84", "#7fcdbb", "#c7e9b4","#ffffd9")) +
  scale_x_continuous(name="\nTime (year)", breaks = c(2000, 2005, 2010, 2015, 2017))+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,100000))+
  theme_light()+
  labs(fill="Income group")+
  theme(legend.position = c(0.8, 0.83),legend.text = element_text(size=10))

p2b 

dev.off()


dat3 <- final %>% 
  group_by(year, region) %>% 
  dplyr::summarise(burdence=sum(burdence), burdenll=sum(burdenll), burdenul=sum(burdenul)) %>% 
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
  geom_line(aes(y=dat3$burdence, x=dat3$year, colour = dat3$region))+
  scale_x_continuous(name="\nTime (year)")+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,60000))+
  theme(legend.position="none")+
  guides(colour=guide_legend(title="WHO Region"))+
  theme_light()

p3


# Stack plot by region

dat3$region <- factor(dat3$region,levels=c("European region", "Region of the Americas", "Eastern Mediterranean region", "Western Pacific region","African region", "South-East Asia region"))


pdf(file="results.v2/Figure_5B.Burden by region stack plot.pdf", height=8, width=12)

p4 <- ggplot(arrange(dat3,region), aes(x=as.numeric(as.character(year)), y=burdence,order=region,fill=region))
p4 <- p4 + geom_area(alpha=0.8)+scale_fill_brewer(palette="Reds", direction = 1,labels=c("European region", "Region of the Americas", "Eastern Mediterranean region", "Western Pacific region","African region", "South-East Asia region")) +
  scale_x_continuous(name="\nTime (year)", breaks = c(2000, 2005, 2010, 2015, 2017))+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,100000))+
  theme_light()+
  labs(fill="WHO region")+
  theme(legend.position = c(0.8, 0.83),legend.text = element_text(size=10))

p4 

dev.off()

# Stack plot by disease

dat5 <- final %>% 
  mutate(subgroup = case_when(
    subgroup %in% c("ari", "tb") ~ "Communicable respiratory disease",
    subgroup %in% c("asthma", "copd") ~ "Chronic respiratory disease",
    subgroup %in% c("ihd", "stroke") ~ "Atherosclerotic cardiovascular disease",
    subgroup=="lungca" ~ "Lung cancer",
    TRUE ~ as.character(subgroup)     
  ))

dat5 <- dat5 %>% 
  group_by(year, subgroup) %>% 
  dplyr::summarise(burdence=sum(burdence), burdenll=sum(burdenll), burdenul=sum(burdenul)) %>% 
  ungroup()


dat5$subgroup <- factor(dat5$subgroup,levels=c("Atherosclerotic cardiovascular disease", "Lung cancer", "Chronic respiratory disease", 
                                           "Communicable respiratory disease"))


pdf(file="results.v2/Figure_4.Burden by cause stack plot.pdf", height=8, width=12)

p5 <- ggplot(arrange(dat5,subgroup), aes(x=as.numeric(as.character(year)), y=burdence,order=subgroup,fill=subgroup))
p5 <- p5 + geom_area(alpha=0.8)+scale_fill_brewer(palette="RdPu", direction = 1, labels=c("Atherosclerotic cardiovascular disease", "Lung cancer", "Chronic respiratory disease", "Communicable respiratory disease")) +
  scale_x_continuous(name="\nTime (year)", breaks = c(2000, 2005, 2010, 2015, 2017))+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,100000))+
  theme_light()+
  labs(fill="Disease")+
  theme(legend.position = c(0.8, 0.83),legend.text = element_text(size=10))


p5 

dev.off()

#2017 data and joining to a map

dat.seventeen <- final %>% 
  group_by(year, iso3) %>% 
  dplyr::summarise(burdence=sum(burdence), burdenll=sum(burdenll), burdenul=sum(burdenul)) %>% 
  filter(year==2017) %>% 
  ungroup()

dat.seventeen$burdencat <- cut2(dat.seventeen$burdence, g=10)


data(dat.seventeen)
sPDF <- joinCountryData2Map(dat.seventeen
                            ,joinCode = "ISO3"
                            ,nameJoinColumn = "iso3")
sPDF <- subset(sPDF, continent != "Antarctica")

pdf("results.v2/burden.cartogram.pdf",
    height=5,
    width=10)

par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
#getting colours
colourPalette <- brewer.pal(10,'OrRd')
#plot map
 #create world map shaped window
mapParams <- mapCountryData(sPDF
                            ,nameColumnToPlot="burdencat"
                            ,addLegend=FALSE
                            ,missingCountryCol = 'dark grey'
                            ,colourPalette=colourPalette, mapTitle=''
                            ,catMethod="categorical")
#adding legend
mapParams$legendText <- c("0", "0.01-0.04", "0.04-1.5", "1.5-5", "5-11", "11-28", "28-103", "103-255", "255-513", "513-17264")

do.call(addMapLegendBoxes
        ,c(mapParams 
           ,ncol=2
           ,x='bottomleft'
           ,title="DALYs (thousands)"))

dev.off()

#descriptive calculations by region
dat3.nought <- dat3[dat3$year==2000,]

dat3.seventeen <- dat3 %>% 
  filter(year==2017) %>% 
  rename_at(c("burdence"), ~c("burdence1"))

x <- cbind(dat3.nought, dat3.seventeen) %>% 
  mutate(percentage=((burdence1-burdence)/burdence)*100) %>%
  mutate(percentage.ll=((burdenll.1-burdenll)/burdenll)*100) %>%
  mutate(percentage.ul=((burdenul.1-burdenul)/burdenul)*100) %>%
  mutate(diff=(burdence1-burdence))


#descriptive calculations by country
dat.nought <- final %>% 
  group_by(year, iso3) %>% 
  dplyr::summarise(burdence1=sum(burdence), burdenll1=sum(burdenll), burdenul1=sum(burdenul)) %>% 
  filter(year==2000) %>% 
  ungroup()

y <- cbind(dat.nought, dat.seventeen) %>% 
  mutate(percentage=((burdence-burdence1)/burdence1)*100) %>% 
  filter(iso3 %in% c("CHN", "IND", "USA", "IDN", "BRA", "PAK", "NGA", "BGD", "RUS", "MEX")) %>% 
  mutate(diff=(burdence-burdence1)) %>% 
  mutate(diff.ll=(burdenll-burdenll1)) %>% 
  mutate(diff.ul=(burdenul-burdenul1))
  
# calculations by cause
dat5.nought <- dat5 %>% filter(year==2000)

dat5.seventeen <- dat5 %>% 
  group_by(year, subgroup) %>% 
  dplyr::summarise(burdence1=sum(burdence), burdenll1=sum(burdenll), burdenul1=sum(burdenul)) %>% 
  ungroup() %>% 
  filter(year==2017)

dat5.diff <- left_join(dat5.nought, dat5.seventeen, by="subgroup") %>% 
  mutate(percentage=((burdence1-burdence)/burdence)*100)

# calculations by year
dat.nought <- dat %>% filter(year==2000)

dat.seventeen1 <- dat %>% 
  filter(year==2017) %>% 
  rename.vars(from=c("burdence", "burdenll", "burdenul"),
              to=c("burdence1", "burdenll1", "burdenul1"))

dat.diff <- cbind(dat.nought, dat.seventeen1) %>% 
  mutate(percentagece=((burdence1-burdence)/burdence)*100) %>% 
  mutate(percentagell=((burdenll1-burdenll)/burdenll)*100) %>% 
  mutate(percentageul=((burdenul1-burdenul)/burdenul)*100) 
  


## supplementary tables

table.dat5 <- dat5 %>% 
  mutate(burden=paste0(round(burdence, digits=0)," [", round(burdenll, digits=0), " - ", round(burdenul, digits=0), "]")) %>% 
  select(year, subgroup, burden) %>% 
  spread(key=subgroup, value=burden)

write.csv(table.dat5, "results.v2/cause_burden.table.csv")


table.dat2 <- dat2 %>% 
  mutate(burden=paste0(round(burdence, digits=0)," [", round(burdenll, digits=0), " - ", round(burdenul, digits=0), "]")) %>% 
  select(year, Income.group, burden) %>% 
  spread(key=Income.group, value=burden) %>% 
  select("year", "Low income", "Lower middle income", "Upper middle income", "High income")

write.csv(table.dat2, "results.v2/income_burden.table.csv")


table.dat3 <- dat3 %>% 
  mutate(burden=paste0(round(burdence, digits=0)," [", round(burdenll, digits=0), " - ", round(burdenul, digits=0), "]")) %>% 
  select(year, region, burden) %>% 
  spread(key=region, value=burden) %>% 
  select("year", "African region", "Region of the Americas", "South-East Asia region", "European region", "Eastern Mediterranean region", 
         "Western Pacific region")

write.csv(table.dat3, "results.v2/regional_burden.table.csv")



iso <- readRDS("data_derived/iso.rds")

table.dat.seventeen <- left_join(dat.seventeen, iso) %>% 
  mutate(burden=paste0(round(burdence, digits=0)," [", round(burdenll, digits=0), " - ", round(burdenul, digits=0), "]")) %>% 
  select(country, burden)

write.csv(table.dat.seventeen, "results.v2/country_burden.table.csv")

## countries in north america

north.america <- read.csv("data_raw/north.america.csv") 

burden.north.america <- left_join(north.america, final) %>% 
  filter(!is.na(burdence))

burden.north.america <- burden.north.america %>% 
  group_by(year) %>% 
  dplyr::summarise(burdence=sum(burdence), burdenll=sum(burdenll), burdenul=sum(burdenul)) %>% 
  ungroup()

burden.north.america.nought <- burden.north.america[burden.north.america$year==2000,]

burden.north.america.seventeen <- burden.north.america %>% 
  filter(year==2017) %>% 
  rename_at(c("burdence"), ~c("burdence1"))

x <- cbind(burden.north.america.nought, burden.north.america.seventeen) %>% 
  mutate(percentage=((burdence1-burdence)/burdence)*100) %>% 
  mutate(diff=(burdence1-burdence))

