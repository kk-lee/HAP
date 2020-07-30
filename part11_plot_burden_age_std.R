
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


final <- read.csv("data_derived.v2/final_age_std.csv")


##==============================================
##                  PLOTS
##==============================================

#2017 data and joining to a map

dat.seventeen <- final %>% 
  group_by(year, iso3) %>% 
  dplyr::summarise(burdence=sum(burdence), burdenll=sum(burdenll), burdenul=sum(burdenul)) %>% 
  filter(year==2017) %>% 
  ungroup() %>% 
  mutate_at(vars(burdence, burdenll, burdenul), funs(.*1000))

dat.seventeen$burdencat <- cut2(dat.seventeen$burdence, g=10)

pdf(file="results.v2/Age-standardised Burden Cartogram 2017.pdf", height=5, width=10)

data(dat.seventeen)
sPDF <- joinCountryData2Map(dat.seventeen
                            ,joinCode = "ISO3"
                            ,nameJoinColumn = "iso3")
sPDF <- subset(sPDF, continent != "Antarctica")

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
mapParams$legendText <- c("0", "7-10", "10-29", "29-67", "67-260", "260-516", "516-1510", "1510-1944", "1944-2566", "2566-5152")

# do.call(addMapLegendBoxes
#         ,c(mapParams 
#            ,ncol=2
#            ,x='bottomleft'
#            ,title="Age-standardised DALYs (thousands)"))

dev.off()


iso <- readRDS("data_derived/iso.rds")

table.dat.seventeen <- left_join(dat.seventeen, iso) %>% 
  mutate(burden=paste0(round(burdence, digits=0)," [", round(burdenll, digits=0), " - ", round(burdenul, digits=0), "]")) %>% 
  select(country, burden)

write.csv(table.dat.seventeen, "results.v2/country_burden_agestd.table.csv")

