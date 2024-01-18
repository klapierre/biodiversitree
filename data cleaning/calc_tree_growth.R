#######################################################################################
##  calc_tree_growth.R: combines tree size data from spring 2022 and 2023, 
##                  then calculates height & RCD growth over summer 2022
##
##  Author: Kelsey McGurrin
##
#######################################################################################

####setup####
library(tidyverse)

# working directory path (add yours if different)
setwd("G:/Shared drives/BDT SEM/2022")

# input
data22<-read_csv("input_data/2022_trophic_tree_final.csv",col_types = cols(Notes2022 = col_skip(), 
                                                                           ht = col_character(), 
                                                                           measDate = col_skip(),
                                                                           crn_ht=col_character()), 
                 na = c("NA","n/a"))
data23<-read_csv("input_data/2023_trophic_tree_final.csv",col_types = cols(Notes2023 = col_skip(), 
                                                                           ht = col_character(), 
                                                                           measDate = col_skip(),
                                                                           crn_ht=col_character()), 
                 na = c("NA","n/a"))

#### clean up formats ####
#comma separate columns and make numeric
#convert each diameter to basal area and then add multiple trunks together

## 2022
data22<-data22 %>%
  separate(DBH,c("DBH1","DBH2","DBH3","DBH4","DBH5"),sep = ",",remove=TRUE,fill= "right",convert = T) %>%
  separate(base,c("base1", "base2", "base3","base4"),sep = ",", remove = TRUE,fill = "right",convert = T) %>%
  separate(ht ,c("ht1", "ht2"), sep = ",", remove = TRUE,fill = "right",convert = T) %>%
  separate(crn_ht ,c("crn_ht1", "crn_ht2"), sep = ",", remove = TRUE,fill = "right",convert = T)

data22[is.na(data22)] <- 0
data22$DBHBasalArea<- apply(data22[,c('DBH1', 'DBH2', 'DBH3', 'DBH4', 'DBH5')], 1, function(x) { sum(0.00007854*x^2) } )
data22$RCDBasalArea<- apply(data22[,c('base1', 'base2', 'base3', 'base4')], 1, function(x) { sum(0.00007854*x^2) } )

## 2023
data23<-data23 %>%
  separate(DBH,c("DBH1","DBH2","DBH3","DBH4","DBH5"),sep = ",",remove=TRUE,fill= "right",convert = T) %>%
  separate(base,c("base1", "base2", "base3","base4"),sep = ",", remove = TRUE,fill = "right",convert = T) %>%
  separate(ht ,c("ht1", "ht2"), sep = ",", remove = TRUE,fill = "right",convert = T) %>%
  separate(crn_ht ,c("crn_ht1", "crn_ht2"), sep = ",", remove = TRUE,fill = "right",convert = T)
data23$DBH1<-as.numeric(data23$DBH1)
data23[is.na(data23)] <- 0
data23$DBHBasalArea<- apply(data23[,c('DBH1', 'DBH2', 'DBH3', 'DBH4', 'DBH5')], 1, function(x) { sum(0.00007854*x^2) } )
data23$RCDBasalArea<- apply(data23[,c('base1', 'base2', 'base3', 'base4')], 1, function(x) { sum(0.00007854*x^2) } )

#### single year size calcs ####

## 2022
data22 <- data22 %>%
  mutate(crn_htAvg=ifelse(data22$crn_ht2<1, data22$crn_ht1, (data22$crn_ht1 +   data22$crn_ht2) /2),
         htAvg = ifelse(data22$ht2<1, data22$ht1, (data22$ht1 + data22$ht2) /2),
         canHt = (htAvg - crn_htAvg),
         radAvg = ((canNE + canSE) / 2),
         ConeCanVol = (1 / 3) * pi * (radAvg ^ 2) * canHt / 1000000,
         RectCanVol = (4 * canNE * canSE * canHt) / 1000000,
         DBHtrunkVol = .42 * DBHBasalArea * htAvg,
         RCDtrunkVol = .42 * RCDBasalArea * htAvg,
         basalRatio = DBHBasalArea / RCDBasalArea
  )
data22<-select(data22,
               "indiv",
               "ht"="htAvg",
               "radAvg",
               "RCDBasalArea",
               "DBHBasalArea",
               "ConeCanVol",
               "RCDtrunkVol",
               "basalRatio")

## 2023
data23 <- data23 %>%
  mutate(crn_htAvg=ifelse(data23$crn_ht2<1, data23$crn_ht1, (data23$crn_ht1 +   data23$crn_ht2) /2),
         htAvg = ifelse(data23$ht2<1, data23$ht1, (data23$ht1 + data23$ht2) /2),
         canHt = (htAvg - crn_htAvg),
         radAvg = ((canNE + canSE) / 2),
         ConeCanVol = (1 / 3) * pi * (radAvg ^ 2) * canHt / 1000000,
         RectCanVol = (4 * canNE * canSE * canHt) / 1000000,
         DBHtrunkVol = .42 * DBHBasalArea * htAvg,
         RCDtrunkVol = .42 * RCDBasalArea * htAvg,
         basalRatio = DBHBasalArea / RCDBasalArea
  )
data23<-select(data23,
               "indiv",
               "ht"="htAvg",
               "radAvg",
               "RCDBasalArea",
               "DBHBasalArea",
               "ConeCanVol",
               "RCDtrunkVol",
               "basalRatio")

#### growth ####
growth<-full_join(data22,data23,by="indiv",suffix = c(".22",".23"))
growth<-growth %>%
  group_by(indiv) %>%
  transmute(ht_growth=ht.23-ht.22,RCD_growth=RCDBasalArea.23-RCDBasalArea.22)%>%
  drop_na()

ggplot(data=growth,aes(x=ht_growth))+geom_histogram()
ggplot(data=growth,aes(x=RCD_growth))+geom_histogram()

## output
write_csv(growth,file = "input_data/tree_growth_22.csv")
