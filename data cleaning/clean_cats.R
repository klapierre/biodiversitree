#######################################################################################
##  clean_cats.R: filters out 2022 from the "clean" overall cat data,
##                makes columns for each sample and totals, abundance and richness
##
##  Author: Kelsey McGurrin
##
#######################################################################################

####setup####
library(tidyverse)

# working directory path (add yours if different)
setwd("G:/Shared drives/BDT SEM/2022")

# input
all<-read_csv("input_data/all lep info.csv")

data22<-all %>%
  filter(year==2022) %>%
  select(-c(time_samp,add_order,Insect_Order,Insect_Family,Host_Families,num_fam,diet,triplicate))

#### clean ####

## abundance: per sample then pivot wider and add
abund<- data22 %>%
  select(-c(sci_name,lep.code)) 
abund_wide<-pivot_wider(abund,names_from = sample,values_from = indvls,values_fn = sum,names_prefix = "S")
abund_wide$S1[is.na(abund_wide$S1)] <- 0
abund_wide$S2[is.na(abund_wide$S2)] <- 0
abund_wide<-abund_wide %>%
  mutate(tot_abund=S1+S2) %>%
  rename(S1_abund=S1,S2_abund=S2)

## richness: per sample then pivot wider
rich<-data22 %>%
  group_by(sample,tree_indiv) %>%
  drop_na() %>%
  summarize(rich=n())
rich_wide<-pivot_wider(rich,names_from = sample,values_from = rich,values_fn = sum,names_prefix = "S")
rich_wide$S1[is.na(rich_wide$S1)] <- 0
rich_wide$S2[is.na(rich_wide$S2)] <- 0
rich_wide<-rename(rich_wide,S1_rich=S1,S2_rich=S2)

## richness: per year then pivot wider
richy<-data22 %>%
  group_by(tree_indiv) %>%
  drop_na() %>%
  summarize(tot_rich=n())

## combine and output
t<-left_join(abund_wide,rich_wide,by="tree_indiv")
out<-left_join(t,richy,by="tree_indiv")
out[is.na(out)] <- 0

write_csv(out,file = "input_data/2022 lep info.csv")
