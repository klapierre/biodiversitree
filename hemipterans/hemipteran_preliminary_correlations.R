################################################################################
##  hemipteran_preliminary_correlations.R: Preliminary data investigation for 
##  biodiversitree multi-trophic project on hemipteran abundance data.
##
##  Authors: Kimberly Komatsu & Kelsey McGurrin
##  Date created: February 21, 2024
################################################################################

library(PerformanceAnalytics)
library(readxl)
library(tidyverse)

setwd('G:\\Shared drives\\BDT SEM')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=26),
             axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=26),
             plot.title = element_text(size=54, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=30))

customBDTcolors<-scale_color_manual(values=c("ACRU"="#798df4","CAAL"="#ff7072","CACA"="#00c5ff","COFL"="#fde9aa",
                                             "FAGR"="#b80073","FRPE"="#1801fd","LIST"="#ffff05","LITU"="#14d600",
                                             "NYSY"="salmon","PLOC"="#0076f4","QUAL"="#791972","QUPA"="#ff0000",
                                             "QURU"="#9c001d","QUVE"="#af804e","ULAM"="#9478ff"))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  




##### read in data #####

data2022 <- read_xlsx('2022\\input_data\\BDT_hemip_data_2022.xlsx')  
 
# may want to drop 3 individuals where no confirmed ID/"unknown" trophic guild for specific analyses- 
# leaving rows in because they should count toward abundance and family diversity
  

##### figures #####

#hemipteran family patterns across tree spp and div treatments
familyTreeAbundance <- data2022 %>% 
  filter(abundance>0, !is.na(family)) %>% 
  group_by(tree_spp, family) %>% 
  summarise(sum_abundance=sum(abundance)) %>% 
  mutate(proportion = round((sum_abundance/sum(sum_abundance)), digits=3)) %>% 
  ungroup() %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

ggplot(familyTreeAbundance, aes(x="", y=proportion, fill=family)) +
  geom_col() +
  coord_polar(theta="y") +
  facet_wrap(~tree_spp) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5), 
        strip.text = element_text(size = 20))

familyDiversityAbundance <- data2022 %>% 
  filter(abundance>0, !is.na(family)) %>% 
  group_by(tree_div, family) %>% 
  summarise(sum_abundance=sum(abundance)) %>% 
  mutate(proportion = round((sum_abundance/sum(sum_abundance)), digits=3)) %>% 
  ungroup() %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

ggplot(familyDiversityAbundance, aes(x="", y=proportion, fill=family)) +
  geom_col() +
  coord_polar(theta="y") +
  facet_wrap(~tree_div) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5), 
        strip.text = element_text(size = 20))

#hemipteran functional group patterns across tree spp and div treatments
functionalTreeAbundance <- data2022 %>% 
  filter(abundance>0, !is.na(trophic_guild)) %>% 
  group_by(tree_spp, trophic_guild) %>% 
  summarise(sum_abundance=sum(abundance)) %>% 
  mutate(proportion = round((sum_abundance/sum(sum_abundance)), digits=3)) %>% 
  ungroup() %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

ggplot(functionalTreeAbundance, aes(x="", y=proportion, fill=trophic_guild)) +
  geom_col() +
  coord_polar(theta="y") +
  facet_wrap(~tree_spp) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5), 
        strip.text = element_text(size = 20))

functionalDiversityAbundance <- data2022 %>% 
  filter(abundance>0, !is.na(trophic_guild)) %>% 
  group_by(tree_div, trophic_guild) %>% 
  summarise(sum_abundance=sum(abundance)) %>% 
  mutate(proportion = round((sum_abundance/sum(sum_abundance)), digits=3)) %>% 
  ungroup() %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

ggplot(functionalDiversityAbundance, aes(x="", y=proportion, fill=trophic_guild)) +
  geom_col() +
  coord_polar(theta="y") +
  facet_wrap(~tree_div) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5), 
        strip.text = element_text(size = 20))


#abundances
ggplot(data=subset(data2022, !is.na(trophic_guild)), aes(x=as.factor(tree_div), y=abundance, color=trophic_guild)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~tree_spp, scales='free')
  
####  hemip abund standard BDT format  ####  
ggplot(data=data2022,aes(x=tree_div, y=abundance, color=tree_spp, group=tree_spp))+
  stat_smooth(method="glm",se=FALSE,linetype=5,linewidth=1)+ 
  stat_smooth(aes(group=1), method="glm",se=T, colour="black", linewidth=2)+
  stat_summary(aes(group=1),colour="black",size=0.85,fun.data = "mean_se")+
  labs(x="Plot Diversity",y="Hemipteran Abundance",title = "2022")+
  customBDTcolors+
  theme_classic()+
  theme(strip.text=element_text(face="bold",size=rel(1.5)))+
  scale_x_continuous(limits=c(0, 13),
                     breaks=c(1,4,12),
                     expand = c(0,0))
#ggsave("2022/hemipteran results/hemip abundance 2022.png",width=4, height=4, units=c("in"))


##### merge on other datasets #####

damage <- read_xlsx('2022\\input_data\\2022_tree_damage.xlsx') %>% 
  mutate(damage=((leaf_1+leaf_2+leaf_3+leaf_4+leaf_5+leaf_6+leaf_7+leaf_8+leaf_9+leaf_10+leaf_11+leaf_12)/12)) %>% 
  select(-(leaf_1:leaf_12), -date) %>% 
  group_by(plot, tree_spp, tree_indiv, metric, pct_deer) %>% 
  summarise(damage_mean=mean(damage)) %>%
  ungroup() %>% 
  pivot_wider(names_from=metric, values_from=damage_mean, values_fill=0)

lep <- read.csv('2022\\input_data\\2022 lep info.csv') %>% 
  group_by(plot, tree_div, tree_spp, tree_indiv) %>% 
  summarise(lep_richness=mean(tot_rich),
            lep_abund=sum(tot_abund)) %>% 
  ungroup()

treeGrowth <- read.csv('2022\\input_data\\tree_growth_22.csv') %>% 
  rename(tree_indiv=indiv)

functionalAbundance <- data2022 %>% 
  filter(abundance>0) %>% 
  group_by(tree_spp, tree_div, plot, tree_indiv, trophic_guild) %>% 
  summarise(abund=sum(abundance)) %>% 
  ungroup() %>% 
  pivot_wider(names_from=trophic_guild, values_from=abund, values_fill=0)
  
allData <- data2022 %>% 
  filter(abundance>0) %>% 
  group_by(tree_spp, tree_div, plot, tree_indiv) %>% 
  summarise(hemip_fam_rich=length(unique(family)),
            hemip_guild_rich=length(unique(trophic_guild)),
            hemip_abund=sum(abundance)) %>% 
  ungroup() %>% 
  left_join(functionalAbundance) %>% 
  left_join(lep) %>% 
  left_join(damage) %>%  
  left_join(treeGrowth) %>% 
  select(-'NA', -unknown)

chartData <- allData %>% 
  select(-tree_spp, -plot, -tree_indiv, -ht_growth, -RCD_growth, -rad_growth, -trunk_growth)

chart.Correlation(chartData, histogram=TRUE, pch=19)


#damage
trt <- data2022 %>% 
  select(tree_spp, tree_div, plot, tree_indiv) %>% 
  unique()

damage2 <- damage %>% 
  pivot_longer(pct_deer:pct_path, names_to='metric', values_to='damage') %>% 
  left_join(trt) %>% 
  filter(!is.na(damage), !is.na(tree_div))

ggplot(data=damage2, aes(x=as.factor(tree_div), y=damage, color=tree_spp)) +
  geom_boxplot() +
  facet_wrap(~metric, scales='free') +
  theme(strip.text = element_text(size = 20))
  
####  damage standard BDT format  ####  

# add tree_div column
damage <- damage %>% 
  left_join(trt) 

## deer
ggplot(data=damage,aes(x=tree_div, y=pct_deer, color=tree_spp, group=tree_spp))+
  stat_smooth(method="glm",se=FALSE,linetype=5,linewidth=1)+ 
  stat_smooth(aes(group=1), method="glm",se=T, colour="black", linewidth=2)+
  stat_summary(aes(group=1),colour="black",size=0.85,fun.data = "mean_se")+
  labs(x="Plot Diversity",y="% of Tree Damaged by Deer",title = "2022")+
  customBDTcolors+
  theme_classic()+
  theme(strip.text=element_text(face="bold",size=rel(1.5)))+
  scale_x_continuous(limits=c(0, 13),
                     breaks=c(1,4,12),
                     expand = c(0,0))
#ggsave("2022/damage results/deer damage 2022.png",width=4, height=4, units=c("in"))

## galls
ggplot(data=damage,aes(x=tree_div, y=num_gall, color=tree_spp, group=tree_spp))+
  stat_smooth(method="glm",se=FALSE,linetype=5,linewidth=1)+ 
  stat_smooth(aes(group=1), method="glm",se=T, colour="black", linewidth=2)+
  stat_summary(aes(group=1),colour="black",size=0.85,fun.data = "mean_se")+
  labs(x="Plot Diversity",y="Average # Galls Per Leaf",title = "2022")+
  customBDTcolors+
  theme_classic()+
  theme(strip.text=element_text(face="bold",size=rel(1.5)))+
  scale_x_continuous(limits=c(0, 13),
                     breaks=c(1,4,12),
                     expand = c(0,0))
#ggsave("2022/damage results/avg num galls 2022.png",width=4, height=4, units=c("in"))

## insect
ggplot(data=damage,aes(x=tree_div, y=pct_insect, color=tree_spp, group=tree_spp))+
  stat_smooth(method="glm",se=FALSE,linetype=5,linewidth=1)+ 
  stat_smooth(aes(group=1), method="glm",se=T, colour="black", linewidth=2)+
  stat_summary(aes(group=1),colour="black",size=0.85,fun.data = "mean_se")+
  labs(x="Plot Diversity",y="Average % Insect Damage Per Leaf",title = "2022")+
  customBDTcolors+
  theme_classic()+
  theme(strip.text=element_text(face="bold",size=rel(1.5)))+
  scale_x_continuous(limits=c(0, 13),
                     breaks=c(1,4,12),
                     expand = c(0,0))
#ggsave("2022/damage results/insect damage 2022.png",width=4, height=4, units=c("in"))

#### hemip richness standard BDT format ####

## by family
ggplot(data=allData,aes(x=tree_div, y=hemip_fam_rich, color=tree_spp, group=tree_spp))+
  stat_smooth(method="glm",se=FALSE,linetype=5,linewidth=1)+ 
  stat_smooth(aes(group=1), method="glm",se=T, colour="black", linewidth=2)+
  stat_summary(aes(group=1),colour="black",size=0.85,fun.data = "mean_se")+
  labs(x="Plot Diversity",y="# Hemipteran Families",title = "2022")+
  customBDTcolors+
  theme_classic()+
  theme(strip.text=element_text(face="bold",size=rel(1.5)))+
  scale_x_continuous(limits=c(0, 13),
                     breaks=c(1,4,12),
                     expand = c(0,0))
#ggsave("2022/hemipteran results/hemip fam rich 2022.png",width=4, height=4, units=c("in"))

## by feeding guild
ggplot(data=allData,aes(x=tree_div, y=hemip_guild_rich, color=tree_spp, group=tree_spp))+
  stat_smooth(method="glm",se=FALSE,linetype=5,linewidth=1)+ 
  stat_smooth(aes(group=1), method="glm",se=T, colour="black", linewidth=2)+
  stat_summary(aes(group=1),colour="black",size=0.85,fun.data = "mean_se")+
  labs(x="Plot Diversity",y="# Hemipteran Guilds",title = "2022")+
  customBDTcolors+
  theme_classic()+
  theme(strip.text=element_text(face="bold",size=rel(1.5)))+
  scale_x_continuous(limits=c(0, 13),
                     breaks=c(1,4,12),
                     expand = c(0,0))
#ggsave("2022/hemipteran results/hemip guild rich 2022.png",width=4, height=4, units=c("in"))