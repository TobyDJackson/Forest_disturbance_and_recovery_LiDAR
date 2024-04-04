

# Written by Toby Jackson for the paper: 
# 'Tall Bornean forests experience higher disturbance rates than eastern Amazonia'
# This script creates figure 2

project_folder= "C:/Users/av23907/Dropbox/3_ALS_summary_rasters/v2_sites_2023/Deposit/"

library(ggridges)
library(tidyverse); library(magrittr)
library(readxl); library(MetBrewer)
library(tidytext); library(patchwork)
library(dplyr)


# Fig 2. Height and change distributions - loop ##########

temp_list=list(); c=1

for (site in c("Tapajos","Ducke","Paracou","Nouragues","Sepilok", "Danum")){
  if (site =="Ducke")    {interval=5.148}
  if (site =="Tapajos")  {interval=4.595}
  if (site =="Paracou")  {interval=4.164}
  if (site =="Nouragues"){interval=4.164}
  if (site =="Sepilok")  {interval=5.315}
  if (site =="Danum")    {interval=5.362}
  if (site =="D13")      {interval=6.258}
  if (site =="D14")      {interval=5.362}
  
  
  
  a_chm     = terra::rast(paste0(project_folder,site,"/a_chm_lspikefree_multi3.1_slope1.75_offset2.1.tif"))
  b_chm     = terra::rast(paste0(project_folder,site,"/b_chm_lspikefree_multi3.1_slope1.75_offset2.1.tif"))
  height_increase = b_chm-a_chm
  whole_stack = terra::rast(list(a_chm,b_chm,height_increase))
  names(whole_stack) = c("a_chm","b_chm","height_increase"  )
  
  temp = whole_stack %>% as.data.frame() %>% dplyr::slice_sample(prop=0.1)
  temp$site = site
  temp$interval = interval
  temp_list[[c]] = temp
  c=c+1
  
}
dt = dplyr::bind_rows(temp_list)


dts = dt  %>% group_by(site) %>% tidyr::drop_na() %>% dplyr::distinct() %>%
  summarize(area_ha = n()/10000,height_increase = mean(height_increase,na.rm=T))

dts %<>% dplyr::mutate(height_increase_py = height_increase/interval)

names(dt)


# __ Plot Fig 2 #################
dt$site %<>% factor(levels = (c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum" )))

p1 = ggplot()+
  geom_density_ridges(data = dl_chm,aes(x=a_chm,y=site_height,fill=site),scale=1.5,color=NA)+
  theme_bw()+theme(legend.position ="none")+#scale_y_discrete(position = "right")+
  scale_fill_manual(values=met.brewer("Klimt", 6)[rev(c(2,1,3,4,5,6))])+
  xlab("Canopy height (m)")+ylab("")+xlim(c(0,80))

p2 = ggplot()+
  geom_density_ridges(data = dl_chm,aes(x=height_increase_py,y=site_height,fill=site),scale=1.5,color=NA)+
  theme_bw()+theme(legend.position ="none")+#scale_y_discrete(position = "right")+
  scale_fill_manual(values=met.brewer("Klimt", 6)[rev(c(2,1,3,4,5,6))])+
  theme(legend.position ="none", axis.text.y = element_blank () , axis.ticks.y = element_blank ())+
  geom_vline(xintercept=0)+
  xlab("Height change (m/yr)")+ylab("")+xlim(c(-5,2))

p2

p3 = p1+p2+plot_annotation(tag_level = "A")
p3

ggsave(plot=p3, filename = paste0(project_folder,"Fig2.Height_change_distributions.png"),
       height = 3, width = 6)
