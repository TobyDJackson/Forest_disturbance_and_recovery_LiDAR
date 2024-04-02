# Working on this December 2023

library(tidyverse); library(magrittr)
library(readxl); library(MetBrewer)
library(tidytext); library(patchwork)
library(corrmorant); library(ggfortify)
library(broom.mixed) # load broom.mixed
library(lme4)
library(ggforestplot)


# Preparing the data ######
dropbox=paste0("C:/Users/av23907/Dropbox/")
pixels = readr::read_csv(file = paste0(dropbox,"3_ALS_summary_rasters/summary_data/allsites_allparams_tile_summary_by_class.csv")) %>% dplyr::filter(percentage_notNA>95)


# Convert to annual change rates
pixels %<>% dplyr::mutate(volume_increase_m3_pha_py = area*height_increase_mean/(interval*tilesize), 
                       proportion_of_area = area/(percentage_notNA*tilesize),
                       height_increase_m_py = height_increase_mean/(interval)) 

pixels %<>% dplyr::filter(site != "D13", site != "D14")

pixels$class %<>% factor(levels = c("Disturbance - canopy","Disturbance - new gap", "Gap - persistent","Gap - recovered","Intact canopy" ))

pixels$site %<>% factor(levels = (c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum" )))
pixels$site_height=NA
pixels$site_height[pixels$site=="Danum"]="Danum (68)"
pixels$site_height[pixels$site=="Sepilok"]="Sepilok (64)"
pixels$site_height[pixels$site=="Nouragues"]="Nouragues (50)"
pixels$site_height[pixels$site=="Tapajos"]="Tapajos (49)"
pixels$site_height[pixels$site=="Paracou"]="Paracou (42)"
pixels$site_height[pixels$site=="Ducke"]="Ducke (39)"
pixels$site_height %<>% factor(levels = (c("Ducke (39)","Paracou (42)","Tapajos (49)","Nouragues (50)","Sepilok (64)","Danum (68)" )))

summary(as.factor(pixels$param))


pixels_per_site_by_class = pixels  %>% 
  dplyr::filter(param == "f") %>%
  dplyr::group_by( raster_algorithm, site,site_height, interval, class) %>% 
  summarize(n=n(),mean_area =mean(area,na.rm=T),
            proportion_of_area =mean(proportion_of_area,na.rm=T),
            height_increase_m_py = mean(height_increase_m_py*area,na.rm=T)/mean_area, 
            volume_increase_m3_pha_py = mean(volume_increase_m3_pha_py,na.rm=T))



# Recurrence times table #########
recurrence_times  = pixels_per_site_by_class %>% 
  dplyr::group_by(site) %>% 
  dplyr::mutate(undisturbed_area = 100-proportion_of_area[class=="Gap - recovered"]-proportion_of_area[class=="Gap - persistent"], 
                recurrence_time  =  interval*undisturbed_area/proportion_of_area[class=="Disturbance - new gap"]) %>%
  summarize(recurrence_time = recurrence_time[1])

intact_growth  = pixels_per_site_by_class %>% dplyr::filter(class=="Intact canopy")



# Fig 3. Disturbance by site ##########

# Volume change 
pixels_per_site_by_class$class %<>% factor(levels = c("Disturbance - canopy","Disturbance - new gap", "Intact canopy","Gap - persistent","Gap - recovered" ))
(p1 = ggplot(arrange(pixels_per_site_by_class,class)
             , aes(y= site_height,x=volume_increase_m3_pha_py ,fill=class))+
    geom_bar(stat="identity")+
    geom_vline(xintercept=0,color="black")+
    theme_bw()+scale_fill_manual(values=met.brewer("Archambault", 7)[c(6,5,1,3,2)])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    #theme(legend.title = element_blank(),legend.position = "none")+
    ylab("")+xlab(expression(~ Canopy ~ volume ~ change ~ m^3 ~ ha^-1 ~ yr^-1))
)



# Area 
pixels_per_site_by_class$class %<>% factor(levels = rev(c("Disturbance - canopy","Disturbance - new gap", "Gap - recovered","Gap - persistent" ,"Intact canopy")))
(p2 = ggplot(arrange(pixels_per_site_by_class,class)
             #%>% dplyr::filter(class != "Intact canopy")
             ,aes(y= site,x=proportion_of_area,fill=class))+geom_bar(stat="identity")+
    theme_bw()+scale_fill_manual(values=met.brewer("Archambault", 7)[rev(c(6,5,2,3,1))])+
    theme(legend.position = "right",legend.title = element_blank())+
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ylab("")+xlab(paste0("Area (%)"))+guides(fill = guide_legend(reverse=TRUE))
)



combined =  p1+p2 +plot_layout(guides = "collect")+plot_annotation(tag_levels = 'A')+
  guides(fill = guide_legend(reverse=TRUE)) & theme(legend.position = "right") 
combined

ggsave(plot=combined, filename = paste0(dropbox,"Rscripts/ChangeDetection/CD_figures/",
                                        "Fig3.Mean_10ha_2panel_lspikefree_fixed_dCHM_5classes.png"),
       height = 3, width = 8)
# NOTE - I couldn't figure out the legend, so I did it manually in powerpoint.


# __ S3: sensitivity to raster algorithm ##############


pixels %<>% mutate(parameters = paste0(tilesize,"ha, ",height_threshold_type,", min ", min_area_threshold,"m2, ",change_layer))
summary(as.factor(pixels$parameters))
pixels_ra = pixels %>% dplyr::filter(parameters =="10ha, fixedH, min 25m2, chm", raster_algorithm !="tspikefree2",raster_algorithm !="thighest2",)
summary(as.factor(pixels_ra$raster_algorithm))
pixels_ra$raster_algorithm %<>% factor(levels = c("lspikefree", "highest","tin" ))

pixels_per_site_by_class = pixels_ra  %>%
  dplyr::group_by(raster_algorithm, site,site_height, interval, class) %>% 
  summarize(n=n(),mean_area =mean(area,na.rm=T),
            proportion_of_area =mean(proportion_of_area,na.rm=T),
            height_increase_m_py = mean(height_increase_m_py*area,na.rm=T)/mean_area, 
            volume_increase_m3_pha_py = mean(volume_increase_m3_pha_py,na.rm=T))

pixels_per_site_by_class$class %<>% factor(levels = c("Disturbance - canopy","Disturbance - new gap", "Intact canopy","Gap - persistent","Gap - recovered" ))



# Volume change 
pixels_per_site_by_class$class %<>% factor(levels = c("Disturbance - canopy","Disturbance - new gap", "Intact canopy","Gap - persistent","Gap - recovered" ))
(p1 = ggplot(arrange(pixels_per_site_by_class,class)
             , aes(y= site,x=volume_increase_m3_pha_py ,fill=class))+
    geom_bar(stat="identity")+
    facet_wrap(~raster_algorithm,nrow=3)+
    geom_vline(xintercept=0,color="black")+
    theme_bw()+scale_fill_manual(values=met.brewer("Archambault", 7)[c(6,5,1,3,2)])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.title = element_blank(),legend.position = "right")+
    ylab("")+xlab(expression(~ Volume ~ change ~ m^3 ~ ha^-1 ~ yr^-1))
)


# Area 
pixels_per_site_by_class$class %<>% factor(levels = rev(c("Disturbance - canopy","Disturbance - new gap", "Gap - recovered","Gap - persistent","Intact canopy" )))
(p2 = ggplot(arrange(pixels_per_site_by_class,class)
             #%>% dplyr::filter(class != "Intact canopy")
             ,aes(y= site,x=proportion_of_area,fill=class))+geom_bar(stat="identity")+
    facet_wrap(~raster_algorithm,nrow=3)+
    theme_bw()+scale_fill_manual(values=met.brewer("Archambault", 7)[rev(c(6,5,2,3,1))])+
    theme(legend.position = "none",legend.title = element_blank())+
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ylab("")+xlab(paste0("Area (%)"))
)



combined =  p1+p2 +plot_layout(guides = "collect")+
  guides(fill = guide_legend(reverse=TRUE)) & theme(legend.position = "right") 
combined

ggsave(plot=combined, filename = paste0(dropbox,"Rscripts/ChangeDetection/CD_figures/",
                                        "FigS3.Volume_change_by_5classes_sensitivity_raster_algorithm.png"),
       height = 8, width = 8)



# __ S4: sensitivity to analysis decisions ##############

pixels %<>% mutate(parameters = paste0(tilesize,"ha, ",height_threshold_type,", min ", min_area_threshold,"m2, ",change_layer))
summary(as.factor(pixels_ad$parameters))

pixels_ad = pixels %>% dplyr::filter(raster_algorithm=="lspikefree")
summary(as.factor(pixels_ad$parameters))


pixels_ad$parameter_changed = NA
pixels_ad$parameter_changed[pixels_ad$parameters=="10ha, fixedH, min 10m2, chm"] = "min area = 10 m2"
pixels_ad$parameter_changed[pixels_ad$parameters=="10ha, fixedH, min 25m2, dsm"] = "height change without ground detection"
pixels_ad$parameter_changed[pixels_ad$parameters=="10ha, meanH, min 25m2, chm"] = "relative height thresholds"
pixels_ad$parameter_changed[pixels_ad$parameters=="1ha, fixedH, min 25m2, chm"] = "1 ha tiles"
pixels_ad$parameter_changed %<>% factor(levels = c("1 ha tiles", "height change without ground detection","relative height thresholds","min area = 10 m2" ))
pixels_ad %<>% dplyr::filter(is.na(parameter_changed)==0)

pixels_per_site_by_class = pixels_ad  %>%
  dplyr::group_by(parameter_changed, site,site_height, interval, class) %>% 
  summarize(n=n(),mean_area =mean(area,na.rm=T),
            proportion_of_area =mean(proportion_of_area,na.rm=T),
            height_increase_m_py = mean(height_increase_m_py*area,na.rm=T)/mean_area, 
            volume_increase_m3_pha_py = mean(volume_increase_m3_pha_py,na.rm=T))

# Volume change 
pixels_per_site_by_class$class %<>% factor(levels = c("Disturbance - canopy","Disturbance - new gap", "Intact canopy","Gap - persistent","Gap - recovered" ))
(p1 = ggplot(arrange(pixels_per_site_by_class,class)
             , aes(y= site,x=volume_increase_m3_pha_py ,fill=class))+
    geom_bar(stat="identity")+
    facet_wrap(~parameter_changed,nrow=4)+
    geom_vline(xintercept=0,color="black")+
    theme_bw()+scale_fill_manual(values=met.brewer("Archambault", 7)[c(6,5,1,3,2)])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.title = element_blank(),legend.position = "right")+
    ylab("")+xlab(expression(~ Volume ~ change ~ m^3 ~ ha^-1 ~ yr^-1))
)


# Area 
pixels_per_site_by_class$class %<>% factor(levels = rev(c("Disturbance - canopy","Disturbance - new gap", "Gap - recovered","Gap - persistent","Intact canopy" )))
(p2 = ggplot(arrange(pixels_per_site_by_class,class)
             ,aes(y= site,x=proportion_of_area,fill=class))+geom_bar(stat="identity")+
    facet_wrap(~parameter_changed,nrow=4)+
    theme_bw()+scale_fill_manual(values=met.brewer("Archambault", 7)[rev(c(6,5,2,3,1))])+
    theme(legend.position = "none",legend.title = element_blank())+
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ylab("")+xlab(paste0("Area (%)"))
)



combined =  p1+p2+plot_layout(guides = "collect")+
  guides(fill = guide_legend(reverse=TRUE)) & theme(legend.position = "right")  
combined

ggsave(plot=combined, filename = paste0(dropbox,"Rscripts/ChangeDetection/CD_figures/",
                                        "FigS4.Volume_change_by_5classes_sensitivity_analysis_decision.png"),
       height = 8, width = 8)


