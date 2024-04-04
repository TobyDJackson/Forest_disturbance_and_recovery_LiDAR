
# Written by Toby Jackson for the paper: 
# 'Tall Bornean forests experience higher disturbance rates than eastern Amazonia'
# This script creates figures 3-6 from the output of the 'Repeat_LiDAR_disturbance_processing' script

# It analyzes canopy height and height change data repeat LiDAR (Danum, Sepilok, Tapajos, Ducke, Paracou, Nouragues)


library(ggridges)
library(tidyverse); library(magrittr)
library(readxl); library(MetBrewer)
library(tidytext); library(patchwork)
library(corrmorant); library(ggfortify)
library(broom.mixed); library(lme4)
library(ggforestplot); library(dplyr)


project_folder= "C:/Users/av23907/Dropbox/3_ALS_summary_rasters/v2_sites_2023/Deposit/"


# Fig 3. Volume and area balance ###########
# This uses all five classes


# __ Preparing the data ######
project_folder= "C:/Users/av23907/Dropbox/3_ALS_summary_rasters/v2_sites_2023/Deposit/"


df = readr::read_csv(file=paste0(project_folder,"combined_results_030424.csv")) %>% dplyr::filter(percentage_notNA>95)
df$tilesize=10
tilesize=10

# Manually add scanning intervals to data frame
df$interval=NA
df$interval[df$site =="Ducke"]=5.148
df$interval[df$site =="Paracou"]=4.164
df$interval[df$site =="Nouragues"]=4.164
df$interval[df$site =="Sepilok"]=5.315
df$interval[df$site =="Danum"]=5.362
df$interval[df$site =="Tapajos"]=4.595


# Convert to annual change rates
df %<>% dplyr::mutate(volume_increase_m3_pha_py = area*height_increase_mean/(interval*tilesize), 
                          proportion_of_area = area/(percentage_notNA*tilesize),
                          height_increase_m_py = height_increase_mean/(interval)) 

# Manually add 99th percentil of canopy height to site names for plotting
df$site %<>% factor(levels = (c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum" )))
df$site_height=NA
df$site_height[df$site=="Danum"]="Danum (68)"
df$site_height[df$site=="Sepilok"]="Sepilok (64)"
df$site_height[df$site=="Nouragues"]="Nouragues (50)"
df$site_height[df$site=="Tapajos"]="Tapajos (49)"
df$site_height[df$site=="Paracou"]="Paracou (42)"
df$site_height[df$site=="Ducke"]="Ducke (39)"
df$site_height %<>% factor(levels = (c("Ducke (39)","Paracou (42)","Tapajos (49)","Nouragues (50)","Sepilok (64)","Danum (68)" )))


# Summarize disturbance and recovery for each site
df_per_site_by_class = df  %>% 
  dplyr::group_by( site,site_height, interval, class) %>% 
  summarize(n=n(),mean_area =mean(area,na.rm=T),
            proportion_of_area =mean(proportion_of_area,na.rm=T),
            height_increase_m_py = mean(height_increase_m_py*area,na.rm=T)/mean_area, 
            volume_increase_m3_pha_py = mean(volume_increase_m3_pha_py,na.rm=T))



# __ Recurrence times table #########
recurrence_times  = df_per_site_by_class %>% 
  dplyr::group_by(site) %>% 
  dplyr::mutate(undisturbed_area = 100-proportion_of_area[class=="Gap - recovered"]-proportion_of_area[class=="Gap - persistent"], 
                recurrence_time  =  interval*undisturbed_area/proportion_of_area[class=="Disturbance - new gap"]) %>%
  summarize(recurrence_time = recurrence_time[1])


# __ plot Fig 3. Disturbance by site ##########

# Volume change 
df_per_site_by_class$class %<>% factor(levels = c("Disturbance - canopy","Disturbance - new gap", "Intact canopy","Gap - persistent","Gap - recovered" ))
(p1 = ggplot(arrange(df_per_site_by_class,class)
             , aes(y= site_height,x=volume_increase_m3_pha_py ,fill=class))+
    geom_bar(stat="identity")+
    geom_vline(xintercept=0,color="black")+
    theme_bw()+scale_fill_manual(values=met.brewer("Archambault", 7)[c(6,5,1,3,2)])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    #theme(legend.title = element_blank(),legend.position = "none")+
    ylab("")+xlab(expression(~ Canopy ~ volume ~ change ~ m^3 ~ ha^-1 ~ yr^-1))
)


# Area 
df_per_site_by_class$class %<>% factor(levels = rev(c("Disturbance - canopy","Disturbance - new gap", "Gap - recovered","Gap - persistent" ,"Intact canopy")))
(p2 = ggplot(arrange(df_per_site_by_class,class)
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

ggsave(plot=combined, filename = paste0(project_folder,"Fig3.Disturbance_vs_recovery.png"),
       height = 3, width = 8)

# Simplify to three classes ######
df3 = df
df3$class2 = "intact"
df3$class2[df3$class =="Disturbance - canopy" | df3$class=="Disturbance - new gap"]="disturbance"
df3$class2[df3$class =="Gap - persistent" | df3$class=="Gap - recovered"]="gaps"
df3 %<>% dplyr::ungroup() %>% dplyr::group_by(site, interval,
                                                    tile_id, percentage_notNA,tile_mean_height_a,tile_max_height_a, tile_mean_height_increase, 
                                                    tile_elevation_mean,tile_elevation_range, tile_pulsedensity_increase, # tile attributes
                                                    class2) %>%
  dplyr::summarize(CHM1_mean = sum(CHM1_mean*area,na.rm=T)/sum(area,na.rm=T), 
                   height_increase_mean = sum(height_increase_mean*area,na.rm=T)/sum(area,na.rm=T), 
                   area=sum(area,na.rm=T))


# Convert to annual change rates
df3 %<>% dplyr::mutate(volume_increase_m3_pha_py = area*height_increase_mean/(interval*tilesize), 
                        proportion_of_area = area/(percentage_notNA*tilesize),
                        height_increase_m_py = height_increase_mean/(interval)) 


# Pivot wider for modelling
df3 %<>% dplyr::select(!c(CHM1_mean,area,height_increase_mean))
df3_wide = df3 %>% pivot_wider(names_from=class2,values_from=c(volume_increase_m3_pha_py, height_increase_m_py,proportion_of_area))
names(df3_wide)


# Fig 4.  Disturbance vs recovery #####
df3_wide$site %<>% factor(levels = rev(c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum" )))

df3_wide %<>% mutate(net_volume_increase = volume_increase_m3_pha_py_intact+volume_increase_m3_pha_py_gaps+volume_increase_m3_pha_py_disturbance)


p1= ggplot(df3_wide , 
           aes((volume_increase_m3_pha_py_disturbance),
               (volume_increase_m3_pha_py_intact+volume_increase_m3_pha_py_gaps),color = site, group = site))+
  geom_point(size=0.3,alpha=0.5)+#geom_smooth(method="lm",se=F)+
  theme_bw()+theme(legend.title = element_blank())+
  geom_abline(intercept=0,slope=-1,color="black")+
  #xlim(c(0,2))+ylim(c(0,0.5))+
  xlim(c(-10000,0))+ylim(c(0,10000))+
  scale_color_manual(values=met.brewer("Klimt", 6)[c(2,1,3,4,5,6)])+
  ylab('Canopy volume increase \n in gaps and intact forest (m3/yr)')+
  xlab('Canopy voume decrease \n from disturbance (m3/yr)')


p2 = ggplot(df3_wide , 
            aes(x=net_volume_increase/10000,y=tile_max_height_a,color = site, group = site))+
  geom_point(size=0.3,alpha=0.5)+#geom_smooth(method="lm",se=F)+
  theme_bw()+theme(legend.title = element_blank())+
  geom_vline(xintercept=0,color="black")+ylim(c(30,90))+
  scale_color_manual(values=met.brewer("Klimt", 6)[c(2,1,3,4,5,6)])+
  #ylim(c(30,100))+
  xlim(c(-2,1))+
  ylab('Initial max canopy height (m)')+
  xlab('Net height change (m/yr)')

pc = p1 + p2 +plot_layout(guides = "collect")  +plot_annotation(tag_levels = 'A') & 
  theme(legend.position = "right", legend.title=element_blank()) &
  guides(colour = guide_legend(override.aes = list(size=5)))
pc

#pc
ggsave(plot=pc, height = 3, width = 8, filename = paste0(project_folder, "Fig4.Disturbance_vs_recovery_vs_height_ha.png"))


# Fig 5.  Disturbance effects + forest #####
df3_wide$site %<>% factor(levels = rev(c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum" )))

gd1 = ggplot(df3_wide , 
             aes(tile_max_height_a,proportion_of_area_disturbance,color = site))+
  geom_point(size=0.3,alpha=0.5)+geom_smooth(method="lm",se=F,linewidth=0.4)+theme_bw()+#scale_y_log10()+
  scale_color_manual(values=met.brewer("Klimt", 6)[c(2,1,3,4,5,6)])+
  ylab('Area disturbed (%)')+xlab("Max canopy height (m)")

gd2 = ggplot(df3_wide, 
             aes(tile_elevation_mean,proportion_of_area_disturbance,color = site))+
  geom_point(size=0.3,alpha=0.5)+geom_smooth(method="lm",se=F,linewidth=0.4)+theme_bw()+
  scale_color_manual(values=met.brewer("Klimt", 6)[c(2,1,3,4,5,6)])+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  ylab('')+xlab("Elevation (m.a.s.l)")


gd3 = ggplot(df3_wide , 
             aes(proportion_of_area_gaps,proportion_of_area_disturbance,color = site))+
  geom_point(size=0.3,alpha=0.5)+geom_smooth(method="lm",se=F,linewidth=0.4)+theme_bw()+
  scale_color_manual(values=met.brewer("Klimt", 6)[c(2,1,3,4,5,6)])+
  ylab('Area disturbed (%)')+xlab("Gap fraction (%)")

gd4 = ggplot(df3_wide , 
             aes(tile_pulsedensity_increase,proportion_of_area_disturbance,color = site))+
  geom_point(size=0.3,alpha=0.5)+geom_smooth(method="lm",se=F,linewidth=0.4)+theme_bw()+
  scale_color_manual(values=met.brewer("Klimt", 6)[c(2,1,3,4,5,6)])+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  ylab('')+xlab("Pulse density increase")


gdc = gd1+gd2+gd3+gd4+plot_annotation(tag_levels = 'A')+
  plot_layout(nrow = 2, guides = "collect") & #guides(color = guide_legend(reverse=TRUE))&
  theme(legend.position = "right",legend.title = element_blank(),panel.grid.minor = element_blank()) &
  guides(colour = guide_legend(override.aes = list(size=5)))
gdc

ggsave(plot=gdc, height = 6, width = 7, filename = paste0(project_folder, "Fig5.Trends_disturbance_10ha_fixedH.png"))


# Run disturbance models ########

temp_models = list();  c=1
for (i in c("Ducke", "Tapajos",  "Nouragues", "Paracou", "Danum", "Sepilok")){

      # Filter to these classes 
      this_df3 = df3_wide %>% ungroup() %>% 
        dplyr::filter(site==i ) %>%  tidyr::drop_na()

      # If there is enough data, run the model
      if (nrow(this_df3)>10){
        
        # Scale the predictor variables
        this_df3 %<>%  mutate(proportion_of_area_gaps=scale(proportion_of_area_gaps), tile_max_height_a=scale(tile_max_height_a), 
                              tile_elevation_mean = scale(tile_elevation_mean),tile_pulsedensity_increase = scale(tile_pulsedensity_increase))
        
        # Disturbance models
        lm4 = lm(proportion_of_area_disturbance~tile_pulsedensity_increase+tile_max_height_a+tile_elevation_mean+proportion_of_area_gaps, data=this_df3);

        lm_final_disturbance=broom::tidy(lm4) %>% mutate(r2=summary(lm4)$adj.r.squared,N=length(lm4$fitted.values))
        lm_final_disturbance %<>% mutate(site = i, vif=c(0,car::vif(lm4)))

        temp_models[[c]] =  lm_final_disturbance
        c=c+1
        
    } # End of if enough data
} # end site loop


dmodels = temp_models %>% dplyr::bind_rows()
dmodels %<>% dplyr::filter(term != "(Intercept)")
dmodels$site     %<>% factor(levels = (c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum")))


# Fig 6 - model forest plot #######

# Rename terms for figure
dmodels$term2 = NA
dmodels$term2[dmodels$term == "tile_pulsedensity_increase"] = "Pulse density increase"
dmodels$term2[dmodels$term == "tile_max_height_a"] = "Max canopy height (m)"
dmodels$term2[dmodels$term == "tile_elevation_mean"] = "Elevation (m.a.s.l)"
dmodels$term2[dmodels$term == "proportion_of_area_gaps"] = "Gap fraction (%)"

ggforestplot::forestplot(
  df = dmodels,
  name=term2,
  estimate = estimate,
  se=std.error,
  colour = site,
  ci=0.95
)+ theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(name = "",values=met.brewer("Klimt", 6,direction=1)[rev(c(2,1,3,4,5,6))])+
  xlim(c(-3,3))+
  xlab("Change in area disturbed (%)")+#theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  guides(color = guide_legend(reverse=TRUE)) #+coord_flip()


ggsave(height = 4, width = 6, filename = paste0(project_folder,"Fig6.Model_effect_sizes.png"))

