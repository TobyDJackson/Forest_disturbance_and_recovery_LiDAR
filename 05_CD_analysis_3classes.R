
library(tidyverse); library(magrittr)
library(readxl); library(MetBrewer)
library(tidytext); library(patchwork)
library(corrmorant); library(ggfortify)
library(broom.mixed); library(lme4)
library(ggforestplot); library(dplyr)

dropbox=paste0("C:/Users/av23907/Dropbox/")
dt = readr::read_csv(file = paste0(dropbox,"3_ALS_summary_rasters/summary_data/allsites_allparams_tile_summary_by_class.csv")) %>% dplyr::filter(percentage_notNA>95)
dt %<>% mutate(parameters = paste0(raster_algorithm, "_",tilesize,"ha, ",height_threshold_type,", min ", min_area_threshold,"m2, ",change_layer))

dt %<>% dplyr::filter(param == "a")
dt %<>% dplyr::filter(site != "D13", site != "D14")

dt %<>% mutate(site_tile = paste(site,tilesize))
summary(as.factor(dt$site_tile))
names(dt)
summary(as.factor(dt$param))
summary(as.factor(dt$class))

# Three classes for models ########
#dt %<>% dplyr::filter(site=="Danum",tile_id==20)
dt$class2 = "intact"
#dt %<>% dplyr::filter(class != "Disturbance - canopy")
dt$class2[dt$class =="Disturbance - canopy" | dt$class=="Disturbance - new gap"]="disturbance"
dt$class2[dt$class =="Gap - persistent" | dt$class=="Gap - recovered"]="gaps"
dt2 = dt %>% dplyr::ungroup() %>% dplyr::group_by(site, parameters,tilesize,interval,height_threshold_type,# Parameters
                                                  tile_id, percentage_notNA,tile_mean_height_a,tile_max_height_a, tile_mean_height_increase, 
                                                  tile_elevation_mean,tile_elevation_range, tile_pd_increase, # tile attributes
                                                  class2) %>%
  dplyr::summarize(CHM1_mean = sum(CHM1_mean*area,na.rm=T)/sum(area,na.rm=T), 
                   height_increase_mean = sum(height_increase_mean*area,na.rm=T)/sum(area,na.rm=T), 
                   area=sum(area,na.rm=T))


# Convert to annual change rates
dt2 %<>% dplyr::mutate(volume_increase_m3_pha_py = area*height_increase_mean/(interval*tilesize), 
                       proportion_of_area = area/(percentage_notNA*tilesize),
                       height_increase_m_py = height_increase_mean/(interval)) 

# Pivot wider
dt2 %<>% dplyr::select(!c(CHM1_mean,area,height_increase_mean))

dt_wide = dt2 %>% pivot_wider(names_from=class2,values_from=c(volume_increase_m3_pha_py, height_increase_m_py,proportion_of_area))
names(dt_wide)


# Fig 4.  Disturbance vs recovery #####
dt_wide$site %<>% factor(levels = rev(c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum" )))

dt_wide %<>% mutate(net_volume_increase = volume_increase_m3_pha_py_intact+volume_increase_m3_pha_py_gaps+volume_increase_m3_pha_py_disturbance)


p1= ggplot(dt_wide , 
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


p2 = ggplot(dt_wide , 
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
ggsave(plot=pc, height = 3, width = 8, filename = paste0(dropbox,"Rscripts/ChangeDetection/CD_figures/", 
                                                         "Fig4.Disturbance_vs_recovery_vs_height_ha.png"))


# Fig 5.  Disturbance effects + forest #####
dt_wide$site %<>% factor(levels = rev(c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum" )))

gd1 = ggplot(dt_wide %>% dplyr::filter(tilesize ==ts), 
             aes(tile_max_height_a,proportion_of_area_disturbance,color = site))+
  geom_point(size=0.3,alpha=0.5)+geom_smooth(method="lm",se=F,linewidth=0.4)+theme_bw()+#scale_y_log10()+
  xlim(c(30,85))+ylim(c(0,35))+
  scale_color_manual(values=met.brewer("Klimt", 6)[c(2,1,3,4,5,6)])+
  ylab('Area disturbed (%)')+xlab("Max canopy height (m)")

gd2 = ggplot(dt_wide %>% dplyr::filter(tilesize ==ts), 
             aes(tile_elevation_mean,proportion_of_area_disturbance,color = site))+
  geom_point(size=0.3,alpha=0.5)+geom_smooth(method="lm",se=F,linewidth=0.4)+theme_bw()+
  ylim(c(0,35))+xlim(c(0,400))+
  scale_color_manual(values=met.brewer("Klimt", 6)[c(2,1,3,4,5,6)])+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  ylab('')+xlab("Elevation (m.a.s.l)")


gd3 = ggplot(dt_wide %>% dplyr::filter(tilesize ==ts), 
             aes(proportion_of_area_gaps,proportion_of_area_disturbance,color = site))+
  geom_point(size=0.3,alpha=0.5)+geom_smooth(method="lm",se=F,linewidth=0.4)+theme_bw()+
  ylim(c(0,35))+xlim(c(0,30))+
  scale_color_manual(values=met.brewer("Klimt", 6)[c(2,1,3,4,5,6)])+
  ylab('Area disturbed (%)')+xlab("Gap fraction (%)")

gd4 = ggplot(dt_wide %>% dplyr::filter(tilesize ==ts), 
             aes(tile_pd_increase,proportion_of_area_disturbance,color = site))+
  geom_point(size=0.3,alpha=0.5)+geom_smooth(method="lm",se=F,linewidth=0.4)+theme_bw()+
  ylim(c(0,35))+xlim(c(0,70))+
  scale_color_manual(values=met.brewer("Klimt", 6)[c(2,1,3,4,5,6)])+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  ylab('')+xlab("Pulse density increase")


gdc = gd1+gd2+gd3+gd4+plot_annotation(tag_levels = 'A')+
  plot_layout(nrow = 2, guides = "collect") & #guides(color = guide_legend(reverse=TRUE))&
  theme(legend.position = "right",legend.title = element_blank(),panel.grid.minor = element_blank()) &
  guides(colour = guide_legend(override.aes = list(size=5)))
gdc

ggsave(plot=gdc, height = 6, width = 7, filename = paste0(dropbox,"Rscripts/ChangeDetection/CD_figures/Models/", 
                                                          "5.Trends_disturbance_10ha_fixedH.png"))


# Run disturbance models ########

site_names = c("Ducke", "Tapajos",  "Nouragues", "Paracou", "Danum", "Sepilok")
temp_models = list(); temp_AIC_models = list(); c=1
for (i in seq(1,length(site_names))){
  (this_site = site_names[i])
  
  # DO THIS FOR DISTURBANCE DATA AS WELL AS PIXELS DATA
  for (ts in c(10)){
    for (htt in c("fixedH")){
      
      
      # Filter to these classes 
      this_dt = dt_wide %>% ungroup() %>% 
        dplyr::filter(tilesize==ts,site==this_site ) %>%
        dplyr::select(disturbance = proportion_of_area_disturbance, 
                      gap_fraction = proportion_of_area_gaps,  
                      max_height = tile_max_height_a,
                      elevation = tile_elevation_mean, 
                      pd_increase = tile_pd_increase) %>%
        tidyr::drop_na()
      
      this_dt <- this_dt[is.finite(rowSums(this_dt)),]
      
      # If there is enough data, run the model
      if (nrow(this_dt)>10){
        
        
        if (1==2){ # Plot data before scaling
          corrmorant::corrmorant(this_dt, style = "blue_red")+
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.position = "none")+
            theme(axis.text=element_text(size=8))+
            ggtitle(paste(this_site,"ALL",ts,"ha",htt))
          ggsave(filename = paste0(dropbox,"Rscripts/ChangeDetection/CD_figures/Models/Corrmorant/Corrmorant2_ALL_",this_site,"_",ts,"_",htt,".png"),
                 height = 6, width = 6)
        } 
        
        # Scale the predictor variables
        this_dt %<>%  mutate(gap_fraction=scale(gap_fraction), max_height=scale(max_height), 
                             elevation = scale(elevation),pd_increase = scale(pd_increase))
        
        # Disturbance models
        lm1 = lm(disturbance~pd_increase, data=this_dt); summary(lm1) 
        lm2 = lm(disturbance~pd_increase+max_height, data=this_dt); summary(lm2) 
        lm3 = lm(disturbance~pd_increase+max_height+elevation, data=this_dt); summary(lm3) 
        lm4 = lm(disturbance~pd_increase+max_height+elevation+gap_fraction, data=this_dt); summary(lm4) 
        
        #plot(lm1)
        
        lm_dt_disturbance = data.frame( variable_added = c( "pd_increase","max_height","elevation" ,"gap_fraction"),
                                        aic = AIC(lm1,lm2,lm3,lm4) ,height_threshold_type=htt,site = this_site,tile_size=ts, 
                                        model = "disturbance")
        
        
        lm_final_disturbance=broom::tidy(lm4) %>% mutate(r2=summary(lm4)$adj.r.squared,N=length(lm4$fitted.values))
        lm_final_disturbance %<>% mutate(height_threshold_type=htt,site = this_site,tile_size=ts, 
                                         aic = stats::AIC(lm4), vif=c(0,car::vif(lm4)), 
                                         model = "disturbance")
    
           
        temp_models[[c]] =  lm_final_disturbance
        temp_AIC_models[[c]] = lm_dt_disturbance
        
        c=c+1
        
        if (1==2){ # Plot the models
          print(ggplot2::autoplot(lm1)+ggtitle(paste(this_site,"ALL",ts,"ha",htt)))
          readline (prompt="Press [enter] to continue") 
          #ggpubr::annotate_figure(g1, top = (paste(this_site,"ALL",ts,htt)))
          #ggpubr::annotate_figure(g1, top = ggpubr::text_grob(paste(this_site,"ALL",ts,htt)))
          #ggsave(plot = g1, filename = paste0(dropbox,"Rscripts/ChangeDetection/CD_figures/Models/Corrmorant/Residuals_ALL_",this_site,"_",ts,"_",htt ,".png"),
          #       height = 6, width = 6)
        }
        
      } # End of if enough data
      #} # End topo loop
    } # end htt loop
  } # end ts loop
} # end site loop


dmodels = temp_models %>% dplyr::bind_rows()
summary(as.factor(dmodels$term))
dmodels_AIC = temp_AIC_models %>% dplyr::bind_rows()

dmodels$r2[dmodels$r2<0]=0
dmodels %<>% mutate(r2_round = paste0("R2 = ",round(r2,digits=1)), 
                    N_r2 = paste0("N = ",N,", R2 = ",round(r2,digits=2)), 
                    site_r2 = paste0(site,", R2 = ",round(r2,digits=2)),
                    size_def = paste(tile_size, "ha,", height_threshold_type))
dmodels %<>% dplyr::filter(term != "(Intercept)")

dmodels$site     %<>% factor(levels = (c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum")))
dmodels_AIC$site %<>% factor(levels = (c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum")))

dmodels_AIC$variable_added    %<>% factor(levels =  rev(c("pd_increase", "max_height","elevation" ,"gap_fraction")))
dmodels$term                  %<>% factor(levels =  rev(c("pd_increase", "max_height","elevation" ,"gap_fraction")))



# Fig 6 - model forest plot #######
# Disturbance - area only
#dmodels2$site_r2     %<>% factor(levels = (c("Ducke, R2 = 0.13","Paracou, R2 = 0.14","Tapajos, R2 = 0.61",
#                                             "Nouragues, R2 = 0.79","Sepilok, R2 = 0.83","Danum, R2 = 0.14")))

dmodels2 =dmodels %>% dplyr::filter(tile_size==10, model == "disturbance" )


# Disturbance - canopy and area 
dmodels2$site_r2     %<>% factor(levels = (c("Ducke, R2 = 0.23","Paracou, R2 = 0.19","Tapajos, R2 = 0.35",
                                             "Nouragues, R2 = 0.53","Sepilok, R2 = 0.55","Danum, R2 = 0.06")))

dmodels2$term    %<>% factor(levels =  (c("pd_increase", "max_height","elevation" ,"gap_fraction")))
summary(as.factor(dmodels2$site_r2))

dmodels2$term2 = NA
dmodels2$term2[dmodels2$term == "pd_increase"] = "Pulse density increase"
dmodels2$term2[dmodels2$term == "max_height"] = "Max canopy height (m)"
dmodels2$term2[dmodels2$term == "elevation"] = "Elevation (m.a.s.l)"
dmodels2$term2[dmodels2$term == "gap_fraction"] = "Gap fraction (%)"

ggforestplot::forestplot(
  df = dmodels2,
  name=term2,
  estimate = estimate,
  se=std.error,
  colour = site_r2,
  ci=0.95
)+ theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(name = "",values=met.brewer("Klimt", 6,direction=1)[rev(c(2,1,3,4,5,6))])+
  xlim(c(-3,3))+
  xlab("Change in area disturbed (%)")+#theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  guides(color = guide_legend(reverse=TRUE)) #+coord_flip()
gdf


ggsave(height = 4, width = 6, filename = paste0(dropbox,"Rscripts/ChangeDetection/CD_figures/Models/", 
                                                          "6.Models_forestplot_disturbance_10ha_fixedH.png"))


# Model fit summary #####

dms = dmodels %>% dplyr::filter(tile_size==ts,  height_threshold_type==gd, 
                                term == "pd_increase")
dms$model %<>% factor(levels = rev(c("disturbance","recovery","centroid","intact" )))
dms$site  %<>% factor(levels = rev(c("Ducke","Tapajos","Paracou","Nouragues","Sepilok","Danum" )))


gdm = ggplot(dms %>% filter(model !="centroid")
             , aes(r2, model, color=site, group=site))+geom_point()+
  scale_color_manual(values=met.brewer("Klimt", 6,direction=1))+theme_bw()+
  xlab("Model R2")+  ylab("") + theme(legend.position = "none")

#clipr::write_clip(dms)

g4 = gdm+gdf+grf+gif+ plot_layout(nrow = 2, guides = "collect") & 
  theme(legend.position = "right",legend.title = element_blank())

ggsave(plot=g4, height = 6, width = 8, filename = paste0(dropbox,"Rscripts/ChangeDetection/CD_figures/Models/", 
                                           "3.Models_with_fit_stats_",ts,"ha_",gd,".png"))


# Filter by percentiles #####

a = dt_wide %>% filter(site=="Danum")
min(a$volume_increase_m3_pha_py_disturbance)

a %<>% mutate(lowp = quantile(volume_increase_m3_pha_py_disturbance, probs = c(0.1)), 
              highp = quantile(volume_increase_m3_pha_py_disturbance, probs = c(0.9)))

percentile = quantile(a$volume_increase_m3_pha_py_disturbance, probs = c(0.25, 0.75))

df_filtered <- dt_wide %>%
  group_by(site) %>%
  mutate(lowp = quantile(volume_increase_m3_pha_py_disturbance, probs = c(0.1)), 
         highp = quantile(volume_increase_m3_pha_py_disturbance, probs = c(0.9))) %>%
  filter(volume_increase_m3_pha_py_disturbance >= lowp & volume_increase_m3_pha_py_disturbance <= highp)



# Cumulative disturbances volume change ######
# How much of the volume change is in small / big gaps 
# NEED THE DISTURBANCE DATA

dt = readr::read_csv(file = paste0(dropbox,"3_ALS_summary_rasters/summary_data/allsites_allparams_gaps.csv")) %>% dplyr::filter(percentage_notNA>95)
dt %<>% dplyr::filter(param == "a")
dt %<>% dplyr::filter(site != "D13", site != "D14")

dt$site_size = NA
dt$site_size[dt$site=="Danum"] = 1509
dt$site_size[dt$site=="Sepilok"] = 1779
dt$site_size[dt$site=="Nouragues"] = 2172
dt$site_size[dt$site=="Tapajos"] = 879
dt$site_size[dt$site=="Paracou"] = 657
dt$site_size[dt$site=="Ducke"] = 1088

dt %<>% dplyr::mutate(volume_increase = area*height_increase_gap) 


# SLOW
# NEED TO GET THESE ALL ON ONE PLOT!
gdcv = dt %>% group_by(site, interval) %>% 
  dplyr::select(gap_or_disturbance_id,area,site,site_size,volume_increase) %>%
  #mutate(dV = 100*dV/sum(dV,rm.na=TRUE)) %>%
  arrange(site,area,volume_increase) %>%
  mutate(cdV=cumsum(abs(volume_increase)), 
         cA = cumsum(area))


gdcv %<>% mutate(cdV_pApy = cdV/(interval*site_size), 
                 cA_pApy = 7*cA/(interval*site_size*100))
gdcv$site %<>% factor(levels = rev(c("Ducke","Paracou","Tapajos","Nouragues","Sepilok","Danum" )))

g1 = ggplot(gdcv,aes(x=area,y=cA_pApy, color=site))+
  geom_point(size = 0.2)+ xlim(c(0,10000))+ theme_bw()+
  ylab(paste0("Cumulative area disturbed (%)"))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  xlab("")+geom_vline(xintercept=1000, color="black")

g2 = ggplot(gdcv,aes(x=area,y=cA_pApy, color=site))+
  geom_point(size = 0.2)+ xlim(c(0,1500))+theme_bw()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  ylab("")+ xlab("")+geom_vline(xintercept=1000, color="black")

g3 = ggplot(gdcv,aes(x=area,y=cdV_pApy, color=site))+
  geom_point(size=0.2)+ xlim(c(0,10000))+  theme_bw()+
  ylab(paste0("Cumulative volume loss (m3)"))+
  xlab(paste0("Disturbance size (m2)"))+geom_vline(xintercept=1000, color="black")

g4 = ggplot(gdcv,aes(x=area,y=cdV_pApy, color=site))+
  geom_point(size=0.2)+   xlim(c(0,1500))+  theme_bw()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  ylab("")+xlab(paste0("Disturbance size (m2)"))+geom_vline(xintercept=1000, color="black")

gc = (g1+g2+g3+g4) + plot_annotation(tag_levels = 'A')+plot_layout(nrow = 2, guides = "collect") & 
  theme(legend.position = "right",legend.title = element_blank(),panel.grid.minor = element_blank())& 
  scale_color_manual(name = "",values=met.brewer("Klimt", 6,direction=1)[(c(2,1,3,4,5,6))])& 
  guides(colour = guide_legend(override.aes = list(size=5)))

#gc

ggsave(plot=gc, filename = paste0(dropbox,"Rscripts/ChangeDetection/CD_figures/", 
                                        "S6.Cumulative_disturbance.png"),
       height=6,width=8)
