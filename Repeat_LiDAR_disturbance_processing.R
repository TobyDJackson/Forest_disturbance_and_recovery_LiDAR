# This script processes canopy height rasters to give disturbance and recovery dynamics
# Written by Toby Jackson for the paper: 
# 'Tall Bornean forests experience higher disturbance rates than eastern Amazonia'

# FUNCTIONS ##########

CD_tile_class_summary = function(tile_stack){
  
  
  # This function summarizes the change in each gap and disturbance based on the rasters only 
  # It relies on the order of cells in the raster, and dplyr summarize. this is much faster than polygonizing 
  tile_tb = tile_stack %>% as.data.frame(xy=TRUE)
  
  # Define classess
  # 0 = background. 1 = Recovering. 2 = Disturbance. 
  # This is a heirarchy. The higher numbers overwrite the others
  tile_tb$class_num = 0 #"Intact canopy"
  tile_tb$class_num[tile_tb$a_gaps>=1]=1 # Recovered gap
  tile_tb$class_num[tile_tb$disturbances>=1]=2 # Disturbance (not resulting in a gap)
  tile_tb$class_num[tile_tb$b_gaps>=1]=3  # New gap
  tile_tb$class_num[tile_tb$b_gaps>=1 & tile_tb$a_gaps>=1]=4  # Persistent gap
  
  #summary(as.factor(tile_tb$class_num))
  
  tile_tb$class=NA
  tile_tb$class[tile_tb$class_num==0] = "Intact canopy"
  tile_tb$class[tile_tb$class_num==1] = "Gap - recovered"
  tile_tb$class[tile_tb$class_num==2] = "Disturbance - canopy"
  tile_tb$class[tile_tb$class_num==3] = "Disturbance - new gap"
  tile_tb$class[tile_tb$class_num==4] = "Gap - persistent"
  
  # Add classs to the rast tile_stack
  #class_rast =tile_tb %>% dplyr::select(x,y,class_num) %>%
  #  terra::rast( type="xyz", crs=terra::crs(tile_stack), digits=6, extent=NULL)
  #plot(class_rast)
  #plot(tile_stack)
  
  # Pixel summary
  df_pixels = tile_tb %>% dplyr::ungroup() %>% dplyr::group_by(class) %>%
    dplyr::summarize(area = n(), 
                     CHM1_mean = mean(a_chm,na.rm=TRUE), 
                     height_increase_mean = mean(height_increase,na.rm=TRUE))
  df_pixels$change_type = "tile_summary_by_class"
  df_pixels$gap_or_disturbance_id =0
  
  
  # combine the different pieces
  tile_summary_tb = df_pixels
  #combined_summary_and_rast = tile_summary_tb
  return(tile_summary_tb)
}


getForestGaps_terra <- function(chm_layer, threshold = 10, size = c(1, 10^4)) {
  chm_layer[chm_layer > threshold] <- NA
  chm_layer[chm_layer <= threshold] <- 1
  gaps <- terra::patches(chm_layer, directions = 8, allowGaps = FALSE)
  rcl <- terra::freq(gaps)
  rcl$layer<-NULL
  rcl[, 2] <- rcl[, 2] * terra::res(chm_layer)[1]^2
  z <- terra::classify(gaps, rcl = rcl, right = FALSE)
  z[is.na(gaps)] <- NA
  gaps[z > size[2]] <- NA
  gaps[z < size[1]] <- NA
  gaps <- terra::patches(gaps, directions = 8, allowGaps = FALSE)
  names(gaps) <- "gaps"
  return(gaps)
}


  
# Packages and paths =========
library(tidyverse); library(magrittr)
library(terra)

# CHANGE THIS PATH
project_folder= "C:/Users/av23907/Dropbox/3_ALS_summary_rasters/v2_sites_2023/Deposit/"

# Set parameters #######
tilesize = 10 # ha
height_threshold_type= "fixedH" 
min_area_threshold = 25

temp_list = list(); c=1
for (site in c("Danum","Sepilok","Tapajos","Ducke","Nouragues","Paracou")){
  

  a_chm = terra::rast(paste0(project_folder,site,"/a_chm_lspikefree_multi3.1_slope1.75_offset2.1.tif"))
  b_chm = terra::rast(paste0(project_folder,site,"/b_chm_lspikefree_multi3.1_slope1.75_offset2.1.tif"))
  height_increase = b_chm - a_chm
  a_pulsedensity  = terra::rast(paste0(project_folder,site,"/a_pulsedensity.tif"))
  b_pulsedensity  = terra::rast(paste0(project_folder,site,"/b_pulsedensity.tif"))
  d_pulsedensity = b_pulsedensity - a_pulsedensity
  dtm = terra::rast(paste0(project_folder,site,"/b_dtm.tif"))
  
  whole_stack = terra::rast(list(a_chm,b_chm,height_increase,dtm,d_pulsedensity))
  names(whole_stack) = c("a_chm","b_chm","height_increase" ,"b_dtm","d_pulsedensity" )
  #whole_stack=terra::rast(list.files(paste0(project_folder,site,"/"),full.names = T))
  #names(whole_stack)
  #names(whole_stack) = c("a_chm","b_chm","b_dtm","height_increase" )
  #terra::plot(whole_stack)
  
  # This should be a set of polygons covering the whole area - I used 10 ha squares.
  chm_agg = whole_stack$a_chm %>% terra::aggregate(round(sqrt(tilesize)*100)) 
  values(chm_agg) = seq(1,nrow(chm_agg)*ncol(chm_agg))
  tile_polygons = terra::as.polygons(chm_agg)
  names(tile_polygons) = "id"
  plot(whole_stack$a_chm)
  plot(tile_polygons,add=T)
  
  # Loop over tiles ####
  for (tile_id in seq(1,nrow(tile_polygons))){
  
    tile_stack = whole_stack %>% terra::crop(tile_polygons[tile_id])
    tile_a_chm_values  = terra::values(tile_stack$a_chm)
    tile_mean_height_a = mean(tile_a_chm_values,na.rm=T)
    tile_dtm_values    = terra::values(tile_stack$b_dtm)
    tile_pulsedensity_increase    = terra::values(tile_stack$d_pulsedensity)
    notNA = sum(is.na(tile_a_chm_values)==F)
    (percentage_notNA = notNA/(tilesize*100))
    
    # __ Define height thresholds #############
    if (height_threshold_type== "fixedH"){ gap_height_threshold = 10; disturbance_height_threshold = 5}
    if (height_threshold_type== "meanH"){ gap_height_threshold = tile_mean_height_a/2; disturbance_height_threshold = tile_mean_height_a/4}
    
    if (is.na(tile_mean_height_a)==0){
      if (percentage_notNA>90){
        
        #  __ Find gaps & disturbance #######
        tile_stack$a_gaps = tile_stack$a_chm %>% getForestGaps_terra(gap_height_threshold, size=c(min_area_threshold,1e6))
        tile_stack$b_gaps = tile_stack$b_chm %>% getForestGaps_terra(gap_height_threshold, size=c(min_area_threshold,1e6)) 
        tile_stack$disturbances = tile_stack$height_increase %>% getForestGaps_terra(-disturbance_height_threshold , size=c(min_area_threshold,1e6)) 
        
        # Log all parameters used
        log = data.frame(site=site,
                         tilesize, tile_id , n_tiles = nrow(tile_polygons), 
                         height_threshold_type , gap_height_threshold,  min_area_threshold , 
                         disturbance_height_threshold, percentage_notNA,
                         tile_pulsedensity_increase = mean(tile_pulsedensity_increase,na.rm=T),
                         tile_mean_height_a = mean(tile_a_chm_values,na.rm=T),
                         tile_max_height_a = max(tile_a_chm_values,na.rm=T),
                         tile_mean_height_increase = mean(terra::values(tile_stack$height_increase),na.rm=T),
                         tile_elevation_mean = mean(tile_dtm_values,na.rm=T),
                         tile_elevation_range = max(tile_dtm_values,na.rm=T) - min(tile_dtm_values,na.rm=T))
 
        # __ Run analysis  #####
        temp_results  = CD_tile_class_summary(tile_stack)
        temp_list[[c]] = cbind(log,temp_results) 
        c=c+1
        
      } # End if notNA statmenet
    } #  End of if mean height != na statement 
    #log = data.frame(datetime_processed = Sys.time(),site = site, tilesize = tilesize, raster_algorithm = raster_algorithm, tile_id = tile_id, height_threshold_type = height_threshold_type)
      print(log)
  } # End tile loop  
}

combined_results = dplyr::bind_rows(temp_list)
# SAVE THESE RESULTS ##########
readr::write_csv(combined_results,file=paste0(project_folder,"combined_results_030424.csv"))






