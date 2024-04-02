# Code started on 20th April 2022 to tidy up the Change Detection (CD) project
# Major update started on 17th August 2023 (!!!) to use Fabian's rasters + look at gap re-growth
# It loops over;
# (1) Loop over the sites; 
#   (2) Loop over tiles (defined by tilesize); 
#       Clip site level rast_stack to 'tile_stack'
#     (3) Loop over gap definition type
#         Find gaps in CHM1 and CHM2, buffer, add ID to tile_stack
#       (4) Loop over change layer (dCHM, dDSM, dDSR) 
#           Find gaps in change_layer, buffer, add ID to tile_stack (these are disturbances)
#           Convert tile_stack to tile_tb (tibble)
#           Classify class between gaps in CHM1, CHM2 and change_layer
#           Summarize heights and changes by overlap for whole tile (base_layer == summary)
#           Summarize heights and changes by gaps in three 'base_layers' 
#           (a) CHM1, (b) CHM2 (c) current change layer (maybe 'gap_definition_layers' is better)
#           Combine, save, next

# The gaps in base_layer == CHM1, 
# with overlap defined using change_layer==dDSM, 
# had changes in dCHM, dDSM and dDSX

# To focus on disturbance always base_layer == change_layer
# change_layer=dCHM, base_layer = CHM1 gives 
#           
# For each of the above, it finds the gaps and the disturbances and the amount of height change in them
# It also determines whether they overlap
# Optionally, it saves the tiled outputs as rasters and polygons


  
# Packages and paths =========
library(tidyverse); library(magrittr)
library(terra)
dropbox   = "C:/Users/av23907/Dropbox/"
temp_files= "C:/Users/av23907/Desktop/temp_files/CD_temp_2023/CD_summary_files/"




params = data.frame(
  param_id ="a",
  tilesize = 10, 
  raster_algorithm="lspikefree_multi3.1_slope1.75_offset2.1", 
  height_threshold_type= "fixedH", 
  min_area_threshold = 25,
  change_layer = "chm")
site="Sepilok"
tile_id=713
param_id=2
params = make_parameter_list(a)
i=6
# 1. Loop over params #######
for (i in seq(1,nrow(params))){  
  
  # Set each parameter as a variable for easy use later.
  tilesize = params$tilesize[i] 
  raster_algorithm = params$raster_algorithm[i]
  height_threshold_type = params$height_threshold_type[i]
  min_area_threshold = params$min_area_threshold[i]
  change_layer = params$change_layer[i]
  
  raster_algorithm_name = raster_algorithm
  if (raster_algorithm =="lspikefree_multi3.1_slope1.75_offset2.1"){raster_algorithm_name = "lspikefree"}
  
  # _ 2. Loop over sites #####
  #for (site in c("D13","D14")){
  for (site in c("D13","D14","Tapajos","Ducke","Paracou","Nouragues","Sepilok", "Danum")){
    if (site =="Ducke")    {interval=5.148}
    if (site =="Tapajos")  {interval=4.595}
    if (site =="Paracou")  {interval=4.164}
    if (site =="Nouragues"){interval=4.164}
    if (site =="Sepilok")  {interval=5.315}
    if (site =="Danum")    {interval=5.362}
    if (site =="D13")      {interval=6.258}
    if (site =="D14")      {interval=5.362}
    
    site_1m_folder = paste0(dropbox,"3_ALS_summary_rasters/v2_sites_2023/",site,'/1m/')
    tile_polygons = terra::vect(paste0(dropbox,"3_ALS_summary_rasters/v2_sites_2023/",site,'/aggregated/',tilesize,'/id.shp'))
    
    whole_b_dtm     = terra::rast(paste0(site_1m_folder,"/b_dtm.tif"))
    whole_a_chm     = terra::rast(paste0(site_1m_folder,"/a_chm_",raster_algorithm,".tif"))
    whole_b_chm     = terra::rast(paste0(site_1m_folder,"/b_chm_",raster_algorithm,".tif"))
    whole_height_increase     = terra::rast(paste0(site_1m_folder,"/d_",change_layer,"_",raster_algorithm,".tif"))
    whole_d_pd     = terra::rast(paste0(site_1m_folder,"/d_pulsedensity.tif"))
    whole_stack = terra::rast(list(whole_b_dtm,whole_height_increase,whole_a_chm,whole_b_chm, whole_d_pd))
    names(whole_stack) = c("b_dtm","height_increase","a_chm", "b_chm","d_pd" )
    
    # __ 3. Loop over tiles ####
    for (tile_id in seq(1,nrow(tile_polygons))){
      
       
      tile_stack = whole_stack %>% terra::crop(tile_polygons[tile_id])
      #plot(tile_stack)
      
      tile_a_chm_values = values(tile_stack$a_chm)
      tile_mean_height_a = mean(tile_a_chm_values,na.rm=T)
      tile_dtm_values = values(tile_stack$b_dtm)
      notNA = sum(is.na(tile_a_chm_values)==F)
      (percentage_notNA = notNA/(tilesize*100))
      if (height_threshold_type== "fixedH"){ gap_height_threshold = 10; disturbance_height_threshold = 5}
      if (height_threshold_type== "meanH"){ gap_height_threshold = tile_mean_height_a/2; disturbance_height_threshold = tile_mean_height_a/4}
      
      if (is.na(tile_mean_height_a)==0){
        if (percentage_notNA>90){
       
   
          #  ___ 4. Find gaps & disturbance #######
          tile_stack$a_gaps = tile_stack$a_chm %>% getForestGaps1(gap_height_threshold, size=c(min_area_threshold,1e6))
          tile_stack$b_gaps = tile_stack$b_chm %>% getForestGaps1(gap_height_threshold, size=c(min_area_threshold,1e6)) 
          tile_stack$disturbances = tile_stack$height_increase %>% getForestGaps1(-disturbance_height_threshold , size=c(min_area_threshold,1e6)) 
  
          # Log all parameters used
          log = data.frame(datetime_processed = Sys.time(),param = params$param_id[i], site , interval, tilesize, tile_id , n_tiles = nrow(tile_polygons), 
                           height_threshold_type , gap_height_threshold,  min_area_threshold , disturbance_height_threshold,
                           raster_algorithm = raster_algorithm_name,  change_layer, percentage_notNA,
                           tile_mean_height_a = mean(tile_a_chm_values,na.rm=T),
                           tile_max_height_a = max(tile_a_chm_values,na.rm=T),
                           tile_mean_height_increase = mean(values(tile_stack$height_increase),na.rm=T),
                           tile_elevation_mean = mean(tile_dtm_values,na.rm=T),
                           tile_elevation_range = max(tile_dtm_values,na.rm=T) - min(tile_dtm_values,na.rm=T),
                           tile_pd_increase = mean(values(tile_stack$d_pd),na.rm=T))
          
        
        
          # ____ 5. MAIN - run & save  #####
          tile_summary_output      = CD_tile_class_summary(tile_stack)
          class_rast = tile_summary_output[[2]]
          tile_summary_output = cbind(log,tile_summary_output[[1]]) 
          
          gaps_output = CD_gap_recovery(tile_stack)
          gaps_combined = cbind(log,gaps_output)
      
          
          # Write class rast 
           #output_folder_path = paste0(tilesize,"ha/",raster_algorithm_name,"/",height_threshold_type, "/",change_layer,"/",site,"/")
          #folder_name = paste0(temp_files,"class_rasters/",output_folder_path);  dir.create(folder_name,recursive=T)
          #if (nrow(temp_tile_summary)>=1){ terra::writeRaster(class_rast, paste0(folder_name,"tile_summary_by_class_",append_summary_files_name,".tif"),
          #                                              names = paste0("class_",append_summary_files_name),overwrite = T)
          #}
          
          append_summary_files_name = paste0(site,"_param_",params$param_id[i],"_",tilesize,"ha_tile",tile_id,"_",raster_algorithm_name,"_",height_threshold_type,"_",change_layer,"_",min_area_threshold)
          
          # Tile summary
          temp_tile_summary = tile_summary_output %>% dplyr::filter(change_type=="tile_summary_by_class")
          #folder_name = paste0(temp_files,"tile_summary_by_class/");  dir.create(folder_name,recursive=T)
          if (nrow(temp_tile_summary)>=1){ readr::write_csv(temp_tile_summary, file = paste0(temp_files,"tile_summary_by_class/tile_summary_by_class_",append_summary_files_name,".csv"))}
          
          # tile_summary_by_class by height
          temp_tile_summary_by_height = tile_summary_output %>% dplyr::filter(change_type=="tile_summary_by_class_by_height")
          #folder_name = paste0(temp_files,"tile_summary_by_class_by_height/",output_folder_path);    dir.create(folder_name,recursive=T)
          if (nrow(temp_tile_summary_by_height)>=1){ readr::write_csv(temp_tile_summary_by_height, file = paste0(temp_files,"tile_summary_by_class_by_height/tile_summary_by_class_by_height_",append_summary_files_name,".csv"))}
          

          # gaps
          #folder_name = paste0(temp_files,"gaps/",output_folder_path);  dir.create(folder_name,recursive=T)
          if (nrow(gaps_combined)>=1){ readr::write_csv(gaps_combined, file = paste0(temp_files,"gaps/gaps_",append_summary_files_name,".csv"))}
          
          # Save logfile
          #readr::write_csv(log, file = paste0(dropbox,"3_ALS_summary_rasters/logfiles/CD_processing_log",Sys.Date(),"_",site,"_",tilesize,"ha.csv"),append=T)
          #print(log)
        
        } # End if notNA statmenet
      } #  End of if mean height != na statement 
      #log = data.frame(datetime_processed = Sys.time(),site = site, tilesize = tilesize, raster_algorithm = raster_algorithm, tile_id = tile_id, height_threshold_type = height_threshold_type)
      print(log)
    } # End tile loop  PARALLEL 
  } # End sites loop
} # End loop over params



# Combine summary files ####

# Single params

tilesize = 10
raster_algorithm="lspikefree_multi3.1_slope1.75_offset2.1"
height_threshold_type= "fixedH"
min_area_threshold = 25
change_layer = "dsm"
for (i in c("tile_summary_by_class", "tile_summary_by_class_by_height", "gaps")){
  
  gd_file_list = list.files(paste0(temp_files,i,"/",tilesize,"ha/",raster_algorithm_name,"/",
                                   height_threshold_type,"/",change_layer,"/"),pattern = ".csv",full.names=T,recursive=T,include.dirs=F)
  gd = dplyr::bind_rows(lapply(gd_file_list, read.csv, header=TRUE))
  readr::write_csv(gd, file = paste0(dropbox,"3_ALS_summary_rasters/summary_data/allsites_",tilesize,"ha_",
                                     raster_algorithm_name,"_",height_threshold_type,"_",change_layer,
                                     "_",i,".csv"))
  print(i)
}

# All params
for (i in c("tile_summary_by_class", "tile_summary_by_class_by_height", "gaps")){
#for (i in c("tile_summary_by_class")){
    
  gd_file_list = list.files(paste0(temp_files,i,"/"),pattern = ".csv",full.names=T,recursive=T,include.dirs=F)
  gd = dplyr::bind_rows(lapply(gd_file_list, read.csv, header=TRUE))
  readr::write_csv(gd, file = paste0(dropbox,"3_ALS_summary_rasters/summary_data/allsites_allparams_",i,".csv"))
  print(i)
}


# OLD - Loop over files to avoid errors
gd_file_list = list.files(paste0(temp_files,"10ha/disturbances/"),pattern = ".csv",full.names=T,recursive=T,include.dirs=F)
c=1; temp_combined=list()
for (i in seq(1,length(gd_file_list))){
  temp = readr::read_csv(gd_file_list[i], show_col_types = FALSE)
  if (nrow(temp>=1)){
    temp_combined[[c]] = temp
    c=c+1
  }
  print(round(100*i/length(gd_file_list)))
}
gd = dplyr::bind_rows(temp_combined)
gd %<>% dplyr::filter(raster_algorithm == "lspikefree_multi3.1_slope1.75_offset2.1", change_layer=="dsm")
gd %<>% dplyr::mutate(site_res = paste0(site, "_",tilesize,"ha"))
summary(as.factor(gd$site_res))
readr::write_csv(gd, file = paste0(dropbox,"3_ALS_summary_rasters/summary_data/All_sites_disturbances_lspikefree_dsm_10ha_all_flips.csv"))



CD_tile_class_summary = function(tile_stack){

  
  # This function summarizes the change in each gap and disturbance based on the rasters only 
  # It relies on the order of cells in the raster, and dplyr summarize. this is much faster than polygonizing 
  tile_tb = tile_stack %>% as.data.frame(xy=TRUE)
  
  # Define classess
  # 0 = background. 1 = Recovering. 2 = Disturbance. 
  # This is a heirarchy. The higher numbers overwrite the others
  tile_tb$class_num = 0 #"Intact canopy"
  tile_tb$class_num[tile_tb$a_gaps>=1]=1 # Recovered gap
  tile_tb$class_num[tile_tb$disturbances>=1]=2 # Disturbance (not resultig in a gap)
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
  class_rast =tile_tb %>% dplyr::select(x,y,class_num) %>%
    terra::rast( type="xyz", crs=crs(tile_stack), digits=6, extent=NULL)
  #plot(class_rast)
  #plot(tile_stack)
  
  # Pixel summary
  df_pixels = tile_tb %>% dplyr::ungroup() %>% dplyr::group_by(class) %>%
    dplyr::summarize(area = n(), 
                     CHM1_mean = mean(a_chm,na.rm=TRUE), 
                     height_increase_mean = mean(height_increase,na.rm=TRUE))
  df_pixels$change_type = "tile_summary_by_class"
  df_pixels$gap_or_disturbance_id =0
  
  
  # Pixel summary with height classes
  tile_tb_by_height = tile_tb %>% dplyr::mutate(a_chm_round = round(a_chm))
  df_pixels_by_height = tile_tb_by_height %>% dplyr::ungroup() %>% 
    dplyr::group_by(a_chm_round,class) %>%
    dplyr::summarize(area = n(), 
                     CHM1_mean = mean(a_chm,na.rm=TRUE), 
                     height_increase_mean = mean(height_increase,na.rm=TRUE))
  df_pixels_by_height$change_type = "tile_summary_by_class_by_height"
  df_pixels_by_height$gap_or_disturbance_id =0
  
  
  # combine the different pieces
  tile_summary_tb = dplyr::bind_rows(df_pixels, df_pixels_by_height)
  combined_summary_and_rast = list(tile_summary_tb,class_rast)
  return(combined_summary_and_rast)
}



CD_gap_recovery = function(tile_stack){
 
  
  # find gap centroids to measure height recovery
  a_gaps_poly = terra::as.polygons(tile_stack$a_gaps)
  if (nrow(a_gaps_poly)>=1){
  a_gaps_centroids = terra::centroids(a_gaps_poly,inside=T)
  df_centroids = terra::extract(tile_stack,a_gaps_centroids)
  df_centroids %<>% dplyr::select(ID, a_chm_centroid= a_chm, height_increase_centroid = height_increase)
  
  
  df_gaps = terra::extract(tile_stack,a_gaps_poly)
  dfs_gaps = df_gaps %>% dplyr::group_by(ID) %>% 
    dplyr::summarize(area = n(), 
                     a_chm_gap = mean(a_chm,na.rm=T), 
                     height_increase_gap = mean(height_increase,na.rm=T))
  
  
  dfs = dplyr::full_join(dfs_gaps, df_centroids, by="ID")
  }
  if (nrow(a_gaps_poly)==0){dfs= data.frame(ID = NA, area=NA, a_chm_gap= NA, height_increase_gap = NA)}
  dfs$change_type = "gaps"
 
  
  disturbances_poly = terra::as.polygons(tile_stack$disturbances)
  if (nrow(disturbances_poly)>=1){
  df_disturbances = terra::extract(tile_stack,disturbances_poly)
  dfs_disturbances = df_disturbances %>% dplyr::group_by(ID) %>% 
    dplyr::summarize(area = n(), 
                     a_chm_gap = mean(a_chm,na.rm=T), 
                     height_increase_gap = mean(height_increase,na.rm=T))
  }
  if (nrow(disturbances_poly)==0){dfs_disturbances= data.frame(ID = NA, area=NA, a_chm_gap= NA, height_increase_gap = NA)}
  dfs_disturbances$change_type = "disturbances"
  
  #ggplot(dfs,aes(height_increase_gap, height_increase_centroid))+geom_point()
  
  dfs2 = dplyr::bind_rows(dfs,dfs_disturbances)
  dfs2 %<>% dplyr::rename(gap_or_disturbance_id = ID) 
  
  if (nrow(dfs2)==0){dfs2[1,]=NA}
  return(dfs2)
}

getForestGaps1 <- function(chm_layer, threshold = 10, size = c(1, 10^4)) {
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





make_parameter_list = function(a){
  parama = data.frame(param_id="a",  tilesize = 10,  raster_algorithm="lspikefree_multi3.1_slope1.75_offset2.1", 
                      height_threshold_type= "fixedH",   min_area_threshold = 25,  change_layer = "chm")
  paramb = data.frame(param_id="b",  tilesize = 1,  raster_algorithm="lspikefree_multi3.1_slope1.75_offset2.1", 
                      height_threshold_type= "fixedH",   min_area_threshold = 25,  change_layer = "chm")
  paramc = data.frame(param_id="c",  tilesize = 10,  raster_algorithm="highest", 
                      height_threshold_type= "fixedH",   min_area_threshold = 25,  change_layer = "chm")
  paramd = data.frame(param_id="d",  tilesize = 10,  raster_algorithm="tin", 
                      height_threshold_type= "fixedH",   min_area_threshold = 25,  change_layer = "chm")
  parame = data.frame(param_id="e",  tilesize = 10,  raster_algorithm="lspikefree_multi3.1_slope1.75_offset2.1", 
                      height_threshold_type= "meanH",   min_area_threshold = 25,  change_layer = "chm")
  paramf = data.frame(param_id="f",  tilesize = 10,  raster_algorithm="lspikefree_multi3.1_slope1.75_offset2.1", 
                      height_threshold_type= "fixedH",   min_area_threshold = 10,  change_layer = "chm")
  paramg = data.frame(param_id="g",  tilesize = 10,  raster_algorithm="lspikefree_multi3.1_slope1.75_offset2.1", 
                      height_threshold_type= "fixedH",   min_area_threshold = 25,  change_layer = "dsm")
  
  params = dplyr::bind_rows(parama,paramb,paramc,paramd,parame,paramf,paramg)
  return(params)
}

