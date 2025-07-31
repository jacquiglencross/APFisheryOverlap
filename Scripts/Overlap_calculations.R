
# ** Overlap with observed penguin tracks **


#load packages
pacman::p_load(tidyverse, stringr, #tidyr
               data.table, imputeTS, #datawrangling
               foreach, doParallel, furrr, progressr, #parallel
               crawl, forecast, sf, trip, moveHMM, #movement analysis
               ggplot2, gganimate, ggnewscale,  #datavis
               beepr, rstudioapi, here, lubridate, zoo)
conflicted::conflicts_prefer(lubridate::year)



# STEP 1: Read in data #####

#clear some space
rm(list = ls())
gc()
invisible(gc())

#directories
output_dir <- ("E:/Chapter 2 HMM/processed_data/")
shapefile_dir <- ("E:/Shapefiles/")
files_dir <- ("E:/Chapter 2 HMM/")
figure_dir <- ("C:/Users/Jacqui Glencross/OneDrive/Documents/PhD/Thesis/Chapter 1 Overlap/Figures/")
setwd(figure_dir)



#shapefiles
coastline <- read_sf(here(shapefile_dir,"gadm36_ZAF_3.shp"))# crs=4326  #read in shapefile of Robben and Dassen Islands

closures <- read.csv(file = here(files_dir,"Closures.csv"), header=T) 
MPA <- read_sf(here(shapefile_dir,"RobbenIslandMPA.shp")) %>%
  mutate(Type = "MPA") %>%
  dplyr::select(geometry, Type)
RI <- read_sf(here(shapefile_dir,"RI20km.shp")) %>%
  mutate(Type = "RI") %>%
  dplyr::select(geometry, Type)
DI<- read_sf(here(shapefile_dir,"DI20km.shp")) %>%
  mutate(Type = "DI") %>%
  dplyr::select(geometry, Type)


closureshp <- MPA %>%
  rbind(.,RI) %>%
  rbind(., DI) %>%
  left_join(., closures)  %>% 
  filter(Year > 2015) %>%
  filter(Year < 2022) %>%
  filter(Type != "MPA") %>%
  filter(!is.na(Closure))


closuressimple <- closures %>%
  filter(Type != "MPA") %>%
  mutate(Colony = ifelse(Type == "RI", "Robben", "Dassen")) %>%
  filter(Year < 2022) %>%
  dplyr::select(-Type)

RIpolygon <- read_sf(here(shapefile_dir,"RI20kmpolygon.shp")) %>%
  mutate(Type = "RI")
DIpolygon <- read_sf(here(shapefile_dir,"DI20kmpolygon.shp")) %>%
  mutate(Type = "DI")
closureareas <- RIpolygon %>%
  rbind(., DIpolygon) 

#data
metadata <- read.csv(here(files_dir,"Metadata_MASTER.csv"), header = T) %>% #df with metadata
  dplyr::select(deployID, Year, Nest_Lat, Nest_Lon, deployID, Colony, Closure, Year) %>%
  mutate(ID = deployID,
         id = deployID)

# bad interpolation
remove_files <- c("01_2019D", 
                  "02_2016D",
                  "03_2015R",
                  "03_2019D",
                  "04_2012R",
                  "04_2014D",
                  "04_2018R",
                  "04_2019R",
                  "05_2017R",
                  "06_2018R",
                  "07_2013R",
                  "07_2015R",
                  "08_2010R",
                  "09_2013D",
                  "11_2016D",
                  "14_2014R",
                  "14_2016R",
                  "15_2013D",
                  "16_2015R",
                  "17_2012D",
                  "18_2015D",
                  "19_2016R",
                  "20_2015R",
                  "20_2017R",
                  "21_2015D",
                  "21_2015R",
                  "25_2017R",
                  "32_2012D") 

gpsdeploys <- readRDS(here(output_dir,"GPS_HMM_complete.RDS")) %>%
  left_join(., (metadata %>% dplyr::select(ID, Nest_Lat, Nest_Lon)), by = "ID") %>% #add nest lat/lon to data
  filter(!(ID %in% remove_files)) %>%
  filter(Year < 2020) %>%
  filter(Year > 2015) %>%
  dplyr::select(ID)
gpsdeploys <- unique(gpsdeploys)
deploys <- gpsdeploys$ID

tracks <- readRDS(here(output_dir, "GPS_final.RDS")) %>%
  rename(., ID = deployID) %>%
  filter((ID %in% deploys)) %>% 
  mutate(Year = year(DateTime)) %>%
  left_join(., metadata) %>%
  st_as_sf(., coords = c("Lon", "Lat"), crs = 4326) %>%  
  mutate(id = ID,
         Week = format(DateTime, format="%W"),
         Month = month(DateTime)) %>%
  dplyr::select(id, geometry, Year, Month, Colony, Closure) %>%  
  mutate(YearMonth = paste(Year,"_",Month)) %>%
  mutate(Type = "African_penguin_track") 
rm(remove_files, gpsdeploys, deploys)


fishingvessels <- readRDS(file = here(output_dir,"GFWR_ALL.RDS")) %>%  
  mutate(Month = month(start),
         id = mmsi) %>%
  dplyr::select(id, geometry, Year, Month) %>%
  mutate(YearMonth = paste(Year,"_",Month)) %>%
  mutate(Type = "Fishing_vessel") %>%
  filter(Year < 2020) %>%
  filter(Year > 2015) %>%
  mutate(Colony = NA,
         Closure = NA) %>%
  st_as_sf(.)
africanpenguin <- readRDS(file = here(output_dir,"foragedivesHMM_all.RDS")) %>%  
  mutate(id = ID,
         Week = format(DateTime, format="%W"),
         Month = month(DateTime)) %>%
  left_join(., metadata) %>%
  dplyr::select(id, geometry, Year, Month, Colony, Closure) %>%  
  mutate(YearMonth = paste(Year,"_",Month)) %>%
  mutate(Type = "African_penguin_ARS") %>%
  filter(Year < 2020) %>%
  filter(Year > 2015)


sf::sf_use_s2 (FALSE)     
#closure area data
fishingvesselsarea <- fishingvessels %>%
  st_intersects(., closureareas, sparse = FALSE) 

fishingclosurearea <- as.data.frame(fishingvesselsarea) %>%
  mutate(intersects = ifelse((V1 == TRUE) | (V2 == TRUE) | (V3 == TRUE) | (V4 == TRUE), TRUE, FALSE)) %>%
  dplyr::select(intersects) %>%
  cbind(., fishingvessels) %>%
  filter(intersects == TRUE) %>%
  dplyr::select(id, geometry, Year, Month, Colony, Closure, YearMonth)
fishingclosurearea <- st_as_sf(fishingclosurearea)
rm(fishingvesselsarea)



africanpenguinarea <- africanpenguin %>%
  st_intersects(., closureareas, sparse = FALSE) 

penguinclosurearea <- as.data.frame(africanpenguinarea) %>%
  mutate(intersects = ifelse((V1 == TRUE) | (V2 == TRUE) | (V3 == TRUE) | (V4 == TRUE), TRUE, FALSE)) %>%
  dplyr::select(intersects) %>%
  cbind(., africanpenguin) %>%
  filter(intersects == TRUE)  %>%  
  left_join(metadata) %>%
  dplyr::select(id, geometry, Year, Month, Colony, Closure, YearMonth)
penguinclosurearea <- st_as_sf(penguinclosurearea)
rm(africanpenguinarea)



tracksarea <- tracks  %>%
  st_intersects(., closureareas, sparse = FALSE) 

tracksclosurearea <- as.data.frame(tracksarea) %>%
  mutate(intersects = ifelse((V1 == TRUE) | (V2 == TRUE) | (V3 == TRUE) | (V4 == TRUE), TRUE, FALSE)) %>%
  dplyr::select(intersects) %>%
  cbind(., tracks) %>%
  filter(intersects == TRUE) %>%
  left_join(metadata) %>%
  dplyr::select(id, geometry, Year, Month, Colony, Closure, YearMonth)
tracksclosurearea <- st_as_sf(tracksclosurearea)
rm(tracksarea)


#merge
allpoints <- rbind(fishingvessels, africanpenguin, tracks)
area_honeycomb_grid = st_make_grid(allpoints, c(0.02, 0.02), what = "polygons", square = FALSE)



# To sf and add grid ID
honeycomb_grid_sf = st_sf(area_honeycomb_grid) %>%
  mutate(grid_id = 1:length(lengths(area_honeycomb_grid)))# add grid ID


# collate points within the grid cells

# STEP 2: Monthly  Points ####
# STEP 2.1: AOI ####
#count number of points in each grid
months <- unique(allpoints$YearMonth)
templist = list()
for (i in months) {
  
  penguinsub_m <- subset(tracks, YearMonth == i)
  ARSsub_m <- subset(africanpenguin, YearMonth == i)
  fishingsub_m <- subset(fishingvessels, YearMonth == i)
  
  
  
  temp<- honeycomb_grid_sf %>%
    mutate(n_penguin = lengths(st_intersects(honeycomb_grid_sf, penguinsub_m)),
           n_dives = lengths(st_intersects(honeycomb_grid_sf, ARSsub_m)),
           n_fishing = lengths(st_intersects(honeycomb_grid_sf, fishingsub_m)))
  templist[[i]] <- temp  %>%
    mutate(YearMonth= i)
}

honeycomb_month <- data.table::rbindlist(templist) 

honeycomb_month  <- separate(honeycomb_month, YearMonth, c("Year", "Month"), sep = " _ ") %>%
  filter(Month > 2) %>%  # remove feb and sep because there are very few points
  filter(Month < 9)
#beep(2)

#STEP 2.2: Closure areas ####
#START HERE#####


colonies <- c("Robben", "Dassen")
months <- unique(allpoints$YearMonth)
templistj <- list()
for (j in colonies){
  ARScolony <- subset(penguinclosurearea, Colony == j)
  trackcolony <- subset(tracksclosurearea, Colony == j)
  #fishingcolony <- subset(fishingclosurearea, Colony == j)
  templisti = list()
  for (i in months) {
    
    
    
    penguinsub_c <- subset(trackcolony, YearMonth == i)
    ARSsub_c <- subset(ARScolony, YearMonth == i)
    fishingsub_c <- subset(fishingclosurearea, YearMonth == i)
    
    
    temp<- honeycomb_grid_sf %>%
      mutate(n_dives = lengths(st_intersects(honeycomb_grid_sf, ARSsub_c)),
             n_penguin = lengths(st_intersects(honeycomb_grid_sf, penguinsub_c)),
             n_fishing = lengths(st_intersects(honeycomb_grid_sf, fishingsub_c)))
    templisti[[i]] <- temp  %>%
      mutate(YearMonth = i)
  }
  monthsdf  <- data.table::rbindlist(templisti) 
  templistj[[j]] <- monthsdf  %>%
    mutate(Colony = j)
  
}


honeycomb_month_closures <- data.table::rbindlist(templistj) 


honeycomb_month_closures  <- separate(honeycomb_month_closures, YearMonth, c("Year", "Month"), sep = " _ ") %>%
  filter(Month > 2) %>%  # remove feb and sep because there are very few points
  filter(Month < 9) %>%
  mutate(Year = as.integer(Year)) %>%
  left_join(., closuressimple)





##### STEP 3: annual overlap ####
#### STEP 3.1: AOI ####
### STEP 3.1.1: equal ####

#remove tiles with no fishing or penguins
honeycomb_count_m = filter(honeycomb_month, n_fishing > 0 | n_dives > 0 | n_penguin > 0) %>%
  mutate(b_fishing = ifelse(n_fishing > 0, 1, 0),
         b_penguin = ifelse(n_penguin > 0, 1, 0),
         b_dives = ifelse(n_dives > 0, 1, 0)) %>% st_sf(.)
## STEP 3.1.1.1: ars ####
# spatial ####
honeycomb_countdifspatial_ars <- honeycomb_count_m %>%
  mutate(overlap = as.factor(ifelse(b_dives == 1 & b_fishing == 1, 1, 0)),
         overlap_c = ifelse(b_dives == 1 & b_fishing == 1, 1, 0))  %>% st_sf(.) #%>%   
#group_by(Year, Month, area_honeycomb_grid) %>%
#summarise(overlap = as.factor(mean(overlap_c)))
arsdifference_spatial_overlap <- honeycomb_countdifspatial_ars %>%
  filter(overlap_c == 1)
arsdifference_spatial_nooverlap <- honeycomb_countdifspatial_ars %>%
  filter(overlap_c == 0)


arsdifference_spatial_plot <- 
  ggplot() + theme_classic() +
  #geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = arsdifference_spatial_nooverlap, alpha = 0.7, fill = "#B0E2FF") +
  geom_sf(data = arsdifference_spatial_overlap, alpha = 0.7, fill = "dodgerblue1") + 
  #scale_fill_manual(labels = c("No overlap", "Overlap"),values=c("#B0E2FF", "dodgerblue1")) +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
arsdifference_spatial_plot
ggsave(arsdifference_spatial_plot, file = paste0(figure_dir,"Spatial overlap_2016-2019 ARS.png"),width = 30, height = 50, units = "cm")

# intensity ####
annualintensity_ars <- honeycomb_count_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_dives),
         Pnf = n_fishing/Nf,
         Pnp = n_dives/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_dives)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity))

annualintensity_rank_ars <- honeycomb_count_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_dives),
         Pnf = n_fishing/Nf,
         Pnp = n_dives/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         overlap1 = Pnp*Pnf,
         overlap2 = overlap1/((n_fishing+n_dives)/(Nf+Np)),
         overlap2 = ifelse(is.na(overlap2),0,overlap2))%>% st_sf(.) %>%
  filter(overlap2 > 0)  %>%
  group_by(Year, grid_id) %>%
  summarise(sumnf = sum(n_fishing),
            sumnp = sum(n_dives)) %>%
  st_drop_geometry() %>%
  mutate(F_rank = min_rank(sumnf),
         P_rank = min_rank(sumnp),
         rankdif = P_rank - F_rank,
         Direction = ifelse(rankdif == 0, "Low", ifelse(rankdif > 1, "Penguin dominated", "Fishery dominated"))) %>%
  dplyr::select(Year, grid_id, Direction)

figannualintensity_ars <- honeycomb_count_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_dives),
         Pnf = n_fishing/Nf,
         Pnp = n_dives/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_dives)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity)) %>%
  left_join(annualintensity_rank_ars) %>%
  mutate(intensitydirection = ifelse(Direction == "Fishery dominated", intensity *-1, intensity))

figannualintensity_arsP1 <- figannualintensity_ars %>% filter(Direction == "Fishery dominated")
figannualintensity_arsP2 <- figannualintensity_ars %>% filter(Direction != "Fishery dominated")
#shows difference on scale -1 to 1
ars_intensity_plot <- 
  ggplot() + theme_classic() +
  geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) +   
  geom_sf(data = figannualintensity_arsP1, aes(fill = intensitydirection*100), alpha = 0.7) + 
  scale_fill_gradient('Overlap intensity (Fishery dominated)', low = "red2", high = "#FFEC8B") +
  new_scale("fill") + 
  geom_sf(data = figannualintensity_arsP2, aes(fill = intensitydirection*100), alpha = 0.7) + 
  scale_fill_gradient('Overlap intensity (Penguin dominated)',low = "lightblue2", high = "#27408B") +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
ars_intensity_plot
ggsave(ars_intensity_plot, file = paste0(figure_dir,"Intensity 2016-2019.png"),width = 30, height = 50, units = "cm")

# table 1: overlap within the AOI, equal overlap, ARS ####


tableannualoverlap_ars <- honeycomb_countdifspatial_ars %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap = sum(overlap_c)/length(overlap)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(annualoverlap = mean(monthlyoverlap),
            "Spatial overlap (%)" = annualoverlap*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap)


annualintensity_ars_table <- honeycomb_count_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_dives),
         Pnf = n_fishing/Nf,
         Pnp = n_dives/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         overlap1 = Pnp*Pnf,
         overlap2 = overlap1/((n_fishing+n_dives)/(Nf+Np)),
         overlap2 = ifelse(is.na(overlap2),0,overlap2))%>% st_sf(.) 

tableannualintensity_ars <- annualintensity_ars_table %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap2 = sum(overlap2)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(annualoverlap2 = mean(monthlyoverlap2),
            "Overlap Intensity (%)" = annualoverlap2*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap2)




table1 <- tableannualoverlap_ars %>%
  left_join(.,tableannualintensity_ars)
write.table(table1, "Annual overlap ars.csv")


## STEP 3.1.1.2: tracks ####
# spatial ####
honeycomb_countdifspatial_track <- honeycomb_count_m %>%
  mutate(overlap = as.factor(ifelse(b_penguin == 1 & b_fishing == 1, 1, 0)),
         overlap_c = ifelse(b_penguin == 1 & b_fishing == 1, 1, 0)) %>% st_sf(.) 
trackdifference_spatial_overlap <- honeycomb_countdifspatial_track %>%
  filter(overlap_c == 1)
trackdifference_spatial_nooverlap <- honeycomb_countdifspatial_track %>%
  filter(overlap_c == 0)


trackdifference_spatial_plot <- 
  ggplot() + theme_classic() +
  #geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = trackdifference_spatial_nooverlap, alpha = 0.7, fill = "#B0E2FF") +
  geom_sf(data = trackdifference_spatial_overlap, alpha = 0.7, fill = "dodgerblue1") + 
  #scale_fill_manual(labels = c("No overlap", "Overlap"),values=c("#B0E2FF", "dodgerblue1")) +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
trackdifference_spatial_plot
ggsave(trackdifference_spatial_plot, file = ("Spatial overlap_2016-2019 Tracks.png"),width = 30, height = 50, units = "cm")

# intensity ####
annualintensity_track <- honeycomb_count_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_penguin)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity))

figannualintensity_track <- honeycomb_count_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_penguin)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity)) %>%
  left_join(annualintensity_rank_ars) %>%
  mutate(intensitydirection = ifelse(Direction == "Fishery dominated", intensity *-1, intensity))


figannualintensity_trackP1 <- figannualintensity_track %>% filter(Direction == "Fishery dominated")
figannualintensity_trackP2 <- figannualintensity_track %>% filter(Direction != "Fishery dominated")
#shows difference on scale -1 to 1
track_intensity_plot <- 
  ggplot() + theme_classic() +
  geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) +   
  geom_sf(data = figannualintensity_trackP1, aes(fill = intensitydirection*100), alpha = 0.7) + 
  scale_fill_gradient('Overlap intensity (Fishery dominated)', low = "red2", high = "#FFEC8B") +
  new_scale("fill") + 
  geom_sf(data = figannualintensity_trackP2, aes(fill = intensitydirection*100), alpha = 0.7) + 
  scale_fill_gradient('Overlap intensity (Penguin dominated)',low = "lightblue2", high = "#27408B") +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
track_intensity_plot

ggsave("Intensity 2016-2019 track.png",width = 30, height = 50, units = "cm")

# table 2: overlap within the AOI, equal overlap, track ####


tableannualoverlap_track <- honeycomb_countdifspatial_track %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap = sum(overlap_c)/length(overlap)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(annualoverlap = mean(monthlyoverlap),
            "Spatial overlap (%)" = annualoverlap*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap)

annualintensity_track <- honeycomb_count_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_penguin)/(Nf+Np)),
         intensity = ifelse(is.na(intensity),0,intensity)) %>% st_sf(.) 


tableannualintensity_track <- annualintensity_track %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap = sum(intensity)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(annualoverlap = mean(monthlyoverlap),
            "Overlap Intensity (%)" = annualoverlap*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap)


table2 <- tableannualoverlap_track %>%
  left_join(.,tableannualintensity_track)
write.table(table2, "Annual overlap track.csv")


# compare ARS and tracks and decide which to use ####
table2a <- table2 %>%
  mutate(Method = "tracks") 
tablecompare <- table1 %>%
  mutate(Method = "ARS") %>%
  rbind(., table2a)

# overlap does not differ substantially.
# continue with just tracks

### STEP 3.1.2: penguin centric ####
#remove tiles with no penguins
honeycomb_count_nofishing = filter(honeycomb_month, n_dives > 0 | n_penguin > 0) %>%
  mutate(b_fishing = ifelse(n_fishing > 0, 1, 0),
         b_dives = ifelse(n_dives > 0, 1, 0),
         b_penguin = ifelse(n_penguin > 0, 1, 0))

## STEP 3.1.2.1: ars ####
# spatial ####
PC_honeycomb_countdifspatial_ars <- honeycomb_count_nofishing %>%
  mutate(overlap = as.factor(ifelse(b_dives == 1 & b_fishing == 1, 1, 0)),
         overlap_c = ifelse(b_dives == 1 & b_fishing == 1, 1, 0)) %>% st_sf(.) 


PC_arsdifference_spatial_overlap <- PC_honeycomb_countdifspatial_ars %>%
  filter(overlap_c == 1)
PC_arsdifference_spatial_nooverlap <- PC_honeycomb_countdifspatial_ars %>%
  filter(overlap_c == 0)


PC_arsdifference_spatial_plot <- 
  ggplot() + theme_classic() +
  #geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = PC_arsdifference_spatial_nooverlap, alpha = 0.7, fill = "#B0E2FF") +
  geom_sf(data = PC_arsdifference_spatial_overlap, alpha = 0.7, fill = "dodgerblue1") + 
  #scale_fill_manual(labels = c("No overlap", "Overlap"),values=c("#B0E2FF", "dodgerblue1")) +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
PC_arsdifference_spatial_plot
ggsave(PC_arsdifference_spatial_overlap, file = ("Penguin centric spatial overlap_2016-2019 ARS.png"),width = 30, height = 50, units = "cm")

# intensity ####
PC_annualintensity_ars <- honeycomb_count_nofishing %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_dives),
         Pnf = n_fishing/Nf,
         Pnp = n_dives/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_dives)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity))

figPCannualintensity_ars <- honeycomb_count_nofishing %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_dives),
         Pnf = n_fishing/Nf,
         Pnp = n_dives/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_dives)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity)) %>%
  left_join(annualintensity_rank_ars) %>%
  mutate(intensitydirection = ifelse(Direction == "Fishery dominated", intensity *-1, intensity))

figPCannualintensity_arsP1 <- figPCannualintensity_ars %>% filter(Direction == "Fishery dominated")
figPCannualintensity_arsP2 <- figPCannualintensity_ars %>% filter(Direction != "Fishery dominated")
#shows difference on scale -1 to 1
PC_ars_intensity_plot <- 
  ggplot() + theme_classic() +
  geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) +   
  geom_sf(data = figPCannualintensity_arsP1, aes(fill = intensitydirection*100), alpha = 0.7) + 
  scale_fill_gradient('Overlap intensity (Fishery dominated)', low = "red2", high = "#FFEC8B") +
  new_scale("fill") + 
  geom_sf(data = figPCannualintensity_arsP2, aes(fill = intensitydirection*100), alpha = 0.7) + 
  scale_fill_gradient('Overlap intensity (Penguin dominated)',low = "lightblue2", high = "#27408B") +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
#shows difference on scale -1 to 1
PC_ars_intensity_plot
ggsave(PC_ars_intensity_plot, file = ("Penguin centric intensity 2016-2019.png"),width = 30, height = 50, units = "cm")

# table 3: overlap within the AOI, penguin centric overlap, ARS ####


PC_tableannualoverlap_ars <- PC_honeycomb_countdifspatial_ars %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap = sum(overlap_c)/length(overlap)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(annualoverlap = mean(monthlyoverlap),
            "Spatial overlap (%)" = annualoverlap*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap)


PC_annualintensity_ars <- honeycomb_count_nofishing %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_dives),
         Pnf = n_fishing/Nf,
         Pnp = n_dives/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         overlap1 = Pnp*Pnf,
         overlap2 = overlap1/((n_fishing+n_dives)/(Nf+Np)),
         overlap2 = ifelse(is.na(overlap2),0,overlap2))%>% st_sf(.) 

PC_tableannualintensity_ars <- PC_annualintensity_ars %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap2 = sum(overlap2)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(annualoverlap2 = mean(monthlyoverlap2),
            "Overlap Intensity (%)" = annualoverlap2*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap2)


table3 <- PC_tableannualoverlap_ars %>%
  left_join(.,PC_tableannualintensity_ars)
write.table(table3, "Penguin centric annual overlap ars.csv")

## STEP 3.1.2.2: tracks ####
# spatial ####
PChoneycomb_countdifspatial_track <- honeycomb_count_nofishing %>%
  mutate(overlap = as.factor(ifelse(b_penguin == 1 & b_fishing == 1, 1, 0)),
         overlap_c = ifelse(b_penguin == 1 & b_fishing == 1, 1, 0)) %>% st_sf(.) 
PCtrackdifference_spatial_overlap <- PChoneycomb_countdifspatial_track %>%
  filter(overlap_c == 1)
PCtrackdifference_spatial_nooverlap <- PChoneycomb_countdifspatial_track %>%
  filter(overlap_c == 0)


PCtrackdifference_spatial_plot <- 
  ggplot() + theme_classic() +
  #geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = PCtrackdifference_spatial_nooverlap, alpha = 0.7,aes(fill = 'No overlap')) +
  geom_sf(data = PCtrackdifference_spatial_overlap, alpha = 0.7, aes(fill = 'Overlap')) + 
  #scale_fill_manual(labels = c("No overlap", "Overlap"),values=c("#B0E2FF", "dodgerblue1")) +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  scale_fill_manual(name='Spatial overlap',
                    breaks=c('No overlap', 'Overlap'),
                    values=c('No overlap'='#B0E2FF', 'Overlap'='dodgerblue1')) +
  facet_wrap(~ Year, ncol = 2)
PCtrackdifference_spatial_plot
ggsave(PCtrackdifference_spatial_plot, file = ("PCSpatial overlap_2016-2019track.png"),width = 30, height = 50, units = "cm")

# intensity ####
PC_annualintensity_track <- honeycomb_count_nofishing %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_penguin)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity))

annualintensity_rank_track <- honeycomb_count_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         overlap1 = Pnp*Pnf,
         overlap2 = overlap1/((n_fishing+n_penguin)/(Nf+Np)),
         overlap2 = ifelse(is.na(overlap2),0,overlap2))%>% st_sf(.) %>%
  filter(overlap2 > 0)  %>%
  group_by(Year, grid_id) %>%
  summarise(sumnf = sum(n_fishing),
            sumnp = sum(n_penguin)) %>%
  st_drop_geometry() %>%
  mutate(F_rank = min_rank(sumnf),
         P_rank = min_rank(sumnp),
         rankdif = P_rank - F_rank,
         Direction = ifelse(rankdif == 0, "Low", ifelse(rankdif > 1, "Penguin dominated", "Fishery dominated"))) %>%
  dplyr::select(Year, grid_id, Direction)


figPCannualintensity_track <- honeycomb_count_nofishing %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_penguin)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity)) %>%
  left_join(annualintensity_rank_track) %>%
  mutate(intensitydirection = ifelse(Direction == "Fishery dominated", intensity *-1, intensity))

figPCannualintensity_trackP1 <- figPCannualintensity_track %>% filter(Direction == "Fishery dominated")
figPCannualintensity_trackP2 <- figPCannualintensity_track %>% filter(Direction != "Fishery dominated")
figPCannualnooverlap_track <- PCtrackdifference_spatial_nooverlap %>% dplyr::select(grid_id, area_honeycomb_grid) %>% 
  distinct(.)
#shows difference on scale -1 to 1
PC_track_intensity_plot <- 
  ggplot() + theme_classic() +
  geom_sf(data = figPCannualnooverlap_track, alpha = 0.7, fill = "#CFCFCF", aes(shape = "No overlap")) +
  labs(shape = NULL) +
  geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "#CFCFCF", alpha = 0.7) +   
  geom_sf(data = figPCannualintensity_trackP1, aes(fill = intensitydirection*-100), alpha = 0.7) + 
  scale_fill_gradient('Overlap intensity \n(Fishery dominated)', low = "#FFEC8B", high = "red2", limits = c(0,70)) +
  new_scale("fill") + 
  geom_sf(data = figPCannualintensity_trackP2, aes(fill = intensitydirection*100), alpha = 0.7) + 
  scale_fill_gradient('Overlap intensity \n(Penguin dominated)',low = "lightblue2", high = "#27408B", limits = c(0,100)) +
  geom_sf(data = coastline, fill = "darkgrey", color = "darkgrey") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.position = "none") + #
  #legend.key.size = unit(1, 'cm'), #change legend key size
  #legend.key.height = unit(1, 'cm'), #change legend key height
  #legend.key.width = unit(1, 'cm'), #change legend key width
  #legend.title = element_text(size=14), #change legend title font size
  #legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
#shows difference on scale -1 to 1
PC_track_intensity_plot
ggsave(PC_track_intensity_plot, file = ("Penguin centric intensity 2016-2019track.png"),width = 30, height = 50, units = "cm")

# table 2b: overlap within the AOI, penguin centric, track ####


PCtableannualoverlap_track <- PChoneycomb_countdifspatial_track %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap = sum(overlap_c)/length(overlap)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(annualoverlap = mean(monthlyoverlap),
            "Spatial overlap (%)" = annualoverlap*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap)

PCannualintensity_track <- honeycomb_count_nofishing %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_penguin)/(Nf+Np)),
         intensity = ifelse(is.na(intensity),0,intensity)) %>% st_sf(.) 


PCtableannualintensity_track <- PCannualintensity_track %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap = sum(intensity)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(annualoverlap = mean(monthlyoverlap),
            "Overlap Intensity (%)" = annualoverlap*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap)


table2b <- PCtableannualoverlap_track %>%
  left_join(.,PCtableannualintensity_track)
write.table(table2b, "PC Annual overlap track.csv")



#### STEP 3.2: Closure areas ####



### STEP 3.2.1: equal ####
#remove tiles with no fishing or penguins
closure_count_m = filter(honeycomb_month_closures, n_fishing > 0 | n_dives > 0 | n_penguin > 0) %>%
  mutate(b_fishing = ifelse(n_fishing > 0, 1, 0),
         b_penguin = ifelse(n_penguin > 0, 1, 0),
         b_dives = ifelse(n_dives > 0, 1, 0)) %>% st_sf(.)

## track
# spatial ####
closure_countdifspatial <- closure_count_m %>%
  mutate(overlap = as.factor(ifelse(b_penguin == 1 & b_fishing == 1, 1, 0)),
         overlap_c = ifelse(b_penguin == 1 & b_fishing == 1, 1, 0))  %>% st_sf(.) #%>%   
#group_by(Year, Month, area_honeycomb_grid) %>%
#summarise(overlap = as.factor(mean(overlap_c)))
closure_spatial_overlap <- closure_countdifspatial %>%
  filter(overlap_c == 1)
closure_spatial_nooverlap <- closure_countdifspatial %>%
  filter(overlap_c == 0)


closure_spatial_plot <- 
  ggplot() + theme_classic() +
  #geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = closure_spatial_nooverlap, alpha = 0.7, fill = "#B0E2FF") +
  geom_sf(data = closure_spatial_overlap, alpha = 0.7, fill = "dodgerblue1") + 
  #scale_fill_manual(labels = c("No overlap", "Overlap"),values=c("#B0E2FF", "dodgerblue1")) +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
closure_spatial_plot
ggsave(closure_spatial_plot, file = ("Spatial overlap_2016-2019 closures.png"),width = 30, height = 50, units = "cm")

# intensity ####

closureannualintensity <- closure_count_m %>%
  group_by(Month, Year, Colony) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_penguin)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity))

#shows difference on scale -1 to 1
closure_intensity_plot <- 
  ggplot() + theme_classic() +
  geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = closureannualintensity, aes(fill = intensity*100), alpha = 0.7) + 
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_fill_viridis_c('Overlap',option = "D", direction = -1) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
closure_intensity_plot
ggsave(closure_intensity_plot, file= ("Intensity 2016-2019 closures.png"),width = 30, height = 50, units = "cm")

# table 4: overlap within the closure areas, equal overlap, tracks ####


tableannualoverlap_closures <- closure_countdifspatial %>%
  group_by(Year, Month,Colony) %>%
  summarise(monthlyoverlap = sum(overlap_c)/length(overlap)) %>%
  ungroup() %>%
  group_by(Year, Colony) %>%
  summarise(annualoverlap = mean(monthlyoverlap),
            "Spatial overlap (%)" = annualoverlap*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap)


closureannualintensity <- closure_count_m %>%
  group_by(Month, Year, Colony) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         overlap1 = Pnp*Pnf,
         overlap2 = overlap1/((n_fishing+n_penguin)/(Nf+Np)),
         overlap2 = ifelse(is.na(overlap2),0,overlap2))%>% st_sf(.) 

closuretableannualintensity <- closureannualintensity %>%
  group_by(Year, Month, Colony) %>%
  summarise(monthlyoverlap2 = sum(overlap2)) %>%
  ungroup() %>%
  group_by(Year, Colony) %>%
  summarise(annualoverlap2 = mean(monthlyoverlap2),
            "Overlap Intensity (%)" = annualoverlap2*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap2)


table4 <- tableannualoverlap_closures %>%
  left_join(.,closuretableannualintensity)
write.table(table4, "Closure annual overlap tracks.csv")


### STEP 3.2.2: penguin centric ####
closure_count_nofishing = filter(honeycomb_month_closures, n_dives > 0 | n_penguin > 0) %>%
  mutate(b_fishing = ifelse(n_fishing > 0, 1, 0),
         b_dives = ifelse(n_dives > 0, 1, 0),
         b_penguin = ifelse(n_penguin > 0, 1, 0))
## tracks
# spatial ####
PC_closure_countdifspatial <- closure_count_nofishing %>%
  mutate(overlap = as.factor(ifelse(b_penguin == 1 & b_fishing == 1, 1, 0)),
         overlap_c = ifelse(b_penguin == 1 & b_fishing == 1, 1, 0)) %>% st_sf(.) 


PC_closure_spatial_overlap <- PC_closure_countdifspatial %>%
  filter(overlap_c == 1)
PC_closure_spatial_nooverlap <- PC_closure_countdifspatial %>%
  filter(overlap_c == 0)


PC_closure_spatial_plot <- 
  ggplot() + theme_classic() +
  #geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = PC_closure_spatial_nooverlap, alpha = 0.7, fill = "#B0E2FF") +
  geom_sf(data = PC_closure_spatial_overlap, alpha = 0.7, fill = "dodgerblue1") + 
  #scale_fill_manual(labels = c("No overlap", "Overlap"),values=c("#B0E2FF", "dodgerblue1")) +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
PC_closure_spatial_plot
ggsave(PC_closure_spatial_plot, file = ("Penguin centric spatial overlap_2016-2019 closures.png"),width = 30, height = 50, units = "cm")

# intensity ####
PC_annualintensity_closure <- closure_count_nofishing %>%
  group_by(Month, Year, Colony) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_penguin)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity))

#shows difference on scale -1 to 1
PC_intensity_plot_closure <- 
  ggplot() + theme_classic() +
  geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = PC_annualintensity_closure, aes(fill = intensity*100), alpha = 0.7) + 
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_fill_viridis_c('Overlap',option = "D", direction = -1) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_wrap(~ Year, ncol = 2)
PC_intensity_plot_closure
ggsave(PC_intensity_plot_closure, file = ("Penguin centric intensity 2016-2019 closure.png"),width = 30, height = 50, units = "cm")

# table 5: overlap within the closure areas, penguin centric overlap, tracks ####


PC_tableannualoverlap_closure <- PC_closure_countdifspatial %>%
  group_by(Year, Month, Colony) %>%
  summarise(monthlyoverlap = sum(overlap_c)/length(overlap)) %>%
  ungroup() %>%
  group_by(Year, Colony) %>%
  summarise(annualoverlap = mean(monthlyoverlap),
            "Spatial overlap (%)" = annualoverlap*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap)


PC_annualintensity_closure <- closure_count_nofishing %>%
  group_by(Month, Year, Colony) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         overlap1 = Pnp*Pnf,
         overlap2 = overlap1/((n_fishing+n_penguin)/(Nf+Np)),
         overlap2 = ifelse(is.na(overlap2),0,overlap2))%>% st_sf(.) 

PC_tableannualintensity_closure <- PC_annualintensity_closure %>%
  group_by(Year, Month, Colony) %>%
  summarise(monthlyoverlap2 = sum(overlap2)) %>%
  ungroup() %>%
  group_by(Year, Colony) %>%
  summarise(annualoverlap2 = mean(monthlyoverlap2),
            "Overlap Intensity (%)" = annualoverlap2*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap2)


table5 <- PC_tableannualoverlap_closure %>%
  left_join(.,PC_tableannualintensity_closure)
write.table(table5, "Penguin centric annual overlap tracks closure.csv")



##### STEP 4: monthly ####
#### STEP 4.1: AOI ####
### STEP 4.1.1: equal ####
## ars
#remove tiles with no fishing or penguins
honeycomb_count_m = filter(honeycomb_month, n_fishing > 0 | n_dives > 0 | n_penguin > 0) %>%
  mutate(b_fishing = ifelse(n_fishing > 0, 1, 0),
         b_penguin = ifelse(n_penguin > 0, 1, 0),
         b_dives = ifelse(n_dives > 0, 1, 0)) %>% st_sf(.) %>%
  filter(as.numeric(Month) > 5)
## STEP 4.1.1.1: tracks ####
# spatial ####
honeycomb_countdifspatial <- honeycomb_count_m %>%
  mutate(overlap = as.factor(ifelse(b_penguin == 1 & b_fishing == 1, 1, 0)),
         overlap_c = ifelse(b_penguin == 1 & b_fishing == 1, 1, 0))  %>% st_sf(.) #%>%   
#group_by(Year, Month, area_honeycomb_grid) %>%
#summarise(overlap = as.factor(mean(overlap_c)))
difference_spatial_overlap <- honeycomb_countdifspatial %>%
  filter(overlap_c == 1)
difference_spatial_nooverlap <- honeycomb_countdifspatial %>%
  filter(overlap_c == 0)


difference_spatial_plot_m <- 
  ggplot() + theme_classic() +
  #geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = difference_spatial_nooverlap, alpha = 0.7, fill = "#B0E2FF") +
  geom_sf(data = difference_spatial_overlap, alpha = 0.7, fill = "dodgerblue1") + 
  #scale_fill_manual(labels = c("No overlap", "Overlap"),values=c("#B0E2FF", "dodgerblue1")) +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_grid(Month ~ Year)
difference_spatial_plot_m
ggsave(difference_spatial_plot_m, file = ("Monthly spatial overlap_2016-2019 tracks.png"),width = 30, height = 50, units = "cm")

# intensity ####
monthlyintensity <- honeycomb_count_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_penguin)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity))

#shows difference on scale -1 to 1
intensity_plot_m <- 
  ggplot() + theme_classic() +
  geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = monthlyintensity, aes(fill = intensity*100), alpha = 0.7) + 
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_fill_viridis_c('Overlap',option = "D", direction = -1) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_grid(Month ~ Year)
intensity_plot_m
ggsave(intensity_plot_m, file = ("Monthly intensity 2016-2019.png"),width = 30, height = 50, units = "cm")

# table 6: overlap within the AOI, equal overlap, ARS ####


tablemonthlyoverlap <- honeycomb_countdifspatial %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap = sum(overlap_c)/length(overlap)) %>%
  ungroup() %>%
  group_by(Month) %>%
  summarise(annualoverlap = mean(monthlyoverlap),
            "Spatial overlap (%)" = annualoverlap*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap)


monthlyintensity <- honeycomb_count_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         overlap1 = Pnp*Pnf,
         overlap2 = overlap1/((n_fishing+n_penguin)/(Nf+Np)),
         overlap2 = ifelse(is.na(overlap2),0,overlap2))%>% st_sf(.) 

tablemonthlyintensity <- monthlyintensity %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap2 = sum(overlap2)) %>%
  ungroup() %>%
  group_by(Month) %>%
  summarise(annualoverlap2 = mean(monthlyoverlap2),
            "Overlap Intensity (%)" = annualoverlap2*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap2)


table6 <- tablemonthlyoverlap %>%
  left_join(.,tablemonthlyintensity)
write.table(table6, "Monthly overlap.csv")





### STEP 4.1.2: penguin centric ####

PC_count_nofishing_m = filter(honeycomb_month, n_dives > 0 | n_penguin > 0) %>%
  mutate(b_fishing = ifelse(n_fishing > 0, 1, 0),
         b_dives = ifelse(n_dives > 0, 1, 0),
         b_penguin = ifelse(n_penguin > 0, 1, 0))%>%
  filter(as.numeric(Month) > 5)
## tracks 
# spatial ####
PC_countdifspatial_m <- PC_count_nofishing_m %>%
  mutate(overlap = as.factor(ifelse(b_penguin == 1 & b_fishing == 1, 1, 0)),
         overlap_c = ifelse(b_penguin == 1 & b_fishing == 1, 1, 0)) %>% st_sf(.) %>%
  filter(as.numeric(Month) > 5)


PC_spatial_overlap_m <- PC_countdifspatial_m %>%
  filter(overlap_c == 1)
PC_spatial_nooverlap_m <- PC_countdifspatial_m %>%
  filter(overlap_c == 0)


PC_spatial_plot_m <- 
  ggplot() + theme_classic() +
  #geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) + 
  geom_sf(data = PC_spatial_nooverlap_m, alpha = 0.7, fill = "#B0E2FF") +
  geom_sf(data = PC_spatial_overlap_m, alpha = 0.7, fill = "dodgerblue1") + 
  #scale_fill_manual(labels = c("No overlap", "Overlap"),values=c("#B0E2FF", "dodgerblue1")) +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_grid(Month~ Year)
PC_spatial_plot_m
ggsave(PC_spatial_plot_m, file = ("Monthly penguin centric spatial overlap_2016-2019 tracks.png"),width = 30, height = 50, units = "cm")

# intensity ####

PC_monthlyintensity <- PC_count_nofishing_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         intensity = (Pnp*Pnf)/((n_fishing+n_penguin)/(Nf+Np)),
         intensity = ifelse(intensity==0,NA,intensity)) %>% st_sf(.) %>%
  filter(!is.na(intensity))


annualintensity_rank_ars <- PC_monthlyintensity %>%
  group_by(Year, grid_id) %>%
  summarise(sumnf = sum(n_fishing),
            sumnp = sum(n_dives)) %>%
  st_drop_geometry() %>%
  mutate(F_rank = min_rank(sumnf),
         P_rank = min_rank(sumnp),
         rankdif = P_rank - F_rank,
         Direction = ifelse(rankdif == 0, "Low", ifelse(rankdif > 1, "Penguin dominated", "Fishery dominated"))) %>%
  dplyr::select(Year, grid_id, Direction)


figmonthintensity_track <- PC_monthlyintensity %>%
  left_join(annualintensity_rank_ars) %>%
  mutate(intensitydirection = ifelse(Direction == "Fishery dominated", intensity *-1, intensity))




figmonthintensity_trackP1 <- figmonthintensity_track %>% filter(Direction == "Fishery dominated")
figmonthintensity_trackP2 <- figmonthintensity_track %>% filter(Direction != "Fishery dominated")
#shows difference on scale -1 to 1
PCtrack_intensity_plot_m <- 
  ggplot() + theme_classic() +
  geom_sf(data = area_honeycomb_grid, fill = NA, lwd = 0.05, color = "lightgrey", alpha = 0.7) +  
  geom_sf(data = PC_spatial_nooverlap_m, alpha = 0.7, fill = "#CFCFCF", aes(shape = "No overlap")) + 
  labs(shape = NULL) +
  geom_sf(data = figmonthintensity_trackP1, aes(fill = intensitydirection*100), alpha = 0.7) + 
  scale_fill_gradient('Overlap intensity \n(Fishery dominated)', low = "red2", high = "#FFEC8B") +
  new_scale("fill") + 
  geom_sf(data = figmonthintensity_trackP2, aes(fill = intensitydirection*100), alpha = 0.7) + 
  scale_fill_gradient('Overlap intensity \n(Penguin dominated)',low = "lightblue2", high = "#27408B") +
  geom_sf(data = coastline, color = "black") +
  geom_sf(data = closureshp, aes(lty = Closure), linewidth = 2) +
  scale_y_continuous(breaks = seq(-34, -33.23,  0.5), limits = c(-34, -33.23)) +
  scale_x_continuous(breaks = seq(17.85, 18.5,  0.5), limits = c(17.85, 18.5)) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size 
  facet_grid(Month ~ Year)
PCtrack_intensity_plot_m


ggsave(PCtrack_intensity_plot_m, file = ("Monthly penguin centric intensity 2016-2019 tracks.png"),width = 30, height = 50, units = "cm")

# table 5: overlap within the AOI, penguin centric overlap, tracks ####


PC_tableoverlap_m <- PC_countdifspatial_m %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap = sum(overlap_c)/length(overlap)) %>%
  ungroup() %>%
  group_by(Month, Year) %>%
  summarise(annualoverlap = mean(monthlyoverlap),
            "Spatial overlap (%)" = annualoverlap*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap)


PC_intensity_m <- PC_count_nofishing_m %>%
  group_by(Month, Year) %>%
  mutate(Nf = sum(n_fishing),
         Np = sum(n_penguin),
         Pnf = n_fishing/Nf,
         Pnp = n_penguin/Np,
         Pnf = ifelse(is.na(Pnf),0,Pnf),
         Pnp = ifelse(is.na(Pnp),0,Pnp),
         overlap1 = Pnp*Pnf,
         overlap2 = overlap1/((n_fishing+n_penguin)/(Nf+Np)),
         overlap2 = ifelse(is.na(overlap2),0,overlap2))%>% st_sf(.) 

PC_tableintensity_m <- PC_intensity_m %>%
  group_by(Year, Month) %>%
  summarise(monthlyoverlap2 = sum(overlap2)) %>%
  ungroup() %>%
  group_by(Month, Year) %>%
  summarise(annualoverlap2 = mean(monthlyoverlap2),
            "Overlap Intensity (%)" = annualoverlap2*100)%>%
  st_drop_geometry() %>%
  dplyr::select(-annualoverlap2)


table6 <- PC_tableoverlap_m %>%
  left_join(.,PC_tableintensity_m)
write.table(table6, "Monthly penguin centric overlap tracks.csv")












