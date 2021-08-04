#DWR Chl data cleaning and joining to zoop data
#EMP Zooplankton Enviornmental Scientist
#arthur.barros@wildlife.ca.gov
#last updated 2/1/2021
rm( list = ls()) #clear env

#Load packages:
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(lubridate)
library(measurements)

#####################################
##Download the zoop data from EDI
#####################################
#https://portal.edirepository.org/nis/mapbrowse?packageid=edi.522.7

#CB_matrix
CB_url<-"https://portal.edirepository.org/nis/dataviewer?packageid=edi.522.7&entityid=c0916b64396edab85b07038e32ff0342"
#desination file
CB_dest<-"Data/CB_matrix.csv"
download.file(CB_url,CB_dest)
CB_matrix<-read.csv("Data/CB_matrix.csv",fileEncoding="UTF-8-BOM") #file incoding to stop import from adding junk characters to first column name

#Pump_matrix
Pump_url<-"https://portal.edirepository.org/nis/dataviewer?packageid=edi.522.7&entityid=0f7ffacf41372643865af053c0b07663"
#desination file
Pump_dest<-"Data/Pump_matrix.csv"
download.file(Pump_url,Pump_dest)
Pump_matrix<-read.csv("Data/Pump_matrix.csv",fileEncoding="UTF-8-BOM") #file incoding to stop import from adding junk characters to first column name

#Mysid_matrix
Mysid_url<-"https://portal.edirepository.org/nis/dataviewer?packageid=edi.522.7&entityid=0080191932b0987243936eff1bb54ee8"
#desination file
Mysid_dest<-"Data/Mysid_matrix.csv"
download.file(Mysid_url,Mysid_dest)
Mysid_matrix<-read.csv("Data/Mysid_matrix.csv",fileEncoding="UTF-8-BOM") #file incoding to stop import from adding junk characters to first column name

#Zoocode lookup table
Zoocode_url<-"https://portal.edirepository.org/nis/dataviewer?packageid=edi.522.7&entityid=a3be37d31cf146c51bc583632dbfcf06"
#desination file
Zoocode_dest<-"Data/Zoocode_lookup.csv"
download.file(Zoocode_url,Zoocode_dest)
Zoocode_lookup<-read.csv("Data/Zoocode_lookup.csv",fileEncoding="UTF-8-BOM") #file incoding to stop import from adding junk characters to first column name

#Station lookup table
Station_url<-"https://portal.edirepository.org/nis/dataviewer?packageid=edi.522.7&entityid=71dd301f30a2bc2e40f5da573dde9f97"
#desination file
Station_dest<-"Data/Station_lookup.csv"
download.file(Station_url,Station_dest)
Station_lookup<-read.csv("Data/Station_lookup.csv")

#Format dates
CB_matrix$SampleDate<-as.Date(CB_matrix$SampleDate,"%m/%d/%Y")
Mysid_matrix$SampleDate<-as.Date(Mysid_matrix$SampleDate,"%m/%d/%Y")
Pump_matrix$SampleDate<-as.Date(Pump_matrix$SampleDate,"%m/%d/%Y")

#remove "adult" from Zooplankton common names
Zoocode_lookup$CommonName<-str_remove(Zoocode_lookup$CommonName," adult")

#Join station info to catch matrices
Station_lookup<-dplyr::select(Station_lookup,-Core,-lat_degrees,-lat_minutes,-lat_seconds,-lon_degrees,-lon_minutes,-lon_seconds,-year_start,-year_end) #remove superfluous columns
CB_matrix<-CB_matrix%>%
  inner_join(Station_lookup)
Mysid_matrix<-Mysid_matrix%>%
  inner_join(Station_lookup)
Pump_matrix<-Pump_matrix%>%
  inner_join(Station_lookup)

#Use pivot_longer to gather catch matrices CPUE data from wide to long format
CB_matrix_long<-pivot_longer(CB_matrix,cols=ACARTELA:CRABZOEA,names_to = "ZooCode",values_to = "CPUE")
Mysid_matrix_long<-pivot_longer(Mysid_matrix,cols=A_aspera:Unidentified_mysid, names_to="ZooCode", values_to="CPUE")
Pump_matrix_long<-pivot_longer(Pump_matrix,cols=LIMNOSINE:BARNNAUP, names_to="ZooCode", values_to="CPUE")

#join in taxonomic information
CB_matrix_long<-CB_matrix_long%>%
  inner_join(Zoocode_lookup)
Mysid_matrix_long<-Mysid_matrix_long%>%
  inner_join(Zoocode_lookup)
Pump_matrix_long<-Pump_matrix_long%>%
  inner_join(Zoocode_lookup)

#select necessary columns
CB_clean<-dplyr::select(CB_matrix_long,SurveyCode,Core,Current,Year,Survey,SurveyRep,SampleDate,StationNZ,DWRStationNo,
                 Secchi,Chl_a,Temperature,ECSurfacePreTow,ECBottomPreTow,
                 lat,lon,km,ZooCode,CPUE,NativeOrIntro,Phylum,Class,Order,Family,CommonName)
Mysid_clean<-dplyr::select(Mysid_matrix_long,SurveyCode,Core,Current,Year,Survey,SurveyRep,SampleDate,StationNZ,DWRStationNo,
                    Secchi,Chl_a,Temperature,ECSurfacePreTow,ECBottomPreTow,
                    lat,lon,km,ZooCode,CPUE,NativeOrIntro,Phylum,Class,Order,Family,CommonName)
Pump_clean<-dplyr::select(Pump_matrix_long,SurveyCode,Core,Current,Year,Survey,SurveyRep,SampleDate,StationNZ,DWRStationNo,
                   Secchi,Chl_a,Temperature,ECSurfacePreTow,ECBottomPreTow,
                   lat,lon,km,ZooCode,CPUE,NativeOrIntro,Phylum,Class,Order,Family,CommonName)

#target specific taxa orders for specific gear types based on gear efficiency
CB_orders<-c("Calanoida","Cladocera")
Pump_orders<-c("Cyclopoida","Rotifer")
Mysid_orders<-c("Mysida")

#filter by target parameters
CB_targets<-CB_clean%>%
  filter(Order%in%CB_orders)
Pump_targets<-Pump_clean%>%
  filter(Order%in%Pump_orders)
Mysid_targets<-Mysid_clean%>%
  filter(Order%in%Mysid_orders)

#add gear type
CB_targets$Gear="CB"
Pump_targets$Gear="Pump"
Mysid_targets$Gear="Mysid"

#Join all catch matrices into one
Zoop_all<-CB_targets%>%
  rbind(Pump_targets)%>%
  rbind(Mysid_targets)


#####################################
##Load and clean the DWR chl data
#####################################

#DWR data downloaded from: https://emp.baydeltalive.com/projects/11282
#8/3/2021

dwr_stations<-read.csv("data/EMP_Discrete_Water_Quality_Stations_1975-2020.csv")
chl_1975_2016<-read_excel("data/Phytoplankton_Algal_Type_Data_1975_-2016.xlsx")
chl_2017<-read_excel("data/EMP_Phytoplankton_2017_2018_Data.xlsx", 
                     sheet = "2017 Algal Type Data")
chl_2018<-read_excel("data/EMP_Phytoplankton_2017_2018_Data.xlsx", 
                     sheet = "2018 Algal Type Data")
chl_2019<-read_excel("data/EMP_Phytoplankton_2019_Data.xlsx", 
                     sheet = "2019 Data by Algal Type")

#Format dates
chl_1975_2016$Date<-as.Date(chl_1975_2016$Date,"%m/%d/%Y")
chl_2017$Date<-as.Date(chl_2017$Date,"%m/%d/%Y")
chl_2018$Date<-as.Date(chl_2018$Date,"%m/%d/%Y")
chl_2019$Date<-as.Date(chl_2019$Date,"%m/%d/%Y")

#select columns
chl_2017<-dplyr::select(chl_2017,Date,"Station Code","Pennate Diatom","Centric Diatom","Chrysophyte","Cyanobacterium","Cryptophyte","Dinoflagellate","Euglenoid","Green Alga")
chl_2018<-dplyr::select(chl_2018,Date,"Station Code","Pennate Diatom","Centric Diatom","Chrysophyte","Cyanobacterium","Cryptophyte","Dinoflagellate","Euglenoid","Green Alga")
chl_2019<-dplyr::select(chl_2019,Date,"Station Code","Pennate Diatom","Centric Diatom","Chrysophyte","Cyanobacterium","Cryptophyte","Dinoflagellate","Euglenoid","Green Alga")
chl_1975_2016<-dplyr::select(chl_1975_2016,Date,"Station Code","Pennate Diatom","Centric Diatom","Chrysophyte","Cyanobacterium","Cryptophyte","Dinoflagellate","Euglenoid","Green Alga")

#join all chl data
DWR_chl<-chl_1975_2016%>%rbind(chl_2017)%>%rbind(chl_2018)%>%rbind(chl_2019)

DWR_chl<-dplyr::rename(DWR_chl,"SampleDate"="Date","DWRStationNo"="Station Code")

#calculate total phytoplankton organisms
DWR_chl$Total_orgs<-rowSums(DWR_chl[3:10])

#####################################
##Join zoop and chl data
#####################################

#There are two ways I'll do this
#1) Join by samples that were towed at the same station same day linking through DWRStationNo

Zoop_chl_station<-Zoop_all%>%inner_join(DWR_chl)

#select target taxa
target_adults<-c("ACARTELA","ACARTIA","EURYTEM","PDIAPFOR","SINOCAL","OITHDAV","LIMNOSPP","LIMNOSINE","LIMNOTET","OTHCYC","KERATELA","OTHROT","POLYARTH","SYNCH","SYNCHBIC","TRICHO","BOSMINA","DAPHNIA","DIAPHAN","OTHCLADO","H_longirostris","N_kadiakensis","N_mercedis")

Zoop_chl_station<-Zoop_chl_station%>%
  filter(ZooCode%in%target_adults)

station_sampling<-unique(dplyr::select(Zoop_chl_station,SampleDate,DWRStationNo))
saveRDS(station_sampling,"outputs/station_sampling.rds")
saveRDS(Zoop_chl_station,"outputs/Zoop_chl_Station.rds")

#2) Join by samples that were towed at the same month in the same region using EDSM strata
library(deltamapr)
library(spacetools)
library(sf)

#fit into deltamapr subregions for Brians adjusted EDSM subregions
WaterMap<-readRDS("Data/WaterMap.rds") #load watermap
stratum<-(R_EDSM_Subregions_Mahardja)%>%
  st_transform(4326)

#get DWR_chl station lat/long
dwr_stations$DWRStationNo<-dwr_stations$Station

#fixed zoop stations
zoop_stations<-unique(dplyr::select(Zoop_all,StationNZ,lon,lat))
zoop_stations<-zoop_stations[complete.cases(zoop_stations$lon),] #remove non-fixed stations
zoop_stations_geom<-st_as_sf(zoop_stations,coords = c("lon", "lat"),remove = TRUE, crs = 4326)

#fixed chl stations
chl_stations<-unique(dplyr::select(dwr_stations,Station,Longitude,Latitude))
chl_stations<-chl_stations[complete.cases(chl_stations$Longitude),] #remove non-fixed stations
chl_stations_geom<-st_as_sf(chl_stations,coords = c("Longitude", "Latitude"),remove = TRUE, crs = 4326)

p<-ggplot()+
  geom_sf(data=stratum,aes(fill=SubRegion))+
  geom_sf(data=spacetools::Delta)+
  geom_sf(data = zoop_stations_geom$geometry, size = 2, 
          shape = 23, fill = "green")+
  geom_sf(data = chl_stations_geom$geometry, size = 2, 
          shape = 24, fill = "blue")+
  theme(legend.position = "none") 
p
ggsave("Figures/station_subregions.png")

#st_join to assign stations to subregions
zoop_station_subregions<-st_join(zoop_stations_geom,stratum,join = st_within)
st_geometry(zoop_station_subregions)<-NULL
zoop_station_subregions<-unique((dplyr::select(zoop_station_subregions,StationNZ,SubRegion)))

chl_station_subregions<-st_join(chl_stations_geom,stratum,join = st_within)
st_geometry(chl_station_subregions)<-NULL
chl_station_subregions<-unique((dplyr::select(chl_station_subregions,Station,SubRegion)))

chl_station_subregions$DWRStationNo<-chl_station_subregions$Station

#do non-fixed ez stations next only for zoop, don't have DWR chl lat/long for ez stations
ez_stations<-c("NZEZ2","NZEZ6","NZES2SJR","NZEZ6SJR","EZ2","EZ6","EZ2-SJR","EZ6-SJR")
zoop_ez_stations<-read.csv("data/EZ_stations.csv")
zoop_ez_stations$SampleDate<-as.Date(zoop_ez_stations$SampleDate,"%m/%d/%Y")

zoop_ez<-Zoop_all%>%inner_join(zoop_ez_stations)
zoop_ez<-unique(dplyr::select(zoop_ez,StationNZ,SampleDate,Lat,Long))
zoop_ez_geom<-st_as_sf(zoop_ez,coords = c("Long", "Lat"),remove = TRUE, crs = 4326)


p<-ggplot()+
  geom_sf(data=stratum,aes(fill=SubRegion))+
  geom_sf(data=spacetools::Delta)+
  geom_sf(data = zoop_ez_geom$geometry, size = 2, 
          shape = 23, fill = "green")+
  theme(legend.position = "none") 
p
ez_station_subregions<-st_join(zoop_ez_geom,stratum,join = st_within)
st_geometry(ez_station_subregions)<-NULL
ez_station_subregions<-unique((dplyr::select(ez_station_subregions,StationNZ,SampleDate,SubRegion)))

#join in all subregion data to chl and zoop datasets
DWR_chl_regions<-DWR_chl%>%left_join(chl_station_subregions)
DWR_chl_regions$DWRSampleDate<-DWR_chl_regions$SampleDate

Zoop_all_regions<-Zoop_all%>%inner_join(zoop_station_subregions)

Zoop_all__ezregions<-Zoop_all%>%
  inner_join(ez_station_subregions)
Zoop_all_regions<-Zoop_all_regions%>%rbind(Zoop_all__ezregions)
Zoop_all_regions$ZoopSampleDate<-Zoop_all_regions$SampleDate

#select for necessary columns
DWR_chl_regions$Survey<-as.numeric(format(as.Date(DWR_chl_regions$DWRSampleDate), "%m"))
DWR_chl_regions$Year<-as.numeric(format(as.Date(DWR_chl_regions$DWRSampleDate), "%Y"))
DWR_chl_regions<-dplyr::select(DWR_chl_regions,DWRSampleDate,DWRStationNo,SubRegion,Year,Survey,3:11)

Zoop_all_regions<-dplyr::select(Zoop_all_regions,ZoopSampleDate,StationNZ,SubRegion,Year,Survey,ZooCode,CommonName,CPUE,Gear)

Zoop_chl_subregions<-Zoop_all_regions%>%inner_join(DWR_chl_regions)
saveRDS(Zoop_chl_subregions,"outputs/Zoop_chl_subregions.rds")
