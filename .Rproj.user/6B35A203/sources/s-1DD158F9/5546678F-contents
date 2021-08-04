#Figure production
#Created by Arthur Barros 
#EMP Zooplankton Enviornmental Scientist
#arthur.barros@wildlife.ca.gov
#last updated 8/3/2021
rm( list = ls()) #clear env
#read in packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(plotrix)
library(purrr)
library(cowplot)
library(plyr)
library(gridExtra)

Zoop_chl_stations<-readRDS("outputs/Zoop_chl_Station.rds")
Zoop_chl_subregions<-readRDS("outputs/Zoop_chl_subregions.rds")

target_adults<-c("ACARTELA","ACARTIA","EURYTEM","PDIAPFOR","SINOCAL","OITHDAV","LIMNOSPP","LIMNOSINE","LIMNOTET","BOSMINA","DAPHNIA","DIAPHAN","H_longirostris","N_kadiakensis","N_mercedis")

######################################
#create sampling coverage heat maps
###################################
samples_stations<-unique(dplyr::select(Zoop_chl_stations,Year,Survey,SampleDate,Gear,DWRStationNo))
samples_stations<-samples_stations%>%
  group_by(Year,Survey,Gear)%>%
  dplyr::summarize(tows=length(DWRStationNo))
samples_stations$Survey<-as.factor(samples_stations$Survey)

gear_types<-c("CB","Pump","Mysid")

plot_list<-list()
for(i in 1:3){
  target_gear<-gear_types[i]
  d<-samples_stations%>%filter(Gear==target_gear)
  p<-ggplot(d,aes(Year,Survey))+
    geom_tile(aes(fill=tows))+
    theme_classic()+
    theme(axis.text.x=element_text(angle=90,size=10),axis.text.y = element_text(size=10),legend.position = "none")+
    scale_fill_gradient(low="white",high="red")+
    geom_text(aes(label=round(tows,1)),size=4)+
    guides(colour=F)+
    ggtitle(paste(target_gear,"Sampling Coverage with DWR Chl data",sep=" "))
  plot_list[[i]]<-p
}
gear_plots<-plot_grid(plotlist=plot_list,labels = "AUTO", align = "v",ncol=1)
gear_plots
saveRDS(gear_plots,"Data/sampling_coverage_DWRStations.rds")
save_plot("Figures/chl_sampling_coverage_DWRStations.png", gear_plots,ncol=1,base_height = 14,base_width = 10)

#for regions now
samples_subregions<-unique(dplyr::select(Zoop_chl_subregions,Year,Survey,SubRegion,Gear))
samples_subregions<-samples_subregions%>%
  group_by(Year,Survey,Gear)%>%
  dplyr::summarize(tows=length(SubRegion))
samples_subregions$Survey<-as.factor(samples_subregions$Survey)


gear_types<-c("CB","Pump","Mysid")

plot_list<-list()
for(i in 1:3){
  target_gear<-gear_types[i]
  d<-samples_subregions%>%filter(Gear==target_gear)
  p<-ggplot(d,aes(Year,Survey))+
    geom_tile(aes(fill=tows))+
    theme_classic()+
    theme(axis.text.x=element_text(angle=90,size=10),axis.text.y = element_text(size=10),legend.position = "none")+
    scale_fill_gradient(low="white",high="red")+
    geom_text(aes(label=round(tows,1)),size=4)+
    guides(colour=F)+
    ggtitle(paste(target_gear,"Sampling Coverage with DWR Chl data",sep=" "))
  plot_list[[i]]<-p
}
gear_plots<-plot_grid(plotlist=plot_list,labels = "AUTO", align = "v",ncol=1)
gear_plots
saveRDS(gear_plots,"Data/sampling_coverage_SubRegions.rds")
save_plot("Figures/chl_sampling_coverage_SubRegions.png", gear_plots,ncol=1,base_height = 14,base_width = 10)

#############################
#Lets do some linear modeling
#############################
#this function creates the model from the dataset (d), then creates a character string with the equation, r-squared, and p-value. This is then used in the plotting
lmp <- function (modelobject) {
  coeffs<- coef(summary(modelobject))
  p<-if (dim(coeffs)[[1]]==1) # if only one row of pvalues then NA
    p.value <- NA
  else 
    p.value <- coeffs[2,4]
  return(p)
}

lm_eqn <- function(df){
  m <- lm(log(month_CPUE+.01) ~ log(month_chl+.01), df);
  eq <- substitute(~~italic(r)^2~"="~r2*","~~italic(p-value)*"="~pvalue,
                   list(intercept = format(coef(m)[1], digits = 2), 
                        slope = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3),
                        pvalue=format(lmp(m),digits=4)))
  as.character(as.expression(eq));
}


#let's run the data grouped by survey and year
#Run for loops for gams by taxa and chl species
winter_months<-c(12,1,2)

d<-dplyr::select(Zoop_chl_subregions,Year,Survey,DWRSampleDate,DWRStationNo,ZoopSampleDate,StationNZ,SubRegion,ZooCode,CommonName,CPUE,12:20)
d<-d%>%filter(!Survey%in%winter_months) #remove winter months
chl_species<-colnames(d[11:19])

for(i in 1:length(target_adults)){
  z<-target_adults[i]
  d1<-d%>%
    filter(ZooCode==z)
  for(s in 1:length(chl_species)){
    chl_taxa<-chl_species[s]
    d2<-dplyr::select(d1,Year,Survey,SubRegion,DWRSampleDate,DWRStationNo,ZoopSampleDate,StationNZ,CommonName,CPUE,matches(chl_taxa))
    d2<-dplyr::rename(d2,"chl_taxa"=10)
    
    #start grouping and averaging similar to FLOAT and DROUGHT averaging
    d3<-d2%>%
      group_by(Year,Survey,SubRegion)%>%
      dplyr::summarise(region_CPUE=mean(CPUE,rm.na=T),region_chl=mean(chl_taxa,rm.na=T))
    d4<-d3%>%
      group_by(Year,Survey)%>%
      dplyr::summarise(month_CPUE=mean(region_CPUE,rm.na=T),month_chl=mean(region_chl,rm.na=T))
    
    eq<-lm_eqn(d4)#apply linear model function
    p <- ggplot(data = d4, aes(x = log(month_chl+.01), y = log(month_CPUE+.01))) +
      geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
      geom_point(aes())+
      geom_smooth(se = T, method = "lm") +
      ggtitle(paste(z,chl_taxa,sep=" "))+
      theme_classic()+
      geom_text(x=-Inf,y=Inf,vjust="top",hjust="left",label=lm_eqn(d4),parse=T)
    p
    plot_list[[s]]<-p
  }
  taxa_plots<-plot_grid(plotlist=plot_list,align = "h",nrow=3)
  taxa_plots
  save_plot(paste("figures/",z,".png",sep=""),taxa_plots,base_height = 12,base_width = 12)
}
