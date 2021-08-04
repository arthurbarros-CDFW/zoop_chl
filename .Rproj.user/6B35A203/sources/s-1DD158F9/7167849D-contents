#taxa and chl gams
rm( list = ls()) #clear env

library(mgcv)
library(hobbsutils)
library(tidyverse)
library(MASS)
library(corrplot)
library(readxl)
library(voxel)
library(lubridate)
library(tidymv)
library(cowplot)

Zoop_chl<-readRDS("outputs/Zoop_chl.rds")

#Run for loops for gams by taxa and chl species
gam_data<-dplyr::select(Zoop_chl,Year,Survey,DWRStationNo,SampleDate,ZooCode,CommonName,CPUE,27:35)

target_adults<-c("ACARTELA","ACARTIA","EURYTEM","PDIAPFOR","SINOCAL","OITHDAV","LIMNOSPP","LIMNOSINE","LIMNOTET","BOSMINA","DAPHNIA","DIAPHAN","H_longirostris","N_kadiakensis","N_mercedis")

chl_species<-colnames(gam_data[8:16])

plot_list<-list()
for(i in 1:length(target_adults)){
  z<-target_adults[i]
  d1<-gam_data%>%
    filter(ZooCode==z)
  for(s in 1:length(chl_species)){
    chl_taxa<-chl_species[s]
    d2<-dplyr::select(d1,Year,Survey,DWRStationNo,SampleDate,CommonName,CPUE,matches(chl_taxa))
    d2<-rename(d2,"chl_taxa"=7)
    gam_m1<-gam(round(CPUE)~s(chl_taxa),data=d2)
    model_p<-predict_gam(gam_m1)
    sb=summary(gam_m1)
    xrng <- range(model_p$chl_taxa)
    yrng <- range(model_p$fit)
    p<-model_p%>%
      ggplot(aes(chl_taxa,fit))+geom_smooth_ci()+
      ggtitle(paste(z,chl_taxa,sep=" "))+
      annotate("text",x=xrng[1],y=seq(yrng[2],yrng[1],length.out = 1),label=paste("p = ",signif(c(sb$s.table[4])),hjust = 0))
    p
    plot_list[[s]]<-p
  }
  taxa_plots<-plot_grid(plotlist=plot_list,align = "h",nrow=1)
  taxa_plots
  save_plot(paste("Figures/GAMS/",z,chl_taxa,"ziP.png",sep="_"),taxa_plots,base_height = 6,base_width = 24)
}
