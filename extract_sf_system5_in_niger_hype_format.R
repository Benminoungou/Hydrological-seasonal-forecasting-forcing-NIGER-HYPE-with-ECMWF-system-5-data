# Script to to extract ECMWF Seasonal forecast system 5 into Niger-HYPE V 2.22 format 
# Authors: Dr Abdou Ali, Bernard Minoungou, Mohamed Hamatan

#--------------
# 1 - General
  rm(list=ls())
  # Load librairies
  library(rgdal)
  library(raster)
  library(ncdf4)
  # Get Niger Sub-bassin shapefile and create centroids
  subs<-readOGR(dsn="C:/Niger_2.22/GIS", layer="Current_Niger_onedelta_subbasins_2012-10-17")
  centroid<-getSpPPolygonsLabptSlots(subs)
  centroid<-as.data.frame(centroid)
  colnames(centroid)<-c("x","y")
  centroid<-as.data.frame(centroid)
  coordinates(centroid)<-~x+y
  projection(centroid)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  
#-------------- 
# 2 - Conversion
  #--------------
  #Define SF file dir
  months<-vector()
  for (i in 1:12) {
      if (i<=9) {
      months<-c(months, paste("0", i, sep=""))
    } 
    else {
      months<-c(months, paste(i))
    }
  }
  
  datadir<-paste("I:/ECMWF_system5/D03_Niger_month", months, 
                 "/nobackup/smhid12/sm_jorro/GLORIOUS/Data/PreProcessSF/ECMWFSyst5/D03_Niger/m",
                 months, sep="")
  #--------------
  #Define SF ensembles index
  ens<-c(0:24)
  enss<-vector()
  for (i in 1:25) {
    if (i<=9) {
      enss<-c(enss, paste("00", i, sep=""))
    } 
    else {
      enss<-c(enss, paste("0", i, sep=""))
    }
  }
  #--------------
  #Define variables names, variables attributes and years
  var.name<-c("PObs", "Tobs", "TMAXobs", "TMINobs")
  var.num<-c(228, 55, 51, 52)
  years<-c(1993:2015)
  rm(i)
  #--------------
  #Get Data and save file into Niger-HYPE format
  for (i in 1:length(months)) {
    for (j in 1:length(ens)) {
      for (l in 1:length(years)) {
       for (k in 1:length(var.num)) {
          fname<-paste(var.num[k], "_D03_Niger_ECMWF-sys5_season_03deg_orig_reg_day_r", 
                       ens[j], "_", years[l],"-",months[i], ".nc",sep="")
          rr<-brick(paste(datadir[i], fname, sep="/"))
      
          #--------------
          #Extract and convert cumilated precipitation data into daily precipitation
          if (var.name[k]=="PObs") {
            get_Data.cum<-t(round(extract(rr, centroid)*1000, 5))
            get_Data<-matrix(NA, nrow=nrow(get_Data.cum), ncol=ncol(get_Data.cum))
            get_Data[1,]<-get_Data.cum[1,]
            for (m in 2:nrow(get_Data)) {
              get_Data[m,]<-round((get_Data.cum[m,]-get_Data.cum[m-1,]),2)
            }
            get_Data[which(get_Data<0)]<-0
            colnames(get_Data)<-slot(subs, "data")[, "SUBID"]
            Date<-paste(substr(rownames(get_Data.cum),2,5), substr(rownames(get_Data.cum),7,8),
                        substr(rownames(get_Data.cum),10,11), sep="-")
          } else {
            #--------------
            #Extract Temperature data and convert from Kelvin into °C
            get_Data<-t(round((extract(rr, centroid)-273.15), 2))
            colnames(get_Data)<-slot(subs, "data")[, "SUBID"]
            Date<-paste(substr(rownames(get_Data),2,5), substr(rownames(get_Data),7,8),
                        substr(rownames(get_Data),10,11), sep="-")
          }
          get_Data<-cbind(as.character(Date),get_Data)
          colnames(get_Data)[1]<-"DATE"
          if (k==1) {
            dir.create(paste("I:/ecmwf_syst5_hype_format/", "D03_Niger_month", months[i], sep=""))
          }
          dir.create(paste("I:/ecmwf_syst5_hype_format/", "D03_Niger_month", months[i], "/",years[l],sep=""))
          file.name<-paste(var.name[k], "_", enss[j], ".txt", sep="")
          current_file.dir<-paste("I:/ecmwf_syst5_hype_format/", "D03_Niger_month", months[i], "/",years[l],sep="")
          #--------------
          #Save data into Niger-HYPE format
          write.table(get_Data, paste(current_file.dir, file.name, sep="/"), col.names=T,
                      row.names=F, quote=F, sep="\t")
          process.progress<-paste(months[i], ens[j], years[l], var.name[k], sep="_")
          print(process.progress)
        }
      }
    }
  }
