### Script to calculate and plot annual parameters of ECMWF SF system5
### Authors: Dr Abdou Ali, Bernard Minoungou, Mohamed Hamatan, Jafet Andersson


  # -----------
  # 1 - General
  #Remove workspace
  rm(list=ls())
  library(raster)
  library(rgdal)
  require(lattice)
  #--------------
  #Get Niger Sub-bassin
  subs<-readOGR(dsn="C:/Niger_2.22/GIS", layer="Current_Niger_onedelta_subbasins_2012-10-17")
  subs1<-subs
  #--------------
  #Magnification for plots
  cexx<-1.5
  #--------------
  #Define SF file and outputs dir
  months<-vector()
  for (i in 1:12) {
    if (i<=9) {
      months<-c(months, paste("0", i, sep=""))
    } 
    else {
      months<-c(months, paste(i))
    }
  }

  datadir<-paste("I:/ecmwf_syst5_hype_format/D03_Niger_month", months, sep="")
  outdir<-paste("I:/analysis_ecmwf_syst5/D03_Niger_month",months, sep="")
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
  #Define variables names, variables attributes, years, plots lab, units
  var.name<-c("PObs", "Tobs", "TMAXobs", "TMINobs")
  var.long.name<-c("Seasonal Precipitation", "Seasonal Mean Temperature", 
                   "Seasonal Max Temperature", "Seasonal Min Temperature")
  units<-c("(mm)", rep("(°C)",3))
  var.num<-c(228, 55, 51, 52)
  years<-c(1993:2015)
  rm(i)
  # ----------------
  # 2 - Calculations and plots
  for (i in 1:length(months)) {
    for (j in 1:length(enss)) { 
      for (k in 1:length(years)) {
        for (l in 1:length(var.num)) {
          # 2.1 - Calculations
          fdir<-paste(datadir[i], years[k], sep="/")
          fname<-paste(var.name[l],"_",enss[j], ".txt", sep="")
          current.file<-read.table(paste(fdir, fname, sep="/"),h=T, sep="\t")
          if (var.name[l]=="PObs") {
            var.mean<-colSums(current.file[,-1])
          } else {
            var.mean<-colMeans(current.file[,-1])
          }
          a1<-slot(subs1, "data")
          a1<-cbind(a1,var.mean)
          slot(subs, "data")<-a1
          classes<-pretty(var.mean,100)
          my.colfun<-colorRampPalette(c("red","yellow","green","deepskyblue","blue"))
          if (j==1) {
            dir.create(outdir[i])
          }
          if( k==1) {
            dir.create(paste(outdir[i], "/ens", enss[j], sep="")) 
          }
          main1<-paste(var.long.name[l], years[k], "Ensemble", enss[j], units[l], sep=" ")
          # ----------------
          # 3.2. Maps (sub-basins)
          png(paste(outdir[i], "/ens", enss[j], "/",
                    paste("map", "_",var.name[l], "_",years[k],".png", sep=""),sep=""),
                    800,800)
          trellis.par.set(par.main.text=list(cex=cexx))
          print(
            spplot(subs, zcol="var.mean", 
                   at=classes, 
                   col.regions=my.colfun(length(classes)-1),
                   main=main1,
                   colorkey=list(labels=list(cex=cexx)),
                   scales=list(draw = TRUE))
          )
          dev.off()
          print(main1)
        }
      }
    }
  }
