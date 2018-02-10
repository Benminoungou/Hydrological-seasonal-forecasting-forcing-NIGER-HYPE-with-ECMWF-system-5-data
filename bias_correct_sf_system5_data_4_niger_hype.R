###---------------------------------------------------------------------------
###Script to bias correct the ECMWF SF System 5 with the WFDEI reanalys data  
###Authors: Bernard Minoungou, Abdou Ali, Mohamed Hamatan                     
###---------------------------------------------------------------------------

#-----------
#1 - General 
  
  #-------------
  # 1.1. Clear work directory and load packages
    rm(list=ls())
    library(raster)
    library(rgdal)
    require(lattice)
    library(qmap)

  #--------------
  # 1.2. Magnification for plots
    cexx<-1.5
    
    
  #--------------
  # 1.3. Define SF system 5 dir
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
    
  #-------------- 
  #1.5. Défine corrected data, bias parameters and Q-Q plot dir
    biaisdir<-paste("I:/ecmwf_syst5_hype_format/Biais/D03_Niger_month", 
                    months, sep="")
    plotdir<-paste("I:/analysis_ecmwf_syst5/Biais/plot/D03_Niger_month",
                   months, sep="")
    biasparamdir<-paste("I:/analysis_ecmwf_syst5/Biais/Bias_parameters/D03_Niger_month",
                        months, sep="")
  #--------------
  #1.6. Define SF ensembles index
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
  # 1.7. Define variables names, variables attributes and years
    var.name<-c("PObs", "Tobs", "TMAXobs", "TMINobs")
    var.name<-c("PObs", "Tobs", "TMAXobs", "TMINobs")
    var.long.name<-c("Precipitation", "Mean Temperature", 
                 "Max Temperature", "Min Temperature")
    units<-c("[mm]", rep("[°C]",3))
    var.num<-c(228, 55, 51, 52)
    years<-c(1993:2015)
    rm(i)

#-----------
#2 - Get WFDEI Data (Pobs, TObs, TMAXObs, TMINObs) 
  wfdei<-list()
  pcp<-read.table("C:/Niger_2.22/Pobs.txt", sep="\t",h=T)
  tobs<-read.table("C:/Niger_2.22/Tobs.txt", sep="\t",h=T)
  tmaxobs<-read.table("C:/Niger_2.22/TMAXobs.txt", sep="\t",h=T)
  tminobs<-read.table("C:/Niger_2.22/TMINobs.txt", sep="\t",h=T)
  wfdei[[1]]<-pcp
  wfdei[[2]]<-tobs
  wfdei[[3]]<-tmaxobs
  wfdei[[4]]<-tminobs
  
#-----------
#3 - Bias correction and save corrected data, bias parameters and Q-Q Plots
for (i in 1:length(months)) {
  dir.create(biaisdir[i]) # Bias output dir
  dir.create(plotdir[i]) # Bias correction QQ Plot dir by initialization month 
  dir.create(biasparamdir[i]) # Bias parameter sub-dir by initialization month
  for (j in 1:length(var.num)) {
    for (k in 1:length(enss)) {
      value.ens.all.year<-NULL
      #-----------
      #3.1 - loop through years
      for (l in 1:length(years)) {
        if (j==2&k==1) {
          dir.create(paste(biaisdir[i], years[l], sep="/"))
        }
        fdir<-paste(datadir[i], years[l], sep="/")
        fname<-paste(var.name[j],"_",enss[k], ".txt", sep="")
        current.file<-read.table(paste(fdir, fname, sep="/"),h=T, sep="\t")
        value.ens.all.year<-rbind(value.ens.all.year,current.file)
      }
      date.common<-intersect(value.ens.all.year[,1],wfdei[[j]][,1])
      mymatch1<-match(date.common,value.ens.all.year[,1])
      mymatch2<-match(date.common,wfdei[[j]][,1])
      value.ens.intersect<-value.ens.all.year[mymatch1,]
      value.wfdei.intersect<-wfdei[[j]][mymatch2,]
      
      #3.2 - Compute bias parameter for precipitation data
      if (j==1) {
        current.fit<-fitQmap(obs = value.wfdei.intersect[,-1],
                             mod =value.ens.intersect[,-1],
                             method ="QUANT", 
                             type="linear", 
                             wet.day=TRUE)
      } else {
        #3.3 - Compute bias parameter for temperature data
        current.fit<-fitQmap(obs = value.wfdei.intersect[,-1],
                             mod =value.ens.intersect[,-1],
                             method ="QUANT", 
                             type="linear")
      }
      
      
      #3.4 - Apply bias correction to historical data
      value.ens.corr<-round(doQmap(value.ens.intersect[,-1], current.fit),2)
      
      
      ###3.5 Save bias parameter for operationnal correction
      save(current.fit, 
           file=paste(biasparamdir[i],"/",var.name[j], "_ens_", enss[k], ".RData",sep=""))
      
      ###3.6 Plot quantile-quantile map to evaluate the correction
      png(paste(plotdir[i], "/QQplot_", var.name[j], "_ens_", enss[k], ".png", sep=""),  800,800)
      xymin<-min(min(rowMeans(value.wfdei.intersect[,-1])),
                min(rowMeans(value.ens.corr, na.rm=T)),
                min(rowMeans(value.ens.intersect[,-1], na.rm=T)))
      xymax<-max(max(rowMeans(value.wfdei.intersect[,-1])),
                max(rowMeans(value.ens.corr, na.rm=T)),
                max(rowMeans(value.ens.intersect[,-1], na.rm=T)))
      qqplot(rowMeans(value.wfdei.intersect[,-1]), 
             rowMeans(value.ens.corr, na.rm=T),
             xlim=c(xymin,xymax), ylim=c(xymin,xymax),
             col="green",
             xlab=paste("Observed",  var.long.name[j], units[j], sep=" "), 
             ylab=paste("Forecasted",  var.long.name[j], units[j], sep=" "), 
             main=paste("Q-Q Plot", var.long.name[j], "Ensemble", enss[k], sep=" "),
             cex.lab=cexx,cex.main=cexx,cex.axis=cexx
             )
      points(sort(rowMeans(value.wfdei.intersect[,-1], na.rm=T)), 
             sort(rowMeans(value.ens.intersect[,-1], na.rm=T)), 
             col = "red")
      abline(0,1, col="grey", lwd=2)
      legend("topleft", c("Corrected", "Uncorrected"), 
             col=c("green", "red"), pch=c(1,1), cex=1.5)
      dev.off()
      
      ###3.7 Save corrected data year by year
      value.ens.all.year.cor<-round(doQmap(value.ens.all.year[,-1], current.fit),2)
      value.ens.all.year.cor<-cbind(value.ens.all.year[,1],value.ens.all.year.cor)
      colnames(value.ens.all.year.cor)<-c("DATE", sub("X","",colnames(value.ens.all.year)[-1]))
      for (m in 1:length(years)) {
        current.save.file<-value.ens.all.year.cor[which(substr(value.ens.all.year.cor[,1],1,4)==years[m]),]
        write.table(current.save.file,
                    paste(biaisdir[i], years[m], paste(var.name[j], "_", enss[k], ".txt", sep=""), sep="/"),
                    sep="\t", col.names=T, row.names=F, quote=F)
      }
      print(paste(months[i],"_", var.name[j], "_ens_", enss[k], sep=""))
    }
  }
}
