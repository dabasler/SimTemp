rm(list = ls(all = TRUE))
while (length(dev.list())>0) dev.off()
op<-par();nop<-names(op);op<-op[nop!="cin" & nop!="cra" & nop!="csi" & nop!="cxy" &nop!="din"];rm("nop")
par(op)
####################################################################################
daylength<-function(lat,day,year){
  #Formula according to H.Glarner (http://herbert.gandraxa.com/herbert/lod.asp)
  # input: nothern latitude in degrees;julian date(s);year
  # output: Daylengthduration in hours
  pp<-NA
  #setconstants
  lat<- (pi/180)*lat				# deg/rad
  day<- day + 11 					# Correct for winter solistice
  day<-day +year%%4*0.25				# Correct for leapyears
  j<-pi/182.625
  axis<-(pi/180)*23.439				# earths ecliptic
  #calculate daylength
  for (i in 1:length(day)){
    m<-1-tan(lat)*tan(axis*cos(j*day[i])) 	#Exposed radius part between sun's zenith and sun's circle
    if (m<0){m<-0} 					# sun never appears
    if (m>2){m<-2}					# sun never disappears
    b<-acos(1-m)/pi					# Exposed fraction of the sun's circle
    pp[i]<-b*24						# Daylength (lat,day)
  }
  return (pp)
}


sm<-function(x,a,b,c){a*sin(((x-c)/182.5)*pi)+b} # Function Yearly slope

getTsim<-function(SM,SSD,dSD,sT=0){
  #SM Mittlerer Jahresverlauf der Temperatur
  #SDD Mittlerer Jahresverlauf der Temperatur-standardabweichung
  dT<-abs(rnorm(365,sd=dSD)) # Verteilung der Änderung von Tag zu Tag
  aT<-rnorm(1,sT,sd=SSD[1]) # Zufällige Start temperature um 0°C
  Tsim<-NA
  for (i in 1:365){
    swtch<--(rbinom(1,1,pnorm(aT,mean=0,sd=SSD[i]))*2-1) # Wahrscheinlichkeit zum vorzeichenwechsel von dT ist abhängig von der Temperatur-standardabweichung für den betreffenden Tag
    aT=aT+dT[i]*swtch
    Tsim<-c(Tsim,aT)
  }
  Tsim[-1]+SM
}

getTsimV<-function(SMV,SSD,dSD,sT=0){
  SM<-sm(c(1:365),rnorm(1,mean=SMV[1],sd=SMV[4]),rnorm(1,mean=SMV[2],sd=SMV[5]),rnorm(1,mean=SMV[3],sd=SMV[6]))
  return(getTsim(SM,SSD,dSD,sT=0))
}

getSM<-function(TY,DOY,plt=FALSE){ ### FIT THE YEARLY SEASONAL SLOPE
  nlf1<-nls(TY~sm(DOY,a,b,c),start = list(a = sd(TY), b = mean(TY), c = 110))
  SM<-sm(c(1:365),coef(nlf1)[1],coef(nlf1)[2],coef(nlf1)[3])
  if (plt==TRUE){
    plot(DOY,TY)
    lines(c(1:365),SM,col=2)
    mtext(side=3,bquote(.(round(coef(nlf1)[1],2))*sin(((x-.(round(coef(nlf1)[3],2))/183)*pi)+.(round(coef(nlf1)[2],2)))))}
  return(SM)
}

getSSD<-function(TY,DOY,fit,plt=FALSE){ # Get Seasonal SD Variation
  rp<-30
  if (!is.null(ncol(TY))){ # Allows to estimate mean SD across different Stations (TYwith multiple columns)
    TSD<-as.numeric(tapply(TY[,1],DOY,sd))
    for (i in 2:ncol(TY)){TSD<-cbind(TSD,as.numeric(tapply(TY[,i],DOY,sd)))}
    TSD<-apply(TSD,1,mean)
  }else{TSD<-as.numeric(tapply(TY,DOY,sd))}
  TSD<-c(TSD[(366-rp):365],TSD,TSD[1:rp])
  ssd<-lm(TSD ~ poly(c(1:length(TSD))-rp, fit))
  SSD<-as.numeric(predict(ssd))
  SSD<-SSD[(rp+1):(rp+365)]
  if (plt==TRUE){
    plot(c(1:365),TSD[(rp+1):(rp+365)],ylim=c(0,max(TSD)),las=1,xlab="DOY",ylab="SD(T)")
    lines(c(1:365),SSD,col=2)}
  return(SSD)
}

getSSDpar<-function(TY,DOY,fit,plt=FALSE){ # Get Seasonal SD Variation
  rp<-30
  if (!is.null(ncol(TY))){ # Allows to estimate mean SD across different Stations (TYwith multiple columns)
    TSD<-as.numeric(tapply(TY[,1],DOY,sd))
    for (i in 2:ncol(TY)){TSD<-cbind(TSD,as.numeric(tapply(TY[,i],DOY,sd,na.rm=TRUE)))}
    TSD<-apply(TSD,1,mean,na.rm=TRUE)
  }else{TSD<-as.numeric(tapply(TY,DOY,sd))}
  TSD<-c(TSD[(366-rp):365],TSD,TSD[1:rp])
  ssd<-lm(TSD ~ poly(c(1:length(TSD))-rp, fit,raw=TRUE))
  cf<-coef(ssd)
  names(cf)[2:(fit+1)]<-letters[1:fit]
  return (cf)
}

polynom<-function(x,cf){
  fit<-length(cf)-1
  y<-0
  for(i in 1:fit){y<-y+cf[i+1]*x^i}
  y<-y+cf[1]
  return(y)
}

getdT_SD<-function(x,plt=FALSE){ # Get SD of Delta T to next day
  dT<-c(c(x,NA)-c(NA,x))[2:length(x)]
  if (plt==TRUE) {DOY<-rep(c(1:365),length(x)/365);plot(DOY[-1],dT);lines(c(0,365),rep(sd(dT),2),col=2);lines(c(0,365),-rep(sd(dT),2),col=2)}
  sd(dT,na.rm=TRUE)
}

getSMpar<-function(TY,DOY,plt=FALSE){ ### FIT THE YEARLY SEASONAL SLOPE  
  nlf1<-nls(TY~sm(DOY,a,b,c),start = list(a = sd(TY), b = mean(TY), c = 110))
  return(coef(nlf1))
}

getSMV<-function(t,plt=FALSE){ # get Seasonal Slope Parameter and Variation across years
  smp<-rep(NA,3)
  for (i in 1:(length(t)/365)){smp<-rbind(smp,getSMpar(t[(1+(i-1)*365):(365+(i-1)*365)],c(1:365)))}
  psd<-apply(smp[-1,],2,sd)
  SMV<-c(getSMpar(t,rep(c(1:365),length(t)/365)),psd)
  names(SMV)<-c(letters[1:3],paste(letters[1:3],"sd",sep="_"))  
    if (plt==TRUE){
      xmx<-SMV[3]+182.5/2
      mx<-sm(round(SMV[3]+182.5/2),SMV[1],SMV[2],SMV[3])
      px<-seq(xmx-30,xmx+30,0.1)
      py<-seq((round(SMV[2])-5),(round(SMV[2])+5),0.1)
      lines(px,100*dnorm(px,xmx,SMV[6])+20,col=4) # Shift
      lines(-100*dnorm(py,SMV[2],SMV[5])+xmx+1,py,col=4)
      points(xmx,mx,col=4,pch=16)
      points(xmx,SMV[2],col=4,pch=16)
      points(xmx,2*SMV[2]-mx,col=4,pch=16)
      segments(xmx,mx,xmx,2*SMV[2]-mx,col=4)
      lines(50*dnorm(py,SMV[2],SMV[4]/2)+xmx+1,  py+mx/2,col=4)
      lines(50*dnorm(py,SMV[2],SMV[4]/2)+xmx+1,  py-mx/2,col=4)
    }
  return(SMV)
}


MonthlyMeans<-function(TD){
  ny<-ceiling(length(TD)/365)
  md<-c(1,cumsum(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),ny)))
  md<-md[md<=length(TD)]
  MM<-NA
  for (i in 1:(length(md)-1)){MM<-c(MM,mean(TD[md[i]:md[i+1]]))}
  MM<-MM[-1]
  if (length(TD)>md[length(md)]){MM<-c(MM,NA)}
  names(MM)<-rep(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","sep","Oct","Nov","Dec"),ny)[1:(length(MM))]
  return(MM)
}


################################

# Examples:
SimTemp.EXAMPLES<-function(){
  
  getTempDOY<-function(filename){
    T<-read.table(paste(inpath,filename,sep=""),header=TRUE,sep=";",skip=2)
    Tt<-strptime(T$time,"%Y%m%d")
    T$Year<-as.numeric(strftime(Tt,"%Y"))
    T$DOY<-as.numeric(strftime(Tt,"%j"))
    T<-T[,c(4,5,3)]
    names(T)[3]<-"TM"
    T<-T[T$DOY<=365,] # leap year day
    attr(T,"id")<-substr(filename,1,nchar(filename)-4)
    return (T)}
  
  # Get Some Temperatures
  inpath<-"C:/PHD/PhenologyModels/Data/ProcessDATA/PM_input/temperature/"
  Tfiles<-list.files(inpath,pattern="/*TG.txt")
  Tf<-getTempDOY(Tfiles[1])
  DOY<-Tf$DOY
  TY<-Tf[,3]
  # Get Some Temperatures
  x<-c(1:365)
  SSD<-getSSD(TY,DOY,6,plt=TRUE) # Get the SD~DOY function
  dSD<-getdT_SD(TY,TRUE)
  SM<-getSM(TY,DOY,TRUE) # Fixed estimate with similar seasonal slope across years
  SMV<-getSMV(TY,TRUE) # variable estimate with variable seasonal slope across years
  
  plot(DOY,TY,col="grey75",ylim=c(-20,30))
  lines(x,SM,"l")
  lines(x,SM,col="blue")
  lines(x,SM-SSD,col="red")
  lines(x,SM+SSD,col="red")
  lines(x,SM+1.65*SSD,col="orange")
  lines(x,SM-1.65*SSD,col="orange")
  lines(DOY[c(1:365)],getTsim(SM,SSD,dSD),"l")
  
  # Generate Years 
  for (i in 1:41){points(DOY[c(1:365)],getTsim(SM,SSD,dSD),col="blue")}
  for (i in 1:41){points(DOY[c(1:365)],getTsimV(SMV,SSD,dSD),col="red")}
  # Example Generate Data:
  n<-10
  TX<-NA;for (i in 1:n){TX<-c(TX,getTsim(SM,SSD,dSD))};TX<-TX[-1] # Similar seasonal course
  TV<-NA;for (i in 1:n){TX<-c(TX,getTsimV(SMV,SSD,dSD))};TV<-TV[-1] # Variable seasonal course
  TX<-NA;sT=0;for (i in 1:n){TX<-c(TX,getTsim(SM,SSD,dSD,sT));sT=TX[365]};TX<-TX[-1] # Similar seasonal course (continous throughout years)
  plot (TX,type="l")
  plot (MonthlyMeans(TX),type="l",ylim=c(-10,25))
  lines (MonthlyMeans(TY[1:(n*365)]),col=2,ylim=c(-10,25))
  
  #
  sT=0
  while (TRUE){
    Temperature<-getTsimV(SMV,SSD,dSD,sT)
    sT<-Temperature[365]
    plot(1:365,Temperature,ylim=c(-20,30),col="blue",type="l",las=1)
    lines(c(0,365),c(0,0))
    Sys.sleep(0.05)
  }
  ############################
  # Test if overall sd is similar to testset
  mean(getTsim(SM,SSD,dSD))
  sd(tapply(Tall[,3],Tall[,1],mean))
  ts<-NA;for (i in 1:41){ts<-c(ts,mean(getTsim(SM,SSD,dSD)))};sd(ts[-1]) # Too low
  ts<-NA;for (i in 1:41){ts<-c(ts,mean(getTsimV(SMV,SSD,dSD)))};sd(ts[-1]) # slightl too high low
  ############################
  # Show how good the normal distributions fits through out the year
  TY<-TM
  d<-2
  for (i in 1:round(365/d)){
    t<-TY[DOY<i*d &DOY>(i-1)*d]
    tt<-t-mean(t)
    dtt<-density(tt)
    plot(dtt$x,dtt$y/max(dtt$y),type="l",ylim=c(0,1),xlim=c(-15,15),main=i)
    dn<-dnorm(seq(-20,20,0.1),mean=0,sd=sd(t))
    lines(seq(-20,20,0.1),dn/max(dn),col="green")
    Sys.sleep(0.05)
  }
}
################################ FIND MEAN ESTIMATES FOR PARAMETERS
# inpath<-"C:/PHD/PhenologyModels/Data/ProcessDATA/PM_input/temperature/"
# Tfiles<-list.files(inpath,pattern="/*TG.txt")
# Tall<-getTempDOY(Tfiles[1])
# for (i in 1:length(Tfiles)){print (i);Tall<-cbind(Tall,getTempDOY(Tfiles[i])[,3])}
# names(Tall)[c(3:ncol(Tall))]<-c(1:(ncol(Tall)-2))
# DOY<-Tall$DOY
# Tlong<-rep(0,(ncol(Tall)-2)*nrow(Tall));DOYlong<-rep(0,(ncol(Tall)-2)*nrow(Tall));l<-nrow(Tall)
# for (i in 3:ncol(Tall)){Tlong[(1+(i-3)*l):(l+(i-3)*l)]<-Tall[,i];DOYlong[(1+(i-3)*l):(l+(i-3)*l)]<-Tall$DOY}
# #
# SMV<-getSMV(Tall[,3]);for (s in 4:ncol(Tall)){SMV<-rbind(SMV,getSMV(Tall[,s]))};SMV<-apply(SMV,2,mean)#mean Estimate across all sites
# #SSD<-getSSD(Tall[,3:ncol(Tall)],DOY,6,plt=TRUE) #mean Estimate across all sites
# SSDpar<-getSSDpar(Tall[,3:ncol(Tall)],DOY,6,plt=TRUE) #mean Estimate across all sites
# SSD<-polynom(1:365,SSDpar)
# dSD<-getdT_SD(Tlong,TRUE) #mean Estimate across all sites
# print (SMV)
# print (SSDpar)
# print (dSD)
# ################################ # RUN FROM SAVED VALUES
# SMV<-c(9.0118765,8.7910412,108.9776101,0.8927081,0.7626841,4.9142875)
# SSDpar<-c(4.705127e+00,5.265827e-03,-7.161727e-04,9.833857e-06,-5.584483e-08,1.413592e-10,-1.308046e-13)
# SSD<-polynom(1:365,SSDpar)
# dSD<-2.327802
# 
# TX<-getTsimV(SMV,SSD,dSD)
# plot(TX,type="l",ylim=c(-15,30),las=1,xlab="DOY",ylab="Daily Mean Temperature")
# lines (cumsum(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))-15,MonthlyMeans(TX),col=2)

# DETRENDED TEMEPRATURE
# 
# Tym<-tapply(Tf$TM,Tf$Year,mean)
# plot(c(1969:2009),Tym,"l")
# lm1<-lm(Tym~c(1969:2009))
# abline(lm1)
# Tyms<-Tym-(coef(lm1)[2]*c(1969:2009)+coef(lm1)[1])
# plot(c(1969:2009),Tyms,"l")
# sd(Tyms)
# lines(c(1969:2009),rnorm(41,0,sd(Tyms)),"l")
# 
# Tf$DOYf<-Tf$Year+Tf$DOY/365
# plot(Tf$DOYf,Tf$TM,"l",ylim=c(-20,30))
# lm1<-lm(Tf$TM~Tf$DOYf)
# abline(lm1)
# Tf$TMs<-Tf$TM-(coef(lm1)[2]*c(1969:2009)+coef(lm1)[1])
# plot(Tf$DOYf,Tf$TMs+ mean(Tf$TM),"l",ylim=c(-20,30))
# lines(Tf$DOYf,Tf$TMs+ mean(Tf$TM),"l",ylim=c(-20,30),col=2)
