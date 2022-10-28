# to process raw USGS flow discharge
library(zoo)
library(data.table)
library(ggplot2)

dv<-read.table('dv_new',header=TRUE,fill=TRUE) # dv_new is daily data downloaded from USGS


# 1969-01-01 to 2018-12-31 # complete 18262?
siteid0<-(unique(dv$site_no))

siteid<-(unique(dv$site_no))
xchar<-nchar(as.character(siteid))
siteid<-siteid[xchar>=8]
siteid<-siteid[!grepl('_cd',siteid)]
siteid<-siteid[!is.na(siteid)] # site id
unique(nchar(as.character(siteid)))

siteid<-siteid[nchar(as.character(siteid))<=15] # site id



splitdv<-split(dv,dv$site_no)
splitdv1<-splitdv[names(splitdv)%in%siteid]
rm(splitdv)


### daily noise

allnoise<-lapply(splitdv1, allusgsnoise_day)

pval<-unlist(allnoise)[seq(2,length(unlist(allnoise)),64)]

noise<--1*unlist(allnoise)[seq(1,length(unlist(allnoise)),64)]
noise0<-noise[noise<10]
site0<-unlist(siteid[noise<10])

xx<-as.data.frame(matrix(unlist(allnoise),ncol=64,byrow = T))
xx<-xx[noise<10,]
names(xx)<-c('noise','p_val_global','noise_7_365','p_val_7_365',
             'noise_7_30','p_val_7_30','noise_7_90','p_val_7_90',
             'noise_30_365','p_val_30_365','noise_30_180','p_val_30_180','noise_365_3650','p_val_365_3650',
             'noise_30_90','pval_30_90','noise_90_180','pval_90_180','noise_180_365', 'pval_180_365',
             'noiser','p_val_globalr','noise_7_365r','p_val_7_365r',
             'noise_7_30r','p_val_7_30r','noise_7_90r','p_val_7_90r',
             'noise_30_365r','p_val_30_365r','noise_30_180r','p_val_30_180r','noise_365_3650r','p_val_365_3650r',
             'noise_30_90r','pval_30_90r','noise_90_180r','pval_90_180r','noise_180_365r', 'pval_180_365r',
             'noisenot','p_val_globalnot','noise_7_365not','p_val_7_365not',
             'noise_7_30not','p_val_7_30not','noise_7_90not','p_val_7_90not',
             'noise_30_365not','p_val_30_365not','noise_30_180not','p_val_30_180not','noise_365_3650not','p_val_365_3650not',
             'noise_30_90not','pval_30_90not','noise_90_180not','pval_90_180not','noise_180_365not', 'pval_180_365not',
             'data lenght_global','missing rate','startdate','enddate')
xx$noise<-xx$noise*(-1)
xx$noise_7_365<-xx$noise_7_365*(-1)
xx$noise_7_30<-xx$noise_7_30*(-1)
xx$noise_7_90<-xx$noise_7_90*(-1)
xx$noise_30_365<-xx$noise_30_365*(-1)
xx$noise_30_180<-xx$noise_30_180*(-1)
xx$noise_365_3650<-xx$noise_365_3650*(-1)
xx$noise_30_90<-xx$noise_30_90*(-1)
xx$noise_90_180<-xx$noise_90_180*(-1)
xx$noise_180_365<-xx$noise_180_365*(-1)

xx$noiser<-xx$noiser*(-1)
xx$noise_7_365r<-xx$noise_7_365r*(-1)
xx$noise_7_30r<-xx$noise_7_30r*(-1)
xx$noise_7_90r<-xx$noise_7_90r*(-1)
xx$noise_30_365r<-xx$noise_30_365r*(-1)
xx$noise_30_180r<-xx$noise_30_180r*(-1)
xx$noise_365_3650r<-xx$noise_365_3650r*(-1)
xx$noise_30_90r<-xx$noise_30_90r*(-1)
xx$noise_90_180r<-xx$noise_90_180r*(-1)
xx$noise_180_365r<-xx$noise_180_365r*(-1)

xx$noisenot<-xx$noisenot*(-1)
xx$noise_7_365not<-xx$noise_7_365not*(-1)
xx$noise_7_30not<-xx$noise_7_30not*(-1)
xx$noise_7_90not<-xx$noise_7_90not*(-1)
xx$noise_30_365not<-xx$noise_30_365not*(-1)
xx$noise_30_180not<-xx$noise_30_180not*(-1)
xx$noise_365_3650not<-xx$noise_365_3650not*(-1)
xx$noise_30_90not<-xx$noise_30_90not*(-1)
xx$noise_90_180not<-xx$noise_90_180not*(-1)
xx$noise_180_365not<-xx$noise_180_365not*(-1)

xx$startdate<-as.Date(xx$startdate)
xx$enddate<-as.Date(xx$enddate)
xx$siteid<-site0

xxx<-nchar(as.character(site0))
xx1<-xx[xxx==8,]
xx2<-xx[xxx==9,]
xx3<-xx[xxx==10,]
xx4<-xx[xxx==15,]
xx0<-rbind(xx1,xx2,xx3,xx4)

write.csv(xx0,'usgs7504.csv',row.names = F)



allusgsnoise_yr<-function(sample){
  #  sample<-splitdv[[2]]
  da2<-sample[,3:4]
  names(da2)=c('time','flow')
  da2$time=as.Date(da2$time,'%Y-%m-%d')
  da2$flow<-as.numeric(as.character(da2$flow))
  da2<-da2[complete.cases(da2),]
  
  ddd=as.numeric(da2$time)
  ddd0=abs(ddd[length(ddd)]-ddd[length(ddd)-1])
  if(ddd0<1){
    da2<-aggregate(flow~time,da2,FUN=mean)  ## convert subdaily to daily
  }
  NL=(as.numeric(da2$time[length(da2$time)])-as.numeric(da2$time[1]))
  NL0<-nrow(da2)
  missingratio<-1-NL0/NL
  aaa=as.Date(0:NL,origin=da2$time[1])
  if(missingratio>=0.05| NL0<=50*365){
    xnoise=c(-99,-99)
    xnoise_max=c(-99,-99)
    xnoise_min=c(-99,-99)
  }else{
    A=as.data.frame(aaa)
    names(A)=c('time')
    da2=data.frame(da2)
    da3=merge(A,da2,by='time',all=TRUE)
    da3$flow=na.approx(da3$flow) # daily flow
    
    da3$yr<-format(da3$time,format='%Y')
    da3$yr<-as.numeric(da3$yr)
    da4<-aggregate(flow~yr,da3,mean) # aggregate to annual mean flow
    
    xnoise<-cal_noise_yr(da4$flow) # in case the flow records are not ending or start in a regular date
    
    da5<-aggregate(flow~yr,da3,max) # aggregate to annual max flow
    xnoise_max<-cal_noise_yr(da5$flow)
    
    da6<-aggregate(flow~yr,da3,min) # aggregate to annual min flow
    xnoise_min<-cal_noise_yr(da6$flow)
    
  }
  result<-unname(cbind(xnoise[1],xnoise[2],xnoise_max[1],xnoise_max[2],xnoise_min[1],xnoise_min[2],
                       NL0,missingratio,
                       as.numeric(da2$time[1]),as.numeric(da2$time[length(da2$time)])))
  return(result)
}


allusgsnoise_day<-function(sample){
  #  sample<-splitdv[[2]]
  da2<-sample[,3:4]
  names(da2)=c('time','flow')
  da2$time=as.Date(da2$time,'%Y-%m-%d')
  da2$flow<-as.numeric(as.character(da2$flow))
  da2<-da2[complete.cases(da2),]
  
  ddd=as.numeric(da2$time)
  ddd0=abs(ddd[length(ddd)]-ddd[length(ddd)-1])
  if(ddd0<1){
    da2<-aggregate(flow~time,da2,FUN=mean)  ## convert subdaily to daily
  }
  NL=(as.numeric(da2$time[length(da2$time)])-as.numeric(da2$time[1]))
  NL0<-nrow(da2)
  missingratio<-1-NL0/NL
  aaa=as.Date(0:NL,origin=da2$time[1])
  if(missingratio>=0.05| NL0<=15*365){
    xnoise=c(-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99)
    xnoise_raw=c(-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99)
    xnoise_notrend=c(-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99)
    
  }else{
    A=as.data.frame(aaa)
    names(A)=c('time')
    da2=data.frame(da2)
    da3=merge(A,da2,by='time',all=TRUE)
    da3$flow=na.approx(da3$flow)
    xnoise<-cal_noise_day(da3$flow)
    xnoise_raw<-cal_noise_day_raw(da3$flow)
    xnoise_notrend<-cal_noise_day_notrend(da3$flow)
  }
  result<-(c(xnoise,xnoise_raw,xnoise_notrend,NL0,missingratio,
             as.numeric(da2$time[1]),as.numeric(da2$time[length(da2$time)])))
  return(result)
}

## remove linear trend & periodic seasonal component
getrandom_stl<-function(x){
  tx<-ts(x,frequency=365.25)
  yyy<-(stl(tx,s.window="periodic",t.window = length(x)))$time.series[,3]
  return(yyy)
}

library(pracma)
##### calculate power spectrum of annual flows
cal_noise_yr<-function(x){
  x=na.approx(x)
  
  x<-x[2:(length(x)-1)] # in case the flow does not end or start as expected.
  
  cutoff=c(0,1)

  if(length(unique(x))==1){ # in case all flows are 0
    slope=-99
    pval=-99
  }else{
    data=spectrum(detrend(x),log='no',plot=FALSE)
    ind0=floor(length(data$spec)*cutoff)
    ind=max(ind0[1],1):ind0[2]
    
    regfit=lm(log10(data$spec)~log10(data$freq))
    slope=summary(regfit)$coefficients[2,1]
    pval=summary(regfit)$coefficients[2,4]
  }
  return(c(slope,pval))
}


##### calculate power spectrum of daily flows (removing seasonal,trend)
cal_noise_day<-function(x){
  x=na.approx(x)
  data=spectrum(getrandom_stl(x),log='no',plot=FALSE)
  
  regfit_global=lm(log10(data$spec)~log10(data$freq))
  slope_global=summary(regfit_global)$coefficients[2,1]
  pval_global=summary(regfit_global)$coefficients[2,4]
  
  x4<-as.data.frame(cbind(data$freq,data$spec))
  names(x4)<-c('freq','spec')
  x4$day<-1/x4$freq # it converted to year after stl
  x5<-x4[x4$day<=365/365 & x4$day>=7/365,]
  regfit_7_365=lm(log10(spec)~log10(freq),data=x5)
  slope_7_365=summary(regfit_7_365)$coefficients[2,1]
  pval_7_365=summary(regfit_7_365)$coefficients[2,4]
  
  x6<-x4[x4$day<=30/365 & x4$day>=7/365,]
  regfit_7_30=lm(log10(spec)~log10(freq),data=x6)
  slope_7_30=summary(regfit_7_30)$coefficients[2,1]
  pval_7_30=summary(regfit_7_30)$coefficients[2,4]
  
  x7<-x4[x4$day<=90/365 & x4$day>=7/365,]
  regfit_7_90=lm(log10(spec)~log10(freq),data=x7)
  slope_7_90=summary(regfit_7_90)$coefficients[2,1]
  pval_7_90=summary(regfit_7_90)$coefficients[2,4]
  
  x8<-x4[x4$day<=365/365 & x4$day>=30/365,]
  regfit_30_365=lm(log10(spec)~log10(freq),data=x8)
  slope_30_365=summary(regfit_30_365)$coefficients[2,1]
  pval_30_365=summary(regfit_30_365)$coefficients[2,4]
  
  x9<-x4[x4$day<=180/365 & x4$day>=30/365,]
  regfit_30_180=lm(log10(spec)~log10(freq),data=x9)
  slope_30_180=summary(regfit_30_180)$coefficients[2,1]
  pval_30_180=summary(regfit_30_180)$coefficients[2,4]
  
  x10<-x4[x4$day<=365*10/365 & x4$day>=365/365,]
  regfit_365_3650=lm(log10(spec)~log10(freq),data=x10)
  slope_365_3650=summary(regfit_365_3650)$coefficients[2,1]
  pval_365_3650=summary(regfit_365_3650)$coefficients[2,4]
  
  x11<-x4[x4$day<=90/365 & x4$day>=30/365,]
  regfit_30_90=lm(log10(spec)~log10(freq),data=x11)
  slope_30_90=summary(regfit_30_90)$coefficients[2,1]
  pval_30_90=summary(regfit_30_90)$coefficients[2,4]
  
  x12<-x4[x4$day<=180/365 & x4$day>=90/365,]
  regfit_90_180=lm(log10(spec)~log10(freq),data=x12)
  slope_90_180=summary(regfit_90_180)$coefficients[2,1]
  pval_90_180=summary(regfit_90_180)$coefficients[2,4]
  
  x13<-x4[x4$day<=365/365 & x4$day>=180/365,]
  regfit_180_365=lm(log10(spec)~log10(freq),data=x13)
  slope_180_365=summary(regfit_180_365)$coefficients[2,1]
  pval_180_365=summary(regfit_180_365)$coefficients[2,4]
  
  return(c(slope_global,pval_global,slope_7_365,pval_7_365,
           slope_7_30,pval_7_30,slope_7_90,pval_7_90,
           slope_30_365,pval_30_365,slope_30_180,pval_30_180,
           slope_365_3650,pval_365_3650,
           slope_30_90,pval_30_90,slope_90_180,pval_90_180,slope_180_365, pval_180_365));
}


### raw data with removing trend, but keeping seasonal signal
cal_noise_day_notrend<-function(x){
  x=na.approx(x)
  data=spectrum(x,log='no',plot=FALSE,detrend=TRUE)
  regfit_global=lm(log10(data$spec)~log10(data$freq))
  slope_global=summary(regfit_global)$coefficients[2,1]
  pval_global=summary(regfit_global)$coefficients[2,4]
  
  x4<-as.data.frame(cbind(data$freq,data$spec))
  names(x4)<-c('freq','spec')
  x4$day<-1/x4$freq # it converted to day
  x5<-x4[x4$day<=365 & x4$day>=7,]
  regfit_7_365=lm(log10(spec)~log10(freq),data=x5)
  slope_7_365=summary(regfit_7_365)$coefficients[2,1]
  pval_7_365=summary(regfit_7_365)$coefficients[2,4]
  
  x6<-x4[x4$day<=30 & x4$day>=7,]
  regfit_7_30=lm(log10(spec)~log10(freq),data=x6)
  slope_7_30=summary(regfit_7_30)$coefficients[2,1]
  pval_7_30=summary(regfit_7_30)$coefficients[2,4]
  
  x7<-x4[x4$day<=90 & x4$day>=7,]
  regfit_7_90=lm(log10(spec)~log10(freq),data=x7)
  slope_7_90=summary(regfit_7_90)$coefficients[2,1]
  pval_7_90=summary(regfit_7_90)$coefficients[2,4]
  
  x8<-x4[x4$day<=365 & x4$day>=30,]
  regfit_30_365=lm(log10(spec)~log10(freq),data=x8)
  slope_30_365=summary(regfit_30_365)$coefficients[2,1]
  pval_30_365=summary(regfit_30_365)$coefficients[2,4]
  
  x9<-x4[x4$day<=180 & x4$day>=30,]
  regfit_30_180=lm(log10(spec)~log10(freq),data=x9)
  slope_30_180=summary(regfit_30_180)$coefficients[2,1]
  pval_30_180=summary(regfit_30_180)$coefficients[2,4]
  
  x10<-x4[x4$day<=365*10 & x4$day>=365,]
  regfit_365_3650=lm(log10(spec)~log10(freq),data=x10)
  slope_365_3650=summary(regfit_365_3650)$coefficients[2,1]
  pval_365_3650=summary(regfit_365_3650)$coefficients[2,4]
  
  x11<-x4[x4$day<=90 & x4$day>=30,]
  regfit_30_90=lm(log10(spec)~log10(freq),data=x11)
  slope_30_90=summary(regfit_30_90)$coefficients[2,1]
  pval_30_90=summary(regfit_30_90)$coefficients[2,4]
  
  x12<-x4[x4$day<=180 & x4$day>=90,]
  regfit_90_180=lm(log10(spec)~log10(freq),data=x12)
  slope_90_180=summary(regfit_90_180)$coefficients[2,1]
  pval_90_180=summary(regfit_90_180)$coefficients[2,4]
  
  x13<-x4[x4$day<=365 & x4$day>=180,]
  regfit_180_365=lm(log10(spec)~log10(freq),data=x13)
  slope_180_365=summary(regfit_180_365)$coefficients[2,1]
  pval_180_365=summary(regfit_180_365)$coefficients[2,4]
  
  
  return(c(slope_global,pval_global,slope_7_365,pval_7_365,
           slope_7_30,pval_7_30,slope_7_90,pval_7_90,
           slope_30_365,pval_30_365,slope_30_180,pval_30_180,
           slope_365_3650,pval_365_3650,
           slope_30_90,pval_30_90,slope_90_180,pval_90_180,slope_180_365, pval_180_365));
}

### raw data without removing trend/mean/seasonal cycle
cal_noise_day_raw<-function(x){
  x=na.approx(x)
  data=spectrum(x,log='no',plot=FALSE,detrend=FALSE,demean=FALSE)
  
  regfit_global=lm(log10(data$spec)~log10(data$freq))
  slope_global=summary(regfit_global)$coefficients[2,1]
  pval_global=summary(regfit_global)$coefficients[2,4]
  
  x4<-as.data.frame(cbind(data$freq,data$spec))
  names(x4)<-c('freq','spec')
  x4$day<-1/x4$freq # it converted to year after stl
  x5<-x4[x4$day<=365 & x4$day>=7,]
  regfit_7_365=lm(log10(spec)~log10(freq),data=x5)
  slope_7_365=summary(regfit_7_365)$coefficients[2,1]
  pval_7_365=summary(regfit_7_365)$coefficients[2,4]
  
  x6<-x4[x4$day<=30 & x4$day>=7,]
  regfit_7_30=lm(log10(spec)~log10(freq),data=x6)
  slope_7_30=summary(regfit_7_30)$coefficients[2,1]
  pval_7_30=summary(regfit_7_30)$coefficients[2,4]
  
  x7<-x4[x4$day<=90 & x4$day>=7,]
  regfit_7_90=lm(log10(spec)~log10(freq),data=x7)
  slope_7_90=summary(regfit_7_90)$coefficients[2,1]
  pval_7_90=summary(regfit_7_90)$coefficients[2,4]
  
  x8<-x4[x4$day<=365 & x4$day>=30,]
  regfit_30_365=lm(log10(spec)~log10(freq),data=x8)
  slope_30_365=summary(regfit_30_365)$coefficients[2,1]
  pval_30_365=summary(regfit_30_365)$coefficients[2,4]
  
  x9<-x4[x4$day<=180 & x4$day>=30,]
  regfit_30_180=lm(log10(spec)~log10(freq),data=x9)
  slope_30_180=summary(regfit_30_180)$coefficients[2,1]
  pval_30_180=summary(regfit_30_180)$coefficients[2,4]
  
  x10<-x4[x4$day<=365*10 & x4$day>=365,]
  regfit_365_3650=lm(log10(spec)~log10(freq),data=x10)
  slope_365_3650=summary(regfit_365_3650)$coefficients[2,1]
  pval_365_3650=summary(regfit_365_3650)$coefficients[2,4]
  
  x11<-x4[x4$day<=90 & x4$day>=30,]
  regfit_30_90=lm(log10(spec)~log10(freq),data=x11)
  slope_30_90=summary(regfit_30_90)$coefficients[2,1]
  pval_30_90=summary(regfit_30_90)$coefficients[2,4]
  
  x12<-x4[x4$day<=180 & x4$day>=90,]
  regfit_90_180=lm(log10(spec)~log10(freq),data=x12)
  slope_90_180=summary(regfit_90_180)$coefficients[2,1]
  pval_90_180=summary(regfit_90_180)$coefficients[2,4]
  
  x13<-x4[x4$day<=365 & x4$day>=180,]
  regfit_180_365=lm(log10(spec)~log10(freq),data=x13)
  slope_180_365=summary(regfit_180_365)$coefficients[2,1]
  pval_180_365=summary(regfit_180_365)$coefficients[2,4]
  
  
  return(c(slope_global,pval_global,slope_7_365,pval_7_365,
           slope_7_30,pval_7_30,slope_7_90,pval_7_90,
           slope_30_365,pval_30_365,slope_30_180,pval_30_180,
           slope_365_3650,pval_365_3650,
           slope_30_90,pval_30_90,slope_90_180,pval_90_180,slope_180_365, pval_180_365));
}




library(randomForest)
library(hydroGOF)
library(ggplot2)
library(dplyr)
library(scales)
library(ggmap)
library(ppcor)
library(scales)
library(grid)
library(gridExtra)

rfnoisedata_final1<-read.csv('new7504.csv')  


noise_new<-xx0[,c('noise','siteid')]
names(noise_new)<-c('noise_new','STAID')

noise_new$STAID<-as.numeric(as.character(noise_new$STAID))

rfnoise1<-merge(noise_new,rfnoisedata_final1,by='STAID')
names(rfnoise1)[2:3]<-c('noise','noise_old')

features<-c('DOF','DOR', 'ORD_CLAS', 'urban2006', 'forest2006', 'wetland2006',
            'wateruse2005', 'tmean_mean', 'tmean_cv', 'ppt_mean', 'ppt_cv',
            'DRAIN_SQKM', 'Elev_m')


rfnoise2<-rfnoise1[,c(features,'noise')]


### after the cross-validation
set.seed(23)
train3=sample(1:nrow(rfnoise2),floor(0.8*nrow(rfnoise2)))

cv_rf<-rfcv(rfnoise2[train3,-which(names(rfnoise2)%in%'noise')],rfnoise2$noise[train3],cv.fold = 10)
with(cv_rf, plot(n.var, error.cv))
bestmtry <- tuneRF(rfnoise2[train3,-which(names(rfnoise2)%in%'noise')],rfnoise2$noise[train3], stepFactor=1.5, improve=1e-4, ntree=2000)
print(bestmtry)

rfmodel31_new<-randomForest(noise~., data=rfnoise2,subset=train3,importance=TRUE,
                            ntree=1000,na.omit=na.omit,mtry=6)
rfmodel31_new
plot(rfmodel31_new)
varImpPlot(rfmodel31_new)
rfmodel3<-rfmodel31_new
plot(rfmodel31_new$predicted,rfmodel31_new$y)

test3<-seq(1,nrow(rfnoise2),1)[-train3]
rfmodel31_newb<-predict(rfmodel31_new,rfnoise2[test3,])
cor(rfnoise2$noise[test3],rfmodel31_newb)

rftrain1<-as.data.frame(cbind(rfmodel31_new$predicted,rfmodel31_new$y))
names(rftrain1)<-c('Predict','Obs')
rftest1<-as.data.frame(cbind(rfnoise2$noise[test3],rfmodel31_newb))
names(rftest1)<-c('Predict','Obs')

cor(rftrain1$Obs,rftrain1$Predict)

ggplot(rftrain1,aes(x=Obs,y=Predict))+
  geom_point()+theme_minimal()+
  geom_abline(intercept = 0,slope=1,color='red',size=2,linetype='dashed')+
  xlab('Noise color (Observed)')+
  ylab('Noise color (Predicted)')+
  xlim(c(0,3))+ylim(c(0,3))+ggtitle('')+
  theme(text=element_text(family='serif',size=14))+
  theme(axis.line = element_line(colour = "black"))

cor(rftest1$Obs,rftest1$Predict)

ggplot(rftest1,aes(x=Obs,y=Predict))+
  geom_point()+theme_minimal()+
  geom_abline(intercept = 0,slope=1,color='red',size=2,linetype='dashed')+
  xlab('Noise color (Observed)')+
  ylab('Noise color (Predicted)')+
  xlim(c(0,3))+ylim(c(0,3))+ggtitle('')+
  theme(text=element_text(family='serif',size=14))+
  theme(axis.line = element_line(colour = "black"))


rftrain1$Diff<-rftrain1$Predict-rftrain1$Obs

ggplot(rftrain1,aes(x=Diff))+theme_minimal()+
  geom_histogram(colour='grey',fill='white')+
  scale_x_continuous(breaks=seq(-3,3,0.5))+
  xlab('Prediction error')+ylab('Number of gauges')+
  theme(axis.line = element_line(colour = "black"))+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank()) +theme(axis.text=element_text(size=14),axis.title=element_text(size=14))+
  theme(text=element_text(family='serif',size=16))

rftest1$Diff<-rftest1$Predict-rftest1$Obs

ggplot(rftest1,aes(x=Diff))+theme_minimal()+
  geom_histogram(colour='grey',fill='white')+
  scale_x_continuous(breaks=seq(-3,3,0.5))+
  xlab('Prediction error')+ylab('Number of gauges')+
  theme(axis.line = element_line(colour = "black"))+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank()) +theme(axis.text=element_text(size=14),axis.title=element_text(size=14))+
  theme(text=element_text(family='serif',size=16))




rmse(rfmodel31_new$y,rfmodel31_new$predicted)
rmse(rfmodel31_newb,rfnoise2$noise[test3])

pbias(rfmodel31_new$y,rfmodel31_new$predicted)
pbias(rfmodel31_newb,rfnoise2$noise[test3])

cor(rfmodel31_new$predicted,rfmodel31_new$y)
plot(rfnoise2$noise[test3],rfmodel31_newb,xlim=c(0,3.5),ylim=c(0,3.5))
abline(a=0,b=1)



rf_train<-rfnoise1[train3,]
rf_train$noise_predict<-rfmodel31_new$predicted


rf_test<-rfnoise1[-train3,]
rf_test$noise_predict<-rfmodel31_newb


rfdaily_1<-data.frame(cbind(rfmodel31_new$y,rfmodel31_new$predicted))
names(rfdaily_1)<-c('Daily_train_obs','Daily_train_pred')
rfdaily_1$Obs<-'obs'
rfdaily_1$pred<-'pred'
x1<-rfdaily_1[,c('Daily_train_obs','Obs')]
x2<-rfdaily_1[,c('Daily_train_pred','pred')]
names(x1)<-c('noise','Type')
names(x2)<-c('noise','Type')
x12<-as.data.frame(rbind(x1,x2))
ggplot(data=x12,aes(x=noise,fill=Type))+
  geom_histogram(position='identity',alpha=0.2,binwidth = 0.1)

rfdaily_2<-data.frame(rfmodel31_newb,rfnoise2$noise[test3])
names(rfdaily_2)<-c('Daily_test_obs','Daily_test_pred')
rfdaily_2$Obs<-'obs'
rfdaily_2$pred<-'pred'
x11<-rfdaily_2[,c('Daily_test_obs','Obs')]
x22<-rfdaily_2[,c('Daily_test_pred','pred')]
names(x11)<-c('noise','Type')
names(x22)<-c('noise','Type')
x122<-as.data.frame(rbind(x11,x22))
ggplot(data=x122,aes(x=noise,fill=Type))+
  geom_histogram(position='identity',alpha=0.2,binwidth = 0.1)

############################# 
daily_obs<-(c(rfmodel31_new$y,rfnoise2$noise[test3]))
daily_pred<-(c(rfmodel31_new$predicted,rfmodel31_newb))
daily_all<-data.frame(cbind(daily_obs,daily_pred))

daily_train<-data.frame(cbind(rfmodel31_new$y,rfmodel31_new$predicted))
daily_test<-data.frame(cbind(rfnoise2$noise[test3],rfmodel31_newb))

names(daily_train)<-c('obs','pred')
names(daily_test)<-c('obs','pred')
daily_train$Type<-'Train'
daily_test$Type<-'Test'
daily_all0<-(rbind(daily_train,daily_test))
daily_all0$obs_cat<-ifelse(daily_all0$obs<(-0.5),1,
                           ifelse(daily_all0$obs<0.5,2,
                                  ifelse(daily_all0$obs<1.5,3,
                                         ifelse(daily_all0$obs<2.5,4,5))))

daily_all0$pred_cat<-ifelse(daily_all0$pred<(-0.5),1,
                            ifelse(daily_all0$pred<0.5,2,
                                   ifelse(daily_all0$pred<1.5,3,
                                          ifelse(daily_all0$pred<2.5,4,5))))

daily_all0$match<-ifelse(daily_all0$pred_cat==daily_all0$obs_cat,1,0)

sum(daily_all0$match[daily_all0$Type=='Train'])/sum(daily_all0$Type=='Train')
sum(daily_all0$match[daily_all0$Type=='Test'])/sum(daily_all0$Type=='Test')


hist(daily_all0$obs_cat)
hist(daily_all0$pred_cat)

daily_bar_train<-melt(daily_all0[daily_all0$Type='Train',c('obs_cat','pred_cat')])


ggplot(data=daily_bar_train)+
  geom_bar(data=daily_bar_train[daily_bar_train$variable=='obs_cat',],
           aes(x=as.factor(value),fill=as.factor(value)),stat='Count')+
  geom_bar(data=daily_bar_train[daily_bar_train$variable=='pred_cat',],aes(x=as.factor(value),
                                                                           fill=as.factor(value)),stat='Count',position = position_dodge(width=0.1))+
  scale_fill_manual(values=c('2'='grey85','3'='pink','4'='red','5'='black'))


ggplot(data=daily_all)+
  geom_histogram(aes(daily_obs,fill='white'),alpha=0.2,binwidth = 0.25,colour=1)+
  geom_histogram(aes(daily_pred,fill='grey'),alpha=0.2,binwidth=0.25)+
  xlab('Daily noise color')+ylab('Number of gauges')+
  scale_fill_manual(name='',values=c('white'='white','grey'='grey'),labels=c('white'='Observed','grey'='Predicted'))+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   legend.position = c(0.1,0.8),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank()) +theme(axis.text=element_text(size=12,family='serif'),axis.title=element_text(size=14,family='serif'))






features<-c('DOF','DOR', 'ORD_CLAS', 'urban2006', 'forest2006', 'wetland2006',
            'wateruse2005', 'tmean_mean', 'tmean_cv', 'ppt_mean', 'ppt_cv',
            'DRAIN_SQKM', 'Elev_m')
figure_data1<-rfnoise2[,(features)]
rf_imp<-varImpPlot(rfmodel31_new) # mtry=6
feature_imp<-rf_imp
figure_data1$noise<-rfnoise2$noise

names(figure_data1)[4:6]<-c('LULC_urban','LULC_forest','LULC_wetland')
names(figure_data1)[7:11]<-c('Anthro_wateruse','Hydro_tmp','Hydro_tmpCV','Hydro_ppt','Hydro_pptCV')
names(figure_data1)[1:3]<-c('Anthro_DOF','Anthro_DOR','Hydro_stream')
names(figure_data1)[12:13]<-c('Geo_area','Geo_elevation')

varcor<-cor(figure_data1)
dat_cor1<-cbind(varcor[1:(nrow(varcor)-1),nrow(varcor)],feature_imp)
dat_cor1<-as.data.frame(dat_cor1)

dat_cor1$`%IncMSE`<-feature_imp[,1] # IncMSE
dat_cor1$gini<-feature_imp[,2] # IncNodePurity


dat_cor1$Type<-c('Anthro','Anthro','Hydro','LULC',
                 'LULC','LULC',
                 'Anthro','Hydro','Hydro','Hydro','Hydro',
                 'Geo','Geo')


data_pie<-dat_cor1[which(names(dat_cor1)%in%c('gini','Type'))]
names(data_pie)<-c('Importance','group')

df<-aggregate(.~group,data_pie,sum)
df <- data.frame(value = df$Importance,
                 Group = df$group) %>%
  mutate(Group = factor(Group, levels = sort(df$group,decreasing = T)),
         cumulative = cumsum(value),
         midpoint = cumulative - value / 2,
         label=paste0(" ",c(round(value[1]/ sum(value) * 100, 0),
                            #         label = paste0(Group, " ", round(value / sum(value) * 100, 1), "%"))
                            #        label=paste0(Group," ",c(round(value[1]/ sum(value) * 100, 0),
                            100-round(sum(value[1])/ sum(value) * 100, 0)-round(sum(value[3])/ sum(value) * 100, 0)-
                              round(sum(value[4])/ sum(value) * 100, 0),round(value[3]/ sum(value) * 100, 0),
                            round(value[4]/ sum(value) * 100, 0)),"%"))

pieplot1<-ggplot(df, aes(x = 1, weight = value, fill = Group)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y") +
  scale_fill_manual(name =" ", values = c('Anthro'='aliceblue', "Geo" = "azure2",'Hydro'='lavender','LULC'='lightcyan2'))+
  #    scale_fill_brewer("Blues") +
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  geom_text(aes(x = 1.2, y = midpoint, label = label),size=8,family='serif')+theme_nothing()




names(dat_cor1)[1]<-'correlation'
dat_cor<-dat_cor1[(order(dat_cor1$gini,decreasing=TRUE)),]
color_cor1<-colorRampPalette(c("purple", "grey", "blue"))
color_cor<-color_cor1(nrow(dat_cor))

colo1<-as.data.frame(cbind(dat_cor1$correlation,dat_cor1$gini))
names(colo1)[1:2]<-c('correlation','mse')
colo1$mse<-colo1$mse/sum(colo1$mse)


### if we use partial correlation in the plot
colo1$correlation<-pcor(rfnoise2)$estimate[14,1:13]  # partial correlation coefficient between the variable and noise

############################
colo2<-colo1[(order(colo1$correlation,decreasing=TRUE)),]
colo2$color_cor<-color_cor
colo2<-colo2[(order(colo2$mse,decreasing=TRUE)),]
name1<-(row.names(dat_cor1))[as.numeric(row.names(colo2))]

names(colo2)[1]<-'PartialCorr'

#figure_rf<-ggplot(data=colo2,aes(x=reorder(row.names(colo2),mse),mse,fill=correlation))+
figure_rf<-ggplot(data=colo2,aes(x=reorder(row.names(colo2),mse),mse,fill=PartialCorr))+
  geom_bar(stat='identity')+scale_x_discrete(labels = rev(name1))+
  #  scale_x_discrete(labels = rev(name1),position='top')+
  scale_fill_gradient2(low='blue',mid='lightgoldenrodyellow',high='red',midpoint = 0)+coord_flip()+
  #  scale_y_reverse(limits=c(0.24,0))+
  ylab('Relative importance')+xlab('')+ylim(c(0,0.25))+
  theme(axis.line.x = element_line(color="black", size = 0.5))+
  theme(axis.text=element_text(size=12,colour='black'))+
  theme(legend.position = c(0.8,0.75))+  #position of color bar
  theme(line = element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.key.width=unit(0.5,"line"),legend.key.height=unit(0.8,"line"))



figure_rf<-figure_rf+annotation_custom(ggplotGrob(pieplot1), xmin = 0, xmax =8, 
                                       ymin = 0.0, ymax = 0.3)


parfeature<-features[c(12,8,6,4,1,2)]

par_drainage1<-partialPlot(rfmodel3,rfnoise2,parfeature[1])
par_temp1<-partialPlot(rfmodel3,rfnoise2,parfeature[2])
par_wetland1<-partialPlot(rfmodel3,rfnoise2,parfeature[3])
par_urban1<-partialPlot(rfmodel3,rfnoise2,parfeature[4])
par_DOF1<-partialPlot(rfmodel3,rfnoise2,parfeature[5])
par_DOR1<-partialPlot(rfmodel3,rfnoise2,parfeature[6])
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

par_drainage2<-ggplot() +
  geom_line(data=as.data.frame(par_drainage),aes(x=range01(x), y=y))


colorspar <- c('Drainage','Temperature','Wetland','Urban','DOR','DOF')
parnoise<-ggplot() +
  geom_line(data=as.data.frame(par_drainage1),aes(x=range01(x), y=y,color='c1'),size=3)+
  geom_line(data=as.data.frame(par_temp1),aes(x=range01(x), y=y,color='c2'),size=2.5)+
  geom_line(data=as.data.frame(par_wetland1),aes(x=range01(x), y=y,color='c3'),linetype = "dotdash",size=2)+
  geom_line(data=as.data.frame(par_urban1),aes(x=range01(x), y=y,color='c4'),linetype = "dashed",size=2)+
  geom_line(data=as.data.frame(par_DOR1),aes(x=range01(x), y=y,color='c5'),size=4)+
  geom_line(data=as.data.frame(par_DOF1),aes(x=range01(x), y=y,color='c6'),size=2)+
  scale_color_manual('',values=c('c1'="azure2",'c2'='lavender','c3'='lightcyan2','c4'='aliceblue','c5'='aliceblue','c6'='aliceblue'),
                     labels = c('Area','Tmp','Wet','Urb','DOR','DOF'))+
  #  scale_color_manual('',values=c('c1'='red','c2'='blue','c3'='black','c4'='lightblue','c5'='grey85','c6'='pink'),
  #                     labels = c('Area','Tmp','Wet','Urb','DOR','DOF'))+
  xlab('')+ylab('Noise color')+  guides(color = guide_legend(nrow = 1))+
  theme(legend.key.height = unit(0.05,'cm' ))+
  theme(legend.key.width =unit(0.02,'cm' ))+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank()) +theme(axis.text=element_text(size=12),axis.title=element_text(size=12)) +
  theme(legend.position ='none')+
  annotate(geom='text',x=0.08,y=1.85,label='Area',color='red')+
  annotate(geom='text',x=0.08,y=1.65,label='Tmp',color='blue')+
  annotate(geom='text',x=0.33,y=1.7,label='Wetland',color='black')+
  annotate(geom='text',x=0.25,y=1.52,label='DOR',color='grey85')+
  annotate(geom='text',x=0.75,y=1.35,label='Urban',color='lightblue')+
  annotate(geom='text',x=0.09,y=1.35,label='DOF',color='pink')
##########################################################################################
################## annual noise
new2597<-read.csv('J:/usgs_flow/noise2597_all.csv')

new2597<-new2597[,c(features,'noise')] # mean noise

set.seed(0)
train33=sample(1:nrow(new2597),floor(0.8*nrow(new2597)))
### mean annual noise
cv_rf2<-rfcv(new2597[train33,-which(names(new2597)%in%'noise')],new2597$noise[train33],cv.fold = 10)
with(cv_rf2, plot(n.var, error.cv))
bestmtry2 <- tuneRF(new2597[train33,-which(names(new2597)%in%'noise')],new2597$noise[train33], stepFactor=1.5, improve=1e-4, ntree=1000)
print(bestmtry2)

rfmodel33<-randomForest(noise~., data=new2597,subset=train33,importance=TRUE,ntree=1000,na.omit=na.omit,mtry=6,set.seed(623))
rfmodel33
varImpPlot(rfmodel33)
features<-c('DOF','DOR', 'ORD_CLAS', 'urban2006', 'forest2006', 'wetland2006',
            'wateruse2005', 'tmean_mean', 'tmean_cv', 'ppt_mean', 'ppt_cv',
            'DRAIN_SQKM', 'Elev_m')
figure_data2<-new2597[,(features)]
rf_imp_yr<-varImpPlot(rfmodel33) # 
feature_imp_yr<-rf_imp_yr
figure_data2$noise<-new2597$noise
names(figure_data2)[4:6]<-c('LULC_Urban','LULC_forest','LULC_wetland')
names(figure_data2)[7:11]<-c('Anthro_wateruse','Hydro_tmp','Hydro_tmpCV','Hydro_ppt','Hydro_pptCV')
names(figure_data2)[1:3]<-c('Anthro_DOF','Anthro_DOR','Hydro_stream')
names(figure_data2)[12:13]<-c('Geo_area','Geo_elevation')
varcor_yr<-cor(figure_data2)
dat_cor_yr<-cbind(varcor_yr[1:(nrow(varcor_yr)-1),nrow(varcor_yr)],feature_imp_yr)
dat_cor_yr<-as.data.frame(dat_cor_yr)
#dat_cor_yr$`%IncMSE`<-feature_imp # this is feature importance extracted from Python 
dat_cor_yr$`%IncMSE`<-feature_imp_yr[,1] # IncMSE
dat_cor_yr$gini<-feature_imp_yr[,2] # IncNodePurity

dat_cor_yr$Type<-c('Anthro','Anthro','Hydro','LULC',
                   'LULC','LULC',
                   'Anthro','Hydro','Hydro','Hydro','Hydro',
                   'Geo','Geo')

data_pie_yr<-dat_cor_yr[which(names(dat_cor_yr)%in%c('gini','Type'))]
names(data_pie_yr)<-c('Importance','group')

df2<-aggregate(.~group,data_pie_yr,sum)
df2 <- data.frame(value = df2$Importance,
                  Group = df2$group) %>%
  mutate(Group = factor(Group, levels = sort(df2$group,decreasing = T)),
         cumulative = cumsum(value),
         midpoint = cumulative - value / 2,
         label = paste0(round(value / sum(value) * 100, 0), "%"))
#df2$label[3]<-'49%'  # summation to 1 (rounding error)


names(dat_cor_yr)[1]<-'correlation'
#dat_cor<-dat_cor_yr[(order(dat_cor_yr$`%IncMSE`,decreasing=TRUE)),]
dat_cor<-dat_cor_yr[(order(dat_cor_yr$gini,decreasing=TRUE)),]
color_cor1<-colorRampPalette(c("purple", "grey", "blue"))
color_cor<-color_cor1(nrow(dat_cor))
#colo1_yr<-as.data.frame(cbind(dat_cor_yr$correlation,dat_cor_yr$`%IncMSE`))
#names(colo1_yr)[1:2]<-c('correlation','mse')
colo1_yr<-as.data.frame(cbind(dat_cor_yr$correlation,dat_cor_yr$gini))
names(colo1_yr)[1:2]<-c('correlation','mse')
colo1_yr$mse<-colo1_yr$mse/sum(colo1_yr$mse)


############### partial correlation
colo1_yr$correlation<-pcor(figure_data2)$estimate[14,1:13]  # partial correlation coefficient between the variable and noise
###############################################

colo2_yr<-colo1_yr[(order(colo1_yr$correlation,decreasing=TRUE)),]
colo2_yr$color_cor<-color_cor
colo2_yr<-colo2_yr[(order(colo2_yr$mse,decreasing=TRUE)),]
name1_yr<-(row.names(dat_cor_yr))[as.numeric(row.names(colo2_yr))]

names(colo2_yr)[1]<-'PartialCorr'
colo2_yr$var<-name1_yr
names(colo2_yr)<-c('Partial2','mse2','color_cor2','var')




DATA<-colo2
DATA$var<-name1
DATA$type<-c('G','H','L','H','G','H','H','L','L','A','A','A','H') # urban is LULC

DATA2<-DATA[,c('mse','var','type')]
colo2_yr2<-colo2_yr[,c('var','mse2')]
colo2_yr2$var[9]<-'Anthro_urban'  ### check which one 9 or 8 to rename
# #DATA2$var[8]<-'Anthro_urban'
# #DATA2$type[8]<-'A'  
DATA2$var[8]<-'Anthro_urban' 


DATA<-merge(DATA2,colo2_yr2,all.x=TRUE,by='var',sort=FALSE)


#DATA$var[1]<-'Geo_Area'
#DATA$var[5]<-'Geo_Elevation'

# DATA$var<-c('Area','Temperature','Temperature CV','Wetland Cover',
#             'Elevation','Precipitation','Precipitation CV','Dam Regulation',
#             'Urban Cover','Forest Cover','Water Use','Fragmentation','Stream Order')
DATA$var<-c('Area','Temperature','Wetland Cover','Temperature CV',
            'Elevation','Precipitation','Precipitation CV','Urban Cover',
            'Forest Cover','Dam Regulation','Water Use','Fragmentation','Stream Order')

DATA$var<-factor(DATA$var,levels=DATA[order(DATA$mse),'var'])


g.mid<-ggplot(DATA,aes(x=1,y=var))+geom_text(aes(label=var),family='serif',size=5)+

  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.935,1.065))+xlab('')+
  theme(axis.title.y=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"))+theme(axis.text=element_text(size=12,family='serif'),axis.title=element_text(size=14,family = 'serif'))

g1 <- ggplot(data = DATA, aes(x = var, y = mse,fill=type)) +
  scale_fill_manual(name =" ", values = c('A'='antiquewhite', "G" = "azure2",'H'='lavender','L'='lightsalmon2'))+
  geom_bar(stat = "identity") +ggtitle('(a)')+
  theme(#axis.title.x = element_blank(), 
    plot.title=element_text(hjust=0,size=14,family = 'serif'),
    axis.title.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    legend.position = 'none',
    plot.margin = unit(c(1,-1,1,0), "mm")) +
  scale_y_reverse(limits = c(0.25,0)) + coord_flip()+ylab('Predictor importance')+
  theme(axis.line.x = element_line(colour = 'black'),
        panel.background = element_blank())+theme(axis.text=element_text(size=14,family='serif'),axis.title=element_text(size=14,family='serif'))

g2 <- ggplot(data = DATA, aes(x = var, y = mse2,fill=type)) +
  scale_fill_manual(name =" ", values = c('A'='antiquewhite', "G" = "azure2",'H'='lavender','L'='lightsalmon2'))+
  xlab(NULL)+
  geom_bar(stat = "identity") + ggtitle("(b)") +
  theme(plot.title=element_text(hjust=0,size=14,family = 'serif'),
        legend.position = 'none',#axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,0,1,-1), "mm")) +
  coord_flip()+ylab('Predictor importance')+ylim(c(0,0.25))+
  theme(axis.line.x = element_line(colour = 'black'),
        panel.background = element_blank())+theme(axis.text=element_text(size=14,family='serif'),axis.title=element_text(size=14,family='serif'))



pieplot_daily<-ggplot(df, aes(x = 1, weight = value, fill = Group)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y") +
  scale_fill_manual(name =" ", values = c('Anthro'='antiquewhite', "Geo" = "azure2",'Hydro'='lavender','LULC'='lightsalmon2'))+
  #    scale_fill_brewer("Blues") +
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  geom_text(aes(x = 1.2, y = midpoint, label = label),size=5,family='serif')+theme_nothing()+
  geom_text(aes(x = 1.88, y = midpoint, label = c('Water use','Geo','Hydro','Land use')),size=5,family='serif')

pieplot_yr<-ggplot(df2, aes(x = 1, weight = value, fill = Group)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y") +
  scale_fill_manual(name =" ", values = c('Anthro'='antiquewhite', "Geo" = "azure2",'Hydro'='lavender','LULC'='lightsalmon2'))+
  #  scale_fill_manual(name =" ", values = c('Anthro'='antiquewhite', "Geo" = "azure2",'Hydro'='lavender','LULC'='lightsalmon2'))+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  geom_text(aes(x = 1.2, y = midpoint, label = label),size=5,family='serif')+theme_nothing()+
  geom_text(aes(x = 1.88, y = midpoint, label = c('Water use','Geo','Hydro','Land use')),size=5,family='serif')

g01<-g1+annotation_custom(ggplotGrob(pieplot_daily), xmin = 0, xmax =6, 
                          ymin = -0.27, ymax = -0.12)


g02<-g2+annotation_custom(ggplotGrob(pieplot_yr), xmin = 0, xmax =6, 
                          ymin = 0.12, ymax = 0.27)


gg1 <- ggplot_gtable(ggplot_build(g01))
gg2 <- ggplot_gtable(ggplot_build(g02))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9))



parnoise_daily<-ggplot() +ggtitle('(c)')+
  theme(plot.title=element_text(hjust=0,size=14,family = 'serif'))+
  geom_line(data=as.data.frame(par_drainage1),aes(x=range01(x), y=y,color='c1'),size=3)+
  geom_line(data=as.data.frame(par_temp1),aes(x=range01(x), y=y,color='c2'),size=2.5)+
  geom_line(data=as.data.frame(par_wetland1),aes(x=range01(x), y=y,color='c3'),linetype = "dotdash",size=2)+
  geom_line(data=as.data.frame(par_urban1),aes(x=range01(x), y=y,color='c4'),linetype = "longdash",size=2.5)+
  geom_line(data=as.data.frame(par_DOR1),aes(x=range01(x), y=y,color='c5'),size=3.5)+
  #  geom_line(data=as.data.frame(par_DOF1),aes(x=range01(x), y=y,color='c6'),size=2)+
  #  scale_color_manual('',values=c('c1'="azure2",'c2'='lavender','c3'='lightsalmon2','c4'='antiquewhite','c5'='antiquewhite','c6'='antiquewhite'),
  #                     labels = c('Area','Tmp','Wet','Urb','DOR','DOF'))+
  scale_color_manual('',values=c('c1'="azure2",'c2'='lavender','c3'='lightsalmon2','c4'='antiquewhite','c5'='antiquewhite'),
                     labels = c('Area','Tmp','Wet','Urb','DOR'))+
  #  scale_color_manual('',values=c('c1'='red','c2'='blue','c3'='black','c4'='lightblue','c5'='grey85','c6'='pink'),
  #                     labels = c('Area','Tmp','Wet','Urb','DOR','DOF'))+
  xlab('')+ylab('Noise color')+  guides(color = guide_legend(nrow = 1))+
  theme(legend.key.height = unit(0.05,'cm' ))+
  theme(legend.key.width =unit(0.02,'cm' ))+xlab('Predictor')+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank()) +theme(axis.text=element_text(size=12,family='serif'),axis.title=element_text(size=14,family='serif')) +
  theme(legend.position ='none')+
  annotate(geom='text',x=0.08,y=1.85,label='Area',color="black",family='serif')+
  annotate(geom='text',x=0.35,y=1.65,label='Temperature',color='black',family='serif')+
  annotate(geom='text',x=0.8,y=1.75,label='Wetland',color='black',family='serif')+
  annotate(geom='text',x=0.85,y=1.52,label='Dam Regulation',color='black',family='serif')+
  annotate(geom='text',x=0.75,y=1.35,label='Urban',color='black',family='serif')


rfnoisedata_final2<-figure_data2
names(rfnoisedata_final2)<-names(rfnoise2)
par_drainage<-partialPlot(rfmodel33,rfnoisedata_final2,parfeature[1])
par_temp<-partialPlot(rfmodel33,rfnoisedata_final2,parfeature[2])
par_wetland<-partialPlot(rfmodel33,rfnoisedata_final2,parfeature[3])
par_urban<-partialPlot(rfmodel33,rfnoisedata_final2,parfeature[4])
par_DOF<-partialPlot(rfmodel33,rfnoisedata_final2,parfeature[5])
par_DOR<-partialPlot(rfmodel33,rfnoisedata_final2,parfeature[6])

parnoise_annual<-ggplot() +ggtitle('(d)')+
  theme(plot.title=element_text(hjust=0,size=14,family = 'serif'))+
  geom_line(data=as.data.frame(par_drainage),aes(x=range01(x), y=y,color='c1'),size=3)+
  geom_line(data=as.data.frame(par_temp),aes(x=range01(x), y=y,color='c2'),size=2.5)+
  geom_line(data=as.data.frame(par_wetland),aes(x=range01(x), y=y,color='c3'),linetype = "dotdash",size=2)+
  geom_line(data=as.data.frame(par_urban),aes(x=range01(x), y=y,color='c4'),linetype = "longdash",size=2.5)+
  geom_line(data=as.data.frame(par_DOR),aes(x=range01(x), y=y,color='c5'),size=3.5)+
  scale_color_manual('',values=c('c1'="azure2",'c2'='lavender','c3'='lightsalmon2','c4'='antiquewhite','c5'='antiquewhite'),
                     labels = c('Area','Tmp','Wet','Urb','DOR'))+
  xlab('')+ylab('Noise color')+  guides(color = guide_legend(nrow = 1))+
  theme(legend.key.height = unit(0.05,'cm' ))+
  theme(legend.key.width =unit(0.02,'cm' ))+xlab('Predictor')+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank()) +theme(axis.text=element_text(size=12,family='serif'),axis.title=element_text(size=14,family='serif')) +
  theme(legend.position ='none')+  
  annotate(geom='text',x=0.02,y=0.3,label='Area',color='black',family='serif')+
  annotate(geom='text',x=0.88,y=0.38,label='Temperature',color='black',family='serif')+
  annotate(geom='text',x=0.13,y=0.28,label='Wetland',color='black',family='serif')+
  annotate(geom='text',x=0.25,y=0.24,label='Dam Regulation',color='black',family='serif')+
  annotate(geom='text',x=0.6,y=0.3,label='Urban',color='black',family='serif')

p01<-grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9))
p02<-grid.arrange(parnoise_daily,parnoise_annual,ncol=2)

grid.arrange(p01,p02)


