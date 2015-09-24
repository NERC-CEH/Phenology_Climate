Climate_Windows <- function(phendat,climdat=NULL,rand=FALSE) {
     
    require(mgcv) 

    if(is.null(climdat)){
       temp_dat=get_chess(phendat)
    }else{
       temp_dat=climdat
    }

    tm <- strptime(temp_dat$Date,"%Y-%m-%d")
    temp_dat$Year = tm$year+1900
    temp_dat$DOY = tm$yday+1
    temp_dat$Month = tm$mon+1
    temp_dat$DayMonth= tm$mday
    temp_dat$DateNonYear = paste(temp_dat$Month,temp_dat$DayMonth,sep="_")


    #################################
    #################################
    #################################

    phendat$Tmean_l = phendat$Tmean_u = phendat$Precip_l = phendat$Precip_u = NA

    Tmean_l_u <- find_window(param="Tmean",tempphen=phendat,tempclim=temp_dat,rand=rand)
    phendat$Tmean_l = Tmean_l_u$Lower[match(phendat$Year,names(Tmean_l_u$Lower))]
    phendat$Tmean_u = Tmean_l_u$Upper[match(phendat$Year,names(Tmean_l_u$Upper))]

   
    Precip_l_u <- find_window(param="Precip",tempphen=phendat,tempclim=temp_dat,rand=rand)
    phendat$Precip_l = Precip_l_u$Lower[match(phendat$Year,names(Precip_l_u$Lower))]
    phendat$Precip_u = Precip_l_u$Upper[match(phendat$Year,names(Precip_l_u$Upper))]

    mod=glm((Day.of.year+365)~Tmean_l + Tmean_u + Precip_l + Precip_u,family=Gamma,data=phendat)
    nullmod = glm((Day.of.year+365)~Year,family=Gamma,data=phendat)


    mt_idx <- match(row.names(summary(mod)$coefficients)[-1],c("Tmean_l","Tmean_u","Precip_l","Precip_u"))
    coef_est <- rep(NA,4) ; coef_est[mt_idx] <-summary(mod)$coefficients[-1,1]
    sterr_est <- rep(NA,4) ; sterr_est[mt_idx] <-summary(mod)$coefficients[-1,2]
    pval_est <- rep(NA,4) ; pval_est[mt_idx] <-summary(mod)$coefficients[-1,4]
    intercept <- summary(mod)$coefficients[1,1]
    resid_var <- var((predict(mod,type="response")-365)-(mod$y-365))
    #resid_var <- var((predict(mod,type="response")-365)-phendat$Day.of.year)

    ## add column names and things like site id, metric, species etc.
    Out_Info=c(coef_est,sterr_est,pval_est,intercept,resid_var,100*(1-(summary(mod)$deviance/summary(mod)$null.deviance)),Tmean_l_u$DOY_l[1],Tmean_l_u$DOY_l[2],Tmean_l_u$DOY_u[1],Tmean_l_u$DOY_u[2],Precip_l_u$DOY_l[1],Precip_l_u$DOY_l[2],Precip_l_u$DOY_u[1],Precip_l_u$DOY_u[2],AIC(mod)-AIC(nullmod),summary(nullmod)$coefficients[2,1],var(phendat$Day.of.year), var(phendat$Tmean_l,na.rm=TRUE), var(phendat$Tmean_u,na.rm=TRUE), var(phendat$Precip_l,na.rm=TRUE), var(phendat$Precip_u,na.rm=TRUE), c(as.matrix(phendat[1,])))

    nm_idx=match(c("Year","Day.of.year","Species.common.name","Species.latin.name", "Monitoring.scheme.name"
     ,"Site.name"
     ,"Single.site..1...or.Area..2."
     ,"Latitude.of.site.or.top.left.hand.corner.if.an.area..decimal.degrees."
     ,"Longitude.of.site.or.top.left.hand.corner.if.an.area..decimal.degrees."
     ,"Latitude.of.bottom.right.hand.corner.if.an.area..decimal.degrees."
     ,"Longitude.of.bottom.right.hand.corner.if.an.area.decimal.degrees."
     ,"Altitude..metres.above.sea.level."
     ,"Freshwater..F...Marine..M...Terrestrial..T."
     ,"Trophic.level"
     ,"Environment"
     ,"Habitat"
     ,"Habitat.specialism"
     ,"Generation.time..under.optimal.conditions."
     ,"Metric.class"
     ,"Specific.metric"
     ,"Event"
     ,"Taxonomic.Class"
     ,"Thermal.Category.1"
     ,"Thermal.Category.2"
     ,"Data.resolution"
     ,"Movement.Class"
     ,"Length.class"
     ,"Mass.class"
     ,"Voltinism"
     ,"Consumer.trophic.specialism"
     ,"UniVal"
     ,"Precip_u"
     ,"Precip_l"
    , "Tmean_u"
    ,"Tmean_l" ),names(phendat))

    return(list(Summary=Out_Info,Mod_Data=phendat[,nm_idx]))

}

find_window = function(param,tempphen,tempclim,rand=FALSE){

  tempphen<- tempphen[tempphen$Year>=1960,]

  all_dates=rep(unique(tempclim$DateNonYear),2)

  store_dat <- data.frame(Date=NA,Rsq=NA,Coef=NA)

  tempclim$modparam = tempclim[,grep(param,names(tempclim))]
  tempclim$match_id = paste(tempclim$Year,tempclim$DateNonYear,sep="_")

  ###loop over days
  for(i in unique(tempclim$DateNonYear)){

    if(i=="2_29"){tempphen1=tempphen[is.element(tempphen$Year,seq(1960,2020,by=4)),]}else{tempphen1=tempphen}

    tempphen1$match_id = paste(tempphen1$Year,i,sep="_")
    tempphen1$mean_temp <- tempclim$modparam[match(tempphen1$match_id,tempclim$match_id)]

    if(!rand){
        mod <- lm(Day.of.year~mean_temp,data=tempphen1)
        store_dat = rbind(store_dat,data.frame(Date=i,Rsq=summary(mod)$r.squared,Coef=summary(mod)$coefficients[2,1]))
    }
    else{
        if(length(tempphen1[,1])>4){
            mod <- lme(Day.of.year~mean_temp,random=list(Year=~1),data=tempphen1,na.action=na.exclude)
            store_dat = rbind(store_dat,data.frame(Date=i,Rsq=summary(lm(tempphen1$Day.of.year ~ mod$fitted[,2]))$r.squared,Coef=summary(mod)$coefficients$fixed[2]))
        }else{
            store_dat = rbind(store_dat,data.frame(Date=i,Rsq=NA,Coef=NA))
        }
    }

  }

  #make sure the previous year is captured too
  store_dat_prev <- data.frame(Date=NA,Rsq=NA,Coef=NA)

  for(i in unique(tempclim$DateNonYear)){

    tempphen$match_id = paste(tempphen$Year-1,i,sep="_")

    ### is the first year / last year a leap year or not?
    if(i=="2_29"){

        tempphen1=tempphen[is.element(tempphen$Year-1,seq(1960,2020,by=4)),]
        tempphen1$mean_temp <- c(tempclim$modparam[match(tempphen1$match_id,tempclim$match_id)])

    }

    else{

        tempphen1=tempphen
        tempphen1$mean_temp <- c(tempclim$modparam[match(tempphen1$match_id,tempclim$match_id)])

    }

    if(!rand){
        mod <- lm(Day.of.year~mean_temp,data=tempphen1)
        store_dat_prev = rbind(store_dat_prev,data.frame(Date=i,Rsq=summary(mod)$r.squared,Coef=summary(mod)$coefficients[2,1]))
    }
    else{
        if(length(tempphen1[,1])>4){
            mod <- lme(Day.of.year~mean_temp,random=list(Year=~1),data=tempphen1,na.action=na.exclude)
            store_dat_prev = rbind(store_dat_prev,data.frame(Date=i,Rsq=summary(lm(tempphen1$Day.of.year ~ mod$fitted[,2]))$r.squared,Coef=summary(mod)$coefficients$fixed[2]))
        }else{
            store_dat_prev = rbind(store_dat_prev,data.frame(Date=i,Rsq=NA,Coef=NA))
        }
    }


  }


  store_dat=store_dat[-1,]
  store_dat_prev=store_dat_prev[-1,]


  MonthName=c("January","February","March","April","May","June","July","August","September","October","November","December")



  par(mfrow=c(2,1))
    plot(c(store_dat_prev[,3],store_dat[,3]),type="l",ylab="Coefficients",xaxt="n",xlab="Month",main=paste(param," Fit :-  Species = ",unique(tempphen$Species.latin.name),"; Metric = ",unique(tempphen$Specific.metric),sep=""),cex.main=0.75)
    abline(v=c(quantile(tempphen$Day.of.year,0.95),quantile(tempphen$Day.of.year,0.95)+365),col="green",lwd=2)
    axis(1,at=c(c(1,which(diff(tempclim$Month)[1:366]!=0)[1:11]),(c(1,which(diff(tempclim$Month)[1:366]!=0)[1:11])+366)),labels=rep(MonthName,2),cex.axis=0.6,las=2)


    x=c(store_dat_prev[,3],store_dat[,3])+100
    x[x<=(0)]=0.001
    y=(1:732)
    sm_mod=gam(x~s(y,k=18,fx=TRUE),family=Gamma(link=log),method="REML",na.action=na.exclude)
    pred_C=predict(sm_mod,se.fit=TRUE,type="response")
    lines(pred_C$fit-100,col="red",lwd=3.5)
    lines((pred_C$fit+1.96*pred_C$se.fit)-100,lty=2,col="red",lwd=3)
    lines((pred_C$fit-1.96*pred_C$se.fit)-100,lty=2,col="red",lwd=3)

    if(mean(tempphen$Day.of.year)<=0){
      pred_C$fit=pred_C$fit[1:366]
      pred_C$se.fit=pred_C$se.fit[1:366]

    }else{
      pred_C$fit=pred_C$fit[ceiling(quantile(tempphen$Day.of.year,0.95)):floor(quantile(tempphen$Day.of.year,0.95)+366)]
      pred_C$se.fit=pred_C$se.fit[ceiling(quantile(tempphen$Day.of.year,0.95)):floor(quantile(tempphen$Day.of.year,0.95)+366)]
    }

    abline(h=0,col="cyan",lwd=2)
    abline(h=quantile(pred_C$fit-100,0.975),col="blue",lwd=3)
    abline(h=quantile(pred_C$fit-100,0.025),col="blue",lwd=3)

  #beta <- coef(sm_mod);Vb <- vcov(sm_mod)

  ## simulate replicate beta vectors from posterior...
  #Cv <- chol(Vb)
  #n.rep=10000;nb <- length(beta)
  #br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep) + beta

  ## turn these into replicate linear predictors...
  #xp <- 0:(length(y))#/length(y)
  #Xp <- predict(sm_mod,newdata=data.frame(y=xp),type="lpmatrix")
  #lp <- Xp%*%br
  #fv <- exp(lp) ## ... finally, replicate expected value vectors

  ## now simulate from Gamma deviates with mean as in fv
  ## and estimated scale...
  #yr <- matrix(rgamma(fv*0,shape=1/sm_mod$scale,scale=fv*sm_mod$scale),nrow(fv),ncol(fv))

  ## compute 95% prediction interval...
  #PI_C <- apply(yr,1,quantile,prob=c(.025,0.975))
  ## and plot it...
  #lines(xp,PI_C[1,]-10,col=2,lwd=2,lty=3);lines(xp,PI_C[2,]-10,col=2,lwd=2,lty=3)

    plot(c(store_dat_prev[,2],store_dat[,2]),type="l",ylab="R Squared",xaxt="n",xlab="Month")
    abline(v=c(quantile(tempphen$Day.of.year,0.95),quantile(tempphen$Day.of.year,0.95)+365),col="green",lwd=2)
    axis(1,at=c(c(1,which(diff(tempclim$Month)[1:366]!=0)[1:11]),(c(1,which(diff(tempclim$Month)[1:366]!=0)[1:11])+366)),labels=rep(MonthName,2),cex.axis=0.6,las=2)


    x=c(store_dat_prev[,2],store_dat[,2])+100
    x[x<=(0)]=0.001
    y=(1:732)
    sm_mod=gam(x~s(y,k=18,fx=TRUE),family=Gamma(link=log),method="REML",na.action=na.exclude)
    pred_R=predict(sm_mod,se.fit=TRUE,type="response")
    lines(pred_R$fit-100,col="red",lwd=3.5)
    lines((pred_R$fit+1.96*pred_R$se.fit)-100,lty=2,col="red",lwd=3)
    lines((pred_R$fit-1.96*pred_R$se.fit)-100,lty=2,col="red",lwd=3)

    if(mean(tempphen$Day.of.year)<=0){
      pred_R$fit=pred_R$fit[1:366]
      pred_R$se.fit=pred_R$se.fit[1:366]
    }else{
      pred_R$fit=pred_R$fit[ceiling(quantile(tempphen$Day.of.year,0.95)):floor(quantile(tempphen$Day.of.year,0.95)+366)]
      pred_R$se.fit=pred_R$se.fit[ceiling(quantile(tempphen$Day.of.year,0.95)):floor(quantile(tempphen$Day.of.year,0.95)+366)]
    }


    abline(h=0,col="cyan",lwd=2)
    abline(h=quantile(pred_R$fit-100,0.975),col="blue",lwd=3)
    abline(h=quantile(pred_R$fit-100,0.025),col="blue",lwd=3)

  #beta <- coef(sm_mod);Vb <- vcov(sm_mod)

  ## simulate replicate beta vectors from posterior...
  #Cv <- chol(Vb)
  #n.rep=10000;nb <- length(beta)
  #br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep) + beta

  ## turn these into replicate linear predictors...
  #xp <- 0:(length(y))#/length(y)
  #Xp <- predict(sm_mod,newdata=data.frame(y=xp),type="lpmatrix")
  #lp <- Xp%*%br
  #fv <- exp(lp) ## ... finally, replicate expected value vectors

  ## now simulate from Gamma deviates with mean as in fv
  ## and estimated scale...
  #yr <- matrix(rgamma(fv*0,shape=1/sm_mod$scale,scale=fv*sm_mod$scale),nrow(fv),ncol(fv))

  ## compute 95% prediction interval...
  #PI_R <- apply(yr,1,quantile,prob=c(.025,0.975))
  ## and plot it...
  #lines(xp,PI_R[1,]-10,col=2,lwd=2,lty=3);lines(xp,PI_R[2,]-10,col=2,lwd=2,lty=3)


  ## get window dates
  if(mean(tempphen$Day.of.year)<=0){
  poss_dates <- all_dates[1:366]
  }else{
  poss_dates <- all_dates[ceiling(quantile(tempphen$Day.of.year,0.95)):floor(quantile(tempphen$Day.of.year,0.95)+366)]
  }

  coef_range_l <- which(((pred_C$fit-1.96*pred_C$se.fit)-100)<=quantile(pred_C$fit-100,0.025))
  coef_range_u <- which(((pred_C$fit+1.96*pred_C$se.fit)-100)>=quantile(pred_C$fit-100,0.975))

  r_sq_range <- which(((pred_R$fit+1.96*pred_R$se.fit)-100)>=quantile(pred_R$fit-100,0.975))

  temp_window_l <- coef_range_l[is.element(coef_range_l,r_sq_range)]
  temp_window_u <- coef_range_u[is.element(coef_range_u,r_sq_range)]

  if(length(temp_window_l)==0){temp_window_l = coef_range_l}
  if(length(temp_window_u)==0){temp_window_u = coef_range_u}


  temp_window_l = find_concurrent_period(temp_window_l,pred_C=pred_C)
  temp_window_u = find_concurrent_period(temp_window_u,pred_C=pred_C)

  tmp_wind_l <- poss_dates[temp_window_l]
  tmp_wind_u <- poss_dates[temp_window_u]

  tmp_slim_l <- tempclim[is.element(tempclim$DateNonYear,tmp_wind_l),]
  tmp_slim_u <- tempclim[is.element(tempclim$DateNonYear,tmp_wind_u),]

  doysl <- tempclim$DOY[match(tmp_wind_l,tempclim$DateNonYear)]
  doysu <- tempclim$DOY[match(tmp_wind_u,tempclim$DateNonYear)]

  if(length(which(doysl>quantile(tempphen$Day.of.year,0.95)))>0){
    idx=which(doysl>quantile(tempphen$Day.of.year,0.95))
    tmp_slim_l$Year[is.element(tmp_slim_l$DateNonYear,tmp_wind_l[idx])]=tmp_slim_l$Year[is.element(tmp_slim_l$DateNonYear,tmp_wind_l[idx])]+1
  }

  if(length(which(doysu>quantile(tempphen$Day.of.year,0.95)))>0){
    idx=which(doysu>quantile(tempphen$Day.of.year,0.95))
    tmp_slim_u$Year[is.element(tmp_slim_u$DateNonYear,tmp_wind_u[idx])]=tmp_slim_u$Year[is.element(tmp_slim_u$DateNonYear,tmp_wind_u[idx])]+1
  }

  mn_temps_l <- tapply(tmp_slim_l$modparam,tmp_slim_l$Year,mean)
  mn_temps_u <- tapply(tmp_slim_u$modparam,tmp_slim_u$Year,mean)

  Temp_matched_l=mn_temps_l[match(tempphen$Year,names(mn_temps_l))]
  Temp_matched_u=mn_temps_u[match(tempphen$Year,names(mn_temps_u))]


  return(list(Lower=Temp_matched_l,Upper=Temp_matched_u,DOY_l=doysl[c(1,length(doysl))],DOY_u=doysu[c(1,length(doysu))]))

}


##############################################
##############################################


find_concurrent_period=function(temp_window,pred_C){

    idx=1
    diff_id = c(1,which(diff(temp_window)>1))

    if(length(diff_id)>1){if((which(diff(temp_window)>1))[1]==1){diff_id=diff_id[-1]}}
    if(length(temp_window)>1){
    for(rv in 1:length(diff_id)){

      if(rv==length(diff_id)){
        idx=c(idx,rep(rv,length((diff_id[rv]+1):length(temp_window))))
      }else{
        idx=c(idx,rep(rv,length((diff_id[rv]+1):(diff_id[rv+1]))))
      }

    }
    }

    mx_coef = which.max(tapply(abs(pred_C$fit[temp_window]),idx,mean))

    return(temp_window[idx==mx_coef])



}

get_chess <- function(phenseries){

    patheast=paste("http://thredds.nerc-lancaster.ac.uk/thredds/dodsC/chess/driving_data/aggregation/precip_aggregation.ncml.ascii?x[",0,":1:",655,"]",sep="")
    east=read.table(patheast,skip=5,header=FALSE,sep=",")[1,]
    pathnorth=paste("http://thredds.nerc-lancaster.ac.uk/thredds/dodsC/chess/driving_data/aggregation/precip_aggregation.ncml.ascii?y[",0,":1:",1056,"]",sep="")
    north=read.table(pathnorth,skip=5,header=FALSE,sep=",")[1,]

    lonid=which(names(phenseries)=="Longitude.of.site.or.top.left.hand.corner.if.an.area..decimal.degrees.")
    latid=which(names(phenseries)=="Latitude.of.site.or.top.left.hand.corner.if.an.area..decimal.degrees.")
    coords <- LatLon_2_OS(phenseries[1,lonid],phenseries[1,latid])

    wy=(which.min((north-coordinates(coords)[2])^2)-1)
    wx=(which.min((east-coordinates(coords)[1])^2)-1)

    wt1=0
    wt2=19357
    wt3=19357

    patht=paste("http://thredds.nerc-lancaster.ac.uk/thredds/dodsC/chess/driving_data/aggregation/tas_aggregation.ncml.ascii?tas[",wt1,":1:",wt2,"][",wy,":1:",wy,"][",wx,":1:",wx,"]",sep="")
    pathp=paste("http://thredds.nerc-lancaster.ac.uk/thredds/dodsC/chess/driving_data/aggregation/precip_aggregation.ncml.ascii?precip[",wt1,":1:",wt3,"][",wy,":1:",wy,"][",wx,":1:",wx,"]",sep="")

    xt=try(read.table(patht,skip=11,header=TRUE,nrows=(1+(wt2-wt1))))
    if(class(xt)=="try-error"){xt=try(read.table(patht,skip=11,header=TRUE,nrows=(1+(wt2-wt1))))}

    xp=try(read.table(pathp,skip=11,header=TRUE,nrows=(1+(wt3-wt1))))
    if(class(xp)=="try-error"){xp=try(read.table(pathp,skip=11,header=TRUE,nrows=(1+(wt3-wt1))))}

    yt=try(read.table(patht,skip=11+(wt2-wt1)+4,header=FALSE,nrows=1,sep=","))
    chess_dates = unlist(strsplit(as.character(unlist(strsplit(as.character(as.POSIXct((min(yt):max(yt))*60*60*24,origin="1961-01-01")[1:length(yt)]),"GMT"))),"[ ]"))[seq(1,length(yt)*2,by=2)]

    temp_dat <- data.frame(Date=chess_dates,Tmean=xt[,1]-273,Precip=xp[,1])

    return(temp_dat)

}



OS_2_LatLon <- function(xloc,yloc){

  require(rgdal);require(sp)
  pt = data.frame(x=xloc,y=yloc)
  coordinates(pt)=~x+y

  proj4string(pt)=CRS("+init=epsg:27700")

  return(spTransform(pt,CRS("+init=epsg:4326")))

}


LatLon_2_OS <- function(xloc,yloc){

  require(rgdal);require(sp)
  pt = data.frame(x=xloc,y=yloc)
  coordinates(pt)=~x+y

  proj4string(pt)=CRS("+init=epsg:4326")

  return(spTransform(pt,CRS("+init=epsg:27700")))

}

split_data <- function(phen_in){

  phen_in$UniVal <- paste(phen_in$Monitoring.scheme.name,phen_in$Site.name,phen_in$Species.latin.name,phen_in$Specific.metric,phen_in$Voltinism,sep="_")

  if(!file.exists(paste(getwd(),"//Split_Files",sep=""))){
     dir.create(file.path(getwd(),"Split_Files//"))
  }
  
  for(k in 1:length(unique(phen_in$UniVal))){

     write.csv(phen_in[phen_in$UniVal==((unique(phen_in$UniVal))[k]),],file=paste("Split_Files//phen_in",k,".csv",sep="_"),row.names=FALSE)

  }

}



##read in all data in output folder

run_model <- function(allphen,Interaction="None"){


allphen$Taxa="Other"
allphen$Taxa[grep("algal",as.character(allphen$Taxonomic.Class.corr))]="Algae"
allphen$Taxa[grep("insect",as.character(allphen$Taxonomic.Class.corr))]="Insect"
allphen$Taxa[grep("fish",as.character(allphen$Taxonomic.Class.corr))]="Fish"
allphen$Taxa[grep("plants",as.character(allphen$Taxonomic.Class.corr))]="Plants"
allphen$Taxa[grep("conifers",as.character(allphen$Taxonomic.Class.corr))]="Plants"
allphen$Taxa[grep("birds",as.character(allphen$Taxonomic.Class.corr))]="Birds"
allphen$Taxa[grep("amphibians",as.character(allphen$Taxonomic.Class.corr))]="Amphibians"
allphen$Taxa[grep("algae",as.character(allphen$Taxonomic.Class.corr))]="Algae"
allphen$Taxa[grep("cladocera",as.character(allphen$Taxonomic.Class.corr))]="Crustacea"
allphen$Taxa[grep("molluscs",as.character(allphen$Taxonomic.Class.corr))]="Molluscs"
allphen$Taxa[grep("barnacles",as.character(allphen$Taxonomic.Class.corr))]="Crustacea"
allphen$Taxa[grep("lobsters",as.character(allphen$Taxonomic.Class.corr))]="Crustacea"
allphen$Taxa[grep("Mammal",as.character(allphen$Taxonomic.Class.corr))]="Mammals"


allphen$Taxa[allphen$Species.latin.name=="Chlorophyll a" & allphen$Taxonomic.Class.corr=="Other"]="Algae"
allphen$Taxa[allphen$Species.latin.name=="Chlorella" & allphen$Taxonomic.Class.corr=="Other"]="Algae"


phendates = tapply(allphen$Day.of.year,allphen$UniVal,median,na.rm=TRUE)

idx=match(allphen$UniVal,names(phendates))
allphen$MeanDOY=NA
allphen$MeanDOY=phendates[idx]

allphen$EarlyLate="Early"
allphen$EarlyLate[allphen$MeanDOY>181]="Late"

allphen$PhenSeas="WinterSpring"
allphen$PhenSeas[allphen$MeanDOY>152 & allphen$MeanDOY<245]="Summer"
allphen$PhenSeas[allphen$MeanDOY>244]="AutumnWinter"



##run model
library(lme4)

allphen$VoltMetric=paste(allphen$Voltinism,allphen$Metric.class,sep="")

switch(Interaction,
None={form=as.formula(Day.of.year~Tmean_l+Tmean_u+Precip_l+Precip_u +(1+Tmean_l+Tmean_u+Precip_l+Precip_u|Species.latin.name/VoltMetric))},
Taxa={form=as.formula(Day.of.year~Taxonomic.Class*Tmean_l+Taxonomic.Class*Tmean_u+Taxonomic.Class*Precip_l+Taxonomic.Class*Precip_u +(1+Tmean_l+Tmean_u+Precip_l+Precip_u|Species.latin.name/VoltMetric))},
Environment={form=as.formula(Day.of.year~Environment*Tmean_l+Environment*Tmean_u+Environment*Precip_l+Environment*Precip_u +(1+Tmean_l+Tmean_u+Precip_l+Precip_u|Species.latin.name/VoltMetric))},
TrophicLevel={form=as.formula(Day.of.year~Trophic.level*Tmean_l+Trophic.level*Tmean_u+Trophic.level*Precip_l+Trophic.level*Precip_u +(1+Tmean_l+Tmean_u+Precip_l+Precip_u|Species.latin.name/VoltMetric))}
)

mod=lmer(form,data=allphen)


###save model
save(mod,file="Model_Out")

return(list(OutMod=mod,ModData=allphen))

#save(allphen,mod,file="Model_Out.Rdata")

}



##########
##########


plot_mod_res <- function(mod,type="Dens"){
library(lme4)
if(type=="Dens"){
par(mfrow=c(1,2))

vars=c("Tmean_l","Tmean_u")
coli = c("black","red","blue")

estt=list();est=c()
for(tx_count in 1:2){


  k=grep(vars[tx_count],names(fixef(mod)))
  est[tx_count]=(fixef(mod)[k])
  estt[[tx_count]]=((fixef(mod)[k]+ranef(mod)[[2]][,tx_count+1]))

}

try(plot(density(estt[[1]]),xlim=range(c(density(estt[[1]])$x,density(estt[[2]])$x)),main="",xlab=expression("Days per "^o*"C"),col=coli[1],ylim=range(c(density(estt[[1]])$y,density(estt[[2]])$y))))
try(lines(density(estt[[2]]),col=coli[2]))

abline(v=est[1],lwd=2,col=coli[1])
abline(v=est[2],lwd=2,col=coli[2])
abline(v=0,lwd=2,lty=2,col="grey")

legend("topright",legend=as.character(vars),text.col = coli,bty="n",lty=1,col=coli,lwd=2)




vars=c("Precip_l","Precip_u")
coli = c("black","red","blue")

estt=list();est=c()
for(tx_count in 1:2){


  k=grep(vars[tx_count],names(fixef(mod)))
  est[tx_count]=(fixef(mod)[k])
  estt[[tx_count]]=((fixef(mod)[k]+ranef(mod)[[2]][,tx_count+3]))

}

try(plot(density(estt[[1]]),xlim=range(c(density(estt[[1]])$x,density(estt[[2]])$x)),main="",xlab="Days per mm",col=coli[1],ylim=range(c(density(estt[[1]])$y,density(estt[[2]])$y))))
try(lines(density(estt[[2]]),col=coli[2]))

abline(v=est[1],lwd=2,col=coli[1])
abline(v=est[2],lwd=2,col=coli[2])
abline(v=0,lwd=2,lty=2,col="grey")

legend("topright",legend=as.character(vars),text.col = coli,bty="n",lty=1,col=coli,lwd=2)



}else{

  if(type=="Slopes"){
  
    par(mfrow=c(2,2))
  
    xc=seq(0,12,len=100)
    xc=seq(0,30,len=100)
    
    coli = c("black","red","blue")
    
    k=grep("Tmean_l",names(fixef(mod)))
    
    plot(c(0,12),c(0,365),ylim=c(0,365),type="n",col="red",lwd=2,xlab="Temperature (Degrees C)",ylab="Day of Year",main="Temperature Lower")
    
    est=fixef(mod)[1]+xc*(fixef(mod)[k])
    
    lines(xc,(est),lwd=2,lty=1)
    
    for (i in 1:length(ranef(mod)[[2]][,1])){
    
      estt=fixef(mod)[1]+ranef(mod)[[2]][i,1] + (((fixef(mod)[k]+ranef(mod)[[2]][i,2]))*xc)
      lines(xc,(estt),col="grey",lwd=1,lty=1)
    
    }
    
    lines(xc,(est),lwd=2,lty=1)
    
    k=grep("Tmean_u",names(fixef(mod)))
    
    plot(c(0,12),c(0,365),ylim=c(0,365),type="n",col="red",lwd=2,xlab="Temperature (Degrees C)",ylab="Day of Year",main="Temperature Upper")
    
    est=fixef(mod)[1]+xc*(fixef(mod)[k])
    
    lines(xc,(est),lwd=2,lty=1)
    
    for (i in 1:length(ranef(mod)[[2]][,1])){
    
      estt=fixef(mod)[1]+ranef(mod)[[2]][i,1] + (((fixef(mod)[k]+ranef(mod)[[2]][i,3]))*xc)
      lines(xc,(estt),col="grey",lwd=1,lty=1)
    
    }
    
    lines(xc,(est),lwd=2,lty=1)
    xc=seq(0,30,len=100)
    coli = c("black","red","blue")
    
    k=grep("Precip_l",names(fixef(mod)))
    
    plot(c(0,30),c(0,365),ylim=c(0,365),type="n",col="red",lwd=2,xlab="Precipitation (mm)",ylab="Day of Year",main="Precipitation Lower")
    
    est=fixef(mod)[1]+xc*(fixef(mod)[k])
    
    lines(xc,(est),lwd=2,lty=1)
    
    for (i in 1:length(ranef(mod)[[2]][,1])){
    
      estt=fixef(mod)[1]+ranef(mod)[[2]][i,1] + (((fixef(mod)[k]+ranef(mod)[[2]][i,4]))*xc)
      lines(xc,(estt),col="grey",lwd=1,lty=1)
    
    }
    
    lines(xc,(est),lwd=2,lty=1)
  
    k=grep("Precip_u",names(fixef(mod)))
    
    plot(c(0,30),c(0,365),ylim=c(0,365),type="n",col="red",lwd=2,xlab="Precipitation (mm)",ylab="Day of Year",main="Precipitation Upper")
    
    est=fixef(mod)[1]+xc*(fixef(mod)[k])
    
    lines(xc,(est),lwd=2,lty=1)
    
    for (i in 1:length(ranef(mod)[[2]][,1])){
    
      estt=fixef(mod)[1]+ranef(mod)[[2]][i,1] + (((fixef(mod)[k]+ranef(mod)[[2]][i,5]))*xc)
      lines(xc,(estt),col="grey",lwd=1,lty=1)
    
    }
    
    lines(xc,(est),lwd=2,lty=1)
  }else(print("Type must be one of Dens or Slopes"))
}


}
####################
####################


plot_output <- function(mod,plot_type="Slopes",Interaction="None"){

  if(Interaction=="None"){
  
      plot_mod_res(mod=mod,type=plot_type)
  
  }else{
  
      allphen=read.csv("AllphenData_from_mod.csv")
  
      switch(Interaction,
      "TrophicLevel"={colid=14},
      "Taxa"={colid=36},
      "Environment"={colid=15}
      )

      levs=unique(allphen[,colid])
      
      coli = rainbow(length(levs))
      
      
      par(mfrow=c(2,2))
      
      #############
      #############
      
      tx_count=1
      k=grep("Tmean_l",names(fixef(mod)))
      est=c()
      all_estt=list()
      
      for(tx in levs){
      
      less_ests=c()
      less_0ests=c()
      
      
      tx_id=list()
      slopes=c()
      if(tx_count==1){
      est[tx_count]=(fixef(mod)[k[tx_count]])
      }else{
      est[tx_count]=(fixef(mod)[k[tx_count]] + fixef(mod)[k[1]])
      }
      estt=c()
      count=1
      for (sp in unique(as.character(allphen[as.character(allphen[,colid])==tx,4]))){
      
        spl=strsplit(as.character(sp),"\\(")[[1]]
        if(length(spl)==1){
        tx_id[[count]]=grep(spl,row.names(ranef(mod)[[2]]))
        }else{
        sp1=strsplit(as.character(sp),"\\(")[[1]][1]
        sp2=strsplit(as.character(sp),"\\)")[[1]][2]
        tx_id[[count]] = intersect(grep(sp1,row.names(ranef(mod)[[2]])),grep(sp2,row.names(ranef(mod)[[2]])))
        }
        if(tx_count==1){
        estt=c(estt, ((fixef(mod)[k[tx_count]]+ranef(mod)[[2]][tx_id[[count]],2])))
        }else{
        estt=c(estt, ((fixef(mod)[k[tx_count]]+fixef(mod)[k[1]]+ranef(mod)[[2]][tx_id[[count]],2])))
        }
        count=count+1
      
      }
      #if(length(estt[!is.na(estt)])>2){
      estt=na.omit(estt)
      
      all_estt[[tx_count]]=estt
      tx_count=tx_count+1
      
      }
      
      
      for(tx_count in 1:length(levs)){
      
      if(tx_count==1){
      try(plot(density(all_estt[[tx_count]]),xlim=c(-10,5),main="",xlab=expression("Days per "^o*"C"),col=coli[tx_count],ylim=c(0,1)))
      #try(plot(density(all_estt[[tx_count]]),xlim=c(-15,10),main="",xlab="Days per mm",col=coli[tx_count],ylim=c(0,1.5)))
      
      }else{
      try(lines(density(all_estt[[tx_count]]),col=coli[tx_count]))
      }
      
      
      abline(v=est[tx_count],lwd=2,col=coli[tx_count])
      abline(v=0,lwd=2,lty=2,col="grey")
      
      
      }
      
      legend("topright",legend=as.character(levs),text.col = coli,bty="n",lty=1,col=coli,lwd=2)
      
      #############
      #############
      
      
      #############
      #############
      
      tx_count=1
      k=grep("Tmean_u",names(fixef(mod)))
      est=c()
      all_estt=list()
      
      for(tx in levs){
      
      less_ests=c()
      less_0ests=c()
      
      
      tx_id=list()
      slopes=c()
      if(tx_count==1){
      est[tx_count]=(fixef(mod)[k[tx_count]])
      }else{
      est[tx_count]=(fixef(mod)[k[tx_count]] + fixef(mod)[k[1]])
      }
      estt=c()
      count=1
      for (sp in unique(as.character(allphen[as.character(allphen[,colid])==tx,4]))){
      
        spl=strsplit(as.character(sp),"\\(")[[1]]
        if(length(spl)==1){
        tx_id[[count]]=grep(spl,row.names(ranef(mod)[[2]]))
        }else{
        sp1=strsplit(as.character(sp),"\\(")[[1]][1]
        sp2=strsplit(as.character(sp),"\\)")[[1]][2]
        tx_id[[count]] = intersect(grep(sp1,row.names(ranef(mod)[[2]])),grep(sp2,row.names(ranef(mod)[[2]])))
        }
        if(tx_count==1){
        estt=c(estt, ((fixef(mod)[k[tx_count]]+ranef(mod)[[2]][tx_id[[count]],3])))
        }else{
        estt=c(estt, ((fixef(mod)[k[tx_count]]+fixef(mod)[k[1]]+ranef(mod)[[2]][tx_id[[count]],3])))
        }
        count=count+1
      
      }
      #if(length(estt[!is.na(estt)])>2){
      estt=na.omit(estt)
      
      all_estt[[tx_count]]=estt
      tx_count=tx_count+1
      
      }
      
      
      for(tx_count in 1:length(levs)){
      
      if(tx_count==1){
      try(plot(density(all_estt[[tx_count]]),xlim=c(-10,5),main="",xlab=expression("Days per "^o*"C"),col=coli[tx_count],ylim=c(0,1)))
      #try(plot(density(all_estt[[tx_count]]),xlim=c(-15,10),main="",xlab="Days per mm",col=coli[tx_count],ylim=c(0,1.5)))
      
      }else{
      try(lines(density(all_estt[[tx_count]]),col=coli[tx_count]))
      }
      
      
      abline(v=est[tx_count],lwd=2,col=coli[tx_count])
      abline(v=0,lwd=2,lty=2,col="grey")
      
      
      }
      
      legend("topright",legend=as.character(levs),text.col = coli,bty="n",lty=1,col=coli,lwd=2)
      
      #############
      #############
      
      
      
      #############
      #############
      
      tx_count=1
      k=grep("Precip_l",names(fixef(mod)))
      est=c()
      all_estt=list()
      
      for(tx in levs){
      
      less_ests=c()
      less_0ests=c()
      
      
      tx_id=list()
      slopes=c()
      if(tx_count==1){
      est[tx_count]=(fixef(mod)[k[tx_count]])
      }else{
      est[tx_count]=(fixef(mod)[k[tx_count]] + fixef(mod)[k[1]])
      }
      estt=c()
      count=1
      for (sp in unique(as.character(allphen[as.character(allphen[,colid])==tx,4]))){
      
        spl=strsplit(as.character(sp),"\\(")[[1]]
        if(length(spl)==1){
        tx_id[[count]]=grep(spl,row.names(ranef(mod)[[2]]))
        }else{
        sp1=strsplit(as.character(sp),"\\(")[[1]][1]
        sp2=strsplit(as.character(sp),"\\)")[[1]][2]
        tx_id[[count]] = intersect(grep(sp1,row.names(ranef(mod)[[2]])),grep(sp2,row.names(ranef(mod)[[2]])))
        }
        if(tx_count==1){
        estt=c(estt, ((fixef(mod)[k[tx_count]]+ranef(mod)[[2]][tx_id[[count]],4])))
        }else{
        estt=c(estt, ((fixef(mod)[k[tx_count]]+fixef(mod)[k[1]]+ranef(mod)[[2]][tx_id[[count]],4])))
        }
        count=count+1
      
      }
      #if(length(estt[!is.na(estt)])>2){
      estt=na.omit(estt)
      
      all_estt[[tx_count]]=estt
      tx_count=tx_count+1
      
      }
      
      
      for(tx_count in 1:length(levs)){
      
      if(tx_count==1){
      #try(plot(density(estt),xlim=c(-15,10),main="",xlab=expression("Days per "^o*"C"),col=coli[tx_count],ylim=c(0,1)))
      try(plot(density(all_estt[[tx_count]]),xlim=c(-4,6),main="",xlab="Days per mm",col=coli[tx_count],ylim=c(0,1.5)))
      
      }else{
      try(lines(density(all_estt[[tx_count]]),col=coli[tx_count]))
      }
      
      
      abline(v=est[tx_count],lwd=2,col=coli[tx_count])
      abline(v=0,lwd=2,lty=2,col="grey")
      
      
      }
      
      legend("topright",legend=as.character(levs),text.col = coli,bty="n",lty=1,col=coli,lwd=2)
      
      #############
      #############
      
      
      #############
      #############
      
      tx_count=1
      k=grep("Precip_u",names(fixef(mod)))
      est=c()
      
      all_estt=list()
      
      for(tx in levs){
      
      less_ests=c()
      less_0ests=c()
      
      
      tx_id=list()
      slopes=c()
      if(tx_count==1){
      est[tx_count]=(fixef(mod)[k[tx_count]])
      }else{
      est[tx_count]=(fixef(mod)[k[tx_count]] + fixef(mod)[k[1]])
      }
      estt=c()
      count=1
      for (sp in unique(as.character(allphen[as.character(allphen[,colid])==tx,4]))){
      
        spl=strsplit(as.character(sp),"\\(")[[1]]
        if(length(spl)==1){
        tx_id[[count]]=grep(spl,row.names(ranef(mod)[[2]]))
        }else{
        sp1=strsplit(as.character(sp),"\\(")[[1]][1]
        sp2=strsplit(as.character(sp),"\\)")[[1]][2]
        tx_id[[count]] = intersect(grep(sp1,row.names(ranef(mod)[[2]])),grep(sp2,row.names(ranef(mod)[[2]])))
        }
        if(tx_count==1){
        estt=c(estt, ((fixef(mod)[k[tx_count]]+ranef(mod)[[2]][tx_id[[count]],5])))
        }else{
        estt=c(estt, ((fixef(mod)[k[tx_count]]+fixef(mod)[k[1]]+ranef(mod)[[2]][tx_id[[count]],5])))
        }
        count=count+1
      
      }
      #if(length(estt[!is.na(estt)])>2){
      estt=na.omit(estt)
      
      all_estt[[tx_count]]=estt
      tx_count=tx_count+1
      
      }
      
      
      for(tx_count in 1:length(levs)){
      
      if(tx_count==1){
      #try(plot(density(estt),xlim=c(-15,10),main="",xlab=expression("Days per "^o*"C"),col=coli[tx_count],ylim=c(0,1)))
      try(plot(density(all_estt[[tx_count]]),xlim=c(-4,6),main="",xlab="Days per mm",col=coli[tx_count],ylim=c(0,1.5)))
      
      }else{
      try(lines(density(all_estt[[tx_count]]),col=coli[tx_count]))
      }
      
      
      abline(v=est[tx_count],lwd=2,col=coli[tx_count])
      abline(v=0,lwd=2,lty=2,col="grey")
      
      
      }
      
      legend("topright",legend=as.character(levs),text.col = coli,bty="n",lty=1,col=coli,lwd=2)
      
      #############
      #############
  
  }


}





#####################
#####################


mismatch_run <- function(phen_in,clim_in,rand=FALSE,plot_out=FALSE,model=FALSE,model_interaction="None",plot_type="Slopes"){

split_data(phen_in)


if(!file.exists(paste(getwd(),"//IndivOut",sep=""))){
     dir.create(file.path(getwd(),"IndivOut//"))
}


fls=list.files("Split_Files//")


pdf(file=paste("IndivOut//OutputPlots",sep="_"))

   for(k in 1:length(fls)){

      phen_indiv=read.csv(paste("Split_Files//",fls[k],sep=""))

      Out_Info <- Climate_Windows(phendat=phen_indiv,climdat=clim_in,rand=rand)

      write.csv(Out_Info$Mod_Data,file=paste("IndivOut//OutputData1_",k,sep=""),row.names=FALSE)
      write.csv(Out_Info$Summary,file=paste("IndivOut//OutputData2_",k,sep=""),row.names=FALSE)

   }

dev.off()


if(model){

fls <- list.files("IndivOut//")
idx <- grep("OutputData1",fls)

## collate all data

for(j in 1:length(idx)){

  curr_phen <- read.csv(paste("IndivOut//",fls[idx[j]],sep=""))

  if(j==1){
     allphen <- curr_phen
  }else{
     allphen <- rbind(allphen,curr_phen)
  }

}

mod_out <- run_model(allphen=allphen,Interaction=model_interaction)

}

fls <- list.files("IndivOut//")
idx <- grep("OutputData2",fls)

## collate all data

for(j in 1:length(idx)){

  curr_phen <- read.csv(paste("IndivOut//",fls[idx[j]],sep=""))

  if(j==1){
     allphen <- curr_phen
  }else{
     allphen <- cbind(allphen,curr_phen)
  }

}
allphen=t(allphen)


if(model & plot_out){

   plot_output(mod=mod_out$OutMod,plot_type=plot_type,Interaction=model_interaction)

}

return(mod_out)


}











