ClimProject <- function(TempScenarioFile,PrecipScenarioFile,Model,AllPhen){

    require(shapefiles)
    require(splancs)

    TempScenario <- read.shapefile(TempScenarioFile)
    PrecipScenario <- read.shapefile(PrecipScenarioFile)
    
    
    for(colm in c(1:4,13,16:23,38:39,62:65)){
    
    AllPhen[,colm] = as.numeric(as.character(AllPhen[,colm]))
    
    }
    
    spat_id_tm=rep(NA,length(AllPhen[,1]))
    
    for(k in 1:length(TempScenario$shp$shp)){
    
      idx=which(inout(AllPhen[,39:38],TempScenario$shp$shp[[k]]$points))
      spat_id_tm[idx]=k
    
    }
    
    spat_id_pr=rep(NA,length(AllPhen[,1]))
    
    for(k in 1:length(PrecipScenario$shp$shp)){
    
      idx=which(inout(AllPhen[,39:38],PrecipScenario$shp$shp[[k]]$points))
      spat_id_pr[idx]=k
    
    }
    
    
    slopes=troph=taxa=absval=c()
    
    cnt=1 ; v=c()
    
    for(k in unique(AllPhen$UniVal)){
    
        idx = which(AllPhen$UniVal==k)
        
        id=idx[which.max(AllPhen$Year[idx])]
        
        
        dts=as.character(unlist(strsplit(as.character(as.POSIXct(as.integer(AllPhen[id,16]:AllPhen[id,17])*60*60*24,origin="2012-01-01")[1:length(AllPhen[id,16]:AllPhen[id,17])]),"GMT")))
        mon1=(strptime(dts,"%Y-%m-%d")$mon+1)
        
        dts=as.character(unlist(strsplit(as.character(as.POSIXct(as.integer(AllPhen[id,18]:AllPhen[id,19])*60*60*24,origin="2012-01-01")[1:length(AllPhen[id,18]:AllPhen[id,19])]),"GMT")))
        mon2=(strptime(dts,"%Y-%m-%d")$mon+1)
        
        if(AllPhen$Monitoring.scheme.name[id]!="SAHFOS"){
        
            dts=as.character(unlist(strsplit(as.character(as.POSIXct(as.integer(AllPhen[id,20]:AllPhen[id,21])*60*60*24,origin="2012-01-01")[1:length(AllPhen[id,20]:AllPhen[id,21])]),"GMT")))
            mon3=(strptime(dts,"%Y-%m-%d")$mon+1)
            
            dts=as.character(unlist(strsplit(as.character(as.POSIXct(as.integer(AllPhen[id,22]:AllPhen[id,23])*60*60*24,origin="2012-01-01")[1:length(AllPhen[id,22]:AllPhen[id,23])]),"GMT")))
            mon4=(strptime(dts,"%Y-%m-%d")$mon+1)
        
        }
        
        
        
        mult1=mean(as.numeric(TempScenario$dbf$dbf[spat_id_tm[cnt],(mon1+5)]))
        mult2=mean(as.numeric(TempScenario$dbf$dbf[spat_id_tm[cnt],(mon2+5)]))
        mult3=mean(as.numeric(PrecipScenario$dbf$dbf[spat_id_pr[cnt],(mon3+5)]))/30.5
        mult4=mean(as.numeric(PrecipScenario$dbf$dbf[spat_id_pr[cnt],(mon4+5)]))/30.5
        
        
        pred1 = sum(c(
        AllPhen[id,1]*AllPhen[id,65] ,
        AllPhen[id,2]*AllPhen[id,64] ,
        AllPhen[id,3]*AllPhen[id,63] ,
        AllPhen[id,4]*AllPhen[id,62] ,
        AllPhen[id,13]),na.rm=TRUE)
        
        pr1=(1/pred1)-365
        
        pred2 = sum(c(
        AllPhen[id,1]*(mult1+AllPhen[id,65]) ,
        AllPhen[id,2]*(mult2+AllPhen[id,64]) ,
        AllPhen[id,3]*(mult3+AllPhen[id,63]) ,
        AllPhen[id,4]*(mult4+AllPhen[id,62]) ,
        AllPhen[id,13]),na.rm=TRUE)
        
        pr2=(1/pred2)-365
        
        slopes[cnt]=(pr2-pr1)/(2050-as.numeric(as.character(AllPhen$Year[id])))
        absval[cnt]=pr2-pr1
        troph[cnt]=AllPhen$Trophic.level[id]
        taxa[cnt]=AllPhen$Taxa[id]
        
        v[cnt]=mult1
        cnt=cnt+1
    
    
    }
    
    plot_vals=function(x){return(c(quantile(x,0.05,na.rm=TRUE),median(x,na.rm=TRUE),quantile(x,0.95,na.rm=TRUE)))}
    res_tr=tapply(slopes,troph,plot_vals)
    res_tx=tapply(slopes,taxa,plot_vals)


    return(list(TaxaProj=res_tx,TrophicProj=res_tr))


}




