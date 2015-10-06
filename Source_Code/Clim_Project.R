
## function to evaluate differences in the estimated day of year for phenological events based on current climate and climate change scenarios 
## changes are split into absolute change and rate fo change for trophic levels and taxonomic groups


ClimProject <- function(TempScenarioFile,PrecipScenarioFile,Model,AllPhen){

    #load required libraries
    require(shapefiles)
    require(splancs)

    #read in shapefiles corresponding to temperature and precipitation scenarios
    TempScenario <- read.shapefile(TempScenarioFile)
    PrecipScenario <- read.shapefile(PrecipScenarioFile)
    
    #ensure columns in data file are coded correctly
    for(colm in c(1:4,13,16:23,38:39,62:65)){
    
        AllPhen[,colm] = as.numeric(as.character(AllPhen[,colm]))
    
    }
    
    #create an empty vector to store matching spatial units in
    spat_id_tm=rep(NA,length(AllPhen[,1]))
    
    #loop over every spatial unit in the scenario data - eg 25km squares
    for(k in 1:length(TempScenario$shp$shp)){
    
      #find which data in the phenology  records is in this unit and asign it a correspoinding id
      idx=which(inout(AllPhen[,39:38],TempScenario$shp$shp[[k]]$points))
      spat_id_tm[idx]=k
    
    }
    
    #repeat spatial matching for the precipitation scenario data
    spat_id_pr=rep(NA,length(AllPhen[,1]))
    
    for(k in 1:length(PrecipScenario$shp$shp)){
    
      idx=which(inout(AllPhen[,39:38],PrecipScenario$shp$shp[[k]]$points))
      spat_id_pr[idx]=k
    
    }
    
    #create empty vecotrs to store results
    slopes=troph=taxa=absval=c()
    
    cnt=1 ; v=c()
    
    
    #loop through all uniquely defined time series in the data set. 
    for(k in unique(AllPhen$UniVal)){
    
        #find all records corresponding to the same unique time series
        idx = which(AllPhen$UniVal==k)
        
        #take the row that shows the latest year in the record
        id=idx[which.max(AllPhen$Year[idx])]
        
        #extract the dates from the selected windows for each record and pull out the corresponding month
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
        
        
        #define multiplyers as the scenario data form corresponding spatial unit from the months as defined by the seasonal windows  
        mult1=mean(as.numeric(TempScenario$dbf$dbf[spat_id_tm[cnt],(mon1+5)]))
        mult2=mean(as.numeric(TempScenario$dbf$dbf[spat_id_tm[cnt],(mon2+5)]))
        mult3=mean(as.numeric(PrecipScenario$dbf$dbf[spat_id_pr[cnt],(mon3+5)]))/30.5
        mult4=mean(as.numeric(PrecipScenario$dbf$dbf[spat_id_pr[cnt],(mon4+5)]))/30.5
        
        #the baseline prediction for the latest year is obtained by multiplying the extracted climate data by the saved coefficients for each series
        pred1 = sum(c(
        AllPhen[id,1]*AllPhen[id,65] ,
        AllPhen[id,2]*AllPhen[id,64] ,
        AllPhen[id,3]*AllPhen[id,63] ,
        AllPhen[id,4]*AllPhen[id,62] ,
        AllPhen[id,13]),na.rm=TRUE)
        #back transform to get a day of year
        pr1=(1/pred1)-365
        
        #repeat prediction but using the new climate data (eg with the scenario change - 3degree increase for example - added on)
        pred2 = sum(c(
        AllPhen[id,1]*(mult1+AllPhen[id,65]) ,
        AllPhen[id,2]*(mult2+AllPhen[id,64]) ,
        AllPhen[id,3]*(mult3+AllPhen[id,63]) ,
        AllPhen[id,4]*(mult4+AllPhen[id,62]) ,
        AllPhen[id,13]),na.rm=TRUE)
        #back transform to get a day of year
        pr2=(1/pred2)-365
        
        #work out the rate of change based on the year of the latest record in the data and the climate scenario for 2050
        slopes[cnt]=(pr2-pr1)/(2050-as.numeric(as.character(AllPhen$Year[id])))
        #work out the absolute change 
        absval[cnt]=pr2-pr1
        #store the corresponding trophic level and taxonomic group for each time series.
        troph[cnt]=AllPhen$Trophic.level[id]
        taxa[cnt]=AllPhen$Taxa[id]
        
        v[cnt]=mult1
        cnt=cnt+1
    
    
    }
    
    #function to evaluate 95 percentiles and median
    plot_vals=function(x){return(c(quantile(x,0.05,na.rm=TRUE),median(x,na.rm=TRUE),quantile(x,0.95,na.rm=TRUE)))}
    #evaluate summary metrics by each grouping (trophic level and taxonomic group
    res_tr=tapply(slopes,troph,plot_vals)
    res_tx=tapply(slopes,taxa,plot_vals)

    #return the summaries for the different goruping classes 
    return(list(TaxaProj=res_tx,TrophicProj=res_tr))


}




