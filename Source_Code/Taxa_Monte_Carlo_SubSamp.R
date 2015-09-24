
colnm="Trophic.level.corr"
colid=which(names(allphen)==colnm)
levs=unique(allphen[,colid])
  
 
    
idx=list()      

for(isim in 1:100){

idx[[isim]]=999
      
for(tx in unique(allphen$Taxa)){

idx[[isim]]=c(idx[[isim]],which(is.element(allphen$UniVal,sample(unique(allphen$UniVal[allphen$Taxa==tx]),100,replace=TRUE))))


}      

idx[[isim]]=idx[[isim]][-1]
      
}





MCmod <- function(ik){

    
    plot_vals=function(x){return(c(quantile(x,0.05,na.rm=TRUE),median(x,na.rm=TRUE),quantile(x,0.95,na.rm=TRUE)))}
    
    
    mod=lmer(Day.of.year~Trophic.level.corr*Tmean_l+Trophic.level.corr*Tmean_u+Trophic.level.corr*Precip_l+Trophic.level.corr*Precip_u + (1+Tmean_l+Tmean_u+Precip_l+Precip_u|Site.name/Species.latin.name/VoltMetric),data=allphen[idx[[ik]],])
    
    
    
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
    for (sp in unique(allphen$Species.latin.name[allphen[,colid]==tx])){
    
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
    
    tl=c()
    
    for(tx_count in 1:length(levs)){
    
    tl=c(tl,plot_vals(all_estt[[tx_count]]))
    
    }
    
    
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
    for (sp in unique(allphen$Species.latin.name[allphen[,colid]==tx])){
    
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
    
    
    pl=c()
    
    for(tx_count in 1:length(levs)){
    
    pl=c(pl,plot_vals(all_estt[[tx_count]]))
    
    }
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
    for (sp in unique(allphen$Species.latin.name[allphen[,colid]==tx])){
    
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
    
    
    tu=c()
    
    for(tx_count in 1:length(levs)){
    
    tu=c(tu,plot_vals(all_estt[[tx_count]]))
    
    }
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
    for (sp in unique(allphen$Species.latin.name[allphen[,colid]==tx])){
    
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
    
    
    pu=c()
    
    for(tx_count in 1:length(levs)){
    
    pu=c(pu,plot_vals(all_estt[[tx_count]]))
    
    }
    
    return(c(tl,tu,pl,pu))

}



##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################


library(snowfall)

sfInit(parallel=TRUE, cpus=16, type="SOCK")

sfLibrary(lme4)

nsim <- 100

sfExport("allphen", "levs","idx","colid")

res <- sfClusterApplyLB(1:nsim, MCmod)         

sfStop()

save(res,file="N://MC100_res.RData")





