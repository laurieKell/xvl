rfs=c("k",               "virgin",           "virgin.",          "r0",    
      "ssbtar",          "sprtar",           "ftar",             "biomasstar", 
      "ssbsprtar",       "fsprtar",          "catchsprtar",      "bmsy", 
      "sprmsy",          "fmsy",             "msy",                "retainmsy") 
names(rfs)=
  c("SSB_Unfished",    "TotBio_Unfished",  "SmryBio_Unfished", "Recr_Unfished", 
    "SSB_Btgt",        "SPR_Btgt",         "Fstd_Btgt",        "TotYield_Btgt", 
    "SSB_SPRtgt",      "Fstd_SPRtgt",      "TotYield_SPRtgt",  "SSB_MSY",  
    "SPR_MSY",         "Fstd_MSY",         "TotYield_MSY",     "RetYield_MSY") 

getTs<-function(x) {
  ts=merge(x$timeseries[,c(1:5,7:8)],ddply(x$catch,.(Yr), with,sum(Obs)))
  if (length(names(ts))==8)
    names(ts)=c("year","area","era","season","biomass","ssb","rec","catch")
  else
    names(ts)=c("year","era","season","biomass","ssb","rec","catch")
  
  ts}

getRf<-function(x){
  names(x$derived_quants)=tolower(names(x$derived_quants))
  rf=subset(x$derived_quants,
            label%in%c("SSB_Unfished","SSB_MSY","Fstd_MSY","TotYield_MSY"))[,1:5]
  rf[,1]=rfs[rf[,1]]
  dimnames(rf)[[1]]=rf[,1]
  
  cbind("quant"=c("hat","var"),as.data.frame(t(rf[,2:3])))}

smrySS<-function(x,covar=TRUE,forecast=TRUE,ncols=250){
  
  ss=llply(x, function(x) SS_output(x, verbose=FALSE,printstats=FALSE,covar=covar,forecast=forecast,ncols=ncols))
  ts=ldply(ss, getTs) 
  pf=ldply(ss, getPellaT)
  rf=ldply(ss, getRf)

  if (!is.null(attributes(x)$split_labels)){
    ts=cbind(attributes(x)$split_labels[ts$.id,],ts)
    pf=cbind(attributes(x)$split_labels[pf$.id,],pf)
    rf=cbind(attributes(x)$split_labels[rf$.id,],rf)
    }
    
  return(list(timeseries=ts,refpts=rf,pfunc=transform(pf,m=1+p)))}

getRefpts<-function(object,value=TRUE){
  names(object$derived_quants)=tolower(names(object$derived_quants))
  rf=subset(object$derived_quants,label%in%rfs[c(1,12,15)])[,2:3]
  names(rf)=c("value","var")
  dimnames(rf)[[1]]=c("k","bmsy","msy")
  rf=t(as.matrix(rf))
  rf["var",]=rf["var",]^2
  rf=as.data.frame(rf)
  
  rf[ifelse(value,1,2),]}

getPellaT<-function(object){
  rfs=c("SSB_Unfished",    "TotBio_Unfished",  "SmryBio_Unfished", "Recr_Unfished",    
        "SSB_Btgt",        "SPR_Btgt",         "Fstd_Btgt",        "TotYield_Btgt",   
        "SSB_SPRtgt",      "Fstd_SPRtgt",      "TotYield_SPRtgt",  "SSB_MSY",         
        "SPR_MSY",         "Fstd_MSY",         "TotYield_MSY",     "RetYield_MSY",
        "SSB_unfished","Dead_Catch_MSY") 
  
  names(object$derived_quants)=tolower(names(object$derived_quants))
  rf=subset(object$derived_quants,label%in%rfs[c(1,17,12,15,18)])[,2:3]
  names(rf)=c("value","var")
  dimnames(rf)[[1]]=c("k","bmsy","msy")
  rf=t(as.matrix(rf))
  rf["var",]=rf["var",]^2
  rf=as.data.frame(rf)
  
  pt=transmute(rf["value",],shape=bmsy/k,
               k    =k,
               p    =optimise(function(x,y) (y-(1/(1+x))^(1/x))^2,
                              c(-0.9999,10),y=shape)$minimum,
               r    =(1+p)*(msy/bmsy))
  
  pt[,c("k","r","p","shape")]}

productionFn<-function(b,r,k,p){
  t1=b*r/p
  t3=(b/k)^p
  t1*(1-t3)}

hat<-function(ssb,r,k,p) tMSE:::productionFn(ssb,r,k,p)
pe <-function(year,ssb,catch,hat){
      data.frame(year=year[-length(year)],
                 pe=log((ssb[-1]+catch[-length(ssb)]-hat[-length(ssb)])/ssb[-length(ssb)]))}

