require(r4ss)
require(stringr)
require(plyr)

require(doParallel)
require(foreach)

## Sets up the files for the jackknife
## Sets up the files for the jackknife
setJK<-function(x){
  
  ## process file
  dirNow=getwd()
  
  #get rid of comments in data file
  dfl =str_trim(scan(x,sep="\n",what=as.character()))
  dfl=dfl[substr(dfl,1,1)!="#"]
  dfl=dfl[nchar(dfl)!=0]
  
  # function to get count number of rows
  rw<-function(x) as.numeric(strsplit(x,"\\s+")[[1]][1])
  
  # number of fleets and surveys  
  nFlt=rw(dfl[6])
  nSry=rw(dfl[7])
  
  # rows with data
  rCtc=seq(rw(dfl[17]))+17
  rU  =max(rCtc)+1+nFlt+nSry+seq(rw(dfl[max(rCtc)+1]))
  nDsc=rw(dfl[max(rU)+1])
  rDsc=max(rU)+nDsc+2+seq(rw(dfl[max(rU)+nDsc+2]))
  nLnc=rw(dfl[max(rDsc)+9])
  rLnc=(max(rDsc)+12):(length(dfl)-13)
  
  ## key for obs, i.e. indices and length comps
  u  =mdply(rU,  function(i) as.numeric(strsplit(dfl[i],"\\s+")[[1]])[1:5])[,-1]
  
  ## No length comps for now
  #lnc=mdply(rLnc,function(i) as.numeric(strsplit(dfl[i],"\\s")[[1]])[1:75])
  
  names(u)=c("year","season","fleet","obs","cv")
  
  u=cbind(u,row=rU)
  
  list(dfl=dfl,u=u)}

setJKNew<-function(x){
  
  ## process file
  dirNow=getwd()
  
  #get rid of comments in data file
  dfl =str_trim(scan(x,sep="\n",what=as.character()))
  dfl=dfl[substr(dfl,1,1)!="#"]
  dfl=dfl[nchar(dfl)!=0]
  
  # function to get count number of rows
  rw<-function(x) as.numeric(strsplit(x,"\\s+")[[1]][1])
  
  nSeas=rw(dfl[3])
  minYr=rw(dfl[1])
  maxYr=rw(dfl[2])
  
  # number of fleets and surveys  
  nFlt=rw(dfl[10])
 
  rU=(seq(length(dfl))[maply(dfl,rw)==-9999][1]+nFlt+1):
     (seq(length(dfl))[maply(dfl,rw)==-9999][2]-1)

  ## key for obs, i.e. indices and length comps
  u  =mdply(rU,  function(i) as.numeric(strsplit(dfl[i],"\\s+")[[1]])[1:5])[,-1]
  
  names(u)=c("year","season","fleet","obs","cv")
  
  u=cbind(u,row=rU)
  
  list(dfl=dfl,u=u)}

mkTmp<-function(){
  dr=tempfile()
  dir.create(dr)
  
  dr}

jkU<-function(i,u,tfl,dat,newVer=FALSE){
  ## copy files from target 
  dirNow=getwd()
  dirTmp=mkTmp()
  setwd(dirTmp)

  file.copy(file.path(dirname(dat),"."),dirTmp,r=T)
  
  #leave out obs
  u[,"fleet"]=-u[,"fleet"]
  for (j in i)
    tfl[j]=paste(unlist(subset(u[,-5],j==row)),sep=" ",collapse=" ")
  cat(tfl,sep="\n",file=file.path(dirTmp,substr(dat,nchar(dirname(dat))+2,nchar(dat))))
  
  #run
  #if (newVer)
  #   system2("wine",args="ss.exe -nohess",stdout=NULL)
  #else  
  #   system2("./ss3_3.24z",args="-nohess",stdout=NULL)
  
  # Linux
  if (R.version$os=='linux-gnu') {
    exe=paste(system.file('bin', 'linux', package="xvl", mustWork=TRUE),
                    ifelse(newVer,"ss_opt","ss3_3.24z"), sep='/')
    file.copy(exe, dirTmp)
    system2(ifelse(newVer,"./ss_opt","./ss3_3.24z"),args="-nohess",stdout=NULL)
  # Windows
  } else if (.Platform$OS.type=='windows') {
    exe = paste(system.file('bin', 'windows', package="xvl", mustWork=TRUE), 
                ifelse(newVer,"SS.exe","ss3.exe"), sep='/')
    
    file.copy(exe, dirTmp)
    
    system2(ifelse(newVer,"ss.exe","ss3.exe"),args="-nohess",stdout=NULL)
  }else 
    stop()

  #get results
  ssf=SS_output(getwd(), 
                forecast  =FALSE, 
                covar     =FALSE,
                checkcor  =FALSE,
                verbose   =FALSE, 
                printstats=FALSE, 
                hidewarn  =TRUE, 
                NoCompOK  =TRUE)
  
  rfs=c("SSB_Unfished",    "TotBio_Unfished", "SmryBio_Unfished", "Recr_Unfished",  "SSB_Btgt",        
        "SPR_Btgt",        "Fstd_Btgt",       "TotYield_Btgt",    "SSB_SPRtgt",     "Fstd_SPRtgt",     
        "TotYield_SPRtgt", "SSB_MSY",         "SPR_MSY",          "Fstd_MSY",       "TotYield_MSY",    
        "RetYield_MSY") 
  
  names(ssf$derived_quants)=tolower(names(ssf$derived_quants))
  rf=subset(ssf$derived_quants,label%in%rfs)
  ts=ssf$timeseries[,1:8]
  
  #clean up  
  setwd(dirNow)
  file.remove(dirTmp,r=T)
  
  list(u=ssf$cpue,rf=rf,ts=ts)}

runHcst<-function(x,n=10,newVer=FALSE){

  ## process files
  if (newVer)
     fls=setJKNew(x)
  else   
     fls=setJK(x)
 
  key=ddply(fls$u,.(fleet), transform, maxyr=max(year))
  key=subset(key,year>=(maxyr-n)&year<maxyr)[,-7]
  
  dir=dirname(x)
  dir.create(file.path(dir,"hcast"))
  
  hRsd=foreach(i=seq(dim(key)[1]),
       .multicombine=TRUE,
       .combine     =rbind.fill,
       .packages    =c("xvl","r4ss")) %dopar%{

       iRw=subset(fls$u,fleet==key[i,"fleet"]&year>key[i,"year"])[,"row"]
       res=xvl:::jkU(iRw,fls$u,fls$dfl,x,newVer)
       rtn=cbind(tail =key[i,"year"],
                 naive=fls$u[as.numeric(dimnames(subset(fls$u,fleet==key[i,"fleet"]&year==key[i,"year"]))[[1]]),"obs"],
                 subset(res$u,Fleet==key[i,"fleet"]&Yr>=key[i,"year"]))[,1:13]

       write.table(res[[1]],file=file.path(dir, "hcast",paste("rsd",i,".csv",sep="")))
       write.table(res[[2]],file=file.path(dir, "hcast",paste("ref",i,".csv",sep="")))
       write.table(res[[3]],file=file.path(dir, "hcast",paste("ts" ,i,".csv",sep="")))

       names(rtn)[3:13]=c("fleet","name","year","season","year.","vuln","obs","hat","q","eff","se")

       rtn}
  
  rsdl=mdply(data.frame(i=seq(dim(key)[1])),function(i)
    read.csv(file.path(dir,"hcast",paste("rsd",i,".csv",sep="")),header=T,sep=" "))
  if (newVer)
    names(rsdl)=c("key","fleet","name","area","year","season","year.","vulnerable","obs","hat","q","q.","se",
                  "dev","like","like.","sp","use")
  else
    names(rsdl)=c("key","fleet","name","year","season","year.","vulnerable","obs","hat","q","q.","se",
                  "dev","like","like.","sp","use")
  
  ts  =mdply(data.frame(i=seq(dim(key)[1])),function(i) 
    read.csv(file.path(dir,"hcast",paste("ts",i,".csv",sep="")),header=T,sep=" "))
  
  names(ts)=c("key","area","year","era","season","biomass","biomass.","ssb","rec")
  
  rf  =mdply(data.frame(i=seq(dim(key)[1])),function(i)   
    read.csv(file.path(dir,"hcast",paste("ref",i,".csv",sep="")),header=T,sep=" "))
  rf=rf[,1:3]
  names(rf)=c("key","variable","value")
  
  return(list(hindcast  =hRsd,
              residuals =rsdl,
              timeseries=ts,
              refpts    =rf,
              key       =data.frame("key"=seq(dim(key)[1]),key)))}

runHcstYr<-function(x,n=5,newVer=FALSE){
  
  ## process files
  if (newVer)
    fls=setJKNew(x)
  else
    fls=setJK(x)
  
  key=subset(fls$u,year>=max(year)-n+1)
  yrs=rev(sort(unique(key$year)))
  
  dir=dirname(x)
  dir.create(file.path(dir,"hyrs"))
  
  hRsd=NULL
  hRsd=foreach(i=yrs[seq(n)],
               .multicombine=TRUE,
               .combine     =rbind.fill,
               .packages    =c("xvl","r4ss")) %dopar%{
  #for(i in yrs[seq(n)]){
  
     iRw=subset(fls$u,year>=i)[,"row"]
     res=xvl:::jkU(iRw,fls$u,fls$dfl,x,newVer)
    
     ##bug
     #ggplot(subset(h1[[1]],year==tail+1))+
     #   geom_point(aes(year,obs))+
     #   geom_point(aes(year,hat),col="red")+
     #   geom_point(aes(year,naive),col="green")+
     #   geom_line(aes(year,obs))+
     #   geom_line(aes(year,hat),col="red")+
     #   geom_line(aes(year,naive),col="green")+
     #   facet_grid(name~.,scale="free")
     naive=subset(fls$u,year==i)[,c("fleet","obs")]
     names(naive)[2]="naive"
     
     
     if (newVer){
         rtn=cbind(tail=i,subset(res$u,Yr>=i))[,1:13]
         names(rtn)[2:13]=c("fleet","name","area","year","season","year.","vuln","obs","hat","q","eff","se")
     } else {
         rtn=cbind(tail=i, subset(res$u,Yr>=i))[,1:12]
         names(rtn)[2:12]=c("fleet","name","year","season","year.","vuln","obs","hat","q","eff","se")
     }
     
     rtn=merge(rtn,naive,by="fleet")
     
     write.table(res[[1]],file=file.path(dir, "hyrs",paste("rsd",i,".csv",sep="")))
     write.table(res[[2]],file=file.path(dir, "hyrs",paste("ref",i,".csv",sep="")))
     write.table(res[[3]],file=file.path(dir, "hyrs",paste("ts" ,i,".csv",sep="")))
                 
     
     rtn
     #hRsd=rbind.fill(hRsd,rtn)
     }
  
  rsdl=mdply(data.frame(i=yrs[seq(n)]),function(i)
        read.csv(file.path(dir,"hyrs",paste("rsd",i,".csv",sep="")),header=T,sep=" "))
  if (dim(rsdl)[2]==20)
    names(rsdl)=c("tail","fleet","name","area","year","season","year.","vulnerable","obs","hat","q","q.","se",
                  "dev","like","like.","sp","use")
  else
    names(rsdl)=c("tail","fleet","name","year","season","year.","vulnerable","obs","hat","q","q.","se",
                 "dev","like","like.","sp","use")
  
  ts  =mdply(data.frame(i=yrs[seq(n)]),function(i) 
        read.csv(file.path(dir,"hyrs",paste("ts",i,".csv",sep="")),header=T,sep=" "))
  names(ts)=c("tail","area","year","era","season","biomass","biomass.","ssb","rec")
  
  rf  =mdply(data.frame(i=yrs[seq(n)]),function(i)   
        read.csv(file.path(dir,"hyrs",paste("ref",i,".csv",sep="")),header=T,sep=" "))
  rf=rf[,1:3]
  names(rf)=c("tail","variable","value")
  
  return(list(hindcast  =hRsd,
              residuals =rsdl,
              timeseries=ts,
              refpts    =rf,
              key       =data.frame("key"=seq(dim(key)[1]),key)))}

runJK<-function(x){ 
  
  ## process files
  dirNow=getwd()
  fls   =setJK(x)
  dir   =dirname(x)
  dat   =substr(x,nchar(dir)+2,nchar(x))
  
  dirX=file.path(dir, "xval")
  dir.create(dirX)
  
  pRsd=NULL
  pRsd=foreach(i=seq(length(fls$u$row)),
     .multicombine=TRUE,
     .combine     =rbind,
     .packages    =c("r4ss","xvl")) %dopar%{
                
     ## copy files from target 
     dirTmp=xvl:::mkTmp()
     setwd(dirTmp)
                
     #file.copy(file.path(dirname(dir),"/"),dirTmp,r=T)
     file.copy(file.path(dirname(dir),"."),,dirTmp,r=T)
     
     iRw=fls$u[i,"row"]
     res=xvl:::jkU(iRw,fls$u,fls$dfl,file.path(dirTmp,dat))
                
     write.table(res[[1]],file=file.path(dir, "xval",paste("prd",i,".csv",sep="")))
     write.table(res[[2]],file=file.path(dir, "xval",paste("ref",i,".csv",sep="")))
     write.table(res[[3]],file=file.path(dir, "xval",paste("ts" ,i,".csv",sep="")))
                
     #clean up  
     setwd(dirNow)
     file.remove(dirTmp,r=TRUE)
                
     res[[1]][i,]}
  
  pRsd=pRsd[,1:11]
  names(pRsd)=c("fleet","name","year","season","year.","vulnerable","obs","hat","q","q.","se")
  
  #rsdl=mdply(data.frame(i=seq(length(fls$u$row))),function(i)
  #  read.csv(file.path(dirX,paste("prd",i,".csv",sep="")),header=T,sep=" ")[i,])

  ts  =mdply(data.frame(i=seq(length(fls$u$row))),function(i) 
    read.csv(file.path(dirX,paste("ts",i,".csv",sep="")),header=T,sep=" "))
  names(ts)=c("i","area","year","era","season","biomass","biomass.","ssb")
  
  rf  =mdply(data.frame(i=seq(length(fls$u$row))),function(i)   
    read.csv(file.path(dirX,paste("ref",i,".csv",sep="")),header=T,sep=" "))
  rf=rf[,1:3]
  names(rf)=c("i","variable","value")
  
  return(list(prediction=pRsd,timeseries=ts,refpts=rf))}

runJKBlock<-function(x,n=5){ 
  
  ## process files
  dirNow=getwd()
  fls   =setJK(x)
  dir   =dirname(x)
  dat   =substr(x,nchar(dir)+2,nchar(x))
  
  dirX=file.path(dir, "xvalBlock")
  dir.create(dirX)
  
  key=ddply(subset(fls$u,year>0),.(fleet), with, 
            {x=c(seq(max(year),min(year),-5),min(year));
            x=x[!duplicated(x)]
            data.frame(max=x[-length(x)],min=x[-1])})
  
  pRsd=NULL
  pRsd=foreach(i=seq(dim(key)[1]),
               .multicombine=TRUE,
               .combine     =rbind,
               .packages    =c("r4ss","xvl")) %dopar%{
                 
                ## copy files from target 
                dirTmp=mkTmp()
                setwd(dirTmp)
                 
                #file.copy(file.path(dirname(dat),"/"),dirTmp,r=T)
                file.copy(file.path(dirname(dat),"."),,dirTmp,r=T)
                
                iRw=subset(fls$u,fleet==key[i,"fleet"]&year>key[i,"min"]&year<=key[i,"max"])[,"row"]
                res=xvl:::jkU(iRw,fls$u,fls$dfl,file.path(dirTmp,dat))
                 
                write.table(res[[1]],file=file.path(dir, "xvalBlock",paste("prd",i,".csv",sep="")))
                write.table(res[[2]],file=file.path(dir, "xvalBlock",paste("ref",i,".csv",sep="")))
                write.table(res[[3]],file=file.path(dir, "xvalBlock",paste("ts" ,i,".csv",sep="")))
                 
                #clean up  
                setwd(dirNow)
                file.remove(dirTmp,r=TRUE)
               
                names(res[[1]])[1:11]=c("fleet","name","year","season","year.","vuln","obs","hat","q","eff","se")

                subset(res[[1]][,1:11],fleet==key[i,"fleet"]&year>key[i,"min"]&year<=key[i,"max"])}
  
  #rsdl=mdply(data.frame(i=seq(length(fls$u$row))),function(i)
  #  read.csv(file.path(dirX,paste("prd",i,".csv",sep="")),header=T,sep=" ")[i,])
  
  ts  =mdply(data.frame(i=seq(length(fls$u$row))),function(i) 
    read.csv(file.path(dirX,paste("ts",i,".csv",sep="")),header=T,sep=" "))
  names(ts)=c("i","area","year","era","season","biomass","biomass.","ssb")
  
  rf  =mdply(data.frame(i=seq(length(fls$u$row))),function(i)   
    read.csv(file.path(dirX,paste("ref",i,".csv",sep="")),header=T,sep=" "))
  rf=rf[,1:3]
  names(rf)=c("i","variable","value")
  
  return(list(prediction=pRsd,timeseries=ts,refpts=rf))}
