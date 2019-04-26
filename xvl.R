library(r4ss)
library(stringr)
library(plyr)
library(doParallel)
library(foreach)
library(dplyr)

library(devtools)
install_github("lauriekell/xvl",force=TRUE)

library(xvl)

cl = makeCluster(3)
registerDoParallel(cl)
      
x        ="/home/laurence/Desktop/mak/run1/data.ss"
ss1      =SS_output(dirname(x))
pf1      =smrySS(dirname(x))
h1       =runHcstYr(x,1)
hindcast1=runHcst(x,1)
j1       =runJK(x)

## baltic cod
x="/home/laurence/Desktop/Dropbox/baltic-ss/EBC2019/EBcod_wgbfas19_dat.ss"
ss1      =SS_output(dirname(x))
pf1      =smrySS(dirname(x))
h1       =runHcstYr(x,1,TRUE)
hindcast1=runHcst(x,1,TRUE)

names(h1[[2]])[5:18]=names(h1[[2]])[4:17]
names(h1[[2]])[4]="what"

lbl=merge(ddply(transform(ss1$cpue,name=Fleet_name),.(name),with,data.frame(y=max(Obs))),
          ddply(subset(h1[[1]],year==tail+1),.(name),
                with,as.character(paste("MASE =",signif(maseSS(obs,hat),4)))))
lbl=cbind(lbl,x=min(ss1$cpue$Yr)+5)

ggplot(subset(h1[[2]],year<=tail+2))+
  geom_point(aes(Yr,Obs),data=transform(ss1$cpue,name=Fleet_name),col="red")+
  facet_grid(name~.,scale="free")+theme_bw()+
  ylab("CPUE")+xlab("Year")+
  geom_text(data=lbl,aes(x=x,y=y*.9,label=V1))+
  geom_line(aes(year,hat,group=tail),lwd=.2,col="grey25")+
  geom_point(aes(year,hat,group=tail),data=subset(h1[[2]],year==tail+2))
  