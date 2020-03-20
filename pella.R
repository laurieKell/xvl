fnlim<-function(p=1.01,r=0.5,B=0,K=1){
  fd<-deriv(~r/p*B*(1-(B/K)^p),"B")
  
  attributes(eval(fd))$gradient}

fnmsy<-function(p=1.01,r=0.5)
  r*(1/(1+p))

fn<-function(par=c(0.5,1.01),fmsy=0.25,flim=0.5) {
  r=par[1]
  p=par[2]
  (flim-fnlim(p,r))^2+(fmsy-fnmsy(p,r))^2}

optim(c(r=0.5,p=1.01),fn,flim=0.5,fmsy=0.25)$par
  