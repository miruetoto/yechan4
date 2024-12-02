ebt2<-function(f,tau)
{
  fsave<-f
  f<-mirror(f)
  len<-length(f)
  t=1:len
  sampled.index<-list()
  missing.index<-list()
  band<-rep(0,len*tau); dim(band)<-c(len,tau)
  
  for (eta in 1:tau){
   sampled.index[[eta]]<-seq(from=eta,to=len,by=tau)
   if(sampled.index[[eta]][length(sampled.index[[eta]])]!=len) sampled.index[[eta]]<-c(sampled.index[[eta]],len)
   missing.index[[eta]]<-(1:len)[-sampled.index[[eta]]]
   band[,eta]<-lin.impute(t,f,missing.index[[eta]])
  }
  L<-apply(band,1,min)
  U<-apply(band,1,max)
  M<-apply(band,1,mean)
  V=U-L
  ltemp<-ceiling(length(fsave)/2)
  index<-ltemp:(len-ltemp-1)
  list(t=1:length(index),f=fsave,V=V[index],L=L[index],U=U[index],M=M[index],tau=tau,band=band[index,],sampled.index=sampled.index)
}

decomp<-function(f,k=4,maxtau=floor(length(f)/2))
{
  len<-length(f);
  power<-sum(f^2)
  modes<-rep(0,maxtau*len); dim(modes)<-c(maxtau,len)
  ampofmodes<-rep(0,maxtau*len); dim(ampofmodes)<-c(maxtau,len)
  modes[1,]<-f
  if(ext)
  modes[2,]<-ext(f=f,tau=2,eps=0.000000000000000000000000005,iter=k)
  ampofmodes[2,]<-ebt(f=modes[2,],tau=2)$V
  res<-modes[1,]-modes[2,]
  for(tau in 3:maxtau){
      modes[tau,]<-ext(f=res,tau=tau,eps=0.000000000000000000000000005,iter=k)
      res<-res-modes[tau,]
      ampofmodes[tau,]<-ebt(f=modes[tau,],tau=tau)$V
  }
  modes[1,]<-res
  ampofmodes[1,]<-res
  list(modes=t(modes),ampofmodes=t(ampofmodes))
}

ext<-function(f,tau,iter=100000,eps=0.05,M="mean",V="volume")
{
  m.extracted.old<-f
  m.extracted<-ebt(f,tau=tau,M=M,V=V)$M
  res<-f-m.extracted
  for(j in 2:iter)
  {
    m.extracted.old<-m.extracted
    m.extracted<-ebt(res,tau=tau,M=M,V=V)$M
    res<-res-m.extracted
  }
  highest.freq<-res
}

ebt<-function(f,t=0,noise=0,tau,M="mean",V="var",interpolation="linear",bound="fix",sampling="equalspace")
{
  fsave<-f
  tsave<-t+(length(t)==1)*1:length(f)
  f<-c(rep(f[1],tau*2),f,rep(f[length(f)],tau*2))
  t=1:length(f)
  if(sum(noise)==0) noise=0*f
  len<-length(f)
  sampled.index<-list()
  missing.index<-list()
  if (interpolation=="const")
  {
    lband<-rep(0,len*tau); dim(lband)<-c(len,tau)
    rband<-rep(0,len*tau); dim(rband)<-c(len,tau)
    band<-rep(0,len*tau*2); dim(band)<-c(len,tau*2)
  }else{
    band<-rep(0,len*tau); dim(band)<-c(len,tau)
  }

  if (interpolation=="cubic")
  {
    for (eta in 1:tau){
      if (sampling=="random") sampled.index[[eta]]<-sort(sample(1:len,1/tau*1000,replace=F))
      else if(sampling=="equalspace") sampled.index[[eta]]<-seq(from=eta,to=len,by=tau)
      if (sampled.index[[eta]][1]!=1) sampled.index[[eta]]=c(1,sampled.index[[eta]])
      if (sampled.index[[eta]][length(sampled.index[[eta]])]!=len) sampled.index[[eta]]=c(sampled.index[[eta]],len)
      missing.index[[eta]]<-(1:len)[-sampled.index[[eta]]]
      result<-smooth.spline(t[sampled.index[[eta]]],f[sampled.index[[eta]]]+noise[sampled.index[[eta]]],spar=0,all.knots=TRUE)
      band[,eta][sampled.index[[eta]]]<-result$y
      band[,eta][missing.index[[eta]]]<-predict(result,t[missing.index[[eta]]])$y
    }
  }else if (interpolation=="linear"){
    for (eta in 1:tau){
      if (sampling=="random") sampled.index[[eta]]<-sort(sample(1:len,1/tau*1000,replace=F))
      else if(sampling=="equalspace") sampled.index[[eta]]<-seq(from=eta,to=len,by=tau)
      if (sampled.index[[eta]][1]!=1) sampled.index[[eta]]<-c(1,sampled.index[[eta]])
      if (sampled.index[[eta]][length(sampled.index[[eta]])]!=len) sampled.index[[eta]]<-c(sampled.index[[eta]],len)
      missing.index[[eta]]<-(1:len)[-sampled.index[[eta]]]
      band[,eta]<-lin.impute(t,f+noise,missing.index[[eta]])
    }
  }else if (interpolation=="const"){
    for (eta in 1:tau){
      if (sampling=="random") sampled.index[[eta]]<-sort(sample(1:len,1/tau*1000,replace=F))
      else if(sampling=="equalspace") sampled.index[[eta]]<-seq(from=eta,to=len,by=tau)
      if (sampled.index[[eta]][1]!=1) sampled.index[[eta]]=c(1,sampled.index[[eta]])
      if (sampled.index[[eta]][length(sampled.index[[eta]])]!=len) sampled.index[[eta]]=c(sampled.index[[eta]],len)
      missing.index[[eta]]<-(1:len)[-sampled.index[[eta]]]
      lband[,eta]<-left.const(f+noise,missing.index[[eta]])
      rband[,eta]<-right.const(f+noise,missing.index[[eta]])
    }
    band<-cbind(lband,rband)
  }
  L<-apply(band,1,min)
  U<-apply(band,1,max)

  M1<-apply(band,1,mean)
  M2<-apply(band,1,median)
  M3<-(L+U)/2
  V1=U-L
  V2=apply(band,1,var) ; V2=V2*(tau-1)/tau;

  if (M=="mean") M<-M1
  else if (M=="median") M<-M2
  else if (M=="volume") M<-M3

  if (V=="volume"){
    V<-V1
    L<-L
    U<-U
  }
  else if (V=="var"){
    V<-V2
    L<-M-V
    U<-M+V
  }

  index<-(2*tau+1):(length(f)-2*tau)
  list(t=tsave,f=fsave,V=V[index],L=L[index],U=U[index],M=M[index],tau=tau,band=band[index,],sampled.index=sampled.index)
}

vis<-function(ebt,band=0,bandcol=0,bandlwd=1,bandpch=0,V=T,M=T,mse=F,main="",xlim=F,ylim=F,obs="point",size=0.5)
{
  if(length(ylim)==1) ylim<-c(min(ebt$L),max(ebt$U)) else ylim<-ylim 
  if(length(xlim)==1) xlim<-range(ebt$t) else xlim<-xlim
  sampled.index<-ebt$sampled.index
  
  if (mse==T) plot(ebt$t,ebt$f,ylim=ylim,xlim=xlim,type='l',main=paste("tau=",ebt$tau, "mse=",ebt$mse,main),xlab="",ylab="",lwd=4,col="gray60")
  else {
    if (obs=="point")  plot(ebt$t,ebt$f,ylim=ylim,xlim=xlim,xlab="",ylab="",main=main,cex=0.5,col="gray60")
    else  plot(ebt$t,ebt$f,ylim=ylim,xlim=xlim,xlab="",ylab="",main=main,col=1,type='l')
  }
  
  if (V==T){
    lines(ebt$t,ebt$U,col=2,lwd=4,lty=1)
    lines(ebt$t,ebt$L,col=4,lwd=4,lty=1)
  }
  if (M==T) lines(ebt$t,ebt$M,col=2,lwd=4)
  if (band[1]!=0)
  {
    bandindex<-band
    if (bandcol[1]==0) bandcol=(1:length(bandindex))+1
    if (bandpch[1]==0) bandpch=(1:length(bandindex))
    for (i in bandindex)
    {
      #lines(ebt$t,ebt$band[,i],col=bandcol[which(bandindex==i)],lwd=2)
      lines(ebt$t,ebt$band[,i],col="gray60",lwd=bandlwd)
      if (obs=="point") points(ebt$t[sampled.index[[i]]],ebt$band[,i][sampled.index[[i]]],col=bandcol[which(bandindex==i)],cex=size)    
    }       
  }
  
  if (mse==T) lines(ebt$t,ebt$f,ylim=ylim,xlim=xlim,xlab="",ylab="",lwd=4,col="gray60")
  else {
    if (obs=="point")  lines(ebt$t,ebt$f,ylim=ylim,xlim=xlim,xlab="",ylab="",cex=0.5,col="gray60")
    else  lines(ebt$t,ebt$f,ylim=ylim,xlim=xlim,xlab="",ylab="",col=1,type='l',lwd=1)
  }
  
}

swave<-function(freq=1,len)
{
  unit<-c(rep(-1,freq),rep(1,freq))
  sqwave<-rep(unit,ceiling(len/(freq*2)))
  sqwave<-sqwave[1:len]
  sqwave
}

linter<-function(p1,p2,x)
{
  x1<-p1[1]; y1<-p1[2]
  x2<-p2[1]; y2<-p2[2]
  m<-abs(x1-x)
  n<-abs(x2-x)
  y<-(m*y2+n*y1)/(m+n)
  y
}

lextra<-function(p1,p2,x)
{
  x1<-p1[1]; y1<-p1[2]
  x2<-p2[1]; y2<-p2[2]
  m<-abs(x1-x)
  n<-abs(x2-x)
  y<-(m*y2-n*y1)/(m-n)
  y
}

lin.impute<-function(t,y,mindex)
{
  len<- length(y)
  index<-1:len 
  mindex <- sort(mindex)
  oindex <- (1:len)[-mindex]
  yin<-rep(NA,len)
  yin[oindex]<-y[oindex]
  
  if(mindex[1]==1)
  {
    index.imputed<-(1:(oindex[1]-1))
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-lextra(p1=c(t[oindex[1]],y[oindex[1]]),p2=c(t[oindex[2]],y[oindex[2]]),x=t[index.imputed[j]])
    }
  }
  
  for (i in 1:(length(oindex)-1))
  {
    if(oindex[i+1]-oindex[i]>1) 
    {
      index.imputed1<-(index[oindex[i]]:index[oindex[i+1]]) ## index which need to be interpolated + bound
      index.imputed2<-index.imputed1[-c(1,length(index.imputed1))] ## index which need to be interpolated
      len2<-length(index.imputed2)
      for (j in 1:len2)
      {
        yin[index.imputed2[j]]<-linter(p1=c(t[index.imputed1[1]],y[index.imputed1[1]]),p2=c(t[index.imputed1[length(index.imputed1)]],y[index.imputed1[length(index.imputed1)]]),x=t[index.imputed2[j]])
      }
    }
  } 
  if((len-oindex[length(oindex)])>0)
  {
    index.imputed<-(oindex[length(oindex)]+1):len
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-lextra(p1=c(t[oindex[length(oindex)-1]],y[oindex[length(oindex)]-1]),p2=c(t[oindex[length(oindex)]],y[oindex[length(oindex)]]),x=t[index.imputed[j]])
    }
  }
  yin
}

left.const<-function(y,mindex)
{
  len<- length(y)
  index<-1:len 
  mindex <- sort(mindex)
  oindex <- (1:len)[-mindex]
  yin<-rep(NA,len)
  yin[oindex]<-y[oindex]
  
  if(mindex[1]==1)
  {
    index.imputed<-(1:(oindex[1]-1))
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-y[oindex[1]]
    }
  }
  
  for (i in 1:(length(oindex)-1))
  {
    if(oindex[i+1]-oindex[i]>1) 
    {
      index.imputed<-(index[oindex[i]]:index[oindex[i+1]])
      index.imputed<-index.imputed[-c(1,length(index.imputed))] ## index which need interpolated
      len2<-length(index.imputed)
      for (j in 1:len2)
      {
        yin[index.imputed[j]]<-y[oindex[i]]
      }
    }
  }
  
  if((len-oindex[length(oindex)])>0)
  {
    index.imputed<-(oindex[length(oindex)]+1):len
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-y[oindex[length(oindex)]]
    }
  }
  yin
}


right.const<-function(y,mindex)
{
  len<- length(y)
  index<-1:len 
  mindex <- sort(mindex)
  oindex <- (1:len)[-mindex]
  yin<-rep(NA,len)
  yin[oindex]<-y[oindex]
  
  if(mindex[1]==1)
  {
    index.imputed<-(1:(oindex[1]-1))
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-y[oindex[1]]
    }
  }
  
  
  for (i in 1:(length(oindex)-1))
  {
    if(oindex[i+1]-oindex[i]>1) 
    {
      index.imputed<-(index[oindex[i]]:index[oindex[i+1]])
      index.imputed<-index.imputed[-c(1,length(index.imputed))] ## index which need interpolated
      len2<-length(index.imputed)
      for (j in 1:len2)
      {
        yin[index.imputed[j]]<-y[oindex[i+1]]
      }
    }
  }
  
  if((len-oindex[length(oindex)])>0)
  {
    index.imputed<-(oindex[length(oindex)]+1):len
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-y[oindex[length(oindex)]]
    }
  }
  yin
}

extract.hfreq<-function(t,f,tau,M="mean",V="sd",tol=0.05,iter=100,interpolation="linear")
{
  m.extracted.old<-f
  m.extracted<-ebt(f,tau=tau,M=M,V=V,interpolation=interpolation)$M
  res<-f-m.extracted
  j<-0
  while((max(abs(m.extracted.old-m.extracted))>tol)&(j<iter))
  {
    m.extracted.old<-m.extracted
    m.extracted<-ebt(res,tau=tau,M=M,V=V,interpolation=interpolation)$M
    res<-(res-m.extracted)
    j<-j+1
  }
  highest.freq<-res
}



extract.hfreq2<-function(t,f,noise=0,tau,M="mean",V="sd",interpolation="linear",tol=0.05,iter=100)
{
  m.extracted.old<-f
  m.extracted<-ebt(f,t,noise=noise,tau=tau,M=M,V=V,interpolation=interpolation)$M
  m.extracted<-ebt(f,t,noise=noise,tau=2,M=M,V=V,interpolation=interpolation)$M ## 2-step!!
  res<-f+noise-m.extracted
  j<-0
  while((max(abs(m.extracted.old-m.extracted))>tol)&(j<iter))
  {
    m.extracted.old<-m.extracted
    m.extracted<-ebt(res,t,noise=noise,tau=tau,M=M,V=V,interpolation=interpolation)$M
    res<-res-m.extracted
    j<-j+1
    plot(res[500:700],ylim=c(-2,2),main=j,type='l')
  }
  highest.freq<-res
  highest.freq
}

mimpute<-function(t,fm,mindex,tau,M="mean",V="sd",interpolation="cubic",bound="extra")
{
  j<-0
  fm_temp<-fm
  ebt.result<-ebt(fm,t,noise=0,tau=tau,M=M,V=V,interpolation=interpolation,bound=bound)
  fm[mindex]<-ebt.result$M[mindex]
  while((abs(sum(fm_temp-fm))>0.00001)&(j<100))
  {
    fm_temp<-fm
    ebt.result<-ebt(fm,t,noise=0,tau=tau,M=M,V=V,interpolation=interpolation,bound=bound)
    fm[mindex]<-ebt.result$M[mindex]
    j<-j+1
  }
  fm
}

zcr<-function(f)
{
  count<-f*0;
  for(i in 2:length(f))
  {
    if((f[i]>0)&(f[i-1]<0)) count[i]<-1
  }
  sum(count)
}

# method 1 
test.var1<-function(t,f,tau,alpha,plot=F,v)
{
  chiV<-ebt(f,t,tau=tau,interpolation="linear",bound="extra",V="chi")$V
  lb1<-qchisq(alpha/2,tau)
  ub1<-qchisq(1-alpha/2,tau)
  if(plot==T)
  {
    plot(chiV)
    abline(h=lb1,col=2)
    abline(h=ub1,col=2)
    abline(v=v[1],col=4)
    abline(v=v[2],col=4)
    abline(v=v[1]-tau,col=4)
    abline(v=v[2]+tau,col=4)
  }
  mis1<-f*0; mis2<-f*0
  mis1[tau:(len-tau)]<-chiV[tau:(len-tau)]<lb1
  mis2[tau:(len-tau)]<-chiV[tau:(len-tau)]>ub1
  mis1rate<-sum(mis1)/length(tau:(len-tau))
  mis2rate<-sum(mis2)/length(tau:(len-tau))
  boundmis<-which(mis1+mis2==1)
  boundmis<-c(boundmis[((v[1]-tau)<boundmis)&(boundmis<v[1])],boundmis[(v[2]<boundmis)&(boundmis<v[2]+tau)])
  out=list(mis1=mis1*1,mis2=mis2*1,mis1rate=mis1rate,mis2rate=mis2rate,boundmis=boundmis)                                                      
}

test.var2<-function(t,f,tau,alpha,plot=F,v)
{
  wv<-f*0
  for(i in (tau+1):(len-tau))
  {
    wv[i]<-sum((f[(i-tau):(i+tau)])^2)
  }
  lb2<-qchisq(alpha/2,tau*2+1)
  ub2<-qchisq(1-alpha/2,tau*2+1)
  if(plot==T)
  {
    plot(wv)
    abline(h=lb2,col=2)
    abline(h=ub2,col=2)
    abline(v=v[1],col=4)
    abline(v=v[2],col=4)
    abline(v=v[1]-tau,col=4)
    abline(v=v[2]+tau,col=4)
  }
  mis1<-f*0; mis2<-f*0
  mis1[tau:(len-tau)]<-wv[tau:(len-tau)]<lb2
  mis2[tau:(len-tau)]<-wv[tau:(len-tau)]>ub2
  mis1rate<-sum(mis1)/length(tau:(len-tau))
  mis2rate<-sum(mis2)/length(tau:(len-tau))
  boundmis<-which(mis1+mis2==1)
  boundmis<-c(boundmis[((v[1]-tau)<boundmis)&(boundmis<v[1])],boundmis[(v[2]<boundmis)&(boundmis<v[2]+tau)])
  out=list(mis1=mis1*1,mis2=mis2*1,mis1rate=mis1rate,mis2rate=mis2rate,boundmis=boundmis)
  out=list(mis1=mis1*1,mis2=mis2*1,misrate=misrate,boundmis=boundmis)                                                      
}

Vmap<-function(f,maxtau,M="mean",V="var")
{
  VM<-rep(0,length(f)*maxtau); dim(VM)<-c(length(f),maxtau);
  
  for(tau in 2:maxtau)
  {
    VM[,tau]<-ebt(f,tau=tau,M=M,V=V)$V
  }  
  VM
}

Mmap<-function(f,maxtau,M="mean",V="var")
{
  MM<-rep(0,length(f)*maxtau); dim(MM)<-c(length(f),maxtau);
  
  for(tau in 2:maxtau)
  {
    MM[,tau]<-ebt(f,tau=tau,M=M,V=V)$M
  }  
  MM
}

Qmap<-function(f,maxtau)
{
  QM<-rep(0,(length(f)-1)*maxtau); dim(QM)<-c(length(f)-1,maxtau);
  len<-length(f)
  ebtresult<-list()
  for(tau in 2:maxtau)
  {
    ebtresult[[tau-1]]<-ebt(f,tau=tau)
  }
  
  quan<-list()
  for (i in 1:(maxtau-1))
  {
    quan[[i]]<-ebtresult[[i]]$band*0
    for (j in 1:len)
    {
      quan[[i]][j,]<-abs(order(ebtresult[[i]]$band[j,]))/length(ebtresult[[i]]$band[j,])      
    }
  }
  
  for(tau in 2:maxtau)
  {
    QM[,tau]<-apply(abs(diff(quan[[tau-1]])),1,sum)/tau*2
  }
  QM
}

## denoising - piecsewise periodic signal 

mirror<-function(signal)
{
  ###mirroing###
  # Period and sampling frequency of input signal
  T = length(signal);
  f_mirror<-c()
  f_mirror[1:floor(T/2)] = signal[floor(T/2):1];
  f_mirror[floor(T/2+1):floor(3*T/2)] = signal;
  f_mirror[floor(3*T/2+1):floor(2*T)] = signal[T:floor(T/2+1)];
  f = f_mirror;
  f
}  

# remove mirroring
demirror<-function(signal)
{
  T = length(signal)/2
  signal <- signal[ceiling((T/2 + 1):(T/2 + T))]
}
# 
# find.tau1<-function(f,t,maxtau) #len(f) must be even number
# {
#   len = length(f);
#   f_mirror<-mirror(f)
#   index=floor(len/2+1):floor(3*len/2)
#   ss<-rep(0,maxtau*len); dim(ss)<-c(maxtau,len) ##ss[tau,i]
#   opt.tau<-f*0;
#   for(i in 1:len)
#   {
#     for(tau in 3:maxtau)
#     {
#       temp_index1<-(index[i]-tau+1):(index[i])
#       temp_index2<-(index[i]):(index[i]+tau-1)
#       ss[tau,i]<-sum(f_mirror[temp_index1])^2+sum(f_mirror[temp_index2])^2
#     }
#     library(features)
#     ss_temp<-features(1:length(ss[,i]),ss[,i])
#     fget.value<-fget(ss_temp)
#     temp1<-fget.value$crit.pts[(fget.value$crit.pts>5)&(fget.value$curvature>0)]
#         
#     if(is.na(temp1[1])){
#       opt.tau[i]<-maxtau
#     }else if((temp1==0)){
#       opt.tau[i]<-maxtau
#     }  else {
#       opt.tau[i]<-round(min(temp1[temp1>0])) ## first local minimum. 
#     }
#   }  
#   opt.tau
# }
# 
# 
# find.tau2<-function(f,t,maxtau) #len(f) must be even number
# {
#   len = length(f);
#   f_mirror<-mirror(f)
#   index=floor(len/2+1):floor(3*len/2)
#   ss<-rep(0,maxtau*len); dim(ss)<-c(maxtau,len) ##ss[tau,i]
#   opt.tau<-f*0;
#   for(tau in 3:maxtau)
#   {
#     ebt(t,f)
#   }
#   for(i in 1:len)
#   {
#     for(tau in 3:maxtau)
#     {
#       ebt.rslt<-ebt(f,t,tau=tau)
#       temp_index1<-(index[i]-tau+1):(index[i])
#       temp_index2<-(index[i]):(index[i]+tau-1)
#       temp_index<-(index[i]-tau+1):(index[i]+tau-1)
#       ss[tau,i]<-(ebt.rslt$V[i]-(f_mirror[temp_index1]-ebt.rslt$M[i])^2/tau)^2+
#         (ebt.rslt$V[i]-(f_mirror[temp_index2]-ebt.rslt$M[i])^2/tau)^2
#     }
#     library(features)
#     ss_temp<-features(1:length(ss[,i]),ss[,i])
#     fget.value<-fget(ss_temp)
#     temp1<-fget.value$crit.pts[(fget.value$crit.pts>5)&(fget.value$curvature>0)]
#     
#     if(is.na(temp1[1])){
#       opt.tau[i]<-maxtau
#     }else if((temp1==0)){
#       opt.tau[i]<-maxtau
#     }  else {
#       opt.tau[i]<-round(min(temp1[temp1>0])) ## first local minimum. 
#     }
#   }  
#   opt.tau
# }

## denoise
denoise<-function(f,taudf)
{
  f-extract.hfreq(f=f,tau=taudf)
}

denoise3<-function(t=length(f),f,tau.list=50:2,tol=0.05,lambda=0,iter=1,msr="V")
{  
  ## define initial 
  f.denoised<-f*0
  ebt.rslt<-list()
  index2<-rep(0,length(f)*length(tau.list)) ; dim(index2)<-c(length(f),length(tau.list))
  index3<-f*0
  maxtau<-tau.list[1]
  
  ## define thresh function 
  thresh<-function(f,lambda=lambda)
  {
    f[f<lambda]<-0
    f
  }
  
  ## select lambda with clustering
  if(lambda==0)
  {
    
  }
  
  ## run ebt 
  len_tau<-length(tau.list)
  for(i in 1:len_tau)
  {
    ebt.rslt[[i]]<-ebt(f,t,tau=tau.list[i])
  }
  
  ## define D
  D<-list()
  for(i in 1:len_tau)
  {
    if(msr=="D") D[[i]]<-ebt.rslt[[i]]$V+(f-ebt.rslt[[i]]$M)^2 ############# V?Ǵ? D?? ?????ؾ???.
    else D[[i]]<-ebt.rslt[[i]]$V
  }
  
  ## define time index such that D[i]<lambda 
  index1<-list()
  if(lambda==0){
    library(stats)
    for(i in 1:len_tau)
    {
      kout<-kmeans(abs(diff(D[[i]])),centers=2)
      cluster.index<-which(kout$centers==min(kout$centers))
      index1[[i]]<-which(kout$cluster==cluster.index)
    }
  }else{
    for (i in 1:len_tau)
    {
      index1[[i]]<-which(thresh(D[[i]],lambda=lambda)==0)
    } 
  }
  
  ## i=1
  i=1
  index2[index1[[i]],i]<-1
  f.denoised<-(f-extract.hfreq(t,f,tau=tau.list[i],tol=tol,iter=iter))*index2[,i]
  
  ## i=2 
  i=2
  index2[index1[[i]],i]<-1
  index3<-index2[,i]*(1-index2[,1])
  f.denoised<-(f-extract.hfreq(t,f,tau=tau.list[i],tol=tol,iter=iter))*index3+f.denoised   
  
  ## i=3~(len-1)
  for(i in 3:(len_tau-1))
  {
    index2[index1[[i]],i]<-1
    index3<-apply((1-index2)[,1:(i-1)],1,prod)*index2[,i]  ### this line is important!! 
    f.denoised<-(f-extract.hfreq(t,f,tau=tau.list[i],tol=tol,iter=iter))*index3+f.denoised       
  }
  
  ## i=len
  i=len_tau
  index2[index1[[i]],i]<-1
  index3<-apply((1-index2)[,1:(i-1)],1,prod)
  f.denoised<-(f-extract.hfreq(t,f,tau=tau.list[i],tol=tol,iter=iter))*index3+f.denoised    
  f.denoised    
}

denoise4<-function(t,f,tau.list)
{  
  library(EbayesThresh)
  f.denoised<-f*0
  ebt.rslt<-list()
  V<-rep(0,length(f)*length(tau.list)) ; dim(V)<-c(length(f),length(tau.list))
  V.inv=V*0
  V.inter<-f*0
  i=1
  ebt.rslt[[i]]<-ebt(f,t,tau=tau.list[i],interpolation="cubic")
  if (sum(ebayesthresh(ebt.rslt[[i]]$V)>0)==0) V[,i]<-0 else V[which(ebayesthresh(ebt.rslt[[i]]$V)>0),i]<-1
  
  V.inv[,i]<-abs(V[,i]-1)
  f.denoised<-(ebt.rslt[[i]]$M)*V.inv[,i]
  
  i=2
  ebt.rslt[[i]]<-ebt(f,t,tau=tau.list[i])
  if (sum(ebayesthresh(ebt.rslt[[i]]$V)>0)==0) V[,i]<-0 else V[which(ebayesthresh(ebt.rslt[[i]]$V)>0),i]<-1
  V.inv[,i]<-abs(V[,i]-1)
  V.inter<-V[,1]*V.inv[,i]  
  f.denoised<-(ebt.rslt[[i]]$M)*V.inter+f.denoised   
  
  for(i in 3:length(tau.list))
  {
    ebt.rslt[[i]]<-ebt(f,t,tau=tau.list[i])
    if (sum(ebayesthresh(ebt.rslt[[i]]$V)>0)==0) V[,i]<-0 else V[which(ebayesthresh(ebt.rslt[[i]]$V)>0),i]<-1
    V.inv[,i]<-abs(V[,i]-1)
    V.inter<-apply(V[,1:(i-1)],1,prod)*V.inv[,i]  
    f.denoised<-(ebt.rslt[[i]]$M)*V.inter+f.denoised    
  }
  f.denoised  
}

denoise5<-function(t,f,tau.list=50:2,tol=0.05,lambda)
{  
  ## define initial 
  f.denoised<-f*0
  ebt.rslt<-list()
  index2<-rep(0,length(f)*length(tau.list)) ; dim(index2)<-c(length(f),length(tau.list))
  index3<-f*0
  maxtau<-tau.list[1]
  
  ## define thresh function 
  thresh<-function(f,lambda=lambda)
  {
    f[f<lambda]<-0
    f
  }
  
  ## run ebt 
  len_tau<-length(tau.list)
  for(i in 1:len_tau)
  {
    ebt.rslt[[i]]<-ebt(f,t,tau=tau.list[i])
  }
  
  ## define D
  D<-list()
  for(i in 1:len_tau)
  {
    D[[i]]<-ebt.rslt[[i]]$V ############# V?Ǵ? D?? ?????ؾ???.
  }
  
  ## define time index such that D[i]<lambda 
  index1<-list()
  for (i in 1:len_tau)
  {
    index1[[i]]<-which(thresh(D[[i]],lambda=lambda)==0)
  }
  
  ## i=1
  i=1
  index2[index1[[i]],i]<-1
  f.denoised<-(f-extract.hfreq(t,f,tau=tau.list[i],tol=tol))*index2[,i]
  
  ## i=2 
  i=2
  index2[index1[[i]],i]<-1
  index3<-index2[,i]*(1-index2[,1])
  f.denoised<-(f-extract.hfreq(t,f,tau=tau.list[i],tol=tol))*index3+f.denoised   
  
  ## i=3~(len-1)
  for(i in 3:(len_tau-1))
  {
    index2[index1[[i]],i]<-1
    index3<-apply((1-index2)[,1:(i-1)],1,prod)*index2[,i]  ### this line is important!! 
    f.denoised<-(f-extract.hfreq(t,f,tau=tau.list[i],tol=tol))*index3+f.denoised       
  }
  
  ## i=len
  i=len_tau
  index2[index1[[i]],i]<-1
  index3<-apply((1-index2)[,1:(i-1)],1,prod)
  f.denoised<-(f-extract.hfreq(t,f,tau=tau.list[i],tol=tol,iter=1))*index3+f.denoised    
  f.denoised    
}


find.tau<-function(t=1:length(f),f,taulist=2:100)
{
  ebtresult<-list()
  corr<-c()
  corr[1]<-1
  for(i in 1:length(taulist))
  {
    ebtresult[[i+1]]<-ebt(f,t,tau=taulist[i])
    corr[i+1]<-(sum(cor(ebtresult[[i+1]]$band))-i-1)/((i+1)^2-(i+1))
  }
  
  quan<-list()
  for (i in (1:(length(ebtresult)-1)))
  {
    quan[[i+1]]<-ebtresult[[i+1]]$band*0
    for (j in 1:length(ebtresult[[i+1]]$band[,1]))
    {
      quan[[i+1]][j,]<-abs(order(ebtresult[[i+1]]$band[j,]))/length(ebtresult[[i+1]]$band[j,])      
    }
  }
  
  diffquan<-list()
  rank<-c()
  rank[1]<-0
  for(i in 1:length(taulist))
  {
    diffquan[[i+1]]<-diff(quan[[i+1]])
    rank[i+1]<-sum(abs(diffquan[[i+1]]))/(i+1)
  }
  list(corr=corr, rank=rank,diffquan=diffquan)
}

get.h<-function(f,tau,iter=1,plot=F,main="",xlim=c(0,pi),log=F,col=1)
{
  len=length(f)
  ker<-c(1:tau,(tau-1):1)
  len_0<-(len-1)/2-(length(ker)-1)/2 
  #### 0000k k k0000 ####
  h<-c(ker,rep(0,len_0),rep(0,len_0))/tau^2
  if(plot==T){
    if(log==F){
      omega_hat<-seq(from=-pi,to=pi,len=len)
      H_iter<-c(1-(1-abs(fft(h))[((len+1)/2):len])^iter,1-(1-abs(fft(h))[1:((len+1)/2-1)])^iter)
      plot(omega_hat,H_iter,type='l',xlab="",ylab="",main=main,xlim=xlim)
      #abline(v=2*pi/tau,lty=2,lwd=0.5);abline(v=-2*pi/tau,lty=2,lwd=0.5)
    }else{
      omega_hat<-seq(from=-pi,to=pi,len=len)
      H_iter<-c(1-(1-abs(fft(h))[((len+1)/2):len])^iter,1-(1-abs(fft(h))[1:((len+1)/2-1)])^iter)
      plot(omega_hat[omega_hat>0],H_iter[omega_hat>0],type='l',xlab="",ylab="",main=main,xlim=xlim,log="y",yaxt="n",ylim=c(0.0000000001,1),col=col)
      # y-axis
      max<- floor(log10(1))
      min<- floor(log10(0.0000000001))
      by<-floor((max-min)/2)
      temp<-seq(min,max,by=by)
      labels <- sapply(temp, function(i) as.expression(bquote(10^ .(i))))
      axis(2, at=10^temp, labels=labels)
      
      #abline(v=2*pi/tau,lty=2,lwd=0.5);abline(v=-2*pi/tau,lty=2,lwd=0.5)
    }
  }
  h_iter<-Re(fft(1-(1-fft(h))^iter,inverse=T))
  hf<-Re(fft(fft(f)*fft(h_iter),inverse=T))/length(h)
  list(omega_hat=omega_hat,H_iter=H_iter,h=h,hf=hf)
}

FB<-function(x)
{
  k<-1:10000
  1+2*sum((1-4*(k^2)*(x^2))*exp(-2*(k^2)*(x^2)))
}

st.test<-function(Kt)
{
  M<-round(log(length(Kt)),0)
  Kccf<-ccf(Kt,Kt,type="covariance",plot=F)$acf
  argmax_ccf<-which(Kccf==max(Kccf))
  sigma12<-sqrt(Kccf[argmax_ccf]+2*sum(Kccf[argmax_ccf:(argmax_ccf+M)]))
  
  Y<-c()
  Z<-c()
  
  for(i in 1:length(Kt))
  {
    Y[i]<-1/sigma12/sqrt(length(Kt))*(sum(Kt[1:i])-mean(Kt))
  }
  
  for(i in 1:length(Kt))
  {
    Z[i]<-Y[i]-(i/length(Kt))*Y[length(Kt)]
  }
  R<-max(Z)-min(Z)
  pval<-round(FB(R),3)
  pval
}

eye.linear<-function(x=1:length(y),y,tau)
{
  f<-y
  ## x??????
  length(f)
  
  ## y??????
  max(f)-min(f)
  
  ## scaling
  c<-length(f)/(max(f)-min(f))
  cf<-f*c
  
  ## length 
  len<-c()
  for (i in 1:(length(cf)-1))
  {
    len[i]<-sqrt(1+(cf[i+1]-cf[i])^2)
  }
  
  ## total length
  sum(len)
  
  ## ?ո????? ??��????
  sum(len)/length(cf)
  
  ## ??��???? -> ??ȯ?ϰ??? ?ϴ? ??ǥ???ϱ? (linear interpolation)
  ## ??ȯ?ϰ??? ?ϴ? x??ǥ
  unit.length<-sum(len)/length(cf) #??��????
  len.cumulated<-len*0 #????????
  
  len.cumulated[1]<-len[1]
  for(i in 2:length(len))
  {
    len.cumulated[i]<-len.cumulated[i-1]+len[i]
  }
  
  n1<-sum(len)/unit.length #??��???̿? ???? ?? sampling ??
  
  unit.axis<-unit.length*(1:n1) #??��??????
  x.axis<-c()
  for(i in 1:n1)
  {
    if(unit.axis[i]>=len.cumulated[length(len.cumulated)]){
      kkk<-n1 
    }else{
      kkk<-min(which(unit.axis[i]<len.cumulated))
    }
    
    if(kkk==1)
    {
      x.axis[i]<-unit.axis[i]/len.cumulated[kkk]
    }else if(kkk==n1){
      x.axis[i]<-n1
    }else{
      x.axis[i]<-kkk+(unit.axis[i]-len.cumulated[kkk-1])/(len.cumulated[kkk]-len.cumulated[kkk-1])
    }
  }
  x.axis # ??ȯ?ϰ??? ?ϴ? x??ǥ 
  ## ??ȯ?ϰ??? ?ϴ? x??ǥ?κ??? ??ȯ?ϰ??? ?ϴ? y??ǥ???ϱ? 
  x.axis.floor<-floor(x.axis) #???? 
  x.axis.ceiling<-ceiling(x.axis) #?ø? 
  y.axis<-c()
  for(i in 1:length(x.axis))
  {
    if(x.axis.floor[i]==0)
    {
      y.axis[i]<-cf[x.axis.ceiling[i]]*x.axis[i]
    }else if(x.axis.floor[i]>=length(y)){
      y.axis[i]<-cf[length(cf)]
    }else {
      y.axis[i]<-
        cf[x.axis.floor[i]]+
        (cf[x.axis.ceiling[i]]-cf[x.axis.floor[i]])*(x.axis[i]-x.axis.floor[i])
    }
  }
  list(x=x.axis,y=ebt(y.axis/c,tau=tau)$M)
}


eye.cubic<-function(x=1:length(y),y,spar=spar,all.knots=TRUE,df=3)
{
  f<-y
  length(f)
  max(f)-min(f)
  
  ## scaling
  c<-length(f)/(max(f)-min(f))
  cf<-f*c
  
  ## length 
  len<-c()
  for (i in 1:(length(cf)-1))
  {
    len[i]<-sqrt(1+(cf[i+1]-cf[i])^2)
  }
  
  ## total length
  sum(len)
  
  sum(len)/length(cf)
  
  unit.length<-sum(len)/length(cf)
  len.cumulated<-len*0
  
  len.cumulated[1]<-len[1]
  for(i in 2:length(len))
  {
    len.cumulated[i]<-len.cumulated[i-1]+len[i]
  }
  
  n1<-sum(len)/unit.length 
  
  unit.axis<-unit.length*(1:n1)
  x.axis<-c()
  for(i in 1:n1)
  {
    if(unit.axis[i]>=len.cumulated[length(len.cumulated)]){
      kkk<-n1 
    }else{
      kkk<-min(which(unit.axis[i]<len.cumulated))
    }
    
    if(kkk==1)
    {
      x.axis[i]<-unit.axis[i]/len.cumulated[kkk]
    }else if(kkk==n1){
      x.axis[i]<-n1
    }else{
      x.axis[i]<-kkk+(unit.axis[i]-len.cumulated[kkk-1])/(len.cumulated[kkk]-len.cumulated[kkk-1])
    }
  }
  x.axis 
  
  sp<-smooth.spline(x,cf,spar=spar,all.knots = all.knots)
  y.axis<-predict(sp,x.axis)$y
  list(x=x.axis,y=y.axis/c)
}

sdtresh<-function(f,sdtresh)
{
  for (i in 1:length(f))
  {
    if(f[i]<3*sdtresh) f[i]<-sdtresh
  }
  f
}

# linear imputation 
fun.initial <- function(x, y, mindex) {
  
  ndata <- length(y)    
  mindex <- sort(mindex)
  oindex <- (1:ndata)[-mindex]
  
  if (mindex[1] == 1) {
    y[1] <- y[min(oindex)]; mindex <- mindex[-1]; oindex <- c(1, oindex)
  }
  if (mindex[length(mindex)] == ndata) {
    y[ndata] <- y[max(oindex)]; mindex <- mindex[-length(mindex)]; oindex <- c(oindex, ndata)
  }
  
  nmindex <- length(mindex)
  noindex <- length(oindex)
  
  tmp <- NULL
  for(i in 1:nmindex) {
    tmpoindex <- seq(1:(noindex-1))[oindex[1:(noindex-1)] < mindex[i] & mindex[i] <= oindex[2:noindex]]
    tmpoindex <- c(tmpoindex, tmpoindex+1)
    tmp <- c(tmp, diff(y[oindex[tmpoindex]]) / diff(x[oindex[tmpoindex]]) * 
               (x[mindex[i]]- x[oindex[tmpoindex[2]]]) + y[oindex[tmpoindex[2]]])
  }
  y[mindex] <- tmp
  
  list(yin=y)
}


covar<-function(t1,t2)
{ 
  u1<-t1[1]; v1<-t1[2]; u2<-t2[1]; v2<-t2[2];
  v1*v2/(u1*u2)*min(u1,u2)-v2/u2*min(v1,u2)-v1/u1*min(u1,v2)+min(v1,v2)
}

generate.gaussian<-function(Tx,nsim=Tx)
{
  sim<-list()
  for(k in 1:nsim)
  {
    distH<-rep(0,Tx*Tx); dim(distH)<-c(Tx,Tx)
    for(i in 1:(Tx-1))
    {
      t1<-c(i/Tx,1/Tx); t2<-c((i+1)/Tx,1/Tx);
      Sigma_11<-abs(covar(t1,t1));Sigma_12<-abs(covar(t1,t2));Sigma_21<-abs(covar(t2,t1)); Sigma_22<-abs(covar(t2,t2));
      
      if(Sigma_22==0) mu_cond<-0 else mu_cond<-Sigma_12/Sigma_22*distH[t1[1]*Tx,t1[2]*Tx]
      if(Sigma_11==0) Sigma_cond<-Sigma_22 else Sigma_cond<-abs(Sigma_22-Sigma_12^2/Sigma_11)
      distH[t2[1]*Tx,t2[2]*Tx]<-rnorm(1,mean=mu_cond,sd=sqrt(Sigma_cond));
      
      for(j in 1:(Tx-1))
      {
        t1<-c(i/Tx,j/Tx); t2<-c(i/Tx,(j+1)/Tx);
        Sigma_11<-abs(covar(t1,t1));Sigma_12<-abs(covar(t1,t2));Sigma_21<-abs(covar(t2,t1)); Sigma_22<-abs(covar(t2,t2)); 
        
        if(Sigma_22==0) mu_cond<-0 else mu_cond<-Sigma_12/Sigma_22*distH[t1[1]*Tx,t1[2]*Tx]
        if(Sigma_11==0) Sigma_cond<-Sigma_22 else Sigma_cond<-abs(Sigma_22-Sigma_12^2/Sigma_11)
        distH[t2[1]*Tx,t2[2]*Tx]<-rnorm(1,mean=mu_cond,sd=sqrt(Sigma_cond));
        if(i>j) distH[t2[1]*Tx,t2[2]*Tx]<-0
      }
    }
    sim[[k]]<-distH
    sim[[k]][,1]<-0
  }
  sim
}

qsc<-function(uT,alpha=0.05,sim.result)
{
  qdist<-c()
  for(i in 1:100) qdist[i]<-max(sim.result[[i]][1:uT,1:uT])
  quantile(qdist,prob=0.95)
}


Dhat<-function(x,uT,vT)
{
  Tx<-length(x)
  Dhat<- sum(x[1:vT])/Tx - sum(x[1:uT])/Tx * (vT/uT)
  Dhat
}

calDhat<-function(x,uT)
{
  calDhat<-c()
  for(j in 1:uT)
  {
    calDhat[j]<-abs(Dhat(x,uT,j))
  }
  max(calDhat)
}

lrvar<-function(x,meanx){
  mm<-round(log(length(x)),0)
  xccf<-ccf(x-meanx,x-meanx,type="covariance",plot=F)$acf #xccf<-ccf(x-smooth.spline(x,spar=0.2)$y,x-smooth.spline(x,spar=0.2)$y,type="covariance",plot=F)$acf
  ii<-which(xccf==max(xccf))
  lrvariance<-xccf[ii]+2*sum(xccf[ii:(ii+mm)])
  lrvariance
}

cp.find<-function(x,sim.result)
{
  Tx<-length(x)
  indicator<-c()
  q<-c()
  calDhatsc<-c()
  calHhatsc<-c()
  
  mm<-round(log(length(x)),0)
  xccf<-ccf(x-smooth.spline(x,spar=0.2)$y,x-smooth.spline(x,spar=0.2)$y,type="covariance",plot=F)$acf
  ii<-which(xccf==max(xccf))
  lrvariance<-xccf[ii]+2*sum(xccf[ii:(ii+mm)])

  for(i in 1:Tx)
  {
    q[i]<-qsc(uT=i,alpha=0.05,sim.result)
    calDhatsc[i]<-calDhat(x,uT=i)/sqrt(lrvariance)
    calHhatsc[i]<-sqrt(Tx)*(calDhatsc[i]-0)
    indicator[i]<-(calHhatsc[i]<q[i])
  }
  Tx<-sum(indicator)
  
  if(Tx<10){
    t_0hat<-Tx
  }else{
    xTx<-x[1:Tx]
    indicator<-c()
    calDhatsc<-c()
    calHhatsc<-c()
    library(sandwich)
    mm<-round(log(length(xTx)),0)  
    xccf<-ccf(xTx-smooth.spline(xTx,spar=0.2)$y,xTx-smooth.spline(xTx,spar=0.2)$y,type="covariance",plot=F)$acf
    ii<-which(xccf==max(xccf))
    lrvariance<-xccf[ii]+2*sum(xccf[ii:(ii+mm)])
    for(i in 1:Tx)
    {
      calDhatsc[i]<-calDhat(xTx,uT=i)/lrvariance
      calHhatsc[i]<-sqrt(Tx)*(calDhatsc[i]-0)
      indicator[i]<-(calHhatsc[i]<q[i])
    }
    t_0hat<-sum(indicator)
  }
  t_0hat
}


cp.find2<-function(x,sim.result,calDhat_mean_tau)
{
  Tx<-length(x)
  indicator<-c()
  q<-c()
  calDhatsc<-c()
  calHhatsc<-c()
  calDhat_sample<-c()
  for(i in 1:Tx)
  {
    q[i]<-qsc(uT=i,alpha=0.05,sim.result)
    calDhat_sample[i]<-calDhat(x,uT=i)
  }
  calDhatsc<-(calDhat_sample)/sd(calDhat_sample-calDhat_mean_tau)
  calHhatsc<-sqrt(Tx)*(calDhatsc)
  indicator<-(calHhatsc<q)
  Tx<-sum(indicator)
  
  if(Tx<10){
    t_0hat<-Tx
  }else{
    calDhat_sample2<-c()
    xTx<-x[1:Tx]
    indicator<-c()
    calDhatsc<-c()
    calHhatsc<-c()
    for(i in 1:Tx)
    {
      calDhat_sample2[i]<-calDhat(x,uT=i)
    }
    calDhatsc<-(calDhat_sample2)/sd(calDhat_sample-calDhat_mean_tau)
    calHhatsc<-sqrt(Tx)*(calDhatsc)
    indicator<-(calHhatsc<q[1:Tx])
    t_0hat<-sum(indicator)
  }
  t_0hat
}

eye.ebt<-function(x=1:length(y),y,tau_df=15,tau_var=15)
{
  f<-y
  # ## variace estimation 
  # tau<-tau_var
  # k<-2*sum((1:(tau-1))^2)/(tau^4)+(1-(1/tau))^2
  # noise2<-extract.hfreq(t,f,tau=tau,iter=0)
  # renoise2<-noise2*sqrt(1/k)
  # 
  ## denoising
  tau<-tau_df
  len<-ebt(f,tau=tau,V="volume")$V
  sum(len)
  # reasonable unit length
  sum(len)/length(f)
  
  unit.length<-sum(len)/length(f) 
  len.cumulated<-len*0     
  
  len.cumulated[1]<-len[1]
  for(i in 2:length(len))
  {
    len.cumulated[i]<-len.cumulated[i-1]+len[i]
  }
  
  n1<-sum(len)/unit.length 
  
  unit.axis<-unit.length*(1:n1) 
  x.axis<-c()
  for(i in 1:n1)
  {
    if(unit.axis[i]>=len.cumulated[length(len.cumulated)]){
      kkk<-n1 
    }else{
      kkk<-min(which(unit.axis[i]<len.cumulated))
    }
    
    if(kkk==1)
    {
      x.axis[i]<-unit.axis[i]/len.cumulated[kkk]
    }else if(kkk==n1){
      x.axis[i]<-n1
    }else{
      x.axis[i]<-kkk+(unit.axis[i]-len.cumulated[kkk-1])/(len.cumulated[kkk]-len.cumulated[kkk-1])
    }
  }
  x.axis 
  x.axis.floor<-floor(x.axis) 
  x.axis.ceiling<-ceiling(x.axis)
  y.axis<-c()
  for(i in 1:length(x.axis))
  {
    if(x.axis.floor[i]==0)
    {
      y.axis[i]<-f[x.axis.ceiling[i]]*x.axis[i]
    }else if(x.axis.floor[i]>=length(y)){
      y.axis[i]<-f[length(f)]
    }else {
      y.axis[i]<-
        f[x.axis.floor[i]]+
        (f[x.axis.ceiling[i]]-f[x.axis.floor[i]])*(x.axis[i]-x.axis.floor[i])
    }
  }
  ebtresult<-ebt(y.axis,tau=tau,V="volume",M="mean")
  list(x=x.axis,yL=ebtresult$L,yU=ebtresult$U,yM=ebtresult$M,y=y.axis)
}

df.last<-function(y,tau_df=10,tau_var=2,LUM="M")
{
  yiter<-y
  xiter<-1:length(y)
  temp.eye<-eye.ebt(y=yiter,tau_df=tau_df,tau_var=tau_var)
  if(LUM=="L") yiter<-temp.eye$yL
  else if(LUM=="U") yiter<-temp.eye$yU
  else yiter<-temp.eye$yM
  true<-temp.eye$x
  down<-floor(temp.eye$x)
  down[which(down==0)]=1
  up<-ceiling(temp.eye$x)
  up[which(up>length(yiter))]<-length(yiter)
  xiter<-xiter[down]+(true-down)*(xiter[up]-xiter[down])
  newxiter<-c(xiter,1:length(xiter))
  newyiter<-c(yiter,rep(NA,length(yiter)))
  newdf<-data.frame(newxiter,newyiter)
  newdf2<-newdf[order(newdf$newxiter),]
  mindex<-which(is.na(newdf2$newyiter)==TRUE)
  predy<-fun.initial(newdf2$newxiter,newdf2$newyiter,mindex=mindex)$yin[mindex]
  predy
}

ads<-function(y,taudf,iter=1)
{
  j<-0
  dy<-df.last(y=y,tau_df=taudf)
  if(iter==0)
  {
    while((max(abs(dyold-dy))>0.01)&(j<20))
    {
      dyold<-dy
      dy<-df.last((dyold+(taudf-1)*y)/taudf,tau_df=taudf)
      #dy<-df.last(-(dyold-(taudf-1)*y)/taudf+y,tau_df=taudf)
      j<-j+1
    }
  }
  for(l in 1:iter)
  {
    dyold<-dy
    dy<-df.last((dyold+(taudf-1)*y)/taudf,tau_df=taudf)
  }
  dy
}

ads2<-function(y,taudf)
{
  j<-0
  dy<-ebt(y,tau=taudf)$M 
  #dy<-y-extract.hfreq(f=y,tau=taudf)
  dyold<-y
  while((max(abs(dyold-dy))>0.01)&(j<10))
  {
    dyold<-dy
    dy<-ebt((dyold+(taudf-1)*y)/taudf,tau=taudf)$M
    #dy<-dyold-extract.hfreq(f=(dyold+(taudf-1)*y)/taudf,tau=taudf)
    j<-j+1
  }
  dy
}


cvlambda<-function(y,lambda=5:25)
{
  cv<-rep(0,length(lambda))
  for(l in 1:length(lambda))
  {
    taudf<-lambda[l]
    index<-1:length(y)
    mindex1<-seq(2,N,by=2)
    mindex2<-index[-mindex1]
    
    yin1<-fun.initial(index,y=y,mindex=mindex1)$yin
    yin2<-fun.initial(index,y=y,mindex=mindex2)$yin
    
    dyin1<-ads(yin1,taudf=taudf)
    dyin2<-ads(yin2,taudf=taudf)
    cv[l]<-(sum((y[mindex1]-dyin1[mindex1])^2)+sum((y[mindex2]-dyin2[mindex2])^2))/N
    cat(round(l/length(lambda)*100,2),"%\n")
  }
  cv
}


cvlambda2<-function(y,lambda=5:25)
{
  cv<-rep(0,length(lambda))
  for(l in 1:length(lambda))
  {
    taudf<-lambda[l]
    index<-1:length(y)
    mindex1<-seq(2,N,by=2)
    mindex2<-index[-mindex1]
    yeven<-y[mindex1]
    yodd<-y[mindex2]
    dyin1<-y*0; dyin2<-y*0;
    dyin1[mindex1]<-ads(yeven,taudf=taudf)
    dyin2[mindex2]<-ads(yodd,taudf=taudf)
    
    yin1<-fun.initial(index,y=dyin1,mindex=mindex2)$yin
    yin2<-fun.initial(index,y=dyin2,mindex=mindex1)$yin

    cv[l]<-(sum((y[mindex1]-yin1[mindex2])^2)+sum((y[mindex2]-yin2[mindex1])^2))/N
    #cat(round(l/length(lambda)*100,2),"%\n")
  }
  cv
}


tauest<-function(f,taulist=15:25,V="volume")
{
  resultOfEBT.V1<-rep(0,length(f)*length(taulist)); dim(resultOfEBT.V1)<-c(length(f),length(taulist))
  resultOfEBT.V<-rep(0,length(f)*length(taulist)); dim(resultOfEBT.V)<-c(length(f),length(taulist))   
  for(tl in 1:length(taulist))
  {
    resultOfEBT.V1[,tl]<-ebt(f,tau=taulist[tl],V=V)$V
    resultOfEBT.V[,tl]<-resultOfEBT.V1[,tl]
    #resultOfEBT.V[taulist[tl]:(length(resultOfEBT.V1[,tl])-taulist[tl]),tl]<-filter(resultOfEBT.V1[,tl],c(1:taulist[tl],(taulist[tl]-1):1)/taulist[tl]^2, sides=2)[taulist[tl]:(length(resultOfEBT.V1[,tl])-taulist[tl])]
    resultOfEBT.V[taulist[tl]:(length(resultOfEBT.V1[,tl])-taulist[tl]),tl]<-filter(resultOfEBT.V1[,tl],c(0,0,1,0,0), sides=2)[taulist[tl]:(length(resultOfEBT.V1[,tl])-taulist[tl])]
  }
  tau.est<-c()
  for(tindex in 1:length(f))
  {
    # x2<-taulist[min(which(resultOfEBT.V[tindex,]==max(resultOfEBT.V[tindex,])))]
    # if (x2==min(taulist)) {
    #   tau.est[tindex]<-x2
    # } else if (x2==max(taulist)) {
    #   tau.est[tindex]<-x2
    # } else{
    # x1<-x2-1
    # x3<-x2+1
    # y2<-resultOfEBT.V[tindex,x2-min(taulist)+1]
    # y1<-resultOfEBT.V[tindex,x1-min(taulist)+1]
    # y3<-resultOfEBT.V[tindex,x3-min(taulist)+1]
    # tau.est[tindex]<-minconvex(p1=c(x1,y1),p2=c(x2,y2),p3=c(x3,y3))
    # }
    tau.est[tindex]<-taulist[min(which(resultOfEBT.V[tindex,]==max(resultOfEBT.V[tindex,])))]
  }
  tau.est
}

decomp.convex<-function(f)
{
  fp<-f*0 #fp is f+
  fp[1]<-f[1]
  for(i in 2:length(f))
  {
    cond<-f[i]-f[i-1]<0
    fp[i]<-(1-cond)*(fp[i-1]+f[i]-f[i-1])+cond*fp[i-1]
  }
  fm<-f-fp #fm is f-
  
  fpp<-f*0 #fpp is f++ 
  fpp[1]<-fp[1]
  fpp[2]<-fp[2]
  for(i in 3:length(f))
  {
    if(fp[i-1]-fp[i-2]==0) cond<-0 else cond<-(fp[i]-fp[i-1])/(fp[i-1]-fp[i-2])<1
    fpp[i]<-(1-cond)*(fpp[i-1]+fpp[i-1]-fpp[i-2]+fp[i]-fp[i-1])+cond*(fpp[i-1]+fpp[i-1]-fpp[i-2])
  }
  fpm<-fp-fpp
  
  fm<- -fm
  fmp<-f*0 #fpp is f++ 
  fmp[1]<- fm[1]
  fmp[2]<- fm[2]
  for(i in 3:length(f))
  {
    if(fm[i-1]-fm[i-2]==0) cond<-0 else cond<-(fm[i]-fm[i-1])/(fm[i-1]-fm[i-2])<1
    fmp[i]<-(1-cond)*(fmp[i-1]+fmp[i-1]-fmp[i-2]+fm[i]-fm[i-1])+cond*(fmp[i-1]+fmp[i-1]-fmp[i-2])
  }
  fmm<-fm-fmp
  fmm<- -fmm #convex
  fmp<- -fmp #concave
  convex=fmm+fpp
}

minconvex<-function(p1,p2,p3)
{
  x1<-p1[1]; y1<-p1[2];
  x2<-p2[1]; y2<-p2[2];
  x3<-p3[1]; y3<-p3[2];
  sol<-solve(cbind(c(x1^2,x2^2,x3^2),c(x1,x2,x3),c(1,1,1)))%*%c(y1,y2,y3)
  a<-sol[1]
  b<-sol[2]
  c<-sol[3]
  -b/2/a
}

rhomap<-function(f,g,maxtau,M="mvolume",V="volume",mom=1)
{
  f<-(f-mean(f))/sd(f)
  g<-(g-mean(g))/sd(g)
  rhoM<-rep(0,length(f)*maxtau); dim(rhoM)<-c(length(f),maxtau);
  febt<-list()
  gebt<-list()
  for(tau in 10:maxtau)
  {
    febt[[tau]]<-ebt(f,tau=tau,V="volume")
    gebt[[tau]]<-ebt(g,tau=tau,V="volume")
    for(i in 1:length(f))
    {
      set1<-febt[[tau]]$band[i,]
      set2<-gebt[[tau]]$band[i,]  
      DV <- c(set1,set2); DV2<-c((set1-mean(set1))^2,(set2-mean(set2))^2)
      IV <- factor(rep(c("A", "B"), c(length(set1), length(set2))))
      library(coin)  # for oneway_test(), pvalue()
      if(mom==1){rhoM[i,tau]<-pvalue(oneway_test(DV ~ IV, alternative="greater",distribution=approximate(B=9999)))}else if(mom==2){
      rhoM[i,tau]<-pvalue(oneway_test(DV ~ IV, alternative="greater",distribution=approximate(B=9999)))*
        pvalue(oneway_test(DV2 ~ IV, alternative="greater",distribution=approximate(B=9999)))}
    }
  }  
  rhoM
}

stand<-function(f)
{
  (f-mean(f))/sd(f)
}


Tmap<-function(f,g,maxtau,M="mvolume",V="volume",mom=1)
{
  f<-(f-mean(f))/sd(f)
  g<-(g-mean(g))/sd(g)
  rhoM<-rep(0,length(f)*maxtau); dim(rhoM)<-c(length(f),maxtau);
  febt<-list()
  gebt<-list()
  for(tau in 10:maxtau)
  {
    febt[[tau]]<-ebt(f,tau=tau,V="volume")
    gebt[[tau]]<-ebt(g,tau=tau,V="volume")
    for(i in 1:length(f))
    {
      set1<-febt[[tau]]$band[i,]
      set2<-gebt[[tau]]$band[i,]  
      if(mom==1){
        rhoM[i,tau]<-t.test(x=set1,y=set2)$p.value
      }else if(mom==2){
        rhoM[i,tau]<-t.test(x=set1,y=set2)$p.value*t.test(x=(set1-mean(set1))^2,y=(set2-mean(set2))^2)$p.value
      }
      
    }
  }  
  rhoM
}

ampanalty<-function(f,tau,lambda=10,k=0.5)
{
  len<-length(mirror(f))
  ltemp<-ceiling(length(f)/2)
  index<-(ltemp-2):(len-ltemp-1)
  mf<-mirror(f)
  V<-ebt(mf,tau=tau)$V[index]
  V2prime<-abs(diff(diff(V))) #V2prime : changes of volume
  index<-exp(-V2prime/lambda)>k 
  # Large V2prime => Large index
  # small k => Large index 
  # small lambda => Large index 
  B<-index
  B
}
  
extfast<-function(f,tau)
{
  index<-seq(from=1,to=length(f),by=length(f)/tau)
  index2<-sort(c(floor(index)[-1],ceiling(index)))
  fft_f<-fft(f)
  fft_f[-index]<-0
  fhat<-fft(fft_f,inverse = T)/Len
  Re(fhat)
}

bext<-function(t,f,tau,M="mean",V="sd",interpolation="linear",tol=0.05,iter=100,fzerotol=0.05)
{
  m.extracted.old<-f
  m.extracted<-ebt(f,tau=tau)$M
  res<-f-m.extracted
  j<-0
  Aindex_1<-rep(1,length(f))
  VV<-ebt(res,tau=tau,V="volume")$V
  V2prime<-c(0,diff(c(0,diff(VV))))
  Aindex_2<-abs(V2prime)<0.05
  Zindex<-(VV<0.0005)*(abs(V2prime)<0.05)
  while((max(abs(m.extracted.old-m.extracted))>tol)&(j<iter))
  {
    if(sum(Aindex_2)==0){
      res<-res
    }else{
      Aindex_2<-Aindex_2*Aindex_1*(1-Zindex)
      Aindex<-which(Aindex_2==1)
      m.extracted.old<-m.extracted
      m.extracted<-ebt(res,tau=tau)$M
      res[Aindex]<-(res-m.extracted)[Aindex]
      j<-j+1
      Aindex_1<-Aindex_2
      VV<-ebt(res,tau=tau,V="volume")$V
      V2prime<-c(0,diff(c(0,diff(VV))))
      Aindex_2<-abs((V2prime))<0.05
      Zindex<-(VV<0.0005)*(abs(V2prime)<0.05)
    } 
  }
  hfreq<-res
  list(hfreq=hfreq,Aindex=Aindex)
}
