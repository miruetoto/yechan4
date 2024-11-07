lin_impute<-function(t,y,mindex)
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


left_const<-function(y,mindex)
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


right_const<-function(y,mindex)
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
#---#

ebt<-function(t,f=NULL,tau,mfunc="mean",vfunc="var",inter_method="linear")
{
  if (is.null(f)){
    f <- t 
    t <- 1:length(f)
  }
  tsave<-t
  fsave<-f
  f<-c(rep(f[1],tau*2),f,rep(f[length(f)],tau*2))
  t<-1:length(f)
  len<-length(f)
  sampled_index<-list()
  missing_index<-list()
  band<-rep(0,len*tau); dim(band)<-c(len,tau)
  if (inter_method=="cubic")
  {
    for (eta in 1:tau){
      sampled_index[[eta]]<-seq(from=eta,to=len,by=tau)
      if (sampled_index[[eta]][1]!=1) sampled_index[[eta]]=c(1,sampled_index[[eta]])
      if (sampled_index[[eta]][length(sampled_index[[eta]])]!=len) sampled_index[[eta]]=c(sampled_index[[eta]],len)
      missing_index[[eta]]<-(1:len)[-sampled_index[[eta]]]
      result<-smooth.spline(t[sampled_index[[eta]]],f[sampled_index[[eta]]],spar=0,all.knots=TRUE)
      band[,eta][sampled_index[[eta]]]<-result$y
      band[,eta][missing_index[[eta]]]<-predict(result,t[missing_index[[eta]]])$y
    }
  }else if (inter_method=="linear"){
    for (eta in 1:tau){
      sampled_index[[eta]]<-seq(from=eta,to=len,by=tau)
      if (sampled_index[[eta]][1]!=1) sampled_index[[eta]]<-c(1,sampled_index[[eta]])
      if (sampled_index[[eta]][length(sampled_index[[eta]])]!=len) sampled_index[[eta]]<-c(sampled_index[[eta]],len)
      missing_index[[eta]]<-(1:len)[-sampled_index[[eta]]]
      band[,eta]<-lin_impute(t,f,missing_index[[eta]])
    }
  }
  U<-apply(band,1,max)    
  L<-apply(band,1,min)


  M1<-apply(band,1,mean)
  M2<-apply(band,1,median)
  M3<-(L+U)/2
  V1=U-L
  V2=apply(band,1,var) ; V2=V2*(tau-1)/tau;

  if (mfunc=="mean") M<-M1
  else if (mfunc=="median") M<-M2
  else if (mfunc=="volume") M<-M3

  if (vfunc=="volume"){
    V<-V1
    L<-L
    U<-U
  }
  else if (vfunc=="var"){
    V<-V2
    L<-M-V
    U<-M+V
  }

  index<-(2*tau+1):(length(f)-2*tau)
  sampled_index <- lapply(sampled_index, function(x) {
  filtered <- x[x <= length(fsave)]
  return(filtered)
  })    
    
  list(t=tsave,f=fsave,V=V[index],L=L[index],U=U[index],M=M[index],tau=tau,band=band[index,],sampled_index=sampled_index)
}

mvmap<-function(t,f=NULL,maxtau,M="mean",V="var",inter_method="linear")
{
  VM<-rep(0,length(f)*maxtau); dim(VM)<-c(length(f),maxtau);
  MM<-VM 
  for(tau in 2:maxtau)
  {
    out <- ebt(t,f,tau=tau,M,V,inter_method)
    VM[,tau]<-out$V 
    MM[,tau]<-out$M 
  }  
  list(t=t,vmap=VM, mmap=MM)
}
