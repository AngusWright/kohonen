generate.kohgroups<-function(som,factor.nbins=Inf,ncores=12) { #{{{
  somclust<-som$unit.classif
  som.dim<-c(som$grid$xdim,som$grid$ydim)
  #Convert the IDs from data-to-SOMcell into data-to-cluster {{{
  if (factor.nbins<som.dim[1]*som.dim[2]) { 
    som.hc = cutree(hclust(dist(som$codes[[1]])), factor.nbins)
    require(foreach)
    require(doParallel)
    registerDoParallel(cores=ncores)
    somind<-foreach(i=1:factor.nbins,.combine=rbind,.inorder=FALSE)%dopar%{ 
      #Assign the data to the cluster if it's SOMcell belongs to the cluster
      ind<-which(som$unit.classif%in%which(som.hc==i))
      return=cbind(ind,rep(i,length(ind)))
    } 
    somclust[somind[,1]]<-somind[,2]
  } else { 
    som.hc<-1:(som.dim[1]*som.dim[2])
  }
  #}}}
  return=list(data.cl=somclust,som.cl=som.hc)
}#}}}

generate.kohgroup.property<-function(property,group,factor.nbins,na.rm=TRUE) { #{{{
  lev<-1:factor.nbins
  mean.prop<-sd.prop<-med.prop<-mad.prop<-range.prop<-cl.num<-rep(NA,length(lev))
  pb<-txtProgressBar(min=0,max=length(lev),style=3) 
  for (i in lev) { 
    setTxtProgressBar(pb,which(lev==i)) 
    ind<-which(group==i)
    cl.num[which(lev==i)]<-length(ind)
    mean.prop[which(lev==i)]<-mean(property[ind],na.rm=na.rm)
    sd.prop[which(lev==i)]<-sd(property[ind],na.rm=na.rm)
    med.prop[which(lev==i)]<-med(property[ind],na.rm=na.rm)
    mad.prop[which(lev==i)]<-mad(property[ind],na.rm=na.rm)
    range.prop[which(lev==i)]<-diff(range(property[ind],na.rm=na.rm))
  }
  close(pb)
  return=data.frame(group=lev,mean=mean.prop,sd=sd.prop,median=med.prop,mad=mad.prop,
                    range.prop=range.prop,num=cl.num)
}#}}}
