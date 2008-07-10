seqcmpinertia<-function(distmatrix,group){
  ret<-list()
  ret$stat<-list()
  #Notation comme pour l'ANOVA, SC=Inertia dans le sens du criète de Ward
  ret$stat$SCtot<-.Call("tmrsubmatrixinertia",distmatrix,1:nrow(distmatrix),PACKAGE="TraMineR")
  ret$stat$SCres=0
  ind<-1:nrow(distmatrix)
  grp<-as.factor(group)
  lgrp<-levels(grp)
  #pour chaque valeur du groupe
  ret$variance<-vector("list",length(lgrp))
  #ret$contrib<-vector("list",length(lgrp))
  for(i in 1:length(lgrp)){
    #on crée le groupe en question
    grpindiv<-ind[grp==lgrp[i]]
    if(length(grpindiv)>1){
      r<-.Call("tmrsubmatrixinertia", distmatrix, as.integer(grpindiv), PACKAGE="TraMineR")
      ret$variance[i]<-r/length(grpindiv)
      ret$stat$SCres<-ret$stat$SCres+r
    }
    else{
      ret$variance[i]<-0
    }



  }
  ret$stat$SCexp<-ret$stat$SCtot-ret$stat$SCres
  ret$stat$PseudoR2<-ret$stat$SCexp/ret$stat$SCtot
  n<-length(ind)
  k<-length(lgrp)
  ret$stat$n<-n
  ret$stat$numgrp<-k
  ret$stat$PseudoBIC<-n*log(ret$stat$SCres/n)+k*log(n)
  ret$stat$PseudoAICu<-log(ret$stat$SCres/(n-k))+(n+k)/(n-k-2)
  ret$stat$PseudoF<-(ret$stat$SCexp/(k-1))/(ret$stat$SCres/(n-k))
  ret$stat$Variance<-ret$stat$SCtot/n
  ret$stat<-as.data.frame(ret$stat)
  ret$variance<-as.data.frame(ret$variance)
#  ret$contrib<-as.data.frame(ret$contrib)
  colnames(ret$variance)<-lgrp
#  colnames(ret$contrib)<-lgrp
  return(ret)
}