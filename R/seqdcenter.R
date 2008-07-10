############################
## Compute distance to center for a group
############################

seqdcenter<-function(distmatrix,group=NULL){
  if(is.null(group)){
    group<-integer(nrow(distmatrix))
    group[]<-1
  }
  grptbl=table(group)
  ret<-numeric(nrow(distmatrix))
  ind<-1:nrow(distmatrix)
  grp<-as.factor(group)
  lgrp<-levels(grp)
  #pour chaque valeur du groupe
  for(i in 1:length(lgrp)){
    #on crée le groupe en question
    grpindiv<-ind[grp==lgrp[i]]
    #on calcul l'inertie intraclasse
    ret[grpindiv]<-.Call("tmrinertiacontrib", distmatrix, as.integer(grpindiv),PACKAGE="TraMineR")
  }
  return(ret)
}