############################
## Compute distance to center for a group
############################

disscenter<-function(diss,group=NULL, medoids.index=FALSE, max.iter=20){
  trim <- 0
  if(inherits(diss, "dist")){
    diss<-as.matrix(diss)
  }
  if(is.null(group)){
    group<-integer(nrow(diss))
    group[]<-1
  }
  grptbl=table(group)
  ret<-numeric(nrow(diss))
  ind<-1:nrow(diss)
  grp<-factor(group)
  lgrp<-levels(grp)
  medoids<-numeric(length(lgrp))
  keep=1-trim
  #pour chaque valeur du groupe
  for(i in 1:length(lgrp)){
    #on crée le groupe en question
    cond<-grp==lgrp[i]
    grpindiv<-sort(ind[cond])
    #on calcul l'inertie intraclasse
    r<-(.Call("tmrsubmatrixinertiaCindividuals", diss, as.integer(grpindiv-1),PACKAGE="TraMineR")/length(grpindiv))
    dc <- .Call("tmrinertiacontrib", diss, as.integer(grpindiv),PACKAGE="TraMineR")-r
    ret[grpindiv]<-dc
    medoids[i] <-sort(which(ret==min(dc)&cond))[1]
    if(trim>0){
      m2<- -1
      iter<-1
      while(iter<max.iter&&m2!=medoids[i]){
        dcm<- diss[grpindiv,medoids[i]]
        maxdist<-quantile(dcm, probs=keep)
        m2 <- medoids[i]
        trimmedgrpindiv<-sort(grpindiv[dcm<=maxdist])
        r<-(.Call("tmrsubmatrixinertiaCindividuals", diss, as.integer(trimmedgrpindiv-1),PACKAGE="TraMineR")/length(trimmedgrpindiv))
        dc <- .Call("tmrinertiacontrib", diss, as.integer(trimmedgrpindiv),PACKAGE="TraMineR")-r
        ret[grpindiv]<-NA
        ret[trimmedgrpindiv]<-dc
        medoids[i] <-sort(which(!is.na(ret)&ret==min(dc)&cond))[1]
        iter<- iter+1
      }
    }
  }
  names(medoids)<-lgrp
  if(medoids.index)return(medoids)
  return(ret)
}