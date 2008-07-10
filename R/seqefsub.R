## ========================================
## Find frequent subsequences
## ========================================

seqefsub<-function(seq,minSupport=NULL, pMinSupport=NULL,maxGap=-1, windowSize=-1, ageMin=-1, ageMax=-1,ageMaxEnd=-1, maxK=-1){
#  message("Event sequence analysis module is still experimental")
  if(!is.seqelist(seq))stop("seq should be a seqelist. See help on seqecreate.")
  if(is.null(minSupport)){
    if(is.null(pMinSupport)){
      stop("You should specify a minimum support through minSupport or pMinSupport argument")
    }
    minSupport<-round(pMinSupport*length(seq))
  }
  classname<-c("seqe")
  subseq<-.Call("tmrfindsubsequences",unlist(list(seq)),as.double(c(maxGap)),as.double(c(windowSize)),
    as.double(c(ageMin)),as.double(c(ageMax)),as.double(c(ageMaxEnd)),as.integer(c(minSupport)),as.integer(c(maxK)),classname,PACKAGE="TraMineR")
  ord<-order(unlist(subseq[1]),decreasing=TRUE)
  ret<-list();
  ret$subseq<-unlist(subseq[2])[ord]
  class(ret$subseq)<-c("seqelist","list")
#  attr(ret$subseq,"dictionnary")<-attr(seq,"dictionnary")
  ret$support<-unlist(subseq[1])[ord]
  return(ret)
}