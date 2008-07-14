## ========================================
## Create events objects
## ========================================


#tmrsequence<-function(id,timestamp,event){
#  .Call("tmrsequence",as.integer(id),as.double(timestamp),as.integer(event), PACKAGE="TraMineR")
#}
seqecreate<-function(id,timestamp, event, endEvent=NULL){
#  warning("Event sequence analysis module is still experimental", call.=FALSE)
  classname<-c("seqe")
  intEvent=NULL
  if(is.factor(event)){
    dictionnary<-levels(event)
    if(!is.null(endEvent)){

      for(i in 1:length(dictionnary)){
        if(dictionnary[i]==endEvent){
          intEvent=i

        }
      }
      if(is.null(intEvent)){
        stop("endEvent not found in event dictionary")
        return(invisible())
      }
    }
  } else {
    dictionnary<-c()
  }
  ret<-.Call("tmrsequenceseveral",as.integer(id),as.double(timestamp),as.integer(event),as.integer(c(intEvent)),classname,as.character(dictionnary), PACKAGE="TraMineR")
  class(ret)<-c("seqelist","list")
  return(ret)
}
#SEXP tmrsequence(SEXP idpers, SEXP time, SEXP event, SEXP classname, SEXP seq)
#seqecreatesub<-function(seq,timestamp, event){
#  warning("Event sequence analysis module is still experimental", call.=FALSE)
#  if(!is.seqelist(seq))stop("seq should be a seqelist. See help on seqecreate.")
#  classname<-c("seqe")
#  e<-factor(event,levels=levels(seq))
#  ret<-list(.Call("tmrsequence",as.integer(-1),as.double(timestamp),as.integer(e),classname,seq, PACKAGE="TraMineR"))
#  class(ret)<-c("seqelist","list")
#  return(ret)
#}