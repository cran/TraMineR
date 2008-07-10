## ========================================
## Get and set length of seqe
## ========================================

seqelength<-function(s){
  seqelength.internal<-function(s){
    if(is.seqe(s))
      return(.Call("tmrsequencegetlength",s, PACKAGE="TraMineR"))

    return(-1)
  }
  if(is.seqelist(s)){
      as.numeric(sapply(unlist(s),seqelength.internal))
  }else if(is.seqe(s)){
    as.numeric(seqelength.internal(s))
  }else{
    stop("s should be a seqelist. See help on seqecreate.")
  }
}
seqesetlength<-function(s, len){
  if(!is.seqelist(s))stop("s should be a seqelist. See help on seqecreate.")
  if(length(s)!=length(len))stop("s and len should be of the same size.")
  return(invisible(.Call("tmrsequencesetlength",s,as.double(len), PACKAGE="TraMineR")))
}

