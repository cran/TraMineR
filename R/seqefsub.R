## ========================================
## Find frequent subsequences
## ========================================

seqefsub<-function(seq,strsubseq=NULL,minSupport=NULL, pMinSupport=NULL,constraint=seqeconstraint(), maxK=-1){
#  message("Event sequence analysis module is still experimental")
  if(!is.seqelist(seq))stop("seq should be a seqelist. See help on seqecreate.")
    if(!inherits(constraint,"seqeconstraint")){
    constraint=seqeconstraint()
    warning("constraint argument should be set using seqeconstraint function. No constraint where used.")
  }

  #We are not looking for frequent but specific subsequences
  if(!is.null(strsubseq)){
    subseq<-seqecreatesub(strsubseq,seq)
    ret<-createsubseqelist(seq,constraint,subseq,data.frame(),type="user")
    support<-colSums(seqeapplysub(ret, method="presence"))
    ord<-order(unlist(support),decreasing=TRUE)
    ret<-createsubseqelist(seq,constraint,subseq[ord],data.frame(Support=(support[ord]/length(seq))),type="user")
    return(ret)
  }else{
    if(is.null(minSupport)){
      if(is.null(pMinSupport)){
        stop("You should specify a minimum support through minSupport or pMinSupport argument")
      }
      minSupport<-round(pMinSupport*length(seq))
    }
    classname<-c("seqe")
    subseq<-.Call("tmrfindsubsequences",unlist(list(seq)),as.double(c(constraint$maxGap)),as.double(c(constraint$windowSize)),
      as.double(c(constraint$ageMin)),as.double(c(constraint$ageMax)),as.double(c(constraint$ageMaxEnd)),as.double(c(constraint$countMethod)), as.integer(c(minSupport)),as.integer(c(maxK)),classname,PACKAGE="TraMineR")
    ord<-order(unlist(subseq[1]),decreasing=TRUE)

	
    ret<-createsubseqelist(seq,constraint,unlist(subseq[2])[ord],data.frame(Support=(unlist(subseq[1])[ord]/length(seq)), Count=(unlist(subseq[1])[ord])))
  }
  return(ret)
}