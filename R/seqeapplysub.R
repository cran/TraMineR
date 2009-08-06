## ========================================
## Count number of time a subsequence appear in a sequence
## ========================================

seqeapplysub<-function(subseq, method="count",constraint=NULL, rules=FALSE){
#message("Event sequence analysis module is still experimental")
if(!is.subseqelist(subseq))stop("subseq should be a subseqelist. See help on seqefsub.")
  lastparam<-as.integer(c(1))
  if(method=="count")
    lastparam<-as.integer(c(1))
  else if(method=="presence")
    lastparam<-as.integer(c(2))
  else if(method=="age")
    lastparam<-as.integer(c(3))
  if(is.null(constraint)){
    constraint<-subseq$constraint
  }
  if(!inherits(constraint,"seqeconstraint")){
    constraint=seqeconstraint()
    warning("constraint argument should be set using seqeconstraint function. No constraint where used.")
  }
  if(!rules) {
  return(.Call("tmrmatrixsubseqinseq",unlist(list(subseq$subseq)),unlist(list(subseq$seqe)),
    as.double(c(constraint$maxGap)),as.double(c(constraint$windowSize)),
    as.double(c(constraint$ageMin)),as.double(c(constraint$ageMax)),as.double(c(constraint$ageMaxEnd)),
    lastparam,PACKAGE="TraMineR"))
  }
  else {
	  return(.Call("tmrmatrixsubseqinseq",unlist(list(subseq$subseq)),unlist(list(subseq$subseq)),
					  as.double(c(constraint$maxGap)),as.double(c(constraint$windowSize)),
					  as.double(c(constraint$ageMin)),as.double(c(constraint$ageMax)),as.double(c(constraint$ageMaxEnd)),
					  lastparam,PACKAGE="TraMineR"))
  }	
}
