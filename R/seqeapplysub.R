## ========================================
## Count number of time a subsequence appear in a sequence
## ========================================

seqeapplysub<-function(subseq, seq, method="count", maxGap=-1, windowSize=-1, ageMin=-1, ageMax=-1, ageMaxEnd=-1){
#message("Event sequence analysis module is still experimental")
if(!is.seqelist(seq))stop("seq should be a seqelist. See help on seqecreate.")
if(!is.seqelist(subseq))stop("subseq should be a seqelist. See help on seqecreate.")
  lastparam<-as.integer(c(1))
  if(method=="count")
    lastparam<-as.integer(c(1))
  else if(method=="presence")
    lastparam<-as.integer(c(2))
  else if(method=="age")
    lastparam<-as.integer(c(3))

  return(.Call("tmrmatrixsubseqinseq",unlist(list(subseq)),unlist(list(seq)),
    as.double(c(maxGap)),as.double(c(windowSize)),
    as.double(c(ageMin)),as.double(c(ageMax)),as.double(c(ageMaxEnd)),
    lastparam,PACKAGE="TraMineR"))
}
