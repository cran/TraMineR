## ============================================
## Retrieve the alphabet from a sequence object
## ============================================

alphabet <- function(seqdata) {

	if (inherits(seqdata,c("stslist","PSTf"))){
    statl <- attr(seqdata,"alphabet")
  }
  else if (inherits(seqdata,"seqelist")){
    statl <- levels(seqdata)
  }
  else {
		stop("seqdata should be a state sequence object, an event sequence object, or a suffix tree. Use seqdef or seqecreate.")
  }

return(statl)
}
	
