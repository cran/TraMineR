## ========================================
## Find frequent subsequences
## ========================================

seqefsub <- function(seq, strsubseq=NULL, minSupport=NULL, pMinSupport=NULL, constraint=seqeconstraint(), maxK=-1, weighted=TRUE){
	if(!weighted) {
		ww <- seqeweight(seq)
		totseq <- length(seq)
		seqeweight(seq) <- as.double(rep(1, totseq))
	} else {
		totseq <- sum(seqeweight(seq))
	}
#  message("Event sequence analysis module is still experimental")
	if (!is.seqelist(seq)) {
		stop(" [!] seq should be a seqelist. See help on seqecreate.")
	}
    if (!inherits(constraint,"seqeconstraint")) {
		constraint=seqeconstraint()
		warning(" [!] constraint argument should be set using seqeconstraint function. No constraint where used.")
	}

  #We are not looking for frequent but specific subsequences
	if(!is.null(strsubseq)){
		subseq <- seqecreatesub(strsubseq, seq)
		ret <- createsubseqelist(seq, constraint, subseq, data.frame(), type="user")
		ww <- seqeweight(seq)
		## Taking weights into account
		support <- colSums(ww*seqeapplysub(ret, method="presence"))
		ord <- order(unlist(support), decreasing=TRUE)
		support <- support[ord]
		ret <- createsubseqelist(seq, constraint, subseq[ord], data.frame(Support=(support/totseq), Count=support),type="user")
		return(ret)
	}else{
		if(is.null(minSupport)){
			if(is.null(pMinSupport)){
				stop(" [!] You should specify a minimum support through minSupport or pMinSupport argument")
			}
			minSupport<-pMinSupport*sum(seqeweight(seq))
		}
		classname <- c("seqe")
		subseq <- .Call("tmrfindsubsequences", unlist(list(seq)), as.double(c(constraint$maxGap)), 
						as.double(c(constraint$windowSize)), as.double(c(constraint$ageMin)), 
						as.double(c(constraint$ageMax)), as.double(c(constraint$ageMaxEnd)), 
						as.double(c(constraint$countMethod)), as.double(c(minSupport)), 
						as.integer(c(maxK)), classname, PACKAGE="TraMineR")
		
		ord <- order(unlist(subseq[1]), decreasing=TRUE)
		support <- unlist(subseq[1])[ord]
		ret <- createsubseqelist(seq,constraint,unlist(subseq[2])[ord],data.frame(Support=(support/totseq), Count=support))
	}
	if(!weighted) {
		seqeweight(seq) <- ww
	}
	return(ret)
}