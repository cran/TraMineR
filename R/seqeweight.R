## ========================================
## Get and set weight of seqe
## ========================================

seqeweight<-function(s){
	seqeweight.internal<-function(s){
		if(is.seqe(s)) {
			return(.Call("tmrsequencegetweight",s, PACKAGE="TraMineR"))
		}
		return(-1)
	}
	
	if (is.seqelist(s)) {
		as.numeric(sapply(unlist(s),seqeweight.internal))
	}else if(is.seqe(s)) {
		as.numeric(seqeweight.internal(s))
	} else {
		stop(" [!] s should be a seqelist. See help on seqecreate.")
	}
}
"seqeweight<-" <- function(s, value){
	if(!is.seqelist(s)) {
		stop(" [!] s should be a seqelist. See help on seqecreate.")
	}
	if(length(s)!=length(value)) {
		stop(" [!] s and weights should be of the same size.")
	}
	.Call("tmrsequencesetweight", s, as.double(value), PACKAGE="TraMineR")
	return(s)
}

seqeisweighted <- function(s) {
	if(!is.seqelist(s)) {
		stop(" [!] s should be a seqelist. See help on seqecreate.")
	}
	weights <- seqeweight(s)
	return(any(weights!=1))
}