TraMineRInternalNodeInit <- function(...){
	return(DTNInit(...))
}

TraMineRInternalSplitInit <- function(...){
	return(DTNsplit(...))
}

TraMineRInternalLayout <- function(...){
	return(TraMineR.setlayout(...))
}

TraMineRInternalSeqeage <- function(...){
	return(seqeage(...))
}

TraMineRInternalLegend <- function(...){
	return(TraMineR.legend(...))
}

TraMineRInternalSeqgbar <- function(...){
	return(seqgbar(...))
}

TraMineRInternalWeightedInertiaDist <- function(diss, diss.size, is.dist, individuals, sweights, var) {
  return(.Call(C_tmrWeightedInertiaDist, diss, diss.size, is.dist, individuals, sweights, var))
}
