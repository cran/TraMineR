
debuglevel <- function(level=NULL) {
	if(is.null(level)){
		return(.Call("getTraMineRDebugLevel",PACKAGE="TraMineR"))
	}
	.Call("setTraMineRDebugLevel",as.integer(level),PACKAGE="TraMineR")
	return(level)
	
}
