## ===================================
## Returns the length of the sequences 
## ===================================

seqlength <- function(seqdata) {

	if (!inherits(seqdata,"stslist")) {
		stop("data is not a sequence object, use 'seqdef' function to create one")
		}
	
	## missing <- attr(seqdata,"missing")
	## if (!is.na(missing)) seqdata[seqdata==missing] <- NA

	void <- attr(seqdata,"void")

	sl <- apply(seqdata,1, function(x) TraMineR.length(x,void))

	return(sl)
}
